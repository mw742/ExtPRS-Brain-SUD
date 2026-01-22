#!/usr/bin/env bash
#SBATCH --job-name=extPRS_on_IDPs_v7
#SBATCH --output=extPRS_on_IDPs_v7_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --partition=highmem

# ----------------------- SAFETY -----------------------
set -euo pipefail
IFS=$'\n\t'
trap 'echo "[fatal] line $LINENO failed" >&2; exit 2' ERR

# ----------------------- ROOTS ------------------------
ROOT="/gpfs/group/ehlers/mwei/ABCD_2/ABCD"
WORKDIR="${ROOT}/work_provided_assets"

# ------------------------ INPUTS ----------------------
IDP_PARQUET="${WORKDIR}/gwas_output/ALL_BRAIN_BASELINE_NORM.parquet"
COVAR_TSV="${WORKDIR}/final_covariates.plusSES.tsv"
ICV_FILE="/gpfs/group/ehlers/mwei/Tables/core/imaging/icv_baseline.tsv"

# PRS input
PRSICE_FILE="/gpfs/group/ehlers/mwei/ABCD/QC_results_round3/PRS_results_perm_all.best"
PRSCS_FILE="/gpfs/group/ehlers/mwei/ABCD/QC_results_round3/PRS-CS/PRScs_individual_scores.profile"
PRS_KIND="${PRS_KIND:-PRSICE}"   # PRSICE or PRSCS

# ------------------------ OUTPUTS ---------------------
OUTDIR_BASE="${WORKDIR}/gwas_output/extPRS_on_IDPs"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTTAG="${OUTTAG:-modaware_v7_${STAMP}}"
OUTDIR="${OUTDIR_BASE}/${OUTTAG}"
mkdir -p "${OUTDIR}/per_idp" "${OUTDIR}/qc"

# ------------------- OPTIONS / FILTERS ----------------
WINSOR="${WINSOR:-0.005}"           # two-sided winsor p (0.005 => 0.5%)
MIN_N="${MIN_N:-200}"               # minimum N per IDP model
ROBUST_SE="${ROBUST_SE:-HC3}"       # HC3 or 'cluster' (set CLUSTER_BY to enable)
CLUSTER_BY="${CLUSTER_BY:-}"        # "", "site", "site_id", etc.
MOTION_COL="${MOTION_COL:-motion}"  # ABCD motion column
GLOBAL_COL="${GLOBAL_COL:-}"        # optional global (e.g., mean_thickness) -> sMRI only
KEEP_TXT="${KEEP_TXT:-}"            # keep-list (FID/IID or single IID)

# Optional IDP selection (default = all numeric)
IDP_LIST="${IDP_LIST:-}"            # file with one IDP/line
IDP_REGEX="${IDP_REGEX:-}"          # regex to include (case-insensitive OK)
IDP_EXCLUDE_REGEX="${IDP_EXCLUDE_REGEX:-^Unnamed}"  # default exclusions

# ---------------- SUPERVISOR-REQUESTED KNOBS ----------
SITE_COL="${SITE_COL:-site}"                             # site column in COVAR_TSV
SCANNER_MANU_COL="${SCANNER_MANU_COL:-scanner_manufacturer}"
SCANNER_SERIAL_COL="${SCANNER_SERIAL_COL:-scanner_device_serial}"
SES_COLS="${SES_COLS:-}"                                 # sensitivity follow-up (comma list)
QC_REGEX="${QC_REGEX:-mean_fd|fd|motion|qc}"             # extra dMRI/fMRI QC fields if present

# ------------------- DIAGNOSTICS ----------------------
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
echo "[diag] job_id=${SLURM_JOB_ID:-NA} cpus=$OMP_NUM_THREADS"
echo "[diag] IDP_PARQUET=$IDP_PARQUET"
echo "[diag] COVAR_TSV=$COVAR_TSV"
echo "[diag] ICV_FILE=$ICV_FILE"
echo "[diag] PRS_KIND=$PRS_KIND"
echo "[diag] OUTDIR=$OUTDIR"
echo "[diag] WINSOR=$WINSOR MIN_N=$MIN_N ROBUST_SE=$ROBUST_SE CLUSTER_BY=${CLUSTER_BY:-<none>}"
echo "[diag] MOTION_COL=${MOTION_COL:-<none>} GLOBAL_COL=${GLOBAL_COL:-<none>} KEEP_TXT=${KEEP_TXT:-<none>}"
echo "[diag] IDP_LIST=${IDP_LIST:-<none>} IDP_REGEX=${IDP_REGEX:-<none>} IDP_EXCLUDE_REGEX=${IDP_EXCLUDE_REGEX:-<none>}"
echo "[diag] SITE_COL=${SITE_COL:-<none>} MANU=${SCANNER_MANU_COL:-<none>} SERIAL=${SCANNER_SERIAL_COL:-<none>} QC_REGEX=${QC_REGEX:-<none>}"
echo "[diag] SES_COLS=${SES_COLS:-<none>}"

# ------------------- INPUT CHECKS ---------------------
[[ -f "$IDP_PARQUET" ]] || { echo "[fatal] Missing IDP_PARQUET"; exit 2; }
[[ -f "$COVAR_TSV"   ]] || { echo "[fatal] Missing COVAR_TSV"; exit 2; }
[[ -f "$ICV_FILE"    ]] || { echo "[fatal] Missing ICV_FILE"; exit 2; }

case "$PRS_KIND" in
  PRSICE) [[ -f "$PRSICE_FILE" ]] || { echo "[fatal] PRSICE_FILE not found"; exit 2; } ;;
  PRSCS)  [[ -f "$PRSCS_FILE"  ]] || { echo "[fatal] PRSCS_FILE not found"; exit 2; } ;;
  *) echo "[fatal] PRS_KIND must be PRSICE or PRSCS"; exit 2;;
esac

if [[ -n "$KEEP_TXT" && ! -f "$KEEP_TXT" ]]; then
  echo "[fatal] KEEP_TXT specified but file not found: ${KEEP_TXT}"; exit 2
fi
if [[ -n "$IDP_LIST" && ! -f "$IDP_LIST" ]]; then
  echo "[fatal] IDP_LIST specified but file not found"; exit 2
fi

# ----------------- ENV / PYTHON CHOICE ----------------
unset PYTHONHOME PYTHONPATH PYTHONNOUSERSITE

# ------------------- EXPORT FOR PY --------------------
export IDP_PARQUET COVAR_TSV ICV_FILE OUTDIR
export PRS_KIND PRSICE_FILE PRSCS_FILE
export WINSOR MIN_N ROBUST_SE CLUSTER_BY MOTION_COL GLOBAL_COL KEEP_TXT
export IDP_LIST IDP_REGEX IDP_EXCLUDE_REGEX
export SITE_COL SCANNER_MANU_COL SCANNER_SERIAL_COL SES_COLS QC_REGEX

# --------------------- PYTHON BLOCK -------------------
micromamba run -n my_env python - <<'PY'
import os, re
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

# ------------ helpers ------------
def need(k):
    v=os.environ.get(k)
    if not v: raise SystemExit(f"[fatal] Missing env: {k}")
    return v

def canon_iid(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().str.replace(r"^SUB-","",regex=True)

def winsor(s, p):
    if p<=0: return pd.to_numeric(s, errors="coerce")
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.quantile(p), s.quantile(1-p)
    return s.clip(lo, hi)

def has_var(s):
    s = pd.to_numeric(s, errors="coerce").dropna()
    return (len(s) > 1) and (s.std(ddof=0) > 0)

def bh(p):
    p = np.asarray(p, float); m = len(p)
    if m == 0: return p
    o = np.argsort(p); r = p[o]
    q = r*m/np.arange(1,m+1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(p); out[o] = np.clip(q,0,1)
    return out

def ensure_numeric_or_dummies(df: pd.DataFrame, col: str, prefix: str):
    """Return (df2, added_columns). If col is non-numeric, one-hot encode with k-1 dummies."""
    if (not col) or (col not in df.columns):
        return df, []
    if pd.api.types.is_numeric_dtype(df[col]):
        return df, [col]
    tmp = pd.to_numeric(df[col], errors="coerce")
    if tmp.notna().all():
        df[col] = tmp
        return df, [col]
    cat = df[col].astype("category")
    d = pd.get_dummies(cat, prefix=prefix, drop_first=True, dummy_na=False)
    df = pd.concat([df.drop(columns=[col]), d], axis=1)
    return df, list(d.columns)

# ------------ env ------------
IDP_PARQUET = need("IDP_PARQUET")
COVAR_TSV   = need("COVAR_TSV")
ICV_FILE    = need("ICV_FILE")
OUTDIR      = need("OUTDIR")
PRS_KIND    = need("PRS_KIND").upper()
PRSICE_FILE = need("PRSICE_FILE")
PRSCS_FILE  = need("PRSCS_FILE")

WINSOR      = float(os.environ.get("WINSOR","0"))
MIN_N       = int(os.environ.get("MIN_N","200"))
ROBUST_SE   = os.environ.get("ROBUST_SE","HC3").upper()
CLUSTER_BY  = os.environ.get("CLUSTER_BY","").strip()
MOTION_COL  = os.environ.get("MOTION_COL","").strip()
GLOBAL_COL  = os.environ.get("GLOBAL_COL","").strip()
KEEP_TXT    = os.environ.get("KEEP_TXT","").strip()

IDP_LIST    = os.environ.get("IDP_LIST","").strip()
IDP_REGEX   = os.environ.get("IDP_REGEX","").strip()
IDP_EXCLUDE_REGEX = os.environ.get("IDP_EXCLUDE_REGEX","").strip()

SITE_COL           = os.environ.get("SITE_COL","").strip()
SCANNER_MANU_COL   = os.environ.get("SCANNER_MANU_COL","").strip()
SCANNER_SERIAL_COL = os.environ.get("SCANNER_SERIAL_COL","").strip()
SES_COLS           = [c.strip() for c in os.environ.get("SES_COLS","").split(",") if c.strip()]
QC_REGEX           = os.environ.get("QC_REGEX","").strip()

OUTDIR_PATH = Path(OUTDIR); (OUTDIR_PATH/"per_idp").mkdir(parents=True, exist_ok=True); (OUTDIR_PATH/"qc").mkdir(exist_ok=True)

# ------------ load tables ------------
idps = pd.read_parquet(IDP_PARQUET)
pd.DataFrame({"col": idps.columns, "dtype": [str(t) for t in idps.dtypes]}).to_csv(OUTDIR_PATH/"qc"/"columns_snapshot.tsv", sep="\t", index=False)

# IID unify
for cand in ["IID","src_subject_id","subjectkey","iid","#IID","id","SUBJECT","participant_id","FID"]:
    if cand in idps.columns:
        idps = idps.rename(columns={cand:"IID"})
        break
if "IID" not in idps.columns: raise SystemExit("[fatal] No IID column in IDP_PARQUET")
idps["IID"] = canon_iid(idps["IID"])

cov = pd.read_csv(COVAR_TSV, sep="\t")
for cand in ["IID","src_subject_id","subjectkey","iid","#IID","id","SUBJECT","participant_id","FID"]:
    if cand in cov.columns:
        cov = cov.rename(columns={cand:"IID"})
        break
cov["IID"] = canon_iid(cov["IID"])

# normalize sex to {0,1}
if "sex" in cov.columns:
    if pd.api.types.is_numeric_dtype(cov["sex"]):
        cov["sex"] = cov["sex"].replace({1:0,2:1,0:np.nan})
    else:
        cov["sex"] = cov["sex"].map({"M":0,"F":1,"MALE":0,"FEMALE":1,"male":0,"female":1})
for req in ["age","sex"]:
    if req not in cov.columns: raise SystemExit(f"[fatal] Covariates missing {req}")

pcs = [c for c in cov.columns if c.upper().startswith("PC")]

# ICV
icv = pd.read_csv(ICV_FILE, sep=None, engine="python")
for cand in ["IID","src_subject_id","subjectkey","iid","#IID","id","SUBJECT","participant_id","FID"]:
    if cand in icv.columns:
        icv = icv.rename(columns={cand:"IID"})
        break
icv["IID"] = canon_iid(icv["IID"])
icv_col=None
for c in ["ICV","eTIV","E_TIV","E_TIV_mm3","intracranial_volume","icv","etiv"]:
    if c in icv.columns: icv_col=c; break
if icv_col is None:
    icv_col = icv.columns[1] if icv.shape[1]>=2 else None
if icv_col is None: raise SystemExit("[fatal] Could not detect ICV column")
icv = icv[["IID", icv_col]].rename(columns={icv_col:"ICV"})

# PRS
prs_path = PRSICE_FILE if PRS_KIND=="PRSICE" else PRSCS_FILE
prs = pd.read_csv(prs_path, sep=r"\s+|\t|,", engine="python")
if "IID" in prs.columns:
    iid_series = prs["IID"].astype(str).str.strip()
else:
    patt = re.compile(r"^(SUB-)?NDAR_[A-Z0-9]{8,}$")
    candidates = [(c, prs[c].astype(str).str.match(patt).sum()) for c in prs.columns if prs[c].dtype==object]
    candidates.sort(key=lambda x:x[1], reverse=True)
    iid_series = prs[candidates[0][0]] if candidates and candidates[0][1]>0 else prs.iloc[:,0]
prs["IID"] = canon_iid(iid_series)
score_col = "SCORESUM" if "SCORESUM" in prs.columns else (
    next((c for c in ["PRS","PGS","SCORE","SCORE1","SCORE1_SUM","PRS_SCORE"] if c in prs.columns), None)
    or prs.select_dtypes("number").columns[-1]
)
prs = prs[["IID", score_col]].rename(columns={score_col:"ext_PRS"})
prs["ext_PRS"] = (prs["ext_PRS"] - prs["ext_PRS"].mean()) / (prs["ext_PRS"].std(ddof=0) or 1.0)

# keep-list
KEEP_TXT = os.environ.get("KEEP_TXT","").strip()
if KEEP_TXT:
    keep = pd.read_csv(KEEP_TXT, sep=r"\s+|\t|,", header=None, engine="python")
    keep = (keep.iloc[:,1] if keep.shape[1]>=2 else keep.iloc[:,0]).astype(str)
    keep = set(canon_iid(keep))
    idps = idps[idps["IID"].isin(keep)]
    cov  =  cov[cov ["IID"].isin(keep)]
    icv  =  icv[icv ["IID"].isin(keep)]
    prs  =  prs[prs ["IID"].isin(keep)]

# intersect
S = set(idps["IID"]) & set(cov["IID"]) & set(icv["IID"]) & set(prs["IID"])
if len(S) < MIN_N: raise SystemExit(f"[fatal] Too few overlapping IIDs: {len(S)}")
idps = idps[idps["IID"].isin(S)].copy()
cov  =  cov[cov["IID"].isin(S)].copy()
icv  =  icv[icv["IID"].isin(S)].copy()
prs  =  prs[prs["IID"].isin(S)].copy()

# ------------ base covariate block ------------
base_cov_cols = ["IID","age","sex"] + pcs
for extra in (SITE_COL, SCANNER_MANU_COL, SCANNER_SERIAL_COL):
    if extra and extra in cov.columns and extra not in base_cov_cols:
        base_cov_cols.append(extra)

# SES (sensitivity)
for col in [c.strip() for c in os.environ.get("SES_COLS","").split(",") if c.strip()]:
    if col in cov.columns and col not in base_cov_cols:
        base_cov_cols.append(col)

covblock = cov[base_cov_cols].merge(prs[["IID","ext_PRS"]], on="IID", how="inner")
covblock = covblock.merge(icv, on="IID", how="left")

# ---- categorical → dummies ----
qc_dir = OUTDIR_PATH / "qc"
added_all = []
covblock, added_site   = ensure_numeric_or_dummies(covblock, SITE_COL, "site");           added_all += [("site", c) for c in added_site]
covblock, added_manu   = ensure_numeric_or_dummies(covblock, SCANNER_MANU_COL, "scanner_manu"); added_all += [("scanner_manu", c) for c in added_manu]
covblock, added_serial = ensure_numeric_or_dummies(covblock, SCANNER_SERIAL_COL, "scanner_serial"); added_all += [("scanner_serial", c) for c in added_serial]

# QC snapshots (fixed unequal-length lists)
pd.DataFrame({"col": covblock.columns}).to_csv(qc_dir/"covblock_columns.tsv", sep="\t", index=False)
pd.DataFrame(added_all, columns=["source","dummy_col"]).to_csv(qc_dir/"categorical_to_dummies.tsv", sep="\t", index=False)

# Cluster group counts
if os.environ.get("CLUSTER_BY","") and (os.environ["CLUSTER_BY"] in cov.columns):
    cov.groupby(os.environ["CLUSTER_BY"], dropna=True).size().to_csv(
        qc_dir/f"cluster_group_counts_{os.environ['CLUSTER_BY']}.tsv", sep="\t"
    )

# ------------ IDP selection ------------
def is_numeric_col(df, c):
    try:
        return pd.api.types.is_numeric_dtype(df[c])
    except Exception:
        return False

numeric_cols = [c for c in idps.columns if c!="IID" and is_numeric_col(idps, c)]
wl = numeric_cols[:]
IDP_LIST = os.environ.get("IDP_LIST","").strip()
IDP_REGEX = os.environ.get("IDP_REGEX","").strip()
IDP_EXCLUDE_REGEX = os.environ.get("IDP_EXCLUDE_REGEX","").strip()

if IDP_LIST:
    f = pd.read_csv(IDP_LIST, header=None)[0].astype(str).str.strip().tolist()
    wl = [w for w in f if w in numeric_cols]
if IDP_REGEX:
    rx = re.compile(IDP_REGEX, flags=re.IGNORECASE)
    wl = [c for c in wl if rx.search(c)]
if IDP_EXCLUDE_REGEX:
    ex = re.compile(IDP_EXCLUDE_REGEX)
    wl = [c for c in wl if not ex.search(c)]
pd.Series(wl, name="idp").to_csv(qc_dir/"idp_working_set.txt", index=False)

# ------------ covariate chooser per IDP ------------
SES_COLS = [c.strip() for c in os.environ.get("SES_COLS","").split(",") if c.strip()]
def covars_for(colname: str):
    name = colname.lower()
    pcs_cols = [c for c in covblock.columns if c.upper().startswith("PC")]
    base = ["age","sex"] + pcs_cols
    base += [c for c in covblock.columns if c.startswith("site_")]
    base += [c for c in covblock.columns if c.startswith("scanner_manu_")]
    base += [c for c in covblock.columns if c.startswith("scanner_serial_")]
    for s in SES_COLS:
        if s in covblock.columns: base.append(s)

    is_non_t1    = any(k in name for k in ["rsi","_dt","fmri","dti","dmri"])
    is_thickness = ("thick" in name) or ("thickness" in name)
    is_sulc      = ("sulc" in name) or ("sulcal" in name)
    is_area      = ("area" in name)
    is_volume    = ("vol" in name) or ("volume" in name)

    if is_non_t1:
        MOTION_COL = os.environ.get("MOTION_COL","").strip()
        if MOTION_COL and (MOTION_COL in covblock.columns): base.append(MOTION_COL)
        QC_REGEX = os.environ.get("QC_REGEX","").strip()
        if QC_REGEX:
            rx = re.compile(QC_REGEX, flags=re.IGNORECASE)
            for q in [c for c in covblock.columns if rx.search(c)]:
                if q not in base: base.append(q)
    else:
        if is_area or is_volume:
            if "ICV" in covblock.columns: base.append("ICV")
        GLOBAL_COL = os.environ.get("GLOBAL_COL","").strip()
        if GLOBAL_COL and (GLOBAL_COL in covblock.columns): base.append(GLOBAL_COL)

    return base

# ------------ modeling ------------
rows, skipped = [], []
for col in wl:
    if col not in idps.columns:
        skipped.append((col, "missing_from_parquet"))
        continue

    model_covs = covars_for(col)
    X_cols = ["ext_PRS"] + model_covs

    need = ["IID", col] + [c for c in X_cols if c in covblock.columns]
    CLUSTER_BY = os.environ.get("CLUSTER_BY","").strip()
    want_cluster = (ROBUST_SE == "CLUSTER") and bool(CLUSTER_BY)

    dat = idps[["IID", col]].merge(covblock[need[2:]+["IID"]], on="IID", how="inner")
    if want_cluster and (CLUSTER_BY in cov.columns):
        dat = dat.merge(cov[["IID", CLUSTER_BY]], on="IID", how="left")

    dat = dat.dropna()
    if len(dat) < MIN_N:
        skipped.append((col, f"too_small_N({len(dat)})"))
        continue

    y = pd.to_numeric(dat[col], errors="coerce")
    WINSOR = float(os.environ.get("WINSOR","0"))
    if WINSOR>0: y = winsor(y, WINSOR)
    y = (y - y.mean()) / (y.std(ddof=0) or 1.0)
    if not has_var(y):
        skipped.append((col, "zero_var_y"))
        continue

    Xf = dat[X_cols].apply(pd.to_numeric, errors="coerce")
    nunique = Xf.nunique(dropna=True)
    keep_covs = ["ext_PRS"] + [c for c in X_cols[1:] if nunique.get(c,2)>1]
    Xf = Xf[keep_covs]
    if not has_var(Xf["ext_PRS"]):
        skipped.append((col, "zero_var_prs"))
        continue

    Xf_const = sm.add_constant(Xf, has_constant="add")

    groups = None
    cov_type = "HC3"
    cov_kwds = None
    if want_cluster:
        if CLUSTER_BY not in dat.columns:
            print(f"[warn] CLUSTER_BY='{CLUSTER_BY}' not found; fallback HC3.")
        else:
            g = pd.Series(dat[CLUSTER_BY]).astype("category")
            n_groups = g.cat.categories.size if g.notna().any() else 0
            if n_groups < 2:
                print(f"[warn] CLUSTER_BY='{CLUSTER_BY}' <2 groups; fallback HC3.")
            else:
                groups = g.cat.codes.values
                cov_type = "cluster"
                cov_kwds = {"groups": groups}

    try:
        fitf = sm.OLS(y.values, Xf_const.values).fit(cov_type=cov_type, cov_kwds=cov_kwds)
    except Exception as e:
        print(f"[warn] {cov_type} fit failed ({e}); fallback HC3.")
        fitf = sm.OLS(y.values, Xf_const.values).fit(cov_type="HC3")

    Xr = dat[[c for c in Xf.columns if c!="ext_PRS"]]
    Xr_const = sm.add_constant(Xr, has_constant="add")
    try:
        fitr = sm.OLS(y.values, Xr_const.values).fit(cov_type=cov_type, cov_kwds=cov_kwds)
    except Exception:
        fitr = sm.OLS(y.values, Xr_const.values).fit(cov_type="HC3")

    names = ["const"] + list(Xf.columns)
    pmap  = dict(zip(names, fitf.pvalues))
    bmap  = dict(zip(names, fitf.params))
    semap = dict(zip(names, fitf.bse))
    if "ext_PRS" not in pmap:
        skipped.append((col, "no_ext_PRS_in_model"))
        continue

    R2f, R2r = float(fitf.rsquared), float(fitr.rsquared)
    partial = max(0.0, (R2f - R2r) / (1.0 - R2r + 1e-12))

    rows.append({
        "idp": col, "n": int(len(dat)),
        "beta": float(bmap["ext_PRS"]),
        "se":   float(semap["ext_PRS"]),
        "p":    float(pmap["ext_PRS"]),
        "r2":   R2f, "r2_reduced": R2r, "partial_r2": partial,
        "used_motion": int(os.environ.get("MOTION_COL","").strip() in keep_covs),
        "used_icv":    int("ICV" in keep_covs),
        "covars": "|".join([c for c in keep_covs if c!="ext_PRS"])
    })

    with open(OUTDIR_PATH/"per_idp"/f"{col}.coef.txt","w") as f:
        f.write(f"idp\t{col}\nn\t{len(dat)}\n")
        f.write(f"beta\t{bmap['ext_PRS']}\nse\t{semap['ext_PRS']}\np\t{pmap['ext_PRS']}\n")
        f.write(f"r2\t{R2f}\nr2_reduced\t{R2r}\npartial_r2\t{partial}\n")
        f.write(f"covars\t{'|'.join(keep_covs)}\n")

pd.DataFrame(skipped, columns=["idp","reason"]).to_csv(OUTDIR_PATH/"qc"/"skipped_idps.tsv", sep="\t", index=False)

# ------------ aggregate + BH/Bonf ------------
res_cols = ["idp","n","beta","se","p","r2","r2_reduced","partial_r2","used_motion","used_icv","covars"]
res = pd.DataFrame(rows, columns=res_cols) if rows else pd.DataFrame(columns=res_cols)
if res.empty:
    print("[warn] No models ran. See qc/columns_snapshot.tsv and qc/skipped_idps.tsv.")
    raise SystemExit(0)

res = res.sort_values("p").reset_index(drop=True)
res["p_BH"] = bh(res["p"].values)
m = res.shape[0]
res["p_Bonf"] = np.minimum(res["p"] * max(m,1), 1.0)

res.to_csv(OUTDIR_PATH/"extPRS_on_IDPs_summary.tsv", sep="\t", index=False)
print(f"[ok] wrote {OUTDIR_PATH/'extPRS_on_IDPs_summary.tsv'} rows={len(res)}")
PY

