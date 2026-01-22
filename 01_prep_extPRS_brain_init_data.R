#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

canon_iid <- function(x) {
  x %>% as.character() %>% str_trim() %>% toupper() %>% str_replace("^SUB-", "")
}

# ---------------- args ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if (length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

ROOT   <- get_arg("--ROOT", "/gpfs/group/ehlers/mwei/ABCD_2/ABCD")
WORK   <- file.path(ROOT, "work_provided_assets")

IDP_CSV    <- get_arg("--IDP_CSV", file.path(WORK, "gwas_output/ALL_BRAIN_BASELINE_NORM.csv.gz"))
COVAR_TSV  <- get_arg("--COVAR_TSV", file.path(WORK, "final_covariates.plusSES.tsv"))
PRS_FILE   <- get_arg("--PRS_FILE", "/gpfs/group/ehlers/mwei/ABCD/QC_results_round3/PRS_results_perm_all.best")
INIT_DIR   <- get_arg("--INIT_DIR", "/gpfs/group/ehlers/mwei/ABCD/init_phenos")

HITS_FILE  <- get_arg("--HITS_FILE", NA_character_)
OUTDIR     <- get_arg("--OUTDIR", NA_character_)

if (is.na(HITS_FILE) || !file.exists(HITS_FILE)) stop("[fatal] HITS_FILE missing or not found: ", HITS_FILE)
if (is.na(OUTDIR) || OUTDIR=="") stop("[fatal] OUTDIR missing")

dir.create(file.path(OUTDIR, "data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR, "qc"), recursive=TRUE, showWarnings=FALSE)

cand_iid <- c("IID","src_subject_id","subjectkey","iid","#IID",
              "id","SUBJECT","participant_id","FID")

# ---------------- 1) read significant IDPs ----------------
sig_tbl <- fread(HITS_FILE)
if (!"idp" %in% names(sig_tbl)) setnames(sig_tbl, 1, "idp")
sig_tbl[, idp := str_trim(as.character(idp))]
sig_idps_new <- unique(sig_tbl$idp)
cat("[diag] significant IDPs from HITS_FILE:", length(sig_idps_new), "\n")

# ---------------- 2) read IDP header only, build mapping ----------------
cat("[diag] reading IDP header:", IDP_CSV, "\n")
idp_header <- fread(IDP_CSV, nrows = 0)
idp_cols_all <- names(idp_header)

iid_col <- NULL
for (c in cand_iid) if (c %in% idp_cols_all) { iid_col <- c; break }
if (is.null(iid_col)) stop("[fatal] No IID-like column found in IDP_CSV header")

old_cols <- setdiff(idp_cols_all, iid_col)

clean_fun <- function(x) sub(".*__", "", as.character(x))
clean_names <- clean_fun(old_cols)

if (any(duplicated(clean_names))) {
  dup <- unique(clean_names[duplicated(clean_names)])
  cat("[fatal] duplicate cleaned IDP names after '__' stripping. Example:\n")
  print(head(dup, 20))
  stop("Cleaned IDP names are not unique. Need stricter mapping rule.")
}

idp_map <- data.table(idp_old = old_cols, idp_clean = clean_names)

map_merged <- merge(
  data.table(idp = sig_idps_new),
  idp_map,
  by.x = "idp",
  by.y = "idp_clean",
  all.x = TRUE,
  sort = FALSE
)

matched <- map_merged[!is.na(idp_old)]
cat("[diag] matched LiJi IDs:", nrow(matched), "of", length(sig_idps_new), "\n")

if (nrow(matched) == 0) {
  cat("[fatal] example LiJi idp names:\n"); print(head(sig_idps_new, 10))
  cat("[fatal] example cleaned IDP column names:\n"); print(head(idp_map$idp_clean, 20))
  stop("No overlapping IDPs after cleaning rule.")
}

use_old <- unique(matched$idp_old)
use_new <- matched$idp[match(use_old, matched$idp_old)]

# ---------------- 3) read only needed IDP columns ----------------
select_cols <- c(iid_col, use_old)
cat("[diag] reading IDP subset columns:", length(select_cols), "\n")
idps <- fread(IDP_CSV, select = select_cols, showProgress = FALSE)
setnames(idps, iid_col, "IID")
idps[, IID := canon_iid(IID)]
setnames(idps, old = use_old, new = use_new)

# ---------------- 4) read covariates ----------------
cov <- fread(COVAR_TSV)
for (c in cand_iid) if (c %in% names(cov)) { setnames(cov, c, "IID"); break }
if (!"IID" %in% names(cov)) stop("[fatal] No IID column in covariates")
cov[, IID := canon_iid(IID)]

# ---------------- 5) read PRS ----------------
prs <- fread(PRS_FILE)
for (c in cand_iid) if (c %in% names(prs)) { setnames(prs, c, "IID"); break }
if (!"IID" %in% names(prs)) stop("[fatal] No IID column in PRS")
prs[, IID := canon_iid(IID)]

prs_cols <- setdiff(names(prs), c("FID","IID"))
if (length(prs_cols) < 1) stop("[fatal] No PRS columns found")

numeric_score_info <- rbindlist(lapply(prs_cols, function(col) {
  v <- suppressWarnings(as.numeric(prs[[col]]))
  data.table(col = col, n_numeric = sum(!is.na(v)), n_total = length(v))
}))
setorder(numeric_score_info, -n_numeric)
score_col <- numeric_score_info$col[1]
cat("[diag] using PRS column:", score_col, "\n")

setnames(prs, score_col, "ext_PRS")
prs[, ext_PRS := suppressWarnings(as.numeric(ext_PRS))]
mu <- mean(prs$ext_PRS, na.rm=TRUE)
sdv <- sd(prs$ext_PRS, na.rm=TRUE); if (is.na(sdv) || sdv==0) sdv <- 1
prs[, ext_PRS := (ext_PRS - mu)/sdv]
prs <- prs[, .(IID, ext_PRS)]

# ---------------- 6) prep one outcome ----------------
prep_one <- function(outcome_name) {
  surv_file <- file.path(INIT_DIR, paste0(outcome_name, "_survival.tsv"))
  if (!file.exists(surv_file)) stop("[fatal] missing survival file: ", surv_file)

  cat("[diag] reading survival:", surv_file, "\n")
  surv <- fread(surv_file)
  for (c in cand_iid) if (c %in% names(surv)) { setnames(surv, c, "IID"); break }
  if (!"IID" %in% names(surv)) stop("[fatal] no IID column in ", surv_file)
  surv[, IID := canon_iid(IID)]

  common <- Reduce(intersect, list(idps$IID, cov$IID, prs$IID, surv$IID))
  cat("[diag]", outcome_name, "common N =", length(common), "\n")

  dat <- surv[IID %in% common] %>%
    left_join(prs, by="IID") %>%
    left_join(cov, by="IID") %>%
    left_join(idps, by="IID") %>%
    dplyr::filter(!is.na(time), !is.na(event), !is.na(ext_PRS))

  out_file <- file.path(OUTDIR, "data", paste0("extPRS_brain_", outcome_name, ".tsv"))
  fwrite(as.data.table(dat), out_file, sep="\t")
  cat("[ok] wrote", out_file, "N=", nrow(dat), "\n")
}

for (nm in c("alcohol","nicotine","cannabis","any")) prep_one(nm)

fwrite(data.table(idp = use_new), file.path(OUTDIR, "qc", "sig_idps_used.tsv"), sep="\t")
cat("[ok] prep finished. OUTDIR=", OUTDIR, "\n")

