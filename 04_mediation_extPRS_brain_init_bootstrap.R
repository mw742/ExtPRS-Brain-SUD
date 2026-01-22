#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mediation)
})

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)

# -------------------------
# robust arg parser: supports
#   --FLAG=value   and   --FLAG value
# -------------------------
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  w <- which(args == flag)
  if (length(w) > 0) {
    i <- w[1] + 1
    if (i <= length(args)) return(args[i])
  }
  default
}

OUTDIR   <- get_arg("--OUTDIR", NA_character_)
SE_MODE  <- tolower(get_arg("--SE_MODE", "cluster"))  # reporting only
SES_MODE <- tolower(get_arg("--SES_MODE", "off"))
SIMS     <- as.integer(get_arg("--SIMS", "5000"))

if (is.na(OUTDIR) || OUTDIR == "") stop("[fatal] OUTDIR missing")

# -------------------------
# paths
# -------------------------
DATA_DIR   <- file.path(OUTDIR, "data")
STAGE4_DIR <- file.path(OUTDIR, "stage4")
dir.create(STAGE4_DIR, recursive = TRUE, showWarnings = FALSE)

CAND_FILE    <- file.path(OUTDIR, "stage3", "extPRS_IDPs_on_initiation_signif.tsv")
RAW_FILE     <- file.path(STAGE4_DIR, "mediation_extPRS_brain_initiation_raw.tsv")
SUMMARY_FILE <- file.path(STAGE4_DIR, "mediation_extPRS_brain_initiation_summary.tsv")
DONE_FILE    <- file.path(STAGE4_DIR, "_STAGE4_DONE.txt")

# SES covariates (only used if present in data)
ses_covars <- c(
  "income_num_z","parent_edu_years_z","parent_married_or_cohab_z",
  "income_num_z_missing","parent_edu_years_z_missing","parent_married_or_cohab_z_missing"
)

write_empty_outputs <- function(msg) {
  if (!file.exists(RAW_FILE)) {
    empty_raw <- data.table(
      outcome=character(), idp=character(), n=integer(), n_case=integer(), n_ctrl=integer(),
      acme=double(), acme_p=double(), ade=double(), ade_p=double(),
      prop_med=double(), prop_med_p=double(),
      se_mode=character(), ses_mode=character(), sims=integer()
    )
    fwrite(empty_raw, RAW_FILE, sep = "\t")
  }
  fwrite(data.table(message = msg, se_mode = SE_MODE, ses_mode = SES_MODE, sims = SIMS),
         SUMMARY_FILE, sep = "\t")
  writeLines("OK", DONE_FILE)
  cat("[ok] stage4 wrote empty outputs:", msg, "\n")
}

# -------------------------
# candidate pairs
# -------------------------
if (!file.exists(CAND_FILE)) stop("[fatal] missing candidate file: ", CAND_FILE)
cands <- tryCatch(fread(CAND_FILE), error = function(e) NULL)
if (is.null(cands) || nrow(cands) == 0) {
  write_empty_outputs("Candidate file exists but is empty")
  quit(save = "no", status = 0)
}
if (!all(c("outcome", "idp") %in% names(cands))) stop("[fatal] candidate file needs columns: outcome,idp")

cands <- unique(cands[, .(outcome, idp)])
cands[, key := paste(outcome, idp, sep = "__")]

# -------------------------
# resume: skip already-done pairs
# -------------------------
done_keys <- character(0)
if (file.exists(RAW_FILE) && file.info(RAW_FILE)$size > 0) {
  done <- tryCatch(fread(RAW_FILE, select = c("outcome", "idp")), error = function(e) NULL)
  if (!is.null(done) && nrow(done) > 0) {
    done <- unique(done)
    done_keys <- paste(done$outcome, done$idp, sep = "__")
  }
}
todo <- cands[!(key %in% done_keys)]
cat("[diag] total pairs:", nrow(cands), " already done:", length(done_keys), " todo:", nrow(todo), "\n")

if (nrow(todo) == 0L) {
  fwrite(
    data.table(message = "Nothing to do (all pairs already completed)",
               n_total_pairs = nrow(cands), n_done_pairs_preexisting = length(done_keys),
               n_todo_pairs = 0, se_mode = SE_MODE, ses_mode = SES_MODE, sims = SIMS),
    SUMMARY_FILE, sep = "\t"
  )
  writeLines("OK", DONE_FILE)
  cat("[ok] nothing to do; done marker written:", DONE_FILE, "\n")
  quit(save = "no", status = 0)
}

# -------------------------
# mediate wrapper with fallback
# -------------------------
safe_mediate <- function(med.fit, out.fit, treat, mediator, sims_n) {
  attempt <- function(boot_flag, sims_try, seed) {
    set.seed(seed)
    tryCatch(
      mediate(med.fit, out.fit, treat = treat, mediator = mediator, boot = boot_flag, sims = sims_try),
      error = function(e) NULL
    )
  }
  res <- attempt(TRUE, sims_n, 2025)
  if (!is.null(res)) return(res)
  res <- attempt(TRUE, min(1000, sims_n), 2026)
  if (!is.null(res)) return(res)
  attempt(FALSE, min(1000, sims_n), 2027)
}

# -------------------------
# diagnostics counters
# -------------------------
diag <- list(
  n_pairs = nrow(todo),
  missing_dat_file = 0L,
  header_read_fail = 0L,
  missing_event_or_prs = 0L,
  idp_missing_in_data = 0L,
  fread_fail = 0L,
  qc_fail = 0L,
  lm_fail = 0L,
  glm_fail = 0L,
  mediate_fail = 0L,
  written_ok = 0L
)

# -------------------------
# covariates (Stage4)
# -------------------------
covars_base <- c("age", "sex", paste0("PC", 1:20), "site", "scanner_manufacturer")
if (SES_MODE == "on") covars_base <- c(covars_base, ses_covars)

# -------------------------
# main loop: process by outcome (read each outcome file ONCE)
# -------------------------
first_write <- !(file.exists(RAW_FILE) && file.info(RAW_FILE)$size > 0)

outcomes <- sort(unique(todo$outcome))
for (oc in outcomes) {
  dat_file <- file.path(DATA_DIR, paste0("extPRS_brain_", oc, ".tsv"))
  if (!file.exists(dat_file)) {
    diag$missing_dat_file <- diag$missing_dat_file + sum(todo$outcome == oc)
    next
  }

  # read header only to figure out available columns
  hdr <- tryCatch(fread(dat_file, nrows = 0), error = function(e) NULL)
  if (is.null(hdr)) {
    diag$header_read_fail <- diag$header_read_fail + sum(todo$outcome == oc)
    next
  }
  cols_all <- names(hdr)

  event_col <- "event"
  prs_col   <- "ext_PRS"
  if (!(event_col %in% cols_all) || !(prs_col %in% cols_all)) {
    diag$missing_event_or_prs <- diag$missing_event_or_prs + sum(todo$outcome == oc)
    next
  }

  # IDPs to test for this outcome
  idps <- todo[outcome == oc, idp]
  idps_present <- idps[idps %in% cols_all]
  diag$idp_missing_in_data <- diag$idp_missing_in_data + (length(idps) - length(idps_present))
  if (length(idps_present) == 0) next

  covars_present <- covars_base[covars_base %in% cols_all]
  select_cols <- unique(c(event_col, prs_col, covars_present, idps_present))

  dat <- tryCatch(fread(dat_file, select = select_cols), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) {
    diag$fread_fail <- diag$fread_fail + length(idps_present)
    next
  }

  # prepare covars for formula (only those actually loaded)
  covars <- covars_present

  # enforce factors for categorical vars if present
  if ("site" %in% names(dat)) dat[, site := as.factor(site)]
  if ("scanner_manufacturer" %in% names(dat)) dat[, scanner_manufacturer := as.factor(scanner_manufacturer)]

  # loop IDPs within this outcome
  for (idp_name in idps_present) {
    cat("[diag] outcome", oc, "idp", idp_name, "\n")

    # subset minimal columns to reduce overhead
    needed <- c(event_col, prs_col, covars, idp_name)
    dt <- dat[, ..needed]

    dt2 <- dt %>%
      filter(!is.na(.data[[idp_name]]), !is.na(.data[[prs_col]]), !is.na(.data[[event_col]])) %>%
      mutate(
        event     = as.integer(.data[[event_col]]),
        ext_raw   = as.numeric(.data[[prs_col]]),
        ext_PRS_z = as.numeric(scale(ext_raw)[, 1]),
        idp_z     = as.numeric(scale(as.numeric(.data[[idp_name]]))[, 1])
      )

    n <- nrow(dt2)
    if (n < 200) { diag$qc_fail <- diag$qc_fail + 1L; next }

    n_case <- sum(dt2$event == 1, na.rm = TRUE)
    n_ctrl <- sum(dt2$event == 0, na.rm = TRUE)
    if (n_case < 20 || n_ctrl < 20) { diag$qc_fail <- diag$qc_fail + 1L; next }

    if (sd(dt2$ext_PRS_z, na.rm = TRUE) < 1e-6 || sd(dt2$idp_z, na.rm = TRUE) < 1e-6) {
      diag$qc_fail <- diag$qc_fail + 1L; next
    }

    covar_str <- if (length(covars) > 0) paste("+", paste(covars, collapse = " + ")) else ""
    f_med <- as.formula(paste("idp_z ~ ext_PRS_z", covar_str))
    f_out <- as.formula(paste("event ~ ext_PRS_z + idp_z", covar_str))

    med.fit <- try(lm(f_med, data = dt2), silent = TRUE)
    if (inherits(med.fit, "try-error")) { diag$lm_fail <- diag$lm_fail + 1L; next }

    out.fit <- try(glm(f_out, data = dt2, family = binomial), silent = TRUE)
    if (inherits(out.fit, "try-error")) { diag$glm_fail <- diag$glm_fail + 1L; next }

    med_res <- safe_mediate(med.fit, out.fit, "ext_PRS_z", "idp_z", SIMS)
    if (is.null(med_res)) { diag$mediate_fail <- diag$mediate_fail + 1L; next }

    row <- data.table(
      outcome = oc, idp = idp_name,
      n = n, n_case = n_case, n_ctrl = n_ctrl,
      acme = med_res$d0, acme_p = med_res$d0.p,
      ade  = med_res$z0, ade_p  = med_res$z0.p,
      prop_med = med_res$n0, prop_med_p = med_res$n0.p,
      se_mode = SE_MODE, ses_mode = SES_MODE, sims = SIMS
    )

    if (first_write) {
      fwrite(row, RAW_FILE, sep = "\t")
      first_write <- FALSE
    } else {
      fwrite(row, RAW_FILE, sep = "\t", append = TRUE, col.names = FALSE)
    }
    diag$written_ok <- diag$written_ok + 1L
  }
}

# If raw does not exist for some reason, create header-only file
if (!file.exists(RAW_FILE)) {
  empty_raw <- data.table(
    outcome=character(), idp=character(), n=integer(), n_case=integer(), n_ctrl=integer(),
    acme=double(), acme_p=double(), ade=double(), ade_p=double(),
    prop_med=double(), prop_med_p=double(),
    se_mode=character(), ses_mode=character(), sims=integer()
  )
  fwrite(empty_raw, RAW_FILE, sep = "\t")
}

# count rows in raw
raw_rows <- 0L
tmp <- tryCatch(fread(RAW_FILE), error = function(e) NULL)
if (!is.null(tmp)) raw_rows <- nrow(tmp)

# -------------------------
# write summary + done marker
# -------------------------
fwrite(
  data.table(
    se_mode = SE_MODE, ses_mode = SES_MODE, sims = SIMS,
    n_total_pairs = nrow(cands),
    n_done_pairs_preexisting = length(done_keys),
    n_todo_pairs = nrow(todo),
    n_successfully_written_this_run = diag$written_ok,
    n_rows_in_raw_after_run = raw_rows,
    # diagnostics:
    missing_dat_file = diag$missing_dat_file,
    header_read_fail = diag$header_read_fail,
    missing_event_or_prs = diag$missing_event_or_prs,
    idp_missing_in_data = diag$idp_missing_in_data,
    fread_fail = diag$fread_fail,
    qc_fail = diag$qc_fail,
    lm_fail = diag$lm_fail,
    glm_fail = diag$glm_fail,
    mediate_fail = diag$mediate_fail
  ),
  SUMMARY_FILE, sep = "\t"
)

writeLines("OK", DONE_FILE)

cat("[ok] stage4 raw written:", RAW_FILE, "\n")
cat("[ok] stage4 summary written:", SUMMARY_FILE, "\n")
cat("[ok] stage4 done marker:", DONE_FILE, "\n")

