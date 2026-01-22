#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(broom)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if (length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

OUTDIR    <- get_arg("--OUTDIR", NA_character_)
SE_MODE   <- tolower(get_arg("--SE_MODE", "cluster"))  # cluster | hc3 (means no cluster)
SES_MODE  <- tolower(get_arg("--SES_MODE", "off"))     # on | off
if (is.na(OUTDIR) || OUTDIR=="") stop("[fatal] OUTDIR missing")

DATA_DIR <- file.path(OUTDIR, "data")
dir.create(file.path(OUTDIR, "stage2"), recursive=TRUE, showWarnings=FALSE)

ses_covars <- c("income_num_z","parent_edu_years_z","parent_married_or_cohab_z",
                "income_num_z_missing","parent_edu_years_z_missing","parent_married_or_cohab_z_missing")

do_direct_prs <- function(outcome_name) {
  dat_file <- file.path(DATA_DIR, paste0("extPRS_brain_", outcome_name, ".tsv"))
  cat("[diag] reading:", dat_file, "\n")
  dat <- fread(dat_file)

  covars <- c("age","sex", paste0("PC",1:20),
              # keep raw site/scanner fields if present (not one-hot here)
              "site","scanner_manufacturer","scanner_device_serial")

  # SES on/off
  if (SES_MODE == "on") covars <- c(covars, ses_covars)
  covars <- covars[covars %in% names(dat)]

  fml <- as.formula(paste("Surv(time, event) ~ ext_PRS +", paste(covars, collapse=" + ")))

  fit <- if (SE_MODE=="cluster" && ("site" %in% names(dat))) {
    coxph(fml, data=dat, ties="efron", cluster=dat$site)
  } else {
    coxph(fml, data=dat, ties="efron")
  }

  td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE)
  out <- td %>%
    filter(term=="ext_PRS") %>%
    mutate(outcome=outcome_name, se_mode=SE_MODE, ses_mode=SES_MODE, n=fit$n) %>%
    select(outcome, se_mode, ses_mode, n, estimate, conf.low, conf.high, p.value)

  out
}

res_direct <- bind_rows(lapply(c("alcohol","nicotine","cannabis","any"), do_direct_prs))

out_file <- file.path(OUTDIR, "stage2", "extPRS_direct_on_initiation.tsv")
fwrite(as.data.table(res_direct), out_file, sep="\t")
cat("[ok] wrote", out_file, "\n")

