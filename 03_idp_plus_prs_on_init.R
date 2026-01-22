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
HITS_FILE <- get_arg("--HITS_FILE", NA_character_)
SE_MODE   <- tolower(get_arg("--SE_MODE", "cluster"))
SES_MODE  <- tolower(get_arg("--SES_MODE", "off"))
if (is.na(OUTDIR) || OUTDIR=="") stop("[fatal] OUTDIR missing")
if (is.na(HITS_FILE) || !file.exists(HITS_FILE)) stop("[fatal] HITS_FILE missing: ", HITS_FILE)

DATA_DIR <- file.path(OUTDIR, "data")
dir.create(file.path(OUTDIR, "stage3"), recursive=TRUE, showWarnings=FALSE)

ses_covars <- c("income_num_z","parent_edu_years_z","parent_married_or_cohab_z",
                "income_num_z_missing","parent_edu_years_z_missing","parent_married_or_cohab_z_missing")

sig_tbl <- fread(HITS_FILE)
if (!"idp" %in% names(sig_tbl)) setnames(sig_tbl, 1, "idp")
sig_idps <- unique(sig_tbl$idp)

run_idp_models <- function(outcome_name) {
  dat_file <- file.path(DATA_DIR, paste0("extPRS_brain_", outcome_name, ".tsv"))
  cat("[diag] reading:", dat_file, "\n")
  dat <- fread(dat_file)

  covars <- c("age","sex", paste0("PC",1:20),
              "site","scanner_manufacturer","scanner_device_serial")
  if (SES_MODE=="on") covars <- c(covars, ses_covars)
  covars <- covars[covars %in% names(dat)]

  results <- list()

  for (idp in sig_idps) {
    if (!idp %in% names(dat)) next

    dat2 <- dat %>%
      filter(!is.na(.data[[idp]]), !is.na(ext_PRS), !is.na(time), !is.na(event)) %>%
      mutate(
        ext_PRS_z = as.numeric(scale(ext_PRS)[,1]),
        idp_z     = as.numeric(scale(.data[[idp]])[,1])
      )

    if (nrow(dat2) < 200) next

    covar_str <- if (length(covars)>0) paste("+", paste(covars, collapse=" + ")) else ""
    fml <- as.formula(paste("Surv(time, event) ~ ext_PRS_z + idp_z", covar_str))

    fit <- try({
      if (SE_MODE=="cluster" && ("site" %in% names(dat2))) {
        coxph(fml, data=dat2, ties="efron", cluster=dat2$site)
      } else {
        coxph(fml, data=dat2, ties="efron")
      }
    }, silent=TRUE)

    if (inherits(fit,"try-error")) next

    td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE)
    for (term in c("ext_PRS_z","idp_z")) {
      if (!(term %in% td$term)) next
      row <- td[td$term==term, ]
      results[[length(results)+1]] <- data.frame(
        outcome=outcome_name, idp=idp, term=term,
        hr=row$estimate, hr_lwr=row$conf.low, hr_upr=row$conf.high,
        p=row$p.value, n=fit$n,
        se_mode=SE_MODE, ses_mode=SES_MODE
      )
    }
  }

  if (length(results)==0) return(NULL)
  bind_rows(results)
}

res_all <- bind_rows(lapply(c("alcohol","nicotine","cannabis","any"), run_idp_models))
if (is.null(res_all) || nrow(res_all)==0) stop("[fatal] no Stage 3 models ran")

res_all <- res_all %>%
  group_by(outcome, term) %>%
  mutate(p_BH = p.adjust(p, method="BH")) %>%
  ungroup()

out_all <- file.path(OUTDIR, "stage3", "extPRS_IDPs_on_initiation_all.tsv")
fwrite(as.data.table(res_all), out_all, sep="\t")
cat("[ok] wrote", out_all, "\n")

cand <- res_all %>% filter(term=="idp_z", p_BH < 0.05)
out_cand <- file.path(OUTDIR, "stage3", "extPRS_IDPs_on_initiation_signif.tsv")
fwrite(as.data.table(cand), out_cand, sep="\t")
cat("[ok] wrote candidates:", out_cand, "\n")

