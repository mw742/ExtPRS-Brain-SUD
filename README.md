# ExtPRS-Brain-SUD


This repository contains code for a **four-stage analysis pipeline** linking **externalizing polygenic risk scores (extPRS)**, **baseline neuroimaging-derived phenotypes (IDPs)**, and **time-to-initiation of substance use** in the ABCD Study®.

The pipeline integrates large-scale IDP association screening, survival analysis, and mediation analysis to identify neurobiological pathways linking genetic liability to adolescent substance use initiation.


---

## Overview

A **four-stage analysis pipeline**:

1. **Stage 1 – extPRS → IDP screening**
   00extPRS_on_IDPs.slurm
   
   IDP-by-IDP linear regression to identify neuroimaging phenotypes associated with extPRS.
   

3. **Stage 2 – extPRS → initiation (direct effect)**
   01_prep_extPRS_brain_init_data.R
   02_extPRS_to_init.R
   
   Cox proportional hazards models for time-to-initiation of alcohol, nicotine, cannabis, and any substance.

5. **Stage 3 – extPRS + IDP → initiation (joint models)**
   03_idp_plus_prs_on_init.R
     
   Joint Cox models to prioritize IDPs associated with initiation beyond extPRS.

7. **Stage 4 – Mediation analysis**
   04_mediation_extPRS_brain_init_bootstrap.R
     
   Bootstrap-based mediation analysis to estimate indirect (ACME) and direct (ADE) effects of extPRS via IDPs.


---

## Data availability

This repository does **not** include individual-level ABCD data, which are controlled-access.

To run the pipeline, users must obtain authorized access to the ABCD Study® and prepare related data. 



---

## Citation

If you use this work, please cite:

Wei, M., & Peng, Q. (2026). *Externalizing Polygenic Liability, Brain Imaging Phenotypes, and Adolescent Substance Use Initiations: A Multistage Association and Mediation Analysis in ABCD*. bioRxiv. https://doi.org/10.64898/2026.04.08

---

## Data access and responsible use

This repository distributes analysis code and does not include participant-level ABCD Study data. Users must obtain input data through authorized providers and comply with applicable data-use agreements and institutional requirements. The software license does not grant permission to redistribute external datasets or restricted derived outputs.

Generated files may contain participant identifiers, polygenic scores, imaging measurements, covariates, and substance-use outcomes. In particular, the prepared `data/extPRS_brain_*.tsv` files contain participant-level information. Do not upload these files, or other outputs containing restricted information, to a public repository without appropriate authorization.

For reproducibility, document the cohort and inclusion criteria, data release/version, accession or DOI where available, PRS source and construction method, phenotype and follow-up definitions, preprocessing steps, software versions, analysis settings, and exact code commit used. Include the acknowledgments required by the relevant data providers.

## Validation status and updates

This repository contains research code that may include undetected errors or limitations. Review of commit `db13a34` included source-code inspection, Bash and embedded Python syntax checks, and limited synthetic checks of covariate handling. It did not include execution of the R analyses, end-to-end pipeline testing, or independent reproduction of the associated study results.

Issues identified during review concern covariate inclusion, PRS-column selection, reporting of inference settings, and reuse of existing results across runs. These issues require correction and validation. Their effects on the associated study results have not been established.

The mediation stage currently models binary substance-use initiation status and does not incorporate follow-up time or censoring. Its estimates should therefore be distinguished from the time-to-event associations estimated in the preceding Cox models. Observational mediation estimates do not, by themselves, establish causal mechanisms.

Users should independently assess whether the implementation and its assumptions are suitable for their intended analyses. Confirmed issues, corrections, and subsequent validation will be documented in the commit history or release notes, including affected versions and known effects on outputs.


