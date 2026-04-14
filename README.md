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



