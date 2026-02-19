# Age-Adjusted Cardiovascular Risk and Hippocampal-Cognitive Coupling Along the Alzheimer's Disease Continuum

**Fernandez-Lozano S, Villeneuve S, Collins LD.**
Submitted to *Alzheimer's & Dementia*.

## Overview

This repository contains the analysis pipeline and manuscript source for a
study examining how cardiovascular risk (CVR) mediates hippocampal
volume-cognition coupling in preclinical Alzheimer's disease. Two CVR
operationalizations are compared:

1. **Framingham Risk Score (FRS)** — a widely used clinical composite
2. **CVR-MIMIC** — a latent factor from a Multiple Indicators, Multiple
   Causes model, age-adjusted at the latent level

Longitudinal mixed-effects models (LME) and latent growth curve models
(LGCM) are used to test cross-sectional associations and mediation
pathways, respectively.

## Data Access

This study uses data from the Alzheimer's Disease Neuroimaging Initiative
(ADNI) and the ADSP Phenotype Harmonization Consortium (ADSP-PHC).

- **ADNI**: Apply at <https://adni.loni.usc.edu/>
- **ADSP-PHC**: Apply at <https://www.niagads.org/adsp/content/phc>

Due to data use agreements, raw data cannot be shared. The `data/`
directory is gitignored; users must place their own data files there after
obtaining access.

## Repository Structure

```
├── R/
│   ├── scripts/          # Numbered pipeline scripts (01–15b)
│   └── utils/            # Shared utility functions
├── config/
│   └── pipeline_config.yaml
├── reports-src/
│   ├── manuscript_ad_submission.qmd
│   ├── references.bib
│   ├── title.tex
│   └── _quarto.yml
├── renv/                 # renv bootstrap
├── renv.lock             # Package versions
├── environment.yml       # Conda environment specification
├── data/                 # (gitignored) Raw and processed data
├── models/               # (gitignored) Model fits and results
├── outputs/              # (gitignored) Tables, figures, reports
└── logs/                 # (gitignored) Pipeline logs
```

## Environment Setup

### Option A: Conda + renv (recommended)

```bash
conda env create -f environment.yml
conda activate ukb_adni_lgcm
# R packages are managed by renv; on first launch:
Rscript -e 'renv::restore()'
```

### Option B: renv only

Requires R >= 4.5 and system dependencies for compiled packages.

```bash
Rscript -e 'renv::restore()'
```

## Pipeline

Scripts are numbered by execution order. Run them sequentially:

```
01  preprocess_volumes.R           Raw volume preprocessing
02  compute_zscores.R              GAMLSS-based z-scores
03  prepare_analysis_data.R        Merge clinical + imaging
04  frs_problem_analysis.R         FRS ceiling/age-weight diagnostics
05  cvr_mimic_model.R              MIMIC CFA model fitting
06  cvr_measurement_invariance.R   Multi-group MI across sex
07  cvr_score_extraction.R         Age-adjusted factor scores
08  lme_hvr_z.R                    LME: HVR-z ~ CVR
09  lme_hc_z.R                     LME: hippocampal cognition
10  lme_sensitivity.R              Sensitivity analyses
11  prepare_lgcm_data.R            Reshape for LGCM
12  lgcm_parallel.R                Parallel-process LGCM
13  lgcm_mediation.R               Mediation bootstrap
14  lgcm_simulation.R              Power / simulation
15b prepare_manuscript_env.R       Assemble manuscript data
```

Dependency chain: `01 → 02 → 03 → 04 → 05 → 06 → 07 → {08, 09, 10} → 11 → {12, 13} → 14 → 15b`

## Rendering the Manuscript

After running the full pipeline:

```bash
cd reports-src
quarto render manuscript_ad_submission.qmd --to pdf
```

Output is written to `outputs/reports/`.

## Configuration

All paths, seeds, and analysis parameters are centralized in
`config/pipeline_config.yaml`. Scripts read paths via `get_data_path()`
and settings via `get_script_setting()` — no hardcoded paths.

## Citation

If you use this code, please cite:

> Fernandez-Lozano S, Villeneuve S, Collins LD. Age-Adjusted
> Cardiovascular Risk and Hippocampal-Cognitive Coupling Along the
> Alzheimer's Disease Continuum. *Alzheimer's & Dementia* (submitted).

## License

Apache License 2.0. See [LICENSE](LICENSE).
