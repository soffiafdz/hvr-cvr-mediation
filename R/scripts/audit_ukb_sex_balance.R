#!/usr/bin/env Rscript
# =============================================================================
# audit_ukb_sex_balance.R
# =============================================================================
# Reproducible audit of the UKB normative-reference sex imbalance.
#
# Background:
#   The GAMLSS normative reference (data/external/gamlss.rds) used to
#   z-score ADNI HVR/HC values is 56.2% female (15,412 F vs 12,022 M).
#   This script recomputes the per-sex counts at each upstream pipeline
#   stage so the imbalance can be defended in a reviewer response.
#
# Stages audited (matches upstream 02_quality_control.R + 05_adjust_headsize.R):
#   1. Raw UKB imaging cohort (any ses-2 or ses-3 MRI date)
#   2. DARQ QC failures (DARQ < 0.25)
#   3. AssemblyNet QC failures (grade != "A")
#   4. ICD-10 exclusions: F (psychiatric), G (neurological), Q0 (CNS malf.)
#   5. Final GAMLSS sample (all splits + train-only split)
#
# Caveats (not reproduced here):
#   - Outlier removal (|Z| > 3 on volume residuals): requires refitting the
#     upstream lm(VAL ~ SEX + AGE + ICC) and is sex-adjusted by design.
#   - Manually-curated bad-EID exclusion list used in 05_adjust_headsize.R
#     (loaded via exclude_ids_file) is small and not included.
#
# Inputs (upstream UKB data; absolute paths):
#   - /ipl/ipl27/sfernandez/data_phd/ukb/fst/ukb_covars.fst
#   - /ipl/ipl27/sfernandez/data_phd/ukb/derivatives/qc_darq-assemblynet.rds
#   - data/external/gamlss.rds (this project)
#
# Output:
#   - outputs/audit/ukb_sex_balance.rds (per-stage counts)
#   - Console: human-readable summary table
#
# Run (this project's renv lacks `fst`; the upstream UKB project has it):
#   cd /data/ipl/ipl27/sfernandez/headsize_sexeffects_ukb
#   Rscript /ipl/ipl27/sfernandez/hvr_cvr_mediation/R/scripts/\
#     audit_ukb_sex_balance.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("fst", quietly = TRUE)) {
  stop(
    "Package `fst` not available in current R library.\n",
    "Run this script from the upstream UKB project root:\n",
    "  cd /data/ipl/ipl27/sfernandez/headsize_sexeffects_ukb\n",
    "  Rscript ",
    "/ipl/ipl27/sfernandez/hvr_cvr_mediation/R/scripts/",
    "audit_ukb_sex_balance.R",
    call. = FALSE
  )
}
library(fst)

# --- Paths (absolute; portable across project roots) ---
ukb_covars.path <- file.path(
  "/ipl/ipl27/sfernandez/data_phd/ukb",
  "fst/ukb_covars.fst"
)
ukb_qc.path <- file.path(
  "/ipl/ipl27/sfernandez/data_phd/ukb",
  "derivatives/qc_darq-assemblynet.rds"
)
gamlss.path <- file.path(
  "/ipl/ipl27/sfernandez/hvr_cvr_mediation",
  "data/external/gamlss.rds"
)

out.dir <- file.path(
  "/ipl/ipl27/sfernandez/hvr_cvr_mediation",
  "outputs/audit"
)
out.path <- file.path(out.dir, "ukb_sex_balance.rds")

# --- Thresholds (mirroring upstream pipeline_config.yaml) ---
DARQ_THRESHOLD <- 0.25
ASNET_PASS <- "A"
ICD_PRIMARY <- "^F|^G|^Q0"

# --- Load ---
cat("Loading UKB covariates ...\n")
covars.dt <- read_fst(ukb_covars.path, as.data.table = TRUE)

cat("Loading QC scores ...\n")
qc.dt <- readRDS(ukb_qc.path)

cat("Loading GAMLSS reference sample ...\n")
gamlss.lst <- readRDS(gamlss.path)
final.dt <- as.data.table(gamlss.lst$DATA$CRS)
rm(gamlss.lst)

# --- Stage 1: raw imaging cohort ---
imaging.dt <- covars.dt[
  !is.na(DATE_mri_ses2) | !is.na(DATE_mri_ses3),
  .(EID, SEX = SEX_r)
]
stage_imaging.dt <- imaging.dt[, .N, keyby = SEX]

# --- Stage 2: QC outcomes (joined with sex) ---
sex.dt <- covars.dt[, .(EID, SEX = SEX_r)]
qc_sex.dt <- merge(qc.dt, sex.dt, by = "EID")

darq_fail.dt <- qc_sex.dt[
  DARQ < DARQ_THRESHOLD, .N, keyby = SEX
]
asnet_fail.dt <- qc_sex.dt[
  ASBLYNET != ASNET_PASS, .N, keyby = SEX
]

# Failure rates (% within sex)
qc_rates.dt <- qc_sex.dt[, .(
  total = .N,
  darq_fail_pct = round(
    100 * mean(DARQ < DARQ_THRESHOLD, na.rm = TRUE), 2
  ),
  asnet_fail_pct = round(
    100 * mean(ASBLYNET != ASNET_PASS, na.rm = TRUE), 2
  )
), keyby = SEX]

# --- Stage 3: ICD-10 exclusions within imaging cohort ---
imaging_eids.v <- imaging.dt$EID
icd_dt.fn <- function(pattern.s) {
  covars.dt[
    EID %in% imaging_eids.v & ICD_10 %like% pattern.s,
    .N, keyby = .(SEX = SEX_r)
  ]
}
icd_f.dt <- icd_dt.fn("^F")
icd_g.dt <- icd_dt.fn("^G")
icd_q0.dt <- icd_dt.fn("^Q0")
icd_pri.dt <- icd_dt.fn(ICD_PRIMARY)

# --- Stage 4: final GAMLSS sample ---
final_all.dt <- final.dt[!duplicated(EID), .N, keyby = SEX]
final_train.dt <- final.dt[
  SPLIT == "train", .SD[1], by = EID
][, .N, keyby = SEX]

# --- Helper: extract per-sex N, defaulting to 0 ---
get_n.fn <- function(dt, sex.s) {
  if (nrow(dt) == 0) return(0L)
  v <- dt[SEX == sex.s, N]
  if (length(v) == 0) 0L else v
}

# --- Assemble summary ---
stages.v <- c(
  "1. Raw imaging cohort (ses-2 or ses-3 MRI)",
  "2. DARQ QC failures (DARQ < 0.25)",
  "3. AssemblyNet QC failures (grade != A)",
  "4. ICD-10 F-code (any)",
  "5. ICD-10 G-code (any)",
  "6. ICD-10 Q0-code (any)",
  "7. ICD-10 primary exclusion (F|G|Q0)",
  "8. Final GAMLSS sample (all splits)",
  "9. Final GAMLSS sample (train only)"
)
stage_dts.lst <- list(
  stage_imaging.dt,
  darq_fail.dt,
  asnet_fail.dt,
  icd_f.dt,
  icd_g.dt,
  icd_q0.dt,
  icd_pri.dt,
  final_all.dt,
  final_train.dt
)

summary.dt <- data.table(
  stage = stages.v,
  female = vapply(
    stage_dts.lst, get_n.fn, integer(1), sex.s = "Female"
  ),
  male = vapply(
    stage_dts.lst, get_n.fn, integer(1), sex.s = "Male"
  )
)
summary.dt[, total := female + male]
summary.dt[, pct_female := round(
  100 * female / pmax(total, 1L), 1
)]

cat("\n=== UKB sex-balance audit ===\n")
print(summary.dt, row.names = FALSE)

cat("\n=== QC failure rates within sex (%) ===\n")
print(qc_rates.dt, row.names = FALSE)

# --- Save ---
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(stages = summary.dt, qc_rates = qc_rates.dt),
  out.path
)
cat(sprintf("\nSaved: %s\n", out.path))
