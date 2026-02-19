#!/usr/bin/env Rscript

# ============================================================
# 07_cvr_score_extraction.R
# ============================================================
# PURPOSE:
#   Extract age-adjusted CVR latent factor scores from
#   the MIMIC model (script 05) for use in LME models.
#   Merge with analysis cohort and standardize both FRS
#   and CVR_mimic on the common sample.
#
# AGE-ADJUSTMENT METHODOLOGY:
#   The MIMIC model is: CVR ~ Age, with 5 indicators.
#   lavPredict() returns eta (the full latent variable),
#   which includes Age-predicted variance (r~0.49 with
#   Age). For LME comparisons with FRS, we need CVR
#   *independent* of Age.
#
#   We use latent-level residualization:
#     zeta = eta - gamma * AGE_c
#   where:
#     eta   = lavPredict(mimic_fit) factor scores
#     gamma = MIMIC regression coefficient (CVR ~ AGE_c)
#     AGE_c = centered age (same as model input)
#
#   This targets the disturbance term of the MIMIC
#   structural equation, which is orthogonal to Age by
#   construction. It is the pragmatic approximation of
#   Residual SEM (Asparouhov & Muthen, 2023).
#
#   Perplexity second-opinion (2026-02-10) confirmed:
#     - Latent-level residualization is defensible
#     - Equivalent to extracting the disturbance term
#     - Changes the construct from "CVR" to
#       "CVR beyond age expectations"
#     - Appropriate when the goal is to compare with
#       another age-confounded measure (FRS)
#
# LGCM NOTE:
#   For LGCM (scripts 12-13), the MIMIC measurement
#   model is embedded directly in the OpenMx structural
#   model. This is the gold-standard single-step SEM
#   approach (Muthen, 2006). No score extraction is
#   needed for that pathway. The zeta scores computed
#   here are ONLY for LME comparisons.
#
# STANDARDIZATION:
#   Both FRS and CVR_mimic (zeta) are z-scored on the
#   SAME common sample (subjects with both measures).
#   This ensures beta coefficients are directly
#   comparable across predictors.
#
# INPUTS:
#   - models/results/lme/cvr_mimic_specification.rds
#   - data/raw/adsp_phc/ADSP_PHC_CVRF_*.csv
#   - data/raw/All_Subjects_LABDATA_*.csv
#   - data/derivatives/lme/cohort_hvr.rds
#
# OUTPUTS:
#   - data/derivatives/lme/cvr_mimic_scores.rds
#   - data/derivatives/lme/cohort_hvr.rds (updated)
#
# ============================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(lavaan)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))

log_script_start("07_cvr_score_extraction.R")
config <- load_config()
validate_config(config)


# ============================================================
# PART 1: PREPARE DATA
# ============================================================
log_section("Part 1: Prepare Data")

# Load MIMIC specification from script 05
spec_path <- get_data_path(
  "models", "cvr_mimic_specification"
)
check_files_exist(spec_path)
spec.lst <- read_rds_safe(
  spec_path, "MIMIC specification"
)
log_info("  Syntax: %s", trimws(spec.lst$lavaan_syntax))

# Load and merge raw indicator data
cvrf_path <- get_data_path("raw", "adsp_phc_cvrf")
cvrf.dt <- fread(cvrf_path)
cvrf.dt <- cvrf.dt[order(RID, EXAMDATE)]
cvrf.dt <- cvrf.dt[, .SD[1], by = RID]

lab_path <- get_data_path("raw", "adni_labdata")
lab.dt <- fread(lab_path)
for (v in c("RCT11", "RCT20", "RCT392")) {
  lab.dt[, (v) := suppressWarnings(
    as.numeric(get(v))
  )]
}
lab.dt <- lab.dt[order(RID, EXAMDATE)]
lab_bl.dt <- lab.dt[, .(
  GLUCOSE = RCT11[1],
  CHOLESTEROL = RCT20[1],
  CREATININE = RCT392[1]
), by = RID]

analysis.dt <- merge(
  cvrf.dt, lab_bl.dt, by = "RID", all.x = TRUE
)
log_info("Merged data: %d subjects", nrow(analysis.dt))

# Prepare model variables
mod.dt <- data.table(
  RID = analysis.dt$RID,
  SEX = ifelse(
    analysis.dt$PHC_Sex == 1, "Male", "Female"
  ),
  AGE_raw = analysis.dt$PHC_Age_CardiovascularRisk,
  FRS = analysis.dt$PHC_ASCVD_10y_FRS_Simple_Ageover30,
  SBP_z = as.numeric(scale(analysis.dt$PHC_SBP)),
  HTN = as.numeric(analysis.dt$PHC_Hypertension),
  GLUCOSE_z = as.numeric(
    scale(analysis.dt$GLUCOSE)
  ),
  CHOL_z = as.numeric(
    scale(analysis.dt$CHOLESTEROL)
  ),
  CREAT_z = as.numeric(
    scale(analysis.dt$CREATININE)
  ),
  AGE_c = as.numeric(scale(
    analysis.dt$PHC_Age_CardiovascularRisk
  ))
)
log_info(
  "Model data: %d subjects, %d with Age",
  nrow(mod.dt), sum(!is.na(mod.dt$AGE_c))
)

# ============================================================
# PART 2: FIT MIMIC & EXTRACT SCORES
# ============================================================
log_section("Part 2: Fit MIMIC & Extract Scores")

# Fit the final MIMIC model
mimic.fit <- sem(
  spec.lst$lavaan_syntax,
  data = mod.dt,
  estimator = "MLR",
  missing = "fiml"
)
log_info(
  "MIMIC converged: %s",
  lavInspect(mimic.fit, "converged")
)

# Report fit
fi.v <- fitMeasures(mimic.fit, c(
  "cfi.robust", "tli.robust",
  "rmsea.robust", "srmr"
))
log_info(
  "  CFI=%.3f TLI=%.3f RMSEA=%.3f SRMR=%.3f",
  fi.v[1], fi.v[2], fi.v[3], fi.v[4]
)

# Extract eta (full latent scores)
eta.v <- as.numeric(lavPredict(mimic.fit))
log_info(
  "  eta: N=%d, M=%.4f, SD=%.4f",
  sum(!is.na(eta.v)),
  mean(eta.v, na.rm = TRUE),
  sd(eta.v, na.rm = TRUE)
)

# Get gamma (Age regression coefficient)
pe.dt <- as.data.table(
  parameterEstimates(mimic.fit, standardized = TRUE)
)
gamma_row <- pe.dt[op == "~" & rhs == "AGE_c"]
gamma_est <- gamma_row$est[1]
gamma_std <- gamma_row$std.all[1]
log_info(
  "  gamma (CVR ~ Age): est=%.4f, std=%.4f",
  gamma_est, gamma_std
)

# Compute zeta = eta - gamma * AGE_c
# This is the latent disturbance: CVR independent
# of Age. Equivalent to the RSEM disturbance term.
zeta.v <- eta.v - gamma_est * mod.dt$AGE_c
mod.dt[, CVR_mimic := zeta.v]

log_info(
  "  zeta: M=%.4f, SD=%.4f",
  mean(zeta.v, na.rm = TRUE),
  sd(zeta.v, na.rm = TRUE)
)

# ============================================================
# PART 3: VALIDATION
# ============================================================
log_section("Part 3: Validation")

# --- Age independence ---
r_age_eta <- cor(
  eta.v, mod.dt$AGE_c,
  use = "complete.obs"
)
r_age_zeta <- cor(
  zeta.v, mod.dt$AGE_c,
  use = "complete.obs"
)
log_info("Age correlations:")
log_info(
  "  r(eta, Age) = %.3f (before adjustment)",
  r_age_eta
)
log_info(
  "  r(zeta, Age) = %.3f (after adjustment)",
  r_age_zeta
)
if (abs(r_age_zeta) > 0.05) {
  log_warn(
    "  zeta-Age r > 0.05; check methodology"
  )
} else {
  log_info("  PASS: zeta orthogonal to Age")
}

# --- Convergent validity with FRS ---
has_frs <- !is.na(mod.dt$FRS)
r_frs_eta <- cor(
  eta.v[has_frs], mod.dt$FRS[has_frs],
  use = "complete.obs"
)
r_frs_zeta <- cor(
  zeta.v[has_frs], mod.dt$FRS[has_frs],
  use = "complete.obs"
)
log_info("")
log_info("FRS convergent validity:")
log_info(
  "  r(eta, FRS) = %.3f", r_frs_eta
)
log_info(
  "  r(zeta, FRS) = %.3f", r_frs_zeta
)

# By sex
for (sx in c("Male", "Female")) {
  idx <- mod.dt$SEX == sx & has_frs
  r_sx <- cor(
    zeta.v[idx], mod.dt$FRS[idx],
    use = "complete.obs"
  )
  log_info(
    "  r(zeta, FRS) [%s] = %.3f", sx, r_sx
  )
}

# --- FRS-Age for comparison ---
r_frs_age <- cor(
  mod.dt$FRS, mod.dt$AGE_c,
  use = "complete.obs"
)
log_info("")
log_info("Comparison: r(FRS, Age) = %.3f", r_frs_age)
log_info(
  "  FRS confounded; zeta age-orthogonal"
)

# ============================================================
# PART 4: SAVE SCORES
# ============================================================
log_section("Part 4: Save Scores")

scores.dt <- mod.dt[, .(
  RID, SEX, AGE_raw, FRS,
  CVR_mimic, AGE_c,
  SBP_z, HTN, GLUCOSE_z, CHOL_z, CREAT_z
)]
log_info(
  "Score table: %d subjects", nrow(scores.dt)
)

scores_path <- get_data_path(
  "models", "cvr_mimic_scores"
)
write_rds_safe(
  scores.dt, scores_path,
  "CVR MIMIC scores"
)

# ============================================================
# PART 5: MERGE WITH COHORT & STANDARDIZE
# ============================================================
log_section("Part 5: Merge & Standardize")

# Load analysis cohort from script 03
cohort_path <- get_data_path(
  "derivatives", "lme_cohort_hvr"
)
check_files_exist(cohort_path)
cohort.dt <- read_rds_safe(
  cohort_path, "LME analysis cohort"
)
log_info(
  "Loaded cohort: %d obs, %d subjects",
  nrow(cohort.dt),
  cohort.dt[, uniqueN(RID)]
)

# Merge CVR_mimic (one score per subject)
cvr_merge.dt <- scores.dt[
  !is.na(CVR_mimic),
  .(RID, CVR_mimic)
]
# Remove any existing CVR_mimic column
if ("CVR_mimic" %in% names(cohort.dt)) {
  cohort.dt[, CVR_mimic := NULL]
}
cohort.dt <- merge(
  cohort.dt, cvr_merge.dt,
  by = "RID", all.x = TRUE
)

n_total <- cohort.dt[VISIT == 1, uniqueN(RID)]
n_cvr <- cohort.dt[
  VISIT == 1 & !is.na(CVR_mimic),
  uniqueN(RID)
]
log_info(
  "CVR_mimic merged: %d/%d subjects",
  n_cvr, n_total
)

# Filter to common sample (both FRS and CVR_mimic)
cohort.dt <- cohort.dt[!is.na(CVR_mimic)]
n_common <- cohort.dt[, uniqueN(RID)]
log_info(
  "Common sample (FRS + CVR_mimic): %d subjects",
  n_common
)

# --- Z-score BOTH on common sample ---
# This ensures beta comparability
log_info("")
log_info("Standardizing on common sample:")

# FRS
frs_mean <- mean(cohort.dt$FRS, na.rm = TRUE)
frs_sd <- sd(cohort.dt$FRS, na.rm = TRUE)
cohort.dt[, FRS_raw := FRS]
cohort.dt[, FRS := (FRS - frs_mean) / frs_sd]
log_info(
  "  FRS: M=%.2f SD=%.2f -> z-scored",
  frs_mean, frs_sd
)

# CVR_mimic (zeta)
cvr_mean <- mean(
  cohort.dt$CVR_mimic, na.rm = TRUE
)
cvr_sd <- sd(
  cohort.dt$CVR_mimic, na.rm = TRUE
)
cohort.dt[
  , CVR_mimic := (CVR_mimic - cvr_mean) / cvr_sd
]
log_info(
  "  CVR_mimic: M=%.4f SD=%.4f -> z-scored",
  cvr_mean, cvr_sd
)

# Verify
log_info("")
log_info("Post-standardization:")
log_info(
  "  FRS: M=%.4f SD=%.4f",
  mean(cohort.dt$FRS, na.rm = TRUE),
  sd(cohort.dt$FRS, na.rm = TRUE)
)
log_info(
  "  CVR_mimic: M=%.4f SD=%.4f",
  mean(cohort.dt$CVR_mimic, na.rm = TRUE),
  sd(cohort.dt$CVR_mimic, na.rm = TRUE)
)

# Store standardization metadata
attr(cohort.dt, "cvr_standardization") <- list(
  method = "global_zscore_common_sample",
  frs_raw_mean = frs_mean,
  frs_raw_sd = frs_sd,
  cvr_mimic_raw_mean = cvr_mean,
  cvr_mimic_raw_sd = cvr_sd,
  n_common_sample = n_common,
  age_adjustment = paste(
    "Latent-level residualization:",
    "zeta = eta - gamma * AGE_c,",
    "where eta = lavPredict(MIMIC),",
    "gamma = MIMIC regression coef.",
    "Targets the disturbance term,",
    "orthogonal to Age."
  ),
  timestamp = Sys.time()
)

# Final summary
bl.dt <- cohort.dt[VISIT == 1]
log_info("")
log_info("FINAL ANALYSIS COHORT:")
log_info(
  "  N = %d subjects, %d observations",
  nrow(bl.dt), nrow(cohort.dt)
)
log_info(
  "  Sex: %d Males, %d Females",
  sum(bl.dt$SEX == "Male"),
  sum(bl.dt$SEX == "Female")
)
r_final <- cor(
  bl.dt$FRS, bl.dt$CVR_mimic,
  use = "complete.obs"
)
log_info(
  "  r(FRS_z, CVR_mimic_z) = %.3f", r_final
)

# Save updated cohort
write_rds_safe(
  cohort.dt, cohort_path,
  "LME cohort (CVR-standardized)"
)

log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info("07_cvr_score_extraction.R COMPLETE")
log_info(paste(rep("=", 70), collapse = ""))
