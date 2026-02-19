#!/usr/bin/env Rscript

# =============================================================================
# 08_prepare_lgcm_data.R
# =============================================================================
# Prepare wide-format data for LGCM analysis using EDT as time metric
#
# Key features:
#   1. EDT (Estimated Distance to Dementia) as time metric - NOT age-centered
#   2. Sex variable for stratified analyses
#   3. FRS × HVR interaction term pre-computed
#   4. SE² columns for fixing residual variance in lavaan/OpenMx
#
# Source: Analysis cohort from 03/04 pipeline (has EDT + CVR)
#
# Outputs:
#   - data/derivatives/lgcm_input_full.rds (wide format, full Aβ+ sample with EDT)
#   - data/derivatives/lgcm_input_male.rds (male subsample)
#   - data/derivatives/lgcm_input_female.rds (female subsample)
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

# Source utilities
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# Initialize
log_script_start("11_prepare_lgcm_data.R")
config <- load_config()
validate_config(config)  # Validate config structure
set_seed()

# --- Configuration ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "prepare_lgcm", default = FALSE
)

output_dir <- get_data_path("derivatives", "lgcm_dir")
ensure_directory(output_dir)

output_full_path <- get_data_path("derivatives", "lgcm_input_full")
output_male_path <- get_data_path("derivatives", "lgcm_input_male")
output_female_path <- get_data_path("derivatives", "lgcm_input_female")

# Check if regeneration is needed using utility function
input_path <- get_data_path("derivatives", "lme_cohort_hvr")
if (!needs_regeneration(output_full_path, input_path, FORCE_REGENERATE)) {
  log_info("Output files are up-to-date and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE to rerun.")
  log_script_end("08_prepare_lgcm_data.R", success = TRUE)
  quit(status = 0)
}

# =============================================================================
# 1. Load LME Cohort (has EDT)
# =============================================================================
log_section("Loading LME Cohort with EDT")

# The LME cohort was prepared in 02_prepare_lme_data.R
# It has the validated EDT values
lme_cohort_path <- get_data_path("derivatives", "lme_cohort_hvr")

if (!file.exists(lme_cohort_path)) {
  log_error("LME cohort not found at: %s", lme_cohort_path)
  log_error("Run 02_prepare_lme_data.R first to generate cohort with EDT")
  stop("Missing LME cohort")
}

lme.dt <- read_rds_safe(lme_cohort_path, "LME cohort with EDT")

# Load MIMIC indicators for OpenMx embedding
mimic_scores_path <- get_data_path(
  "models", "cvr_mimic_scores"
)
check_files_exist(mimic_scores_path)
mimic.dt <- read_rds_safe(
  mimic_scores_path, "CVR MIMIC scores"
)
mimic_ind_cols.v <- c(
  "SBP_z", "HTN", "GLUCOSE_z",
  "CHOL_z", "CREAT_z"
)
missing_ind.v <- setdiff(
  mimic_ind_cols.v, names(mimic.dt)
)
if (length(missing_ind.v) > 0) {
  stop(sprintf(
    "MIMIC scores missing indicators: %s",
    paste(missing_ind.v, collapse = ", ")
  ))
}
merge_cols.v <- c("RID", mimic_ind_cols.v)
mimic_merge.dt <- unique(
  mimic.dt[, ..merge_cols.v]
)
# Remove existing cols if present
for (col in mimic_ind_cols.v) {
  if (col %in% names(lme.dt))
    lme.dt[, (col) := NULL]
}
lme.dt <- merge(
  lme.dt, mimic_merge.dt,
  by = "RID", all.x = TRUE
)
log_info(
  "MIMIC indicators merged: %d/%d",
  sum(!is.na(lme.dt[[mimic_ind_cols.v[1]]])),
  nrow(lme.dt)
)

log_info("LME cohort loaded:")
log_info("  N observations: %d", nrow(lme.dt))
log_info("  N subjects: %d", uniqueN(lme.dt$PTID))

# Validate using utility functions
required_cols.v <- c(
  "PTID", "VISIT", "YRS_from_bl", "EDT", "EDT_bl_years",
  "SEX", "FRS", "AGE", "EDUC",
  "HVR_z", "HC_z", "HVR_raw", "HC_raw", "TIV",
  "MEM", "LAN", "EXF", "MEM_SE", "LAN_SE", "EXF_SE"
)
validate_columns(lme.dt, required_cols.v, "LME cohort")
validate_not_empty(lme.dt, "LME cohort")
validate_age(lme.dt, "AGE", min_age = 50, max_age = 100)
validate_categorical(
  lme.dt, "SEX", get_factor_levels("SEX"),
  allow_extra = FALSE
)
# FRS is now z-scored upstream (M=0, SD=1); skip range check
# validate_num_range would fail on negative z-scores

# Check if CVR_mimic is available (merged in script 03)
HAS_CVR_MIMIC <- "CVR_mimic" %in% names(lme.dt) &&
  sum(!is.na(lme.dt$CVR_mimic)) > 0
if (HAS_CVR_MIMIC) {
  log_info(
    "CVR_mimic available: %d non-NA observations",
    sum(!is.na(lme.dt$CVR_mimic))
  )
} else {
  stop(
    "CVR_mimic not in cohort. Run script 07 first."
  )
}

log_info("All validation checks passed")

# =============================================================================
# 2. Validate EDT
# =============================================================================
log_section("Validating EDT")

# EDT is the estimated distance to dementia in years
# (Raket 2020, progmod package)
# Negative = pre-onset (before expected AD conversion)
# Positive = post-onset (past expected AD conversion)

edt_summary.dt <- lme.dt[, .(
  n = .N,
  n_na = sum(is.na(EDT)),
  min = min(EDT, na.rm = TRUE),
  max = max(EDT, na.rm = TRUE),
  mean = mean(EDT, na.rm = TRUE),
  sd = sd(EDT, na.rm = TRUE)
)]

log_info("EDT summary:")
log_info("  N total: %d, N missing: %d",
         edt_summary.dt$n, edt_summary.dt$n_na)
log_info("  Range: [%.2f, %.2f] years",
         edt_summary.dt$min, edt_summary.dt$max)
log_info("  Mean (SD): %.2f (%.2f)",
         edt_summary.dt$mean, edt_summary.dt$sd)

# Remove observations without EDT
lme.dt <- lme.dt[!is.na(EDT)]
log_info("After removing NA EDT: %d observations, %d subjects",
         nrow(lme.dt), uniqueN(lme.dt$PTID))

# =============================================================================
# 3. Prepare Variables for LGCM
# =============================================================================
log_section("Preparing Variables")

# Standardize variable names for LGCM
# First validate that required columns exist (fail early if not)
old_names.v <- c(
  "HVR_z", "HC_z", "HVR_raw", "HC_raw", "TIV",
  "MEM", "LAN", "EXF", "MEM_SE", "LAN_SE", "EXF_SE"
)
missing_cols.v <- setdiff(old_names.v, names(lme.dt))
if (length(missing_cols.v) > 0) {
  log_error("Missing required columns for renaming: %s",
            paste(missing_cols.v, collapse = ", "))
  stop("Missing columns in LME cohort")
}

setnames(lme.dt,
         old_names.v,
         c("HVR_Z", "HC_Z", "HVR_RAW", "HC_RAW", "TIV",
           "PHC_MEM", "PHC_LAN", "PHC_EXF",
           "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE"))

# Log brain measures availability
log_info("Brain measures for LGCM: HVR_Z, HC_Z, HVR_RAW, HC_RAW (+TIV)")

# Create timepoint labels (T0, T1, ...)
setorder(lme.dt, PTID, VISIT)
lme.dt[, TIMEPOINT := paste0("T", VISIT - 1)]

# =============================================================================
# FRS Centering Strategy
# =============================================================================
# We CENTER FRS within sex (subtract sex-specific mean) rather than STANDARDIZE
# (divide by SD) for the following reasons:
#
# 1. INTERPRETABILITY: Centered coefficients represent the effect at the
#    sex-specific mean FRS. This is more clinically meaningful than "per SD".
#
# 2. SEX DIFFERENCES: Men and women have different FRS distributions
#    (different means and SDs). Sex-centering removes between-sex mean
#    differences while preserving within-sex variation.
#
# 3. INTERACTION INTERPRETATION: For the FRS × HVR interaction, centering
#    ensures the main effect of HVR is evaluated at average FRS (within sex).
#
# 4. P-VALUES UNAFFECTED: Centering/standardizing changes coefficient values
#    but NOT statistical significance (p-values identical).
#
# NOTE: For standardized betas, we use mxStandardizeRAMpaths() on the fitted
# OpenMx model, which computes proper standardized path coefficients based on
# the observed variance-covariance structure.
# =============================================================================
# Within-sex centering for BOTH CVR measures
# FRS and CVR_mimic arrive globally z-scored (M=0, SD=1)
# from script 03. Here we additionally center within sex
# to remove between-sex mean differences. This ensures
# main effects of HVR are evaluated at the sex-specific
# mean CVR, which is more interpretable for sex-stratified
# LGCM analyses.
lme.dt[, FRS_sex_centered := FRS - mean(
  FRS, na.rm = TRUE
), by = SEX]

# Verify centering
frs_by_sex.dt <- lme.dt[, .(
  n = .N,
  frs_mean = mean(FRS, na.rm = TRUE),
  frs_sd = sd(FRS, na.rm = TRUE),
  frs_centered_mean = mean(
    FRS_sex_centered, na.rm = TRUE
  )
), by = SEX]
log_info("FRS by sex (z-scored, pre-centering):")
log_info("  Male: mean=%.3f, SD=%.3f",
         frs_by_sex.dt[SEX == "Male", frs_mean],
         frs_by_sex.dt[SEX == "Male", frs_sd])
log_info("  Female: mean=%.3f, SD=%.3f",
         frs_by_sex.dt[SEX == "Female", frs_mean],
         frs_by_sex.dt[SEX == "Female", frs_sd])
log_info("  (FRS_sex_centered mean ~ 0 within each sex)")

# CVR_mimic: same within-sex centering as FRS
# This parallel treatment ensures both CVR measures
# receive identical preprocessing for LGCM mediation.
if (HAS_CVR_MIMIC) {
  lme.dt[, CVR_sex_centered := CVR_mimic - mean(
    CVR_mimic, na.rm = TRUE
  ), by = SEX]

  cvr_by_sex.dt <- lme.dt[, .(
    cvr_mean = mean(CVR_mimic, na.rm = TRUE),
    cvr_sd = sd(CVR_mimic, na.rm = TRUE),
    cvr_centered_mean = mean(
      CVR_sex_centered, na.rm = TRUE
    )
  ), by = SEX]
  log_info("CVR_mimic by sex (z-scored, pre-centering):")
  log_info("  Male: mean=%.3f, SD=%.3f",
           cvr_by_sex.dt[SEX == "Male", cvr_mean],
           cvr_by_sex.dt[SEX == "Male", cvr_sd])
  log_info("  Female: mean=%.3f, SD=%.3f",
           cvr_by_sex.dt[SEX == "Female", cvr_mean],
           cvr_by_sex.dt[SEX == "Female", cvr_sd])
  log_info(
    "  (CVR_sex_centered mean ~ 0 within each sex)"
  )
}

# Compute SE^2 (variance) for each domain
for (dom in c("MEM", "LAN", "EXF")) {
  se_col_name <- paste0("PHC_", dom, "_SE")
  var_col_name <- paste0("PHC_", dom, "_VAR")
  if (se_col_name %in% names(lme.dt)) {
    lme.dt[, (var_col_name) := get(se_col_name)^2]
  }
}

# =============================================================================
# 4. Create Wide Format
# =============================================================================
log_section("Creating Wide Format")

# Determine valid timepoints (those with sufficient data)
tp_counts.dt <- lme.dt[, .N, by = TIMEPOINT]
setorder(tp_counts.dt, TIMEPOINT)
log_info("Observations per timepoint:")
for (i in seq_len(nrow(tp_counts.dt))) {
  log_info("  %s: %d", tp_counts.dt$TIMEPOINT[i], tp_counts.dt$N[i])
}

# Use configurable number of timepoints (default: 6 = T0-T5) for LGCM stability
MAX_TIMEPOINTS <- get_script_setting("lgcm", "max_timepoints", default = 6)
valid_tps.v <- paste0("T", 0:(MAX_TIMEPOINTS - 1))
lme.dt <- lme.dt[TIMEPOINT %in% valid_tps.v]
log_info("Using %d timepoints: %s",
         MAX_TIMEPOINTS, paste(valid_tps.v, collapse = ", "))

# Get baseline values (time-invariant covariates)
# Note: HVR_Z is pivoted to HVR_Z_T0, HVR_Z_T1, etc. - don't include here
# Check for APOE4 availability
# Baseline covariates (time-invariant)
# Both FRS and CVR_mimic included when available,
# each with raw (z-scored) and sex-centered versions.
HAS_APOE4 <- "APOE4" %in% names(lme.dt) &&
  sum(!is.na(lme.dt$APOE4)) > 0

# Build baseline column list dynamically
bl_cols.v <- c(
  "PTID", "SEX", "AGE", "EDUC",
  "FRS", "FRS_sex_centered"
)
if (HAS_APOE4) bl_cols.v <- c(bl_cols.v, "APOE4")
if (HAS_CVR_MIMIC) {
  bl_cols.v <- c(
    bl_cols.v, "CVR_mimic", "CVR_sex_centered"
  )
}
# MIMIC indicators for OpenMx embedding
MIMIC_INDICATORS <- c(
  "SBP_z", "HTN", "GLUCOSE_z",
  "CHOL_z", "CREAT_z"
)
avail_mimic.v <- MIMIC_INDICATORS[
  MIMIC_INDICATORS %in% names(lme.dt)
]
HAS_MIMIC_INDICATORS <- length(avail_mimic.v) ==
  length(MIMIC_INDICATORS)
if (HAS_MIMIC_INDICATORS) {
  bl_cols.v <- c(bl_cols.v, avail_mimic.v)
  log_info("MIMIC indicators available for LGCM")
} else {
  stop(sprintf(
    "Missing MIMIC indicators: %s",
    paste(setdiff(
      MIMIC_INDICATORS, avail_mimic.v
    ), collapse = ", ")
  ))
}
bl_cols.v <- c(bl_cols.v, "EDT_bl_years")

baseline.dt <- lme.dt[
  TIMEPOINT == "T0",
  ..bl_cols.v
]
setnames(
  baseline.dt,
  c("AGE", "EDT_bl_years"),
  c("AGE_bl", "EDT_bl")
)

log_info("Baseline covariates: %s",
         paste(names(baseline.dt), collapse = ", "))

# Time-varying columns to pivot (all brain measures for parallel process LGCM)
time_varying.v <- c(
  "PHC_MEM", "PHC_LAN", "PHC_EXF",
  "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE",
  "PHC_MEM_VAR", "PHC_LAN_VAR", "PHC_EXF_VAR",
  "HVR_Z", "HC_Z", "HVR_RAW", "HC_RAW", "TIV",  # All brain measures at each tp
  "YRS_from_bl",
  "EDT"  # EDT as time metric for OpenMx ITVS
)

# Pivot to wide format
wide.dt <- dcast(
  lme.dt,
  PTID ~ TIMEPOINT,
  value.var = time_varying.v,
  sep = "_"
)

# Merge with baseline
wide.dt <- merge(baseline.dt, wide.dt, by = "PTID")

log_info("Wide format: %d subjects, %d columns", nrow(wide.dt), ncol(wide.dt))

# =============================================================================
# 4b. Add EDT² Columns for Quadratic LGCM
# =============================================================================
# Quadratic LGCM requires EDT² as definition variable loadings for the
# quadratic slope factor. See Sterba (2014) for methodology.
log_info("Adding EDT² columns for quadratic LGCM...")

edt_cols.v <- grep("^EDT_T[0-9]+$", names(wide.dt), value = TRUE)
for (edt_col in edt_cols.v) {
  sq_col <- paste0(edt_col, "_sq")
  wide.dt[, (sq_col) := get(edt_col)^2]
}

log_info("  Added %d EDT² columns: %s",
         length(edt_cols.v),
         paste(paste0(edt_cols.v, "_sq"), collapse = ", "))

# =============================================================================
# 5. Compute Derived Variables for LGCM
# =============================================================================
log_section("Computing Derived Variables")
# Count timepoints per subject (cognitive and HVR)
mem_cols.v <- grep("^PHC_MEM_T[0-9]+$", names(wide.dt), value = TRUE)
hvr_cols.v <- grep("^HVR_Z_T[0-9]+$", names(wide.dt), value = TRUE)
wide.dt[, n_cog_timepoints := rowSums(!is.na(.SD)), .SDcols = mem_cols.v]
wide.dt[, n_hvr_timepoints := rowSums(!is.na(.SD)), .SDcols = hvr_cols.v]

log_info("Timepoint counts:")
log_info("  Cognitive: %d subjects with >=2 timepoints",
         sum(wide.dt$n_cog_timepoints >= 2))
log_info("  HVR: %d subjects with >=2 timepoints",
         sum(wide.dt$n_hvr_timepoints >= 2))

# Create FRS × Brain interaction terms for ALL brain measures
# Using sex-centered FRS for interpretability

# HVR_Z interactions
wide.dt[, FRS_x_HVR := FRS_sex_centered * HVR_Z_T0]
wide.dt[, FRS_x_HVR_raw := FRS * HVR_Z_T0]

# HC_Z interactions
wide.dt[, FRS_x_HC := FRS_sex_centered * HC_Z_T0]
wide.dt[, FRS_x_HC_raw := FRS * HC_Z_T0]

# HVR_RAW interactions
wide.dt[, FRS_x_HVR_RAW := FRS_sex_centered * HVR_RAW_T0]

# HC_RAW interactions (also include TIV for covariate)
wide.dt[, FRS_x_HC_RAW := FRS_sex_centered * HC_RAW_T0]

# Standardize brain measures for interaction (optional)
wide.dt[, HVR_Z_T0_std := scale(HVR_Z_T0)[, 1]]
wide.dt[, HC_Z_T0_std := scale(HC_Z_T0)[, 1]]
wide.dt[, HVR_RAW_T0_std := scale(HVR_RAW_T0)[, 1]]
wide.dt[, HC_RAW_T0_std := scale(HC_RAW_T0)[, 1]]
wide.dt[, FRS_std := scale(FRS)[, 1]]
wide.dt[, FRS_x_HVR_std := FRS_std * HVR_Z_T0_std]
wide.dt[, FRS_x_HC_std := FRS_std * HC_Z_T0_std]

# CVR_mimic interaction terms
wide.dt[
  , CVR_x_HVR := CVR_sex_centered * HVR_Z_T0
]
wide.dt[
  , CVR_x_HC := CVR_sex_centered * HC_Z_T0
]
log_info(
  "  CVR_x_HVR: CVR (sex-ctr) x HVR_Z_T0"
)
log_info(
  "  CVR_x_HC: CVR (sex-ctr) x HC_Z_T0"
)

log_info("Derived variables created:")
log_info("  FRS_x_HVR: FRS (sex-centered) x HVR_Z_T0")
log_info("  FRS_x_HC: FRS (sex-centered) x HC_Z_T0")
log_info("  FRS_x_HVR_RAW: FRS (sex-centered) x HVR_RAW_T0")
log_info("  FRS_x_HC_RAW: FRS (sex-centered) x HC_RAW_T0")
log_info("  n_cog/n_hvr_timepoints: parallel process counts")

# Compute average SE^2 per domain (for model comparison)
for (dom in c("MEM", "LAN", "EXF")) {
  var_cols.v <- grep(
    paste0("PHC_", dom, "_VAR_T"), names(wide.dt), value = TRUE
  )
  if (length(var_cols.v) > 0) {
    avg_col_name <- paste0("PHC_", dom, "_VAR_AVG")
    wide.dt[,
            (avg_col_name) := rowMeans(.SD, na.rm = TRUE),
            .SDcols = var_cols.v]
  }
}

# =============================================================================
# 6. Create Analysis Samples
# =============================================================================
log_section("Creating Analysis Samples")

# Full sample: FRS + HVR + at least 2 timepoints (cog AND HVR)
has_complete.v <- !is.na(wide.dt$FRS) & !is.na(wide.dt$HVR_Z_T0) &
                  wide.dt$n_cog_timepoints >= 2 & wide.dt$n_hvr_timepoints >= 2

full_sample.dt <- wide.dt[has_complete.v]
log_info("Full sample (FRS + HVR + >=2 timepoints): %d subjects",
         nrow(full_sample.dt))

# Sex-stratified samples
male_sample.dt <- full_sample.dt[SEX == "Male"]
female_sample.dt <- full_sample.dt[SEX == "Female"]

log_info("Male sample: %d subjects", nrow(male_sample.dt))
log_info("Female sample: %d subjects", nrow(female_sample.dt))

# =============================================================================
# 7. Sample Summaries
# =============================================================================
log_section("Sample Summaries")

summarize_sample.fn <- function(dt, name) {
  log_info("%s (N=%d):", name, nrow(dt))
  log_info("  Age (baseline): %.1f (%.1f)",
           mean(dt$AGE_bl, na.rm = TRUE),
           sd(dt$AGE_bl, na.rm = TRUE))
  if (name == "Full Sample") {
    log_info("  Sex: %d Female (%.1f%%)",
             sum(dt$SEX == "Female", na.rm = TRUE),
             100 * mean(dt$SEX == "Female", na.rm = TRUE))
  }
  log_info("  Education: %.1f (%.1f)",
           mean(dt$EDUC, na.rm = TRUE),
           sd(dt$EDUC, na.rm = TRUE))
  log_info("  FRS: %.1f (%.1f)",
           mean(dt$FRS, na.rm = TRUE),
           sd(dt$FRS, na.rm = TRUE))
  log_info("  HVR z-score (T0): %.2f (%.2f)",
           mean(dt$HVR_Z_T0, na.rm = TRUE),
           sd(dt$HVR_Z_T0, na.rm = TRUE))
  log_info("  EDT (baseline): %.1f (%.1f)",
           mean(dt$EDT_bl, na.rm = TRUE),
           sd(dt$EDT_bl, na.rm = TRUE))
  log_info("  Cog timepoints: %.1f (range %d-%d)",
           mean(dt$n_cog_timepoints),
           min(dt$n_cog_timepoints),
           max(dt$n_cog_timepoints))
  log_info("  HVR timepoints: %.1f (range %d-%d)",
           mean(dt$n_hvr_timepoints),
           min(dt$n_hvr_timepoints),
           max(dt$n_hvr_timepoints))
}

summarize_sample.fn(full_sample.dt, "Full Sample")
summarize_sample.fn(male_sample.dt, "Male Sample")
summarize_sample.fn(female_sample.dt, "Female Sample")

# =============================================================================
# 8. Add Metadata and Save
# =============================================================================
log_section("Saving Results")

add_metadata.fn <- function(dt, name) {
  attr(dt, "creation_date") <- Sys.time()
  attr(dt, "n_subjects") <- nrow(dt)
  attr(dt, "valid_timepoints") <- valid_tps.v
  attr(dt, "time_metric") <- "EDT"
  attr(dt, "time_metric_note") <- paste(
    "EDT in years; negative = pre-onset, positive = post-onset"
  )
  attr(dt, "frs_centering") <- paste(
    "FRS_sex_centered: within-sex centered"
  )
  int_terms <- c("FRS_x_HVR", "FRS_x_HVR_std")
  if (HAS_CVR_MIMIC) {
    attr(dt, "cvr_centering") <- paste(
      "CVR_sex_centered: within-sex centered"
    )
  }
  attr(dt, "interaction_terms") <- int_terms
  attr(dt, "has_cvr_mimic") <- HAS_CVR_MIMIC
  attr(dt, "has_mimic_indicators") <-
    HAS_MIMIC_INDICATORS
  attr(dt, "sample_name") <- name
  dt
}

full_sample.dt <- add_metadata.fn(full_sample.dt, "Full AB+ Sample with EDT")
male_sample.dt <- add_metadata.fn(male_sample.dt, "Male AB+ Sample")
female_sample.dt <- add_metadata.fn(female_sample.dt, "Female AB+ Sample")

write_rds_safe(
  full_sample.dt, output_full_path, description = "LGCM EDT full sample"
)
write_rds_safe(
  male_sample.dt, output_male_path, description = "LGCM EDT male sample"
)
write_rds_safe(
  female_sample.dt, output_female_path, description = "LGCM EDT female sample"
)

log_info("")
log_info("=== OUTPUT FILES ===")
log_info("  Full sample: %s", output_full_path)
log_info("  Male sample: %s", output_male_path)
log_info("  Female sample: %s", output_female_path)

log_script_end("08_prepare_lgcm_data.R", success = TRUE)
