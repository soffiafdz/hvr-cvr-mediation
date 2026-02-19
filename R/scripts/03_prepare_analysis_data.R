#!/usr/bin/env Rscript

# =============================================================================
# 03_prepare_analysis_data.R
# =============================================================================
# Prepare base analysis cohort: z-scores, cognition, biomarkers,
# amyloid filtering, disease time (EDT). FRS kept as raw values.
#
# CVR_mimic merge and global standardization of BOTH CVR
# measures happens in 04_cvr_latent_factor.R (Part 6).
#
# Inputs (paths from config/pipeline_config.yaml):
#   - derivatives/z_scores: z-scores (from 02_compute_zscores.R)
#   - raw/adsp_phc_*: ADSP-PHC data (cognition, amyloid, etc.)
#   - raw/apoe: APOE genotype data
#
# Outputs:
#   - derivatives/lme_cohort_hvr: Base cohort (FRS raw, no CVR)
#   - derivatives/tau_power_info: Tau subsample power info
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(lubridate)
})

# Source utilities
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/disease_time.R"))

# Initialize
log_script_start("03_prepare_analysis_data.R")
config <- load_config()
set_seed()

# --- Configuration ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "prepare_data", default = FALSE
)

output_hvr_path <- get_data_path("derivatives", "lme_cohort_hvr")

# Check if output exists and skip if not forcing
if (!FORCE_REGENERATE && file.exists(output_hvr_path)) {
  log_info("Output file exists and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate = TRUE to rerun.")
  log_script_end("03_prepare_analysis_data.R", success = TRUE)
  quit(status = 0)
}

ensure_directory(dirname(output_hvr_path))

# =============================================================================
# 1. Load Z-Scores
# =============================================================================
log_section("Loading HVR Z-Scores")

zscore_path <- get_data_path("derivatives", "z_scores")
check_files_exist(zscore_path)
z_scores.dt <- read_rds_safe(zscore_path, "HVR z-scores")
log_info("Loaded z-scores: %d rows", nrow(z_scores.dt))

# Filter to HVR, bilateral (summed L+R), no adjustment
# SIDE == "LR" matches corrected z-score computation
# (GAMLSS models trained on summed bilateral volumes)
# ADJ == "NON" because HVR ratio inherently normalizes head-size
hvr.dt <- z_scores.dt[
  ROI == "HVR" & SIDE == "LR" & ADJ == "NON",
  .(PTID, EXAMDATE, MRI_date = EXAMDATE, HVR_z = Z)
]
setkey(hvr.dt, PTID, EXAMDATE)
log_info("HVR z-scores: %d observations from %d subjects",
         nrow(hvr.dt), hvr.dt[, uniqueN(PTID)])

# SENSITIVITY: Extract HC z-scores (bilateral, NON adjusted)
hc.dt <- z_scores.dt[
  ROI == "HC" & SIDE == "LR" & ADJ == "NON",
  .(PTID, EXAMDATE, HC_z = Z)
]
setkey(hc.dt, PTID, EXAMDATE)
log_info("HC z-scores (NON, sensitivity): %d observations", nrow(hc.dt))

# SENSITIVITY: Load raw HVR values (not z-scored)
hchvr_path <- get_data_path("derivatives", "adni_hchvr")
hchvr.dt <- read_rds_safe(hchvr_path, "HC-HVR volumes")

# Compute bilateral HVR from raw volumes (NON adjustment)
hvr_raw.dt <- hchvr.dt[ADJ == "NON", .(
  HVR_raw = sum(HC) / (sum(HC) + sum(VC))
), by = .(PTID, EXAMDATE)]

# Standardize HVR_raw (within-sample z-score) for comparable effect sizes
hvr_raw.dt[, HVR_raw := scale(HVR_raw)[, 1]]
setkey(hvr_raw.dt, PTID, EXAMDATE)
log_info("Raw HVR standardized (sensitivity): %d obs, M=%.2f, SD=%.2f",
         nrow(hvr_raw.dt),
         hvr_raw.dt[, mean(HVR_raw, na.rm = TRUE)],
         hvr_raw.dt[, sd(HVR_raw, na.rm = TRUE)])

# SENSITIVITY S3: Raw HC (unadjusted) with TIV covariate
# Extract bilateral HC (sum L+R) and TIV from NON-adjusted data
hc_raw.dt <- hchvr.dt[ADJ == "NON", .(
  HC_raw = sum(HC),
  TIV = unique(ICC)
), by = .(PTID, EXAMDATE)]

# Standardize HC_raw for comparable effect sizes
hc_raw.dt[, HC_raw := scale(HC_raw)[, 1]]
hc_raw.dt[, TIV := scale(TIV)[, 1]]
setkey(hc_raw.dt, PTID, EXAMDATE)
log_info("Raw HC + TIV (sensitivity): %d obs", nrow(hc_raw.dt))

# Merge sensitivity measures into hvr.dt
hvr.dt <- merge(hvr.dt, hc.dt, by = c("PTID", "EXAMDATE"), all.x = TRUE)
hvr.dt <- merge(hvr.dt, hvr_raw.dt, by = c("PTID", "EXAMDATE"), all.x = TRUE)
hvr.dt <- merge(hvr.dt, hc_raw.dt, by = c("PTID", "EXAMDATE"), all.x = TRUE)
log_info("Merged sensitivity: HC_z=%d, HVR_raw=%d, HC_raw=%d non-NA",
         sum(!is.na(hvr.dt$HC_z)),
         sum(!is.na(hvr.dt$HVR_raw)),
         sum(!is.na(hvr.dt$HC_raw)))

# =============================================================================
# 2. Load ADSP-PHC Data
# =============================================================================
log_section("Loading ADSP-PHC Data")

adsp_phc_dir <- get_data_path("raw", "adsp_phc")

# =============================================================================
# Amyloid Positivity: Combined PET + CSF Definition
# =============================================================================
# Amyloid positivity determined from ADSP-PHC data using two modalities:
# 1. PET imaging: >25 Centiloids (PHC_AMYLOID_STATUS == 1)
# 2. CSF biomarkers: A+ from AT_class (GMM-derived cutoffs, Timsina et al. 2024)
# =============================================================================

# --- PET Amyloid ---
amy_file <- get_data_path("raw", "adsp_phc_amyloid")
check_files_exist(amy_file)
amy.dt <- read_csv_safe(amy_file,
  select = c("PTID", "PHC_SCANDATE", "PHC_Sex", "PHC_AMYLOID_STATUS"),
  description = "Amyloid PET"
)
setkey(amy.dt, PTID, PHC_SCANDATE)
setnames(amy.dt, c("PHC_Sex", "PHC_AMYLOID_STATUS"), c("SEX", "AMY_stat"))

# Validate and recode sex variable
# PHC files use numeric coding (1=Male, 2=Female) while ADNIMERGE uses strings
# Check which coding scheme is present before recoding
sex_values <- unique(amy.dt$SEX)
if (all(sex_values %in% c(1, 2, NA_integer_, NA_real_))) {
  # Numeric coding (PHC format): 1=Male, 2=Female
  amy.dt[, SEX := fifelse(SEX == 1, "Male", "Female")]
  log_info("Sex recoded from numeric (1=Male, 2=Female) to character")
} else if (all(sex_values %in% c("Male", "Female", "M", "F", NA_character_))) {
  # Already character coding
  amy.dt[SEX %in% c("M"), SEX := "Male"]
  amy.dt[SEX %in% c("F"), SEX := "Female"]
  log_info("Sex variable already in character format")
} else {
  # Unexpected coding - log warning and show values
  sex_str <- paste(sex_values, collapse = ", ")
  log_warn("Unexpected SEX coding detected: %s", sex_str)
  log_warn("Please verify sex variable coding in source data")
  log_warn("Assuming 1=Male, 2=Female and proceeding...")
  amy.dt[, SEX := fifelse(SEX == 1, "Male", "Female")]
}

# Final validation
final_sex_values <- unique(amy.dt$SEX)
valid_sex_values <- c("Male", "Female", NA_character_)
invalid_values <- setdiff(final_sex_values, valid_sex_values)
if (length(invalid_values) > 0) {
  stop("Sex variable has unexpected values after recoding: ",
       paste(invalid_values, collapse = ", "),
       ". Please check source data coding.")
}
log_info("Sex validation passed: %d Male, %d Female",
         sum(amy.dt$SEX == "Male", na.rm = TRUE),
         sum(amy.dt$SEX == "Female", na.rm = TRUE))

log_info("Amyloid PET: %d scans, %d subjects",
         nrow(amy.dt), amy.dt[, uniqueN(PTID)])

# --- CSF Amyloid ---
# Load CSF biomarker data with AT_class (GMM-derived A/T classification)
# AT_class uses data-driven cutoffs via GMM (Timsina et al. 2024)
# Categories: A+T+, A+T-, A-T+, A-T-
csf_file <- get_data_path("raw", "adsp_phc_biomarker")
if (file.exists(csf_file)) {
  csf.dt <- read_csv_safe(csf_file,
    select = c("PTID", "PHC_Sex", "AT_class"),
    description = "CSF Biomarkers"
  )

  # Identify CSF A+ subjects (AT_class starts with "A+")
  csf.dt[, CSF_Apos := grepl("^A\\+", AT_class)]
  csf_apos.dt <- unique(csf.dt[CSF_Apos == TRUE, .(PTID)])

  # Recode sex from CSF file
  csf_sex_values <- unique(csf.dt$PHC_Sex)
  if (all(csf_sex_values %in% c(1, 2, NA_integer_, NA_real_))) {
    csf.dt[, SEX := fifelse(PHC_Sex == 1, "Male", "Female")]
  } else {
    csf.dt[, SEX := PHC_Sex]
  }
  csf_sex.dt <- unique(csf.dt[CSF_Apos == TRUE, .(PTID, SEX)])

  log_info("CSF A+ subjects: %d (from AT_class: %s)",
           nrow(csf_apos.dt),
           paste(unique(csf.dt[CSF_Apos == TRUE, AT_class]), collapse = ", "))
} else {
  log_warn("CSF biomarker file not found: %s", csf_file)
  log_warn("Proceeding with PET-only amyloid positivity")
  csf_apos.dt <- data.table(PTID = character(0))
  csf_sex.dt <- data.table(PTID = character(0), SEX = character(0))
}

# --- Tau PET (for power analysis only) ---
tau_file <- get_data_path("raw", "adsp_phc_tau")
tau_available <- file.exists(tau_file)
if (tau_available) {
  tau.dt <- read_csv_safe(tau_file,
    select = c("PTID", "PHC_SCANDATE", "META_TEMPORAL_SUVR"),
    description = "Tau PET"
  )
  setkey(tau.dt, PTID, PHC_SCANDATE)
  setnames(tau.dt, "META_TEMPORAL_SUVR", "TAU")
  tau.dt[, TAU_date := PHC_SCANDATE]
  log_info("Tau PET: %d scans, %d subjects",
           nrow(tau.dt), tau.dt[, uniqueN(PTID)])
} else {
  log_warn("Tau PET file not found: %s", tau_file)
  tau.dt <- data.table(PTID = character(0))
}

# Cognition
cog_file <- get_data_path("raw", "adsp_phc_cognition")
check_files_exist(cog_file)
cog.dt <- read_csv_safe(cog_file,
  select = c("PTID", "EXAMDATE",
             "PHC_Age_Cognition", "PHC_Education",
             "PHC_MEM", "PHC_LAN", "PHC_EXF",
             "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE",
             "PHC_Diagnosis"),
  description = "Cognition"
)
setkey(cog.dt, PTID)
setnames(cog.dt,
         c("PHC_Age_Cognition", "PHC_Education",
           "PHC_MEM", "PHC_LAN", "PHC_EXF",
           "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE",
           "PHC_Diagnosis"),
         c("AGE", "EDUC",
           "MEM", "LAN", "EXF",
           "MEM_SE", "LAN_SE", "EXF_SE",
           "DX"))
n_cog_database <- cog.dt[, uniqueN(PTID)]
log_info("Cognition: %d observations, %d subjects",
         nrow(cog.dt), n_cog_database)

# CVRF (Framingham Risk Score and PC1)
FRS_COLUMN <- get_parameter("frs", "column")
PC1_COLUMN <- get_parameter("cvrf_pc1", "column")
cvrf_file <- get_data_path("raw", "adsp_phc_cvrf")
check_files_exist(cvrf_file)
cvrf.dt <- read_csv_safe(cvrf_file,
  select = c("RID", "PTID", FRS_COLUMN, PC1_COLUMN),
  description = "CVRF"
)
setkey(cvrf.dt, PTID)
setnames(cvrf.dt, c(FRS_COLUMN, PC1_COLUMN), c("FRS", "CVRF_PC1"))
cvrf.dt <- cvrf.dt[!is.na(FRS)]  # Keep if FRS available (PC1 may have some NA)
log_info("CVRF (FRS): %d subjects with FRS data", nrow(cvrf.dt))
log_info("CVRF (PC1): %d subjects with PC1 data", sum(!is.na(cvrf.dt$CVRF_PC1)))

# IMPORTANT: FRS was validated for ages 30-62 (D'Agostino et al., 2008)
# ADNI participants are typically older. Log this limitation.
# Age check will be done after merge with cognition data (which has AGE)

# APOE Genotype
apoe_file <- get_data_path("raw", "apoe")
check_files_exist(apoe_file)
apoe.dt <- read_csv_safe(apoe_file,
  select = c("PTID", "GENOTYPE"),
  description = "APOE genotype"
)
setkey(apoe.dt, PTID)
# Extract APOE4 carrier status from genotype string (e.g., "3/4" -> TRUE)
apoe.dt[, APOE4 := grepl("4", GENOTYPE)]
apoe.dt <- unique(apoe.dt[, .(PTID, APOE4)])  # One row per subject
log_info("APOE: %d subjects with genotype data", nrow(apoe.dt))
log_info("  APOE4 carriers: %d (%.1f%%)",
         sum(apoe.dt$APOE4), 100 * mean(apoe.dt$APOE4))

# =============================================================================
# 3. Filter and Merge
# =============================================================================
log_section("Filtering and Merging")

# =============================================================================
# Combine PET A+ and CSF A+ for amyloid-positive definition
# =============================================================================
# PET A+: PHC_AMYLOID_STATUS == 1 (>25 Centiloids)
# CSF A+: AT_class starts with "A+" (GMM-derived cutoffs)
# Final: Union of PET A+ OR CSF A+

amy_pet_pos.dt <- unique(
  amy.dt[AMY_stat == 1, .(PTID, SEX, AMY_source = "PET")]
)
amy_csf_pos.dt <- csf_sex.dt[, .(PTID, SEX, AMY_source = "CSF")]

# Combine: prefer PET if subject has both, otherwise use CSF
amy_pos.dt <- rbind(amy_pet_pos.dt, amy_csf_pos.dt, fill = TRUE)
amy_pos.dt <- amy_pos.dt[, .SD[1], by = PTID]  # Keep first (PET prioritized)
amy_pos.dt <- unique(amy_pos.dt[, .(PTID, SEX)])

# Log contribution from each source
n_pet_only <- length(setdiff(amy_pet_pos.dt$PTID, amy_csf_pos.dt$PTID))
n_csf_only <- length(setdiff(amy_csf_pos.dt$PTID, amy_pet_pos.dt$PTID))
n_both <- length(intersect(amy_pet_pos.dt$PTID, amy_csf_pos.dt$PTID))

log_info("Amyloid-positive subjects (combined): %d", nrow(amy_pos.dt))
log_info("  PET A+ only: %d", n_pet_only)
log_info("  CSF A+ only: %d", n_csf_only)
log_info("  Both PET and CSF A+: %d", n_both)

# Add FRS and PC1 to cognition data
cog.dt <- cvrf.dt[cog.dt, on = "PTID", nomatch = NULL]
log_info("After FRS merge: %d observations, %d subjects",
         nrow(cog.dt), cog.dt[, uniqueN(PTID)])

# Add APOE to cognition data
# Left join - keep all cog even if no APOE
cog.dt <- apoe.dt[cog.dt, on = "PTID"]
n_apoe_missing <- sum(is.na(cog.dt$APOE4))
if (n_apoe_missing > 0) {
  pct_missing <- 100 * n_apoe_missing / nrow(cog.dt)
  log_warn("APOE missing for %d observations (%.1f%%)",
           n_apoe_missing, pct_missing)
  log_warn("These will be excluded from APOE-stratified analyses")
}
log_info("After APOE merge: %d observations, %d subjects",
         nrow(cog.dt), cog.dt[, uniqueN(PTID)])

# Filter to amyloid-positive
cog.dt <- amy_pos.dt[cog.dt, on = "PTID", nomatch = NULL]
log_info("After Aβ+ filter: %d observations, %d subjects",
         nrow(cog.dt), cog.dt[, uniqueN(PTID)])

# FRS age validity assessment
# Framingham Risk Score validated for ages 30-62 (D'Agostino et al., 2008)
# Record statistics for Limitations section of manuscript
FRS_MIN_AGE <- 30
FRS_MAX_AGE <- 62

frs_validity <- list()
if ("AGE" %in% names(cog.dt)) {
  # Helper: check if age outside validated range
  is_outside_range <- function(age) {
    age < FRS_MIN_AGE | age > FRS_MAX_AGE
  }

  frs_validity <- cog.dt[!is.na(FRS), .(
    n_total = .N,
    n_subjects = uniqueN(PTID),
    n_below_30 = sum(AGE < FRS_MIN_AGE, na.rm = TRUE),
    n_above_62 = sum(AGE > FRS_MAX_AGE, na.rm = TRUE),
    n_outside_range = sum(is_outside_range(AGE), na.rm = TRUE),
    pct_outside = 100 * sum(is_outside_range(AGE), na.rm = TRUE) / .N,
    mean_age = mean(AGE, na.rm = TRUE),
    sd_age = sd(AGE, na.rm = TRUE),
    min_age = min(AGE, na.rm = TRUE),
    max_age = max(AGE, na.rm = TRUE)
  )]

  # Compute at subject level (baseline only) for manuscript reporting
  baseline_ages <- cog.dt[!is.na(FRS), .SD[1], by = PTID][, AGE]
  n_outside <- sum(is_outside_range(baseline_ages), na.rm = TRUE)
  frs_validity$n_subjects_outside <- n_outside
  frs_validity$pct_subjects_outside <- 100 * n_outside / length(baseline_ages)

  log_info("FRS VALIDITY STATISTICS (for Limitations section):")
  log_info("  Validated age range: %d-%d years", FRS_MIN_AGE, FRS_MAX_AGE)
  log_info("  Sample age: %.1f ± %.1f years (range: %.0f-%.0f)",
           frs_validity$mean_age, frs_validity$sd_age,
           frs_validity$min_age, frs_validity$max_age)
  log_info("  Subjects outside validated range: %d/%d (%.1f%%)",
           frs_validity$n_subjects_outside, frs_validity$n_subjects,
           frs_validity$pct_subjects_outside)
  log_info("  Observations outside validated range: %d/%d (%.1f%%)",
           frs_validity$n_outside_range, frs_validity$n_total,
           frs_validity$pct_outside)
}

# Store for downstream use (e.g., manuscript generation)
attr(cog.dt, "frs_validity") <- frs_validity

# Format dates and diagnosis
cog.dt[, `:=`(
  COG_date = ymd(EXAMDATE),
  DX = factor(DX, levels = 1:3, labels = c("CU", "MCI", "AD"))
)]
setkey(cog.dt, PTID, EXAMDATE)

# =============================================================================
# 4. Create Analysis Cohort
# =============================================================================
# Note: cohort_hvr.dt contains ALL brain measures (HVR_z, HC_z, HVR_raw, HC_raw)
# for the same subjects. The "hvr" name is historical - same cohort used for all.
log_section("Creating Analysis Cohort")

# Rolling join: get closest preceding brain measurement
# Using roll = Inf to match OHBM methodology and maximize sample size
# hvr.dt already contains HC_z, HVR_raw, HC_raw merged in Section 1
cohort_hvr.dt <- hvr.dt[cog.dt, roll = Inf, nomatch = NULL]

# Log time gap diagnostics
if (nrow(cohort_hvr.dt) > 0) {
  # MRI_date preserved from hvr.dt; COG_date from cog.dt
  # Gap = cognitive assessment date - MRI scan date (should be >= 0)
  cohort_hvr.dt[, HVR_COG_gap_days := as.numeric(COG_date - MRI_date)]
  gap_summary <- cohort_hvr.dt[, .(
    mean_gap = mean(HVR_COG_gap_days, na.rm = TRUE),
    max_gap = max(HVR_COG_gap_days, na.rm = TRUE),
    n_over_180 = sum(HVR_COG_gap_days > 180, na.rm = TRUE)
  )]
  log_info("HVR-Cognition gap: mean=%.0f days, max=%.0f days, >180d: %d",
           gap_summary$mean_gap, gap_summary$max_gap,
           gap_summary$n_over_180)
}

# Add visit structure
cohort_hvr.dt[order(COG_date), `:=`(
  VISIT = seq_len(.N),
  YRS_from_bl = as.numeric(COG_date - COG_date[1]) / 365.25
), by = PTID]

# Check for subjects with only 1 visit (contribute no slope variance info)
n_single_visit <- cohort_hvr.dt[, .N, by = PTID][N == 1, .N]
if (n_single_visit > 0) {
  log_warn("HVR cohort: %d subjects have only 1 visit", n_single_visit)
  log_warn("These cannot contribute to random slope estimation")
  log_warn("Consider filtering or acknowledging in methods")
}

log_info("ANALYSIS COHORT (all brain measures):")
log_info("  Subjects: %d", cohort_hvr.dt[, uniqueN(PTID)])
log_info("  Observations: %d", nrow(cohort_hvr.dt))
log_info("  Visits/subject: %.1f (range: %d-%d)",
         cohort_hvr.dt[, .N, by = PTID][, mean(N)],
         cohort_hvr.dt[, .N, by = PTID][, min(N)],
         cohort_hvr.dt[, .N, by = PTID][, max(N)])
log_info("  Brain measures: HVR_z, HC_z, HVR_raw, HC_raw")

# =============================================================================
# 5. Tau Subsample Power Analysis
# =============================================================================
# Compute tau subsample size for post-hoc power analysis in 04a-d scripts.
# Same subjects used for all brain measures, so power info applies to all.
# Tau analyses excluded due to insufficient power for interaction effects.
log_section("Tau Subsample Power Analysis")

tau_power_info <- list(
  tau_available = tau_available,
  n_tau_subjects = 0,
  n_tau_observations = 0,
  n_analysis_subjects = cohort_hvr.dt[, uniqueN(PTID)],
  n_analysis_observations = nrow(cohort_hvr.dt),
  exclusion_reason = "Insufficient sample size for interaction effects",
  note = "Same subjects for all brain measures (HVR_z, HC_z, HVR_raw, HC_raw)"
)

if (tau_available && nrow(tau.dt) > 0) {
  # Compute what the tau subsample would be (intersection with analysis cohort)
  tau_subjects.v <- unique(tau.dt$PTID)
  analysis_subjects.v <- unique(cohort_hvr.dt$PTID)
  tau_overlap.v <- intersect(tau_subjects.v, analysis_subjects.v)

  # Join tau to cognition to get observation count
  cog_tau.dt <- tau.dt[cog.dt, roll = Inf, nomatch = NULL]
  cog_tau.dt[, EXAMDATE := as.IDate(COG_date)]
  setkey(cog_tau.dt, PTID, EXAMDATE)

  # Join brain measures to get final tau subsample
  tau_subsample.dt <- hvr.dt[cog_tau.dt, roll = Inf, nomatch = NULL]

  tau_power_info$n_tau_subjects <- length(tau_overlap.v)
  tau_power_info$n_tau_observations <- nrow(tau_subsample.dt)
  tau_power_info$sample_ratio <- tau_power_info$n_tau_subjects /
                                  tau_power_info$n_analysis_subjects

  log_info("TAU SUBSAMPLE (for power analysis):")
  log_info("  Tau PET subjects: %d", length(tau_subjects.v))
  log_info("  Overlap with analysis cohort: %d subjects", length(tau_overlap.v))
  log_info("  Tau subsample observations: %d", nrow(tau_subsample.dt))
  log_info("  Sample ratio (tau/analysis): %.1f%%",
           100 * tau_power_info$sample_ratio)

  # Power calculation note
  log_info("")
  log_info("POWER ANALYSIS NOTE:")
  log_info("  Analysis cohort N=%d provides adequate power for interactions",
           tau_power_info$n_analysis_subjects)
  log_info("  Tau subsample N=%d is %.1f%% of analysis cohort",
           tau_power_info$n_tau_subjects,
           100 * tau_power_info$sample_ratio)
  log_info("  SE scales with sqrt(N): SE_tau ≈ %.1fx SE_analysis",
           sqrt(tau_power_info$n_analysis_subjects /
                tau_power_info$n_tau_subjects))
  log_info("  Post-hoc power analysis in 04a-d scripts will quantify this")

  # Clean up temporary objects
  rm(cog_tau.dt, tau_subsample.dt)
} else {
  log_warn("Tau PET data not available - skipping tau power analysis")
}

# =============================================================================
# 6. Sample Summary
# =============================================================================
log_section("Sample Summary")

baseline.dt <- cohort_hvr.dt[VISIT == 1]
log_info("Analysis Cohort (N=%d):", nrow(baseline.dt))
log_info("  Age: %.1f (%.1f) [%.0f-%.0f]",
         mean(baseline.dt$AGE, na.rm = TRUE),
         sd(baseline.dt$AGE, na.rm = TRUE),
         min(baseline.dt$AGE, na.rm = TRUE),
         max(baseline.dt$AGE, na.rm = TRUE))
log_info("  Sex: %d Female (%.1f%%)",
         sum(baseline.dt$SEX == "Female"),
         100 * mean(baseline.dt$SEX == "Female"))
log_info("  Education: %.1f (%.1f)",
         mean(baseline.dt$EDUC, na.rm = TRUE),
         sd(baseline.dt$EDUC, na.rm = TRUE))
if ("APOE4" %in% names(baseline.dt)) {
  log_info("  APOE4 carriers: %d (%.1f%%)",
           sum(baseline.dt$APOE4, na.rm = TRUE),
           100 * mean(baseline.dt$APOE4, na.rm = TRUE))
}
log_info("  FRS: %.1f (%.1f) [%.0f-%.0f]",
         mean(baseline.dt$FRS, na.rm = TRUE),
         sd(baseline.dt$FRS, na.rm = TRUE),
         min(baseline.dt$FRS, na.rm = TRUE),
         max(baseline.dt$FRS, na.rm = TRUE))
if ("CVRF_PC1" %in% names(baseline.dt)) {
  log_info("  CVRF PC1: %.2f (%.2f)",
           mean(baseline.dt$CVRF_PC1, na.rm = TRUE),
           sd(baseline.dt$CVRF_PC1, na.rm = TRUE))
}
log_info("  HVR z-score: %.2f (%.2f)",
         mean(baseline.dt$HVR_z, na.rm = TRUE),
         sd(baseline.dt$HVR_z, na.rm = TRUE))
log_info("  Diagnosis: CU=%d, MCI=%d, AD=%d",
         sum(baseline.dt$DX == "CU"),
         sum(baseline.dt$DX == "MCI"),
         sum(baseline.dt$DX == "AD"))

# =============================================================================
# 7. Compute Estimated Disease Time (EDT)
# =============================================================================
log_section("Computing Estimated Disease Time (EDT)")

# Compute EDT for analysis cohort using progmod
# Multivariate ADAS13 + MMSE model (Kühnel et al., Stat Med 2021)
# Falls back to univariate ADAS-13 (Raket, Front Big Data 2020) if needed
# This avoids circularity with harmonized cognitive outcomes (MEM, LAN, EXF)
log_info("Computing EDT for analysis cohort...")
cohort_hvr.dt <- compute_disease_time(
  cohort_hvr.dt,
  config = config,
  use_precomputed = TRUE,
  fit_multivariate = TRUE
)

# Validate EDT
edt_valid <- validate_edt(cohort_hvr.dt)

if (edt_valid) {
  log_info("Analysis cohort EDT: range [%.2f, %.2f] years",
           min(cohort_hvr.dt$EDT, na.rm = TRUE),
           max(cohort_hvr.dt$EDT, na.rm = TRUE))
  attr(cohort_hvr.dt, "edt_available") <- TRUE
} else {
  log_warn("Analysis cohort EDT not available or invalid")
  attr(cohort_hvr.dt, "edt_available") <- FALSE
}

# =============================================================================
# 7b. CONSORT Inclusion Counts
# =============================================================================
# Save inclusion counts at each filtering step for CONSORT diagram
# These must come from data — never hardcoded in manuscript
log_section("CONSORT Inclusion Counts")

inclusion_counts <- list(
  n_database = n_cog_database,
  n_amyloid_pos = nrow(amy_pos.dt),
  n_with_mri = length(
    intersect(amy_pos.dt$PTID, unique(hvr.dt$PTID))
  ),
  n_analysis = cohort_hvr.dt[, uniqueN(PTID)]
)
log_info("CONSORT counts:")
log_info("  Database: %d", inclusion_counts$n_database)
log_info("  Amyloid+: %d", inclusion_counts$n_amyloid_pos)
log_info("  With MRI: %d", inclusion_counts$n_with_mri)
log_info("  Analysis: %d", inclusion_counts$n_analysis)

# =============================================================================
# 8. Save Outputs
# =============================================================================
log_section("Saving Results")

# Transfer FRS validity stats to output cohort
attr(cohort_hvr.dt, "frs_validity") <- frs_validity

write_rds_safe(cohort_hvr.dt, output_hvr_path,
               description = "LME analysis cohort (all brain measures)")

# Save tau power info for downstream LME scripts
tau_power_path <- get_data_path("derivatives", "tau_power_info")
write_rds_safe(tau_power_info, tau_power_path,
               description = "Tau subsample power info")

# Save FRS validity separately for easy manuscript access
frs_validity_path <- path_join(dirname(output_hvr_path), "frs_validity_stats.rds")
write_rds_safe(frs_validity, frs_validity_path,
               description = "FRS validity statistics")

# Save CONSORT inclusion counts
inclusion_path <- path_join(
  dirname(output_hvr_path), "inclusion_counts.rds"
)
write_rds_safe(inclusion_counts, inclusion_path,
               description = "CONSORT inclusion counts")

log_info("Saved: %s", output_hvr_path)
log_info("Saved: %s", tau_power_path)
log_info("Saved: %s", frs_validity_path)
log_info("Saved: %s", inclusion_path)

log_script_end("03_prepare_analysis_data.R", success = TRUE)
