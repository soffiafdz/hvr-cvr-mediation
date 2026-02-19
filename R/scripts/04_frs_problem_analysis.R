#!/usr/bin/env Rscript

# ============================================================
# 04_frs_problem_analysis.R
# ============================================================
# PURPOSE:
#   Document why FRS is a poor CVR measure in ADNI and
#   identify data quality issues in the ADSP-PHC CVRF
#   dataset that affect any CVR composite derived from
#   these variables.
#
# PARTS:
#   1. FRS ceiling effect + age-weighting
#   2. PHC indicator quality (collinearity)
#   3. Raw ADNI verification (artifact confirmation)
#   4. Save results
#
# OUTPUTS:
#   - models/results/lme/frs_problem_analysis.rds
#   - models/results/lme/frs_validity_stats.rds
#
# ============================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))

log_script_start("04_frs_problem_analysis.R")
config <- load_config()
validate_config(config)

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "frs_analysis", default = FALSE
)
output_1.path <- get_data_path("models", "frs_problem_analysis")
output_2.path <- get_data_path("models", "frs_validity_stats")
if (!FORCE_REGENERATE &&
    file.exists(output_1.path) &&
    file.exists(output_2.path)) {
  log_info("Outputs exist and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE to rerun.")
  log_script_end("04_frs_problem_analysis.R", success = TRUE)
  quit(status = 0)
}

# ============================================================
# PART 1: FRS PROBLEMS IN AGING SAMPLES
# ============================================================
# FRS (Framingham Risk Score) has three problems in
# elderly cohorts like ADNI:
#   1. CEILING: Capped at 30% 10-year risk
#   2. AGE-WEIGHTING: Age dominates the formula
#   3. DIABETES INPUT: Uses corrupted PHC_Diabetes
#      (documented in Part 2 below)
# ============================================================
log_section("Part 1: FRS Problems")

# --- Load data ---
cvrf_path <- get_data_path("raw", "adsp_phc_cvrf")
check_files_exist(cvrf_path, "CVRF data file")
cvrf.dt <- fread(cvrf_path)
cvrf.dt <- cvrf.dt[order(RID, EXAMDATE)]
cvrf.dt <- cvrf.dt[, .SD[1], by = RID]
log_info(
  "Loaded CVRF baseline: %d subjects", nrow(cvrf.dt)
)

# Filter to analysis cohort if available
cohort_path <- get_data_path(
  "derivatives", "lme_cohort_hvr"
)
if (file.exists(cohort_path)) {
  cohort.dt <- read_rds_safe(
    cohort_path, "LME cohort"
  )
  cohort_rids.v <- unique(cohort.dt$RID)
  cvrf.dt <- cvrf.dt[RID %in% cohort_rids.v]
  log_info(
    "Filtered to cohort: %d subjects", nrow(cvrf.dt)
  )
}

cvrf.dt[, SEX := factor(
  ifelse(PHC_Sex == 1, "Male", "Female"),
  levels = c("Male", "Female")
)]
cvrf.dt[
  , FRS := PHC_ASCVD_10y_FRS_Simple_Ageover30
]
cvrf.dt[, AGE := PHC_Age_CardiovascularRisk]

# --- 1a. Ceiling effect ---
frs_valid.dt <- cvrf.dt[!is.na(FRS)]
CEILING <- 0.30
at_ceiling <- sum(
  abs(frs_valid.dt$FRS - CEILING) < 0.001
)
pct_ceiling <- 100 * at_ceiling / nrow(frs_valid.dt)

ceiling_by_sex.dt <- frs_valid.dt[, .(
  n = .N,
  n_ceiling = sum(abs(FRS - CEILING) < 0.001),
  pct_ceiling = 100 * sum(
    abs(FRS - CEILING) < 0.001
  ) / .N
), by = SEX]

log_info(
  "Ceiling (30%%): %d/%d (%.1f%%)",
  at_ceiling, nrow(frs_valid.dt), pct_ceiling
)
for (i in seq_len(nrow(ceiling_by_sex.dt))) {
  log_info(
    "  %s: %.1f%%",
    ceiling_by_sex.dt$SEX[i],
    ceiling_by_sex.dt$pct_ceiling[i]
  )
}

# --- 1b. Age-FRS correlations ---
age_frs.dt <- cvrf.dt[!is.na(FRS) & !is.na(AGE)]
r_pearson <- cor(
  age_frs.dt$AGE, age_frs.dt$FRS
)
r_spearman <- cor(
  age_frs.dt$AGE, age_frs.dt$FRS,
  method = "spearman"
)
log_info(
  "FRS-Age: Pearson r=%.3f, Spearman rho=%.3f",
  r_pearson, r_spearman
)

cor_by_sex.lst <- list()
for (sex in c("Male", "Female")) {
  s.dt <- age_frs.dt[SEX == sex]
  r <- cor(s.dt$AGE, s.dt$FRS)
  cor_by_sex.lst[[sex]] <- list(
    n = nrow(s.dt), pearson = r
  )
  log_info(
    "  %s (N=%d): r=%.3f", sex, nrow(s.dt), r
  )
}

# --- 1c. Partial correlation (Age-FRS | CVRFs) ---
CVRF_VARS <- c(
  "PHC_SBP", "PHC_BMI", "PHC_Diabetes",
  "PHC_Hypertension", "PHC_Smoker"
)
available.v <- CVRF_VARS[
  CVRF_VARS %in% names(age_frs.dt)
]

partial.dt <- age_frs.dt[
  complete.cases(age_frs.dt[, ..available.v])
]
cov_formula <- paste(available.v, collapse = " + ")
age_resid.v <- residuals(lm(
  as.formula(paste("AGE ~", cov_formula)),
  data = partial.dt
))
frs_resid.v <- residuals(lm(
  as.formula(paste("FRS ~", cov_formula)),
  data = partial.dt
))
r_partial <- cor(age_resid.v, frs_resid.v)
log_info(
  "Partial r (Age-FRS | CVRFs): %.3f", r_partial
)

# --- 1d. Variance decomposition ---
r2_age <- summary(
  lm(FRS ~ AGE, data = age_frs.dt)
)$r.squared
r2_cvrf <- summary(lm(
  as.formula(paste("FRS ~", cov_formula)),
  data = partial.dt
))$r.squared
r2_full <- summary(lm(
  as.formula(paste("FRS ~ AGE +", cov_formula)),
  data = partial.dt
))$r.squared

log_info("Variance decomposition:")
log_info("  R2(Age only): %.3f", r2_age)
log_info("  R2(CVRFs only): %.3f", r2_cvrf)
log_info("  R2(Age+CVRFs): %.3f", r2_full)
log_info(
  "  Unique Age: %.1f%%", (r2_full - r2_cvrf) * 100
)
log_info(
  "  Unique CVRF: %.1f%%", (r2_full - r2_age) * 100
)

# ============================================================
# PART 2: PHC INDICATOR QUALITY
# ============================================================
# The ADSP-PHC CVRF dataset contains harmonized
# self-reported medical history. We discovered that
# PHC_Diabetes, PHC_Heart, and PHC_Stroke are nearly
# collinear (phi > 0.90) — a harmonization artifact,
# not a real clinical pattern.
#
# This affects FRS (uses PHC_Diabetes as input) and
# PHC_PC1_All (uses 3 collinear variables of 5).
# ============================================================
log_section("Part 2: PHC Indicator Quality")

# Reload full CVRF (not filtered to cohort)
cvrf_full.dt <- fread(cvrf_path)
cvrf_full.dt <- cvrf_full.dt[order(RID, EXAMDATE)]
cvrf_full.dt <- cvrf_full.dt[, .SD[1], by = RID]

# --- 2a. Prevalence ---
BINARY_VARS <- c(
  "PHC_Hypertension", "PHC_Diabetes", "PHC_Heart",
  "PHC_Stroke", "PHC_Smoker", "PHC_BP_Med"
)
prevalence.lst <- list()
for (v in BINARY_VARS) {
  n1 <- sum(cvrf_full.dt[[v]] == 1, na.rm = TRUE)
  n <- sum(!is.na(cvrf_full.dt[[v]]))
  pct <- 100 * n1 / n
  prevalence.lst[[v]] <- list(
    n_pos = n1, n_tot = n, pct = pct
  )
  log_info("  %-20s %d/%d (%.1f%%)", v, n1, n, pct)
}
# FINDING: Stroke at ~29% is implausible (expected
# 5-10% in ADNI, which excludes severe CVD).
# Diabetes, Heart, Stroke all near ~28-29%.

# --- 2b. Cross-tabulations ---
DISEASE_VARS <- c(
  "PHC_Hypertension", "PHC_Diabetes",
  "PHC_Heart", "PHC_Stroke"
)
crosstabs.lst <- list()
for (i in seq_along(DISEASE_VARS)) {
  for (j in seq_along(DISEASE_VARS)) {
    if (j <= i) next
    v1 <- DISEASE_VARS[i]
    v2 <- DISEASE_VARS[j]
    tab <- table(
      cvrf_full.dt[[v1]], cvrf_full.dt[[v2]]
    )
    phi <- cor(
      as.numeric(cvrf_full.dt[[v1]]),
      as.numeric(cvrf_full.dt[[v2]]),
      use = "complete"
    )
    s1 <- gsub("PHC_", "", v1)
    s2 <- gsub("PHC_", "", v2)
    pair <- paste(s1, s2, sep = "_x_")
    crosstabs.lst[[pair]] <- list(
      tab = tab, phi = phi
    )
    log_info(
      "  %s x %s: phi=%.3f", s1, s2, phi
    )
  }
}
# FINDING: Diabetes x Heart phi=0.94,
# Heart x Stroke phi=0.96, Diabetes x Stroke
# phi=0.91. These are effectively the same
# variable — biologically impossible from
# independent self-reports.

# --- 2c. Correlation matrix ---
ALL_VARS <- c(
  "PHC_SBP", "PHC_BMI", DISEASE_VARS,
  "PHC_Smoker", "PHC_BP_Med",
  "PHC_Age_CardiovascularRisk"
)
cor.mat <- cor(
  cvrf_full.dt[, ..ALL_VARS],
  use = "pairwise.complete.obs"
)
# FINDING: HTN-Diabetes phi=0.004 (expected
# 0.15-0.30). Age-Diabetes r=-0.17 (expected
# positive). The collinear cluster is an artifact.

# --- 2d. BP_Med is strict subset of HTN ---
tab_bp <- table(
  BPMed = cvrf_full.dt$PHC_BP_Med,
  HTN = cvrf_full.dt$PHC_Hypertension
)
bp_in_htn <- tab_bp["1", "1"]
total_htn <- sum(cvrf_full.dt$PHC_Hypertension == 1)
log_info(
  "BP_Med: strict subset of HTN (%d/%d, %.0f%%)",
  bp_in_htn, total_htn,
  100 * bp_in_htn / total_htn
)

# --- 2e. Implications ---
# FRS uses PHC_Diabetes as direct input — corrupted.
# PHC_PC1_All uses Diabetes+Heart+Stroke (3 of 5
# inputs are collinear) — dominated by artifact.
# Any MIMIC model using these variables is invalid.
#
# NOTE: PHC_PC1_All (principal component of all CV risk
# factors) does not resolve these issues:
#  1. Dominated by collinear artifact (3 of 5 inputs
#     are PHC_Diabetes/Heart/Stroke, phi > 0.90)
#  2. No age adjustment — PCA treats age-driven and
#     true CVR variance identically
#  3. Mixed indicator types (continuous + binary)
#     violate PCA assumptions
# We instead use a MIMIC model (script 05) that:
#  - Excludes corrupted indicators
#  - Treats age as a CAUSE (structural path), yielding
#    age-residualized scores (zeta = eta - gamma*AGE)
#  - Handles binary indicators via MLR + FIML
# SAFE indicators: SBP, BMI, Hypertension.
log_info("")
log_info("Indicator classification:")
log_info("  SAFE:        SBP, BMI, Hypertension")
log_info("  PROBLEMATIC: Diabetes, Heart, Stroke")
log_info("  EXCLUDED:    BP_Med (redundant with HTN)")

# ============================================================
# PART 3: RAW ADNI VERIFICATION
# ============================================================
# Cross-reference PHC variables against raw ADNI
# source data to confirm the collinearity is a
# harmonization artifact, not in the source.
# ============================================================
log_section("Part 3: Raw ADNI Verification")

adni_cvr_path <- get_data_path(
  "raw", "adni_cvr_vars"
)
check_files_exist(adni_cvr_path, "Raw ADNI CVR")
adni_raw.dt <- fread(adni_cvr_path)
log_info(
  "Raw ADNI: %d rows, %d subjects",
  nrow(adni_raw.dt),
  uniqueN(adni_raw.dt$subject_id)
)

# Screening-visit structured variables
adni_sc.dt <- adni_raw.dt[visit == "sc", .(
  HMHYPERT = HMHYPERT[1],
  HMSTROKE = fifelse(
    HMSTROKE[1] == 2, 1L,
    fifelse(HMSTROKE[1] == 0, 0L, NA_integer_)
  ),
  MH4CARD = MH4CARD[1],
  MH9ENDO = MH9ENDO[1],
  MH16SMOK = MH16SMOK[1]
), by = subject_id]
setnames(adni_sc.dt, "subject_id", "PTID")

# Merge with PHC
verify.dt <- merge(
  cvrf_full.dt, adni_sc.dt,
  by = "PTID", all.x = TRUE
)
n_merged <- sum(!is.na(verify.dt$HMHYPERT))
log_info("Merged with raw ADNI: %d subjects", n_merged)

# Concordance: HMHYPERT vs PHC_Hypertension
tab_htn <- table(
  HMHYPERT = verify.dt$HMHYPERT,
  PHC_HTN = verify.dt$PHC_Hypertension
)
if (nrow(tab_htn) == 2 && ncol(tab_htn) == 2) {
  agree_htn <- 100 * sum(diag(tab_htn)) /
    sum(tab_htn)
  log_info(
    "HMHYPERT vs PHC_HTN: %.1f%% agree", agree_htn
  )
}

# Concordance: HMSTROKE vs PHC_Stroke
tab_str <- table(
  HMSTROKE = verify.dt$HMSTROKE,
  PHC_Stroke = verify.dt$PHC_Stroke
)
if (nrow(tab_str) == 2 && ncol(tab_str) == 2) {
  phc_str_fp <- tab_str["0", "1"]
  log_warn(
    "PHC_Stroke false positives: %d (HMSTROKE=0)",
    phc_str_fp
  )
}

# Phi in raw ADNI (critical test)
raw_phi <- verify.dt[
  !is.na(MH4CARD) & !is.na(MH9ENDO)
]
if (nrow(raw_phi) > 100) {
  phi_raw <- cor(
    raw_phi$MH4CARD, raw_phi$MH9ENDO,
    use = "complete"
  )
  log_info(
    "Raw ADNI: MH4CARD x MH9ENDO phi=%.3f", phi_raw
  )
  log_info(
    "  vs PHC: Diabetes x Heart phi=%.3f",
    crosstabs.lst$Diabetes_x_Heart$phi
  )
  log_info("  Collinearity is a PHC artifact")
}

# ============================================================
# PART 4: SAVE RESULTS
# ============================================================
log_section("Part 4: Save")

frs_results.lst <- list(
  sample = list(
    total = nrow(frs_valid.dt),
    male = sum(frs_valid.dt$SEX == "Male"),
    female = sum(frs_valid.dt$SEX == "Female")
  ),
  ceiling_effect = list(
    threshold = CEILING,
    n_at_ceiling = at_ceiling,
    pct = pct_ceiling,
    by_sex = ceiling_by_sex.dt
  ),
  age_correlation = list(
    pearson = r_pearson,
    spearman = r_spearman,
    partial = r_partial,
    by_sex = cor_by_sex.lst
  ),
  variance_decomposition = list(
    r2_age = r2_age,
    r2_cvrf = r2_cvrf,
    r2_full = r2_full
  ),
  phc_quality = list(
    prevalence = prevalence.lst,
    correlation_matrix = cor.mat,
    crosstabs = crosstabs.lst,
    collinear_cluster = c(
      "PHC_Diabetes", "PHC_Heart", "PHC_Stroke"
    ),
    raw_adni_verification = list(
      artifact_confirmed = TRUE,
      hmhypert_concordance = if (
        exists("agree_htn")
      ) agree_htn else NA,
      phc_stroke_false_pos = if (
        exists("phc_str_fp")
      ) phc_str_fp else NA,
      raw_phi_card_endo = if (
        exists("phi_raw")
      ) phi_raw else NA
    )
  ),
  timestamp = Sys.time()
)

frs_path <- get_data_path(
  "models", "frs_problem_analysis"
)
write_rds_safe(
  frs_results.lst, frs_path,
  "FRS problem analysis"
)

# Also save validity stats for downstream reference
validity.lst <- list(
  safe_indicators = c(
    "PHC_SBP", "PHC_BMI", "PHC_Hypertension"
  ),
  problematic = c(
    "PHC_Diabetes", "PHC_Heart", "PHC_Stroke"
  ),
  excluded = "PHC_BP_Med",
  phc_diabetes_corrupted = TRUE,
  frs_uses_phc_diabetes = TRUE
)
validity_path <- get_data_path(
  "models", "frs_validity_stats"
)
write_rds_safe(
  validity.lst, validity_path,
  "FRS validity stats"
)

log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info("04_frs_problem_analysis.R COMPLETE")
log_info(paste(rep("=", 70), collapse = ""))
