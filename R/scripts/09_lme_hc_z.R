#!/usr/bin/env Rscript

# =============================================================================
# 06_lme_hc_z.R - HC z-score LME Analysis (FRS + CVR_mimic)
# =============================================================================
# Longitudinal Mixed-Effects Models: CVR x HC_z Interaction
#
# This script runs LME models for BOTH CVR measures (FRS and CVR_mimic)
# in a single unified pipeline, enabling direct comparison.
#
# BRAIN MEASURE: HC z-score (UK Biobank normative model)
#   - Total hippocampal volume normalized using GAMLSS centiles
#   - Age/sex-adjusted deviation from healthy aging trajectory
#   - Negative z-scores indicate smaller-than-expected hippocampus
#   - Complements HVR_z by testing total volume vs ventricle-
#     adjusted ratio
#
# HYPOTHESES:
#   1. FRS moderates the HC-cognition trajectory (3-way interaction)
#   2. If FRS effects are age-driven, CVR_mimic (age-partialled)
#      should show NULL effects
#
# MODEL SPECIFICATION:
#   Cognition ~ (YRS + YRS^2) x CVR x HC_z + Age_bl + EDUC +
#               APOE4 + (YRS | PTID)
#
#   - YRS: Years from baseline (calendar time)
#   - YRS^2: Quadratic time (captures acceleration of decline)
#   - CVR: FRS or CVR_mimic
#   - HC_z: Hippocampal volume z-score
#   - Age_bl: Age at first observation in analysis cohort
#   - EDUC: Years of education
#   - APOE4: APOE e4 carrier status (0/1)
#   - (YRS | PTID): Random intercepts and slopes by subject
#
# QUADRATIC TIME JUSTIFICATION:
#   Script 06b_lme_yrs_nonlinearity.R tested linear vs quadratic:
#   All 12 models (2 sexes x 2 brain measures x 3 domains) showed
#   significantly better fit with quadratic time (LRT p < 0.05).
#
# CVR MEASURES:
#   FRS: Framingham Risk Score (10-year cardiovascular risk)
#     - Standard clinical composite; correlated with age
#   CVR_mimic: MIMIC model latent factor
#     - Age treated as cause, not indicator; residual CVR
#     - Orthogonal to age by construction
#
# ANALYSIS STRATEGY:
#   1. Sex-stratified models (PRIMARY) - NIH SABV policy
#   2. FRS vs CVR_mimic comparison
#
# WHY BOTH HVR_z AND HC_z?
#   - HVR (ratio) controls for head size and ventricular expansion
#   - HC (volume) is the traditional AD biomarker
#   - If effects replicate: finding is robust to measurement approach
#   - If effects differ: HVR may capture unique variance
#
# METHODOLOGICAL REFERENCES:
#   - SABV Policy: NIH ORWH (2016)
#     https://orwh.od.nih.gov/sex-gender/
#       nih-policy-sex-biological-variable
#   - Clayton & Collins (2014) Nature 509:282-283
#   - Ferretti et al. (2018) J Alzheimers Dis 64:S37-S50
#   - Jack et al. (2018) Alzheimers Dement 14:535-562
#     "NIA-AA Research Framework" - hippocampal atrophy as
#     neurodegeneration marker
#
# MULTIPLE COMPARISON CORRECTION:
#   FDR (Benjamini-Hochberg) applied WITHIN each CVR measure:
#   6 tests each (3 domains x 2 sexes). NOT across both measures.
#
# INPUTS:
#   - data/derivatives/lme/cohort_hvr.rds
#   - data/derivatives/lme/cvr_mimic_scores.rds
#
# OUTPUTS:
#   - models/results/lme/lme_hc_z_results.rds (both measures)
# =============================================================================

# -------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(lme4)
  library(lmerTest)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

log_script_start("09_lme_hc_z.R")
config <- load_config()
validate_config(config)
validate_packages(c("lme4", "lmerTest", "data.table"))
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "lme_hc", default = FALSE
)
output.path <- get_data_path("models", "lme_hc_z")
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE to rerun.")
  log_script_end("09_lme_hc_z.R", success = TRUE)
  quit(status = 0)
}

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------
BRAIN_VAR <- "HC_z"
BRAIN_LABEL <- "HC z-score (UKB normative model)"

CVR_MEASURES <- list(
  FRS = "FRS",
  CVR_mimic = "CVR_mimic"
)

COGNITIVE_DOMAINS <- get_parameter("cognitive_domains")

# LME settings from config
LME_ESTIMATOR <- get_script_setting(
  "lme", "estimator", default = "REML"
)
LME_OPTIMIZER <- get_script_setting(
  "lme", "optimizer", default = "bobyqa"
)
LME_MAXITER <- get_script_setting(
  "lme", "max_iter", default = 100000
)

results_dir.path <- get_data_path("models", "lme_results_dir")
ensure_directory(results_dir.path)

log_info("=== HC_z ANALYSIS: %s ===", BRAIN_LABEL)
log_info("CVR measures: %s", paste(names(CVR_MEASURES), collapse = ", "))
log_info("")

# =================================================================
# 1. Load and Validate Data
# =================================================================
log_section("Loading Data")

cohort.path <- get_data_path("derivatives", "lme_cohort_hvr")
check_files_exist(cohort.path)
cohort.dt <- read_rds_safe(cohort.path, "HVR cohort")

log_info(
  "Loaded cohort: %d subjects, %d observations",
  cohort.dt[, uniqueN(PTID)], nrow(cohort.dt)
)

# --- Load / verify CVR_mimic scores ---
# CVR_mimic is pre-merged by script 07 into cohort.
# If missing, attempt merge from scores file.
has_cvr_mimic <- FALSE

if ("CVR_mimic" %in% names(cohort.dt) &&
    sum(!is.na(cohort.dt$CVR_mimic)) > 0) {
  has_cvr_mimic <- TRUE
  n_with_cvr <- sum(!is.na(cohort.dt$CVR_mimic))
  log_info(
    "CVR_mimic already in cohort: %d obs",
    n_with_cvr
  )
  cvr_age_cor <- cor(
    cohort.dt$CVR_mimic, cohort.dt$AGE,
    use = "complete.obs"
  )
  log_info(
    "CVR_mimic-Age cor: r = %.3f (expect ~0)",
    cvr_age_cor
  )
} else {
  cvr_path <- file.path(
    dirname(cohort.path),
    "cvr_mimic_scores.rds"
  )
  if (!file.exists(cvr_path)) {
    log_warn(
      "CVR_mimic not in cohort and scores "
    )
    log_warn(
      "file not found: %s", cvr_path
    )
    stop("CVR_mimic scores unavailable")
  }
  cvr_scores.dt <- read_rds_safe(
    cvr_path, "CVR_mimic scores"
  )
  if (!"RID" %in% names(cvr_scores.dt)) {
    stop("CVR_mimic scores missing RID column")
  }
  cvr_merge.dt <- cvr_scores.dt[
    , .(RID, CVR_mimic)
  ]
  cohort.dt <- merge(
    cohort.dt, cvr_merge.dt,
    by = "RID", all.x = TRUE
  )
  n_with_cvr <- sum(
    !is.na(cohort.dt$CVR_mimic)
  )
  if (n_with_cvr == 0) {
    stop("No observations with CVR_mimic")
  }
  has_cvr_mimic <- TRUE
  log_info(
    "CVR_mimic merged: %d obs", n_with_cvr
  )
}

# Required columns for FRS analysis
required_cols.v <- c(
  "PTID", "EXAMDATE", "YRS_from_bl", "AGE", "EDUC",
  "SEX", "FRS", BRAIN_VAR, COGNITIVE_DOMAINS
)
validate_columns(cohort.dt, required_cols.v, "HVR cohort")
validate_not_empty(cohort.dt, "HVR cohort")

# Check APOE4 availability
HAS_APOE <- "APOE4" %in% names(cohort.dt) &&
  !all(is.na(cohort.dt$APOE4))
if (HAS_APOE) {
  apoe_n <- sum(
    cohort.dt[, .(apoe = APOE4[1]), by = PTID]$apoe,
    na.rm = TRUE
  )
  apoe_pct <- 100 * apoe_n / cohort.dt[, uniqueN(PTID)]
  log_info(
    "APOE4 available: %d carriers (%.1f%%)",
    apoe_n, apoe_pct
  )
} else {
  log_warn("APOE4 not available - excluding from models")
}

# =================================================================
# 2. Compute Age at Baseline (Age_bl)
# =================================================================
log_section("Computing Baseline Age")

setorder(cohort.dt, PTID, EXAMDATE)
if (!"Age_bl" %in% names(cohort.dt)) {
  baseline_age.dt <- cohort.dt[
    , .(Age_bl = AGE[1]), by = PTID
  ]
  cohort.dt <- merge(
    cohort.dt, baseline_age.dt,
    by = "PTID", all.x = TRUE
  )
}

log_info("Age_bl = age at first observation in analysis cohort")
log_info(
  "  Range: %.1f - %.1f years",
  min(cohort.dt$Age_bl), max(cohort.dt$Age_bl)
)
log_info(
  "  Mean (SD): %.1f (%.1f)",
  mean(cohort.dt$Age_bl), sd(cohort.dt$Age_bl)
)

# Note: Quadratic time modeled via I(YRS_from_bl^2)
# Reference: R-sig-ME discussion on lmer quadratic terms
log_info(
  "YRS range: %.1f - %.1f years",
  min(cohort.dt$YRS_from_bl), max(cohort.dt$YRS_from_bl)
)

# =================================================================
# 3. Define Model Formulas and Helper Functions
# =================================================================
# Model: Cognition ~ YRS x CVR x HC_z + I(YRS^2) x CVR x HC_z
#                   + covariates + (YRS | PTID)
#
# The 3-way interaction YRS x CVR x HC_z tests:
#   Does the CVR x HC effect on cognition CHANGE over time?
#
# The I(YRS^2) x CVR x HC_z terms test:
#   Does the 3-way interaction ACCELERATE over time?
#
# Using I() for squared term follows R-sig-ME recommendations.
# -------------------------------------------------------------------

COVARIATES_STRAT <- if (HAS_APOE) {
  "Age_bl + EDUC + APOE4"
} else {
  "Age_bl + EDUC"
}

fit_lme_growth.fn <- function(formula, data, description) {
  log_info("  Fitting: %s", description)

  result.lst <- tryCatch({
    warnings.lst <- list()
    model.fit <- withCallingHandlers(
      lmer(
        formula, data = data,
        REML = (LME_ESTIMATOR == "REML"),
        control = lmerControl(
          optimizer = LME_OPTIMIZER,
          optCtrl = list(maxfun = LME_MAXITER)
        )
      ),
      warning = function(w) {
        warnings.lst <<- c(
          warnings.lst, conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    )

    is_singular <- isSingular(model.fit)
    converged <- !is_singular

    if (is_singular) {
      log_warn(
        "    Singular fit (random effect variance near zero)"
      )
    }
    if (length(warnings.lst) > 0) {
      log_warn(
        "    Warnings: %s",
        paste(warnings.lst, collapse = "; ")
      )
    }

    coefs.mat <- summary(model.fit)$coefficients
    list(
      fit = model.fit,
      coefficients = coefs.mat,
      converged = converged,
      is_singular = is_singular,
      warnings = warnings.lst,
      n_obs = nrow(data),
      n_subjects = data[, uniqueN(PTID)]
    )
  }, error = function(e) {
    log_warn("    Model failed: %s", e$message)
    list(
      fit = NULL, converged = FALSE,
      error = e$message
    )
  })

  result.lst
}

extract_interaction.fn <- function(
  model.res, cvr_var, brain_var
) {
  if (is.null(model.res$fit) || !model.res$converged) {
    return(list(
      term = NA, beta = NA, se = NA,
      t_value = NA, p_value = NA
    ))
  }

  coefs.mat <- model.res$coefficients

  # Match any ordering of the 3-way interaction
  # YRS_from_bl : cvr_var : brain_var (all permutations)
  parts.v <- c("YRS_from_bl", cvr_var, brain_var)
  perms.v <- apply(
    expand.grid(
      parts.v, parts.v, parts.v
    ), 1, function(x) {
      if (length(unique(x)) == 3) {
        paste(x, collapse = ":")
      } else {
        NA
      }
    }
  )
  perms.v <- perms.v[!is.na(perms.v)]
  pattern <- paste(perms.v, collapse = "|")

  idx <- grep(pattern, rownames(coefs.mat))

  if (length(idx) == 0) {
    return(list(
      term = NA, beta = NA, se = NA,
      t_value = NA, p_value = NA
    ))
  }

  idx <- idx[1]
  list(
    term = rownames(coefs.mat)[idx],
    beta = coefs.mat[idx, "Estimate"],
    se = coefs.mat[idx, "Std. Error"],
    t_value = coefs.mat[idx, "t value"],
    p_value = coefs.mat[idx, "Pr(>|t|)"]
  )
}

# =================================================================
# 4. PRIMARY ANALYSIS: Sex-Stratified Models (SABV)
#    Loop over CVR measures
# =================================================================
log_section("PRIMARY: Sex-Stratified Analysis (SABV)")
log_info("")
log_info(
  "Rationale: NIH SABV policy + known AD sexual dimorphism"
)
log_info(
  "Reference: Clayton & Collins (2014) Nature 509:282-283"
)
log_info("")

all_stratified.lst <- list()

for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]

  # Determine analysis data for this CVR measure
  if (cvr_name == "CVR_mimic") {
    if (!has_cvr_mimic) {
      log_warn("Skipping CVR_mimic (not available)")
      next
    }
    analysis.dt <- cohort.dt[!is.na(CVR_mimic)]
  } else {
    analysis.dt <- cohort.dt
  }

  log_info(
    "========== CVR measure: %s ==========", cvr_name
  )
  log_info(
    "Model: Cognition ~ YRS * %s * %s + I(YRS^2) * %s * %s + %s + (YRS | PTID)",
    CVR_VAR, BRAIN_VAR, CVR_VAR, BRAIN_VAR,
    COVARIATES_STRAT
  )
  log_info("")

  stratified.lst <- list()

  for (sex in c("Male", "Female")) {
    log_info("--- %s ---", toupper(sex))

    sex_cohort.dt <- analysis.dt[SEX == sex]
    n_subj <- sex_cohort.dt[, uniqueN(PTID)]
    n_obs <- nrow(sex_cohort.dt)
    log_info(
      "N = %d subjects, %d observations",
      n_subj, n_obs
    )

    stratified.lst[[sex]] <- list()

    for (domain in COGNITIVE_DOMAINS) {
      formula_str <- sprintf(
        paste0(
          "%s ~ YRS_from_bl * %s * %s + ",
          "I(YRS_from_bl^2) * %s * %s + ",
          "%s + (YRS_from_bl | PTID)"
        ),
        domain, CVR_VAR, BRAIN_VAR,
        CVR_VAR, BRAIN_VAR, COVARIATES_STRAT
      )

      model.res <- fit_lme_growth.fn(
        as.formula(formula_str),
        sex_cohort.dt,
        sprintf(
          "%s: YRS x %s x %s + I(YRS^2)",
          domain, CVR_VAR, BRAIN_VAR
        )
      )

      int.lst <- extract_interaction.fn(
        model.res, CVR_VAR, BRAIN_VAR
      )

      stratified.lst[[sex]][[domain]] <- list(
        model = model.res,
        interaction = int.lst,
        formula = formula_str,
        n_subjects = n_subj,
        n_obs = n_obs
      )

      if (!is.na(int.lst$p_value)) {
        sig <- ifelse(int.lst$p_value < 0.05, " *", "")
        log_info(
          "    %s: b = %.4f (SE = %.4f), p = %.4f%s",
          domain, int.lst$beta, int.lst$se,
          int.lst$p_value, sig
        )
      } else {
        log_info(
          "    %s: Model did not converge", domain
        )
      }
    }
    log_info("")
  }

  all_stratified.lst[[cvr_name]] <- stratified.lst
}

# =================================================================
# 5. Create Summary Table
# =================================================================
log_section("Summary Table")

summary_rows.lst <- list()

# --- Stratified results (primary) for each CVR measure ---
for (cvr_name in names(all_stratified.lst)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]

  for (sex in c("Male", "Female")) {
    for (domain in COGNITIVE_DOMAINS) {
      res <- all_stratified.lst[[cvr_name]][[sex]][[domain]]
      int <- res$interaction

      summary_rows.lst[[length(summary_rows.lst) + 1]] <-
        data.table(
          Brain_Measure = BRAIN_VAR,
          CVR_Measure = cvr_name,
          Analysis = "Primary_Stratified",
          Sex = sex,
          Domain = domain,
          N_subjects = res$n_subjects,
          N_obs = res$n_obs,
          Beta = int$beta,
          SE = int$se,
          t_value = int$t_value,
          p_value = int$p_value
        )
    }
  }
}

summary.dt <- rbindlist(summary_rows.lst)

# -------------------------------------------------------------------
# FDR correction WITHIN each CVR measure (6 tests each)
# 3 domains x 2 sexes = 6 tests per CVR measure
# -------------------------------------------------------------------
for (cvr_name in names(CVR_MEASURES)) {
  mask <- summary.dt$Analysis == "Primary_Stratified" &
    summary.dt$CVR_Measure == cvr_name
  if (any(mask)) {
    summary.dt[
      mask,
      p_fdr := p.adjust(p_value, method = "BH")
    ]
  }
}

summary.dt[, Significant := p_value < 0.05]
summary.dt[, Significant_FDR := p_fdr < 0.05]

# Log summary per CVR measure
for (cvr_name in names(CVR_MEASURES)) {
  strat_mask <- summary.dt$Analysis == "Primary_Stratified" &
    summary.dt$CVR_Measure == cvr_name
  n_tests <- sum(strat_mask)
  n_sig <- sum(
    summary.dt[strat_mask, Significant], na.rm = TRUE
  )
  n_fdr <- sum(
    summary.dt[strat_mask, Significant_FDR], na.rm = TRUE
  )
  log_info(
    "%s: %d tests, %d sig (uncorrected), %d sig (FDR)",
    cvr_name, n_tests, n_sig, n_fdr
  )
}

# =================================================================
# 6. FRS vs CVR_mimic Comparison
# =================================================================
log_section("FRS vs CVR_mimic Comparison (HC_z)")

if (has_cvr_mimic &&
    "CVR_mimic" %in% summary.dt$CVR_Measure) {
  log_info("")
  log_info(paste(rep("=", 55), collapse = ""))
  log_info("COMPARISON: FRS vs CVR_mimic effects on HC_z")
  log_info(paste(rep("=", 55), collapse = ""))

  for (sex in c("Male", "Female")) {
    log_info("")
    log_info("--- %s ---", sex)

    for (domain in COGNITIVE_DOMAINS) {
      frs_row <- summary.dt[
        Analysis == "Primary_Stratified" &
          CVR_Measure == "FRS" &
          Sex == sex &
          Domain == domain
      ]
      cvr_row <- summary.dt[
        Analysis == "Primary_Stratified" &
          CVR_Measure == "CVR_mimic" &
          Sex == sex &
          Domain == domain
      ]

      if (nrow(frs_row) > 0 && nrow(cvr_row) > 0) {
        frs_sig <- ifelse(
          frs_row$p_value < 0.05, "*", ""
        )
        cvr_sig <- ifelse(
          cvr_row$p_value < 0.05, "*", ""
        )

        log_info("  %s:", domain)
        log_info(
          "    FRS:       b = %7.4f, p = %.4f%s",
          frs_row$Beta, frs_row$p_value, frs_sig
        )
        log_info(
          "    CVR_mimic: b = %7.4f, p = %.4f%s",
          cvr_row$Beta, cvr_row$p_value, cvr_sig
        )
      }
    }
  }

  # Interpretation guidance
  log_info("")
  log_info("Interpretation:")
  log_info(
    "  If FRS sig but CVR_mimic null: age-weighting drives effect"
  )
  log_info(
    "  If both sig: genuine CVR effect beyond age confounding"
  )
  log_info(
    "  If neither sig: no CVR x HC_z moderation of decline"
  )
} else {
  log_warn(
    "CVR_mimic not available - skipping comparison"
  )
}

# =================================================================
# 7. Key Findings Summary
# =================================================================
log_section("KEY FINDINGS")
log_info("")
log_info("Brain measure: %s", BRAIN_LABEL)
log_info("Primary analysis: Sex-stratified (SABV)")
log_info("")

for (cvr_name in names(all_stratified.lst)) {
  log_info("--- CVR: %s ---", cvr_name)
  for (sex in c("Male", "Female")) {
    log_info("%s:", sex)
    for (domain in COGNITIVE_DOMAINS) {
      res <- all_stratified.lst[[cvr_name]][[sex]][[domain]]
      int <- res$interaction
      if (!is.na(int$p_value)) {
        sig <- ifelse(int$p_value < 0.05, "***", "")
        log_info(
          "  %s: b = %.4f, p = %.4f %s",
          domain, int$beta, int$p_value, sig
        )
      }
    }
  }
  log_info("")
}

# =================================================================
# 8. Save Results
# =================================================================
log_section("Saving Results")

cvr_validation.lst <- list()
if (has_cvr_mimic) {
  cvr_validation.lst <- list(
    age_correlation = cvr_age_cor,
    expected_null = abs(cvr_age_cor) < 0.1
  )
}

output.lst <- list(
  # Metadata
  brain_measure = BRAIN_VAR,
  brain_label = BRAIN_LABEL,
  cvr_measures = CVR_MEASURES,
  analysis_type = "primary",

  # Results by CVR measure
  stratified = all_stratified.lst,
  summary_table = summary.dt,

  # CVR_mimic validation
  cvr_validation = cvr_validation.lst,

  # Methodology documentation
  methodology = list(
    primary_analysis = paste(
      "Sex-stratified (NIH SABV policy)"
    ),
    model_stratified = sprintf(
      "Cognition ~ YRS * CVR * %s + %s + (YRS | PTID)",
      BRAIN_VAR, COVARIATES_STRAT
    ),
    fdr_correction = paste(
      "BH within each CVR measure",
      "(6 tests: 3 domains x 2 sexes)"
    ),
    has_apoe = HAS_APOE,
    cvr_rationale = paste(
      "CVR_mimic from MIMIC model partials out age.",
      "If FRS effects are age-driven,",
      "CVR_mimic effects should be null."
    ),
    references = list(
      sabv_policy = "NIH ORWH SABV Policy (2016)",
      sabv_rationale = paste(
        "Clayton & Collins (2014)",
        "Nature 509:282-283"
      ),
      ad_sex_diff = paste(
        "Ferretti et al. (2018)",
        "J Alzheimers Dis 64:S37-S50"
      ),
      nia_aa_framework = paste(
        "Jack et al. (2018)",
        "Alzheimers Dement 14:535-562"
      )
    )
  ),

  # Config
  config = list(
    estimator = LME_ESTIMATOR,
    optimizer = LME_OPTIMIZER,
    timestamp = Sys.time()
  )
)

output.path <- file.path(
  results_dir.path, "lme_hc_z_results.rds"
)
write_rds_safe(
  output.lst, output.path, "LME HC_z results (both CVR)"
)
log_info("Saved: %s", output.path)

log_script_end("06_lme_hc_z.R", success = TRUE)
