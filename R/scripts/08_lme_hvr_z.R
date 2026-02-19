#!/usr/bin/env Rscript

# =============================================================================
# 05_lme_hvr_z.R - LME Models: CVR x HVR_z Interaction (FRS + CVR_mimic)
# =============================================================================
# Longitudinal Mixed-Effects Models: CVR x HVR_z Interaction
#
# BRAIN MEASURE: HVR z-score (UK Biobank normative model)
#   - Hippocampal-to-Ventricle Ratio normalized using GAMLSS centiles
#   - Age/sex-adjusted deviation from healthy aging trajectory
#   - Negative z-scores indicate accelerated hippocampal atrophy
#
# CVR MEASURES:
#   1. FRS  - Framingham Risk Score (10-year cardiovascular risk)
#   2. CVR_mimic - Age-adjusted CVR from MIMIC latent factor model
#
# HYPOTHESIS:
#   Cardiovascular risk (FRS) moderates the relationship between
#   hippocampal integrity (HVR) and cognitive trajectories in
#   amyloid-positive adults. Higher FRS weakens the protective
#   association between preserved HVR and maintained cognition.
#
#   If FRS effects are driven by age-weighting, CVR_mimic x HVR
#   interaction should NOT be significant. The key comparison is
#   FRS vs CVR_mimic effects.
#
# MODEL SPECIFICATION:
#   Cognition ~ (YRS + YRS^2) x CVR x HVR_z +
#               Age_bl + EDUC + APOE4 + (YRS | PTID)
#
#   - YRS: Years from baseline (centered for interpretation)
#   - YRS^2: Quadratic time term (captures acceleration)
#   - CVR: FRS or CVR_mimic (looped)
#   - HVR_z: Hippocampal-to-Ventricle Ratio z-score
#   - Age_bl: Age at first observation in analysis cohort
#   - EDUC: Years of education
#   - APOE4: APOE e4 carrier status (0/1)
#   - (YRS | PTID): Random intercepts and slopes by subject
#
# QUADRATIC TIME JUSTIFICATION:
#   Script 06b_lme_yrs_nonlinearity.R tested linear vs quadratic
#   time: All 12 models (2 sexes x 2 brain measures x 3 domains)
#   showed significantly better fit with quadratic time
#   (LRT p < 0.05). This captures acceleration of cognitive
#   decline over follow-up.
#
# MIMIC MODEL RATIONALE:
#   The Multiple Indicators Multiple Causes (MIMIC) model treats
#   age as a CAUSE of CVR, not an indicator. This extracts CVR
#   factor scores representing residual CVR after accounting for
#   age, yielding a measure orthogonal to age by construction.
#
# ANALYSIS STRATEGY:
#   1. Sex-stratified models (PRIMARY) - NIH SABV policy
#   2. FRS vs CVR_mimic comparison
#
# WHY SEX-STRATIFIED (SABV)?
#   NIH policy requires consideration of sex as a biological
#   variable (SABV). AD shows marked sexual dimorphism in risk,
#   progression, and biomarkers. When sex differences are
#   hypothesized, stratified analysis is preferred over
#   interaction terms for interpretability and power.
#
# METHODOLOGICAL REFERENCES:
#   - SABV Policy: NIH ORWH (2016)
#   - Clayton & Collins (2014) Nature 509:282-283
#   - Mazure & Swendsen (2016) Biol Psychiatry 80:e11-e14
#   - Ferretti et al. (2018) J Alzheimers Dis 64:S37-S50
#   - Buckley et al. (2019) Ann Neurol 85:430-442
#   - MIMIC: Joreskog & Goldberger (1975) J Econom 3:313-329
#
# MULTIPLE COMPARISON CORRECTION:
#   FDR (Benjamini-Hochberg) applied WITHIN each CVR measure:
#   6 tests per measure (3 domains x 2 sexes).
#   FRS and CVR_mimic are corrected separately because they
#   address distinct hypotheses.
#
# INPUTS:
#   - data/derivatives/lme/cohort_hvr.rds
#
# OUTPUTS:
#   - models/results/lme/lme_hvr_z_results.rds
# =============================================================================

# -----------------------------------------------------------
# Setup
# -----------------------------------------------------------
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

log_script_start("08_lme_hvr_z.R")
config <- load_config()
validate_config(config)
validate_packages(c("lme4", "lmerTest", "data.table"))
set_seed()

# -----------------------------------------------------------
# Configuration
# -----------------------------------------------------------
BRAIN_VAR <- "HVR_z"
BRAIN_LABEL <- "HVR z-score (UKB normative model)"
COGNITIVE_DOMAINS <- get_parameter("cognitive_domains")

CVR_MEASURES <- list(
  FRS = "FRS",
  CVR_mimic = "CVR_mimic"
)

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

log_info("=== LME ANALYSIS: %s ===", BRAIN_LABEL)
log_info(
  "CVR measures: %s",
  paste(names(CVR_MEASURES), collapse = ", ")
)
log_info("")

# ============================================================
# 1. Load and Validate Data
# ============================================================
log_section("Loading Data")

cohort.path <- get_data_path("derivatives", "lme_cohort_hvr")
check_files_exist(cohort.path)
cohort.dt <- read_rds_safe(cohort.path, "HVR cohort")

# Required columns (common to both CVR measures)
required_cols.v <- c(
  "PTID", "EXAMDATE", "YRS_from_bl", "AGE", "EDUC",
  "SEX", "FRS", "CVR_mimic", BRAIN_VAR,
  COGNITIVE_DOMAINS
)
validate_columns(cohort.dt, required_cols.v, "HVR cohort")
validate_not_empty(cohort.dt, "HVR cohort")

log_info(
  "Loaded: %d subjects, %d observations",
  cohort.dt[, uniqueN(PTID)], nrow(cohort.dt)
)

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

# -----------------------------------------------------------
# Validate CVR_mimic (one-time check)
# -----------------------------------------------------------
# CVR_mimic is pre-merged and pre-standardized (z-scored on
# the common analysis sample) in 03_prepare_lme_data.R S7c.
# Both FRS and CVR_mimic share the same reference population,
# so their LME beta coefficients are directly comparable.
n_cvr_na <- sum(is.na(cohort.dt[["CVR_mimic"]]))
if (n_cvr_na > 0) {
  log_warn(
    "%d observations with NA CVR_mimic - filtering",
    n_cvr_na
  )
  cohort.dt <- cohort.dt[!is.na(CVR_mimic)]
}

# Age-independence is a design feature of the MIMIC model:
# age is modelled as a CAUSE of CVR, so extracted factor
# scores represent CVR BEYOND age expectations.
cvr_age_cor <- cor(
  cohort.dt[["CVR_mimic"]], cohort.dt$AGE,
  use = "complete.obs"
)
log_info(
  "CVR_mimic-Age correlation: r = %.3f (expect ~0)",
  cvr_age_cor
)
if (abs(cvr_age_cor) > 0.1) {
  log_warn(
    "CVR_mimic shows unexpected age correlation (r=%.3f)",
    cvr_age_cor
  )
}

log_info(
  "Sample after CVR_mimic check: %d subjects, %d obs",
  cohort.dt[, uniqueN(PTID)], nrow(cohort.dt)
)

# ============================================================
# 2. Compute Age at Baseline (Age_bl)
# ============================================================
# Age_bl = age at FIRST observation in OUR analysis cohort.
# This is NOT necessarily ADNI study entry - it is entry into
# this specific analysis after all inclusion criteria applied.
# -----------------------------------------------------------
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

log_info("Age_bl = age at first obs in analysis cohort")
log_info(
  "  Range: %.1f - %.1f years",
  min(cohort.dt$Age_bl), max(cohort.dt$Age_bl)
)
log_info(
  "  Mean (SD): %.1f (%.1f)",
  mean(cohort.dt$Age_bl), sd(cohort.dt$Age_bl)
)

# Quadratic time modeled via I(YRS_from_bl^2) in formulas
# Reference: R-sig-ME discussion on lmer quadratic terms
log_info(
  "YRS range: %.1f - %.1f years",
  min(cohort.dt$YRS_from_bl),
  max(cohort.dt$YRS_from_bl)
)

# Report both CVR measures (pre-standardized from script 03)
for (cvr_name in names(CVR_MEASURES)) {
  cvr_col <- CVR_MEASURES[[cvr_name]]
  log_info(
    "%s (pre-standardized): M=%.4f, SD=%.4f",
    cvr_col,
    mean(cohort.dt[[cvr_col]], na.rm = TRUE),
    sd(cohort.dt[[cvr_col]], na.rm = TRUE)
  )
}

# ============================================================
# 3. Define Model Formulas and Helper Functions
# ============================================================
# Model:
#   Cognition ~ YRS x CVR x HVR_z +
#               I(YRS^2) x CVR x HVR_z +
#               covariates + (YRS | PTID)
#
# The 3-way interaction YRS x CVR x HVR_z tests:
#   Does the CVR x HVR effect on cognition CHANGE over time?
#   i.e., Does cardiovascular risk increasingly moderate the
#   brain-cognition relationship as disease progresses?
#
# The I(YRS^2) x CVR x HVR_z terms test:
#   Does the 3-way interaction ACCELERATE over time?
#
# Using I() for squared term follows R-sig-ME recommendations.
# Random slopes (YRS | PTID) account for individual
# differences in cognitive trajectories.
# -----------------------------------------------------------

# Build covariate strings
COVARIATES_STRAT <- if (HAS_APOE) {
  "Age_bl + EDUC + APOE4"
} else {
  "Age_bl + EDUC"
}
COVARIATES_FULL <- if (HAS_APOE) {
  "Age_bl + EDUC + SEX + APOE4"
} else {
  "Age_bl + EDUC + SEX"
}

fit_lme_growth.fn <- function(formula, data, description) {
  # Fit LME model with convergence handling
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
        "    Singular fit (random effect var near zero)"
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
      fit = NULL, converged = FALSE, error = e$message
    )
  })

  result.lst
}

extract_interaction.fn <- function(model.res, cvr_var) {
  # Extract 3-way interaction: YRS:cvr_var:BRAIN_VAR
  # cvr_var: character name of CVR column (e.g. "FRS")
  if (is.null(model.res$fit) || !model.res$converged) {
    return(list(
      term = NA, beta = NA, se = NA,
      t_value = NA, p_value = NA
    ))
  }

  coefs.mat <- model.res$coefficients

  # Match any ordering of the 3-way interaction
  pattern <- paste0(
    "YRS_from_bl:", cvr_var, ":", BRAIN_VAR, "|",
    "YRS_from_bl:", BRAIN_VAR, ":", cvr_var, "|",
    cvr_var, ":YRS_from_bl:", BRAIN_VAR, "|",
    cvr_var, ":", BRAIN_VAR, ":YRS_from_bl|",
    BRAIN_VAR, ":YRS_from_bl:", cvr_var, "|",
    BRAIN_VAR, ":", cvr_var, ":YRS_from_bl"
  )

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

# ============================================================
# 4. PRIMARY ANALYSIS: Sex-Stratified Models (SABV)
# ============================================================
# Following NIH SABV policy: analyze males and females
# separately when sex differences are hypothesized.
# AD shows marked sexual dimorphism.
#
# Model per sex per CVR measure:
#   Cognition ~ YRS x CVR x HVR_z +
#               I(YRS^2) x CVR x HVR_z +
#               covariates + (YRS | PTID)
#
# References:
#   Clayton & Collins (2014) Nature 509:282-283
#   Buckley et al. (2019) Ann Neurol 85:430-442
# -----------------------------------------------------------
log_section("PRIMARY: Sex-Stratified Analysis (SABV)")
log_info("")
log_info("Rationale: NIH SABV + known AD sexual dimorphism")
log_info("")

# Store all stratified results keyed by CVR measure
all_stratified.lst <- list()

for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]
  log_info(
    "====== CVR measure: %s ======", CVR_VAR
  )
  log_info(
    "Model: Cog ~ YRS x %s x %s + I(YRS^2) x %s x %s + %s + (YRS | PTID)",
    CVR_VAR, BRAIN_VAR, CVR_VAR, BRAIN_VAR,
    COVARIATES_STRAT
  )
  log_info("")

  stratified_results.lst <- list()

  for (sex in c("Male", "Female")) {
    log_info("--- %s ---", toupper(sex))

    sex_cohort.dt <- cohort.dt[SEX == sex]
    n_subj <- sex_cohort.dt[, uniqueN(PTID)]
    n_obs <- nrow(sex_cohort.dt)
    log_info(
      "N = %d subjects, %d observations",
      n_subj, n_obs
    )

    stratified_results.lst[[sex]] <- list()

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
        model.res, CVR_VAR
      )

      stratified_results.lst[[sex]][[domain]] <- list(
        model = model.res,
        interaction = int.lst,
        formula = formula_str,
        n_subjects = n_subj,
        n_obs = n_obs
      )

      # Report result
      if (!is.na(int.lst$p_value)) {
        sig <- ifelse(int.lst$p_value < 0.05, " *", "")
        log_info(
          "    %s: B = %.4f (SE = %.4f), p = %.4f%s",
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

  all_stratified.lst[[cvr_name]] <- stratified_results.lst
}

# ============================================================
# 5. Create Summary Table
# ============================================================
log_section("Summary Table")

summary_rows.lst <- list()

# Stratified results (primary) - both CVR measures
for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]
  strat.lst <- all_stratified.lst[[cvr_name]]

  for (sex in c("Male", "Female")) {
    for (domain in COGNITIVE_DOMAINS) {
      res <- strat.lst[[sex]][[domain]]
      int <- res$interaction

      summary_rows.lst[[length(summary_rows.lst) + 1]] <-
        data.table(
          Brain_Measure = BRAIN_VAR,
          CVR_Measure = CVR_VAR,
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

# FDR correction WITHIN each CVR measure (6 tests each:
# 3 domains x 2 sexes). Corrected separately because FRS
# and CVR_mimic address distinct hypotheses.
for (cvr_name in names(CVR_MEASURES)) {
  cvr_col <- CVR_MEASURES[[cvr_name]]
  mask <- summary.dt$Analysis == "Primary_Stratified" &
    summary.dt$CVR_Measure == cvr_col
  summary.dt[mask,
    p_fdr := p.adjust(p_value, method = "BH")
  ]
}

summary.dt[, Significant := p_value < 0.05]
summary.dt[, Significant_FDR := p_fdr < 0.05]

# Report per CVR measure
for (cvr_name in names(CVR_MEASURES)) {
  cvr_col <- CVR_MEASURES[[cvr_name]]
  sub.dt <- summary.dt[
    Analysis == "Primary_Stratified" &
      CVR_Measure == cvr_col
  ]
  log_info(
    "%s primary: %d tests, %d sig (uncorr), %d sig (FDR)",
    cvr_col, nrow(sub.dt),
    sum(sub.dt$Significant, na.rm = TRUE),
    sum(sub.dt$Significant_FDR, na.rm = TRUE)
  )
}

# ============================================================
# 6. FRS vs CVR_mimic Comparison
# ============================================================
# Because both FRS and CVR_mimic are z-scored on the same
# common sample (script 03, section 7c), the beta
# coefficients are on the same scale and directly comparable.
# A 1-SD increase in FRS vs CVR_mimic can be meaningfully
# compared. If CVR_mimic betas are null while FRS betas are
# significant, this suggests FRS effects are driven by the
# age-weighting inherent in the Framingham equation rather
# than pure cardiovascular risk.
# -----------------------------------------------------------
log_section("FRS vs CVR_mimic Comparison")
log_info("")
log_info("COMPARISON: FRS vs CVR_mimic (%s)", BRAIN_VAR)
log_info(paste(rep("=", 50), collapse = ""))

for (sex in c("Male", "Female")) {
  log_info("")
  log_info("--- %s ---", sex)

  for (domain in COGNITIVE_DOMAINS) {
    frs_int <- all_stratified.lst[["FRS"]][[
      sex
    ]][[domain]]$interaction
    cvr_int <- all_stratified.lst[["CVR_mimic"]][[
      sex
    ]][[domain]]$interaction

    frs_sig <- ifelse(
      !is.na(frs_int$p_value) & frs_int$p_value < 0.05,
      "*", ""
    )
    cvr_sig <- ifelse(
      !is.na(cvr_int$p_value) & cvr_int$p_value < 0.05,
      "*", ""
    )

    log_info("  %s:", domain)
    if (!is.na(frs_int$p_value)) {
      log_info(
        "    FRS:       B = %7.4f, p = %.4f%s",
        frs_int$beta, frs_int$p_value, frs_sig
      )
    } else {
      log_info("    FRS:       did not converge")
    }
    if (!is.na(cvr_int$p_value)) {
      log_info(
        "    CVR_mimic: B = %7.4f, p = %.4f%s",
        cvr_int$beta, cvr_int$p_value, cvr_sig
      )
    } else {
      log_info("    CVR_mimic: did not converge")
    }
  }
}

# ============================================================
# 7. Key Findings Summary
# ============================================================
log_section("KEY FINDINGS")
log_info("")
log_info("Brain measure: %s", BRAIN_LABEL)
log_info("Primary analysis: Sex-stratified (SABV)")
log_info("")

for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]
  log_info("CVR: %s", CVR_VAR)

  for (sex in c("Male", "Female")) {
    log_info("  %s:", sex)
    for (domain in COGNITIVE_DOMAINS) {
      res <- all_stratified.lst[[cvr_name]][[
        sex
      ]][[domain]]
      int <- res$interaction
      if (!is.na(int$p_value)) {
        sig <- ifelse(int$p_value < 0.05, "***", "")
        log_info(
          "    %s: B = %.4f, p = %.4f %s",
          domain, int$beta, int$p_value, sig
        )
      }
    }
  }
  log_info("")
}

# ============================================================
# 7b. Bootstrap Comparison: FRS vs CVR_mimic 3-way
# ============================================================
# Clustered bootstrap (resampling subjects) to
# formally compare the 3-way interaction beta
# between FRS and CVR_mimic models on the same
# data, yielding a paired delta_beta with BC CI.
# -----------------------------------------------------------
log_section(
  "Bootstrap: FRS vs CVR_mimic 3-way"
)

N_BOOT_LME <- 1000
BOOT_CI_LEVEL <- 0.95
comparison.lst <- list()

boot_chk_dir <- file.path(
  results_dir.path, "lme_boot_checkpoints"
)
ensure_directory(boot_chk_dir)

for (sex in c("Male", "Female")) {
  sex_cohort.dt <- cohort.dt[SEX == sex]
  ptids.v <- unique(sex_cohort.dt$PTID)
  n_subj <- length(ptids.v)

  for (domain in COGNITIVE_DOMAINS) {
    cmp_key <- paste(sex, domain, sep = "_")

    # Check checkpoint
    boot_chk <- file.path(
      boot_chk_dir,
      paste0(cmp_key, ".rds")
    )
    if (file.exists(boot_chk)) {
      cached <- read_rds_safe(
        boot_chk,
        paste("boot chk", cmp_key)
      )
      if (!is.null(cached) &&
          isTRUE(cached$n_boot == N_BOOT_LME)
      ) {
        log_info(
          "%s %s: boot checkpoint found",
          sex, domain
        )
        comparison.lst[[cmp_key]] <- cached
        next
      }
    }

    log_info(
      "%s %s: bootstrapping %d reps...",
      sex, domain, N_BOOT_LME
    )

    frs_form <- as.formula(sprintf(
      paste0(
        "%s ~ YRS_from_bl * FRS * %s + ",
        "I(YRS_from_bl^2) * FRS * %s + ",
        "%s + (YRS_from_bl | .BOOT_ID)"
      ),
      domain, BRAIN_VAR,
      BRAIN_VAR, COVARIATES_STRAT
    ))
    cvr_form <- as.formula(sprintf(
      paste0(
        "%s ~ YRS_from_bl * CVR_mimic",
        " * %s + ",
        "I(YRS_from_bl^2) * CVR_mimic",
        " * %s + ",
        "%s + (YRS_from_bl | .BOOT_ID)"
      ),
      domain, BRAIN_VAR,
      BRAIN_VAR, COVARIATES_STRAT
    ))

    # Pre-split data by subject for speed
    split_idx.lst <- split(
      seq_len(nrow(sex_cohort.dt)),
      sex_cohort.dt$PTID
    )

    delta_beta.v <- numeric(N_BOOT_LME)
    set.seed(42 + which(
      c("Male", "Female") == sex
    ) * 10 + which(
      COGNITIVE_DOMAINS == domain
    ))

    for (b in seq_len(N_BOOT_LME)) {
      boot_ids <- sample(
        ptids.v, n_subj, replace = TRUE
      )
      # Expand + assign unique boot ID
      boot_rows <- unlist(
        split_idx.lst[boot_ids]
      )
      boot.dt <- sex_cohort.dt[boot_rows]
      boot.dt[, .BOOT_ID := rep(
        seq_along(boot_ids),
        lengths(split_idx.lst[boot_ids])
      )]

      frs_ok <- tryCatch({
        frs_fit <- suppressWarnings(lmer(
          frs_form, data = boot.dt,
          REML = (LME_ESTIMATOR == "REML"),
          control = lmerControl(
            optimizer = LME_OPTIMIZER,
            optCtrl = list(
              maxfun = LME_MAXITER
            )
          )
        ))
        frs_coefs <- fixef(frs_fit)
        TRUE
      }, error = function(e) FALSE)

      cvr_ok <- tryCatch({
        cvr_fit <- suppressWarnings(lmer(
          cvr_form, data = boot.dt,
          REML = (LME_ESTIMATOR == "REML"),
          control = lmerControl(
            optimizer = LME_OPTIMIZER,
            optCtrl = list(
              maxfun = LME_MAXITER
            )
          )
        ))
        cvr_coefs <- fixef(cvr_fit)
        TRUE
      }, error = function(e) FALSE)

      if (isTRUE(frs_ok) &&
          isTRUE(cvr_ok)) {
        # Extract 3-way term
        frs_3w <- frs_coefs[grep(
          paste0(
            "YRS_from_bl:FRS:",
            BRAIN_VAR, "|",
            "YRS_from_bl:",
            BRAIN_VAR, ":FRS|",
            "FRS:YRS_from_bl:",
            BRAIN_VAR
          ),
          names(frs_coefs)
        )]
        cvr_3w <- cvr_coefs[grep(
          paste0(
            "YRS_from_bl:CVR_mimic:",
            BRAIN_VAR, "|",
            "YRS_from_bl:",
            BRAIN_VAR, ":CVR_mimic|",
            "CVR_mimic:YRS_from_bl:",
            BRAIN_VAR
          ),
          names(cvr_coefs)
        )]
        if (length(frs_3w) == 1 &&
            length(cvr_3w) == 1) {
          delta_beta.v[b] <-
            frs_3w - cvr_3w
        } else {
          delta_beta.v[b] <- NA
        }
      } else {
        delta_beta.v[b] <- NA
      }
    }

    # Compute BC CI
    valid.v <- delta_beta.v[
      !is.na(delta_beta.v)
    ]

    # Point estimate from original models
    frs_int <- all_stratified.lst[["FRS"]][[
      sex
    ]][[domain]]$interaction
    cvr_int <- all_stratified.lst[[
      "CVR_mimic"
    ]][[sex]][[domain]]$interaction
    delta_est <- frs_int$beta - cvr_int$beta

    if (length(valid.v) >= 100) {
      alpha <- 1 - BOOT_CI_LEVEL
      pb <- mean(valid.v < delta_est)
      pb <- max(0.001, min(0.999, pb))
      z0 <- qnorm(pb)
      za_lo <- qnorm(alpha / 2)
      za_hi <- qnorm(1 - alpha / 2)
      p_lo <- pnorm(2 * z0 + za_lo)
      p_hi <- pnorm(2 * z0 + za_hi)
      ci_lo <- as.numeric(
        quantile(valid.v, probs = p_lo)
      )
      ci_hi <- as.numeric(
        quantile(valid.v, probs = p_hi)
      )
      sig <- (ci_lo > 0 || ci_hi < 0)

      # Bootstrap p-value
      boot_p <- 2 * min(
        mean(valid.v < 0),
        mean(valid.v > 0)
      )

      comparison.lst[[cmp_key]] <- list(
        delta_beta = delta_est,
        ci_lower = ci_lo,
        ci_upper = ci_hi,
        boot_se = sd(valid.v),
        boot_p = boot_p,
        significant = sig,
        n_boot = N_BOOT_LME,
        n_valid = length(valid.v)
      )

      sig_s <- if (sig) " *" else ""
      log_info(
        "  %s %s: dB=%.4f [%.4f, %.4f]%s",
        sex, domain, delta_est,
        ci_lo, ci_hi, sig_s
      )
    } else {
      log_warn(
        "  %s %s: %d valid (need 100)",
        sex, domain, length(valid.v)
      )
      comparison.lst[[cmp_key]] <- list(
        delta_beta = delta_est,
        ci_lower = NA,
        ci_upper = NA,
        boot_se = NA,
        boot_p = NA,
        significant = NA,
        n_boot = N_BOOT_LME,
        n_valid = length(valid.v)
      )
    }

    write_rds_safe(
      comparison.lst[[cmp_key]],
      boot_chk,
      paste("boot checkpoint", cmp_key)
    )
  }
}

# ============================================================
# 8. Save Results
# ============================================================
log_section("Saving Results")

output.lst <- list(
  # Metadata
  brain_measure = BRAIN_VAR,
  brain_label = BRAIN_LABEL,
  cvr_measures = CVR_MEASURES,
  analysis_type = "primary",

  # Primary stratified results (both CVR measures)
  stratified = all_stratified.lst,

  # Bootstrap comparison (FRS vs CVR_mimic)
  comparison = comparison.lst,

  # Combined summary table
  summary_table = summary.dt,

  # CVR_mimic validation
  cvr_validation = list(
    age_correlation = cvr_age_cor,
    expected_null = abs(cvr_age_cor) < 0.1
  ),

  # Methodology documentation
  methodology = list(
    primary_analysis = paste(
      "Sex-stratified (NIH SABV policy)"
    ),
    model_stratified = sprintf(
      "Cognition ~ YRS x CVR x %s + %s + (YRS | PTID)",
      BRAIN_VAR, COVARIATES_STRAT
    ),
    fdr_correction = paste(
      "BH within each CVR measure",
      "(6 tests: 3 domains x 2 sexes)"
    ),
    has_apoe = HAS_APOE,
    rationale = paste(
      "CVR_mimic from MIMIC model partials out age",
      "from cardiovascular risk. If FRS effects are",
      "age-driven, CVR_mimic effects should be null."
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
      stratification = paste(
        "Buckley et al. (2019)",
        "Ann Neurol 85:430-442"
      ),
      mimic_model = paste(
        "Joreskog & Goldberger (1975)",
        "J Econom 3:313-329"
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
  results_dir.path, "lme_hvr_z_results.rds"
)
write_rds_safe(
  output.lst, output.path,
  "LME HVR_z results (FRS + CVR_mimic)"
)
log_info("Saved: %s", output.path)

log_script_end("05_lme_hvr_z.R", success = TRUE)
