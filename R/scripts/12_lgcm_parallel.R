#!/usr/bin/env Rscript

# =============================================================================
# 09_lgcm_parallel.R
# =============================================================================
# Parallel Process LGCM: HVR and Cognitive Trajectories
#
# METHODOLOGICAL APPROACH:
#   This script follows the principle of "ideal model first, pragmatic fallback."
#
#   1. UNIVARIATE TESTS (scripts 08, 08b) established that BOTH cognitive and
#      HVR trajectories are significantly QUADRATIC (accelerating decline):
#      - MEM: chi-sq diff = 1102.56, p < 10^-241
#      - LAN: chi-sq diff = 680.78, p < 10^-146
#      - EXF: chi-sq diff = 607.48, p < 10^-130
#      - HVR: chi-sq diff = 453.30, p < 10^-97
#
#   2. Therefore, the IDEAL parallel process model would use quadratic growth
#      for both HVR and cognition (6 latent factors).
#
#   3. We ATTEMPT the quadratic parallel process model first (Part 2).
#      If it converges with reasonable estimates: use those results.
#      If it fails (identification issues, huge SEs): document failure, fall back.
#
#   4. If quadratic fails, we FALL BACK to linear parallel process (Part 3).
#      The linear slope captures the average rate of change, which is sufficient
#      for testing whether individual differences in HVR decline correlate with
#      individual differences in cognitive decline.
#
# WHY QUADRATIC PARALLEL PROCESS TYPICALLY FAILS:
#   - 6 latent factors (hvr_i, hvr_s, hvr_q, cog_i, cog_s, cog_q) with full
#     covariance structure = 21 free covariances (plus means, variances)
#   - Sample sizes of N=400-900 are insufficient for identification
#   - Quadratic variance (var_q) is often near-zero, causing boundary issues
#   - References:
#     - Bollen & Curran (2006). Latent Curve
#       Models. Wiley. ISBN: 978-0471455929
#     - Preacher (2010). Latent growth curve
#       models. In Hancock & Mueller (Eds.),
#       Reviewer's Guide to Quantitative
#       Methods. Routledge.
#       DOI: 10.4324/9780203861554-17
#     - Curran, Obeidat & Losardo (2010).
#       Twelve frequently asked questions
#       about growth curve modeling.
#       J Cogn Dev, 11(2), 121-136.
#       DOI: 10.1080/15248371003699969
#
# KEY HYPOTHESIS (tested in Part 3 or Part 2 if successful):
#   Does baseline HVR × FRS interaction predict cognitive slope?
#   (Higher cardiovascular risk weakens protective HVR effect)
#
# COUPLING TEST:
#   hvr_s ~~ cog_s: Do individual differences in HVR decline correlate with
#   individual differences in cognitive decline? (Prerequisite for mediation)
#
# INPUTS:
#   - data/derivatives/lgcm/lgcm_input_full.rds
#   - data/derivatives/lgcm/lgcm_input_male.rds
#   - data/derivatives/lgcm/lgcm_input_female.rds
#
# OUTPUTS:
#   - models/results/lgcm/lgcm_parallel_results.rds (primary results)
#   - models/results/lgcm/lgcm_quadratic_attempt.rds (documents what we tried)
#
# REFERENCES:
#   - Newsom (2015). Longitudinal Structural
#     Equation Modeling. Routledge.
#     DOI: 10.4324/9781315871318
#   - McArdle (2009). Latent variable modeling
#     of differences and changes with
#     longitudinal data. Annu Rev Psychol,
#     60, 577-605.
#     DOI: 10.1146/annurev.psych.60.110707.163612
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(OpenMx)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/disease_time.R"))
source(here("R/utils/openmx_helpers.R"))

log_script_start("12_lgcm_parallel.R")
config <- load_config()
validate_config(config)
validate_packages(c("OpenMx", "data.table"))
set_seed()

# --- Configuration ---
COGNITIVE_DOMAINS <- get_parameter("cognitive_domains")
VALID_TIMEPOINTS <- get_valid_timepoints()
N_TIMEPOINTS <- length(VALID_TIMEPOINTS)

# OpenMx options
mxOption(NULL, "Default optimizer", "SLSQP")
mxOption(NULL, "Number of Threads", 2)

results_dir <- get_data_path("models", "lgcm_results_dir")
ensure_directory(results_dir)

# =============================================================================
# PART 1: DATA PREPARATION
# =============================================================================
log_section("Part 1: Data Preparation")

full.dt <- read_rds_safe(
  get_data_path("derivatives", "lgcm_input_full"), "Full EDT sample"
)
male.dt <- read_rds_safe(
  get_data_path("derivatives", "lgcm_input_male"), "Male EDT sample"
)
female.dt <- read_rds_safe(
  get_data_path("derivatives", "lgcm_input_female"), "Female EDT sample"
)

log_info("Loaded samples:")
log_info("  Full: N = %d", nrow(full.dt))
log_info("  Male: N = %d", nrow(male.dt))
log_info("  Female: N = %d", nrow(female.dt))

# Check MIMIC indicator availability
MIMIC_INDICATORS <- c(
  "SBP_z", "HTN", "GLUCOSE_z",
  "CHOL_z", "CREAT_z"
)
HAS_MIMIC_IND <- all(
  MIMIC_INDICATORS %in% names(full.dt)
)
if (HAS_MIMIC_IND) {
  log_info(
    "MIMIC indicators available for embedding"
  )
} else {
  stop(
    "MIMIC indicators missing. Run script 11."
  )
}

# Variable names
edt_cols.v <- paste0("EDT_", VALID_TIMEPOINTS)
hvr_cols.v <- paste0("HVR_Z_", VALID_TIMEPOINTS)

# -----------------------------------------------------------------------------
# 1a. Handle Missing EDT via Hybrid Extrapolation
# -----------------------------------------------------------------------------
log_info("")
log_info("Extrapolating missing EDT values...")

n_missing_before <- sum(is.na(as.matrix(full.dt[, ..edt_cols.v])))
log_info("  Before: %d missing EDT values", n_missing_before)

full_result.lst <- extrapolate_edt_hybrid(full.dt, valid_timepoints = VALID_TIMEPOINTS)
full.dt <- full_result.lst$data
pop_rate <- full_result.lst$pop_rate

male_result.lst <- extrapolate_edt_hybrid(male.dt, pop_rate = pop_rate,
                                           valid_timepoints = VALID_TIMEPOINTS)
male.dt <- male_result.lst$data

female_result.lst <- extrapolate_edt_hybrid(female.dt, pop_rate = pop_rate,
                                             valid_timepoints = VALID_TIMEPOINTS)
female.dt <- female_result.lst$data

n_missing_after <- sum(is.na(as.matrix(full.dt[, ..edt_cols.v])))
log_info("  After: %d missing EDT values", n_missing_after)

# -----------------------------------------------------------------------------
# 1b. Center EDT for Numerical Stability
# -----------------------------------------------------------------------------
all_edt.v <- unlist(full.dt[, ..edt_cols.v])
edt_mean <- mean(all_edt.v, na.rm = TRUE)
log_info("  EDT grand mean: %.3f years (used for centering)", edt_mean)

for (col in edt_cols.v) {
  full.dt[, (col) := get(col) - edt_mean]
  male.dt[, (col) := get(col) - edt_mean]
  female.dt[, (col) := get(col) - edt_mean]
}

# -----------------------------------------------------------------------------
# 1c. Compute EDT² for Quadratic Models
# -----------------------------------------------------------------------------
edt_sq_cols.v <- paste0(edt_cols.v, "_sq")
for (i in seq_along(edt_cols.v)) {
  col <- edt_cols.v[i]
  sq_col <- edt_sq_cols.v[i]
  full.dt[, (sq_col) := get(col)^2]
  male.dt[, (sq_col) := get(col)^2]
  female.dt[, (sq_col) := get(col)^2]
}

log_info("  EDT centered and squared (for quadratic models)")

# -----------------------------------------------------------------------------
# 1d. Check for APOE4
# -----------------------------------------------------------------------------
HAS_APOE4 <- "APOE4" %in% names(full.dt) && sum(!is.na(full.dt$APOE4)) > 0
if (HAS_APOE4) {
  full.dt[, APOE4_num := as.numeric(APOE4)]
  male.dt[, APOE4_num := as.numeric(APOE4)]
  female.dt[, APOE4_num := as.numeric(APOE4)]
  log_info("  APOE4 available: %.1f%% carriers",
           100 * mean(full.dt$APOE4_num, na.rm = TRUE))
} else {
  log_warn("  APOE4 not available")
}

# =============================================================================
# PART 2: ATTEMPT QUADRATIC PARALLEL PROCESS (IDEAL MODEL)
# =============================================================================
log_section("Part 2: Quadratic Parallel Process (Ideal Model)")

log_info("")
log_info("RATIONALE:")
log_info("  Univariate tests (scripts 08, 08b) showed trajectories are QUADRATIC.")
log_info("  Therefore, we first attempt quadratic parallel process models.")
log_info("")
log_info("MODEL STRUCTURE (6 latent factors):")
log_info("  HVR: hvr_i (intercept) + hvr_s (slope) + hvr_q (quadratic)")
log_info("  COG: cog_i (intercept) + cog_s (slope) + cog_q (quadratic)")
log_info("  Cross-process: hvr_s ~~ cog_s (slope coupling)")
log_info("")

#' Create quadratic parallel process LGCM
#'
#' 6 latent factors with sparse covariance structure to aid identification.
#' Only theoretically-motivated covariances are estimated.
create_quadratic_parallel.fn <- function(data, domain, has_apoe = FALSE) {

  cog_vars.v <- paste0("PHC_", domain, "_", VALID_TIMEPOINTS)
  hvr_vars.v <- paste0("HVR_Z_", VALID_TIMEPOINTS)

  manifests.v <- c(cog_vars.v, hvr_vars.v)
  covariate_vars.v <- c("AGE_bl", "EDUC")
  if (has_apoe) covariate_vars.v <- c(covariate_vars.v, "APOE4_num")
  manifests.v <- c(manifests.v, covariate_vars.v)

  latents.v <- c("hvr_i", "hvr_s", "hvr_q", "cog_i", "cog_s", "cog_q")

  # HVR trajectory paths
  hvr_int <- mxPath(from = "hvr_i", to = hvr_vars.v,
                    free = FALSE, values = rep(1, N_TIMEPOINTS))
  hvr_slope <- mxPath(from = "hvr_s", to = hvr_vars.v,
                      free = FALSE, labels = paste0("data.", edt_cols.v))
  hvr_quad <- mxPath(from = "hvr_q", to = hvr_vars.v,
                     free = FALSE, labels = paste0("data.", edt_sq_cols.v))
  hvr_resid <- mxPath(from = hvr_vars.v, arrows = 2,
                      free = TRUE, values = rep(0.5, N_TIMEPOINTS),
                      labels = paste0("hvr_resid_", VALID_TIMEPOINTS))
  hvr_means <- mxPath(from = "one", to = hvr_vars.v,
                      free = FALSE, values = rep(0, N_TIMEPOINTS))

  # Cognitive trajectory paths
  cog_int <- mxPath(from = "cog_i", to = cog_vars.v,
                    free = FALSE, values = rep(1, N_TIMEPOINTS))
  cog_slope <- mxPath(from = "cog_s", to = cog_vars.v,
                      free = FALSE, labels = paste0("data.", edt_cols.v))
  cog_quad <- mxPath(from = "cog_q", to = cog_vars.v,
                     free = FALSE, labels = paste0("data.", edt_sq_cols.v))
  cog_resid <- mxPath(from = cog_vars.v, arrows = 2,
                      free = TRUE, values = rep(0.3, N_TIMEPOINTS),
                      labels = paste0("cog_resid_", VALID_TIMEPOINTS))
  cog_means <- mxPath(from = "one", to = cog_vars.v,
                      free = FALSE, values = rep(0, N_TIMEPOINTS))

  # Factor means
  factor_means <- mxPath(from = "one", to = latents.v,
                         free = TRUE,
                         values = c(0, -0.05, -0.005, 0, -0.05, -0.01),
                         labels = paste0("mean_", latents.v))

  # Factor variances
  factor_var <- mxPath(from = latents.v, arrows = 2,
                       free = TRUE, values = c(1, 0.1, 0.001, 1, 0.1, 0.001),
                       labels = paste0("var_", latents.v))

  # Within-process covariances (i-s, i-q, s-q for each process)
  hvr_cov_is <- mxPath(from = "hvr_i", to = "hvr_s", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_hvr_is")
  hvr_cov_iq <- mxPath(from = "hvr_i", to = "hvr_q", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_hvr_iq")
  hvr_cov_sq <- mxPath(from = "hvr_s", to = "hvr_q", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_hvr_sq")

  cog_cov_is <- mxPath(from = "cog_i", to = "cog_s", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_cog_is")
  cog_cov_iq <- mxPath(from = "cog_i", to = "cog_q", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_cog_iq")
  cog_cov_sq <- mxPath(from = "cog_s", to = "cog_q", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_cog_sq")

  # CROSS-PROCESS covariances (SPARSE - only key ones)
  cross_int <- mxPath(from = "hvr_i", to = "cog_i", arrows = 2,
                      free = TRUE, values = 0.3, labels = "cov_hvr_i_cog_i")
  cross_slope <- mxPath(from = "hvr_s", to = "cog_s", arrows = 2,
                        free = TRUE, values = 0.01, labels = "cov_hvr_s_cog_s")
  cross_quad <- mxPath(from = "hvr_q", to = "cog_q", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_hvr_q_cog_q")

  # Covariates -> factor means (simplified: only to intercepts)
  age_paths <- mxPath(from = "AGE_bl", to = c("hvr_i", "cog_i"),
                      free = TRUE, values = 0,
                      labels = c("age_hvr_i", "age_cog_i"))
  educ_paths <- mxPath(from = "EDUC", to = c("hvr_i", "cog_i"),
                       free = TRUE, values = 0,
                       labels = c("educ_hvr_i", "educ_cog_i"))

  # Covariate variances/means
  pred_var <- mxPath(from = covariate_vars.v, arrows = 2,
                     free = FALSE, values = rep(1, length(covariate_vars.v)))
  pred_means <- mxPath(from = "one", to = covariate_vars.v,
                       free = FALSE, values = rep(0, length(covariate_vars.v)))

  # Build model
  paths.lst <- list(
    hvr_int, hvr_slope, hvr_quad, hvr_resid, hvr_means,
    cog_int, cog_slope, cog_quad, cog_resid, cog_means,
    factor_means, factor_var,
    hvr_cov_is, hvr_cov_iq, hvr_cov_sq,
    cog_cov_is, cog_cov_iq, cog_cov_sq,
    cross_int, cross_slope, cross_quad,
    pred_var, pred_means, age_paths, educ_paths
  )

  if (has_apoe) {
    apoe_paths <- mxPath(from = "APOE4_num", to = c("hvr_i", "cog_i"),
                         free = TRUE, values = 0,
                         labels = c("apoe_hvr_i", "apoe_cog_i"))
    paths.lst <- c(paths.lst, list(apoe_paths))
  }

  model <- mxModel(
    paste0("QuadParallel_", domain),
    type = "RAM",
    manifestVars = manifests.v,
    latentVars = latents.v,
    mxData(observed = as.data.frame(data), type = "raw")
  )

  for (path in paths.lst) {
    model <- mxModel(model, path)
  }

  return(model)
}

# -----------------------------------------------------------------------------
# Attempt Quadratic Models (Male and Female Samples)
# -----------------------------------------------------------------------------
log_info("Attempting quadratic parallel process (Male/Female)...")
log_info("")

quadratic_results.lst <- list()
quadratic_viable <- FALSE

for (sample_name in c("Male", "Female")) {
  log_info("--- %s sample ---", sample_name)
  sample.dt <- switch(sample_name,
    "Male" = male.dt,
    "Female" = female.dt
  )

  for (domain in COGNITIVE_DOMAINS) {
    key <- paste(sample_name, domain, sep = "_")
    log_info(
      "  %s: Fitting quadratic parallel process...",
      domain
    )

    result <- list(
      sample = sample_name,
      domain = domain,
      attempted = TRUE,
      converged = FALSE,
      viable = FALSE,
      reason = NULL,
      params = NULL
    )

    model_result <- tryCatch({
      list(
        value = create_quadratic_parallel.fn(
          sample.dt, domain, has_apoe = HAS_APOE4
        ),
        error = NULL
      )
    }, error = function(e) {
      list(value = NULL, error = conditionMessage(e))
    })
    model <- model_result$value

    if (!is.null(model_result$error)) {
      result$reason <- paste(
        "Model creation error:", model_result$error
      )
    }

    if (!is.null(model)) {
      fit_result <- tryCatch({
        list(
          value = mxTryHard(
            model, silent = TRUE, extraTries = 10
          ),
          error = NULL
        )
      }, error = function(e) {
        list(
          value = NULL,
          error = conditionMessage(e)
        )
      })
      fit <- fit_result$value

      if (!is.null(fit_result$error)) {
        result$reason <- paste(
          "Fitting error:", fit_result$error
        )
      }

      if (!is.null(fit)) {
        status <- fit$output$status$code
        result$converged <- (status %in% c(0, 1))

        if (result$converged) {
          params <- summary(fit)$parameters
          max_se <- max(
            params$Std.Error, na.rm = TRUE
          )

          if (max_se > 100) {
            result$viable <- FALSE
            result$reason <- sprintf(
              "Unreliable estimates: max SE = %.1f",
              max_se
            )
            log_warn(
              "    Converged but UNRELIABLE: max SE = %.1f",
              max_se
            )
          } else {
            result$viable <- TRUE
            result$params <- params
            quadratic_viable <- TRUE
            log_info(
              "    SUCCESS: Converged with good estimates"
            )
          }
        } else {
          result$reason <- paste(
            "Did not converge:",
            fit$output$status$message
          )
          log_warn("    Did not converge")
        }
      }
    } else {
      log_warn("    Model creation failed")
    }

    quadratic_results.lst[[key]] <- result
  }
}

# -----------------------------------------------------------------------------
# Document Quadratic Attempt Outcome
# -----------------------------------------------------------------------------
log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info("QUADRATIC PARALLEL PROCESS OUTCOME")
log_info(paste(rep("=", 70), collapse = ""))

n_viable <- sum(
  sapply(quadratic_results.lst, function(x) x$viable)
)
n_quad_total <- length(quadratic_results.lst)

if (quadratic_viable) {
  log_info("")
  log_info(
    "  STATUS: %d/%d models converged with viable estimates",
    n_viable, n_quad_total
  )
  log_info("  ACTION: Could use quadratic results (but proceeding to linear for")
  log_info("          comparison and consistency with mediation analysis)")
} else {
  log_info("")
  log_info("  STATUS: Quadratic parallel process FAILED")
  log_info("  REASONS:")
  for (domain in names(quadratic_results.lst)) {
    r <- quadratic_results.lst[[domain]]
    if (!r$viable) {
      log_info("    %s: %s", domain, r$reason)
    }
  }
  log_info("")
  log_info("  This is a known issue with quadratic parallel process models:")
  log_info("    - 6 latent factors require very large N for identification")
  log_info("    - Quadratic variance often near zero (boundary issues)")
  log_info("    - References: Bollen & Curran (2006), Preacher (2010)")
  log_info("")
  log_info("  ACTION: Falling back to LINEAR parallel process (Part 3)")
}

log_info(paste(rep("=", 70), collapse = ""))

# Save quadratic attempt documentation
saveRDS(
  list(
    results = quadratic_results.lst,
    any_viable = quadratic_viable,
    n_viable = n_viable,
    timestamp = Sys.time(),
    note = paste(
      "Documents attempt to fit quadratic parallel process LGCMs.",
      "Univariate tests showed quadratic trajectories, but parallel process",
      "models with 6 latent factors failed due to identification issues."
    )
  ),
  file.path(results_dir, "lgcm_quadratic_attempt.rds")
)
log_info("Saved: lgcm_quadratic_attempt.rds (documents quadratic attempt)")

# =============================================================================
# PART 3: FIXED QUADRATIC PARALLEL PROCESS
# =============================================================================
#
# JUSTIFICATION & REFERENCES:
#   Fully random quadratic parallel process
#   (6 latent factors) failed due to empirical
#   under-identification (Bollen & Curran, 2006,
#   Ch. 4, ISBN: 978-0471455929).
#
#   Fixing growth factor variance to zero when
#   it is non-estimable is standard practice:
#     - Ram & Grimm (2009). Growth mixture
#       modeling. Int J Behav Dev, 33(6),
#       565-576. DOI: 10.1177/0165025409343765
#       ("the offending variance ... was fixed
#       to =0 and the re-estimated model
#       interpreted with caution")
#     - Grimm, Ram & Estabrook (2017). Growth
#       Modeling. Guilford Press.
#       ISBN: 978-1462526062
#
#   The fixed quadratic retains the population-
#   level curvature (free mean) while avoiding
#   the covariance structure expansion that
#   causes identification failure. This adds
#   only 2 free parameters to the linear model.
#
log_section(
  "Part 3: Fixed Quadratic Parallel Process"
)

log_info("")
log_info("RATIONALE:")
log_info(
  "  - Fully random quadratic failed (Part 2)"
)
log_info(
  "  - Fixed quadratic adds population-level"
)
log_info(
  "    curvature (2 free params) to linear base"
)
log_info(
  "  - No random quadratic variance (avoids"
)
log_info(
  "    identification issues)"
)
log_info("")
log_info("MODEL STRUCTURE (4+2 latent factors):")
log_info(
  "  HVR: hvr_i + hvr_s + hvr_q (fixed)"
)
log_info(
  "  COG: cog_i + cog_s + cog_q (fixed)"
)
log_info(
  "  Cross: hvr_i~~cog_i, hvr_s~~cog_s"
)
log_info(
  "  Covariates: FRS, HVR_Z_T0, AGE,"
)
log_info("    EDUC, APOE4")
log_info("")

#' Extract fit indices from a converged OpenMx model
#'
#' Computes CFI, TLI, RMSEA, SRMR via reference models.
#' Returns NULL on failure (e.g., non-positive definite).
#' @param fit A fitted mxModel object
#' @return Named list of fit indices or NULL
extract_fit_indices.fn <- function(fit) {
  # ITVS models (definition variables for EDT)
  # cannot use mxRefModels: it ignores
  # individually-varying loadings, making
  # CFI/TLI/RMSEA/chi-sq meaningless.
  # Only -2LL, AIC, BIC are valid.
  m2ll <- fit$output$fit
  n_param <- length(fit$output$estimate)
  n_obs <- fit$data$numObs
  aic_val <- m2ll + 2 * n_param
  bic_val <- m2ll + log(n_obs) * n_param
  result <- list(
    minus2LL = m2ll,
    n_param = n_param,
    AIC = aic_val,
    BIC = bic_val
  )
  log_info(
    "    Fit: -2LL=%.1f AIC=%.1f BIC=%.1f",
    m2ll, aic_val, bic_val
  )
  result
}

#' Nested model LRT (likelihood ratio test)
#'
#' Fixes a parameter to 0 in the full model,
#' re-fits, and computes LRT. Valid for ITVS
#' because both models share the same definition
#' variable structure.
#'
#' @param full_fit Fitted MxModel (full model)
#' @param fix_param Label of parameter to fix
#' @param fix_value Value to fix to (default 0)
#' @return list(chi_diff, df_diff, p) or NULL
nested_lrt.fn <- function(
    full_fit, fix_param, fix_value = 0) {
  reduced <- tryCatch(
    omxSetParameters(
      full_fit,
      labels = fix_param,
      free = FALSE,
      values = fix_value
    ),
    error = function(e) NULL
  )
  if (is.null(reduced)) return(NULL)
  reduced_fit <- tryCatch(
    mxRun(reduced, silent = TRUE),
    error = function(e) NULL
  )
  if (is.null(reduced_fit)) return(NULL)
  if (reduced_fit$output$status$code > 1) {
    return(NULL)
  }
  chi_diff <- reduced_fit$output$fit -
    full_fit$output$fit
  if (chi_diff < 0) chi_diff <- 0
  p <- pchisq(
    chi_diff, df = 1, lower.tail = FALSE
  )
  list(
    chi_diff = chi_diff,
    df_diff = 1,
    p = p
  )
}

#' Create linear parallel process LGCM
#'
#' Standard 4-factor parallel process model with EDT-based ITVS.
create_linear_parallel.fn <- function(
    data, domain,
    cvr_var = "FRS_sex_centered",
    cvr_label = "frs",
    has_apoe = FALSE) {

  cog_vars.v <- paste0("PHC_", domain, "_", VALID_TIMEPOINTS)
  hvr_vars.v <- paste0("HVR_Z_", VALID_TIMEPOINTS)
  cog_var_cols.v <- paste0("PHC_", domain, "_VAR_", VALID_TIMEPOINTS)

  # Get average SE² for cognitive residual constraints
  avg_cog_se2.v <- sapply(cog_var_cols.v, function(col) {
    if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else 0.04
  })

  manifests.v <- c(cog_vars.v, hvr_vars.v)
  predictor_vars.v <- c(
    cvr_var, "AGE_bl", "EDUC"
  )
  if (has_apoe) predictor_vars.v <- c(predictor_vars.v, "APOE4_num")
  manifests.v <- c(manifests.v, predictor_vars.v)

  latents.v <- c("hvr_i", "hvr_s", "cog_i", "cog_s")

  # HVR trajectory
  hvr_int_paths <- mxPath(from = "hvr_i", to = hvr_vars.v,
                          free = FALSE, values = rep(1, N_TIMEPOINTS))
  hvr_slope_paths <- mxPath(from = "hvr_s", to = hvr_vars.v,
                            free = FALSE, labels = paste0("data.", edt_cols.v))
  hvr_resid <- mxPath(from = hvr_vars.v, arrows = 2,
                      free = TRUE, values = rep(0.5, N_TIMEPOINTS),
                      labels = paste0("hvr_resid_", VALID_TIMEPOINTS))
  hvr_means <- mxPath(from = "one", to = hvr_vars.v,
                      free = FALSE, values = rep(0, N_TIMEPOINTS))

  # Cognitive trajectory
  cog_int_paths <- mxPath(from = "cog_i", to = cog_vars.v,
                          free = FALSE, values = rep(1, N_TIMEPOINTS))
  cog_slope_paths <- mxPath(from = "cog_s", to = cog_vars.v,
                            free = FALSE, labels = paste0("data.", edt_cols.v))
  cog_resid <- mxPath(from = cog_vars.v, arrows = 2,
                      free = FALSE, values = avg_cog_se2.v)  # Fixed to SE²
  cog_means <- mxPath(from = "one", to = cog_vars.v,
                      free = FALSE, values = rep(0, N_TIMEPOINTS))

  # Factor variances
  factor_var <- mxPath(from = latents.v, arrows = 2,
                       free = TRUE, values = c(1, 0.1, 1, 0.1),
                       labels = paste0("var_", latents.v))

  # Within-process covariances
  hvr_cov_is <- mxPath(from = "hvr_i", to = "hvr_s", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_hvr_is")
  cog_cov_is <- mxPath(from = "cog_i", to = "cog_s", arrows = 2,
                       free = TRUE, values = 0, labels = "cov_cog_is")

  # Cross-process covariances
  cross_int <- mxPath(from = "hvr_i", to = "cog_i", arrows = 2,
                      free = TRUE, values = 0.3, labels = "cov_hvr_i_cog_i")
  cross_slope <- mxPath(from = "hvr_s", to = "cog_s", arrows = 2,
                        free = TRUE, values = 0, labels = "cov_hvr_s_cog_s")
  cross_cogi_hvrs <- mxPath(from = "cog_i", to = "hvr_s", arrows = 2,
                            free = TRUE, values = 0, labels = "cov_cog_i_hvr_s")

  # Factor means
  factor_means <- mxPath(from = "one", to = latents.v,
                         free = TRUE, values = c(0, -0.02, 0, -0.05),
                         labels = paste0("mean_", latents.v))

  # Predictors
  pred_var <- mxPath(from = predictor_vars.v, arrows = 2,
                     free = FALSE, values = rep(1, length(predictor_vars.v)))
  pred_means <- mxPath(from = "one", to = predictor_vars.v,
                       free = FALSE, values = rep(0, length(predictor_vars.v)))

  age_paths <- mxPath(from = "AGE_bl", to = latents.v,
                      free = TRUE, values = rep(0, 4),
                      labels = paste0("age_", latents.v))
  educ_paths <- mxPath(from = "EDUC", to = latents.v,
                       free = TRUE, values = rep(0, 4),
                       labels = paste0("educ_", latents.v))

  # CVR and baseline HVR covariate paths
  cvr_paths <- mxPath(
    from = cvr_var,
    to = c("cog_i", "cog_s"),
    free = TRUE, values = c(0, 0),
    labels = c(
      paste0(cvr_label, "_cog_i"),
      paste0(cvr_label, "_cog_s")
    )
  )
  hvr_bl_path <- mxPath(
    from = "HVR_Z_T0", to = "cog_s",
    free = TRUE, values = 0.01,
    labels = "hvr_bl_cog_s"
  )

  # Build model
  paths.lst <- list(
    hvr_int_paths, hvr_slope_paths,
    hvr_resid, hvr_means,
    cog_int_paths, cog_slope_paths,
    cog_resid, cog_means,
    factor_var, hvr_cov_is, cog_cov_is,
    factor_means,
    cross_int, cross_slope, cross_cogi_hvrs,
    pred_var, pred_means,
    age_paths, educ_paths,
    cvr_paths, hvr_bl_path
  )

  if (has_apoe) {
    apoe_paths <- mxPath(from = "APOE4_num", to = latents.v,
                         free = TRUE, values = rep(0, 4),
                         labels = paste0("apoe_", latents.v))
    paths.lst <- c(paths.lst, list(apoe_paths))
  }

  model <- mxModel(
    paste0("LinearParallel_", domain),
    type = "RAM",
    manifestVars = manifests.v,
    latentVars = latents.v,
    mxData(observed = as.data.frame(data), type = "raw")
  )

  for (path in paths.lst) {
    model <- mxModel(model, path)
  }

  return(model)
}

#' Create fixed-quadratic parallel process LGCM
#'
#' Extends the linear model with population-level
#' curvature (no random quadratic variance).
#' Adds 2 phantom latent factors (hvr_q, cog_q)
#' with variance fixed near zero and free means.
#' Only 2 new free parameters vs linear model.
#'
#' @param data Data frame
#' @param domain Cognitive domain
#' @param cvr_var CVR variable name
#' @param cvr_label CVR label prefix
#' @param has_apoe Whether APOE4 is available
#' @return OpenMx RAM model
create_fixedquad_parallel.fn <- function(
    data, domain,
    cvr_var = "FRS_sex_centered",
    cvr_label = "frs",
    has_apoe = FALSE) {

  # Build the full linear model first
  model <- create_linear_parallel.fn(
    data, domain,
    cvr_var = cvr_var,
    cvr_label = cvr_label,
    has_apoe = has_apoe
  )

  # Variable names for loadings
  cog_vars.v <- paste0(
    "PHC_", domain, "_", VALID_TIMEPOINTS
  )
  hvr_vars.v <- paste0(
    "HVR_Z_", VALID_TIMEPOINTS
  )

  # Quadratic loadings = EDT^2 (def vars)
  hvr_q_paths <- mxPath(
    from = "hvr_q", to = hvr_vars.v,
    free = FALSE,
    labels = paste0(
      "data.", edt_sq_cols.v
    )
  )
  cog_q_paths <- mxPath(
    from = "cog_q", to = cog_vars.v,
    free = FALSE,
    labels = paste0(
      "data.", edt_sq_cols.v
    )
  )

  # Variance fixed near 0 (not exactly 0
  # to avoid singularity in mxFactorScores)
  quad_var <- mxPath(
    from = c("hvr_q", "cog_q"),
    arrows = 2, free = FALSE,
    values = c(1e-10, 1e-10)
  )

  # Free means = population-average curvature
  quad_means <- mxPath(
    from = "one",
    to = c("hvr_q", "cog_q"),
    free = TRUE,
    values = c(-0.001, -0.001),
    labels = c("mean_hvr_q", "mean_cog_q")
  )

  # Add new latents + paths to existing model
  # Only pass NEW latent names (OpenMx merges)
  model <- mxModel(
    model,
    name = paste0(
      "FixedQuadParallel_", domain
    ),
    latentVars = c("hvr_q", "cog_q"),
    hvr_q_paths, cog_q_paths,
    quad_var, quad_means
  )
  return(model)
}

#' Create parallel process LGCM with embedded MIMIC
#'
#' Builds a parallel process model where CVR is a
#' latent factor measured by 5 indicators (MIMIC
#' measurement model). Main effects only — no
#' interaction (latent predictor cannot form product).
#'
#' @param data Data frame with indicators + outcomes
#' @param domain Cognitive domain (MEM, LAN, EXF)
#' @param has_apoe Whether APOE4 is available
#' @return OpenMx RAM model
create_mimic_parallel.fn <- function(
    data, domain, has_apoe = FALSE) {

  cog_vars.v <- paste0(
    "PHC_", domain, "_", VALID_TIMEPOINTS
  )
  hvr_vars.v <- paste0(
    "HVR_Z_", VALID_TIMEPOINTS
  )

  # Manifests: outcomes + indicators + covariates
  covariate_vars.v <- c("AGE_bl", "EDUC")
  if (has_apoe) {
    covariate_vars.v <- c(
      covariate_vars.v, "APOE4_num"
    )
  }
  manifests.v <- c(
    cog_vars.v, hvr_vars.v,
    MIMIC_INDICATORS, covariate_vars.v
  )

  latents.v <- c(
    "hvr_i", "hvr_s", "cog_i", "cog_s", "CVR"
  )

  # --- HVR trajectory ---
  hvr_int <- mxPath(
    from = "hvr_i", to = hvr_vars.v,
    free = FALSE,
    values = rep(1, N_TIMEPOINTS)
  )
  hvr_slope <- mxPath(
    from = "hvr_s", to = hvr_vars.v,
    free = FALSE,
    labels = paste0("data.", edt_cols.v)
  )
  hvr_resid <- mxPath(
    from = hvr_vars.v, arrows = 2,
    free = TRUE,
    values = rep(0.5, N_TIMEPOINTS),
    labels = paste0(
      "hvr_resid_", VALID_TIMEPOINTS
    )
  )
  hvr_means <- mxPath(
    from = "one", to = hvr_vars.v,
    free = FALSE,
    values = rep(0, N_TIMEPOINTS)
  )

  # --- Cognitive trajectory ---
  cog_int <- mxPath(
    from = "cog_i", to = cog_vars.v,
    free = FALSE,
    values = rep(1, N_TIMEPOINTS)
  )
  cog_slope <- mxPath(
    from = "cog_s", to = cog_vars.v,
    free = FALSE,
    labels = paste0("data.", edt_cols.v)
  )
  cog_resid <- mxPath(
    from = cog_vars.v, arrows = 2,
    free = TRUE,
    values = rep(0.3, N_TIMEPOINTS),
    labels = paste0(
      "cog_resid_", VALID_TIMEPOINTS
    )
  )
  cog_means <- mxPath(
    from = "one", to = cog_vars.v,
    free = FALSE,
    values = rep(0, N_TIMEPOINTS)
  )

  # --- MIMIC measurement model ---
  # Factor loadings: first fixed to 1 (marker)
  cvr_load_marker <- mxPath(
    from = "CVR", to = "SBP_z",
    free = FALSE, values = 1
  )
  cvr_load_free <- mxPath(
    from = "CVR",
    to = c(
      "HTN", "GLUCOSE_z",
      "CHOL_z", "CREAT_z"
    ),
    free = TRUE,
    values = c(0.5, 0.5, -0.3, 0.3),
    labels = c(
      "lam_HTN", "lam_GLUC",
      "lam_CHOL", "lam_CREAT"
    )
  )

  # CAUSE: Age -> CVR
  age_cvr <- mxPath(
    from = "AGE_bl", to = "CVR",
    free = TRUE, values = 0.1,
    labels = "gamma_age"
  )

  # CVR factor variance (disturbance)
  cvr_var_path <- mxPath(
    from = "CVR", arrows = 2,
    free = TRUE, values = 0.5,
    labels = "var_CVR"
  )

  # Indicator residuals
  ind_resid <- mxPath(
    from = MIMIC_INDICATORS, arrows = 2,
    free = TRUE, values = 0.5,
    labels = paste0(
      "resid_", MIMIC_INDICATORS
    )
  )

  # Residual covariances
  rescov_sbp_chol <- mxPath(
    from = "SBP_z", to = "CHOL_z",
    arrows = 2, free = TRUE, values = 0,
    labels = "rescov_sbp_chol"
  )
  rescov_sbp_creat <- mxPath(
    from = "SBP_z", to = "CREAT_z",
    arrows = 2, free = TRUE, values = 0,
    labels = "rescov_sbp_creat"
  )

  # Indicator means fixed to 0 (z-scored)
  ind_means <- mxPath(
    from = "one", to = MIMIC_INDICATORS,
    free = FALSE, values = 0
  )

  # CVR mean (free)
  cvr_mean <- mxPath(
    from = "one", to = "CVR",
    free = TRUE, values = 0,
    labels = "mean_CVR"
  )

  # --- Growth factor parameters ---
  growth_latents <- c(
    "hvr_i", "hvr_s", "cog_i", "cog_s"
  )
  factor_var <- mxPath(
    from = growth_latents, arrows = 2,
    free = TRUE,
    values = c(1, 0.1, 1, 0.1),
    labels = paste0("var_", growth_latents)
  )
  factor_means <- mxPath(
    from = "one", to = growth_latents,
    free = TRUE,
    values = c(0, -0.02, 0, -0.05),
    labels = paste0("mean_", growth_latents)
  )

  # Within-process covariances
  hvr_cov_is <- mxPath(
    from = "hvr_i", to = "hvr_s",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_hvr_is"
  )
  cog_cov_is <- mxPath(
    from = "cog_i", to = "cog_s",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_cog_is"
  )

  # Cross-process covariances
  cross_int <- mxPath(
    from = "hvr_i", to = "cog_i",
    arrows = 2, free = TRUE,
    values = 0.3,
    labels = "cov_hvr_i_cog_i"
  )
  cross_slope <- mxPath(
    from = "hvr_s", to = "cog_s",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_hvr_s_cog_s"
  )

  # --- Structural paths: CVR main effects ---
  cvr_cog_i <- mxPath(
    from = "CVR", to = "cog_i",
    free = TRUE, values = 0,
    labels = "cvr_cog_i"
  )
  cvr_cog_s <- mxPath(
    from = "CVR", to = "cog_s",
    free = TRUE, values = 0,
    labels = "cvr_cog_s"
  )
  cvr_hvr_s <- mxPath(
    from = "CVR", to = "hvr_s",
    free = TRUE, values = 0,
    labels = "cvr_hvr_s"
  )

  # --- Covariates ---
  pred_var <- mxPath(
    from = covariate_vars.v, arrows = 2,
    free = FALSE,
    values = rep(
      1, length(covariate_vars.v)
    )
  )
  pred_means <- mxPath(
    from = "one", to = covariate_vars.v,
    free = FALSE,
    values = rep(
      0, length(covariate_vars.v)
    )
  )
  age_paths <- mxPath(
    from = "AGE_bl", to = growth_latents,
    free = TRUE, values = rep(0, 4),
    labels = paste0("age_", growth_latents)
  )
  educ_paths <- mxPath(
    from = "EDUC", to = growth_latents,
    free = TRUE, values = rep(0, 4),
    labels = paste0(
      "educ_", growth_latents
    )
  )

  # Build model
  paths.lst <- list(
    hvr_int, hvr_slope, hvr_resid,
    hvr_means,
    cog_int, cog_slope, cog_resid,
    cog_means,
    cvr_load_marker, cvr_load_free,
    age_cvr, cvr_var_path,
    ind_resid,
    rescov_sbp_chol, rescov_sbp_creat,
    ind_means, cvr_mean,
    factor_var, factor_means,
    hvr_cov_is, cog_cov_is,
    cross_int, cross_slope,
    cvr_cog_i, cvr_cog_s, cvr_hvr_s,
    pred_var, pred_means,
    age_paths, educ_paths
  )

  if (has_apoe) {
    apoe_paths <- mxPath(
      from = "APOE4_num",
      to = growth_latents,
      free = TRUE, values = rep(0, 4),
      labels = paste0(
        "apoe_", growth_latents
      )
    )
    paths.lst <- c(
      paths.lst, list(apoe_paths)
    )
  }

  model <- mxModel(
    paste0("MIMICParallel_", domain),
    type = "RAM",
    manifestVars = manifests.v,
    latentVars = latents.v,
    mxData(
      observed = as.data.frame(data),
      type = "raw"
    )
  )

  for (path in paths.lst) {
    model <- mxModel(model, path)
  }

  return(model)
}

# -----------------------------------------------------------------------------
# Fit Fixed Quadratic Models
# -----------------------------------------------------------------------------
# Also fit linear (nested) for LRT comparison.
# All results extracted from fixed-quad model.
# -----------------------------------------------------------------------------
linear_results.lst <- list()

for (sample_name in c("Male", "Female")) {
  log_info("")
  log_info(
    "=== %s SAMPLE ===",
    toupper(sample_name)
  )

  data <- switch(
    sample_name,
    "Male" = male.dt,
    "Female" = female.dt
  )

  for (domain in COGNITIVE_DOMAINS) {
    key <- paste(
      sample_name, domain, sep = "_"
    )
    log_info(
      "  %s (N=%d)", domain, nrow(data)
    )

    result <- list(
      sample = sample_name,
      domain = domain,
      n = nrow(data),
      converged = FALSE,
      model_type = "fixed_quadratic",
      params = NULL,
      coupling = NULL,
      quad_means = NULL,
      quad_comparison = NULL
    )

    tryCatch({
      # --- Fit fixed-quad model ---
      fq_model <- create_fixedquad_parallel.fn(
        data, domain,
        has_apoe = HAS_APOE4
      )
      fq_fit <- mxTryHard(
        fq_model, silent = TRUE,
        extraTries = 5
      )

      if (fq_fit$output$status$code != 0) {
        log_warn("    Did not converge")
        linear_results.lst[[key]] <- result
        next
      }

      result$converged <- TRUE

      # --- Fit indices ---
      result$fit_indices <-
        extract_fit_indices.fn(fq_fit)

      # --- Parameters ---
      s_tmp <- tryCatch(
        summary(fq_fit),
        error = function(e) NULL
      )
      params <- if (!is.null(s_tmp)) {
        as.data.table(s_tmp$parameters)
      } else {
        NULL
      }
      result$params <- params

      # --- Quadratic means ---
      if (!is.null(params)) {
        hvr_q_row <- params[
          name == "mean_hvr_q"
        ]
        cog_q_row <- params[
          name == "mean_cog_q"
        ]
        result$quad_means <- list(
          hvr_q = hvr_q_row$Estimate,
          hvr_q_se = hvr_q_row$Std.Error,
          cog_q = cog_q_row$Estimate,
          cog_q_se = cog_q_row$Std.Error
        )
        log_info(
          "    hvr_q=%.5f (SE=%.5f)",
          hvr_q_row$Estimate,
          hvr_q_row$Std.Error
        )
        log_info(
          "    cog_q=%.5f (SE=%.5f)",
          cog_q_row$Estimate,
          cog_q_row$Std.Error
        )
      }

      # --- Coupling ---
      cov_row <- if (!is.null(params)) {
        params[
          name == "cov_hvr_s_cog_s"
        ]
      } else {
        data.table()
      }
      if (nrow(cov_row) > 0) {
        z_val <- cov_row$Estimate /
          cov_row$Std.Error
        p_val <- 2 * pnorm(-abs(z_val))
        # Slope variances for correlation
        v_hvr <- params[
          name == "var_hvr_s", Estimate
        ]
        v_cog <- params[
          name == "var_cog_s", Estimate
        ]
        cor_val <- cov_row$Estimate /
          sqrt(v_hvr * v_cog)
        result$coupling <- list(
          cov = cov_row$Estimate,
          se = cov_row$Std.Error,
          z = z_val, p = p_val,
          cor = cor_val
        )
        sig <- ifelse(
          p_val < 0.05, "*", ""
        )
        log_info(
          "    hvr_s~~cog_s: %.4f, p=%.4e %s",
          cov_row$Estimate, p_val, sig
        )
      }

      # --- Nested LRT: coupling ---
      lrt_coupling <- nested_lrt.fn(
        fq_fit, "cov_hvr_s_cog_s"
      )
      if (!is.null(lrt_coupling)) {
        result$lrt_coupling <- lrt_coupling
        log_info(
          "    LRT coupling: X2=%.2f, p=%.4f",
          lrt_coupling$chi_diff,
          lrt_coupling$p
        )
      }

      # --- LRT vs linear (fix quad to 0) ---
      lin_nested <- tryCatch({
        reduced <- omxSetParameters(
          fq_fit,
          labels = c(
            "mean_hvr_q", "mean_cog_q"
          ),
          free = FALSE, values = 0
        )
        mxRun(reduced, silent = TRUE)
      }, error = function(e) NULL)

      if (!is.null(lin_nested) &&
          lin_nested$output$status$code == 0
      ) {
        lin_m2ll <- lin_nested$output$fit
        lin_np <- length(
          lin_nested$output$estimate
        )
        lin_nobs <- lin_nested$data$numObs
        lin_aic <- lin_m2ll + 2 * lin_np

        fq_m2ll <- result$fit_indices$minus2LL
        fq_aic <- result$fit_indices$AIC

        chi_diff <- lin_m2ll - fq_m2ll
        if (chi_diff < 0) chi_diff <- 0
        lrt_p <- pchisq(
          chi_diff, df = 2,
          lower.tail = FALSE
        )
        aic_diff <- lin_aic - fq_aic

        result$quad_comparison <- list(
          chi_diff = chi_diff,
          df = 2,
          p = lrt_p,
          aic_linear = lin_aic,
          aic_fixedquad = fq_aic,
          aic_diff = aic_diff
        )
        log_info(
          "    vs linear: X2=%.1f, p=%.2e, dAIC=%.1f",
          chi_diff, lrt_p, aic_diff
        )
      }

      # --- Factor scores ---
      fs_res <- tryCatch({
        mxFactorScores(
          fq_fit,
          type = "Regression",
          minManifests = 2
        )
      }, error = function(e) {
        log_warn(
          "    Factor scores failed: %s",
          e$message
        )
        NULL
      })
      if (!is.null(fs_res)) {
        scores.mat <- fs_res[, , 1]
        se.mat <- fs_res[, , 2]
        fs.dt <- data.table(
          PTID = data$PTID,
          hvr_i = scores.mat[, "hvr_i"],
          hvr_s = scores.mat[, "hvr_s"],
          cog_i = scores.mat[, "cog_i"],
          cog_s = scores.mat[, "cog_s"],
          hvr_i_se = se.mat[, "hvr_i"],
          hvr_s_se = se.mat[, "hvr_s"],
          cog_i_se = se.mat[, "cog_i"],
          cog_s_se = se.mat[, "cog_s"]
        )
        result$factor_scores <- fs.dt
        log_info(
          "    Factor scores: N=%d",
          nrow(fs.dt)
        )
      }
    }, error = function(e) {
      log_warn("    Error: %s", e$message)
    })

    linear_results.lst[[key]] <- result
  }
}

# =============================================================================
# PART 4: SUMMARY
# =============================================================================
log_section("Part 4: Summary")

log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info("PARALLEL PROCESS LGCM RESULTS")
log_info(paste(rep("=", 70), collapse = ""))
log_info("")
log_info("METHODOLOGY:")
log_info(
  "  1. Attempted QUADRATIC parallel process"
)
if (quadratic_viable) {
  log_info(
    "     -> Some converged (see rds)"
  )
} else {
  log_info(
    "     -> FAILED (under-identification)"
  )
}
log_info(
  "  2. Fit FIXED QUADRATIC parallel process"
)
n_conv <- sum(vapply(
  linear_results.lst,
  function(r) r$converged, logical(1)
))
log_info(
  "     -> %d/%d converged",
  n_conv, length(linear_results.lst)
)
log_info("")

log_info("KEY TEST 1: Slope Coupling (hvr_s ~~ cog_s)")
log_info("  Do HVR decline and cognitive decline correlate?")
log_info("")
for (key in names(linear_results.lst)) {
  r <- linear_results.lst[[key]]
  if (!r$converged || is.null(r$coupling)) next
  sig <- ifelse(r$coupling$p < 0.05, "*", "")
  log_info("  %s: cov=%.4f, p=%.2e %s", key, r$coupling$cov, r$coupling$p, sig)
}

log_info("")
log_info(paste(rep("=", 70), collapse = ""))

# =============================================================================
# PART 5: SAVE RESULTS
# =============================================================================
log_section("Part 5: Saving Results")

# Summary table
summary_rows.lst <- list()
for (key in names(linear_results.lst)) {
  r <- linear_results.lst[[key]]
  mt <- if (is.null(r$model_type)) {
    "linear"
  } else {
    r$model_type
  }
  summary_rows.lst[[
    length(summary_rows.lst) + 1
  ]] <- data.table(
    Sample = r$sample,
    Domain = r$domain,
    N = r$n,
    Converged = r$converged,
    Model_type = mt,
    Coupling_cov = if (
      !is.null(r$coupling)
    ) {
      r$coupling$cov
    } else {
      NA
    },
    Coupling_p = if (
      !is.null(r$coupling)
    ) {
      r$coupling$p
    } else {
      NA
    }
  )
}
summary.dt <- rbindlist(summary_rows.lst)

output.lst <- list(
  linear_results = linear_results.lst,
  quadratic_attempted = TRUE,
  quadratic_viable = quadratic_viable,
  summary = summary.dt,
  methodology = list(
    approach = paste(
      "Fixed quadratic parallel process.",
      "Population-level curvature,",
      "no random quadratic variance."
    ),
    full_quadratic_outcome = if (
      quadratic_viable
    ) "Some converged" else "Failed",
    fixedquad_justification = paste(
      "Fully random quadratic (6 latent",
      "factors) failed due to empirical",
      "under-identification. Fixed",
      "quadratic adds 2 free params",
      "(population curvature) to linear",
      "base without covariance structure",
      "expansion."
    )
  ),
  config = list(
    timepoints = VALID_TIMEPOINTS,
    edt_mean = edt_mean,
    has_apoe = HAS_APOE4,
    timestamp = Sys.time()
  ),
  references = c(
    paste(
      "Bollen & Curran (2006).",
      "Latent Curve Models.",
      "Wiley. ISBN: 978-0471455929"
    ),
    paste(
      "Curran, Obeidat & Losardo (2010).",
      "Twelve frequently asked questions",
      "about growth curve modeling.",
      "J Cogn Dev, 11(2), 121-136.",
      "DOI: 10.1080/15248371003699969"
    ),
    paste(
      "Newsom (2015). Longitudinal",
      "Structural Equation Modeling.",
      "Routledge.",
      "DOI: 10.4324/9781315871318"
    ),
    paste(
      "Preacher (2010). Latent growth",
      "curve models. In Hancock &",
      "Mueller (Eds.), Reviewer's Guide",
      "to Quantitative Methods.",
      "Routledge.",
      "DOI: 10.4324/9780203861554-17"
    )
  ),
  fit_assessment = list(
    valid_indices = c("-2LL", "AIC", "BIC"),
    nested_tests = "coupling",
    limitations = c(
      paste(
        "Incremental fit indices",
        "(CFI, TLI, RMSEA) are not",
        "available for LGCM with",
        "individually-varying time",
        "scores (ITVS/definition",
        "variables). The saturated",
        "reference model required by",
        "these indices cannot",
        "accommodate subject-specific",
        "slope loadings (Grimm, Ram &",
        "Estabrook, 2017,",
        "ISBN: 978-1462526062;",
        "Mehta & West, 2000,",
        "DOI: 10.1037/1082-989X.5.1.23)."
      ),
      paste(
        "Model adequacy was assessed",
        "via: (a) information criteria",
        "(AIC, BIC) for relative model",
        "comparison, (b) nested",
        "likelihood ratio tests for",
        "specific structural paths",
        "(coupling),",
        "(c) parameter estimate",
        "plausibility and convergence,",
        "and (d) simulation-based",
        "parameter recovery",
        "(Script 14)."
      ),
      paste(
        "Fixed quadratic parallel",
        "process was adopted after",
        "fully random quadratic",
        "models failed due to",
        "empirical under-identification",
        "(near-zero quadratic variance;",
        "Preacher, 2010,",
        "DOI: 10.4324/9780203861554-17).",
        "Population-level curvature",
        "is captured; individual",
        "quadratic variation is not."
      )
    )
  )
)

saveRDS(output.lst, file.path(results_dir, "lgcm_parallel_results.rds"))
log_info("Saved: lgcm_parallel_results.rds")

log_script_end("12_lgcm_parallel.R")
