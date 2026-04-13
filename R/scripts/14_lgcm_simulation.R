#!/usr/bin/env Rscript

# =============================================================================
# 11_lgcm_simulation.R
# =============================================================================
# Simulation Study: SE-Constrained vs Free Estimation
#
# PURPOSE:
#   METHODOLOGICAL VALIDATION — NOT a power analysis.
#   Demonstrates that fixing residual variance to SE^2
#   (known measurement error from normative models)
#   reduces parameter bias and improves CI coverage
#   compared to freely estimating residual variance.
#
# WHY THIS IS NOT A POWER ANALYSIS:
#   Power analysis would simulate under the actual sample
#   size and plausible effect size range, then report
#   detection rates (% of replications where p < 0.05).
#   This script instead validates the ESTIMATION METHOD
#   by comparing two approaches (constrained vs free)
#   under known truth. If reviewers ask about power,
#   a separate analysis is needed.
#
# WHY THIS MATTERS FOR THE PAPER:
#   The SE-constraint is a key methodological innovation.
#   Normative z-scores have known SEs from the GAMLSS
#   model, so fixing residual variance to SE^2 uses this
#   information rather than discarding it. This simulation
#   shows the approach is unbiased and well-calibrated
#   (correct coverage). Reviewers will ask "why fix
#   residuals?" — this script provides the empirical
#   justification.
#
# DESIGN:
#   1. Generate data with known true parameters + EDT
#   2. Fit models WITH SE constraint (innovation)
#   3. Fit models WITHOUT SE constraint (free estimation)
#   4. Compare: bias, coverage, RMSE
#
# Outputs:
#   - models/results/simulation_results.rds
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(OpenMx)
  library(MASS)
})

# Source utilities
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# Initialize
log_script_start("14_lgcm_simulation.R")
config <- load_config()
set_seed()

# OpenMx settings
mxOption(NULL, "Default optimizer", "SLSQP")

# --- Configuration ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "simulation_study", default = FALSE
)

N_REPS <- get_script_setting("simulation", "n_replications", default = 500)
SAMPLE_SIZES <- get_script_setting("simulation", "sample_sizes",
                                   default = c(200, 500))
N_TIMEPOINTS <- get_script_setting("simulation", "n_timepoints", default = 4)

results_dir <- get_data_path("models", "lgcm_results_dir")
ensure_directory(results_dir)

output_path <- path_join(results_dir, "simulation_results.rds")

# Check if outputs exist
if (!FORCE_REGENERATE && file.exists(output_path)) {
  log_info("Output file exists and force_regenerate=FALSE")
  log_info("Skipping simulation. Set force_regenerate=TRUE to rerun.")
  log_script_end("11_lgcm_simulation.R", success = TRUE)
  quit(status = 0)
}

# =============================================================================
# 1. Define True Parameters
# =============================================================================
log_section("Setting True Parameters")

# True population parameters — calibrated from
# ADNI LGCM parallel process model estimates
# (averaged across 6 sex × domain models).
#
# Derivation:
#   intercept_mean: ~0 (z-scored composites)
#   slope_mean: -0.19 (mean cognitive decline
#     per EDT year across domains)
#   intercept_var: 0.17 (mean latent intercept
#     variance from parallel process models)
#   slope_var: 0.002 (mean latent slope variance;
#     GRID CENTER — varied at 0.5x/1x/1.5x)
#   intercept_slope_cov: 0.0004 (mean i-s cov)
#   hvr_effect: 0.12 (mean b-path coupling;
#     GRID CENTER — varied at 0.5x/1x/1.5x)
#   se_true: 0.29 (mean ADSP-PHC composite SE)
#   se_variation: 0.10 (SD of SEs across obs)
#   edt_mean/sd: 8.5/5.0 (observed EDT dist)
#   edt_change: -1 (approx inter-visit spacing)

# Base (ADNI-calibrated) parameters
BASE_PARAMS <- list(
  intercept_mean = 0,
  slope_mean = -0.19,
  intercept_var = 0.17,
  slope_var = 0.002,
  intercept_slope_cov = 0.0004,
  hvr_effect = 0.12,
  se_true = 0.29,
  se_variation = 0.10,
  edt_mean = 8.5,
  edt_sd = 5.0,
  edt_change = -1
)

# Parameter grid: vary hvr_effect and slope_var
# at 0.5x, 1x, 1.5x observed values to show
# robustness across realistic conditions
GRID_MULTIPLIERS <- c(0.5, 1.0, 1.5)
param_grid.dt <- CJ(
  hvr_mult = GRID_MULTIPLIERS,
  sv_mult = GRID_MULTIPLIERS
)
param_grid.dt[, grid_id := .I]

log_info(
  "Parameter grid: %d conditions",
  nrow(param_grid.dt)
)

# For backward compatibility with single-run code,
# set TRUE_PARAMS to the base (1x/1x) condition
TRUE_PARAMS <- BASE_PARAMS

log_info("True parameters:")
log_info("  Intercept mean: %.2f", TRUE_PARAMS$intercept_mean)
log_info("  Slope mean: %.2f (per year EDT)", TRUE_PARAMS$slope_mean)
log_info("  HVR -> Slope effect: %.2f (KEY PARAMETER)", TRUE_PARAMS$hvr_effect)
log_info("  Measurement SE: %.2f (+/- %.2f)", TRUE_PARAMS$se_true,
         TRUE_PARAMS$se_variation)

# =============================================================================
# 2. Data Generation Function
# =============================================================================
log_section("Setting Up Data Generation")

generate_lgcm_data <- function(n, n_tp, params) {
  # Generate latent growth factors
  Sigma_growth <- matrix(c(
    params$intercept_var, params$intercept_slope_cov,
    params$intercept_slope_cov, params$slope_var
  ), 2, 2)

  growth_factors <- mvrnorm(n, c(params$intercept_mean, params$slope_mean),
                            Sigma_growth)
  intercepts <- growth_factors[, 1]
  slopes <- growth_factors[, 2]

  # Generate covariate (HVR z-score, standardized)
  HVR <- rnorm(n)

  # Add HVR effect to slopes
  slopes <- slopes + params$hvr_effect * HVR

  # Generate EDT for each subject (baseline EDT varies by subject)
  edt_baseline <- rnorm(n, params$edt_mean, params$edt_sd)

  # Generate time-varying data
  data_list <- list()

  for (t in seq_len(n_tp)) {
    # EDT at this timepoint (decreases over time)
    edt_t <- edt_baseline + (t - 1) * params$edt_change

    # True score: intercept + slope * EDT
    # Using EDT relative to baseline for the loading
    edt_loading <- edt_t - edt_baseline
    true_score <- intercepts + slopes * edt_loading

    # Observation-specific SE (varies around true SE)
    se_obs <- pmax(0.1, rnorm(n, params$se_true, params$se_variation))

    # Observed score with measurement error
    obs_score <- true_score + rnorm(n, 0, se_obs)

    data_list[[paste0("Y_T", t - 1)]] <- obs_score
    data_list[[paste0("SE_T", t - 1)]] <- se_obs
    data_list[[paste0("VAR_T", t - 1)]] <- se_obs^2
    data_list[[paste0("EDT_T", t - 1)]] <- edt_loading  # EDT loading for ITVS
  }

  data_list$HVR <- HVR
  data_list$id <- seq_len(n)

  return(as.data.frame(data_list))
}

# =============================================================================
# 3. OpenMx Model Functions
# =============================================================================

create_openmx_model <- function(data, n_tp, use_se_constraint = TRUE,
                                model_name = "LGCM") {
  obs_vars.v <- paste0("Y_T", 0:(n_tp - 1))
  edt_vars.v <- paste0("EDT_T", 0:(n_tp - 1))
  var_cols.v <- paste0("VAR_T", 0:(n_tp - 1))

  # Get SE² values for residuals
  if (use_se_constraint) {
    se2_vals.v <- sapply(var_cols.v, function(col) mean(data[[col]],
                                                         na.rm = TRUE))
  } else {
    se2_vals.v <- rep(NA, n_tp)  # Will be freely estimated
  }

  manifests.v <- c(obs_vars.v, "HVR")
  latents.v <- c("cog_i", "cog_s")

  # Intercept loadings (all = 1)
  intercept_paths <- mxPath(
    from = "cog_i", to = obs_vars.v,
    free = FALSE, values = rep(1, n_tp)
  )

  # Slope loadings - EDT as definition variables (ITVS)
  slope_paths <- mxPath(
    from = "cog_s", to = obs_vars.v,
    free = FALSE, labels = paste0("data.", edt_vars.v)
  )

  # Factor variances
  factor_var <- mxPath(
    from = latents.v, arrows = 2,
    free = TRUE, values = c(0.5, 0.02), labels = c("var_i", "var_s")
  )

  # Factor covariance
  factor_cov <- mxPath(
    from = "cog_i", to = "cog_s", arrows = 2,
    free = TRUE, values = 0, labels = "cov_is"
  )

  # Factor means
  factor_means <- mxPath(
    from = "one", to = latents.v,
    free = TRUE, values = c(0, -0.05), labels = c("mean_i", "mean_s")
  )

  # Residual variances
  if (use_se_constraint) {
    residual_paths <- mxPath(
      from = obs_vars.v, arrows = 2,
      free = FALSE, values = se2_vals.v
    )
  } else {
    residual_paths <- mxPath(
      from = obs_vars.v, arrows = 2,
      free = TRUE, values = rep(0.1, n_tp),
      labels = paste0("resvar_", 0:(n_tp - 1))
    )
  }

  # Zero manifest means
  manifest_means <- mxPath(
    from = "one", to = obs_vars.v,
    free = FALSE, values = rep(0, n_tp)
  )

  # HVR -> Slope (key effect)
  hvr_path <- mxPath(
    from = "HVR", to = "cog_s",
    free = TRUE, values = 0.1, labels = "beta_hvr"
  )

  # HVR exogenous specification
  hvr_var <- mxPath(from = "HVR", arrows = 2, free = FALSE, values = 1)
  hvr_mean <- mxPath(from = "one", to = "HVR", free = FALSE, values = 0)

  # Build model
  model <- mxModel(
    model_name,
    type = "RAM",
    manifestVars = manifests.v,
    latentVars = latents.v,
    mxData(observed = data, type = "raw"),
    intercept_paths, slope_paths,
    factor_var, factor_cov, factor_means,
    residual_paths, manifest_means,
    hvr_path, hvr_var, hvr_mean
  )

  return(model)
}

# =============================================================================
# 4. Helper Function: Extract SE from OpenMx Model
# =============================================================================

extract_param_se <- function(fit, param_name) {

  # Method 1: Extract from summary()$parameters (most reliable)
  se_val <- tryCatch({
    summ <- summary(fit)
    params <- summ$parameters
    row <- params[params$name == param_name, ]
    if (nrow(row) > 0 && !is.na(row$Std.Error)) {
      row$Std.Error
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # Method 2: Fallback to vcov() if summary failed
  if (is.na(se_val)) {
    se_val <- tryCatch({
      vc <- vcov(fit)
      if (!is.null(vc) && param_name %in% rownames(vc)) {
        sqrt(diag(vc))[param_name]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
  }

  return(as.numeric(se_val))
}

# =============================================================================
# 5. Fit One Replication
# =============================================================================

fit_one_replication <- function(n, n_tp, params, rep_id) {
  # Generate data
  dat <- generate_lgcm_data(n, n_tp, params)

  results <- list(rep_id = rep_id, n = n)

  # Model 1: SE-Constrained (innovation)
  model_se <- create_openmx_model(dat, n_tp, use_se_constraint = TRUE,
                                  model_name = "LGCM_SE")

  fit_se <- tryCatch({
    mxTryHard(model_se, silent = TRUE)
  }, error = function(e) NULL)

  if (!is.null(fit_se) && fit_se$output$status$code == 0) {
    results$se_est <- coef(fit_se)["beta_hvr"]
    # Extract SE from summary()$parameters (robust method)
    results$se_se <- extract_param_se(fit_se, "beta_hvr")
    # Compute CI
    results$se_ci_lower <- results$se_est - 1.96 * results$se_se
    results$se_ci_upper <- results$se_est + 1.96 * results$se_se
    results$se_converged <- TRUE
  } else {
    results$se_converged <- FALSE
  }

  # Model 2: Free estimation
  model_free <- create_openmx_model(dat, n_tp, use_se_constraint = FALSE,
                                    model_name = "LGCM_Free")

  fit_free <- tryCatch({
    mxTryHard(model_free, silent = TRUE)
  }, error = function(e) NULL)

  if (!is.null(fit_free) && fit_free$output$status$code == 0) {
    results$free_est <- coef(fit_free)["beta_hvr"]
    # Extract SE from summary()$parameters (robust method)
    results$free_se <- extract_param_se(fit_free, "beta_hvr")
    results$free_ci_lower <- results$free_est - 1.96 * results$free_se
    results$free_ci_upper <- results$free_est + 1.96 * results$free_se
    results$free_converged <- TRUE
  } else {
    results$free_converged <- FALSE
  }

  return(as.data.frame(results))
}

# =============================================================================
# 6. Run Simulation
# =============================================================================
log_section("Running Simulation Study")

all_results <- list()

for (gi in seq_len(nrow(param_grid.dt))) {
  hm <- param_grid.dt$hvr_mult[gi]
  sm <- param_grid.dt$sv_mult[gi]
  gid <- param_grid.dt$grid_id[gi]

  # Build condition-specific params
  cond_params <- BASE_PARAMS
  cond_params$hvr_effect <-
    BASE_PARAMS$hvr_effect * hm
  cond_params$slope_var <-
    BASE_PARAMS$slope_var * sm

  log_info(
    "Grid %d/%d: hvr=%.3f (%.1fx), sv=%.4f (%.1fx)",
    gi, nrow(param_grid.dt),
    cond_params$hvr_effect, hm,
    cond_params$slope_var, sm
  )

  for (n in SAMPLE_SIZES) {
    log_info(
      "  N = %d (%d reps)", n, N_REPS
    )
    results_n <- lapply(
      seq_len(N_REPS), function(rep) {
        if (rep %% 100 == 0) {
          log_info("    Rep %d/%d", rep, N_REPS)
        }
        fit_one_replication(
          n, N_TIMEPOINTS, cond_params, rep
        )
      }
    )
    res.dt <- rbindlist(results_n, fill = TRUE)
    res.dt[, grid_id := gid]
    res.dt[, hvr_mult := hm]
    res.dt[, sv_mult := sm]
    res.dt[, true_hvr := cond_params$hvr_effect]
    res.dt[, true_sv := cond_params$slope_var]
    all_results[[paste(gid, n, sep = "_")]] <-
      res.dt
  }
}

results.dt <- rbindlist(all_results, fill = TRUE)

# =============================================================================
# 7. Compute Summary Statistics
# =============================================================================
log_section("Computing Summary Statistics")

summary_stats <- results.dt[, {
  tv <- true_hvr[1]
  .(
    true_hvr = tv,
    true_sv = true_sv[1],
    hvr_mult = hvr_mult[1],
    sv_mult = sv_mult[1],
    n_converged_se = sum(
      se_converged, na.rm = TRUE
    ),
    n_converged_free = sum(
      free_converged, na.rm = TRUE
    ),
    bias_se = mean(
      se_est - tv, na.rm = TRUE
    ),
    bias_free = mean(
      free_est - tv, na.rm = TRUE
    ),
    rmse_se = sqrt(mean(
      (se_est - tv)^2, na.rm = TRUE
    )),
    rmse_free = sqrt(mean(
      (free_est - tv)^2, na.rm = TRUE
    )),
    avg_se_se = mean(se_se, na.rm = TRUE),
    avg_se_free = mean(
      free_se, na.rm = TRUE
    ),
    coverage_se = mean(
      se_ci_lower <= tv &
        se_ci_upper >= tv,
      na.rm = TRUE
    ),
    coverage_free = mean(
      free_ci_lower <= tv &
        free_ci_upper >= tv,
      na.rm = TRUE
    )
  )
}, by = .(grid_id, n)]

log_info("")
log_info("SIMULATION RESULTS (Grid)")
log_info("")

for (i in seq_len(nrow(summary_stats))) {
  row <- summary_stats[i]
  log_info(
    "Grid %d, N=%d (hvr=%.3f, sv=%.4f):",
    row$grid_id, row$n,
    row$true_hvr, row$true_sv
  )
  log_info(
    "  SE-Con: Bias=%.5f RMSE=%.4f Cov=%.1f%%",
    row$bias_se, row$rmse_se,
    row$coverage_se * 100
  )
  log_info(
    "  Free:   Bias=%.5f RMSE=%.4f Cov=%.1f%%",
    row$bias_free, row$rmse_free,
    row$coverage_free * 100
  )
}

# =============================================================================
# 8. Save Results
# =============================================================================
log_section("Saving Results")

output <- list(
  raw_results = results.dt,
  summary = summary_stats,
  true_params = BASE_PARAMS,
  param_grid = param_grid.dt,
  config = list(
    n_reps = N_REPS,
    sample_sizes = SAMPLE_SIZES,
    n_timepoints = N_TIMEPOINTS,
    grid_multipliers = GRID_MULTIPLIERS,
    methodology = "OpenMx with EDT (ITVS)"
  )
)

write_rds_safe(output, output_path, description = "Simulation results")

log_script_end("11_lgcm_simulation.R", success = TRUE)
