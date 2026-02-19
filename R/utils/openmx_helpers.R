#!/usr/bin/env Rscript

# =============================================================================
# openmx_helpers.R
# =============================================================================
# Utility functions for OpenMx LGCM models with Individually-Varying Time Scores
#
# Purpose:
#   Validate lavaan LGCM results using OpenMx's definition variable capability
#   for true individually-varying time scores (ITVS)
#
# Key Advantage over lavaan:
#   OpenMx supports definition variables, allowing factor loadings to vary
#   by person based on their actual disease time values
#
# Functions:
#   - create_bivariate_lgcm_openmx(): Create bivariate LGCM with ITVS
#   - fit_openmx_model(): Fit model with error handling
#   - extract_openmx_params(): Extract parameter estimates
#   - compare_lavaan_openmx(): Compare results between packages
# =============================================================================

# Check and load OpenMx
if (!requireNamespace("OpenMx", quietly = TRUE)) {
  stop("OpenMx package required for ITVS validation. Install with: install.packages('OpenMx')")
}

suppressPackageStartupMessages({
  library(OpenMx)
})

# =============================================================================
# Helper Operators
# =============================================================================

#' Null coalescing operator
#'
#' Returns y if x is NULL, otherwise returns x.
#' @param x Value to check
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# Model Creation
# =============================================================================

#' Create bivariate LGCM with individually-varying time scores
#'
#' Uses OpenMx definition variables to allow person-specific time loadings
#' based on actual disease time values.
#'
#' @param data Data frame with wide-format data
#' @param cog_vars Character vector of cognitive variable names
#' @param hvr_vars Character vector of HVR variable names
#' @param time_vars Character vector of disease time variable names (definition vars)
#' @param covariates Character vector of covariate names
#' @param name Model name
#' @return OpenMx model object
create_bivariate_lgcm_openmx <- function(data,
                                          cog_vars,
                                          hvr_vars,
                                          time_vars,
                                          covariates = c("AGE_bl", "EDUC"),
                                          name = "Bivariate_LGCM_ITVS") {

  n_tp <- length(cog_vars)
  n_cov <- length(covariates)

  # Manifest variables
  manifest_vars <- c(cog_vars, hvr_vars)

  # Latent variables
  latent_vars <- c("cog_i", "cog_s", "hvr_i", "hvr_s")

  # --- Factor Loadings ---

  # Cognitive intercept: all 1s
  cog_i_loadings <- mxPath(
    from = "cog_i",
    to = cog_vars,
    free = FALSE,
    values = rep(1, n_tp)
  )

  # Cognitive slope: definition variables for ITVS
  # Each person's time values are stored in time_vars columns
  cog_s_loadings <- mxPath(
    from = "cog_s",
    to = cog_vars,
    free = FALSE,
    labels = paste0("data.", time_vars)  # Definition variable syntax
  )

  # HVR intercept: all 1s
  hvr_i_loadings <- mxPath(
    from = "hvr_i",
    to = hvr_vars,
    free = FALSE,
    values = rep(1, n_tp)
  )

  # HVR slope: same definition variables
  hvr_s_loadings <- mxPath(
    from = "hvr_s",
    to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", time_vars)
  )

  # --- Latent Variable Means ---
  latent_means <- mxPath(
    from = "one",
    to = latent_vars,
    free = TRUE,
    values = c(0, 0, 0, 0),
    labels = c("mean_cog_i", "mean_cog_s", "mean_hvr_i", "mean_hvr_s")
  )

  # --- Latent Variable Variances ---
  latent_variances <- mxPath(
    from = latent_vars,
    arrows = 2,
    free = TRUE,
    values = c(1, 0.1, 1, 0.1),
    labels = c("var_cog_i", "var_cog_s", "var_hvr_i", "var_hvr_s")
  )

  # --- Within-Construct Covariances ---
  cog_covariance <- mxPath(
    from = "cog_i",
    to = "cog_s",
    arrows = 2,
    free = TRUE,
    values = 0,
    labels = "cov_cog_is"
  )

  hvr_covariance <- mxPath(
    from = "hvr_i",
    to = "hvr_s",
    arrows = 2,
    free = TRUE,
    values = 0,
    labels = "cov_hvr_is"
  )

  # --- Cross-Domain Covariances ---
  cross_intercept_cov <- mxPath(
    from = "cog_i",
    to = "hvr_i",
    arrows = 2,
    free = TRUE,
    values = 0,
    labels = "cov_cog_hvr_i"
  )

  # --- KEY HYPOTHESIS: Cross-Domain Effects ---
  # HVR intercept and slope predict cognitive slope
  cross_effects <- mxPath(
    from = c("hvr_i", "hvr_s"),
    to = "cog_s",
    arrows = 1,
    free = TRUE,
    values = c(0, 0),
    labels = c("beta_hvr_i", "beta_hvr_s")
  )

  # --- Residual Variances ---
  cog_residuals <- mxPath(
    from = cog_vars,
    arrows = 2,
    free = TRUE,
    values = rep(0.5, n_tp),
    labels = paste0("resid_cog_", seq_len(n_tp))
  )

  hvr_residuals <- mxPath(
    from = hvr_vars,
    arrows = 2,
    free = TRUE,
    values = rep(0.5, n_tp),
    labels = paste0("resid_hvr_", seq_len(n_tp))
  )

  # --- Covariate Effects (if any) ---
  covariate_paths <- list()
  if (length(covariates) > 0) {
    # Add covariates as predictors of all latent variables
    for (cov in covariates) {
      covariate_paths[[cov]] <- mxPath(
        from = cov,
        to = latent_vars,
        arrows = 1,
        free = TRUE,
        values = rep(0, 4),
        labels = paste0(cov, "_on_", latent_vars)
      )
    }

    # Covariate means (fixed to 0, data centered)
    cov_means <- mxPath(
      from = "one",
      to = covariates,
      free = FALSE,
      values = rep(0, n_cov)
    )

    # Covariate variances
    cov_vars <- mxPath(
      from = covariates,
      arrows = 2,
      free = TRUE,
      values = rep(1, n_cov),
      labels = paste0("var_", covariates)
    )

    covariate_paths[["means"]] <- cov_means
    covariate_paths[["vars"]] <- cov_vars
  }

  # --- Combine into Model ---
  model <- mxModel(
    name,
    type = "RAM",
    manifestVars = manifest_vars,
    latentVars = latent_vars,
    mxData(observed = data, type = "raw"),

    # Factor structure
    cog_i_loadings,
    cog_s_loadings,
    hvr_i_loadings,
    hvr_s_loadings,

    # Latent means and variances
    latent_means,
    latent_variances,

    # Covariances
    cog_covariance,
    hvr_covariance,
    cross_intercept_cov,

    # Cross-domain effects
    cross_effects,

    # Residuals
    cog_residuals,
    hvr_residuals
  )

  # Add covariate paths if present
  if (length(covariate_paths) > 0) {
    for (path in covariate_paths) {
      model <- mxModel(model, path)
    }
  }

  model
}


# =============================================================================
# Model Fitting
# =============================================================================

#' Fit OpenMx model with error handling
#'
#' @param model OpenMx model object
#' @param intervals Compute confidence intervals
#' @return List with fit object and diagnostics
fit_openmx_model <- function(model, intervals = FALSE) {

  result <- list(
    fit = NULL,
    converged = FALSE,
    status = NA,
    issues = character(0)
  )

  tryCatch({
    # Fit the model
    fit <- mxRun(model, intervals = intervals, silent = TRUE)

    result$fit <- fit
    result$status <- fit$output$status$code

    # Check convergence (status 0 = OK)
    if (result$status == 0) {
      result$converged <- TRUE
    } else {
      result$issues <- c(result$issues,
                         paste("Optimizer status:", result$status))
    }

    # Check for negative variances
    params <- omxGetParameters(fit)
    var_params <- params[grepl("^var_|^resid_", names(params))]
    if (any(var_params < 0)) {
      neg_vars <- names(var_params)[var_params < 0]
      result$issues <- c(result$issues,
                         paste("Negative variance(s):", paste(neg_vars, collapse = ", ")))
    }

  }, error = function(e) {
    result$issues <- c(result$issues, paste("Fitting error:", e$message))
  })

  result
}


# =============================================================================
# Parameter Extraction
# =============================================================================

#' Extract parameter estimates from OpenMx fit
#'
#' @param fit OpenMx fit object
#' @return Data frame with parameter estimates
extract_openmx_params <- function(fit) {

  if (is.null(fit)) {
    return(data.frame())
  }

  # Get all parameters
  params <- omxGetParameters(fit)

  # Try to get standard errors
  se <- tryCatch({
    summary(fit)$parameters$Std.Error
  }, error = function(e) {
    rep(NA, length(params))
  })

  # Get parameter names from summary
  param_names <- tryCatch({
    summary(fit)$parameters$name
  }, error = function(e) {
    names(params)
  })

  data.frame(
    name = names(params),
    estimate = as.numeric(params),
    se = se[match(names(params), param_names)],
    stringsAsFactors = FALSE
  )
}


#' Extract cross-domain effects from OpenMx fit
#'
#' @param fit OpenMx fit object
#' @return Data frame with cross-domain effects
extract_cross_effects_openmx <- function(fit) {

  params <- extract_openmx_params(fit)

  cross_params <- params[grepl("^beta_hvr", params$name), ]

  if (nrow(cross_params) > 0) {
    cross_params$effect_type <- ifelse(
      grepl("_i$", cross_params$name), "intercept", "slope"
    )
  }

  cross_params
}


# =============================================================================
# Comparison Functions
# =============================================================================

#' Compare lavaan and OpenMx results
#'
#' @param lavaan_fit Fitted lavaan model
#' @param openmx_fit Fitted OpenMx model
#' @return Data frame comparing key parameters
compare_lavaan_openmx <- function(lavaan_fit, openmx_fit) {

  comparison <- data.frame(
    parameter = character(),
    lavaan_est = numeric(),
    lavaan_se = numeric(),
    openmx_est = numeric(),
    openmx_se = numeric(),
    difference = numeric(),
    pct_diff = numeric(),
    stringsAsFactors = FALSE
  )

  if (is.null(lavaan_fit) || is.null(openmx_fit)) {
    return(comparison)
  }

  # Extract lavaan parameters
  lavaan_params <- lavaan::parameterEstimates(lavaan_fit)

  # Focus on cross-domain effects (key hypothesis)
  lavaan_hvr <- lavaan_params[
    lavaan_params$lhs == "cog_s" &
      lavaan_params$op == "~" &
      grepl("hvr", lavaan_params$rhs),
  ]

  # Extract OpenMx parameters
  openmx_params <- extract_openmx_params(openmx_fit)
  openmx_hvr <- openmx_params[grepl("^beta_hvr", openmx_params$name), ]

  # Match and compare
  for (i in seq_len(nrow(lavaan_hvr))) {
    lav_label <- lavaan_hvr$label[i]
    lav_est <- lavaan_hvr$est[i]
    lav_se <- lavaan_hvr$se[i]

    # Find matching OpenMx parameter
    omx_match <- openmx_hvr[openmx_hvr$name == lav_label, ]

    if (nrow(omx_match) > 0) {
      omx_est <- omx_match$estimate[1]
      omx_se <- omx_match$se[1]

      comparison <- rbind(comparison, data.frame(
        parameter = lav_label,
        lavaan_est = lav_est,
        lavaan_se = lav_se,
        openmx_est = omx_est,
        openmx_se = omx_se,
        difference = lav_est - omx_est,
        pct_diff = 100 * (lav_est - omx_est) / abs(lav_est),
        stringsAsFactors = FALSE
      ))
    }
  }

  comparison
}


#' Generate concordance summary
#'
#' @param comparison Data frame from compare_lavaan_openmx
#' @param threshold Acceptable percentage difference
#' @return List with concordance metrics
assess_concordance <- function(comparison, threshold = 10) {

  if (nrow(comparison) == 0) {
    return(list(
      concordant = FALSE,
      message = "No parameters to compare"
    ))
  }

  max_pct_diff <- max(abs(comparison$pct_diff), na.rm = TRUE)
  mean_pct_diff <- mean(abs(comparison$pct_diff), na.rm = TRUE)

  concordant <- max_pct_diff < threshold

  list(
    concordant = concordant,
    max_pct_diff = max_pct_diff,
    mean_pct_diff = mean_pct_diff,
    message = if (concordant) {
      sprintf("Results concordant (max diff: %.1f%%, threshold: %d%%)",
              max_pct_diff, threshold)
    } else {
      sprintf("Results DIVERGENT (max diff: %.1f%% exceeds threshold: %d%%)",
              max_pct_diff, threshold)
    }
  )
}


# =============================================================================
# Model Diagnostics (Heywood Cases, SE Issues, etc.)
# =============================================================================

#' Check OpenMx model for Heywood cases and other diagnostic issues
#'
#' @param model.fit Fitted OpenMx model
#' @param model_name Optional name for logging
#' @param log_fn Logging function (default: message)
#' @return List with diagnostic results
#' @export
check_openmx_diagnostics.fn <- function(model.fit, model_name = NULL,
                                        log_fn = message) {
  name_str <- if (!is.null(model_name)) paste0("[", model_name, "] ") else ""

  result.lst <- list(
    model_name = model_name,
    converged = FALSE,
    status_code = NA,
    status_msg = "",
    heywood_cases = character(0),
    large_se = character(0),
    boundary_estimates = character(0),
    warnings = character(0),
    is_admissible = TRUE
  )

 # Check if fit object is valid
  if (is.null(model.fit)) {
    result.lst$warnings <- c(result.lst$warnings, "Model fit is NULL")
    result.lst$is_admissible <- FALSE
    return(result.lst)
  }

  # Status check
  result.lst$status_code <- model.fit$output$status$code
  result.lst$status_msg <- model.fit$output$status$message %||% ""
  result.lst$converged <- (result.lst$status_code == 0)

  if (!result.lst$converged) {
    result.lst$warnings <- c(result.lst$warnings,
      sprintf("Non-convergence (status=%d): %s",
              result.lst$status_code, result.lst$status_msg))
    result.lst$is_admissible <- FALSE
    return(result.lst)
  }

  # Extract parameters
  summ <- tryCatch(summary(model.fit), error = function(e) NULL)
  if (is.null(summ)) {
    result.lst$warnings <- c(result.lst$warnings, "Could not extract summary")
    return(result.lst)
  }

  params.dt <- as.data.table(summ$parameters)

  # --- Check 1: Heywood Cases (negative variances) ---
  var_params.dt <- params.dt[grepl("^var_|_resid_|Var\\[", name)]
  if (nrow(var_params.dt) > 0) {
    neg_var.dt <- var_params.dt[Estimate < 0]
    if (nrow(neg_var.dt) > 0) {
      for (i in seq_len(nrow(neg_var.dt))) {
        msg <- sprintf("HEYWOOD: %s = %.6f (negative variance)",
                       neg_var.dt$name[i], neg_var.dt$Estimate[i])
        result.lst$heywood_cases <- c(result.lst$heywood_cases, msg)
        log_fn(paste0(name_str, "WARNING: ", msg))
      }
      result.lst$is_admissible <- FALSE
    }
  }

  # --- Check 2: Very large SEs (potential identification issues) ---
  if ("Std.Error" %in% names(params.dt)) {
    # Flag SEs > 10x the estimate (for non-tiny estimates) or SE > 100
    large_se.dt <- params.dt[!is.na(Std.Error) & abs(Estimate) > 0.001 &
                             (Std.Error > abs(Estimate) * 10 | Std.Error > 100)]
    if (nrow(large_se.dt) > 0) {
      for (i in seq_len(nrow(large_se.dt))) {
        msg <- sprintf("LARGE_SE: %s (SE=%.4f, Est=%.4f)",
                       large_se.dt$name[i], large_se.dt$Std.Error[i],
                       large_se.dt$Estimate[i])
        result.lst$large_se <- c(result.lst$large_se, msg)
        log_fn(paste0(name_str, "WARNING: ", msg))
      }
    }
  }

  # --- Check 3: Boundary estimates (variances near zero) ---
  if (nrow(var_params.dt) > 0) {
    boundary.dt <- var_params.dt[Estimate >= 0 & Estimate < 1e-6]
    if (nrow(boundary.dt) > 0) {
      for (i in seq_len(nrow(boundary.dt))) {
        msg <- sprintf("BOUNDARY: %s = %.2e (near-zero variance)",
                       boundary.dt$name[i], boundary.dt$Estimate[i])
        result.lst$boundary_estimates <- c(result.lst$boundary_estimates, msg)
        log_fn(paste0(name_str, "NOTE: ", msg))
      }
    }
  }

  return(result.lst)
}


#' Extract comprehensive fit indices from OpenMx model
#'
#' Computes standard SEM fit indices: CFI, TLI, RMSEA, SRMR
#' Requires computing saturated and independence reference models.
#'
#' @param model.fit Fitted OpenMx model
#' @param compute_refmodels Whether to compute CFI/TLI/RMSEA (slower)
#' @return Named list of fit indices
#' @export
extract_openmx_fit_indices.fn <- function(model.fit, compute_refmodels = TRUE) {
  if (is.null(model.fit)) return(NULL)

  summ <- tryCatch(summary(model.fit), error = function(e) NULL)
  if (is.null(summ)) return(NULL)

  # Basic fit indices (always available)
  fit.lst <- list(
    minus2LL = model.fit$output$minimum,
    AIC = tryCatch(AIC(model.fit), error = function(e) NA),
    BIC = tryCatch(BIC(model.fit), error = function(e) NA),
    df = summ$degreesOfFreedom,
    n_params = summ$estimatedParameters,
    n_obs = summ$numObs,
    status = model.fit$output$status$code,
    CFI = NA,
    TLI = NA,
    RMSEA = NA,
    RMSEA_CI_lower = NA,
    RMSEA_CI_upper = NA,
    SRMR = NA,
    timestamp = Sys.time()
  )

  # Compute reference models for CFI/TLI/RMSEA
  if (compute_refmodels) {
    tryCatch({
      # Get reference models (saturated and independence)
      ref_models <- mxRefModels(model.fit, run = TRUE)

      if (!is.null(ref_models)) {
        # Extract fit from reference models
        sat_fit <- ref_models$Saturated
        ind_fit <- ref_models$Independence

        if (!is.null(sat_fit) && !is.null(ind_fit)) {
          # Model chi-square (vs saturated)
          chi_sq_model <- model.fit$output$minimum - sat_fit$output$minimum
          df_model <- summ$degreesOfFreedom

          # Independence model chi-square
          chi_sq_ind <- ind_fit$output$minimum - sat_fit$output$minimum
          ind_summ <- summary(ind_fit)
          df_ind <- ind_summ$degreesOfFreedom

          # CFI = 1 - (chi_model - df_model) / (chi_ind - df_ind)
          # But use max(chi - df, 0) to handle negative values
          ncp_model <- max(chi_sq_model - df_model, 0)
          ncp_ind <- max(chi_sq_ind - df_ind, 0)
          if (ncp_ind > 0) {
            fit.lst$CFI <- 1 - ncp_model / ncp_ind
          }

          # TLI (NNFI) = [(chi_ind/df_ind) - (chi_model/df_model)] /
          #              [(chi_ind/df_ind) - 1]
          if (df_model > 0 && df_ind > 0) {
            ratio_ind <- chi_sq_ind / df_ind
            ratio_model <- chi_sq_model / df_model
            if (ratio_ind > 1) {
              fit.lst$TLI <- (ratio_ind - ratio_model) / (ratio_ind - 1)
            }
          }

          # RMSEA = sqrt(max(chi_model - df_model, 0) / (df_model * (N - 1)))
          n <- summ$numObs
          if (df_model > 0 && n > 1) {
            rmsea_num <- max(chi_sq_model - df_model, 0)
            fit.lst$RMSEA <- sqrt(rmsea_num / (df_model * (n - 1)))

            # RMSEA 90% CI (using non-central chi-square)
            # Lower bound
            if (chi_sq_model > df_model) {
              lambda_lower <- tryCatch({
                uniroot(function(lambda) {
                  pchisq(chi_sq_model, df_model, ncp = lambda) - 0.95
                }, c(0, chi_sq_model * 3))$root
              }, error = function(e) 0)
              fit.lst$RMSEA_CI_lower <- sqrt(lambda_lower / (df_model * (n - 1)))
            } else {
              fit.lst$RMSEA_CI_lower <- 0
            }

            # Upper bound
            lambda_upper <- tryCatch({
              uniroot(function(lambda) {
                pchisq(chi_sq_model, df_model, ncp = lambda) - 0.05
              }, c(0, chi_sq_model * 5 + 100))$root
            }, error = function(e) NA)
            if (!is.na(lambda_upper)) {
              fit.lst$RMSEA_CI_upper <- sqrt(lambda_upper / (df_model * (n - 1)))
            }
          }
        }
      }

      # SRMR - compute from residual correlations
      # This requires the observed and model-implied covariance matrices
      obs_cov <- tryCatch(model.fit$data$observed, error = function(e) NULL)
      if (!is.null(obs_cov) && is.data.frame(obs_cov)) {
        # Get manifest variables
        manifests <- model.fit$manifestVars
        if (length(manifests) > 0) {
          obs_cov_mat <- cov(obs_cov[, manifests, drop = FALSE], use = "pairwise")
          impl_cov <- tryCatch({
            mxGetExpected(model.fit, "covariance")
          }, error = function(e) NULL)

          if (!is.null(impl_cov) && !is.null(obs_cov_mat)) {
            # Standardize both matrices
            obs_sd <- sqrt(diag(obs_cov_mat))
            impl_sd <- sqrt(diag(impl_cov))

            obs_cor <- obs_cov_mat / outer(obs_sd, obs_sd)
            impl_cor <- impl_cov / outer(impl_sd, impl_sd)

            # SRMR = sqrt(mean of squared residual correlations, lower triangle)
            resid_cor <- obs_cor - impl_cor
            lower_tri <- resid_cor[lower.tri(resid_cor, diag = FALSE)]
            fit.lst$SRMR <- sqrt(mean(lower_tri^2, na.rm = TRUE))
          }
        }
      }

    }, error = function(e) {
      # Silently fail - fit indices will remain NA
    })
  }

  return(fit.lst)
}


#' Extract key parameters with z-values and p-values
#'
#' @param model.fit Fitted OpenMx model
#' @param params_of_interest Character vector of parameter name patterns
#' @return data.table with parameter estimates, SEs, and p-values
#' @export
extract_lgcm_params.fn <- function(model.fit,
                                   params_of_interest = NULL) {
  if (is.null(model.fit)) return(NULL)

  summ <- tryCatch(summary(model.fit), error = function(e) NULL)
  if (is.null(summ)) return(NULL)

  params.dt <- as.data.table(summ$parameters)

  # Filter to parameters of interest if specified
  if (!is.null(params_of_interest)) {
    pattern <- paste(params_of_interest, collapse = "|")
    params.dt <- params.dt[grepl(pattern, name)]
  }

  # Add z-values and p-values
  if (nrow(params.dt) > 0 && "Std.Error" %in% names(params.dt)) {
    params.dt[, z := Estimate / Std.Error]
    params.dt[, p := 2 * pnorm(-abs(z))]
    params.dt[, sig := fcase(
      p < 0.001, "***",
      p < 0.01, "**",
      p < 0.05, "*",
      p < 0.10, ".",
      default = ""
    )]
  }

  return(params.dt)
}


#' Print diagnostic summary for a set of models
#'
#' @param diagnostics.lst List of diagnostic results
#' @param log_fn Logging function
#' @export
print_diagnostic_summary.fn <- function(diagnostics.lst, log_fn = message) {
  n_models <- length(diagnostics.lst)
  n_converged <- sum(sapply(diagnostics.lst, function(x) x$converged))
  n_admissible <- sum(sapply(diagnostics.lst, function(x) x$is_admissible))

  all_heywood <- unlist(lapply(diagnostics.lst, function(x) x$heywood_cases))
  all_large_se <- unlist(lapply(diagnostics.lst, function(x) x$large_se))
  all_boundary <- unlist(lapply(diagnostics.lst, function(x) x$boundary_estimates))
  all_warnings <- unlist(lapply(diagnostics.lst, function(x) x$warnings))

  log_fn("")
  log_fn("=== DIAGNOSTIC SUMMARY ===")
  log_fn(sprintf("  Models: %d total, %d converged, %d admissible",
                 n_models, n_converged, n_admissible))

  if (length(all_heywood) > 0) {
    log_fn(sprintf("  Heywood cases: %d", length(all_heywood)))
    for (h in all_heywood) log_fn(sprintf("    - %s", h))
  } else {
    log_fn("  Heywood cases: None")
  }

  if (length(all_large_se) > 0) {
    log_fn(sprintf("  Large SEs: %d", length(all_large_se)))
    for (s in all_large_se) log_fn(sprintf("    - %s", s))
  } else {
    log_fn("  Large SEs: None")
  }

  if (length(all_boundary) > 0) {
    log_fn(sprintf("  Boundary estimates: %d", length(all_boundary)))
  }

  if (length(all_warnings) > 0) {
    log_fn(sprintf("  Other warnings: %d", length(all_warnings)))
    for (w in all_warnings) log_fn(sprintf("    - %s", w))
  }

  log_fn("")

  # Return summary stats
  invisible(list(
    n_models = n_models,
    n_converged = n_converged,
    n_admissible = n_admissible,
    n_heywood = length(all_heywood),
    n_large_se = length(all_large_se),
    n_boundary = length(all_boundary),
    n_warnings = length(all_warnings)
  ))
}
