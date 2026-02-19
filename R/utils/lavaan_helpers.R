#!/usr/bin/env Rscript

# =============================================================================
# lavaan_helpers.R
# =============================================================================
# Utility functions for lavaan model diagnostics and extraction
#
# Functions:
#   - check_convergence(): Check model convergence with diagnostics
#   - extract_fit_indices(): Extract and format fit indices
# =============================================================================

suppressPackageStartupMessages({
  library(lavaan)
})

# =============================================================================
# Model Checking and Diagnostics
# =============================================================================

#' Check model convergence with detailed diagnostics
#'
#' @param fit.fit Fitted lavaan model
#' @param model_name Optional name for logging
#' @return List with convergence status and diagnostics
#' @export
check_convergence <- function(fit.fit, model_name = NULL) {
  name_str <- if (!is.null(model_name)) {
    paste0(" (", model_name, ")")
  } else {
    ""
  }

  result.lst <- list(
    converged = FALSE,
    admissible = FALSE,
    issues = character(0)
  )

  # Check if fit object is valid
  if (is.null(fit.fit)) {
    result.lst$issues <- c(result.lst$issues, "Model fit is NULL")
    return(result.lst)
  }

  # Basic convergence
  result.lst$converged <- lavInspect(fit.fit, "converged")
  if (!result.lst$converged) {
    result.lst$issues <- c(result.lst$issues, "Model did not converge")
  }

  # Check for inadmissible solutions
  tryCatch({
    # Negative variances
    est.dt <- parameterEstimates(fit.fit)
    variances.dt <- est.dt[est.dt$op == "~~" & est.dt$lhs == est.dt$rhs, ]
    neg_var.dt <- variances.dt[variances.dt$est < 0, ]
    if (nrow(neg_var.dt) > 0) {
      result.lst$issues <- c(
        result.lst$issues,
        sprintf("Negative variance(s): %s",
                paste(neg_var.dt$lhs, collapse = ", "))
      )
    }

    # Correlations > 1
    cors.v <- lavInspect(fit.fit, "cor.lv")
    if (is.matrix(cors.v) && any(abs(cors.v[lower.tri(cors.v)]) > 1)) {
      result.lst$issues <- c(result.lst$issues,
                             "Latent variable correlation > 1")
    }

    # Standard errors
    if (any(is.na(est.dt$se[est.dt$free > 0]))) {
      result.lst$issues <- c(result.lst$issues, "Missing standard errors")
    }

    result.lst$admissible <- length(result.lst$issues) == 0 ||
      (length(result.lst$issues) == 1 && !result.lst$converged)

  }, error = function(e) {
    result.lst$issues <- c(result.lst$issues,
                           paste("Diagnostic error:", e$message))
  })

  result.lst
}

# =============================================================================
# Parameter Extraction
# =============================================================================

#' Extract fit indices
#'
#' @param fit.fit Fitted lavaan model
#' @param indices.v Vector of fit index names
#' @return Named vector of fit indices
#' @export
extract_fit_indices <- function(
    fit.fit,
    indices.v = c("chisq", "df", "pvalue", "cfi", "tli",
                  "rmsea", "srmr", "aic", "bic")) {
  if (is.null(fit.fit) || !lavInspect(fit.fit, "converged")) {
    return(setNames(rep(NA, length(indices.v)), indices.v))
  }

  all_indices.v <- fitmeasures(fit.fit)

  # Handle scaled indices (MLR)
  result.v <- sapply(indices.v, function(idx) {
    scaled_name <- paste0(idx, ".scaled")
    robust_name <- paste0(idx, ".robust")

    if (scaled_name %in% names(all_indices.v)) {
      all_indices.v[scaled_name]
    } else if (robust_name %in% names(all_indices.v)) {
      all_indices.v[robust_name]
    } else if (idx %in% names(all_indices.v)) {
      all_indices.v[idx]
    } else {
      NA
    }
  })

  names(result.v) <- indices.v
  result.v
}
