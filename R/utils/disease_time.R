# =============================================================================
# disease_time.R
# =============================================================================
# Estimated Disease Time (EDT) computation using progmod
#
# References:
#   Univariate (pre-computed fallback):
#     Raket LL. Statistical Disease Progression Modeling in Alzheimer Disease.
#     Front Big Data. 2020;3:24. doi:10.3389/fdata.2020.00024
#
#   Multivariate ADAS + MMSE (preferred):
#     Kühnel L, Berger A-K, Markussen B, Raket LL. Simultaneous modeling of
#     Alzheimer's disease progression via multiple cognitive scales.
#     Stat Med. 2021;40:3251-66. doi:10.1002/sim.8932
#
# IMPORTANT: EDT is computed from RAW cognitive tests (ADAS-Cog, MMSE) to
# avoid circularity when predicting harmonized cognitive composites.
#
# Options:
#   1. Fit multivariate ADAS13 + MMSE model (Kühnel et al., 2021)
#   2. Fall back to pre-computed ADNI_disease_stage_bl from progmod
#      (Raket, 2020)
# =============================================================================

#' Compute Estimated Disease Time using progmod
#'
#' @param data data.table with columns: PTID, YRS_from_bl
#' @param config Pipeline configuration (optional, for cache path and refit
#'   flag)
#' @param use_precomputed Logical; use pre-computed ADNI predictions
#'   (default TRUE)
#' @param fit_multivariate Logical; if precomputed fails, fit ADAS+MMSE model
#' @return data.table with EDT column added
#'
#' @details EDT is computed from raw cognitive tests (ADAS-cog) via progmod.
#' The TRUE simultaneous multivariate model (Kühnel et al., 2021) uses both
#' ADAS-Cog 13 and MMSE to inform a SHARED disease time random effect.
#' For follow-up visits, EDT is propagated as: EDT = baseline_EDT + YRS_from_bl
#'
compute_disease_time <- function(data,
                                  config = NULL,
                                  use_precomputed = TRUE,
                                  fit_multivariate = TRUE) {

  if (!requireNamespace("progmod", quietly = TRUE)) {
    log_warn("progmod package not available")
    log_warn("Install with: renv::install('larslau/progmod')")
    return(add_null_edt(data))
  }

  log_info("Computing Estimated Disease Time (EDT) using progmod")

  # Get cache settings from config
  cache_path <- NULL
  refit <- TRUE
  if (!is.null(config)) {
    cache_path <- config$parameters$disease_time$cache_path
    if (!is.null(cache_path)) {
      cache_path <- file.path(here::here(), cache_path)
    }
    refit_setting <- config$parameters$disease_time$refit
    if (!is.null(refit_setting)) {
      refit <- as.logical(refit_setting)
    }
  }

  # Strategy 1: Try TRUE simultaneous multivariate ADAS13 + MMSE model
  # (preferred)
  # (Kühnel et al., Stat Med 2021)
  if (fit_multivariate) {
    result <- fit_multivariate_edt(data, cache_path = cache_path, refit = refit)
    if (!is.null(result)) {
      return(result)
    }
    log_info("Multivariate fit failed; falling back to pre-computed EDT")
  }

  # Strategy 2: Fall back to pre-computed ADNI baseline EDT from progmod package
  # Based on univariate ADAS-cog 13 trajectories (Raket, 2020)
  if (use_precomputed) {
    result <- try_precomputed_edt(data)
    if (!is.null(result)) {
      return(result)
    }
  }

  log_warn("EDT computation failed - returning NULL EDT")
  return(add_null_edt(data))
}

#' Use pre-computed ADNI EDT from progmod package
#'
#' @param data data.table with PTID and YRS_from_bl columns
#' @return data.table with EDT added, or NULL if insufficient coverage
#'
try_precomputed_edt <- function(data) {

  tryCatch({
    # Load Raket's pre-computed ADNI disease stage (from progmod package)
    # This is based on longitudinal ADAS-cog 13 trajectories
    library(progmod)
    utils::data(
      "ADNI_disease_stage_bl", package = "progmod", envir = environment()
    )

    edt_bl.dt <- data.table::as.data.table(ADNI_disease_stage_bl)
    log_info(
      "Loaded progmod ADNI_disease_stage_bl: %d subjects", nrow(edt_bl.dt)
    )
    log_info("  EDT range: %.1f to %.1f months (%.1f to %.1f years)",
             min(edt_bl.dt$pred_AD_month), max(edt_bl.dt$pred_AD_month),
             min(edt_bl.dt$pred_AD_month)/12, max(edt_bl.dt$pred_AD_month)/12)

    # Extract RID from PTID (format: XXX_S_XXXX, RID is last 4 digits)
    data[, RID := as.integer(sub(".*_S_", "", PTID))]

    # Merge baseline EDT
    data <- merge(data, edt_bl.dt, by = "RID", all.x = TRUE)

    # Convert to years and compute visit-level EDT
    # EDT at visit = baseline EDT + time since baseline
    data[, EDT_bl_years := pred_AD_month / 12]
    data[, EDT := EDT_bl_years + YRS_from_bl]

    # Clean up
    data[, c("pred_AD_month", "RID") := NULL]

    # Check coverage
    n_with_edt <- sum(!is.na(data$EDT))
    n_total <- nrow(data)
    pct_coverage <- 100 * n_with_edt / n_total

    log_info("Pre-computed EDT coverage: %d/%d observations (%.1f%%)",
             n_with_edt, n_total, pct_coverage)
    log_info("  Subjects with EDT: %d/%d",
             data[!is.na(EDT), uniqueN(PTID)],
             data[, uniqueN(PTID)])

    if (pct_coverage < 50) {
      log_warn(
        "Low EDT coverage (%.1f%%) - will try multivariate fit", pct_coverage
      )
      data[, c("EDT", "EDT_bl_years") := NULL]
      return(NULL)
    }

    attr(data, "edt_available") <- TRUE
    attr(data, "edt_method") <- "progmod_precomputed"
    attr(data, "edt_coverage") <- pct_coverage

    return(data)

  }, error = function(e) {
    log_warn("Pre-computed EDT failed: %s", e$message)
    return(NULL)
  })
}

#' Fit multivariate ADAS13 + MMSE progmod model
#'
#' Uses progmod to fit a simultaneous multivariate model where both
#' ADAS-Cog 13 and MMSE inform a shared disease time random effect.
#' This implements the approach from Kühnel et al. (Stat Med 2021).
#'
#' Key features:
#' - Long-format stacked data with scale factor
#' - Shared random effect for disease time shift (s)
#' - Heteroscedastic residuals via varIdent weights
#' - Scale-specific growth parameters (l, g, v)
#'
#' @param data data.table with PTID and YRS_from_bl columns
#' @param cache_path Path to save/load cached model fit (optional)
#' @param refit If TRUE, refit model even if cache exists (default TRUE)
#' @return data.table with EDT added, or NULL if fitting fails
#'
fit_multivariate_edt <- function(data, cache_path = NULL, refit = TRUE) {

  tryCatch({
    if (!requireNamespace("ADNIMERGE", quietly = TRUE)) {
      log_warn("ADNIMERGE package required for multivariate EDT fitting")
      return(NULL)
    }

    library(progmod)
    library(ADNIMERGE)

    log_info("Fitting multivariate ADAS13 + MMSE progmod model")
    log_info("  (Kühnel et al., Stat Med 2021)")

    # Load ADNIMERGE data
    utils::data("adnimerge", package = "ADNIMERGE", envir = environment())
    adni <- data.table::as.data.table(adnimerge)

    # Extract RIDs from our data (for later filtering)
    data[, RID := as.integer(sub(".*_S_", "", PTID))]
    our_rids.v <- unique(data$RID)

    # -------------------------------------------------------------------------
    # Prepare wide-format data for ALL ADNIMERGE subjects
    # Fit model on full dataset to ensure all subjects have random effects
    # -------------------------------------------------------------------------
    pm_wide.df <- adni[!is.na(DX.bl),
                    .(id = RID,
                      time = as.numeric(as.character(Month)),
                      ADAS13 = ADAS13,
                      MMSE = MMSE,
                      DX = as.character(DX.bl))]
    pm_wide.df <- as.data.frame(pm_wide.df)

    # Create diagnosis indicators (MCI includes LMCI, EMCI, SMC)
    pm_wide.df$MCI <- as.numeric(pm_wide.df$DX %in% c("LMCI", "EMCI", "SMC"))
    pm_wide.df$DEM <- as.numeric(pm_wide.df$DX == "AD")

    n_subj_all <- length(unique(pm_wide.df$id))
    log_info("Full ADNIMERGE data prepared: %d subjects", n_subj_all)
    log_info("  Analysis cohort: %d subjects", length(our_rids.v))

    # -------------------------------------------------------------------------
    # Check for cached model fit
    # -------------------------------------------------------------------------
    multi_fit.fit <- NULL
    if (!is.null(cache_path) && file.exists(cache_path) && !refit) {
      log_info(
        "Loading cached multivariate progmod fit from: %s", basename(cache_path)
      )
      cached <- tryCatch({
        readRDS(cache_path)
      }, error = function(e) {
        log_warn("Failed to load cached fit: %s", e$message)
        NULL
      })

      # Verify cache was fit on full ADNIMERGE (not a subset)
      if (!is.null(cached)) {
        re <- nlme::ranef(cached)
        if (nrow(re) >= n_subj_all * 0.9) {  # Allow 10% tolerance
          multi_fit.fit <- cached
          log_info("  Cached fit contains %d subjects", nrow(re))
        } else {
          log_info("  Cached fit has only %d subjects (need ~%d), refitting...",
                   nrow(re), n_subj_all)
        }
      }
    }

    # -------------------------------------------------------------------------
    # Fit model if not loaded from cache
    # -------------------------------------------------------------------------
    if (is.null(multi_fit.fit)) {
      log_info("Stacking to long format for simultaneous fitting...")

      tmp_adas.df <- pm_wide.df[!is.na(pm_wide.df$ADAS13),
                          c("id", "time", "MCI", "DEM", "ADAS13")]
      names(tmp_adas.df)[5] <- "value"
      tmp_adas.df$scale <- "ADAS13"

      tmp_mmse.df <- pm_wide.df[!is.na(pm_wide.df$MMSE),
                          c("id", "time", "MCI", "DEM", "MMSE")]
      names(tmp_mmse.df)[5] <- "value"
      tmp_mmse.df$scale <- "MMSE"

      long_data.df <- rbind(tmp_adas.df, tmp_mmse.df)
      long_data.df$scale <- factor(long_data.df$scale)

      n_obs <- nrow(long_data.df)
      n_adas <- sum(long_data.df$scale == "ADAS13")
      n_mmse <- sum(long_data.df$scale == "MMSE")
      log_info(
        "  Long-format data: %d subjects, %d observations", n_subj_all, n_obs
      )
      log_info("  ADAS13: %d obs, MMSE: %d obs", n_adas, n_mmse)

      log_info("Fitting multivariate progmod model on full ADNIMERGE...")

      # Starting values (empirically determined from sample fits)
      fixed_start_coef.v <- c(
        l.scaleADAS13 = 0.24,
        l.scaleMMSE = -0.02,
        s.MCI = 82,
        s.DEM = 158,
        g.scaleADAS13 = 3.66,
        g.scaleMMSE = 3.43,
        v.scaleADAS13 = 10.56,
        v.scaleMMSE = 28.65
      )

      multi_fit.fit <- tryCatch({
        progmod(
          value ~ exp_model(time, l, s, g, v),
          data = long_data.df,
          fixed = list(
            l ~ scale + 0,      # Scale-specific lower asymptote
            # Shared shift (diagnosis-specific, CN = reference)
            s ~ MCI + DEM + 0,
            g ~ scale + 0,      # Scale-specific growth rate
            v ~ scale + 0       # Scale-specific upper asymptote
          ),
          random = list(
            s ~ 1,              # Shared random shift
            v ~ scale           # Scale-specific random intercept
          ),
          groups = ~ id,
          start = fixed_start_coef.v,
          weights = nlme::varIdent(form = ~ 1 | scale),
          control = nlme::nlmeControl(
            maxIter = 200,
            pnlsTol = 0.01,
            returnObject = TRUE
          )
        )
      }, error = function(e) {
        log_warn("Multivariate progmod failed: %s", e$message)
        NULL
      })

      # Save to cache if successful
      if (!is.null(multi_fit.fit) && !is.null(cache_path)) {
        log_info("Saving multivariate fit to cache: %s", basename(cache_path))
        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(multi_fit.fit, cache_path)
      }
    }

    if (is.null(multi_fit.fit)) {
      log_warn("Multivariate model failed - no fallback available")
      return(NULL)
    }

    # -------------------------------------------------------------------------
    # Extract fixed and random effects
    # -------------------------------------------------------------------------
    fe <- nlme::fixef(multi_fit.fit)
    re <- nlme::ranef(multi_fit.fit)

    log_info("Model converged successfully!")
    log_info("  Fixed effects: s.MCI=%.1f, s.DEM=%.1f months",
             fe["s.MCI"], fe["s.DEM"])

    # -------------------------------------------------------------------------
    # Compute individual EDT estimates
    # EDT = population shift (diagnosis-specific) + individual random shift
    # -------------------------------------------------------------------------
    log_info("Computing individual EDT estimates...")

    # Get unique subject info
    dx_bl.df <- unique(pm_wide.df[, c("id", "DX", "MCI", "DEM")])
    dx_bl.df <- dx_bl.df[!duplicated(dx_bl.df$id), ]

    # Population shift based on diagnosis
    dx_bl.df$pop_shift <- ifelse(
      dx_bl.df$DEM == 1, fe["s.DEM"],
      ifelse(dx_bl.df$MCI == 1, fe["s.MCI"], 0)
    )

    # Individual random shift (SHARED across both scales!)
    # The random effect for 's' is stored as "s.(Intercept)" in the ranef output
    dx_bl.df$ind_shift <- re[as.character(dx_bl.df$id), "s.(Intercept)"]

    # Total EDT = population + individual shift (in months)
    dx_bl.df$EDT_bl_months <- dx_bl.df$pop_shift + dx_bl.df$ind_shift

    n_with_edt <- sum(!is.na(dx_bl.df$EDT_bl_months))
    log_info("  Subjects with EDT: %d/%d", n_with_edt, nrow(dx_bl.df))

    # Create EDT estimates table
    edt_estimates.dt <- data.table::data.table(
      RID = dx_bl.df$id,
      EDT_bl_months = dx_bl.df$EDT_bl_months
    )

    # -------------------------------------------------------------------------
    # Merge EDT estimates with input data
    # -------------------------------------------------------------------------
    data <- merge(data, edt_estimates.dt, by = "RID", all.x = TRUE)

    # Convert to years and add time from baseline
    data[, EDT_bl_years := EDT_bl_months / 12]
    data[, EDT := EDT_bl_years + YRS_from_bl]

    # Clean up (keep RID for downstream use)
    data[, c("EDT_bl_months") := NULL]

    # Check coverage
    n_with_edt <- sum(!is.na(data$EDT))
    pct_coverage <- 100 * n_with_edt / nrow(data)

    log_info("Multivariate EDT coverage: %d/%d (%.1f%%)",
             n_with_edt, nrow(data), pct_coverage)

    # Report EDT range
    edt_range.v <- range(data$EDT, na.rm = TRUE)
    log_info("  EDT range: %.1f to %.1f years", edt_range.v[1], edt_range.v[2])

    attr(data, "edt_available") <- TRUE
    attr(data, "edt_method") <- "progmod_simultaneous_adas_mmse"
    attr(data, "edt_coverage") <- pct_coverage
    attr(data, "edt_reference") <- "Kühnel et al., Stat Med 2021"

    return(data)

  }, error = function(e) {
    log_warn("Multivariate EDT fitting failed: %s", e$message)
    return(NULL)
  })
}

#' Add NULL EDT column (fallback)
#'
#' @param data data.table
#' @return data.table with EDT = NA
#'
add_null_edt <- function(data) {
  data[, EDT := NA_real_]
  data[, EDT_bl_years := NA_real_]
  attr(data, "edt_available") <- FALSE
  attr(data, "edt_method") <- "none"
  return(data)
}

#' Validate EDT estimates
#'
#' @param data data.table with EDT column
#' @return Logical indicating if EDT is valid
#'
validate_edt <- function(data) {
  if (!"EDT" %in% names(data)) return(FALSE)
  if (all(is.na(data$EDT))) return(FALSE)

  # Check for reasonable range (e.g., -10 to +30 years from reference)
  edt_range.v <- range(data$EDT, na.rm = TRUE)
  if (edt_range.v[2] - edt_range.v[1] > 50) {
    log_warn("EDT range suspiciously wide: [%.1f, %.1f] years",
             edt_range.v[1], edt_range.v[2])
    return(FALSE)
  }

  # Check coverage
  pct_valid <- 100 * mean(!is.na(data$EDT))
  if (pct_valid < 50) {
    log_warn("EDT coverage too low: %.1f%%", pct_valid)
    return(FALSE)
  }

  return(TRUE)
}

#' Hybrid EDT extrapolation for LGCM
#'
#' Extrapolates missing EDT values using a hybrid strategy:
#' 1. If subject has >= 2 observed EDT values: use individual linear fit
#' 2. If subject has 1 observed EDT value: use population-average rate
#'
#' @param dt data.table with EDT columns (EDT_T0, EDT_T1, etc.)
#' @param pop_rate Population-average EDT change per visit (computed if NULL)
#' @param valid_timepoints Character vector of timepoint labels
#'   (e.g., c("T0", "T1", ...)). If NULL, uses get_valid_timepoints() from
#'   config
#' @return List with extrapolated data and metadata:
#'   - data: data.table with extrapolated EDT values
#'   - pop_rate: population-average rate used
#'   - n_individual_extrap: count of individual extrapolations
#'   - n_population_extrap: count of population-rate extrapolations
#'
#' @export
extrapolate_edt_hybrid <- function(dt,
                                   pop_rate = NULL,
                                   valid_timepoints = NULL) {
  dt <- data.table::copy(dt)

  # Get timepoints from config if not provided
  if (is.null(valid_timepoints)) {
    valid_timepoints <- get_valid_timepoints()
  }
  n_timepoints <- length(valid_timepoints)

  edt_cols.v <- paste0("EDT_", valid_timepoints)
  visit_idx.v <- 0:(n_timepoints - 1)

  # Track extrapolation for sensitivity analysis
  n_individual_extrap <- 0
  n_population_extrap <- 0

  # Step 1: Compute population-average EDT rate if not provided
  if (is.null(pop_rate)) {
    # Collect all (visit_index, EDT) pairs across all subjects
    all_pairs.dt <- data.table::data.table(
      visit = integer(), edt = numeric(), ptid = character()
    )
    for (i in seq_len(nrow(dt))) {
      edt_vals.v <- as.numeric(dt[i, ..edt_cols.v])
      valid <- !is.na(edt_vals.v)
      if (sum(valid) >= 1) {
        all_pairs.dt <- rbind(all_pairs.dt, data.table::data.table(
          visit = visit_idx.v[valid],
          edt = edt_vals.v[valid],
          ptid = dt$PTID[i]
        ))
      }
    }

    # Fit linear model: EDT ~ visit
    pop_fit.fit <- lm(edt ~ visit, data = all_pairs.dt)
    pop_rate <- coef(pop_fit.fit)["visit"]
    log_info("  Population-average EDT rate: %.4f per visit", pop_rate)
  }

  # Step 2: Extrapolate for each subject
  for (i in seq_len(nrow(dt))) {
    edt_vals.v <- as.numeric(dt[i, ..edt_cols.v])
    valid_idx.v <- which(!is.na(edt_vals.v))
    missing_idx.v <- which(is.na(edt_vals.v))

    if (length(missing_idx.v) == 0) next  # No missing values

    if (length(valid_idx.v) >= 2) {
      # Individual extrapolation: fit EDT ~ visit_index
      subj_fit.fit <- lm(edt_vals.v[valid_idx.v] ~ valid_idx.v)
      intercept <- coef(subj_fit.fit)[1]
      slope <- coef(subj_fit.fit)[2]

      # Extrapolate to missing timepoints
      for (j in missing_idx.v) {
        # j-1 because visit_idx.v is 0-based
        edt_vals.v[j] <- intercept + slope * (j - 1)
        n_individual_extrap <- n_individual_extrap + 1
      }

    } else if (length(valid_idx.v) == 1) {
      # Population-rate extrapolation
      base_visit <- valid_idx.v - 1  # Convert to 0-based
      base_edt <- edt_vals.v[valid_idx.v]

      for (j in missing_idx.v) {
        target_visit <- j - 1  # Convert to 0-based
        edt_vals.v[j] <- base_edt + pop_rate * (target_visit - base_visit)
        n_population_extrap <- n_population_extrap + 1
      }
    }

    # Update data table
    for (j in seq_along(edt_cols.v)) {
      data.table::set(dt, i, edt_cols.v[j], edt_vals.v[j])
    }
  }

  list(
    data = dt,
    pop_rate = pop_rate,
    n_individual_extrap = n_individual_extrap,
    n_population_extrap = n_population_extrap
  )
}
