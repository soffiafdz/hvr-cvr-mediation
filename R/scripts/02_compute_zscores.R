#!/usr/bin/env Rscript

# =============================================================================
# 01_compute_zscores.R
# =============================================================================
# Calculate deviation Z-scores for ADNI subjects using normative GAMLSS models
# fit on UK Biobank data.
#
# CRITICAL: GAMLSS prediction requires the original training data to be
# accessible in the environment with the same structure as when models were fit.
#
# Inputs (paths from config/pipeline_config.yaml):
#   - raw/adni_age: ADNI age data
#   - raw/adni_hchvr: ADNI head-size adjusted volumes
#   - external/gamlss_models: GAMLSS normative models (from UKB)
#
# Outputs:
#   - derivatives/z_scores: Z-scores for ADNI subjects
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(progress)
  library(gamlss)
  library(gamlss.add)
})

# Source utilities
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# Initialize
log_script_start("01_compute_zscores.R")
config <- load_config()
set_seed()

# --- Configuration ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "compute_zscores",
  default = FALSE
)

output.path <- get_data_path("derivatives", "z_scores")

# Check if output exists and skip if not forcing
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output file exists and force_regenerate=FALSE: %s", output.path)
  log_info("Skipping z-score calculation. Set force_regenerate=TRUE to rerun.")
  log_script_end("01_compute_zscores.R", success = TRUE)
  quit(status = 0)
}

# --- Constants ---
log_section("Setting Up Constants")

# ROI mapping: ADNI uses VC (Ventricles), UKB models use LV (Lateral Ventricles)
ROI_MAP <- c(
  HC  = "HC",
  VC  = "LV", # Map ADNI VC to model LV
  HVR = "HVR"
)

ADJS <- get_parameter("adjustment_methods")
log_info("Head-size adjustment methods: %s", paste(ADJS, collapse = ", "))

# --- Load Input Data ---
log_section("Loading Input Data")

# ADNI age data (from legacy preprocessing)
adni_age.path <- get_data_path("derivatives", "adni_age")
check_files_exist(adni_age.path)
adni_age.dt <- read_rds_safe(adni_age.path, "ADNI age data")
# Strip "labelled" class from AGE (haven package artifact) for GAMLSS compat
adni_age.dt[, AGE := as.numeric(AGE)]
log_info("ADNI age data: %d observations", nrow(adni_age.dt))

# ADNI HC-HVR adjusted volumes (from script 01_preprocess_volumes.R)
adni_hchvr.path <- get_data_path("derivatives", "adni_hchvr")
check_files_exist(adni_hchvr.path)
adni_hchvr.dt <- read_rds_safe(adni_hchvr.path, "ADNI HC-HVR adjusted volumes")
log_info("ADNI HC-HVR data: %d observations", nrow(adni_hchvr.dt))

# GAMLSS models
gamlss.path <- get_data_path("external", "gamlss_models")
check_files_exist(gamlss.path)
gamlss.lst <- read_rds_safe(gamlss.path, "GAMLSS models")
log_info("GAMLSS models loaded successfully")

# --- Load ADNIMERGE for Demographics ---
log_section("Loading ADNIMERGE Demographics")

# Check for ADNIMERGE
if (!requireNamespace("ADNIMERGE", quietly = TRUE)) {
  adnimerge_tar <- list.files(
    here("data"),
    pattern = "ADNIMERGE.*\\.tar\\.gz$", full.names = TRUE
  )
  if (length(adnimerge_tar) > 0) {
    log_info("Installing ADNIMERGE from local tar: %s", adnimerge_tar[1])
    install.packages(adnimerge_tar[1], repos = NULL, type = "source")
  } else {
    stop("ADNIMERGE package not found. Please install it.", call. = FALSE)
  }
}

library(ADNIMERGE)
data(adnimerge)
setDT(adnimerge)

adni_demog.dt <- adnimerge[, .(
  PTID,
  SEX = as.character(PTGENDER),
  EDUC_num = as.numeric(PTEDUCAT) # Already in years for ADNI
)] |>
  unique() |>
  setkey(PTID)

log_info("ADNIMERGE demographics: %d subjects", nrow(adni_demog.dt))
log_info(
  "Sex distribution: %d Female, %d Male",
  adni_demog.dt[SEX == "Female", .N],
  adni_demog.dt[SEX == "Male", .N]
)

rm(adnimerge)

# --- Data Preparation ---
log_section("Preparing ADNI Data")

# Merge age and demographics
adni_base.dt <- merge(
  adni_age.dt[, -"VISCODE"], # Remove duplicate VISCODE
  adni_demog.dt,
  by = "PTID",
  all.x = TRUE
)
setkey(adni_base.dt, PTID, EXAMDATE)
log_info("Merged base data: %d observations", nrow(adni_base.dt))

# Check for missing demographics
n_missing_demog <- adni_base.dt[is.na(SEX) | is.na(EDUC_num), .N]
if (n_missing_demog > 0) {
  log_warn("Missing demographics for %d observations", n_missing_demog)
  adni_base.dt <- adni_base.dt[!is.na(SEX) & !is.na(EDUC_num)]
}

# Get unique base info (PTID, EXAMDATE, ICC)
unique_cols.v <- c("PTID", "EXAMDATE", "VISCODE", "DX", "ICC")
adni_unique.dt <- adni_hchvr.dt[, .SD, .SDcols = unique_cols.v] |>
  unique() |>
  setkey(PTID, EXAMDATE)

# Merge with demographics
adni_merged.dt <- merge(
  adni_unique.dt,
  adni_base.dt,
  by = c("PTID", "EXAMDATE"),
  all.x = TRUE
)

# Create long format with ROI and VAL columns
# For HC: use NON (ICC included in GAMLSS model) for sensitivity analysis
# For VC: use RES (residuals) - head-size adjusted via regression
# For HVR: use NON adjustment (ratio inherently normalized)
adni_long.dt <- rbind(
  # HC - NON adjustment (ICC in GAMLSS model, for sensitivity analysis)
  adni_hchvr.dt["NON", on = "ADJ"][
    , .(PTID, EXAMDATE, SIDE, ADJ = "NON", ROI = "HC", VAL = HC)
  ],
  # VC (mapped to LV) - residuals adjusted
  adni_hchvr.dt["RES", on = "ADJ"][
    , .(PTID, EXAMDATE, SIDE, ADJ, ROI = "VC", VAL = VC)
  ],
  # HVR - use NON adjustment (GAMLSS models only have HVR with NON)
  # HVR is inherently head-size normalized as a ratio
  adni_hchvr.dt["NON", on = "ADJ"][
    , .(PTID, EXAMDATE, SIDE, ADJ = "NON", ROI = "HVR", VAL = HVR)
  ],
  use.names = TRUE
)

# Add bilateral values using SUM approach (consistent with UKB GAMLSS training)
# CRITICAL: GAMLSS models were trained on SUMMED bilateral volumes, not averages
# HVR bilateral: (HC_L + HC_R) / ((HC_L + HC_R) + (LV_L + LV_R))
# See: headsize_sexeffects_ukb/R/05_adjust_headsize.R:134
adni_bilateral.dt <- rbind(
  # HC bilateral - sum of left + right (NON for sensitivity analysis)
  adni_hchvr.dt["NON", on = "ADJ"][
    , .(VAL = sum(HC, na.rm = TRUE), SIDE = "LR", ADJ = "NON", ROI = "HC"),
    by = .(PTID, EXAMDATE)
  ],
  # VC bilateral - sum of left + right (residuals adjusted)
  adni_hchvr.dt["RES", on = "ADJ"][
    , .(VAL = sum(VC, na.rm = TRUE), SIDE = "LR", ADJ = "RES", ROI = "VC"),
    by = .(PTID, EXAMDATE)
  ],
  # HVR bilateral - computed from summed raw volumes as HC / (HC + VC)
  # Uses NON adjustment to match GAMLSS model structure
  adni_hchvr.dt["NON", on = "ADJ"][
    , .(
      HC_total = sum(HC, na.rm = TRUE),
      VC_total = sum(VC, na.rm = TRUE)
    ),
    by = .(PTID, EXAMDATE)
  ][, .(
    VAL = HC_total / (HC_total + VC_total),
    SIDE = "LR", ADJ = "NON", ROI = "HVR"
  ),
  by = .(PTID, EXAMDATE)
  ],
  use.names = TRUE
)

adni_long.dt <- rbind(adni_long.dt, adni_bilateral.dt, use.names = TRUE)

# Merge with demographics
adni_data.dt <- merge(
  adni_long.dt,
  adni_merged.dt,
  by = c("PTID", "EXAMDATE"),
  all.x = TRUE
) |>
  na.omit()

setkey(adni_data.dt, ROI, ADJ, SIDE, SEX)
log_info("Final ADNI data: %d observations", nrow(adni_data.dt))

# --- Calculate Z-scores ---
log_section("Calculating Z-scores from GAMLSS Models")

# Create sorting table for iteration
sort.dt <- adni_data.dt[, .(ROI, ADJ, SIDE, SEX)] |> unique()
setkey(sort.dt)
log_info("Processing %d ROI/ADJ/SIDE/SEX combinations", nrow(sort.dt))

# Initialize output
adni_data.dt[, ZSCORE := NA_real_]

pb <- progress_bar$new(
  format = "Z-scores | :what [:bar] :current/:total\n",
  total = nrow(sort.dt), clear = FALSE, width = 100, show_after = 0
)
pb$tick(0)

# Track diagnostics
zscore_diag.lst <- list()

for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  roi_adni <- row.dt$ROI # ADNI ROI name (HC, VC, HVR)
  adj <- row.dt$ADJ
  side_adni <- row.dt$SIDE # ADNI SIDE (L, R, LR)
  sex <- row.dt$SEX

  pb$tick(tokens = list(what = sprintf(
    "%s - %s (%s; %s)", sex, roi_adni, side_adni, adj
  )))

  # Map ADNI ROI to model ROI
  # List indexing returns NULL if key not found (not NA)
  roi_model <- ROI_MAP[[roi_adni]]
  if (is.null(roi_model)) {
    log_warn("No model mapping for ROI: %s", roi_adni)
    next
  }

  # Map ADNI SIDE to model SIDE
  # Bilateral data now uses "LR" directly (matching GAMLSS model naming)
  # Keep mapping for backwards compatibility if "AVG" appears
  side_model <- if (side_adni == "AVG") "LR" else side_adni

  # Get final model
  final_mod.lst <- gamlss.lst$FINAL[[sex]][[roi_model]][[adj]][[side_model]]
  if (is.null(final_mod.lst)) {
    log_debug("No model for: %s %s %s %s", sex, roi_model, adj, side_model)
    next
  }

  fam <- final_mod.lst$FAMILY
  mod.fit <- final_mod.lst$FIT

  # Get ADNI data subset
  data.subdt <- adni_data.dt[row.dt][VAL > 0 & !is.na(VAL)]
  if (!nrow(data.subdt)) {
    log_debug("No data for: %s %s %s %s", sex, roi_adni, side_adni, adj)
    next
  }

  # --- Diagnostic: Check for potential extrapolation ---
  # CRITICAL: GAMLSS predict() looks for the original data object by name.
  # The models were fit with `data = gamlss_final.subdt`, so we must create

  # a variable with this exact name containing the training data.
  gamlss_final.subdt <- final_mod.lst$DATA
  train.dt <- gamlss_final.subdt

  train_age_range.v <- range(train.dt$AGE)
  data_age_range.v <- range(data.subdt$AGE)
  train_icc_range.v <- range(train.dt$ICC)
  data_icc_range.v <- range(data.subdt$ICC)

  zscore_age_extrap <- data_age_range.v[1] < train_age_range.v[1] ||
    data_age_range.v[2] > train_age_range.v[2]
  zscore_icc_extrap <- data_icc_range.v[1] < train_icc_range.v[1] ||
    data_icc_range.v[2] > train_icc_range.v[2]

  if (zscore_age_extrap || zscore_icc_extrap) {
    log_debug(
      paste(
        "%s %s (%s; %s): Z-score data may extrapolate",
        "(Age: [%.1f-%.1f] vs [%.1f-%.1f],",
        "ICC: [%.0f-%.0f] vs [%.0f-%.0f])"
      ),
      sex, roi_adni, side_adni, adj,
      data_age_range.v[1], data_age_range.v[2],
      train_age_range.v[1], train_age_range.v[2],
      data_icc_range.v[1], data_icc_range.v[2],
      train_icc_range.v[1], train_icc_range.v[2]
    )
  }

  # Calculate Z-scores using fitted model
  # CRITICAL: Use .SD approach with `new` parameter (not `newdata`) to ensure

  # GAMLSS can find the training data object in the environment.
  # Transfer models (no SITE) use: AGE, EDUC_num, ICC
  pred_cols.v <- c("AGE", "EDUC_num", "ICC")

  zscore_refit_warning <- FALSE
  withCallingHandlers(
    {
      if (fam == "NO") {
        # Gaussian distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pNO(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols.v, "VAL")
        ]
      } else if (fam == "L_NO") {
        # Logit-Normal distribution
        eps <- 1e-6
        data.subdt[, VAL_logit := qlogis(pmin(pmax(VAL, eps), 1 - eps))]
        data.subdt[
          ,
          ZSCORE := qnorm(
            pNO(
              VAL_logit,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols.v, "VAL_logit")
        ]
        data.subdt[, VAL_logit := NULL]
      } else if (fam == "BE") {
        # Beta distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pBE(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols.v, "VAL")
        ]
      } else if (fam == "BCCG") {
        # Box-Cox Cole & Green distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pBCCG(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response"),
              nu = predict(mod.fit, "nu", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols.v, "VAL")
        ]
      } else if (fam == "BCT") {
        # Box-Cox t distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pBCT(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response"),
              nu = predict(mod.fit, "nu", new = .SD, type = "response"),
              tau = predict(mod.fit, "tau", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols.v, "VAL")
        ]
      } else {
        # Unknown distribution family
        log_warn(
          "Unknown GAMLSS family '%s' for %s %s (%s; %s)",
          fam, sex, roi_adni, side_adni, adj
        )
        data.subdt[, ZSCORE := NA_real_]
      }
    },
    warning = function(w) {
      if (grepl("discrepancy.*re-fit", w$message)) {
        zscore_refit_warning <<- TRUE
        log_debug(
          "%s %s (%s; %s): GAMLSS refit warning (z-scores)",
          sex, roi_adni, side_adni, adj
        )
      }
      invokeRestart("muffleWarning")
    }
  )

  # Store diagnostics
  zscore_diag.lst[[sex]][[roi_adni]][[adj]][[side_adni]] <- list(
    REFIT_WARNING = zscore_refit_warning,
    AGE_EXTRAP = zscore_age_extrap,
    ICC_EXTRAP = zscore_icc_extrap,
    N = nrow(data.subdt),
    MODEL_FAMILY = fam
  )

  # Assign Z-scores back to main data
  adni_data.dt[
    data.subdt,
    on = c("PTID", "EXAMDATE", "SEX", "SIDE", "ADJ", "ROI"),
    ZSCORE := i.ZSCORE
  ]
}

# --- Diagnostic Summary ---
log_section("Z-score Calculation Diagnostics")

# Extract diagnostics
zscore_diag.dt <- sort.dt[,
  {
    diag.lst <- zscore_diag.lst[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(diag.lst)) {
      NULL
    } else {
      list(
        REFIT_WARNING = diag.lst$REFIT_WARNING,
        AGE_EXTRAP = diag.lst$AGE_EXTRAP,
        ICC_EXTRAP = diag.lst$ICC_EXTRAP,
        N = diag.lst$N,
        MODEL_FAMILY = diag.lst$MODEL_FAMILY
      )
    }
  },
  by = .(SEX, ROI, ADJ, SIDE)
]

# Summary statistics
if (nrow(zscore_diag.dt) > 0) {
  n_zscore_models <- nrow(zscore_diag.dt)
  n_zscore_refit <- sum(zscore_diag.dt$REFIT_WARNING, na.rm = TRUE)
  n_zscore_age_extrap <- sum(zscore_diag.dt$AGE_EXTRAP, na.rm = TRUE)
  n_zscore_icc_extrap <- sum(zscore_diag.dt$ICC_EXTRAP, na.rm = TRUE)
  total_observations <- sum(zscore_diag.dt$N, na.rm = TRUE)

  log_info("Total z-score calculations: %d model combinations", n_zscore_models)
  log_info("Total observations scored: %d", total_observations)
  log_info(
    "GAMLSS refit warnings: %d (%.1f%%)",
    n_zscore_refit, 100 * n_zscore_refit / n_zscore_models
  )
  log_info(
    "Potential extrapolation: Age=%d (%.1f%%), ICC=%d (%.1f%%)",
    n_zscore_age_extrap, 100 * n_zscore_age_extrap / n_zscore_models,
    n_zscore_icc_extrap, 100 * n_zscore_icc_extrap / n_zscore_models
  )
}

# Z-score summary
n_valid_zscores <- adni_data.dt[!is.na(ZSCORE), .N]
n_total <- nrow(adni_data.dt)
log_info(
  "Valid Z-scores: %d / %d (%.1f%%)",
  n_valid_zscores, n_total, 100 * n_valid_zscores / n_total
)

if (n_valid_zscores > 0) {
  log_info(
    "Z-score range: [%.2f, %.2f], mean: %.2f, sd: %.2f",
    adni_data.dt[!is.na(ZSCORE), min(ZSCORE)],
    adni_data.dt[!is.na(ZSCORE), max(ZSCORE)],
    adni_data.dt[!is.na(ZSCORE), mean(ZSCORE)],
    adni_data.dt[!is.na(ZSCORE), sd(ZSCORE)]
  )
}

# --- Add Extrapolation Flags to Data ---
log_section("Adding Extrapolation Flags to Output")

# Merge diagnostic flags into main data for downstream sensitivity analyses
# This allows filtering out potentially extrapolated z-scores
if (nrow(zscore_diag.dt) > 0) {
  # Create flag columns from diagnostics (group-level flags)
  flag_cols.dt <- zscore_diag.dt[, .(
    SEX, ROI, ADJ, SIDE,
    GAMLSS_REFIT = REFIT_WARNING,
    AGE_EXTRAPOLATED = AGE_EXTRAP,
    ICC_EXTRAPOLATED = ICC_EXTRAP
  )]

  # Merge flags into main data
  adni_data.dt <- merge(
    adni_data.dt, flag_cols.dt,
    by = c("SEX", "ROI", "ADJ", "SIDE"),
    all.x = TRUE
  )

  # Set NA flags to FALSE
  adni_data.dt[is.na(GAMLSS_REFIT), GAMLSS_REFIT := FALSE]
  adni_data.dt[is.na(AGE_EXTRAPOLATED), AGE_EXTRAPOLATED := FALSE]
  adni_data.dt[is.na(ICC_EXTRAPOLATED), ICC_EXTRAPOLATED := FALSE]

  # Summary of flagged observations
  n_refit <- adni_data.dt[GAMLSS_REFIT == TRUE & !is.na(ZSCORE), .N]
  n_age_extrap <- adni_data.dt[AGE_EXTRAPOLATED == TRUE & !is.na(ZSCORE), .N]
  n_icc_extrap <- adni_data.dt[ICC_EXTRAPOLATED == TRUE & !is.na(ZSCORE), .N]

  log_info("Flagged observations (for sensitivity analyses):")
  log_info("  GAMLSS_REFIT=TRUE: %d observations", n_refit)
  log_info("  AGE_EXTRAPOLATED=TRUE: %d observations", n_age_extrap)
  log_info("  ICC_EXTRAPOLATED=TRUE: %d observations", n_icc_extrap)
}

# =============================================================================
# NORMATIVE MODEL TRANSFER VALIDATION
# =============================================================================
# Validate that UK Biobank normative models transfer appropriately to ADNI
# by checking that A-beta negative cognitively unimpaired controls have
# HVR z-scores centered near zero (expected for healthy reference group).
# =============================================================================
log_section("Normative Model Transfer Validation")

# Load amyloid PET data (PHC_AMYLOID_STATUS: 0=negative, 1=positive)
amy_file <- get_data_path("raw", "adsp_phc_amyloid")
if (file.exists(amy_file)) {
  amy_cols.v <- c("PTID", "PHC_AMYLOID_STATUS", "PHC_Diagnosis")
  amy.dt <- fread(amy_file, select = amy_cols.v)

  # Load cognition data for diagnosis (PHC_Diagnosis: 1=CU, 2=MCI, 3=Dementia)
  cog_file <- get_data_path("raw", "adsp_phc_cognition")
  cog.dt <- fread(cog_file, select = c("PTID", "PHC_Diagnosis"))

  # Identify A-beta negative (status=0) cognitively unimpaired (DX=1)
  # If any positive scan, classify as positive (conservative)
  amy_baseline.dt <- amy.dt[
    order(PTID, -PHC_AMYLOID_STATUS)
  ][, .SD[1], by = PTID]
  cog_baseline.dt <- cog.dt[, .SD[1], by = PTID] # First visit diagnosis

  control_ids.v <- intersect(
    amy_baseline.dt[PHC_AMYLOID_STATUS == 0, unique(PTID)],
    cog_baseline.dt[PHC_Diagnosis == 1, unique(PTID)]
  )

  log_info(
    "A-beta negative CU controls identified: %d subjects",
    length(control_ids.v)
  )

  # UKB GAMLSS training sample age range: 46-82 years
  # Filter to this range for valid normative comparison
  UKB_AGE_MIN <- 46
  UKB_AGE_MAX <- 82

  # =========================================================================
  # COMPARABLE SUBSAMPLE: UKB-matched demographics
  # Criteria: A-beta negative, CU, age 46-82 (same as UKB training sample)
  # If normative models transfer well, z-scores should be ~N(0, 1)
  # =========================================================================
  log_info("COMPARABLE SUBSAMPLE ANALYSIS (UKB-matched demographics):")
  log_info("  Criteria: A-beta neg, CU, age %d-%d", UKB_AGE_MIN, UKB_AGE_MAX)

  # Get baseline z-scores for controls within UKB age range
  # Check both HVR and HC to assess normative model transfer
  control_hvr.dt <- adni_data.dt[
    PTID %in% control_ids.v &
      ROI == "HVR" &
      SIDE == "LR" &
      ADJ == "NON" &
      AGE >= UKB_AGE_MIN & AGE <= UKB_AGE_MAX &
      !is.na(ZSCORE)
  ][order(PTID, EXAMDATE)][, .SD[1], by = PTID]

  control_hc.dt <- adni_data.dt[
    PTID %in% control_ids.v &
      ROI == "HC" &
      SIDE == "LR" &
      ADJ == "NON" &
      AGE >= UKB_AGE_MIN & AGE <= UKB_AGE_MAX &
      !is.na(ZSCORE)
  ][order(PTID, EXAMDATE)][, .SD[1], by = PTID]

  n_comparable_hvr <- nrow(control_hvr.dt)
  n_comparable_hc <- nrow(control_hc.dt)

  # Comparable subsample statistics
  comparable_validation.lst <- list(
    ukb_age_range = c(min = UKB_AGE_MIN, max = UKB_AGE_MAX),
    criteria = "A-beta negative, CU, age 46-82"
  )

  if (n_comparable_hvr >= 10) {
    hvr_mean <- mean(control_hvr.dt$ZSCORE, na.rm = TRUE)
    hvr_sd <- sd(control_hvr.dt$ZSCORE, na.rm = TRUE)
    hvr_se <- hvr_sd / sqrt(n_comparable_hvr)
    hvr_t <- t.test(control_hvr.dt$ZSCORE, mu = 0)

    comparable_validation.lst$HVR <- list(
      n = n_comparable_hvr,
      mean = hvr_mean,
      sd = hvr_sd,
      se = hvr_se,
      ci_95 = c(lower = hvr_mean - 1.96 * hvr_se,
                upper = hvr_mean + 1.96 * hvr_se),
      t_statistic = hvr_t$statistic,
      t_p_value = hvr_t$p.value,
      deviation_from_expected = abs(hvr_mean) > 0.5 || abs(hvr_sd - 1) > 0.3,
      zscores = control_hvr.dt$ZSCORE
    )

    log_info("  HVR z-scores (N=%d): mean = %.3f, SD = %.3f",
             n_comparable_hvr, hvr_mean, hvr_sd)
    log_info("    Expected: mean ~ 0, SD ~ 1")
    if (comparable_validation.lst$HVR$deviation_from_expected) {
      log_warn("    DEVIATION: z-scores deviate from expected N(0,1)")
    } else {
      log_info("    OK: z-scores consistent with N(0,1)")
    }
  } else {
    log_warn("  HVR: Insufficient comparable subjects (N=%d)", n_comparable_hvr)
    comparable_validation.lst$HVR <- list(
      n = n_comparable_hvr,
      insufficient_data = TRUE
    )
  }

  if (n_comparable_hc >= 10) {
    hc_mean <- mean(control_hc.dt$ZSCORE, na.rm = TRUE)
    hc_sd <- sd(control_hc.dt$ZSCORE, na.rm = TRUE)
    hc_se <- hc_sd / sqrt(n_comparable_hc)
    hc_t <- t.test(control_hc.dt$ZSCORE, mu = 0)

    comparable_validation.lst$HC <- list(
      n = n_comparable_hc,
      mean = hc_mean,
      sd = hc_sd,
      se = hc_se,
      ci_95 = c(lower = hc_mean - 1.96 * hc_se,
                upper = hc_mean + 1.96 * hc_se),
      t_statistic = hc_t$statistic,
      t_p_value = hc_t$p.value,
      deviation_from_expected = abs(hc_mean) > 0.5 || abs(hc_sd - 1) > 0.3,
      zscores = control_hc.dt$ZSCORE
    )

    log_info("  HC z-scores (N=%d): mean = %.3f, SD = %.3f",
             n_comparable_hc, hc_mean, hc_sd)
    log_info("    Expected: mean ~ 0, SD ~ 1")
    if (comparable_validation.lst$HC$deviation_from_expected) {
      log_warn("    DEVIATION: z-scores deviate from expected N(0,1)")
    } else {
      log_info("    OK: z-scores consistent with N(0,1)")
    }
  } else {
    log_warn("  HC: Insufficient comparable subjects (N=%d)", n_comparable_hc)
    comparable_validation.lst$HC <- list(
      n = n_comparable_hc,
      insufficient_data = TRUE
    )
  }

  # Interpretation for manuscript
  interpret_transfer <- function(hvr.lst, hc.lst) {
    if (isTRUE(hvr.lst$insufficient_data) || isTRUE(hc.lst$insufficient_data)) {
      "Insufficient data for normative transfer validation"
    } else {
      hvr_ok <- !hvr.lst$deviation_from_expected
      hc_ok <- !hc.lst$deviation_from_expected

      if (hvr_ok && hc_ok) {
        sprintf(
          paste0(
            "UKB normative models transfer well to ADNI: ",
            "comparable subsample (N=%d) shows HC z-scores (M=%.2f, SD=%.2f) ",
            "and HVR z-scores (M=%.2f, SD=%.2f) consistent with expected ",
            "N(0,1) distribution."
          ),
          hc.lst$n, hc.lst$mean, hc.lst$sd, hvr.lst$mean, hvr.lst$sd
        )
      } else if (!hvr_ok && !hc_ok) {
        sprintf(
          paste0(
            "CAUTION: Both HC (M=%.2f, SD=%.2f) and HVR (M=%.2f, SD=%.2f) ",
            "z-scores deviate from expected N(0,1) in comparable subsample, ",
            "suggesting calibration differences between UKB and ADNI. ",
            "Raw HVR sensitivity analysis recommended."
          ),
          hc.lst$mean, hc.lst$sd, hvr.lst$mean, hvr.lst$sd
        )
      } else if (!hvr_ok) {
        sprintf(
          paste0(
            "HVR z-scores (M=%.2f, SD=%.2f) deviate from expected N(0,1), ",
            "while HC z-scores (M=%.2f, SD=%.2f) are calibrated. ",
            "Consider HC-based sensitivity analysis."
          ),
          hvr.lst$mean, hvr.lst$sd, hc.lst$mean, hc.lst$sd
        )
      } else {
        sprintf(
          paste0(
            "HC z-scores (M=%.2f, SD=%.2f) deviate from expected N(0,1), ",
            "while HVR z-scores (M=%.2f, SD=%.2f) are calibrated."
          ),
          hc.lst$mean, hc.lst$sd, hvr.lst$mean, hvr.lst$sd
        )
      }
    }
  }

  comparable_validation.lst$interpretation <- interpret_transfer(
    comparable_validation.lst$HVR, comparable_validation.lst$HC
  )
  log_info("")
  log_info("Interpretation: %s", comparable_validation.lst$interpretation)

  # =========================================================================
  # DEMOGRAPHICS TABLES: UKB Normative vs ADNI Comparable Subsample
  # =========================================================================
  log_info("")
  log_info("DEMOGRAPHICS TABLES:")


  # --- UKB Normative Sample Demographics (GAMLSS training) ---
  # Load from GAMLSS source data (the actual data used for model fitting)
  # This ensures we report demographics for the exact subjects used
  ukb_gamlss.path <- get_data_path(
    "external", "gamlss_models"
  )

  if (file.exists(ukb_gamlss.path)) {
    ukb_gamlss.lst <- readRDS(ukb_gamlss.path)
    ukb_crs.dt <- ukb_gamlss.lst$DATA$CRS

    # Get unique subjects (source data is in long format with multiple ROIs)
    ukb_subj.dt <- ukb_crs.dt[!duplicated(EID), .(EID, SEX, AGE, EDUC_num)]

    n_ukb <- nrow(ukb_subj.dt)
    n_ukb_female <- ukb_subj.dt[SEX == "Female", .N]
    n_ukb_male <- ukb_subj.dt[SEX == "Male", .N]

    ukb_demog.lst <- list(
      sample = "UKB Normative",
      description = "UK Biobank GAMLSS training sample",
      n_total = n_ukb,
      n_female = n_ukb_female,
      n_male = n_ukb_male,
      pct_female = 100 * n_ukb_female / n_ukb,
      age_mean = ukb_subj.dt[, mean(AGE, na.rm = TRUE)],
      age_sd = ukb_subj.dt[, sd(AGE, na.rm = TRUE)],
      age_min = ukb_subj.dt[, min(AGE, na.rm = TRUE)],
      age_max = ukb_subj.dt[, max(AGE, na.rm = TRUE)],
      educ_mean = ukb_subj.dt[, mean(EDUC_num, na.rm = TRUE)],
      educ_sd = ukb_subj.dt[, sd(EDUC_num, na.rm = TRUE)]
    )
    rm(ukb_gamlss.lst, ukb_crs.dt, ukb_subj.dt)
  } else {
    log_warn("UKB GAMLSS file not found: %s", ukb_gamlss.path)
    ukb_demog.lst <- list(
      sample = "UKB Normative",
      description = "UK Biobank (GAMLSS file not found)",
      n_total = NA, n_female = NA, n_male = NA, pct_female = NA,
      age_mean = NA, age_sd = NA, age_min = NA, age_max = NA,
      educ_mean = NA, educ_sd = NA
    )
  }

  log_info("  UKB Normative Sample (N=%d):", ukb_demog.lst$n_total)
  log_info("    Sex: %d F (%.1f%%), %d M",
           ukb_demog.lst$n_female, ukb_demog.lst$pct_female,
           ukb_demog.lst$n_male)
  log_info("    Age: %.1f (SD=%.1f), range %.0f-%.0f",
           ukb_demog.lst$age_mean, ukb_demog.lst$age_sd,
           ukb_demog.lst$age_min, ukb_demog.lst$age_max)
  log_info("    Education: %.1f years (SD=%.1f)",
           ukb_demog.lst$educ_mean, ukb_demog.lst$educ_sd)

  # --- ADNI Comparable Subsample Demographics ---
  # Get full demographics for subjects in comparable subsample
  comparable_ptids.v <- control_hvr.dt$PTID

  # Merge with full demographics
  comparable_demog.dt <- adni_merged.dt[
    PTID %in% comparable_ptids.v
  ][order(PTID, EXAMDATE)][, .SD[1], by = PTID]

  n_comp <- nrow(comparable_demog.dt)
  n_comp_female <- comparable_demog.dt[SEX == "Female", .N]
  n_comp_male <- comparable_demog.dt[SEX == "Male", .N]

  adni_comp_demog.lst <- list(
    sample = "ADNI Comparable",
    description = "ADNI A-beta neg, CU, age 46-82 (matches UKB criteria)",
    n_total = n_comp,
    n_female = n_comp_female,
    n_male = n_comp_male,
    pct_female = 100 * n_comp_female / n_comp,
    age_mean = comparable_demog.dt[, mean(AGE, na.rm = TRUE)],
    age_sd = comparable_demog.dt[, sd(AGE, na.rm = TRUE)],
    age_min = comparable_demog.dt[, min(AGE, na.rm = TRUE)],
    age_max = comparable_demog.dt[, max(AGE, na.rm = TRUE)],
    educ_mean = comparable_demog.dt[, mean(EDUC_num, na.rm = TRUE)],
    educ_sd = comparable_demog.dt[, sd(EDUC_num, na.rm = TRUE)]
  )

  log_info("  ADNI Comparable Subsample (N=%d):", adni_comp_demog.lst$n_total)
  log_info("    Sex: %d F (%.1f%%), %d M",
           adni_comp_demog.lst$n_female, adni_comp_demog.lst$pct_female,
           adni_comp_demog.lst$n_male)
  log_info("    Age: %.1f (SD=%.1f), range %.0f-%.0f",
           adni_comp_demog.lst$age_mean, adni_comp_demog.lst$age_sd,
           adni_comp_demog.lst$age_min, adni_comp_demog.lst$age_max)
  log_info("    Education: %.1f years (SD=%.1f)",
           adni_comp_demog.lst$educ_mean, adni_comp_demog.lst$educ_sd)

  # Combine demographics for output
  comparable_validation.lst$demographics <- list(
    ukb_normative = ukb_demog.lst,
    adni_comparable = adni_comp_demog.lst
  )

  # Create demographics comparison table (data.table format)
  demog_comparison.dt <- data.table(
    Sample = c("UKB Normative", "ADNI Comparable"),
    N = c(ukb_demog.lst$n_total, adni_comp_demog.lst$n_total),
    Female_N = c(ukb_demog.lst$n_female, adni_comp_demog.lst$n_female),
    Female_Pct = c(ukb_demog.lst$pct_female, adni_comp_demog.lst$pct_female),
    Age_Mean = c(ukb_demog.lst$age_mean, adni_comp_demog.lst$age_mean),
    Age_SD = c(ukb_demog.lst$age_sd, adni_comp_demog.lst$age_sd),
    Age_Min = c(ukb_demog.lst$age_min, adni_comp_demog.lst$age_min),
    Age_Max = c(ukb_demog.lst$age_max, adni_comp_demog.lst$age_max),
    Educ_Mean = c(ukb_demog.lst$educ_mean, adni_comp_demog.lst$educ_mean),
    Educ_SD = c(ukb_demog.lst$educ_sd, adni_comp_demog.lst$educ_sd)
  )

  comparable_validation.lst$demographics_table <- demog_comparison.dt

  # --- Generate gt demographic table ---
  if (requireNamespace("gt", quietly = TRUE)) {
    library(gt)

    # Create formatted table for publication
    demog_gt.tbl <- demog_comparison.dt[, .(
      Sample,
      N = format(N, big.mark = ","),
      `Female, N (%)` = sprintf("%s (%.1f%%)",
                                format(Female_N, big.mark = ","), Female_Pct),
      `Age, M (SD)` = sprintf("%.1f (%.1f)", Age_Mean, Age_SD),
      `Age Range` = sprintf("%.0f-%.0f", Age_Min, Age_Max),
      `Education, M (SD)` = sprintf("%.1f (%.1f)", Educ_Mean, Educ_SD)
    )] |>
      gt() |>
      tab_header(
        title = "Normative Transfer Validation: Sample Demographics",
        subtitle = "UKB GAMLSS training sample vs ADNI comparable subsample"
      ) |>
      tab_source_note(
        "ADNI Comparable: A-beta negative, cognitively unimpaired, age 46-82"
      ) |>
      cols_align(align = "right", columns = -Sample) |>
      cols_align(align = "left", columns = Sample)

    # Save table
    demog_table.path <- here(
      "outputs/tables/normative_transfer_demographics.html"
    )
    ensure_directory(dirname(demog_table.path))
    gtsave(demog_gt.tbl, demog_table.path)
    log_info("Saved: %s", demog_table.path)

    # Also save RDS for programmatic access
    demog_rds.path <- here(
      "outputs/tables/normative_transfer_demographics.rds"
    )
    write_rds_safe(demog_comparison.dt, demog_rds.path,
                   description = "Normative transfer demographics")
  } else {
    log_warn("gt package not available, skipping demographics table export")
  }

  # =========================================================================
  # IPW SENSITIVITY: Age-Weighted Validation
  # =========================================================================
  # The unweighted validation shows HVR z-scores deviate from zero.

  # However, ADNI is ~8 years older on average (71.3 vs 63.4) even within
  # the same age range. IPW can test whether this age mismatch explains
  # the offset.
  #
  # Method: Weight ADNI subjects to match UKB age distribution using
  # inverse probability weights derived from age bin proportions.
  # =========================================================================
  log_info("")
  log_info("IPW SENSITIVITY: Age-Weighted Validation")

  # Create age bins (5-year intervals matching typical epidemiological bins)
  age_breaks.v <- seq(45, 85, by = 5)
  control_hvr.dt[, age_bin := cut(AGE, breaks = age_breaks.v,
                                   include.lowest = TRUE, right = FALSE)]

  # UKB age distribution (approximate from normal based on demographics)
  ukb_age_probs.v <- diff(pnorm(age_breaks.v,
                                 mean = ukb_demog.lst$age_mean,
                                 sd = ukb_demog.lst$age_sd))
  ukb_age_probs.v <- ukb_age_probs.v / sum(ukb_age_probs.v)
  names(ukb_age_probs.v) <- levels(control_hvr.dt$age_bin)

  # ADNI age distribution (observed)
  adni_age_counts.v <- table(control_hvr.dt$age_bin)
  adni_age_probs.v <- adni_age_counts.v / sum(adni_age_counts.v)

  # Calculate IPW: w_i = P_UKB(age_bin) / P_ADNI(age_bin)
  ipw_table.dt <- data.table(
    age_bin = names(ukb_age_probs.v),
    ukb_prob = ukb_age_probs.v,
    adni_prob = as.numeric(adni_age_probs.v[names(ukb_age_probs.v)]),
    adni_n = as.numeric(adni_age_counts.v[names(ukb_age_probs.v)])
  )
  ipw_table.dt[is.na(adni_prob) | adni_prob == 0, adni_prob := 0.001]
  ipw_table.dt[, weight := ukb_prob / adni_prob]
  ipw_table.dt[is.infinite(weight), weight := 0]

  log_info("  Age bin weights for IPW:")
  for (i in seq_len(nrow(ipw_table.dt))) {
    row <- ipw_table.dt[i]
    log_info("    %s: UKB=%.3f, ADNI=%.3f (n=%d), weight=%.2f",
             row$age_bin, row$ukb_prob, row$adni_prob, row$adni_n, row$weight)
  }

  # Assign weights to each observation
  control_hvr.dt <- merge(control_hvr.dt,
                           ipw_table.dt[, .(age_bin, ipw = weight)],
                           by = "age_bin", all.x = TRUE)
  control_hvr.dt[is.na(ipw), ipw := 0]

  # Compute IPW-weighted statistics
  w <- control_hvr.dt$ipw
  w_norm <- w / sum(w, na.rm = TRUE)
  ipw_mean <- weighted.mean(control_hvr.dt$ZSCORE, w, na.rm = TRUE)

  # Kish effective sample size
  n_eff <- 1 / sum(w_norm^2, na.rm = TRUE)

  # Weighted SD and SE
  ipw_var <- sum(w_norm * (control_hvr.dt$ZSCORE - ipw_mean)^2, na.rm = TRUE)
  ipw_sd <- sqrt(ipw_var * n_comparable_hvr / (n_comparable_hvr - 1))
  ipw_se <- ipw_sd / sqrt(n_eff)

  # Weighted t-test (approximate)
  ipw_t <- ipw_mean / ipw_se
  ipw_p <- 2 * pt(-abs(ipw_t), df = n_eff - 1)
  ipw_ci <- c(lower = ipw_mean - 1.96 * ipw_se,
              upper = ipw_mean + 1.96 * ipw_se)

  # Does IPW resolve the offset?
  ipw_centered <- abs(ipw_mean) < 0.5 && ipw_ci["lower"] < 0.5

  log_info("  IPW-weighted results:")
  log_info("    N effective: %.1f (of %d)", n_eff, n_comparable_hvr)
  log_info("    Mean HVR z-score: %.3f (SE = %.3f)", ipw_mean, ipw_se)
  log_info("    SD: %.3f", ipw_sd)
  log_info("    95%% CI: [%.3f, %.3f]", ipw_ci["lower"], ipw_ci["upper"])
  log_info("    t = %.2f, p = %.4f", ipw_t, ipw_p)

  if (ipw_centered) {
    log_info("  RESULT: IPW RESOLVES the offset (mean within [-0.5, 0.5])")
    ipw_interpretation <- paste0(
      "IPW correction resolves the validation offset. After weighting ADNI ",
      "to match UKB age distribution, mean HVR z-score = ",
      sprintf("%.2f (95%% CI: %.2f, %.2f)", ipw_mean, ipw_ci["lower"],
              ipw_ci["upper"]),
      ", which is acceptable. The original offset was due to age distribution ",
      "mismatch (ADNI ~8 years older), not normative model miscalibration."
    )
  } else {
    log_info("  RESULT: IPW does NOT fully resolve the offset")
    ipw_interpretation <- paste0(
      "IPW correction does not fully resolve the validation offset. After ",
      "weighting, mean HVR z-score = ",
      sprintf("%.2f (95%% CI: %.2f, %.2f)", ipw_mean, ipw_ci["lower"],
              ipw_ci["upper"]),
      ". The persistent offset may indicate: (1) genuine cohort differences ",
      "in HVR beyond age, (2) scanner/protocol differences, or (3) selection ",
      "effects (ADNI recruits memory-concerned individuals)."
    )
  }

  log_info("  Interpretation: %s", ipw_interpretation)

  # Add IPW results to comparable_validation.lst
  comparable_validation.lst$ipw <- list(
    description = "IPW by age (5-year bins) to match UKB distribution",
    n = n_comparable_hvr,
    n_effective = n_eff,
    mean = ipw_mean,
    sd = ipw_sd,
    se = ipw_se,
    ci_95 = ipw_ci,
    t_statistic = ipw_t,
    p_value = ipw_p,
    centered_near_zero = ipw_centered,
    weight_table = ipw_table.dt,
    interpretation = ipw_interpretation
  )

  # Stratified validation by age bin (additional sensitivity)
  stratified_results.dt <- control_hvr.dt[, .(
    n = .N,
    mean_z = mean(ZSCORE, na.rm = TRUE),
    sd_z = sd(ZSCORE, na.rm = TRUE),
    se_z = sd(ZSCORE, na.rm = TRUE) / sqrt(.N)
  ), by = age_bin][order(age_bin)]
  stratified_results.dt[, ci_lower := mean_z - 1.96 * se_z]
  stratified_results.dt[, ci_upper := mean_z + 1.96 * se_z]
  stratified_results.dt[, includes_zero := ci_lower <= 0 & ci_upper >= 0]

  log_info("  Stratified by age bin:")
  for (i in seq_len(nrow(stratified_results.dt))) {
    row <- stratified_results.dt[i]
    zero_str <- ifelse(row$includes_zero, "includes 0", "EXCLUDES 0")
    log_info("    %s (n=%d): mean=%.2f, 95%% CI [%.2f, %.2f] - %s",
             row$age_bin, row$n, row$mean_z, row$ci_lower, row$ci_upper,
             zero_str)
  }

  comparable_validation.lst$stratified_hvr <- stratified_results.dt

  # Clean up temporary column
  control_hvr.dt[, c("age_bin", "ipw") := NULL]

  # ---------------------------------------------------------
  # HC stratified validation by age bin (parallel to HVR)
  # ---------------------------------------------------------
  control_hc.dt[, age_bin := cut(
    AGE, breaks = age_breaks.v,
    include.lowest = TRUE, right = FALSE
  )]

  stratified_hc.dt <- control_hc.dt[, .(
    n = .N,
    mean_z = mean(ZSCORE, na.rm = TRUE),
    sd_z = sd(ZSCORE, na.rm = TRUE),
    se_z = sd(ZSCORE, na.rm = TRUE) / sqrt(.N)
  ), by = age_bin][order(age_bin)]
  stratified_hc.dt[, ci_lower := mean_z - 1.96 * se_z]
  stratified_hc.dt[, ci_upper := mean_z + 1.96 * se_z]
  stratified_hc.dt[,
    includes_zero := ci_lower <= 0 & ci_upper >= 0
  ]

  log_info("  HC Stratified by age bin:")
  for (i in seq_len(nrow(stratified_hc.dt))) {
    row <- stratified_hc.dt[i]
    zero_str <- ifelse(
      row$includes_zero, "includes 0", "EXCLUDES 0"
    )
    log_info(
      "    %s (n=%d): mean=%.2f, 95%% CI [%.2f, %.2f] - %s",
      row$age_bin, row$n, row$mean_z,
      row$ci_lower, row$ci_upper, zero_str
    )
  }

  comparable_validation.lst$stratified_hc <- stratified_hc.dt

  # Clean up temporary column
  control_hc.dt[, age_bin := NULL]

  # =========================================================================
  # FULL SAMPLE VALIDATION (for backwards compatibility)
  # =========================================================================
  log_info("")
  log_info("FULL SAMPLE VALIDATION (all A-beta neg CU, no age filter):")

  # For backwards compatibility, also keep full sample without age filter
  control_zscores.dt <- adni_data.dt[
    PTID %in% control_ids.v &
      ROI == "HVR" &
      SIDE == "LR" &
      ADJ == "NON" &
      !is.na(ZSCORE)
  ][order(PTID, EXAMDATE)][, .SD[1], by = PTID]

  n_controls_with_zscore <- nrow(control_zscores.dt)

  if (n_controls_with_zscore >= 10) {
    # Compute validation statistics
    mean_z <- mean(control_zscores.dt$ZSCORE, na.rm = TRUE)
    sd_z <- sd(control_zscores.dt$ZSCORE, na.rm = TRUE)
    se_z <- sd_z / sqrt(n_controls_with_zscore)
    ci_lower <- mean_z - 1.96 * se_z
    ci_upper <- mean_z + 1.96 * se_z

    # One-sample t-test against mu=0
    t_result <- t.test(control_zscores.dt$ZSCORE, mu = 0)

    # Determine if centered near zero (mean within 0.5 SD)
    centered_near_zero <- abs(mean_z) < 0.5 &&
      ci_lower < 0.25 && ci_upper > -0.25

    # Build interpretation string
    interp_ok <- paste0(
      "Mean HVR z-score for ADNI A-beta negative CU controls is ",
      "%.3f (95%% CI: %.3f, %.3f), confirming UKB norms transfer."
    )
    interp_warn <- paste0(
      "WARNING: Mean HVR z-score for ADNI CU controls is ",
      "%.3f (95%% CI: %.3f, %.3f), suggesting calibration offset."
    )
    interp_str <- if (centered_near_zero) {
      sprintf(interp_ok, mean_z, ci_lower, ci_upper)
    } else {
      sprintf(interp_warn, mean_z, ci_lower, ci_upper)
    }

    transfer_validation.lst <- list(
      n_controls = length(control_ids.v),
      n_with_zscore = n_controls_with_zscore,
      mean_z = mean_z,
      sd_z = sd_z,
      se_z = se_z,
      ci_95 = c(lower = ci_lower, upper = ci_upper),
      t_statistic = t_result$statistic,
      t_p_value = t_result$p.value,
      centered_near_zero = centered_near_zero,
      interpretation = interp_str
    )

    log_info("NORMATIVE TRANSFER VALIDATION RESULTS:")
    log_info("  A-beta negative CU controls: N = %d", n_controls_with_zscore)
    log_info("  Mean HVR z-score: %.3f (SD = %.3f)", mean_z, sd_z)
    log_info("  95%% CI: [%.3f, %.3f]", ci_lower, ci_upper)
    log_info(
      "  t-test vs mu=0: t = %.3f, p = %.4f",
      t_result$statistic, t_result$p.value
    )
    centered_str <- ifelse(centered_near_zero, "YES", "NO - INVESTIGATE")
    log_info("  Centered near zero: %s", centered_str)

    # Save validation results (combine full sample + comparable subsample)
    combined_validation.lst <- list(
      full_sample = transfer_validation.lst,
      comparable_subsample = comparable_validation.lst
    )

    validation_path <- get_data_path(
      "derivatives", "normative_transfer_validation"
    )
    ensure_directory(dirname(validation_path))
    write_rds_safe(
      combined_validation.lst, validation_path,
      description = "Normative model transfer validation"
    )
  } else {
    log_warn(
      "Insufficient CU controls for validation (N=%d, need >=10)",
      n_controls_with_zscore
    )
    interp_insufficient <- "Insufficient CU controls for robust validation"
    transfer_validation.lst <- list(
      n_controls = length(control_ids.v),
      n_with_zscore = n_controls_with_zscore,
      insufficient_data = TRUE,
      interpretation = interp_insufficient
    )

    # Save validation results (combine full sample + comparable subsample)
    combined_validation.lst <- list(
      full_sample = transfer_validation.lst,
      comparable_subsample = comparable_validation.lst
    )

    validation_path <- get_data_path(
      "derivatives", "normative_transfer_validation"
    )
    ensure_directory(dirname(validation_path))
    write_rds_safe(
      combined_validation.lst, validation_path,
      description = "Normative transfer validation (insufficient)"
    )
  }
} else {
  log_warn("Amyloid PET file not found: %s", amy_file)
  log_warn("Skipping normative transfer validation")
}

# --- Prepare Output ---
log_section("Preparing Output")

# Rename ZSCORE to Z for consistency with reference
output.dt <- copy(adni_data.dt)
setnames(output.dt, "ZSCORE", "Z")

# Add metadata
attr(output.dt, "creation_date") <- Sys.time()
attr(output.dt, "n_subjects") <- uniqueN(output.dt$PTID)
attr(output.dt, "n_observations") <- nrow(output.dt)

log_info(
  "Output dataset: %d observations from %d subjects",
  nrow(output.dt), uniqueN(output.dt$PTID)
)

# --- Save Output ---
log_section("Saving Results")

ensure_directory(dirname(output.path))
write_rds_safe(
  output.dt, output.path,
  description = "ADNI Z-scores from UKB normative models"
)

# Save diagnostics
diag_path <- get_data_path("models", "diagnostics", "zscore")
if (!is.null(diag_path)) {
  ensure_directory(dirname(diag_path))
  write_rds_safe(
    zscore_diag.dt, diag_path,
    description = "Z-score calculation diagnostics"
  )
}

log_script_end("01_compute_zscores.R", success = TRUE)
