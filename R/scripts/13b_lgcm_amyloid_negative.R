#!/usr/bin/env Rscript

# =============================================================================
# 13b_lgcm_amyloid_negative.R
# =============================================================================
# LGCM mediation in amyloid-negative, cognitively unimpaired (CU).
#
# Follows up on LME results (10b) showing CVR_mimic moderates
# the HVR-cognition relationship in A-/CU.  Tests whether
# CVR_mimic also mediates through hippocampal decline (a-path).
#
# SPECIFICATION:
#   - Linear LGCM (no quadratic — smaller sample)
#   - YRS_from_bl as time metric (EDT not meaningful in A-/CU)
#   - Pooled sample, SEX as covariate (no sex stratification)
#   - SE²-constrained cognitive residuals (same as A+ models)
#   - No bootstrap (point estimates + Wald SEs only)
#
# INPUTS:
#   - models/results/lme/lme_amyloid_negative.rds (cohort)
#   - data/derivatives/lme/cvr_mimic_scores.rds (MIMIC indicators)
#
# OUTPUTS:
#   - models/results/lgcm/lgcm_amyloid_negative.rds
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(OpenMx)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

log_script_start("13b_lgcm_amyloid_negative.R")
config <- load_config()
validate_config(config)
set_seed()

mxOption(NULL, "Number of Threads",
         parallel::detectCores() - 1L)

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "lgcm_amy_neg",
  default = FALSE
)
output.path <- file.path(
  get_data_path("models", "lgcm_results_dir"),
  "lgcm_amyloid_negative.rds"
)
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists; skipping.")
  log_script_end(
    "13b_lgcm_amyloid_negative.R", success = TRUE
  )
  quit(status = 0)
}

COG_DOMAINS <- c("MEM", "LAN", "EXF")

# ============================================================
# 1. Load and reshape data
# ============================================================
log_section("Part 1: Data Preparation")

lme_res.lst <- read_rds_safe(
  file.path(
    get_data_path("models", "lme_results_dir"),
    "lme_amyloid_negative.rds"
  ),
  "A-/CU LME results"
)
long.dt <- copy(lme_res.lst$cohort)
log_info(
  "Loaded A-/CU cohort: %d subjects, %d obs",
  long.dt[, uniqueN(PTID)], nrow(long.dt)
)

# --- Merge MIMIC indicators for MIMIC-embedded models ---
mimic.dt <- read_rds_safe(
  get_data_path("models", "cvr_mimic_scores"),
  "CVR MIMIC scores"
)
MIMIC_IND <- c(
  "SBP_z", "HTN", "GLUCOSE_z",
  "CHOL_z", "CREAT_z"
)

# Map RID -> PTID
mimic.dt[, RID_char := sprintf("%04d", RID)]
ptid_map.dt <- unique(long.dt[, .(PTID)])
ptid_map.dt[, RID_char := sub(
  "^\\d+_S_", "", PTID
)]
mimic_merge.dt <- merge(
  ptid_map.dt,
  mimic.dt[, c("RID_char", MIMIC_IND),
           with = FALSE],
  by = "RID_char", all.x = TRUE
)
mimic_merge.dt[, RID_char := NULL]
# One row per subject — take first
mimic_merge.dt <- mimic_merge.dt[
  , .SD[1], by = PTID
]
for (col in MIMIC_IND) {
  if (col %in% names(long.dt))
    long.dt[, (col) := NULL]
}
long.dt <- merge(
  long.dt, mimic_merge.dt,
  by = "PTID", all.x = TRUE
)
log_info(
  "MIMIC indicators merged: %d/%d with SBP_z",
  sum(!is.na(long.dt$SBP_z)), nrow(long.dt)
)

# --- Require >=2 distinct MRI scans ---
# Rolling join can match one MRI to multiple
# cognitive visits; these subjects have no
# longitudinal brain data for slope estimation
n_mri.dt <- long.dt[
  , .(n_mri = uniqueN(MRI_date)), by = PTID
]
mri_keep.v <- n_mri.dt[n_mri >= 2, PTID]
n_dropped <- long.dt[, uniqueN(PTID)] -
  length(mri_keep.v)
long.dt <- long.dt[PTID %in% mri_keep.v]
log_info(
  "After >=2 distinct MRI scans: %d subjects "
  , long.dt[, uniqueN(PTID)]
)
if (n_dropped > 0) {
  log_info(
    "  Dropped %d with single MRI scan",
    n_dropped
  )
}

# --- Rename cognitive columns (if not already) ---
rename_map <- c(
  MEM = "PHC_MEM", LAN = "PHC_LAN",
  EXF = "PHC_EXF",
  MEM_SE = "PHC_MEM_SE", LAN_SE = "PHC_LAN_SE",
  EXF_SE = "PHC_EXF_SE"
)
for (old_nm in names(rename_map)) {
  new_nm <- rename_map[[old_nm]]
  if (old_nm %in% names(long.dt) &&
      !new_nm %in% names(long.dt)) {
    setnames(long.dt, old_nm, new_nm)
  }
}

# --- Timepoint label ---
setorder(long.dt, PTID, VISIT)
long.dt[, TIMEPOINT := paste0("T", VISIT - 1)]

# Determine max timepoints per subject
max_tp <- long.dt[, max(VISIT), by = PTID]
n_tp_use <- min(
  max(max_tp$V1),
  4L  # Cap at T0-T3 (smaller sample than A+)
)
VALID_TP <- paste0("T", 0:(n_tp_use - 1))
long.dt <- long.dt[TIMEPOINT %in% VALID_TP]
log_info("Using timepoints: %s",
         paste(VALID_TP, collapse = ", "))

# --- YRS columns for definition variables ---
YRS_COLS <- paste0("YRS_", VALID_TP)
HVR_COLS <- paste0("HVR_Z_", VALID_TP)

# --- Pivot to wide ---
# Time-varying: HVR_z, cognitive, SE, YRS
tv_cols.v <- c(
  "HVR_z", "PHC_MEM", "PHC_LAN", "PHC_EXF",
  "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE",
  "YRS_from_bl"
)
id_cols.v <- c(
  "PTID", "SEX", "Age_bl", "EDUC", "APOE4",
  "FRS", "FRS_raw", "CVR_mimic",
  MIMIC_IND
)

# Keep only needed columns
keep.v <- c(id_cols.v, tv_cols.v, "TIMEPOINT")
keep.v <- intersect(keep.v, names(long.dt))
long_sub.dt <- long.dt[, ..keep.v]

# Baseline row for time-invariant
baseline.dt <- long_sub.dt[
  TIMEPOINT == "T0",
  ..id_cols.v
]
baseline.dt <- baseline.dt[
  , .SD[1], by = PTID
]

# Pivot time-varying
wide.dt <- dcast(
  long_sub.dt,
  PTID ~ TIMEPOINT,
  value.var = tv_cols.v
)
wide.dt <- merge(
  baseline.dt, wide.dt, by = "PTID"
)
log_info(
  "Wide format: %d subjects, %d columns",
  nrow(wide.dt), ncol(wide.dt)
)

# --- Require >=2 HVR + >=2 cognitive timepoints ---
hvr_tp.v <- grep(
  "^HVR_z_T", names(wide.dt), value = TRUE
)
mem_tp.v <- grep(
  "^PHC_MEM_T", names(wide.dt), value = TRUE
)
wide.dt[, n_hvr := rowSums(
  !is.na(.SD)
), .SDcols = hvr_tp.v]
wide.dt[, n_cog := rowSums(
  !is.na(.SD)
), .SDcols = mem_tp.v]

# Count distinct HVR values (not just non-NA)
# Rolling join can repeat one MRI across visits
wide.dt[, n_hvr_distinct := apply(
  .SD, 1, function(x) {
    length(unique(x[!is.na(x)]))
  }
), .SDcols = hvr_tp.v]

n_pre <- nrow(wide.dt)
wide.dt <- wide.dt[
  n_hvr_distinct >= 2 & n_cog >= 2
]
log_info(
  "After >=2 distinct MRI + >=2 cog: %d/%d",
  nrow(wide.dt), n_pre
)
wide.dt[
  , c("n_hvr", "n_cog", "n_hvr_distinct") := NULL
]

# --- Rename HVR_z columns ---
for (tp in VALID_TP) {
  old_nm <- paste0("HVR_z_", tp)
  new_nm <- paste0("HVR_Z_", tp)
  if (old_nm %in% names(wide.dt)) {
    setnames(wide.dt, old_nm, new_nm)
  }
}
# Rename YRS columns
for (tp in VALID_TP) {
  old_nm <- paste0("YRS_from_bl_", tp)
  new_nm <- paste0("YRS_", tp)
  if (old_nm %in% names(wide.dt)) {
    setnames(wide.dt, old_nm, new_nm)
  }
}

# Fill NA definition variables with 0
# (corresponding manifest vars are also NA so
# these rows contribute no likelihood; OpenMx
# requires non-NA definition variables)
for (tp in VALID_TP) {
  yrs_col <- paste0("YRS_", tp)
  if (yrs_col %in% names(wide.dt)) {
    wide.dt[is.na(get(yrs_col)),
            (yrs_col) := 0]
  }
}

# --- Center predictors ---
age_mean <- mean(wide.dt$Age_bl, na.rm = TRUE)
educ_mean <- mean(wide.dt$EDUC, na.rm = TRUE)
wide.dt[, AGE_c := Age_bl - age_mean]
wide.dt[, EDUC_c := EDUC - educ_mean]
wide.dt[, SEX_num := as.integer(SEX == "Male")]

# FRS: center on sample mean
frs_mean <- mean(wide.dt$FRS, na.rm = TRUE)
wide.dt[, FRS_c := FRS - frs_mean]

# CVR_mimic: center on sample mean
cvr_mean <- mean(
  wide.dt$CVR_mimic, na.rm = TRUE
)
wide.dt[, CVR_c := CVR_mimic - cvr_mean]

# APOE4 numeric
wide.dt[, APOE4_num := as.numeric(APOE4)]

HAS_APOE <- sum(
  !is.na(wide.dt$APOE4_num)
) > 0

log_info("AGE_c mean: %.2f", age_mean)
log_info("EDUC_c mean: %.2f", educ_mean)
log_info("FRS_c mean: %.4f", frs_mean)
log_info("CVR_c mean: %.4f", cvr_mean)
log_info("N final: %d", nrow(wide.dt))
log_info("APOE4 available: %s", HAS_APOE)

# --- Compute SE² for cognitive residual fixing ---
se_cols.lst <- list(
  MEM = paste0("PHC_MEM_SE_", VALID_TP),
  LAN = paste0("PHC_LAN_SE_", VALID_TP),
  EXF = paste0("PHC_EXF_SE_", VALID_TP)
)

# ============================================================
# 2. Model builder (linear LGCM mediation)
# ============================================================
log_section("Part 2: Model Specification")

build_linear_mediation.fn <- function(
    cog_vars, hvr_vars, predictor,
    name = "LinMed", data = NULL) {

  n_tp <- length(cog_vars)
  latents <- c(
    "hvr_i", "hvr_s", "cog_i", "cog_s"
  )

  # HVR trajectory (linear, YRS definition vars)
  hvr_int <- mxPath(
    from = "hvr_i", to = hvr_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  hvr_slope <- mxPath(
    from = "hvr_s", to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", YRS_COLS[1:n_tp])
  )

  # Cognitive trajectory
  cog_int <- mxPath(
    from = "cog_i", to = cog_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  cog_slope <- mxPath(
    from = "cog_s", to = cog_vars,
    free = FALSE,
    labels = paste0("data.", YRS_COLS[1:n_tp])
  )

  # Growth factor means and (co)variances
  # Starting values based on observed data:
  # HVR_z var ~1.5, cog var ~0.3
  lat_means <- mxPath(
    from = "one", to = latents,
    free = TRUE,
    values = c(-0.5, -0.02, 0, -0.01),
    labels = paste0("mean_", latents)
  )
  lat_vars <- mxPath(
    from = latents, arrows = 2,
    free = TRUE,
    values = c(0.8, 0.001, 0.2, 0.001),
    lbound = c(0, 0, 0, 0),
    labels = paste0("var_", latents)
  )
  hvr_cov <- mxPath(
    from = "hvr_i", to = "hvr_s",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_hvr_is"
  )
  cog_cov <- mxPath(
    from = "cog_i", to = "cog_s",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_cog_is"
  )
  cross_int <- mxPath(
    from = "hvr_i", to = "cog_i",
    arrows = 2, free = TRUE, values = 0,
    labels = "cov_hvr_i_cog_i"
  )

  # HVR residuals: free
  hvr_resid <- mxPath(
    from = hvr_vars, arrows = 2,
    free = TRUE, values = rep(0.3, n_tp),
    lbound = rep(1e-6, n_tp)
  )

  # Cognitive residuals: free estimation
  # (SE² constraint used in A+ models requires
  # larger time span for stable decomposition;
  # with shorter YRS_from_bl the constraint
  # causes non-PD implied covariance)
  cog_resid <- mxPath(
    from = cog_vars, arrows = 2,
    free = TRUE, values = rep(0.1, n_tp),
    lbound = rep(1e-6, n_tp)
  )

  # Manifest means fixed to 0
  hvr_mm <- mxPath(
    from = "one", to = hvr_vars,
    free = FALSE, values = rep(0, n_tp)
  )
  cog_mm <- mxPath(
    from = "one", to = cog_vars,
    free = FALSE, values = rep(0, n_tp)
  )

  # Predictors
  preds <- c(
    predictor, "AGE_c", "EDUC_c", "SEX_num"
  )
  if (HAS_APOE) preds <- c(preds, "APOE4_num")

  pred_var <- mxPath(
    from = preds, arrows = 2,
    free = FALSE, values = 1
  )
  pred_mn <- mxPath(
    from = "one", to = preds,
    free = FALSE, values = 0
  )

  # Mediation paths
  a_path <- mxPath(
    from = predictor, to = "hvr_s",
    free = TRUE, values = 0, labels = "a"
  )
  b_path <- mxPath(
    from = "hvr_s", to = "cog_s",
    free = TRUE, values = 0.5, labels = "b"
  )
  cprime_path <- mxPath(
    from = predictor, to = "cog_s",
    free = TRUE, values = 0,
    labels = "cprime"
  )
  pred_int <- mxPath(
    from = predictor,
    to = c("hvr_i", "cog_i"),
    free = TRUE, values = 0,
    labels = c("pred_hvr_i", "pred_cog_i")
  )

  # Covariate paths to all latent factors
  cov_labels.fn <- function(prefix) {
    paste0(prefix, "_", latents)
  }

  age_paths <- mxPath(
    from = "AGE_c", to = latents,
    free = TRUE, values = 0,
    labels = cov_labels.fn("age")
  )
  educ_paths <- mxPath(
    from = "EDUC_c", to = latents,
    free = TRUE, values = 0,
    labels = cov_labels.fn("educ")
  )
  sex_paths <- mxPath(
    from = "SEX_num", to = latents,
    free = TRUE, values = 0,
    labels = cov_labels.fn("sex")
  )

  # Build model
  model <- mxModel(
    name, type = "RAM",
    manifestVars = c(hvr_vars, cog_vars, preds),
    latentVars = latents
  )
  model <- mxModel(
    model,
    hvr_int, hvr_slope, cog_int, cog_slope,
    lat_means, lat_vars,
    hvr_cov, cog_cov, cross_int,
    hvr_resid, cog_resid,
    hvr_mm, cog_mm,
    pred_var, pred_mn,
    a_path, b_path, cprime_path, pred_int,
    age_paths, educ_paths, sex_paths
  )

  if (HAS_APOE) {
    apoe_paths <- mxPath(
      from = "APOE4_num", to = latents,
      free = TRUE, values = 0,
      labels = cov_labels.fn("apoe")
    )
    model <- mxModel(model, apoe_paths)
  }

  # Indirect effect algebra
  indirect <- mxAlgebra(
    a * b, name = "indirect"
  )
  total <- mxAlgebra(
    indirect + cprime, name = "total"
  )
  model <- mxModel(model, indirect, total)

  return(model)
}

# ============================================================
# 3. Fit models
# ============================================================
log_section("Part 3: Fitting Models")

all_results.lst <- list()

PRED_MAP <- list(
  FRS = "FRS_c",
  CVR_mimic = "CVR_c"
)

for (pred_name in names(PRED_MAP)) {
  pred_var <- PRED_MAP[[pred_name]]
  log_info("=== Predictor: %s (%s) ===",
           pred_name, pred_var)

  for (domain in COG_DOMAINS) {
    key <- paste0(domain, "_", pred_name)
    log_info("  --- %s ---", key)

    cog_cols.v <- paste0(
      "PHC_", domain, "_", VALID_TP
    )
    hvr_cols.v <- paste0("HVR_Z_", VALID_TP)

    # Determine available timepoints
    avail_tp.v <- VALID_TP[
      sapply(VALID_TP, function(tp) {
        cog_col <- paste0("PHC_", domain, "_", tp)
        hvr_col <- paste0("HVR_Z_", tp)
        yrs_col <- paste0("YRS_", tp)
        all(c(cog_col, hvr_col, yrs_col) %in%
              names(wide.dt)) &&
          sum(!is.na(wide.dt[[cog_col]])) >= 10 &&
          sum(!is.na(wide.dt[[hvr_col]])) >= 10
      })
    ]
    n_tp <- length(avail_tp.v)
    log_info("    Available timepoints: %d (%s)",
             n_tp, paste(avail_tp.v, collapse=","))

    if (n_tp < 2) {
      log_warn("    <2 timepoints, skipping")
      all_results.lst[[key]] <- list(
        converged = FALSE,
        error = "Insufficient timepoints"
      )
      next
    }

    cog_use.v <- paste0(
      "PHC_", domain, "_", avail_tp.v
    )
    hvr_use.v <- paste0(
      "HVR_Z_", avail_tp.v
    )

    # Temporarily adjust YRS_COLS for this model
    yrs_use.v <- paste0("YRS_", avail_tp.v)

    # Build data subset
    keep_cols.v <- c(
      hvr_use.v, cog_use.v, yrs_use.v,
      pred_var, "AGE_c", "EDUC_c",
      "SEX_num"
    )
    # SE columns for residual fixing
    se_use.v <- paste0(
      "PHC_", domain, "_SE_", avail_tp.v
    )
    keep_cols.v <- c(keep_cols.v, se_use.v)
    if (HAS_APOE) {
      keep_cols.v <- c(keep_cols.v, "APOE4_num")
    }
    keep_cols.v <- intersect(
      keep_cols.v, names(wide.dt)
    )

    model_data.dt <- wide.dt[
      !is.na(get(pred_var)), ..keep_cols.v
    ]
    log_info("    N = %d", nrow(model_data.dt))

    # Compute SE² for cognitive residuals
    avg_se2.v <- sapply(se_use.v, function(col) {
      if (col %in% names(model_data.dt)) {
        mean(
          model_data.dt[[col]]^2, na.rm = TRUE
        )
      } else 0.04
    })

    # Temporarily set YRS_COLS for the builder
    # (the builder reads from the global)
    old_yrs <- YRS_COLS
    YRS_COLS <<- yrs_use.v

    model <- tryCatch({
      build_linear_mediation.fn(
        cog_vars = cog_use.v,
        hvr_vars = hvr_use.v,
        predictor = pred_var,
        name = key,
        data = model_data.dt
      )
    }, error = function(e) {
      log_warn("    Build error: %s", e$message)
      NULL
    })
    YRS_COLS <<- old_yrs

    if (is.null(model)) {
      all_results.lst[[key]] <- list(
        converged = FALSE,
        error = "Model build failed"
      )
      next
    }

    # Attach data
    model <- mxModel(
      model,
      mxData(
        observed = as.data.frame(model_data.dt),
        type = "raw"
      )
    )

    # Fit with mxTryHard (retries with perturbed
    # starting values if initial attempt fails)
    fit <- tryCatch({
      mxTryHard(
        model, extraTries = 20,
        bestInitsOutput = FALSE,
        silent = TRUE
      )
    }, error = function(e) {
      log_warn("    Fit error: %s", e$message)
      NULL
    })

    if (is.null(fit)) {
      all_results.lst[[key]] <- list(
        converged = FALSE,
        error = "mxRun failed"
      )
      next
    }

    status <- fit$output$status$code
    converged <- status == 0
    log_info(
      "    Status: %d (%s)",
      status,
      if (converged) "converged" else "FAILED"
    )

    # Extract parameters
    params <- summary(fit)$parameters
    all_params <- summary(fit)$parameters

    # Key paths
    a_est <- params[
      params$name == "a", "Estimate"
    ]
    a_se <- params[
      params$name == "a", "Std.Error"
    ]
    b_est <- params[
      params$name == "b", "Estimate"
    ]
    b_se <- params[
      params$name == "b", "Std.Error"
    ]
    cp_est <- params[
      params$name == "cprime", "Estimate"
    ]
    cp_se <- params[
      params$name == "cprime", "Std.Error"
    ]

    # Indirect effect
    ind_est <- mxEval(indirect, fit)
    ind_se <- tryCatch({
      se <- mxSE(indirect, fit)
      as.numeric(se)
    }, error = function(e) NA_real_)

    a_p <- if (length(a_se) && !is.na(a_se))
      2 * pnorm(-abs(a_est / a_se)) else NA
    b_p <- if (length(b_se) && !is.na(b_se))
      2 * pnorm(-abs(b_est / b_se)) else NA
    ind_p <- if (!is.na(ind_se) && ind_se > 0)
      2 * pnorm(-abs(ind_est / ind_se)) else NA

    log_info(
      "    a-path: %.4f (SE=%.4f, p=%.4f)",
      a_est, a_se, a_p
    )
    log_info(
      "    b-path: %.4f (SE=%.4f, p=%.4f)",
      b_est, b_se, b_p
    )
    log_info(
      "    c'-path: %.4f (SE=%.4f)",
      cp_est, cp_se
    )
    log_info(
      "    indirect: %.6f (SE=%.6f)",
      as.numeric(ind_est),
      if (!is.na(ind_se)) ind_se else 0
    )

    all_results.lst[[key]] <- list(
      converged = converged,
      status = status,
      n = nrow(model_data.dt),
      n_timepoints = n_tp,
      all_params = all_params,
      a_path = list(
        est = a_est, se = a_se, p = a_p
      ),
      b_path = list(
        est = b_est, se = b_se, p = b_p
      ),
      cprime = list(
        est = cp_est, se = cp_se
      ),
      indirect = list(
        est = as.numeric(ind_est),
        se = ind_se, p = ind_p
      ),
      fit_summary = summary(fit)
    )
  }
}

# ============================================================
# 4. Summary
# ============================================================
log_section("KEY FINDINGS (A-/CU LGCM)")

for (pred_name in names(PRED_MAP)) {
  log_info("=== %s ===", pred_name)
  for (domain in COG_DOMAINS) {
    key <- paste0(domain, "_", pred_name)
    r <- all_results.lst[[key]]
    if (is.null(r) || !r$converged) {
      log_info(
        "  %s: did not converge", domain
      )
      next
    }
    a_sig <- ifelse(
      !is.na(r$a_path$p) & r$a_path$p < 0.05,
      " *", ""
    )
    b_sig <- ifelse(
      !is.na(r$b_path$p) & r$b_path$p < 0.05,
      " *", ""
    )
    log_info(
      "  %s (N=%d): a=%.4f (p=%.4f%s), b=%.4f (p=%.4f%s), ind=%.6f",
      domain, r$n,
      r$a_path$est, r$a_path$p, a_sig,
      r$b_path$est, r$b_path$p, b_sig,
      r$indirect$est
    )
  }
}

# ============================================================
# 5. Save (only if at least one model converged)
# ============================================================
log_section("Saving Results")

n_converged <- sum(sapply(
  all_results.lst,
  function(r) isTRUE(r$converged)
))
log_info("%d/%d models converged",
         n_converged, length(all_results.lst))

results.lst <- list(
  results = all_results.lst,
  sample_info = lme_res.lst$sample_info,
  specification = "linear_mediation_YRS",
  sample = "Amyloid-negative, baseline CU"
)

if (n_converged > 0) {
  ensure_directory(dirname(output.path))
  write_rds_safe(
    results.lst, output.path,
    "A-/CU LGCM results"
  )
} else {
  log_warn("No models converged; not saving.")
}

log_script_end(
  "13b_lgcm_amyloid_negative.R", success = TRUE
)
