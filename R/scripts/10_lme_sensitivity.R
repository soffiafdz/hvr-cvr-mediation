#!/usr/bin/env Rscript

# ============================================================
# 07_lme_sensitivity.R - LME SENSITIVITY ANALYSES
# ============================================================
# Consolidates all LME sensitivity analyses:
#
# PART 1: Raw Brain Measure Models (from 06c + 06d)
#   - HVR_raw x {FRS, CVR_mimic}  (full sample)
#   - HC_raw  x {FRS, CVR_mimic} + TIV  (full sample)
#   - Includes YRS^2 term (consistent with 08/09)
#   - Tests robustness to normative modeling approach
#   Saves: lme_raw_sensitivity_results.rds
#
# PART 2: Attrition/Dropout Analysis (from 08 S1)
#   - Baseline characteristic comparisons
#   - Logistic regression predicting dropout
#   - Model-based attrition sensitivity (completers vs full)
#   Saves: attrition_analysis.rds
#
# PART 3: EDT Linearity Test (from 08 S2)
#   - Linear vs quadratic EDT models via LRT
#   Saves: edt_linearity_test.rds
#
# PART 4: EDT Early/Middle Restriction (from 08 S3)
#   - Excludes late EDT tertile
#   Saves: edt_early_middle_sensitivity.rds
#
# INPUTS:
#   - data/derivatives/lme/cohort_hvr.rds
#   - data/derivatives/lme/cvr_mimic_scores.rds (optional)
#
# OUTPUTS:
#   - models/results/lme/lme_raw_sensitivity_results.rds
#   - models/results/lme/attrition_analysis.rds
#   - models/results/lgcm/edt_linearity_test.rds
#   - models/results/lme/edt_early_middle_sensitivity.rds
# ============================================================

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
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

log_script_start("10_lme_sensitivity.R")
config <- load_config()
validate_config(config)
validate_packages(c("lme4", "lmerTest", "data.table"))
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "lme_sensitivity", default = FALSE
)
output.path <- get_data_path(
  "models", "lme_raw_sensitivity"
)
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE to rerun.")
  log_script_end("10_lme_sensitivity.R", success = TRUE)
  quit(status = 0)
}

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
COGNITIVE_DOMAINS <- get_parameter("cognitive_domains")

LME_ESTIMATOR <- get_script_setting(
  "lme", "estimator", default = "REML"
)
LME_OPTIMIZER <- get_script_setting(
  "lme", "optimizer", default = "bobyqa"
)
LME_MAXITER <- get_script_setting(
  "lme", "max_iter", default = 100000
)

lme_results_dir.path <- get_data_path(
  "models", "lme_results_dir"
)
lgcm_results_dir.path <- get_data_path(
  "models", "lgcm_results_dir"
)
ensure_directory(lme_results_dir.path)
ensure_directory(lgcm_results_dir.path)

# Output paths
raw_sens.path <- file.path(
  lme_results_dir.path,
  "lme_raw_sensitivity_results.rds"
)
attrition.path <- file.path(
  lme_results_dir.path, "attrition_analysis.rds"
)
linearity.path <- file.path(
  lgcm_results_dir.path, "edt_linearity_test.rds"
)
edt_sens.path <- file.path(
  lme_results_dir.path,
  "edt_early_middle_sensitivity.rds"
)

# ============================================================
# Load Data
# ============================================================
log_section("Loading Data")

cohort.path <- get_data_path("derivatives", "lme_cohort_hvr")
check_files_exist(cohort.path)
cohort.dt <- read_rds_safe(cohort.path, "HVR cohort")

log_info(
  "Loaded HVR cohort: %d subjects, %d observations",
  cohort.dt[, uniqueN(PTID)], nrow(cohort.dt)
)

# --- Check APOE4 availability ---
HAS_APOE <- "APOE4" %in% names(cohort.dt) &&
  !all(is.na(cohort.dt$APOE4))
if (HAS_APOE) {
  log_info("APOE4 available")
} else {
  log_warn("APOE4 not available - excluding from models")
}

# --- Compute baseline age ---
log_section("Computing Baseline Age")
setorder(cohort.dt, PTID, EXAMDATE)
baseline_age.dt <- cohort.dt[, .(Age_bl = AGE[1]), by = PTID]
cohort.dt <- merge(
  cohort.dt, baseline_age.dt,
  by = "PTID", all.x = TRUE
)
log_info(
  "Age_bl range: %.1f - %.1f years",
  min(cohort.dt$Age_bl), max(cohort.dt$Age_bl)
)

# --- Try loading CVR_mimic scores ---
HAS_CVR_MIMIC <- FALSE
cvr_mimic.path <- tryCatch(
  get_data_path("models", "cvr_mimic_scores"),
  error = function(e) NULL
)
if (!is.null(cvr_mimic.path) && file.exists(cvr_mimic.path)) {
  cvr_scores.dt <- tryCatch({
    scores <- read_rds_safe(
      cvr_mimic.path, "CVR MIMIC scores"
    )
    if (is.data.table(scores) &&
        "CVR_mimic" %in% names(scores) &&
        "PTID" %in% names(scores)) {
      scores
    } else if (is.list(scores) &&
               "scores" %in% names(scores)) {
      as.data.table(scores$scores)
    } else {
      NULL
    }
  }, error = function(e) {
    log_warn("Failed to load CVR_mimic: %s", e$message)
    NULL
  })

  if (!is.null(cvr_scores.dt) &&
      "CVR_mimic" %in% names(cvr_scores.dt)) {
    if (!"CVR_mimic" %in% names(cohort.dt)) {
      # Scores file may use RID (from script 07)
      # or PTID (if pre-merged upstream)
      if ("PTID" %in% names(cvr_scores.dt)) {
        merge_key <- "PTID"
      } else if ("RID" %in% names(cvr_scores.dt) &&
                 "RID" %in% names(cohort.dt)) {
        merge_key <- "RID"
      } else {
        merge_key <- NULL
        log_warn("No common key to merge CVR_mimic")
      }
      if (!is.null(merge_key)) {
        cohort.dt <- merge(
          cohort.dt,
          cvr_scores.dt[
            , .SD, .SDcols = c(merge_key, "CVR_mimic")
          ],
          by = merge_key, all.x = TRUE
        )
      }
    }
    n_cvr <- sum(!is.na(cohort.dt$CVR_mimic))
    if (n_cvr >= 50) {
      HAS_CVR_MIMIC <- TRUE
      log_info(
        "CVR_mimic available: %d non-missing obs", n_cvr
      )
    } else {
      log_warn(
        "CVR_mimic has only %d non-missing obs, skipping",
        n_cvr
      )
    }
  }
} else {
  log_warn(
    "CVR_mimic scores file not found - FRS only"
  )
}

# ============================================================
# Helper Functions
# ============================================================

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

    if (is_singular) log_warn("    Singular fit")
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

extract_3way_interaction.fn <- function(
    model.res, cvr_var, brain_var
) {
  if (is.null(model.res$fit) || !model.res$converged) {
    return(list(
      term = NA, beta = NA, se = NA,
      t_value = NA, p_value = NA
    ))
  }

  coefs.mat <- model.res$coefficients
  # All orderings of the 3-way interaction
  vars.v <- c("YRS_from_bl", cvr_var, brain_var)
  perms.lst <- list(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  patterns.v <- vapply(perms.lst, function(p) {
    paste(vars.v[p], collapse = ":")
  }, character(1))
  pattern <- paste(patterns.v, collapse = "|")

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
# PART 1: RAW BRAIN MEASURE MODELS
# ============================================================
# Loop: BRAIN_MEASURES x CVR_MEASURES
# HVR_raw x {FRS, CVR_mimic}  (full sample)
# HC_raw  x {FRS, CVR_mimic} + TIV  (full sample)
# No YRS^2 term. No sex stratification.
# ============================================================
log_section("PART 1: Raw Brain Measure Models")

# Brain measures configuration
BRAIN_MEASURES.lst <- list(
  HVR_raw = list(
    var = "HVR_raw",
    label = "HVR raw (sample-standardized)",
    needs_tiv = FALSE
  ),
  HC_raw = list(
    var = "HC_raw",
    label = "HC raw + TIV covariate",
    needs_tiv = TRUE
  )
)

# CVR measures: FRS always; CVR_mimic if available
CVR_MEASURES.lst <- list(FRS = "FRS")
if (HAS_CVR_MIMIC) {
  CVR_MEASURES.lst[["CVR_mimic"]] <- "CVR_mimic"
}
log_info(
  "CVR measures: %s",
  paste(names(CVR_MEASURES.lst), collapse = ", ")
)

# Build covariate strings
BASE_COV <- if (HAS_APOE) {
  "Age_bl + EDUC + SEX + APOE4"
} else {
  "Age_bl + EDUC + SEX"
}
TIV_COV <- if (HAS_APOE) {
  "Age_bl + EDUC + SEX + TIV + APOE4"
} else {
  "Age_bl + EDUC + SEX + TIV"
}

raw_results.lst <- list()
raw_summary_rows.lst <- list()

for (brain_name in names(BRAIN_MEASURES.lst)) {
  measure <- BRAIN_MEASURES.lst[[brain_name]]
  brain_var <- measure$var
  needs_tiv <- measure$needs_tiv

  # Check brain variable availability
  if (!brain_var %in% names(cohort.dt) ||
      sum(!is.na(cohort.dt[[brain_var]])) < 50) {
    log_warn(
      "%s not available in data, skipping", brain_var
    )
    next
  }

  # TIV check for HC models
  if (needs_tiv) {
    if (!"TIV" %in% names(cohort.dt) ||
        all(is.na(cohort.dt$TIV))) {
      log_warn("TIV not available for %s, skipping", brain_var)
      next
    }
  }

  covariates_str <- if (needs_tiv) TIV_COV else BASE_COV

  for (cvr_name in names(CVR_MEASURES.lst)) {
    cvr_var <- CVR_MEASURES.lst[[cvr_name]]

    log_info("")
    log_info(
      "=== %s x %s (raw sensitivity) ===",
      brain_name, cvr_name
    )
    log_info(
      "Model: Cog ~ YRS * %s * %s + YRS^2 * %s * %s + %s + (YRS | PTID)",
      cvr_var, brain_var, cvr_var, brain_var,
      covariates_str
    )

    # Subset to complete cases for this combination
    model_data.dt <- cohort.dt[
      !is.na(get(brain_var)) &
        !is.na(get(cvr_var)) &
        !is.na(YRS_from_bl)
    ]
    if (needs_tiv) {
      model_data.dt <- model_data.dt[!is.na(TIV)]
    }

    combo_key <- paste(brain_name, cvr_name, sep = "_")
    raw_results.lst[[combo_key]] <- list()

    for (domain in COGNITIVE_DOMAINS) {
      domain_data.dt <- model_data.dt[
        !is.na(get(domain))
      ]
      if (nrow(domain_data.dt) < 50) {
        log_warn(
          "  %s: insufficient data (N=%d), skip",
          domain, nrow(domain_data.dt)
        )
        next
      }

      formula_str <- sprintf(
        paste0(
          "%s ~ YRS_from_bl * %s * %s",
          " + I(YRS_from_bl^2) * %s * %s",
          " + %s + (YRS_from_bl | PTID)"
        ),
        domain, cvr_var, brain_var,
        cvr_var, brain_var,
        covariates_str
      )

      desc_str <- sprintf(
        "%s: YRS x %s x %s",
        domain, cvr_name, brain_name
      )
      if (needs_tiv) {
        desc_str <- paste0(desc_str, " + TIV")
      }

      model.res <- fit_lme_growth.fn(
        as.formula(formula_str),
        domain_data.dt, desc_str
      )

      int.lst <- extract_3way_interaction.fn(
        model.res, cvr_var, brain_var
      )

      raw_results.lst[[combo_key]][[domain]] <- list(
        model = model.res,
        interaction = int.lst,
        formula = formula_str,
        n_subjects = domain_data.dt[, uniqueN(PTID)],
        n_obs = nrow(domain_data.dt)
      )

      if (!is.na(int.lst$p_value)) {
        sig <- ifelse(int.lst$p_value < 0.05, " *", "")
        log_info(
          "  %s: B=%.4f (SE=%.4f), p=%.4f%s",
          domain, int.lst$beta, int.lst$se,
          int.lst$p_value, sig
        )
      }

      # Summary row
      raw_summary_rows.lst[[
        length(raw_summary_rows.lst) + 1
      ]] <- data.table(
        Brain_Measure = brain_name,
        CVR_Measure = cvr_name,
        Analysis = "raw_sensitivity",
        Domain = domain,
        N_subjects = domain_data.dt[, uniqueN(PTID)],
        N_obs = nrow(domain_data.dt),
        Beta = int.lst$beta,
        SE = int.lst$se,
        t_value = int.lst$t_value,
        p_value = int.lst$p_value
      )
    }
  }
}

# Build summary table
if (length(raw_summary_rows.lst) > 0) {
  raw_summary.dt <- rbindlist(raw_summary_rows.lst)
  raw_summary.dt[, Significant := p_value < 0.05]

  log_section("Raw Sensitivity Summary")
  log_info(
    "Total models fitted: %d", nrow(raw_summary.dt)
  )
  log_info(
    "Significant 3-way interactions: %d",
    sum(raw_summary.dt$Significant, na.rm = TRUE)
  )
} else {
  raw_summary.dt <- data.table()
  log_warn("No raw sensitivity models were fitted")
}

# Key findings
log_section("PART 1 KEY FINDINGS")
for (combo_key in names(raw_results.lst)) {
  log_info("")
  log_info("--- %s ---", combo_key)
  for (domain in COGNITIVE_DOMAINS) {
    res <- raw_results.lst[[combo_key]][[domain]]
    if (!is.null(res) && !is.na(res$interaction$p_value)) {
      int <- res$interaction
      sig <- ifelse(int$p_value < 0.05, "***", "")
      log_info(
        "  %s: B=%.4f, p=%.4f %s",
        domain, int$beta, int$p_value, sig
      )
    }
  }
}

# Save Part 1 results
log_section("Saving Part 1 Results")

raw_output.lst <- list(
  results = raw_results.lst,
  summary_table = raw_summary.dt,
  cvr_measures = names(CVR_MEASURES.lst),
  brain_measures = names(BRAIN_MEASURES.lst),

  methodology = list(
    purpose = paste0(
      "Sensitivity: raw brain measures, ",
      "test robustness to normative modeling"
    ),
    model_template = paste0(
      "Cog ~ YRS * CVR * Brain",
      " + YRS^2 * CVR * Brain",
      " + covariates + (YRS | PTID)"
    ),
    note = paste0(
      "Full sample only, no sex stratification.",
      " HC models include TIV."
    ),
    has_apoe = HAS_APOE,
    has_cvr_mimic = HAS_CVR_MIMIC
  ),

  config = list(
    estimator = LME_ESTIMATOR,
    optimizer = LME_OPTIMIZER,
    timestamp = Sys.time()
  )
)

write_rds_safe(
  raw_output.lst, raw_sens.path,
  "LME raw sensitivity results"
)

# ============================================================
# PART 2: ATTRITION / DROPOUT ANALYSIS
# ============================================================
# Completers (>=4 visits) vs Dropouts (<4 visits)
# 2a. Baseline characteristic comparison
# 2b. Logistic regression predicting dropout
# 2c. Model-based attrition sensitivity
# ============================================================
log_section("PART 2: Attrition/Dropout Analysis")

# Count visits per subject
visit_counts.dt <- cohort.dt[, .N, by = PTID]
setnames(visit_counts.dt, "N", "n_visits")

COMPLETER_THRESHOLD <- 4
visit_counts.dt[
  , completer := fifelse(
    n_visits >= COMPLETER_THRESHOLD, 1L, 0L
  )
]

log_info(
  "Completer definition: >= %d visits",
  COMPLETER_THRESHOLD
)
log_info(
  "  Completers: %d subjects (%.1f%%)",
  visit_counts.dt[completer == 1, .N],
  100 * visit_counts.dt[completer == 1, .N] /
    nrow(visit_counts.dt)
)
log_info(
  "  Dropouts: %d subjects (%.1f%%)",
  visit_counts.dt[completer == 0, .N],
  100 * visit_counts.dt[completer == 0, .N] /
    nrow(visit_counts.dt)
)

# Merge completer status
cohort.dt <- merge(
  cohort.dt,
  visit_counts.dt[, .(PTID, n_visits, completer)],
  by = "PTID"
)

# Baseline characteristics (first visit per subject)
baseline.dt <- cohort.dt[
  order(PTID, YRS_from_bl)
][, .SD[1], by = PTID]

log_info("Baseline observations: %d", nrow(baseline.dt))

# --- 2a. Compare baseline characteristics ---
log_info("")
log_info("--- Comparing Baseline Characteristics ---")

# Brain measures for attrition (z-scores)
BRAIN_Z.lst <- list(
  HVR_z = list(
    var = "HVR_z", cov_extra = "",
    label = "HVR z-score"
  ),
  HC_z = list(
    var = "HC_z", cov_extra = "",
    label = "HC z-score"
  )
)

brain_vars.v <- sapply(BRAIN_Z.lst, function(m) m$var)
continuous_vars.v <- c(
  "AGE", "EDUC", "FRS", brain_vars.v
)
categorical_vars.v <- "SEX"
if (HAS_APOE) {
  categorical_vars.v <- c(
    categorical_vars.v, "APOE4"
  )
}

# Filter to existing columns
continuous_vars.v <- continuous_vars.v[
  continuous_vars.v %in% names(baseline.dt)
]
categorical_vars.v <- categorical_vars.v[
  categorical_vars.v %in% names(baseline.dt)
]

comparison_results.lst <- list()

for (v in continuous_vars.v) {
  comp_vals.v <- baseline.dt[completer == 1, get(v)]
  drop_vals.v <- baseline.dt[completer == 0, get(v)]

  comp_vals.v <- comp_vals.v[!is.na(comp_vals.v)]
  drop_vals.v <- drop_vals.v[!is.na(drop_vals.v)]

  if (length(comp_vals.v) >= 5 &&
      length(drop_vals.v) >= 5) {
    tt.res <- t.test(comp_vals.v, drop_vals.v)

    comparison_results.lst[[v]] <- data.table(
      Variable = v,
      Completers_mean = mean(comp_vals.v),
      Completers_sd = sd(comp_vals.v),
      Dropouts_mean = mean(drop_vals.v),
      Dropouts_sd = sd(drop_vals.v),
      Diff = mean(comp_vals.v) - mean(drop_vals.v),
      t_statistic = tt.res$statistic,
      p_value = tt.res$p.value,
      Test = "t-test"
    )

    sig <- ifelse(tt.res$p.value < 0.05, "*", "")
    log_info(
      "  %s: Comp=%.2f (%.2f), Drop=%.2f (%.2f), p=%.4f %s",
      v,
      mean(comp_vals.v), sd(comp_vals.v),
      mean(drop_vals.v), sd(drop_vals.v),
      tt.res$p.value, sig
    )
  } else {
    log_warn(
      "  %s: Insufficient data for comparison", v
    )
  }
}

# Categorical: chi-square
for (v in categorical_vars.v) {
  if (v %in% names(baseline.dt)) {
    tbl <- table(baseline.dt[[v]], baseline.dt$completer)
    if (all(dim(tbl) >= 2)) {
      chi.res <- chisq.test(tbl)

      if (v == "SEX") {
        completers_pct <- 100 * sum(
          baseline.dt[completer == 1, SEX == "Male"]
        ) / baseline.dt[completer == 1, .N]
        dropouts_pct <- 100 * sum(
          baseline.dt[completer == 0, SEX == "Male"]
        ) / baseline.dt[completer == 0, .N]
      } else if (v == "APOE4") {
        completers_pct <- 100 * sum(
          baseline.dt[completer == 1, APOE4 == TRUE]
        ) / baseline.dt[completer == 1, .N]
        dropouts_pct <- 100 * sum(
          baseline.dt[completer == 0, APOE4 == TRUE]
        ) / baseline.dt[completer == 0, .N]
      } else {
        completers_pct <- NA
        dropouts_pct <- NA
      }

      comparison_results.lst[[v]] <- data.table(
        Variable = v,
        Completers_mean = completers_pct,
        Completers_sd = NA,
        Dropouts_mean = dropouts_pct,
        Dropouts_sd = NA,
        Diff = if (!is.na(completers_pct)) {
          completers_pct - dropouts_pct
        } else {
          NA
        },
        t_statistic = chi.res$statistic,
        p_value = chi.res$p.value,
        Test = "chi-square"
      )

      sig <- ifelse(chi.res$p.value < 0.05, "*", "")
      log_info(
        paste0(
          "  %s: chi-sq=%.2f, p=%.4f %s",
          " (Comp: %.1f%%, Drop: %.1f%%)"
        ),
        v, chi.res$statistic, chi.res$p.value, sig,
        completers_pct, dropouts_pct
      )
    }
  }
}

comparison.dt <- rbindlist(
  comparison_results.lst, fill = TRUE
)

# --- 2b. Logistic regression: dropout prediction ---
log_info("")
log_info("--- Logistic Regression: Dropout Prediction ---")

baseline.dt[, dropout := 1L - completer]

available_predictors.v <- c(
  "FRS", brain_vars.v, "AGE", "SEX"
)
if (HAS_APOE) {
  available_predictors.v <- c(
    available_predictors.v, "APOE4"
  )
}
available_predictors.v <- available_predictors.v[
  available_predictors.v %in% names(baseline.dt)
]

# Remove predictors with too few non-NA values
for (pred in available_predictors.v) {
  if (sum(!is.na(baseline.dt[[pred]])) < 50) {
    available_predictors.v <- setdiff(
      available_predictors.v, pred
    )
  }
}

logistic_results.lst <- list()

if (length(available_predictors.v) >= 2) {
  logit_formula <- as.formula(paste(
    "dropout ~",
    paste(available_predictors.v, collapse = " + ")
  ))

  dropout.fit <- tryCatch({
    glm(
      logit_formula, data = baseline.dt,
      family = binomial()
    )
  }, error = function(e) {
    log_warn("Logistic regression failed: %s", e$message)
    NULL
  })

  if (!is.null(dropout.fit)) {
    model_summary <- summary(dropout.fit)
    coef.dt <- as.data.table(
      coef(model_summary), keep.rownames = "term"
    )
    setnames(coef.dt, c(
      "term", "estimate", "std_error",
      "z_value", "p_value"
    ))

    or_ci.mat <- tryCatch({
      exp(cbind(
        OR = coef(dropout.fit),
        confint(dropout.fit)
      ))
    }, error = function(e) {
      log_warn("CI calculation failed: %s", e$message)
      exp(cbind(OR = coef(dropout.fit)))
    })

    or.dt <- as.data.table(
      or_ci.mat, keep.rownames = "term"
    )

    log_info(
      "Logistic regression results (predicting dropout):"
    )
    for (i in seq_len(nrow(coef.dt))) {
      term <- coef.dt$term[i]
      est <- coef.dt$estimate[i]
      p <- coef.dt$p_value[i]
      or <- exp(est)
      sig <- ifelse(p < 0.05, "*", "")
      log_info(
        "  %s: OR=%.2f, p=%.4f %s",
        term, or, p, sig
      )
    }

    any_significant <- any(
      coef.dt$p_value[-1] < 0.05
    )

    interpretation <- if (any_significant) {
      sig_terms.v <- coef.dt$term[
        coef.dt$p_value < 0.05 &
          coef.dt$term != "(Intercept)"
      ]
      paste0(
        sprintf(
          paste0(
            "CAUTION: Significant predictors",
            " of dropout: %s. "
          ),
          paste(sig_terms.v, collapse = ", ")
        ),
        "Findings may be affected by attrition bias."
      )
    } else {
      paste0(
        "No significant predictors of dropout. ",
        "Attrition bias unlikely to affect findings."
      )
    }

    log_info("")
    log_info("INTERPRETATION: %s", interpretation)

    logistic_results.lst <- list(
      formula = as.character(logit_formula),
      coefficients = coef.dt,
      odds_ratios = or.dt,
      aic = AIC(dropout.fit),
      n_completers = baseline.dt[completer == 1, .N],
      n_dropouts = baseline.dt[completer == 0, .N],
      any_significant = any_significant,
      interpretation = interpretation
    )
  } else {
    logistic_results.lst <- list(
      error = "Model fitting failed",
      n_completers = baseline.dt[completer == 1, .N],
      n_dropouts = baseline.dt[completer == 0, .N]
    )
  }
} else {
  log_warn(
    "Insufficient predictors for logistic regression"
  )
  logistic_results.lst <- list(
    error = "Insufficient predictors"
  )
}

# --- 2c. Model-based attrition sensitivity ---
log_info("")
log_info("--- Model-Based Attrition Sensitivity ---")
log_info(
  "Comparing 3-way models: Full sample vs Completers"
)

COV_LONG <- if (HAS_APOE) {
  "Age_bl + EDUC + SEX + APOE4"
} else {
  "Age_bl + EDUC + SEX"
}

completers.dt <- cohort.dt[
  PTID %in% visit_counts.dt[completer == 1, PTID]
]
log_info(
  "Completers subsample: %d subjects, %d obs",
  completers.dt[, uniqueN(PTID)],
  nrow(completers.dt)
)

model_attrition_results.lst <- list()

for (measure_name in names(BRAIN_Z.lst)) {
  measure <- BRAIN_Z.lst[[measure_name]]
  brain_var <- measure$var

  if (!brain_var %in% names(cohort.dt) ||
      sum(!is.na(cohort.dt[[brain_var]])) < 100) {
    log_warn("%s: Not available, skipping", measure_name)
    next
  }

  log_info("")
  log_info(
    "=== Attrition Sensitivity: %s ===",
    measure$label
  )

  for (domain in COGNITIVE_DOMAINS) {
    log_info("  --- %s ---", domain)

    full_data.dt <- cohort.dt[
      !is.na(get(domain)) & !is.na(FRS) &
        !is.na(get(brain_var)) & !is.na(YRS_from_bl)
    ]
    comp_data.dt <- completers.dt[
      !is.na(get(domain)) & !is.na(FRS) &
        !is.na(get(brain_var)) & !is.na(YRS_from_bl)
    ]

    n_full <- nrow(full_data.dt)
    n_comp <- nrow(comp_data.dt)
    n_subj_full <- full_data.dt[, uniqueN(PTID)]
    n_subj_comp <- comp_data.dt[, uniqueN(PTID)]

    if (n_comp < 100) {
      log_warn(
        "    Insufficient completers (N=%d)", n_comp
      )
      next
    }

    # Model with quadratic time (matches 06a/06b)
    f_str <- sprintf(
      paste0(
        "%s ~ YRS_from_bl * FRS * %s",
        " + I(YRS_from_bl^2) * FRS * %s",
        " + %s + (YRS_from_bl | PTID)"
      ),
      domain, brain_var, brain_var, COV_LONG
    )
    f_model <- as.formula(f_str)

    m_full.fit <- tryCatch({
      lmer(
        f_model, data = full_data.dt,
        REML = TRUE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 100000)
        )
      )
    }, error = function(e) {
      log_warn("    Full model failed: %s", e$message)
      NULL
    })

    m_comp.fit <- tryCatch({
      lmer(
        f_model, data = comp_data.dt,
        REML = TRUE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 100000)
        )
      )
    }, error = function(e) {
      log_warn(
        "    Completers model failed: %s", e$message
      )
      NULL
    })

    if (!is.null(m_full.fit) && !is.null(m_comp.fit)) {
      coef_full.mat <- summary(m_full.fit)$coefficients
      coef_comp.mat <- summary(m_comp.fit)$coefficients

      int_pattern <- paste0(
        "YRS_from_bl:FRS:", brain_var, "|",
        "YRS_from_bl:", brain_var, ":FRS|",
        "FRS:YRS_from_bl:", brain_var, "|",
        "FRS:", brain_var, ":YRS_from_bl|",
        brain_var, ":YRS_from_bl:FRS|",
        brain_var, ":FRS:YRS_from_bl"
      )
      int_term <- grep(
        int_pattern, rownames(coef_full.mat),
        value = TRUE
      )[1]

      if (!is.na(int_term) &&
          int_term %in% rownames(coef_comp.mat)) {
        beta_full <- coef_full.mat[int_term, "Estimate"]
        se_full <- coef_full.mat[int_term, "Std. Error"]
        p_full <- coef_full.mat[int_term, "Pr(>|t|)"]

        beta_comp <- coef_comp.mat[int_term, "Estimate"]
        se_comp <- coef_comp.mat[int_term, "Std. Error"]
        p_comp <- coef_comp.mat[int_term, "Pr(>|t|)"]

        same_direction <- sign(beta_full) ==
          sign(beta_comp)
        same_significance <- (p_full < 0.05) ==
          (p_comp < 0.05)
        beta_change_pct <- (beta_comp - beta_full) /
          abs(beta_full) * 100

        result_key <- paste(
          measure_name, domain, sep = "_"
        )
        model_attrition_results.lst[[
          result_key
        ]] <- data.table(
          Measure = measure_name,
          Domain = domain,
          N_obs_full = n_full,
          N_subj_full = n_subj_full,
          N_obs_comp = n_comp,
          N_subj_comp = n_subj_comp,
          Beta_full = beta_full,
          SE_full = se_full,
          P_full = p_full,
          Beta_completers = beta_comp,
          SE_completers = se_comp,
          P_completers = p_comp,
          Beta_change_pct = beta_change_pct,
          Same_direction = same_direction,
          Same_significance = same_significance,
          Consistent = same_direction &
            same_significance
        )

        sig_f <- ifelse(p_full < 0.05, "*", "")
        sig_c <- ifelse(p_comp < 0.05, "*", "")
        consist <- ifelse(
          same_direction & same_significance,
          "CONSISTENT", "DIFFERS"
        )

        log_info(
          "    Full:       B=%.6f (SE=%.6f), p=%.4f%s",
          beta_full, se_full, p_full, sig_f
        )
        log_info(
          "    Completers: B=%.6f (SE=%.6f), p=%.4f%s",
          beta_comp, se_comp, p_comp, sig_c
        )
        log_info(
          "    Change: %.1f%%, %s",
          beta_change_pct, consist
        )
      } else {
        log_warn("    3-way interaction term not found")
      }
    }
  }
}

# Summarize model-based attrition
if (length(model_attrition_results.lst) > 0) {
  model_attrition.dt <- rbindlist(
    model_attrition_results.lst
  )
  all_consistent <- all(model_attrition.dt$Consistent)

  model_attrition_interpretation <- if (all_consistent) {
    paste0(
      "Model-based sensitivity: FRS x Brain x YRS ",
      "effects CONSISTENT between full sample and ",
      "completers. Attrition bias unlikely."
    )
  } else {
    inconsistent.dt <- model_attrition.dt[
      Consistent == FALSE
    ]
    incon_str <- paste(
      sprintf(
        "%s/%s",
        inconsistent.dt$Measure,
        inconsistent.dt$Domain
      ),
      collapse = ", "
    )
    sprintf(
      paste0(
        "Model-based sensitivity: INCONSISTENT ",
        "for %s. Consider potential attrition bias."
      ),
      incon_str
    )
  }

  log_info("")
  log_info(
    "MODEL-BASED ATTRITION: %s",
    model_attrition_interpretation
  )
} else {
  model_attrition.dt <- NULL
  model_attrition_interpretation <- paste0(
    "Model-based sensitivity: Could not be computed"
  )
  all_consistent <- NA
}

# Save attrition results
attrition_results.lst <- list(
  comparison_table = comparison.dt,
  logistic_model = logistic_results.lst,
  model_sensitivity = model_attrition.dt,
  model_sensitivity_consistent = all_consistent,
  model_sensitivity_interpretation =
    model_attrition_interpretation,
  n_completers = baseline.dt[completer == 1, .N],
  n_dropouts = baseline.dt[completer == 0, .N],
  completer_threshold = COMPLETER_THRESHOLD,
  baseline_interpretation =
    logistic_results.lst$interpretation,
  model_interpretation =
    model_attrition_interpretation
)

write_rds_safe(
  attrition_results.lst, attrition.path,
  description = "Attrition analysis results"
)
log_info("Saved attrition analysis to %s", attrition.path)

# ============================================================
# PART 3: EDT LINEARITY SENSITIVITY TEST
# ============================================================
# Linear vs quadratic EDT models via LRT
# ============================================================
log_section("PART 3: EDT Linearity Sensitivity Test")

COVARIATES_CROSS <- if (HAS_APOE) {
  "AGE + EDUC + SEX + APOE4"
} else {
  "AGE + EDUC + SEX"
}

HAS_EDT <- "EDT" %in% names(cohort.dt) &&
  sum(!is.na(cohort.dt$EDT)) > 100

if (HAS_EDT) {
  log_info(
    "EDT available: %d obs with EDT values",
    sum(!is.na(cohort.dt$EDT))
  )

  edt_mean <- mean(cohort.dt$EDT, na.rm = TRUE)
  edt_sd <- sd(cohort.dt$EDT, na.rm = TRUE)
  cohort.dt[
    , EDT_scaled := (EDT - edt_mean) / edt_sd
  ]
  cohort.dt[, EDT_scaled_sq := EDT_scaled^2]

  log_info(
    "EDT scaling: mean=%.2f, SD=%.2f",
    edt_mean, edt_sd
  )

  linearity_results.lst <- list()

  for (domain in COGNITIVE_DOMAINS) {
    log_info("")
    log_info(
      "--- %s: Linear vs Quadratic EDT ---", domain
    )

    domain_data.dt <- cohort.dt[
      !is.na(get(domain)) & !is.na(EDT_scaled)
    ]
    n_obs <- nrow(domain_data.dt)
    n_subj <- domain_data.dt[, uniqueN(PTID)]

    if (n_obs < 100) {
      log_warn(
        "Insufficient data for %s (N=%d)", domain,
        n_obs
      )
      next
    }

    log_info(
      "  Data: %d observations, %d subjects",
      n_obs, n_subj
    )

    f_linear <- as.formula(sprintf(
      "%s ~ EDT_scaled + %s + (1 | PTID)",
      domain, COVARIATES_CROSS
    ))
    f_quad <- as.formula(sprintf(
      "%s ~ EDT_scaled + EDT_scaled_sq + %s + (1 | PTID)",
      domain, COVARIATES_CROSS
    ))

    m_linear.fit <- tryCatch({
      lmer(
        f_linear, data = domain_data.dt,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
    }, error = function(e) {
      log_warn("Linear model failed: %s", e$message)
      NULL
    })

    m_quad.fit <- tryCatch({
      lmer(
        f_quad, data = domain_data.dt,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
    }, error = function(e) {
      log_warn("Quadratic model failed: %s", e$message)
      NULL
    })

    if (!is.null(m_linear.fit) && !is.null(m_quad.fit)) {
      lrt.res <- anova(m_linear.fit, m_quad.fit)

      linear_aic <- AIC(m_linear.fit)
      quad_aic <- AIC(m_quad.fit)
      lrt_chisq <- lrt.res$Chisq[2]
      lrt_df <- lrt.res$Df[2]
      lrt_p <- lrt.res$`Pr(>Chisq)`[2]

      quad_coef <- fixef(m_quad.fit)["EDT_scaled_sq"]
      quad_se <- sqrt(
        diag(vcov(m_quad.fit))
      )["EDT_scaled_sq"]

      quad_significant <- lrt_p < 0.05

      linearity_results.lst[[domain]] <- data.table(
        Domain = domain,
        N_obs = n_obs,
        N_subjects = n_subj,
        Linear_AIC = linear_aic,
        Quadratic_AIC = quad_aic,
        AIC_diff = linear_aic - quad_aic,
        LRT_chisq = lrt_chisq,
        LRT_df = lrt_df,
        LRT_p = lrt_p,
        Quad_coef = quad_coef,
        Quad_SE = quad_se,
        Quadratic_significant = quad_significant,
        Interpretation = ifelse(
          quad_significant,
          paste0(
            "Non-linear EDT effect detected",
            " - consider quadratic term"
          ),
          "Linear EDT assumption supported"
        )
      )

      sig <- ifelse(quad_significant, "***", "")
      log_info(
        "  Linear AIC: %.1f, Quadratic AIC: %.1f",
        linear_aic, quad_aic
      )
      log_info(
        "  LRT: chi-sq=%.2f, df=%d, p=%.4f %s",
        lrt_chisq, lrt_df, lrt_p, sig
      )
      log_info(
        "  Quadratic coef: %.4f (SE=%.4f)",
        quad_coef, quad_se
      )
      log_info(
        "  Conclusion: %s",
        linearity_results.lst[[domain]]$Interpretation
      )
    }
  }

  if (length(linearity_results.lst) > 0) {
    linearity.dt <- rbindlist(linearity_results.lst)
    any_nonlinear <- any(
      linearity.dt$Quadratic_significant
    )

    overall_interpretation <- if (any_nonlinear) {
      nonlinear_domains.v <- linearity.dt[
        Quadratic_significant == TRUE, Domain
      ]
      paste0(
        sprintf(
          "CAUTION: Non-linear EDT for %s. ",
          paste(nonlinear_domains.v, collapse = ", ")
        ),
        "Consider quadratic term sensitivity."
      )
    } else {
      paste0(
        "Linear EDT assumption supported across ",
        "all cognitive domains."
      )
    }

    log_info("")
    log_info("OVERALL: %s", overall_interpretation)

    linearity_output.lst <- list(
      results_table = linearity.dt,
      any_nonlinear = any_nonlinear,
      edt_mean = edt_mean,
      edt_sd = edt_sd,
      interpretation = overall_interpretation
    )

    write_rds_safe(
      linearity_output.lst, linearity.path,
      description = "EDT linearity test results"
    )
    log_info(
      "Saved linearity test to %s", linearity.path
    )
  } else {
    log_warn("No linearity results to save")
    write_rds_safe(
      list(
        error = "No models converged",
        edt_mean = edt_mean, edt_sd = edt_sd
      ),
      linearity.path,
      description = "EDT linearity test (failed)"
    )
  }
} else {
  log_warn(
    "EDT not available in cohort - skipping linearity"
  )
  write_rds_safe(
    list(error = "EDT not available in data"),
    linearity.path,
    description = "EDT linearity test (EDT unavailable)"
  )
}

# ============================================================
# PART 4: EDT EARLY/MIDDLE RESTRICTION
# ============================================================
# Exclude late EDT tertile; re-run LME for z-score
# brain measures. Tests whether findings hold in
# early/middle disease stages only.
# ============================================================
log_section("PART 4: EDT Early/Middle Sensitivity")

if (HAS_EDT) {
  baseline_edt.dt <- cohort.dt[
    order(PTID, YRS_from_bl)
  ][, .SD[1], by = PTID]
  edt_tertiles.v <- quantile(
    baseline_edt.dt$EDT,
    probs = c(1 / 3, 2 / 3), na.rm = TRUE
  )

  cohort.dt[, EDT_tertile := fcase(
    EDT <= edt_tertiles.v[1], "Early",
    EDT <= edt_tertiles.v[2], "Middle",
    default = "Late"
  )]

  log_info(
    paste0(
      "EDT tertile boundaries: ",
      "Early <= %.2f, Middle <= %.2f, Late > %.2f"
    ),
    edt_tertiles.v[1], edt_tertiles.v[2],
    edt_tertiles.v[2]
  )
  log_info(
    "  Early: %d obs, Middle: %d obs, Late: %d obs",
    cohort.dt[EDT_tertile == "Early", .N],
    cohort.dt[EDT_tertile == "Middle", .N],
    cohort.dt[EDT_tertile == "Late", .N]
  )

  cohort_em.dt <- cohort.dt[
    EDT_tertile %in% c("Early", "Middle")
  ]

  log_info(
    "Early+Middle subsample: %d subjects, %d obs",
    cohort_em.dt[, uniqueN(PTID)],
    nrow(cohort_em.dt)
  )

  edt_sensitivity_results.lst <- list()

  for (measure_name in names(BRAIN_Z.lst)) {
    measure <- BRAIN_Z.lst[[measure_name]]
    brain_var <- measure$var
    cov_extra <- measure$cov_extra

    if (!brain_var %in% names(cohort.dt) ||
        sum(!is.na(cohort.dt[[brain_var]])) < 100) {
      log_warn(
        "%s: Not available, skipping", measure_name
      )
      next
    }

    log_info("")
    log_info(
      "=== EDT Sensitivity: %s ===", measure$label
    )

    for (domain in COGNITIVE_DOMAINS) {
      log_info("  --- %s ---", domain)

      domain_data.dt <- cohort_em.dt[
        !is.na(get(domain)) & !is.na(FRS) &
          !is.na(get(brain_var))
      ]
      n_obs <- nrow(domain_data.dt)
      n_subj <- domain_data.dt[, uniqueN(PTID)]

      if (n_obs < 50) {
        log_warn(
          "    Insufficient data (N=%d)", n_obs
        )
        next
      }

      full_data.dt <- cohort.dt[
        !is.na(get(domain)) & !is.na(FRS) &
          !is.na(get(brain_var))
      ]

      f_str <- sprintf(
        "%s ~ %s * FRS + %s %s + (1 | PTID)",
        domain, brain_var, COVARIATES_CROSS,
        cov_extra
      )
      f_m1 <- as.formula(f_str)

      m1_em.fit <- tryCatch({
        lmer(
          f_m1, data = domain_data.dt,
          REML = TRUE,
          control = lmerControl(
            optimizer = "bobyqa"
          )
        )
      }, error = function(e) NULL)

      m1_full.fit <- tryCatch({
        lmer(
          f_m1, data = full_data.dt,
          REML = TRUE,
          control = lmerControl(
            optimizer = "bobyqa"
          )
        )
      }, error = function(e) NULL)

      if (!is.null(m1_em.fit) &&
          !is.null(m1_full.fit)) {
        coef_em.mat <- summary(
          m1_em.fit
        )$coefficients
        coef_full.mat <- summary(
          m1_full.fit
        )$coefficients

        int_pattern <- sprintf(
          "%s:FRS|FRS:%s", brain_var, brain_var
        )
        int_term <- grep(
          int_pattern, rownames(coef_em.mat),
          value = TRUE
        )[1]

        if (!is.na(int_term)) {
          beta_em <- coef_em.mat[
            int_term, "Estimate"
          ]
          se_em <- coef_em.mat[
            int_term, "Std. Error"
          ]
          p_em <- coef_em.mat[
            int_term, "Pr(>|t|)"
          ]

          beta_full <- coef_full.mat[
            int_term, "Estimate"
          ]
          se_full <- coef_full.mat[
            int_term, "Std. Error"
          ]
          p_full <- coef_full.mat[
            int_term, "Pr(>|t|)"
          ]

          result_key <- paste(
            measure_name, domain, sep = "_"
          )
          edt_sensitivity_results.lst[[
            result_key
          ]] <- data.table(
            Measure = measure_name,
            Domain = domain,
            N_obs_em = n_obs,
            N_subj_em = n_subj,
            N_obs_full = nrow(full_data.dt),
            Beta_early_middle = beta_em,
            SE_early_middle = se_em,
            P_early_middle = p_em,
            Beta_full = beta_full,
            SE_full = se_full,
            P_full = p_full,
            Beta_change_pct =
              (beta_em - beta_full) /
              abs(beta_full) * 100,
            Consistent =
              (p_em < 0.05) == (p_full < 0.05)
          )

          sig_em <- ifelse(p_em < 0.05, "*", "")
          sig_full <- ifelse(
            p_full < 0.05, "*", ""
          )
          log_info(
            "    Full: B=%.4f, p=%.4f%s",
            beta_full, p_full, sig_full
          )
          log_info(
            "    E+M:  B=%.4f, p=%.4f%s",
            beta_em, p_em, sig_em
          )
        }
      }
    }
  }

  if (length(edt_sensitivity_results.lst) > 0) {
    edt_sens.dt <- rbindlist(
      edt_sensitivity_results.lst
    )
    all_edt_consistent <- all(edt_sens.dt$Consistent)

    edt_interpretation <- if (all_edt_consistent) {
      paste0(
        "FRS x Brain findings CONSISTENT ",
        "when excluding late EDT."
      )
    } else {
      inconsistent.dt <- edt_sens.dt[
        Consistent == FALSE
      ]
      incon_str <- paste(
        sprintf(
          "%s/%s",
          inconsistent.dt$Measure,
          inconsistent.dt$Domain
        ),
        collapse = ", "
      )
      sprintf(
        paste0(
          "WARNING: Inconsistent for %s ",
          "when excluding late EDT."
        ),
        incon_str
      )
    }

    log_info("")
    log_info("EDT SENSITIVITY: %s", edt_interpretation)

    edt_sens_output.lst <- list(
      results_table = edt_sens.dt,
      tertile_cutoffs = edt_tertiles.v,
      all_consistent = all_edt_consistent,
      interpretation = edt_interpretation
    )

    write_rds_safe(
      edt_sens_output.lst, edt_sens.path,
      description = "EDT Early+Middle sensitivity"
    )
    log_info(
      "Saved EDT sensitivity to %s", edt_sens.path
    )
  }
} else {
  log_warn(
    "EDT not available - skipping Early+Middle test"
  )
}

# ============================================================
# Final Summary
# ============================================================
log_section("Summary")
log_info("SENSITIVITY ANALYSES COMPLETE")
log_info("")

log_info("1. RAW BRAIN MEASURE MODELS:")
if (nrow(raw_summary.dt) > 0) {
  log_info(
    "   - %d models across %s x %s",
    nrow(raw_summary.dt),
    paste(names(BRAIN_MEASURES.lst), collapse = "/"),
    paste(names(CVR_MEASURES.lst), collapse = "/")
  )
  log_info(
    "   - Significant: %d / %d",
    sum(raw_summary.dt$Significant, na.rm = TRUE),
    nrow(raw_summary.dt)
  )
} else {
  log_info("   - No models fitted")
}
log_info("")

log_info("2. ATTRITION ANALYSIS:")
log_info(
  "   - Completers: %d, Dropouts: %d",
  attrition_results.lst$n_completers,
  attrition_results.lst$n_dropouts
)
if (!is.null(attrition_results.lst$baseline_interpretation)) {
  log_info(
    "   - Baseline: %s",
    attrition_results.lst$baseline_interpretation
  )
}
if (!is.null(attrition_results.lst$model_interpretation)) {
  log_info(
    "   - Models: %s",
    attrition_results.lst$model_interpretation
  )
}
log_info("")

if (HAS_EDT && exists("linearity_output.lst") &&
    !is.null(linearity_output.lst$results_table)) {
  log_info("3. EDT LINEARITY TEST:")
  for (i in seq_len(
    nrow(linearity_output.lst$results_table)
  )) {
    row <- linearity_output.lst$results_table[i]
    log_info(
      "   - %s: %s (p=%.4f)",
      row$Domain, row$Interpretation, row$LRT_p
    )
  }
} else {
  log_info(
    "3. EDT LINEARITY TEST: Not performed (no EDT)"
  )
}
log_info("")

if (HAS_EDT && exists("edt_sens_output.lst") &&
    !is.null(edt_sens_output.lst$results_table)) {
  log_info("4. EDT EARLY/MIDDLE SENSITIVITY:")
  log_info(
    "   - %s", edt_sens_output.lst$interpretation
  )
} else {
  log_info("4. EDT EARLY/MIDDLE SENSITIVITY: Not performed")
}

log_script_end("07_lme_sensitivity.R", success = TRUE)
