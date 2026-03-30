#!/usr/bin/env Rscript

# =============================================================================
# 10b_lme_amyloid_negative.R
# =============================================================================
# Sensitivity analysis: LME 3-way interaction in amyloid-negative,
# cognitively unimpaired (CU) participants.
#
# Addresses reviewer concern that AD pathology may dominate
# HVR-cognitive associations, masking true CVR effects.
# Tests whether FRS and CVR_mimic moderate the HVR-cognition
# relationship outside the AD context.
#
# SAMPLE: Amyloid-negative, baseline CU, no sex stratification.
#
# MODEL:
#   Cognition ~ (YRS + YRS^2) x CVR x HVR_z +
#               Age_bl + EDUC + SEX + APOE4 +
#               (YRS | PTID)
#
# INPUTS:
#   - data/derivatives/adni_z-scores.rds
#   - data/derivatives/adni_hc-hvr_adj.rds
#   - data/derivatives/lme/cvr_mimic_scores.rds
#   - data/raw/adsp_phc/ADSP_PHC_PET_Amyloid_Simple*
#   - data/raw/adsp_phc/ADSP_PHC_COGN*
#   - data/raw/adsp_phc/ADSP_PHC_CVRF*
#
# OUTPUTS:
#   - models/results/lme/lme_amyloid_negative.rds
# =============================================================================

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

log_script_start("10b_lme_amyloid_negative.R")
config <- load_config()
validate_config(config)
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "lme_amy_neg",
  default = FALSE
)
output.path <- file.path(
  get_data_path("models", "lme_results_dir"),
  "lme_amyloid_negative.rds"
)
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE.")
  log_script_end(
    "10b_lme_amyloid_negative.R", success = TRUE
  )
  quit(status = 0)
}

# --- Configuration ---
BRAIN_VAR <- "HVR_z"
COGNITIVE_DOMAINS <- get_parameter(
  "cognitive_domains"
)
CVR_MEASURES <- list(
  FRS = "FRS", CVR_mimic = "CVR_mimic"
)

LME_ESTIMATOR <- get_script_setting(
  "lme", "estimator", default = "REML"
)
LME_OPTIMIZER <- get_script_setting(
  "lme", "optimizer", default = "bobyqa"
)
LME_MAXITER <- get_script_setting(
  "lme", "max_iter", default = 100000
)

# ============================================================
# 1. Identify Amyloid-Negative CU Subjects
# ============================================================
log_section("Identifying A-/CU Sample")

amy.dt <- fread(
  get_data_path("raw", "adsp_phc_amyloid"),
  select = c("PTID", "PHC_AMYLOID_STATUS")
)
# Conservative: if any scan positive, classify as A+
amy_subj.dt <- amy.dt[
  order(PTID, -PHC_AMYLOID_STATUS)
][, .SD[1], by = PTID]
amy_neg_ids.v <- amy_subj.dt[
  PHC_AMYLOID_STATUS == 0, unique(PTID)
]
log_info("Amyloid-negative subjects: %d",
         length(amy_neg_ids.v))

# Baseline CU (PHC_Diagnosis == 1)
cog_raw.dt <- fread(
  get_data_path("raw", "adsp_phc_cognition"),
  select = c(
    "PTID", "EXAMDATE", "PHC_Diagnosis",
    "PHC_Age_Cognition", "PHC_Education",
    "PHC_MEM", "PHC_LAN", "PHC_EXF",
    "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE"
  )
)
bl_dx.dt <- cog_raw.dt[
  order(PTID, EXAMDATE), .SD[1], by = PTID
]
cu_ids.v <- bl_dx.dt[
  PHC_Diagnosis == 1, unique(PTID)
]
amy_neg_cu_ids.v <- intersect(
  amy_neg_ids.v, cu_ids.v
)
log_info("Amyloid-negative CU: %d subjects",
         length(amy_neg_cu_ids.v))

# ============================================================
# 2. Build Longitudinal Cohort
# ============================================================
log_section("Building A-/CU Cohort")

# --- HVR z-scores ---
z.dt <- read_rds_safe(
  get_data_path("derivatives", "z_scores"),
  "HVR z-scores"
)
hvr.dt <- z.dt[
  ROI == "HVR" & SIDE == "LR" & ADJ == "NON" &
    PTID %in% amy_neg_cu_ids.v,
  .(PTID, EXAMDATE, MRI_date = EXAMDATE,
    HVR_z = Z)
]
setkey(hvr.dt, PTID, EXAMDATE)
log_info("HVR z-scores: %d obs, %d subjects",
         nrow(hvr.dt), hvr.dt[, uniqueN(PTID)])

# --- Cognitive data (all visits, A-/CU only) ---
cog.dt <- cog_raw.dt[PTID %in% amy_neg_cu_ids.v]
setnames(cog.dt,
  c("PHC_Age_Cognition", "PHC_Education",
    "PHC_MEM", "PHC_LAN", "PHC_EXF",
    "PHC_MEM_SE", "PHC_LAN_SE", "PHC_EXF_SE",
    "PHC_Diagnosis"),
  c("AGE", "EDUC", "MEM", "LAN", "EXF",
    "MEM_SE", "LAN_SE", "EXF_SE", "DX"))
cog.dt[, COG_date := as.Date(EXAMDATE)]
cog.dt[, DX := factor(
  DX, levels = 1:3,
  labels = c("CU", "MCI", "AD")
)]

# --- APOE ---
apoe_file <- get_data_path("raw", "apoe")
apoe.dt <- fread(apoe_file, select = c(
  "PTID", "GENOTYPE"
))
apoe.dt <- unique(apoe.dt[, .(PTID, GENOTYPE)])
apoe.dt[, APOE4 := as.integer(grepl("4", GENOTYPE))]
apoe.dt[, GENOTYPE := NULL]
apoe.dt <- apoe.dt[, .SD[1], by = PTID]
cog.dt <- apoe.dt[cog.dt, on = "PTID"]

# --- SEX ---
cvrf_file <- get_data_path("raw", "adsp_phc_cvrf")
sex_src.dt <- fread(cvrf_file, select = c(
  "PTID", "PHC_Sex"
))
sex_src.dt <- sex_src.dt[, .SD[1], by = PTID]
sex_src.dt[, SEX := fifelse(
  PHC_Sex == 1, "Male", "Female"
)]
sex_src.dt[, PHC_Sex := NULL]
cog.dt <- sex_src.dt[cog.dt, on = "PTID"]

setkey(cog.dt, PTID, EXAMDATE)
log_info("Cognitive data: %d obs, %d subjects",
         nrow(cog.dt), cog.dt[, uniqueN(PTID)])

# --- Rolling join: brain + cognition ---
cohort.dt <- hvr.dt[cog.dt, roll = Inf,
                     nomatch = NULL]

# Visit structure
cohort.dt[order(COG_date), `:=`(
  VISIT = seq_len(.N),
  YRS_from_bl = as.numeric(
    COG_date - COG_date[1]
  ) / 365.25
), by = PTID]

# Require >=2 visits
visit_counts.dt <- cohort.dt[, .N, by = PTID]
keep_ids.v <- visit_counts.dt[N >= 2, PTID]
cohort.dt <- cohort.dt[PTID %in% keep_ids.v]
log_info(
  "After >=2 visits filter: %d subjects, %d obs",
  cohort.dt[, uniqueN(PTID)], nrow(cohort.dt)
)

# --- FRS (raw from CVRF, then z-score) ---
frs.dt <- fread(cvrf_file, select = c(
  "PTID",
  "PHC_ASCVD_10y_FRS_Simple_Ageover30"
))
frs.dt <- frs.dt[, .SD[1], by = PTID]
setnames(frs.dt,
  "PHC_ASCVD_10y_FRS_Simple_Ageover30",
  "FRS_raw"
)
cohort.dt <- merge(
  cohort.dt, frs.dt, by = "PTID", all.x = TRUE
)

# --- CVR_mimic (from script 07 output) ---
cvr.dt <- read_rds_safe(
  get_data_path("models", "cvr_mimic_scores"),
  "CVR_mimic scores"
)
# Map RID -> PTID via ADNI convention
# RID is numeric suffix of PTID (XXX_S_RRRR)
cvr.dt[, RID_char := sprintf("%04d", RID)]
# Build PTID mapping from cohort
ptid_map.dt <- unique(cohort.dt[, .(PTID)])
ptid_map.dt[, RID_char := sub(
  "^\\d+_S_", "", PTID
)]
cvr_merge.dt <- merge(
  ptid_map.dt, cvr.dt[, .(RID_char, CVR_mimic)],
  by = "RID_char", all.x = TRUE
)
cvr_merge.dt[, RID_char := NULL]
cohort.dt <- merge(
  cohort.dt, cvr_merge.dt,
  by = "PTID", all.x = TRUE
)

# Filter to common sample (both FRS and CVR_mimic)
n_pre <- cohort.dt[, uniqueN(PTID)]
cohort.dt <- cohort.dt[
  !is.na(FRS_raw) & !is.na(CVR_mimic)
]
n_post <- cohort.dt[, uniqueN(PTID)]
log_info(
  "Common sample (FRS + CVR_mimic): %d/%d",
  n_post, n_pre
)

# Z-score both on this sample
cohort.dt[, FRS := scale(FRS_raw)[, 1]]
cohort.dt[, CVR_mimic := scale(CVR_mimic)[, 1]]

# Baseline age
setorder(cohort.dt, PTID, EXAMDATE)
bl_age.dt <- cohort.dt[
  , .(Age_bl = AGE[1]), by = PTID
]
cohort.dt <- merge(
  cohort.dt, bl_age.dt,
  by = "PTID", all.x = TRUE
)

# --- Sample summary ---
log_section("A-/CU Cohort Summary")
n_subj <- cohort.dt[, uniqueN(PTID)]
n_obs <- nrow(cohort.dt)
bl.dt <- cohort.dt[VISIT == 1]
log_info("N = %d subjects, %d observations",
         n_subj, n_obs)
log_info(
  "Sex: %d Males, %d Females",
  sum(bl.dt$SEX == "Male"),
  sum(bl.dt$SEX == "Female")
)
log_info(
  "Age: M=%.1f (SD=%.1f), range %.0f-%.0f",
  mean(bl.dt$Age_bl), sd(bl.dt$Age_bl),
  min(bl.dt$Age_bl), max(bl.dt$Age_bl)
)
log_info(
  "Visits: median=%.0f, range %d-%d",
  median(cohort.dt[, .N, by = PTID]$N),
  min(cohort.dt[, .N, by = PTID]$N),
  max(cohort.dt[, .N, by = PTID]$N)
)
log_info(
  "Follow-up: M=%.1f yrs, max=%.1f yrs",
  mean(cohort.dt[
    , max(YRS_from_bl), by = PTID
  ]$V1),
  max(cohort.dt$YRS_from_bl)
)

# CVR-Age correlation check
r_cvr_age <- cor(
  bl.dt$CVR_mimic, bl.dt$Age_bl,
  use = "complete.obs"
)
log_info(
  "r(CVR_mimic, Age) = %.3f (expect ~0)",
  r_cvr_age
)

HAS_APOE <- "APOE4" %in% names(cohort.dt) &&
  !all(is.na(cohort.dt$APOE4))

# ============================================================
# 3. Fit LME Models
# ============================================================
log_section("LME 3-Way Interaction (A-/CU)")

COVARIATES <- if (HAS_APOE) {
  "Age_bl + EDUC + SEX + APOE4"
} else {
  "Age_bl + EDUC + SEX"
}

fit_lme.fn <- function(formula, data, desc) {
  log_info("  Fitting: %s", desc)
  tryCatch({
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
    coefs.mat <- summary(model.fit)$coefficients
    list(
      fit = model.fit,
      coefficients = coefs.mat,
      converged = !isSingular(model.fit),
      is_singular = isSingular(model.fit),
      warnings = warnings.lst,
      n_obs = nrow(data),
      n_subjects = data[, uniqueN(PTID)]
    )
  }, error = function(e) {
    log_warn("    Failed: %s", e$message)
    list(fit = NULL, converged = FALSE,
         error = e$message)
  })
}

extract_3way.fn <- function(model.res, cvr_var) {
  if (is.null(model.res$fit) ||
      !model.res$converged) {
    return(list(
      term = NA, beta = NA, se = NA,
      t_value = NA, p_value = NA
    ))
  }
  coefs.mat <- model.res$coefficients
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

# --- Run models ---
all_results.lst <- list()
summary_rows.lst <- list()

for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]
  log_info("=== CVR: %s ===", CVR_VAR)

  all_results.lst[[cvr_name]] <- list()

  for (domain in COGNITIVE_DOMAINS) {
    formula_str <- sprintf(
      paste0(
        "%s ~ YRS_from_bl * %s * %s + ",
        "I(YRS_from_bl^2) * %s * %s + ",
        "%s + (YRS_from_bl | PTID)"
      ),
      domain, CVR_VAR, BRAIN_VAR,
      CVR_VAR, BRAIN_VAR, COVARIATES
    )

    model.res <- fit_lme.fn(
      as.formula(formula_str), cohort.dt,
      sprintf("%s: YRS x %s x %s",
              domain, CVR_VAR, BRAIN_VAR)
    )

    int.lst <- extract_3way.fn(
      model.res, CVR_VAR
    )

    all_results.lst[[cvr_name]][[domain]] <- list(
      model = model.res,
      interaction = int.lst,
      formula = formula_str,
      n_subjects = model.res$n_subjects,
      n_obs = model.res$n_obs
    )

    if (!is.na(int.lst$p_value)) {
      sig <- ifelse(
        int.lst$p_value < 0.05, " *", ""
      )
      log_info(
        "  %s: B=%.4f (SE=%.4f), p=%.4f%s",
        domain, int.lst$beta, int.lst$se,
        int.lst$p_value, sig
      )
    } else {
      log_info(
        "  %s: did not converge", domain
      )
    }

    summary_rows.lst[[
      length(summary_rows.lst) + 1
    ]] <- data.table(
      CVR_Measure = CVR_VAR,
      Domain = domain,
      N_subjects = model.res$n_subjects,
      N_obs = model.res$n_obs,
      Beta = int.lst$beta,
      SE = int.lst$se,
      t_value = int.lst$t_value,
      p_value = int.lst$p_value
    )
  }
}

summary.dt <- rbindlist(summary_rows.lst)

# FDR within each CVR measure (3 tests each)
for (cvr_name in names(CVR_MEASURES)) {
  cvr_col <- CVR_MEASURES[[cvr_name]]
  mask <- summary.dt$CVR_Measure == cvr_col
  summary.dt[mask,
    p_fdr := p.adjust(p_value, method = "BH")
  ]
}
summary.dt[, Significant := p_value < 0.05]
summary.dt[, Significant_FDR := p_fdr < 0.05]

# ============================================================
# 4. Summary
# ============================================================
log_section("KEY FINDINGS (A-/CU)")
log_info("")
for (cvr_name in names(CVR_MEASURES)) {
  CVR_VAR <- CVR_MEASURES[[cvr_name]]
  sub.dt <- summary.dt[CVR_Measure == CVR_VAR]
  log_info(
    "%s: %d/3 sig (uncorr), %d/3 sig (FDR)",
    CVR_VAR,
    sum(sub.dt$Significant, na.rm = TRUE),
    sum(sub.dt$Significant_FDR, na.rm = TRUE)
  )
  for (i in seq_len(nrow(sub.dt))) {
    r <- sub.dt[i]
    sig <- ifelse(r$Significant, " *", "")
    log_info(
      "  %s: B=%.4f, p=%.4f, p_fdr=%.4f%s",
      r$Domain, r$Beta, r$p_value,
      r$p_fdr, sig
    )
  }
}

# ============================================================
# 5. Save
# ============================================================
log_section("Saving Results")

results.lst <- list(
  summary_table = summary.dt,
  results = all_results.lst,
  sample_info = list(
    n_subjects = n_subj,
    n_observations = n_obs,
    n_male = sum(bl.dt$SEX == "Male"),
    n_female = sum(bl.dt$SEX == "Female"),
    age_mean = mean(bl.dt$Age_bl),
    age_sd = sd(bl.dt$Age_bl),
    sample = "Amyloid-negative, baseline CU"
  ),
  cohort = cohort.dt
)

ensure_directory(dirname(output.path))
write_rds_safe(
  results.lst, output.path,
  "A-/CU LME results"
)

log_info("")
log_script_end(
  "10b_lme_amyloid_negative.R", success = TRUE
)
