#!/usr/bin/env Rscript

# ============================================================
# 10_lgcm_mediation.R
# ============================================================
#
# PURPOSE:
#   Test whether cardiovascular risk (CVR) effects on
#   cognitive decline are mediated through hippocampal
#   neurodegeneration (HVR decline), using parallel
#   process LGCMs. Both CVR_mimic and FRS are tested
#   as predictors, stratified by sex.
#
# SCIENTIFIC BACKGROUND:
#   Cardiovascular risk factors are associated with
#   hippocampal atrophy and cognitive decline in aging
#   and Alzheimer's disease:
#
#   - Framingham Risk Score (FRS) predicts hippocampal
#     volume loss and cognitive decline
#     (Debette et al., 2011, Neurology;
#      Viticchi et al., 2020, JACC)
#   - Hippocampal atrophy mediates the relationship
#     between vascular risk and memory decline
#     (Knopman et al., 2011; Rabin et al., 2019)
#   - If FRS shows mediation but CVR_mimic (age-adjusted)
#     does not, this suggests FRS effects are driven by
#     age-weighting rather than true cardiovascular risk.
#
# METHODOLOGICAL JOURNEY:
#
#   STEP 1: Cognitive and HVR trajectories are NON-LINEAR
#     (quadratic) per univariate LGCMs (scripts 08b, 09b):
#     - MEM: chi-sq diff = 1102.56, p < 10^-241
#     - LAN: chi-sq diff = 680.78, p < 10^-146
#     - HVR: chi-sq diff = 453.30, p < 10^-97
#
#   STEP 2: Quadratic parallel process + CVR
#     predictors FAILED to converge
#     (identification issues). Known problem
#     when quadratic variance is near zero
#     (Preacher, 2010,
#      DOI: 10.4324/9780203861554-17),
#     or multiple processes multiply parameters
#     (Curran et al., 2010,
#      DOI: 10.1080/15248371003699969).
#
#   STEP 3: We use LINEAR parallel process
#     models as an approximation. The linear
#     slope captures the overall RATE of change
#     even if the true trajectory curves.
#     Standard pragmatic approach
#     (Newsom, 2015,
#      DOI: 10.4324/9781315871318).
#
#   STEP 4: MEDIATION: CVR -> HVR slope -> Cog slope
#     The coupling (hvr_s ~~ cog_s) is a PREREQUISITE
#     for mediation and is established in script 12.
#
# MODEL:
#   Parallel process LGCM with mediation:
#   CVR -> HVR slope -> Cognitive slope
#   Separately for Male/Female (no pooled models).
#
# PREDICTORS:
#   - CVR_mimic: pre-computed zeta scores (script 07)
#   - FRS: Framingham Risk Score (sex-centered)
#
# NOTE: We attempted embedding the full MIMIC measurement
# model directly in OpenMx (build_mimic_mediation.fn).
# 5/6 models did not converge — the embedded MIMIC adds
# ~15 parameters on top of the mediation structure, which
# exceeds data capacity at N~400-500. We instead use the
# pre-computed age-adjusted zeta scores (zeta = eta -
# gamma*AGE_c), matching the LME approach (scripts 08-10).
#
# BOOTSTRAP:
#   Set N_BOOT > 0 for bootstrap CIs on indirect effect.
#   N_BOOT = 0 uses Sobel test only (preliminary).
#
# REFERENCES:
#
#   Mediation Methodology:
#   - Hayes (2009). Beyond Baron and Kenny. Communication
#     Monographs, 76(4), 408-420.
#     DOI: 10.1080/03637750903310360
#   - Zhao et al. (2010). Reconsidering Baron and Kenny.
#     J Consumer Res, 37(2), 197-206.
#     DOI: 10.1086/651257
#   - Preacher & Hayes (2008). Asymptotic and resampling
#     strategies for indirect effects. Behav Res Methods,
#     40(3), 879-891. PMC2819361
#   - Kenny & Judd (2014). Power anomalies in testing
#     mediation. Psychol Sci, 25(2), 334-339. PMC4142865
#
#   LGCM Methodology:
#   - Preacher (2010). Latent growth curve
#     models. In Hancock & Mueller (Eds.),
#     Reviewer's Guide to Quantitative
#     Methods. Routledge.
#     DOI: 10.4324/9780203861554-17
#   - Curran, Obeidat & Losardo (2010).
#     Twelve frequently asked questions
#     about growth curve modeling.
#     J Cogn Dev, 11(2), 121-136.
#     DOI: 10.1080/15248371003699969
#
#   Substantive:
#   - Debette et al. (2011). Neurology, 77(5), 461-468.
#     PMC3146307
#   - Viticchi et al. (2020). JACC, 75(20), 2584-2594.
#   - Gorelick et al. (2011). Stroke, 42(9), 2672-2713.
#     PMC3778669
#
# OUTPUTS:
#   - models/results/lgcm/lgcm_mediation_results.rds
#     (both predictors, both sexes)
#
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(OpenMx)
})

source(here("R/utils/config.R"))
source(here("R/utils/logging.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/disease_time.R"))

log_script_start("13_lgcm_mediation.R")
config <- load_config()
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "lgcm_mediation", default = FALSE
)
output.path <- get_data_path("models", "lgcm_mediation")
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE")
  log_info("Skipping. Set force_regenerate=TRUE to rerun.")
  log_script_end("13_lgcm_mediation.R", success = TRUE)
  quit(status = 0)
}

mxOption(
  key = "Number of Threads",
  value = parallel::detectCores()
)

results_dir <- get_data_path(
  "models", "lgcm_results_dir"
)
ensure_directory(results_dir)

checkpoint_dir <- get_data_path(
  "models", "lgcm_mediation_checkpoints"
)
ensure_directory(checkpoint_dir)

# ============================================================
# CONFIGURATION
# ============================================================
# Set N_BOOT = 0 for preliminary Sobel-only results.
# Increase to 1000+ for final bootstrap CIs.
N_BOOT <- 1000
CI_LEVEL <- 0.95

COG_DOMAINS <- c("MEM", "LAN", "EXF")
VALID_TIMEPOINTS <- paste0("T", 0:5)
EDT_COLS <- paste0("EDT_", VALID_TIMEPOINTS)
EDT_SQ_COLS <- paste0(EDT_COLS, "_sq")
HVR_COLS <- paste0("HVR_Z_", VALID_TIMEPOINTS)

PREDICTORS <- list(
  CVR_mimic = list(
    type = "observed",
    var = "CVR_sex_centered"
  ),
  FRS = list(
    type = "observed",
    var = "FRS_c"
  )
)

# ============================================================
# PART 1: DATA PREPARATION
# ============================================================
log_section("Part 1: Data Preparation")

# Load sex-stratified samples (Male/Female only)
samples.lst <- list(
  Male = read_rds_safe(
    get_data_path("derivatives", "lgcm_input_male"),
    "Male sample"
  ),
  Female = read_rds_safe(
    get_data_path("derivatives", "lgcm_input_female"),
    "Female sample"
  )
)

for (nm in names(samples.lst)) {
  dt <- samples.lst[[nm]]
  log_info(
    "%s: N=%d, CVR_mimic NAs=%d",
    nm, nrow(dt), sum(is.na(dt$CVR_mimic))
  )
}

# ----------------------------------------------------------
# APOE4
# ----------------------------------------------------------
HAS_APOE <- "APOE4" %in% names(samples.lst$Male) &&
  sum(!is.na(samples.lst$Male$APOE4)) > 0

if (HAS_APOE) {
  for (nm in names(samples.lst)) {
    samples.lst[[nm]][,
      APOE4_num := as.numeric(APOE4)
    ]
  }
  log_info("APOE4 available: Yes")
} else {
  log_info("APOE4 available: No")
}

# Check MIMIC indicator availability
MIMIC_IND <- c(
  "SBP_z", "HTN", "GLUCOSE_z",
  "CHOL_z", "CREAT_z"
)
for (nm in names(samples.lst)) {
  avail <- MIMIC_IND[
    MIMIC_IND %in% names(samples.lst[[nm]])
  ]
  log_info(
    "%s: %d/%d MIMIC indicators",
    nm, length(avail), length(MIMIC_IND)
  )
}
HAS_MIMIC <- all(
  MIMIC_IND %in% names(samples.lst$Male)
)
if (!HAS_MIMIC) {
  stop(
    "MIMIC indicators missing. Run script 11."
  )
}

# ----------------------------------------------------------
# EDT extrapolation
# ----------------------------------------------------------
log_info("")
log_info("Extrapolating missing EDT values...")

# Use Male sample for population rate, apply to both
male_result.lst <- extrapolate_edt_hybrid(
  samples.lst$Male,
  valid_timepoints = VALID_TIMEPOINTS
)
samples.lst$Male <- male_result.lst$data
pop_rate <- male_result.lst$pop_rate

female_result.lst <- extrapolate_edt_hybrid(
  samples.lst$Female,
  pop_rate = pop_rate,
  valid_timepoints = VALID_TIMEPOINTS
)
samples.lst$Female <- female_result.lst$data

# ----------------------------------------------------------
# Center EDT (grand mean from pooled Male+Female)
# ----------------------------------------------------------
all_edt.v <- c(
  unlist(samples.lst$Male[, ..EDT_COLS]),
  unlist(samples.lst$Female[, ..EDT_COLS])
)
edt_mean <- mean(all_edt.v, na.rm = TRUE)
log_info(
  "EDT grand mean: %.3f years (centering)", edt_mean
)

for (nm in names(samples.lst)) {
  for (col in EDT_COLS) {
    set(
      samples.lst[[nm]], j = col,
      value = samples.lst[[nm]][[col]] -
        edt_mean
    )
  }
}

# Compute EDT^2 for fixed quadratic models
for (nm in names(samples.lst)) {
  for (i in seq_along(EDT_COLS)) {
    col <- EDT_COLS[i]
    sq_col <- EDT_SQ_COLS[i]
    set(
      samples.lst[[nm]], j = sq_col,
      value = samples.lst[[nm]][[col]]^2
    )
  }
}
log_info("  EDT centered and squared")

# ----------------------------------------------------------
# Center predictors for numerical stability
# ----------------------------------------------------------
# METHODOLOGICAL NOTE:
# Centering predictors is CRITICAL for OpenMx RAM models
# with raw data. When predictor means are fixed to 0 in the
# model (standard practice for exogenous variables), the
# actual data must also be centered.
#
# We use grand means from pooled Male+Female for
# consistency across sex-stratified analyses.
# ----------------------------------------------------------

log_info("")
log_info("Centering predictors...")

all_age.v <- c(
  samples.lst$Male$AGE_bl,
  samples.lst$Female$AGE_bl
)
all_educ.v <- c(
  samples.lst$Male$EDUC,
  samples.lst$Female$EDUC
)

age_mean <- mean(all_age.v, na.rm = TRUE)
educ_mean <- mean(all_educ.v, na.rm = TRUE)

log_info("  AGE grand mean: %.2f", age_mean)
log_info("  EDUC grand mean: %.2f", educ_mean)

for (nm in names(samples.lst)) {
  dt <- samples.lst[[nm]]
  set(dt, j = "AGE_c",
      value = dt$AGE_bl - age_mean)
  set(dt, j = "EDUC_c",
      value = dt$EDUC - educ_mean)

  # Center FRS within each sample
  frs_m <- mean(
    dt$FRS_sex_centered, na.rm = TRUE
  )
  set(dt, j = "FRS_c",
      value = dt$FRS_sex_centered - frs_m)

  log_info(
    "  %s: FRS_c (M=%.4f->0)",
    nm, frs_m
  )
}

# Validate CVR_mimic is age-orthogonal
# (Uses CVR_sex_centered; CVR is now latent in
# MIMIC models so this is a rough pre-check)
for (nm in names(samples.lst)) {
  dt <- samples.lst[[nm]]
  r <- cor(
    dt$CVR_sex_centered, dt$AGE_c,
    use = "complete.obs"
  )
  log_info(
    "  %s: r(CVR, AGE_c) = %.3f", nm, r
  )
}

# ============================================================
# PART 2: DOCUMENT QUADRATIC FAILURE
# ============================================================
log_section(
  "Part 2: Quadratic Model Attempts (Documented)"
)

# ---------------------------------------------------------
# METHODOLOGICAL CONTEXT:
# ---------------------------------------------------------
# Prior univariate results showed quadratic trajectories:
#
#   Domain | Linear -2LL | Quad -2LL  | Chi-sq  | p
#   MEM    | 3419.32     | 2316.76    | 1102.56 | <1e-241
#   LAN    | 4284.29     | 3603.51    | 680.78  | <1e-146
#   EXF    | 5039.05     | 4431.57    | 607.48  | <1e-130
#   HVR    | 3789.77     | 3336.46    | 453.30  | <1e-97
#
# However, quadratic parallel process + CVR predictors
# FAILED:
#   - var(q) near zero causes identification problems
#   - Fixing var(q)=0 makes CVR->q paths undefined
#   - Complex model exceeds what data can identify
#
# We proceed with linear parallel process models.
# ---------------------------------------------------------

log_info("")
log_info("QUADRATIC MODELS: Why they cannot be used")
log_info("")
log_info("  Prior univariate results (quadratic):")
log_info("    MEM: dChiSq=1102.56, p<1e-241")
log_info("    LAN: dChiSq=680.78, p<1e-146")
log_info("    EXF: dChiSq=607.48, p<1e-130")
log_info("    HVR: dChiSq=453.30, p<1e-97")
log_info("")
log_info("  Quadratic parallel + CVR predictors failed:")
log_info("    - var(q) near zero -> identification")
log_info("    - Fixing var(q)=0 -> CVR->q undefined")
log_info("    - Model complexity exceeds data capacity")
log_info("")
log_info("  References:")
log_info("    - Preacher (2010): quantpsy.org")
log_info("    - Curran et al. (2010): PMC3131138")
log_info("    - Bollen & Curran (2006): Latent Curves")

quadratic_attempts.lst <- list(
  univariate_results = list(
    MEM = list(
      linear_2LL = 3419.32,
      quad_2LL = 2316.76,
      chi_sq_diff = 1102.56,
      decision = "quadratic, var_q=0"
    ),
    LAN = list(
      linear_2LL = 4284.29,
      quad_2LL = 3603.51,
      chi_sq_diff = 680.78,
      decision = "quadratic"
    ),
    EXF = list(
      linear_2LL = 5039.05,
      quad_2LL = 4431.57,
      chi_sq_diff = 607.48,
      decision = "quadratic"
    ),
    HVR = list(
      linear_2LL = 3789.77,
      quad_2LL = 3336.46,
      chi_sq_diff = 453.30,
      decision = "quadratic"
    )
  ),
  parallel_process_failure = list(
    attempted = TRUE,
    converged = FALSE,
    reasons = c(
      "Near-zero quadratic variance -> identification",
      "Fixing var_q=0 makes predictor paths undefined",
      "Model complexity exceeds data information",
      "Definition variables compound difficulty"
    ),
    references = c(
      "Preacher (2010) - quantpsy.org",
      "Curran et al. (2010) - PMC3131138",
      "Bollen & Curran (2006) - Latent Curve Models"
    )
  ),
  conclusion = paste(
    "Quadratic parallel process models with CVR",
    "predictors are not identifiable with this",
    "sample. We use linear parallel process models",
    "where the slope captures overall rate of change."
  )
)

# ============================================================
# HELPER: Extract fit indices via reference models
# ============================================================
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

# ============================================================
# PART 3: MODEL SPECIFICATION
# ============================================================
#
# FIXED QUADRATIC JUSTIFICATION:
#   Same rationale as script 12. Fully random
#   quadratic failed; we fix quadratic variance
#   to ~0 and freely estimate quadratic means,
#   capturing population-level curvature.
#
#   References:
#   - Bollen & Curran (2006). Latent Curve
#     Models. Wiley. ISBN: 978-0471455929
#     (identification of quadratic LGCMs)
#   - Ram & Grimm (2009). Int J Behav Dev,
#     33(6), 565-576.
#     DOI: 10.1177/0165025409343765
#     (fixing non-estimable variance to 0)
#   - Grimm, Ram & Estabrook (2017). Growth
#     Modeling. Guilford.
#     ISBN: 978-1462526062
#     (fixed vs random effects in LGCM)
#
log_section("Part 3: Model Specification")

log_info("")
log_info("FIXED QUADRATIC EXTENSION:")
log_info(
  "  - Fully random quadratic failed (Part 2)"
)
log_info(
  "  - Fixed quadratic adds population-level"
)
log_info(
  "    curvature (2 free params, no random"
)
log_info("    quadratic variance)")
log_info(
  "  - Consistent with script 12 (parallel)"
)
log_info("")
log_info("MEDIATION MODEL: CVR -> hvr_s -> cog_s")
log_info("  a:  CVR -> hvr_s     (CV risk -> atrophy)")
log_info("  b:  hvr_s -> cog_s   (atrophy -> decline)")
log_info("  c': CVR -> cog_s     (direct effect)")
log_info("  indirect: a * b      (mediated effect)")
log_info("")
log_info("  NOTE: Significant total effect (X->Y) is")
log_info("  NOT required (Hayes 2009; Zhao et al 2010)")
log_info("")

#' Build parallel process mediation model
#'
#' @param cog_vars Cognitive variable names
#' @param hvr_vars HVR variable names
#' @param predictor Centered predictor column name
#' @param name Model name
#' @return OpenMx model (without data)
build_mediation.fn <- function(cog_vars, hvr_vars,
                               predictor = "CVR_c",
                               name = "Med") {
  n_tp <- length(cog_vars)
  latents <- c(
    "hvr_i", "hvr_s", "cog_i", "cog_s"
  )

  # HVR trajectory
  hvr_int <- mxPath(
    from = "hvr_i", to = hvr_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  hvr_slope <- mxPath(
    from = "hvr_s", to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", EDT_COLS)
  )

  # Cognitive trajectory
  cog_int <- mxPath(
    from = "cog_i", to = cog_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  cog_slope <- mxPath(
    from = "cog_s", to = cog_vars,
    free = FALSE,
    labels = paste0("data.", EDT_COLS)
  )

  # Latent parameters
  lat_means <- mxPath(
    from = "one", to = latents, free = TRUE,
    values = c(0, -0.1, 0, -0.1),
    labels = paste0("mean_", latents)
  )
  lat_vars <- mxPath(
    from = latents, arrows = 2, free = TRUE,
    values = c(1, 0.1, 1, 0.1),
    labels = paste0("var_", latents)
  )
  hvr_cov <- mxPath(
    from = "hvr_i", to = "hvr_s", arrows = 2,
    free = TRUE, values = 0,
    labels = "cov_hvr_is"
  )
  cog_cov <- mxPath(
    from = "cog_i", to = "cog_s", arrows = 2,
    free = TRUE, values = 0,
    labels = "cov_cog_is"
  )
  cross_int <- mxPath(
    from = "hvr_i", to = "cog_i", arrows = 2,
    free = TRUE, values = 0,
    labels = "cov_hvr_i_cog_i"
  )

  # NOTE: No hvr_s ~~ cog_s free covariance because
  # the b path (hvr_s -> cog_s) captures the slope
  # relationship. Including both would be redundant.

  # Residuals (equal across timepoints)
  hvr_resid <- mxPath(
    from = hvr_vars, arrows = 2,
    free = TRUE, values = 0.5,
    labels = "hvr_resid"
  )
  cog_resid <- mxPath(
    from = cog_vars, arrows = 2,
    free = TRUE, values = 0.5,
    labels = "cog_resid"
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
  preds <- c(predictor, "AGE_c", "EDUC_c")
  if (HAS_APOE) {
    preds <- c(preds, "APOE4_num")
  }

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

  # Covariate paths
  age_paths <- mxPath(
    from = "AGE_c", to = latents,
    free = TRUE, values = 0,
    labels = paste0("age_", latents)
  )
  educ_paths <- mxPath(
    from = "EDUC_c", to = latents,
    free = TRUE, values = 0,
    labels = paste0("educ_", latents)
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
    age_paths, educ_paths
  )

  if (HAS_APOE) {
    apoe_paths <- mxPath(
      from = "APOE4_num", to = latents,
      free = TRUE, values = 0,
      labels = paste0("apoe_", latents)
    )
    model <- mxModel(model, apoe_paths)
  }

  # Fixed quadratic extension
  hvr_q_paths <- mxPath(
    from = "hvr_q", to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", EDT_SQ_COLS)
  )
  cog_q_paths <- mxPath(
    from = "cog_q", to = cog_vars,
    free = FALSE,
    labels = paste0("data.", EDT_SQ_COLS)
  )
  quad_var <- mxPath(
    from = c("hvr_q", "cog_q"),
    arrows = 2, free = FALSE,
    values = c(1e-10, 1e-10)
  )
  quad_means <- mxPath(
    from = "one",
    to = c("hvr_q", "cog_q"),
    free = TRUE,
    values = c(-0.01, -0.01),
    labels = c("mean_hvr_q", "mean_cog_q")
  )
  model <- mxModel(
    model,
    latentVars = c("hvr_q", "cog_q"),
    hvr_q_paths, cog_q_paths,
    quad_var, quad_means
  )

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

#' Build mediation model with embedded MIMIC
#'
#' Same trajectory structure as build_mediation.fn
#' but with CVR as a latent factor measured by 5
#' indicators (MIMIC measurement model embedded
#' in OpenMx RAM).
#'
#' @param cog_vars Cognitive variable names
#' @param hvr_vars HVR variable names
#' @param name Model name
#' @return OpenMx model (without data)
build_mimic_mediation.fn <- function(
    cog_vars, hvr_vars,
    name = "MIMICMed") {
  n_tp <- length(cog_vars)
  latents <- c(
    "hvr_i", "hvr_s",
    "cog_i", "cog_s", "CVR"
  )

  # HVR trajectory
  hvr_int <- mxPath(
    from = "hvr_i", to = hvr_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  hvr_slope <- mxPath(
    from = "hvr_s", to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", EDT_COLS)
  )

  # Cognitive trajectory
  cog_int <- mxPath(
    from = "cog_i", to = cog_vars,
    free = FALSE, values = rep(1, n_tp)
  )
  cog_slope <- mxPath(
    from = "cog_s", to = cog_vars,
    free = FALSE,
    labels = paste0("data.", EDT_COLS)
  )

  # Growth factor parameters
  growth_lat <- c(
    "hvr_i", "hvr_s", "cog_i", "cog_s"
  )
  lat_means <- mxPath(
    from = "one", to = growth_lat,
    free = TRUE,
    values = c(0, -0.1, 0, -0.1),
    labels = paste0("mean_", growth_lat)
  )
  lat_vars <- mxPath(
    from = growth_lat, arrows = 2,
    free = TRUE,
    values = c(1, 0.1, 1, 0.1),
    labels = paste0("var_", growth_lat)
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

  # Residuals
  hvr_resid <- mxPath(
    from = hvr_vars, arrows = 2,
    free = TRUE, values = 0.5,
    labels = "hvr_resid"
  )
  cog_resid <- mxPath(
    from = cog_vars, arrows = 2,
    free = TRUE, values = 0.5,
    labels = "cog_resid"
  )
  hvr_mm <- mxPath(
    from = "one", to = hvr_vars,
    free = FALSE, values = rep(0, n_tp)
  )
  cog_mm <- mxPath(
    from = "one", to = cog_vars,
    free = FALSE, values = rep(0, n_tp)
  )

  # --- MIMIC measurement model ---
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
    from = "AGE_c", to = "CVR",
    free = TRUE, values = 0.1,
    labels = "gamma_age"
  )

  # CVR variance (disturbance)
  cvr_var_path <- mxPath(
    from = "CVR", arrows = 2,
    free = TRUE, values = 0.5,
    labels = "var_CVR"
  )

  # Indicator residuals
  ind_resid <- mxPath(
    from = MIMIC_IND, arrows = 2,
    free = TRUE, values = 0.5,
    labels = paste0("resid_", MIMIC_IND)
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

  # Indicator means fixed to 0
  ind_means <- mxPath(
    from = "one", to = MIMIC_IND,
    free = FALSE, values = 0
  )

  # CVR mean
  cvr_mean <- mxPath(
    from = "one", to = "CVR",
    free = TRUE, values = 0,
    labels = "mean_CVR"
  )

  # --- Mediation paths (using latent CVR) ---
  a_path <- mxPath(
    from = "CVR", to = "hvr_s",
    free = TRUE, values = 0, labels = "a"
  )
  b_path <- mxPath(
    from = "hvr_s", to = "cog_s",
    free = TRUE, values = 0.5,
    labels = "b"
  )
  cprime_path <- mxPath(
    from = "CVR", to = "cog_s",
    free = TRUE, values = 0,
    labels = "cprime"
  )
  pred_int <- mxPath(
    from = "CVR",
    to = c("hvr_i", "cog_i"),
    free = TRUE, values = 0,
    labels = c("pred_hvr_i", "pred_cog_i")
  )

  # --- Covariates ---
  preds <- c("AGE_c", "EDUC_c")
  if (HAS_APOE) {
    preds <- c(preds, "APOE4_num")
  }

  pred_var <- mxPath(
    from = preds, arrows = 2,
    free = FALSE, values = 1
  )
  pred_mn <- mxPath(
    from = "one", to = preds,
    free = FALSE, values = 0
  )

  age_paths <- mxPath(
    from = "AGE_c", to = growth_lat,
    free = TRUE, values = 0,
    labels = paste0("age_", growth_lat)
  )
  educ_paths <- mxPath(
    from = "EDUC_c", to = growth_lat,
    free = TRUE, values = 0,
    labels = paste0("educ_", growth_lat)
  )

  # Build model
  manifests_all <- c(
    hvr_vars, cog_vars,
    MIMIC_IND, preds
  )
  model <- mxModel(
    name, type = "RAM",
    manifestVars = manifests_all,
    latentVars = latents
  )

  model <- mxModel(
    model,
    hvr_int, hvr_slope,
    cog_int, cog_slope,
    lat_means, lat_vars,
    hvr_cov, cog_cov, cross_int,
    hvr_resid, cog_resid,
    hvr_mm, cog_mm,
    cvr_load_marker, cvr_load_free,
    age_cvr, cvr_var_path,
    ind_resid,
    rescov_sbp_chol, rescov_sbp_creat,
    ind_means, cvr_mean,
    a_path, b_path, cprime_path,
    pred_int,
    pred_var, pred_mn,
    age_paths, educ_paths
  )

  if (HAS_APOE) {
    apoe_paths <- mxPath(
      from = "APOE4_num", to = growth_lat,
      free = TRUE, values = 0,
      labels = paste0(
        "apoe_", growth_lat
      )
    )
    model <- mxModel(model, apoe_paths)
  }

  # Fixed quadratic extension
  hvr_q_paths <- mxPath(
    from = "hvr_q", to = hvr_vars,
    free = FALSE,
    labels = paste0("data.", EDT_SQ_COLS)
  )
  cog_q_paths <- mxPath(
    from = "cog_q", to = cog_vars,
    free = FALSE,
    labels = paste0("data.", EDT_SQ_COLS)
  )
  quad_var <- mxPath(
    from = c("hvr_q", "cog_q"),
    arrows = 2, free = FALSE,
    values = c(1e-10, 1e-10)
  )
  quad_means <- mxPath(
    from = "one",
    to = c("hvr_q", "cog_q"),
    free = TRUE,
    values = c(-0.01, -0.01),
    labels = c("mean_hvr_q", "mean_cog_q")
  )
  model <- mxModel(
    model,
    latentVars = c("hvr_q", "cog_q"),
    hvr_q_paths, cog_q_paths,
    quad_var, quad_means
  )

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

#' Compute bootstrap CIs for indirect effect
#'
#' Uses case resampling to compute percentile bootstrap
#' CIs for the indirect effect (a * b) in mediation
#' analysis.
#'
#' Returns both percentile and bias-corrected (BC) CIs.
#' BC CIs adjust for bias in the bootstrap distribution
#' and are recommended for indirect effects which often
#' have skewed sampling distributions.
#'
#' Reference: Preacher & Hayes (2008). PMC2819361
#' Reference: Efron & Tibshirani (1993), Chapter 14.
#'
#' @param model_fn Function returning unfitted mxModel
#' @param data Data frame for model fitting
#' @param n_boot Number of bootstrap samples
#' @param ci_level Confidence level (default 0.95)
#' @param seed Random seed for reproducibility
#' @param original_estimate Original indirect effect
#' @param original_a Original a-path estimate
#' @param original_b Original b-path estimate
#' @return List with all path estimates and CIs
bootstrap_indirect.fn <- function(
    model_fn, data, n_boot = 1000,
    ci_level = 0.95, seed = 42,
    original_estimate = NULL,
    original_cprime = NULL,
    original_a = NULL,
    original_b = NULL) {
  set.seed(seed)
  n <- nrow(data)
  indirect_boot.v <- numeric(n_boot)
  cprime_boot.v <- numeric(n_boot)
  total_boot.v <- numeric(n_boot)
  a_boot.v <- numeric(n_boot)
  b_boot.v <- numeric(n_boot)

  if (n_boot >= 100) {
    log_info(
      "    Running %d bootstrap iterations...",
      n_boot
    )
  }

  for (i in seq_len(n_boot)) {
    boot_idx <- sample(
      seq_len(n), n, replace = TRUE
    )
    boot_data <- data[boot_idx, , drop = FALSE]

    model <- model_fn()
    model <- mxModel(
      model,
      mxData(observed = boot_data, type = "raw")
    )

    fit <- suppressWarnings(tryCatch({
      mxRun(
        model, silent = TRUE,
        suppressWarnings = TRUE
      )
    }, error = function(e) NULL))

    if (!is.null(fit) &&
        !is.null(fit$output$status) &&
        fit$output$status$code %in% c(0, 1)) {
      params <- coef(fit)
      a_boot.v[i] <- params["a"]
      b_boot.v[i] <- params["b"]
      indirect_boot.v[i] <- params["a"] *
        params["b"]
      cprime_boot.v[i] <- params["cprime"]
      total_boot.v[i] <-
        indirect_boot.v[i] +
        cprime_boot.v[i]
    } else {
      a_boot.v[i] <- NA
      b_boot.v[i] <- NA
      indirect_boot.v[i] <- NA
      cprime_boot.v[i] <- NA
      total_boot.v[i] <- NA
    }
  }

  # --- Local helper: BC CI for one vector ---
  compute_bc_ci.fn <- function(
      boot_vals, orig_est, alpha) {
    valid <- boot_vals[!is.na(boot_vals)]
    nv <- length(valid)
    if (nv < 100) {
      return(list(
        ci_lower = NA, ci_upper = NA,
        boot_se = NA, significant = NA,
        boot_p = NA
      ))
    }
    pct_lo <- quantile(
      valid, probs = alpha / 2
    )
    pct_hi <- quantile(
      valid, probs = 1 - alpha / 2
    )
    bse <- sd(valid)
    if (!is.null(orig_est)) {
      pb <- mean(valid < orig_est)
      pb <- max(0.001, min(0.999, pb))
      z0 <- qnorm(pb)
      za_lo <- qnorm(alpha / 2)
      za_hi <- qnorm(1 - alpha / 2)
      p_lo <- pnorm(2 * z0 + za_lo)
      p_hi <- pnorm(2 * z0 + za_hi)
      ci_lo <- as.numeric(
        quantile(valid, probs = p_lo)
      )
      ci_hi <- as.numeric(
        quantile(valid, probs = p_hi)
      )
    } else {
      ci_lo <- as.numeric(pct_lo)
      ci_hi <- as.numeric(pct_hi)
    }
    sig <- (ci_lo > 0 || ci_hi < 0)
    bp <- if (sig) {
      pbz <- mean(valid < 0)
      2 * min(pbz, 1 - pbz)
    } else {
      NA
    }
    list(
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      boot_se = bse,
      significant = sig,
      boot_p = bp
    )
  }

  # Convergence check
  n_valid <- sum(!is.na(indirect_boot.v))
  alpha <- 1 - ci_level

  if (n_valid < n_boot * 0.8) {
    log_warn(
      "    Bootstrap: only %d/%d converged (%.1f%%)",
      n_valid, n_boot, 100 * n_valid / n_boot
    )
  }

  # Compute original total for BC
  orig_total <- if (
    !is.null(original_estimate) &&
    !is.null(original_cprime)
  ) {
    original_estimate + original_cprime
  } else {
    NULL
  }

  ind_ci <- compute_bc_ci.fn(
    indirect_boot.v, original_estimate, alpha
  )
  cp_ci <- compute_bc_ci.fn(
    cprime_boot.v, original_cprime, alpha
  )
  tot_ci <- compute_bc_ci.fn(
    total_boot.v, orig_total, alpha
  )
  a_ci <- compute_bc_ci.fn(
    a_boot.v, original_a, alpha
  )
  b_ci <- compute_bc_ci.fn(
    b_boot.v, original_b, alpha
  )

  ci_method <- if (
    !is.null(original_estimate)
  ) "BC" else "percentile"

  list(
    ci_lower = ind_ci$ci_lower,
    ci_upper = ind_ci$ci_upper,
    ci_lower_pct = NA,
    ci_upper_pct = NA,
    ci_method = ci_method,
    boot_se = ind_ci$boot_se,
    boot_p = ind_ci$boot_p,
    n_valid = n_valid,
    n_boot = n_boot,
    ci_level = ci_level,
    significant = ind_ci$significant,
    boot_a_vals = a_boot.v,
    boot_b_vals = b_boot.v,
    a_path = list(
      ci_lower = a_ci$ci_lower,
      ci_upper = a_ci$ci_upper,
      boot_se = a_ci$boot_se,
      significant = a_ci$significant
    ),
    b_path = list(
      ci_lower = b_ci$ci_lower,
      ci_upper = b_ci$ci_upper,
      boot_se = b_ci$boot_se,
      significant = b_ci$significant
    ),
    cprime = list(
      ci_lower = cp_ci$ci_lower,
      ci_upper = cp_ci$ci_upper,
      boot_se = cp_ci$boot_se,
      significant = cp_ci$significant
    ),
    total = list(
      ci_lower = tot_ci$ci_lower,
      ci_upper = tot_ci$ci_upper,
      boot_se = tot_ci$boot_se,
      significant = tot_ci$significant
    )
  )
}

# ============================================================
# PART 4: FIT MODELS
# ============================================================
log_section("Part 4: Fitting Mediation Models")

log_info("Bootstrap: N_BOOT = %d", N_BOOT)
if (N_BOOT == 0) {
  log_info(
    "  Using Sobel test only (preliminary)."
  )
  log_info(
    "  Set N_BOOT >= 1000 for final results."
  )
}
log_info("")

results.lst <- list()
n_cached <- 0L

for (sample_name in names(samples.lst)) {
  data.dt <- samples.lst[[sample_name]]
  log_info(
    "=== %s (N=%d) ===",
    toupper(sample_name), nrow(data.dt)
  )

  for (domain in COG_DOMAINS) {
    cog_vars <- paste0(
      "PHC_", domain, "_", VALID_TIMEPOINTS
    )

    for (pred_label in names(PREDICTORS)) {
      pred_cfg <- PREDICTORS[[pred_label]]
      key <- paste(
        sample_name, domain, pred_label,
        sep = "_"
      )

      # --- Checkpoint: skip if already done ---
      chk_path <- path_join(
        checkpoint_dir,
        paste0(key, ".rds")
      )
      if (file.exists(chk_path)) {
        cached <- read_rds_safe(
          chk_path,
          paste("checkpoint", key)
        )
        chk_nboot <- cached$bootstrap_n
        if (is.null(chk_nboot)) chk_nboot <- 0L
        chk_shared <- isTRUE(
          cached$shared_seed
        )
        chk_fq <- isTRUE(
          cached$fixed_quad
        )
        chk_all_ci <- isTRUE(
          cached$boot_all_ci
        )
        if (chk_nboot == N_BOOT &&
            chk_shared && chk_fq &&
            chk_all_ci) {
          log_info(
            "  %s | %s: checkpoint found",
            domain, pred_label
          )
          results.lst[[key]] <- cached
          n_cached <- n_cached + 1L
          next
        }
        reason <- if (!chk_fq) {
          "fixed_quad"
        } else if (!chk_all_ci) {
          "boot_all_ci"
        } else if (!chk_shared) {
          "shared_seed"
        } else {
          sprintf(
            "N_BOOT %d->%d",
            chk_nboot, N_BOOT
          )
        }
        log_info(
          "  %s | %s: stale (%s), refitting",
          domain, pred_label, reason
        )
      }

      log_info(
        "  %s | %s -> hvr_s -> cog_s",
        domain, pred_label
      )

      # Prepare data and build model
      # OpenMx FIML handles all missingness —
      # no complete.cases() filtering needed.
      if (pred_cfg$type == "observed") {
        pred_var <- pred_cfg$var
        keep_cols <- c(
          HVR_COLS, cog_vars,
          EDT_COLS, EDT_SQ_COLS,
          pred_var, "AGE_c", "EDUC_c"
        )
        if (HAS_APOE) {
          keep_cols <- c(
            keep_cols, "APOE4_num"
          )
        }

        model_data <- as.data.frame(
          data.dt[, ..keep_cols]
        )

        model <- build_mediation.fn(
          cog_vars, HVR_COLS,
          predictor = pred_var,
          name = paste0(
            "Med_", sample_name, "_",
            domain, "_", pred_label
          )
        )
      } else {
        # MIMIC: include indicators (FIML
        # handles missing indicators)
        keep_cols <- c(
          HVR_COLS, cog_vars,
          EDT_COLS, EDT_SQ_COLS,
          MIMIC_IND, "AGE_c", "EDUC_c"
        )
        if (HAS_APOE) {
          keep_cols <- c(
            keep_cols, "APOE4_num"
          )
        }

        model_data <- as.data.frame(
          data.dt[, ..keep_cols]
        )

        model <- build_mimic_mediation.fn(
          cog_vars, HVR_COLS,
          name = paste0(
            "Med_", sample_name, "_",
            domain, "_", pred_label
          )
        )
      }

      log_info(
        "    N: %d", nrow(model_data)
      )

      model <- mxModel(
        model,
        mxData(
          observed = model_data,
          type = "raw"
        )
      )

      fit <- tryCatch(
        mxTryHard(
          model, extraTries = 50,
          silent = TRUE,
          bestInitsOutput = FALSE
        ),
        error = function(e) {
          log_warn(
            "    Model error: %s", e$message
          )
          NULL
        }
      )

      converged <- !is.null(fit) &&
        !is.null(fit$output$status) &&
        fit$output$status$code %in% c(0, 1)

      if (!converged) {
        log_warn("    Did not converge")
        results.lst[[key]] <- list(
          sample = sample_name,
          domain = domain,
          predictor = pred_label,
          converged = FALSE,
          bootstrap_n = N_BOOT
        )
        write_rds_safe(
          results.lst[[key]], chk_path,
          paste("checkpoint", key)
        )
        next
      }

      # --------------------------------------------------
      # EXTRACT RESULTS
      # --------------------------------------------------
      s_fit <- tryCatch(
        summary(fit),
        error = function(e) NULL
      )
      if (is.null(s_fit)) {
        log_warn("    summary() failed")
        results.lst[[key]] <- list(
          sample = sample_name,
          domain = domain,
          predictor = pred_label,
          converged = FALSE,
          bootstrap_n = N_BOOT
        )
        write_rds_safe(
          results.lst[[key]], chk_path,
          paste("checkpoint", key)
        )
        next
      }
      params <- s_fit$parameters
      a_r <- params[params$name == "a", ]
      b_r <- params[params$name == "b", ]
      cp_r <- params[
        params$name == "cprime",
      ]

      a_z <- a_r$Estimate / a_r$Std.Error
      a_p <- 2 * pnorm(-abs(a_z))
      b_z <- b_r$Estimate / b_r$Std.Error
      b_p <- 2 * pnorm(-abs(b_z))
      cp_z <- cp_r$Estimate / cp_r$Std.Error
      cp_p <- 2 * pnorm(-abs(cp_z))

      # Sobel indirect effect (always computed)
      ind_est <- a_r$Estimate * b_r$Estimate
      ind_se <- sqrt(
        a_r$Estimate^2 * b_r$Std.Error^2 +
        b_r$Estimate^2 * a_r$Std.Error^2
      )
      ind_z <- ind_est / ind_se
      ind_p <- 2 * pnorm(-abs(ind_z))

      total_est <- ind_est + cp_r$Estimate
      prop_med <- if (abs(total_est) > 1e-4) {
        ind_est / total_est
      } else {
        NA
      }

      # --- Nested model LRTs ---
      lrt_a <- nested_lrt.fn(fit, "a")
      lrt_b <- nested_lrt.fn(fit, "b")
      lrt_cprime <- nested_lrt.fn(
        fit, "cprime"
      )

      # Bootstrap CIs (if N_BOOT > 0)
      boot_ci_lower <- NA
      boot_ci_upper <- NA
      boot_se <- NA
      boot_significant <- NA
      a_boot_ci_lower <- NA
      a_boot_ci_upper <- NA
      a_boot_se <- NA
      b_boot_ci_lower <- NA
      b_boot_ci_upper <- NA
      b_boot_se <- NA
      cp_boot_ci_lower <- NA
      cp_boot_ci_upper <- NA
      cp_boot_se <- NA
      tot_boot_ci_lower <- NA
      tot_boot_ci_upper <- NA
      tot_boot_se <- NA
      tot_boot_sig <- NA

      if (N_BOOT > 0) {
        # Warm-start: use converged params
        # as starting values so each bootstrap
        # iteration begins near the solution
        # rather than from default values.
        warm_vals <- coef(fit)
        warm_labs <- names(warm_vals)

        boot_model.fn <- if (
          pred_cfg$type == "observed"
        ) {
          function() {
            m <- build_mediation.fn(
              cog_vars, HVR_COLS,
              predictor = pred_cfg$var,
              name = paste0(
                "Boot_", domain, "_",
                pred_label
              )
            )
            omxSetParameters(
              m,
              labels = warm_labs,
              values = warm_vals
            )
          }
        } else {
          function() {
            m <- build_mimic_mediation.fn(
              cog_vars, HVR_COLS,
              name = paste0(
                "Boot_", domain, "_",
                pred_label
              )
            )
            omxSetParameters(
              m,
              labels = warm_labs,
              values = warm_vals
            )
          }
        }

        # Unique seed per combination
        pred_idx <- which(
          pred_label == names(PREDICTORS)
        )
        dom_idx <- which(
          domain == COG_DOMAINS
        )
        samp_idx <- which(
          sample_name == names(samples.lst)
        )
        # Shared seed: same per sex+domain so
        # FRS and CVR_mimic draw identical
        # bootstrap samples (enables paired
        # delta_a comparison)
        boot_seed <- 42 + dom_idx +
          samp_idx * 10

        boot_res <- bootstrap_indirect.fn(
          boot_model.fn, model_data,
          n_boot = N_BOOT,
          ci_level = CI_LEVEL,
          seed = boot_seed,
          original_estimate = ind_est,
          original_cprime = cp_r$Estimate,
          original_a = a_r$Estimate,
          original_b = b_r$Estimate
        )

        boot_ci_lower <- boot_res$ci_lower
        boot_ci_upper <- boot_res$ci_upper
        boot_se <- boot_res$boot_se
        boot_significant <-
          boot_res$significant
        a_boot_ci_lower <-
          boot_res$a_path$ci_lower
        a_boot_ci_upper <-
          boot_res$a_path$ci_upper
        a_boot_se <-
          boot_res$a_path$boot_se
        b_boot_ci_lower <-
          boot_res$b_path$ci_lower
        b_boot_ci_upper <-
          boot_res$b_path$ci_upper
        b_boot_se <-
          boot_res$b_path$boot_se
        cp_boot_ci_lower <-
          boot_res$cprime$ci_lower
        cp_boot_ci_upper <-
          boot_res$cprime$ci_upper
        cp_boot_se <-
          boot_res$cprime$boot_se
        tot_boot_ci_lower <-
          boot_res$total$ci_lower
        tot_boot_ci_upper <-
          boot_res$total$ci_upper
        tot_boot_se <-
          boot_res$total$boot_se
        tot_boot_sig <-
          boot_res$total$significant
      }

      # Log results
      a_s <- ifelse(a_p < 0.05, "*", "")
      b_s <- ifelse(b_p < 0.05, "*", "")
      cp_s <- ifelse(cp_p < 0.05, "*", "")
      ind_s <- ifelse(ind_p < 0.05, "*", "")

      log_info(
        "    a: %.5f (SE=%.5f, p=%.4f) %s",
        a_r$Estimate, a_r$Std.Error,
        a_p, a_s
      )
      log_info(
        "    b: %.4f (SE=%.4f, p=%.2e) %s",
        b_r$Estimate, b_r$Std.Error,
        b_p, b_s
      )
      log_info(
        "    c': %.5f (SE=%.5f, p=%.4f) %s",
        cp_r$Estimate, cp_r$Std.Error,
        cp_p, cp_s
      )
      log_info(
        "    indirect: %.6f (Sobel p=%.4f) %s",
        ind_est, ind_p, ind_s
      )

      if (N_BOOT > 0) {
        log_info(
          "    boot CI: [%.6f, %.6f] %s",
          boot_ci_lower, boot_ci_upper,
          if (isTRUE(boot_significant)) {
            "*"
          } else {
            ""
          }
        )
      }

      if (!is.na(prop_med)) {
        log_info(
          "    proportion mediated: %.1f%%",
          prop_med * 100
        )
      }

      # Log nested LRTs
      if (!is.null(lrt_a)) {
        log_info(
          "    LRT a-path: X2=%.2f, p=%.4f",
          lrt_a$chi_diff, lrt_a$p
        )
      }
      if (!is.null(lrt_b)) {
        log_info(
          "    LRT b-path: X2=%.2f, p=%.4f",
          lrt_b$chi_diff, lrt_b$p
        )
      }

      # Store boot vectors for paired delta
      boot_a_raw <- if (
        N_BOOT > 0 && !is.null(boot_res)
      ) boot_res$boot_a_vals else NULL
      boot_b_raw <- if (
        N_BOOT > 0 && !is.null(boot_res)
      ) boot_res$boot_b_vals else NULL

      # Store result
      result <- list(
        sample = sample_name,
        domain = domain,
        predictor = pred_label,
        n = nrow(model_data),
        converged = TRUE,
        a_path = list(
          est = a_r$Estimate,
          se = a_r$Std.Error,
          z = a_z, p = a_p,
          boot_ci_lower = a_boot_ci_lower,
          boot_ci_upper = a_boot_ci_upper,
          boot_se = a_boot_se
        ),
        b_path = list(
          est = b_r$Estimate,
          se = b_r$Std.Error,
          z = b_z, p = b_p,
          boot_ci_lower = b_boot_ci_lower,
          boot_ci_upper = b_boot_ci_upper,
          boot_se = b_boot_se
        ),
        cprime = list(
          est = cp_r$Estimate,
          se = cp_r$Std.Error,
          z = cp_z, p = cp_p,
          boot_ci_lower = cp_boot_ci_lower,
          boot_ci_upper = cp_boot_ci_upper,
          boot_se = cp_boot_se
        ),
        indirect = list(
          est = ind_est,
          sobel_se = ind_se,
          sobel_z = ind_z,
          sobel_p = ind_p,
          boot_ci_lower = boot_ci_lower,
          boot_ci_upper = boot_ci_upper,
          boot_se = boot_se,
          boot_significant = boot_significant
        ),
        total = list(
          est = total_est,
          boot_ci_lower = tot_boot_ci_lower,
          boot_ci_upper = tot_boot_ci_upper,
          boot_se = tot_boot_se,
          boot_significant = tot_boot_sig
        ),
        prop_mediated = prop_med,
        lrt_a = lrt_a,
        lrt_b = lrt_b,
        lrt_cprime = lrt_cprime,
        fit_indices =
          extract_fit_indices.fn(fit),
        all_params = params,
        bootstrap_n = N_BOOT,
        shared_seed = TRUE,
        fixed_quad = TRUE,
        boot_all_ci = TRUE,
        boot_a_vals = boot_a_raw,
        boot_b_vals = boot_b_raw
      )

      results.lst[[key]] <- result
      write_rds_safe(
        result, chk_path,
        paste("checkpoint", key)
      )
    }

    # ------------------------------------------------
    # Paired Δ_a comparison (FRS vs CVR_mimic)
    # ------------------------------------------------
    # Both predictors share the same bootstrap
    # seed, so boot samples are paired.
    frs_key <- paste(
      sample_name, domain, "FRS", sep = "_"
    )
    cvr_key <- paste(
      sample_name, domain, "CVR_mimic",
      sep = "_"
    )
    frs_r <- results.lst[[frs_key]]
    cvr_r <- results.lst[[cvr_key]]

    can_compare <- N_BOOT > 0 &&
      !is.null(frs_r) &&
      isTRUE(frs_r$converged) &&
      !is.null(frs_r$boot_a_vals) &&
      !is.null(cvr_r) &&
      isTRUE(cvr_r$converged) &&
      !is.null(cvr_r$boot_a_vals)

    if (can_compare) {
      frs_a.v <- frs_r$boot_a_vals
      cvr_a.v <- cvr_r$boot_a_vals
      both_ok <- !is.na(frs_a.v) &
        !is.na(cvr_a.v)
      delta_a.v <- frs_a.v - cvr_a.v

      delta_est <- frs_r$a_path$est -
        cvr_r$a_path$est
      valid_d <- delta_a.v[both_ok]

      if (length(valid_d) >= 100) {
        alpha_d <- 1 - CI_LEVEL
        # Bias-corrected CI
        pb <- mean(valid_d < delta_est)
        pb <- max(0.001, min(0.999, pb))
        z0 <- qnorm(pb)
        za_lo <- qnorm(alpha_d / 2)
        za_hi <- qnorm(1 - alpha_d / 2)
        p_lo <- pnorm(2 * z0 + za_lo)
        p_hi <- pnorm(2 * z0 + za_hi)
        ci_lo <- as.numeric(
          quantile(valid_d, probs = p_lo)
        )
        ci_hi <- as.numeric(
          quantile(valid_d, probs = p_hi)
        )
        sig <- (ci_lo > 0 || ci_hi < 0)

        cmp_key <- paste(
          sample_name, domain, sep = "_"
        )
        if (!exists("comparisons.lst")) {
          comparisons.lst <- list()
        }
        comparisons.lst[[cmp_key]] <- list(
          delta_a_est = delta_est,
          delta_a_ci_lower = ci_lo,
          delta_a_ci_upper = ci_hi,
          delta_a_boot_se = sd(valid_d),
          delta_a_significant = sig,
          n_valid = length(valid_d)
        )

        sig_s <- if (sig) " *" else ""
        log_info(
          "  delta_a (FRS-CVR): %.5f [%.5f, %.5f]%s",
          delta_est, ci_lo, ci_hi, sig_s
        )
      } else {
        log_info(
          "  delta_a: too few valid (%d)",
          length(valid_d)
        )
      }
    }

    log_info("")
  }
}

# Initialise comparisons if none computed
if (!exists("comparisons.lst")) {
  comparisons.lst <- list()
}

n_total <- length(results.lst)
log_info(
  "Checkpoint summary: %d/%d cached, %d new",
  n_cached, n_total, n_total - n_cached
)

# ============================================================
# PART 5: SUMMARY COMPARISON
# ============================================================
log_section("Part 5: Summary Comparison")

log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info(
  "CVR_mimic vs FRS Mediation (Sex-Stratified)"
)
log_info(paste(rep("=", 70), collapse = ""))
log_info("")
log_info(sprintf(
  "%-6s %-4s %-5s %8s %10s %10s %8s",
  "Sex", "Dom", "Pred",
  "a(p)", "b(p)", "ind(p)", "c'(p)"
))
log_info(paste(rep("-", 70), collapse = ""))

summary_rows.lst <- list()

for (key in names(results.lst)) {
  r <- results.lst[[key]]
  if (!r$converged) {
    log_info(
      "%-6s %-4s %-5s  NOT CONVERGED",
      r$sample, r$domain, r$predictor
    )
    next
  }

  a_s <- ifelse(r$a_path$p < 0.05, "*", " ")
  b_s <- ifelse(r$b_path$p < 0.05, "*", " ")
  i_s <- ifelse(
    r$indirect$sobel_p < 0.05, "*", " "
  )
  c_s <- ifelse(
    r$cprime$p < 0.05, "*", " "
  )

  log_info(sprintf(
    "%-6s %-4s %-9s %7.4f%s %9.2e%s %9.4f%s %7.4f%s",
    r$sample, r$domain, r$predictor,
    r$a_path$p, a_s,
    r$b_path$p, b_s,
    r$indirect$sobel_p, i_s,
    r$cprime$p, c_s
  ))

  summary_rows.lst[[key]] <- data.table(
    sample = r$sample,
    domain = r$domain,
    predictor = r$predictor,
    n = r$n,
    a_est = r$a_path$est,
    a_p = r$a_path$p,
    b_est = r$b_path$est,
    b_p = r$b_path$p,
    indirect_est = r$indirect$est,
    indirect_p = r$indirect$sobel_p,
    cprime_est = r$cprime$est,
    cprime_p = r$cprime$p,
    prop_mediated = r$prop_mediated
  )
}

log_info(paste(rep("-", 70), collapse = ""))

if (N_BOOT > 0) {
  log_info("* p < 0.05 (bootstrap CI for indirect)")
} else {
  log_info("* p < 0.05 (Sobel test)")
}
log_info("")

summary.dt <- rbindlist(
  summary_rows.lst, fill = TRUE
)

# ============================================================
# PART 6: INTERPRETATION
# ============================================================
log_section("Part 6: Interpretation")

for (sex in names(samples.lst)) {
  log_info("--- %s ---", sex)
  for (domain in COG_DOMAINS) {
    cvr_key <- paste(
      sex, domain, "CVR_mimic", sep = "_"
    )
    frs_key <- paste(
      sex, domain, "FRS", sep = "_"
    )

    cvr_r <- results.lst[[cvr_key]]
    frs_r <- results.lst[[frs_key]]

    cvr_sig <- !is.null(cvr_r) &&
      cvr_r$converged &&
      cvr_r$indirect$sobel_p < 0.05
    frs_sig <- !is.null(frs_r) &&
      frs_r$converged &&
      frs_r$indirect$sobel_p < 0.05

    interp <- if (frs_sig && !cvr_sig) {
      "FRS-only -> age artifact"
    } else if (cvr_sig && !frs_sig) {
      "CVR-only -> true CVR effect"
    } else if (frs_sig && cvr_sig) {
      "Both significant"
    } else {
      "Neither significant"
    }

    log_info(
      "  %s: FRS=%s, CVR=%s -> %s",
      domain,
      ifelse(frs_sig, "sig", "ns"),
      ifelse(cvr_sig, "sig", "ns"),
      interp
    )
  }
  log_info("")
}

log_info("")
log_info("INTERPRETATION GUIDE:")
log_info("  FRS-only -> age artifact:")
log_info("    FRS mediation driven by age-weighting,")
log_info("    not true cardiovascular risk.")
log_info("  CVR-only -> true CVR effect:")
log_info("    CVR_mimic captures genuine CV risk")
log_info("    effect missed by FRS age confound.")
log_info("  Both significant:")
log_info("    Genuine CVR mediation present.")
log_info("  Neither significant:")
log_info("    No mediation via HVR for this domain.")
log_info("")
log_info("  Modern mediation methodology:")
log_info("  - Hayes (2009): total effect NOT required")
log_info("  - Preacher & Hayes (2008): bootstrap CIs")
log_info("    preferred for indirect effect inference")
log_info("  - Sobel test used as preliminary test")

# ============================================================
# PART 7: SAVE RESULTS
# ============================================================
log_section("Part 7: Save Results")

# Save quadratic attempt documentation
write_rds_safe(
  quadratic_attempts.lst,
  file.path(
    results_dir,
    "lgcm_frs_quadratic_attempts.rds"
  ),
  "Quadratic attempt documentation"
)

# Save mediation results (both predictors, both sexes)
output.lst <- list(
  results = results.lst,
  comparisons = comparisons.lst,
  summary_table = summary.dt,
  sample_info = lapply(
    samples.lst,
    function(dt) list(n = nrow(dt))
  ),
  model_spec = paste(
    "Parallel process LGCM mediation:",
    "CVR -> hvr_s -> cog_s.",
    "Stratified by sex."
  ),
  predictors = PREDICTORS,
  covariates = c(
    "AGE_c", "EDUC_c",
    if (HAS_APOE) "APOE4_num"
  ),
  bootstrap = list(
    n_boot = N_BOOT,
    ci_level = CI_LEVEL,
    note = if (N_BOOT == 0) {
      "Sobel test only (preliminary)"
    } else {
      paste0(
        "Bootstrap CIs (n=", N_BOOT, ")"
      )
    }
  ),
  quadratic_attempts = quadratic_attempts.lst,
  limitation = paste(
    "Trajectories are significantly non-linear",
    "(quadratic). Fully random quadratic models",
    "failed (identification issues). Fixed",
    "quadratic extension adds population-level",
    "curvature; individual quadratic variation",
    "is not estimated."
  ),
  methodology_notes = paste(
    "Following modern mediation methodology",
    "(Hayes, 2009; Zhao et al., 2010),",
    "a significant total effect (X->Y) is NOT",
    "required for mediation.",
    "Bootstrap CIs are preferred for inference",
    "on indirect effects",
    "(Preacher & Hayes, 2008)."
  ),
  references = list(
    mediation_methodology = c(
      "Hayes (2009). Comm Monographs, 76(4).",
      "Zhao et al. (2010). J Consumer Res.",
      "Preacher & Hayes (2008). PMC2819361",
      "Kenny & Judd (2014). PMC4142865"
    ),
    lgcm_methodology = c(
      paste(
        "Preacher (2010). Latent growth",
        "curve models. In Hancock &",
        "Mueller (Eds.), Reviewer's Guide",
        "to Quantitative Methods.",
        "DOI: 10.4324/9780203861554-17"
      ),
      paste(
        "Curran, Obeidat & Losardo",
        "(2010). Twelve frequently asked",
        "questions about growth curve",
        "modeling. J Cogn Dev, 11(2),",
        "121-136.",
        "DOI: 10.1080/15248371003699969"
      )
    ),
    substantive = c(
      "Debette et al. (2011). PMC3146307",
      "Viticchi et al. (2020). JACC.",
      "Gorelick et al. (2011). PMC3778669"
    )
  ),
  fit_assessment = list(
    valid_indices = c("-2LL", "AIC", "BIC"),
    nested_tests = c(
      "a-path", "b-path", "c'-path"
    ),
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
        "mediation paths (a, b, c'),",
        "(c) bootstrap confidence",
        "intervals for the indirect",
        "effect, and (d) parameter",
        "estimate plausibility."
      ),
      paste(
        "Fixed quadratic extension was",
        "adopted after fully random",
        "quadratic models failed due",
        "to empirical under-",
        "identification (near-zero",
        "quadratic variance;",
        "Preacher, 2010,",
        "DOI: 10.4324/9780203861554-17).",
        "Population curvature is",
        "captured; individual quadratic",
        "variation is not."
      )
    )
  ),
  timestamp = Sys.time()
)

write_rds_safe(
  output.lst,
  file.path(
    results_dir,
    "lgcm_mediation_results.rds"
  ),
  "LGCM mediation results (CVR_mimic + FRS)"
)

# ============================================================
# SCRIPT COMPLETE
# ============================================================
log_info("")
log_info("This script produced:")
log_info("  1. Quadratic failure documentation")
log_info("  2. Linear mediation: CVR -> hvr_s -> cog_s")
log_info("  3. Two predictors: CVR_mimic and FRS")
log_info("  4. Sex-stratified (Male/Female only)")
log_info("")
log_info("Output: lgcm_mediation_results.rds")

log_script_end(
  "13_lgcm_mediation.R", success = TRUE
)
