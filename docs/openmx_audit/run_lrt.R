#!/usr/bin/env Rscript

# =============================================================================
# Isolated LRT: homoscedastic vs heteroscedastic HVR residuals
# =============================================================================
# Read-only against existing pipeline outputs.
# Writes ONLY to docs/openmx_audit/.
# Does not modify R/scripts/ or any cached RDS/checkpoints.
#
# Cell: Male, MEM domain, CVR_mimic predictor (primary manuscript cell).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(OpenMx)
})

source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/disease_time.R"))
source(here("R/utils/logging.R"))

mxOption(NULL, "Default optimizer", "SLSQP")
mxOption(NULL, "Number of Threads", 2)

OUT_DIR <- here("docs/openmx_audit")
OUT_FILE <- file.path(OUT_DIR, "lrt_Male_MEM_CVR_mimic.rds")
LOG_FILE <- file.path(OUT_DIR, "lrt_Male_MEM_CVR_mimic.log")

# Tee output to log file
sink(LOG_FILE, split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("=== LRT: HVR residual structure ===\n")
cat("Timestamp:", format(Sys.time()), "\n")
cat("Cell: Male, MEM, CVR_mimic\n\n")

# -----------------------------------------------------------------------------
# Load data (read-only)
# -----------------------------------------------------------------------------
male.dt <- readRDS(here("data/derivatives/lgcm/lgcm_input_male.rds"))
cat("Loaded male sample: N =", nrow(male.dt), "\n")

VALID_TIMEPOINTS <- paste0("T", 0:5)
EDT_COLS <- paste0("EDT_", VALID_TIMEPOINTS)
EDT_SQ_COLS <- paste0(EDT_COLS, "_sq")
HVR_COLS <- paste0("HVR_Z_", VALID_TIMEPOINTS)
COG_COLS <- paste0("PHC_MEM_", VALID_TIMEPOINTS)
COG_VAR_COLS <- paste0("PHC_MEM_VAR_", VALID_TIMEPOINTS)

# -----------------------------------------------------------------------------
# Reproduce script 13 prep (inline; no writes)
# -----------------------------------------------------------------------------
male.dt <- copy(male.dt)  # do not modify the loaded object in place

# EDT extrapolation
extrap <- extrapolate_edt_hybrid(male.dt, valid_timepoints = VALID_TIMEPOINTS)
male.dt <- extrap$data

n_edt_na <- sum(is.na(as.matrix(male.dt[, ..EDT_COLS])))
cat("After EDT extrapolation, remaining NAs:", n_edt_na, "\n")
if (n_edt_na > 0) stop("EDT NAs remain after extrapolation; OpenMx will fail.")

# Center EDT (script 13 uses pooled male+female grand mean; here we use male-only
# because we are doing a within-cell LRT, where the chi-square difference is
# invariant to additive shifts in the time scale. The two models share centering.)
edt_mean <- mean(unlist(male.dt[, ..EDT_COLS]), na.rm = TRUE)
cat("EDT centering mean (male-only, for LRT internal consistency):",
    round(edt_mean, 3), "\n")
for (col in EDT_COLS) {
  set(male.dt, j = col, value = male.dt[[col]] - edt_mean)
}
for (i in seq_along(EDT_COLS)) {
  set(male.dt, j = EDT_SQ_COLS[i], value = male.dt[[EDT_COLS[i]]]^2)
}

# Centered predictors
age_mean <- mean(male.dt$AGE_bl, na.rm = TRUE)
educ_mean <- mean(male.dt$EDUC, na.rm = TRUE)
set(male.dt, j = "AGE_c", value = male.dt$AGE_bl - age_mean)
set(male.dt, j = "EDUC_c", value = male.dt$EDUC - educ_mean)

# APOE
HAS_APOE <- "APOE4" %in% names(male.dt) && sum(!is.na(male.dt$APOE4)) > 0
if (HAS_APOE) male.dt[, APOE4_num := as.numeric(APOE4)]
cat("APOE4 available:", HAS_APOE, "\n\n")

# -----------------------------------------------------------------------------
# Local copy of build_mediation.fn with a switch for residual labels
# Mirrors R/scripts/13_lgcm_mediation.R::build_mediation.fn
# -----------------------------------------------------------------------------
build_mediation.fn <- function(cog_vars, hvr_vars,
                               predictor = "CVR_c",
                               name = "Med",
                               data = NULL,
                               skip_age_cov = FALSE,
                               hvr_resid_mode = c("homoscedastic",
                                                  "heteroscedastic")) {
  hvr_resid_mode <- match.arg(hvr_resid_mode)
  n_tp <- length(cog_vars)
  latents <- c("hvr_i", "hvr_s", "cog_i", "cog_s")

  hvr_int <- mxPath(from = "hvr_i", to = hvr_vars,
                    free = FALSE, values = rep(1, n_tp))
  hvr_slope <- mxPath(from = "hvr_s", to = hvr_vars,
                      free = FALSE,
                      labels = paste0("data.", EDT_COLS))
  cog_int <- mxPath(from = "cog_i", to = cog_vars,
                    free = FALSE, values = rep(1, n_tp))
  cog_slope <- mxPath(from = "cog_s", to = cog_vars,
                      free = FALSE,
                      labels = paste0("data.", EDT_COLS))

  lat_means <- mxPath(from = "one", to = latents, free = TRUE,
                      values = c(0, -0.1, 0, -0.1),
                      labels = paste0("mean_", latents))
  lat_vars <- mxPath(from = latents, arrows = 2, free = TRUE,
                     values = c(1, 0.1, 1, 0.1),
                     labels = paste0("var_", latents))
  hvr_cov <- mxPath(from = "hvr_i", to = "hvr_s", arrows = 2,
                    free = TRUE, values = 0, labels = "cov_hvr_is")
  cog_cov <- mxPath(from = "cog_i", to = "cog_s", arrows = 2,
                    free = TRUE, values = 0, labels = "cov_cog_is")
  cross_int <- mxPath(from = "hvr_i", to = "cog_i", arrows = 2,
                      free = TRUE, values = 0, labels = "cov_hvr_i_cog_i")

  # HVR residuals: the variant under test
  hvr_resid_labels <- if (hvr_resid_mode == "homoscedastic") {
    "hvr_resid"
  } else {
    paste0("hvr_resid_", VALID_TIMEPOINTS)
  }
  hvr_resid <- mxPath(from = hvr_vars, arrows = 2,
                      free = TRUE, values = 0.5,
                      labels = hvr_resid_labels)

  # Cognitive residuals: fixed to SE² from CFA (matches script 13)
  cog_var_cols.v <- sub("_(T\\d+)$", "_VAR_\\1", cog_vars)
  avg_cog_se2.v <- sapply(cog_var_cols.v, function(col) {
    if (!is.null(data) && col %in% names(data)) {
      mean(data[[col]], na.rm = TRUE)
    } else {
      0.04
    }
  })
  cog_resid <- mxPath(from = cog_vars, arrows = 2,
                      free = FALSE, values = avg_cog_se2.v)

  hvr_mm <- mxPath(from = "one", to = hvr_vars,
                   free = FALSE, values = rep(0, n_tp))
  cog_mm <- mxPath(from = "one", to = cog_vars,
                   free = FALSE, values = rep(0, n_tp))

  preds <- if (skip_age_cov) {
    c(predictor, "EDUC_c")
  } else {
    c(predictor, "AGE_c", "EDUC_c")
  }
  if (HAS_APOE) preds <- c(preds, "APOE4_num")

  pred_var <- mxPath(from = preds, arrows = 2,
                     free = FALSE, values = 1)
  pred_mn <- mxPath(from = "one", to = preds,
                    free = FALSE, values = 0)

  a_path <- mxPath(from = predictor, to = "hvr_s",
                   free = TRUE, values = 0, labels = "a")
  b_path <- mxPath(from = "hvr_s", to = "cog_s",
                   free = TRUE, values = 0.5, labels = "b")
  cprime_path <- mxPath(from = predictor, to = "cog_s",
                        free = TRUE, values = 0, labels = "cprime")
  pred_int <- mxPath(from = predictor,
                     to = c("hvr_i", "cog_i"),
                     free = TRUE, values = 0,
                     labels = c("pred_hvr_i", "pred_cog_i"))

  educ_paths <- mxPath(from = "EDUC_c", to = latents,
                       free = TRUE, values = 0,
                       labels = paste0("educ_", latents))

  model <- mxModel(name, type = "RAM",
                   manifestVars = c(hvr_vars, cog_vars, preds),
                   latentVars = latents)
  model <- mxModel(model,
                   hvr_int, hvr_slope, cog_int, cog_slope,
                   lat_means, lat_vars,
                   hvr_cov, cog_cov, cross_int,
                   hvr_resid, cog_resid,
                   hvr_mm, cog_mm,
                   pred_var, pred_mn,
                   a_path, b_path, cprime_path, pred_int,
                   educ_paths)

  if (!skip_age_cov) {
    age_paths <- mxPath(from = "AGE_c", to = latents,
                        free = TRUE, values = 0,
                        labels = paste0("age_", latents))
    model <- mxModel(model, age_paths)
  }
  if (HAS_APOE) {
    apoe_paths <- mxPath(from = "APOE4_num", to = latents,
                         free = TRUE, values = 0,
                         labels = paste0("apoe_", latents))
    model <- mxModel(model, apoe_paths)
  }

  # Fixed quadratic extension
  hvr_q_paths <- mxPath(from = "hvr_q", to = hvr_vars,
                        free = FALSE,
                        labels = paste0("data.", EDT_SQ_COLS))
  cog_q_paths <- mxPath(from = "cog_q", to = cog_vars,
                        free = FALSE,
                        labels = paste0("data.", EDT_SQ_COLS))
  quad_var <- mxPath(from = c("hvr_q", "cog_q"),
                     arrows = 2, free = FALSE,
                     values = c(1e-10, 1e-10))
  quad_means <- mxPath(from = "one",
                       to = c("hvr_q", "cog_q"),
                       free = TRUE,
                       values = c(-0.01, -0.01),
                       labels = c("mean_hvr_q", "mean_cog_q"))
  model <- mxModel(model,
                   latentVars = c("hvr_q", "cog_q"),
                   hvr_q_paths, cog_q_paths,
                   quad_var, quad_means)

  indirect <- mxAlgebra(a * b, name = "indirect")
  total <- mxAlgebra(indirect + cprime, name = "total")
  model <- mxModel(model, indirect, total)

  model
}

# -----------------------------------------------------------------------------
# Fit both variants
# -----------------------------------------------------------------------------
keep_cols <- unique(c(HVR_COLS, COG_COLS, EDT_COLS, EDT_SQ_COLS,
                      "CVR_sex_centered", "AGE_c", "EDUC_c"))
if (HAS_APOE) keep_cols <- c(keep_cols, "APOE4_num")
model_data <- as.data.frame(male.dt[, ..keep_cols])

# Match script 13's predictor centering: CVR is centered within sample
frs_m <- mean(male.dt$CVR_sex_centered, na.rm = TRUE)
model_data$CVR_c <- model_data$CVR_sex_centered - frs_m
model_data$CVR_sex_centered <- NULL  # not used in model

fit_variant <- function(mode) {
  cat(sprintf("\n--- Fitting variant: %s ---\n", mode))
  t0 <- Sys.time()
  model <- build_mediation.fn(
    cog_vars = COG_COLS,
    hvr_vars = HVR_COLS,
    predictor = "CVR_c",
    name = paste0("Med_Male_MEM_CVR_mimic_", mode),
    data = male.dt,
    skip_age_cov = FALSE,
    hvr_resid_mode = mode
  )
  model <- mxModel(model, mxData(observed = model_data, type = "raw"))
  fit <- mxTryHard(model, extraTries = 50, silent = TRUE,
                   bestInitsOutput = FALSE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  status: %d, time: %.1fs\n",
              fit$output$status$code, dt))
  fit
}

fit_homo <- fit_variant("homoscedastic")
fit_hete <- fit_variant("heteroscedastic")

# -----------------------------------------------------------------------------
# LRT and summary
# -----------------------------------------------------------------------------
summarize <- function(fit, label) {
  s <- summary(fit)
  list(
    label = label,
    status = fit$output$status$code,
    minus2LL = fit$output$fit,
    n_param = length(fit$output$estimate),
    n_obs = fit$data$numObs,
    AIC = fit$output$fit + 2 * length(fit$output$estimate),
    BIC = fit$output$fit + log(fit$data$numObs) *
      length(fit$output$estimate),
    a = s$parameters[s$parameters$name == "a",
                     c("Estimate", "Std.Error")],
    b = s$parameters[s$parameters$name == "b",
                     c("Estimate", "Std.Error")],
    cprime = s$parameters[s$parameters$name == "cprime",
                          c("Estimate", "Std.Error")]
  )
}

summ_homo <- summarize(fit_homo, "homoscedastic")
summ_hete <- summarize(fit_hete, "heteroscedastic")

chi_diff <- summ_homo$minus2LL - summ_hete$minus2LL
df_diff <- summ_hete$n_param - summ_homo$n_param
p_value <- pchisq(chi_diff, df = df_diff, lower.tail = FALSE)
aic_diff <- summ_homo$AIC - summ_hete$AIC
bic_diff <- summ_homo$BIC - summ_hete$BIC

cat("\n========================================\n")
cat("LRT RESULTS\n")
cat("========================================\n")
cat(sprintf("Homoscedastic  -2LL = %.3f  (n_param = %d)\n",
            summ_homo$minus2LL, summ_homo$n_param))
cat(sprintf("Heteroscedastic -2LL = %.3f  (n_param = %d)\n",
            summ_hete$minus2LL, summ_hete$n_param))
cat(sprintf("Chi-square diff = %.3f, df = %d, p = %.4g\n",
            chi_diff, df_diff, p_value))
cat(sprintf("dAIC (homo - hete) = %.3f\n", aic_diff))
cat(sprintf("dBIC (homo - hete) = %.3f\n", bic_diff))
cat(sprintf("Decision: %s\n",
            ifelse(p_value < 0.05,
                   "Heteroscedastic preferred (constraint rejected)",
                   "Homoscedastic defensible (constraint not rejected)")))

cat("\n--- Path estimates: homoscedastic ---\n")
print(summ_homo$a); print(summ_homo$b); print(summ_homo$cprime)
cat("\n--- Path estimates: heteroscedastic ---\n")
print(summ_hete$a); print(summ_hete$b); print(summ_hete$cprime)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
saveRDS(list(
  cell = "Male_MEM_CVR_mimic",
  homoscedastic = summ_homo,
  heteroscedastic = summ_hete,
  lrt = list(chi_diff = chi_diff, df = df_diff, p = p_value,
             aic_diff = aic_diff, bic_diff = bic_diff),
  fit_homo = fit_homo,
  fit_hete = fit_hete,
  timestamp = Sys.time()
), OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
