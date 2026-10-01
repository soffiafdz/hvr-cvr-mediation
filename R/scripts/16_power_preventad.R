#!/usr/bin/env Rscript

# =============================================================================
# 16_power_preventad.R
# =============================================================================
# Monte Carlo power analysis for a CIHR fellowship application.
#
# QUESTION:
#   What is the power to detect the CVR_mimic x HVR_z x YRS
#   three-way interaction (LME, cognition as outcome) in a
#   cohort with PREVENT-AD's design, given the effect sizes
#   and variance components estimated in this repo's ADNI
#   models? Secondary: power for the sex-by-moderation
#   contrast (SEX x CVR_mimic x HVR_z x YRS) in the pooled model.
#
# TARGET DESIGN (PREVENT-AD; Villeneuve et al. 2025,
# Alzheimers Dement, doi 10.1002/alz.70653):
#   N = 387 (279 women, 108 men), annual visits,
#   95% with >= 1 follow-up, 79% still followed beyond
#   10 years. Retention is interpolated geometrically between
#   year 1 (0.95) and year 10 (0.79) and extended beyond.
#
# EFFECT-SIZE ANCHORS (all parameters loaded from fits):
#   PRIMARY:   A-/CU pooled-sex LME, Memory, CVR_mimic 3-way
#              (models/results/lme/lme_amyloid_negative.rds).
#              This is the scientifically relevant anchor:
#              PREVENT-AD is younger, CU and pre-amyloid.
#   SECONDARY: A+ females LME, Executive Function, CVR_mimic
#              3-way (models/results/lme/lme_hvr_z_results.rds),
#              reported only as an upper bound.
#   SEX CONTRAST: the 4-way YRS:CVR_mimic:HVR_z:SEX term from a
#              sex-interaction refit of the A-/CU Memory model on
#              its saved analysis frame (fit in this script); an
#              A+ EXF refit (both sexes) serves as upper bound.
#
# DESIGNS: 5 / 8 / 10 / 12 annual waves with dropout, 8 waves
#   complete. The 12-wave design (years 0-11) is taken as the
#   design matching PREVENT-AD's real follow-up. Plus a
#   predictor-spread sensitivity and MDE grid scans.
#
# GENERATING MODEL (per anchor):
#   All fixed effects, random-effect covariance (intercept,
#   YRS slope), and residual variance from the fitted lmer
#   object. Baseline predictors (CVR_mimic, HVR_z) drawn from
#   a bivariate normal with the ADNI baseline means, SDs and
#   correlation. The target coefficient is set to the anchor
#   value (or a grid value for minimum-detectable-effect).
#
# FITTED MODEL (same as manuscript, scripts 08 / 10b):
#   Cog ~ YRS * CVR * HVR_z + I(YRS^2) * CVR * HVR_z +
#         Age_bl + EDUC + APOE4 [+ SEX] + (YRS | PTID)
#   Sex contrast adds "* SEX" to both interaction blocks.
#   REML, bobyqa, Satterthwaite p-value (lmerTest), alpha 0.05.
#
# simr is NOT installed in the renv; simulation is hand-rolled
# with lme4 / lmerTest (no packages installed).
#
# INPUTS:
#   - models/results/lme/lme_amyloid_negative.rds
#   - models/results/lme/lme_hvr_z_results.rds
#
# OUTPUTS:
#   - outputs/power_preventad.rds
#   - outputs/power_preventad.md
#
# ENV OVERRIDES (for smoke tests only):
#   POWER_N_REPS  replicates per condition (default 1000)
#   POWER_N_CORES parallel workers (default 20)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))

log_script_start("16_power_preventad.R")
config <- load_config()
SEED <- get_seed()

# -----------------------------------------------------------
# Constants
# -----------------------------------------------------------
N_REPS <- as.integer(Sys.getenv("POWER_N_REPS", "1000"))
N_CORES <- min(
  as.integer(Sys.getenv("POWER_N_CORES", "20")),
  parallel::detectCores()
)
ALPHA <- 0.05
TARGET_POWER <- 0.80
TERM_3WAY <- "YRS_from_bl:CVR_mimic:HVR_z"
TERM_4WAY <- "YRS_from_bl:CVR_mimic:HVR_z:SEXMale"

# PREVENT-AD design (Villeneuve et al. 2025)
N_WOMEN <- 279L
N_MEN <- 108L
RETENTION_YEAR1 <- 0.95   # 95% with at least one follow-up
RETENTION_YEAR10 <- 0.79  # 79% still followed beyond 10 years
# Wave w is at year w - 1; geometric decay between the two
# known points, extended at the same rate beyond year 10.
DROPOUT_RATE <- 1 - (RETENTION_YEAR10 / RETENTION_YEAR1)^(1 / 9)

DESIGNS <- list(
  conservative = list(
    label = "5 annual waves, dropout",
    n_waves = 5L, dropout = TRUE
  ),
  primary = list(
    label = "8 annual waves, dropout",
    n_waves = 8L, dropout = TRUE
  ),
  ten = list(
    label = "10 annual waves, dropout",
    n_waves = 10L, dropout = TRUE
  ),
  twelve = list(
    label = "12 annual waves, dropout",
    n_waves = 12L, dropout = TRUE
  ),
  complete = list(
    label = "8 annual waves, complete",
    n_waves = 8L, dropout = FALSE
  )
)
# Design taken to match PREVENT-AD's real follow-up
MATCHED_DESIGN <- "twelve"
WAVE_DESIGNS <- c("conservative", "primary", "ten", "twelve")

SAMPLES <- list(
  women = list(label = "Women only", n_w = N_WOMEN, n_m = 0L),
  men = list(label = "Men only", n_w = 0L, n_m = N_MEN),
  pooled = list(
    label = "Pooled (SEX covariate)", n_w = N_WOMEN, n_m = N_MEN
  )
)

# Manuscript-reported anchor values (verified below)
REPORTED_PRIMARY <- list(beta = -0.0138, se = 0.0054, p = 0.011)
REPORTED_SECONDARY <- list(beta = -0.0226, p = 0.008)

# |beta| grids for minimum detectable effect (0 = type I check)
MDE_GRID <- c(
  0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060
)
MDE_GRID_SEX <- c(
  0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080, 0.100, 0.120
)
# MDE scan only for the scientifically relevant (A-/CU) anchor
MDE_ANCHORS <- "primary"

# Predictor-spread sensitivity (primary anchor). PREVENT-AD is
# younger and nearer the UKB norm, so HVR_z may be less dispersed
# than in ADNI; CVR_mimic may be re-standardised on PREVENT-AD.
# Only SDs matter: the highest-order interaction coefficient and
# its SE are invariant to centring of the constituent variables.
SPREAD_SCENARIOS <- list(
  hvr_sd_1.0 = list(
    label = "HVR_z SD = 1.0 (CVR SD as ADNI)",
    override = list(hvr_sd = 1.0)
  ),
  unit_sd = list(
    label = "HVR_z SD = 1.0, CVR SD = 1.0",
    override = list(hvr_sd = 1.0, cvr_sd = 1.0)
  ),
  tight_hvr = list(
    label = "HVR_z SD = 0.8, CVR SD = 1.0",
    override = list(hvr_sd = 0.8, cvr_sd = 1.0)
  )
)
SPREAD_DESIGNS <- c("primary", "twelve")

LME_OPTIMIZER <- get_script_setting(
  "lme", "optimizer", default = "bobyqa"
)
LME_MAXITER <- get_script_setting(
  "lme", "max_iter", default = 100000
)
LME_CONTROL <- lmerControl(
  optimizer = LME_OPTIMIZER,
  optCtrl = list(maxfun = LME_MAXITER)
)

ANEG_PATH <- file.path(
  get_data_path("models", "lme_results_dir"),
  "lme_amyloid_negative.rds"
)
APOS_PATH <- get_data_path("models", "lme_hvr_z")
OUT_RDS <- here("outputs/power_preventad.rds")
OUT_MD <- here("outputs/power_preventad.md")

rel_path.fn <- function(p) {
  # Repo-relative path for reporting
  sub(paste0(here(), "/"), "", p, fixed = TRUE)
}
ANEG_REL <- rel_path.fn(ANEG_PATH)
APOS_REL <- rel_path.fn(APOS_PATH)

log_info("Replicates per condition: %d", N_REPS)
log_info("Cores: %d", N_CORES)
log_info("Seed: %d", SEED)

# ============================================================
# 1. Load fitted models and extract parameters
# ============================================================
log_section("Loading fitted LME objects")

aneg.res <- read_rds_safe(ANEG_PATH, "A-/CU LME results")
apos.res <- read_rds_safe(APOS_PATH, "A+ LME results")

extract_params.fn <- function(fit, label, source_path,
                              term = TERM_3WAY) {
  # Pull everything the simulation needs from one lmer fit.
  # `term` is the fixed-effect coefficient under test.
  stopifnot(inherits(fit, "merMod"))
  frame.dt <- as.data.table(fit@frame)
  setorder(frame.dt, PTID, YRS_from_bl)
  bl.dt <- frame.dt[, .SD[1], by = PTID]

  vc.dt <- as.data.table(VarCorr(fit))
  get_vc.fn <- function(g, v1, v2 = NA) {
    if (is.na(v2)) {
      vc.dt[grp == g & var1 == v1 & is.na(var2), vcov]
    } else {
      vc.dt[grp == g & var1 == v1 & var2 == v2, vcov]
    }
  }
  var_int <- get_vc.fn("PTID", "(Intercept)")
  var_slp <- get_vc.fn("PTID", "YRS_from_bl")
  cov_is <- get_vc.fn("PTID", "(Intercept)", "YRS_from_bl")
  var_res <- vc.dt[grp == "Residual", vcov]

  coefs.mat <- summary(fit)$coefficients
  stopifnot(term %in% rownames(coefs.mat))
  b3 <- coefs.mat[term, "Estimate"]
  se3 <- coefs.mat[term, "Std. Error"]

  nvis.v <- frame.dt[, .N, by = PTID]$N
  fu.v <- frame.dt[, max(YRS_from_bl), by = PTID]$V1

  has_sex <- "SEX" %in% names(frame.dt)
  apoe.v <- as.integer(bl.dt$APOE4)

  list(
    label = label,
    source_path = source_path,
    term = term,
    formula = deparse1(formula(fit)),
    fixed_formula = lme4::nobars(formula(fit))[-2],
    response = as.character(formula(fit)[[2]]),
    fixef = fixef(fit),
    coef_table = coefs.mat,
    beta3 = b3,
    se3 = se3,
    p3 = coefs.mat[term, "Pr(>|t|)"],
    beta3_ci = c(
      lower = b3 - qnorm(0.975) * se3,
      upper = b3 + qnorm(0.975) * se3
    ),
    var_int = var_int,
    var_slp = var_slp,
    cov_int_slp = cov_is,
    var_res = var_res,
    g.mat = matrix(
      c(var_int, cov_is, cov_is, var_slp), 2, 2
    ),
    is_singular = isSingular(fit),
    is_reml = isREML(fit),
    n_obs = nrow(frame.dt),
    n_subj = nrow(bl.dt),
    has_sex = has_sex,
    apoe_logical = is.logical(frame.dt$APOE4),
    pred = list(
      cvr_mean = mean(bl.dt$CVR_mimic),
      cvr_sd = sd(bl.dt$CVR_mimic),
      hvr_mean = mean(bl.dt$HVR_z),
      hvr_sd = sd(bl.dt$HVR_z),
      cvr_hvr_cor = cor(bl.dt$CVR_mimic, bl.dt$HVR_z),
      age_mean = mean(bl.dt$Age_bl),
      age_sd = sd(bl.dt$Age_bl),
      educ_mean = mean(bl.dt$EDUC),
      educ_sd = sd(bl.dt$EDUC),
      apoe_rate = mean(apoe.v, na.rm = TRUE),
      prop_female = if (has_sex) {
        mean(bl.dt$SEX == "Female")
      } else {
        NA_real_
      }
    ),
    visits = list(
      n_visits_median = median(nvis.v),
      n_visits_q1 = as.numeric(quantile(nvis.v, 0.25)),
      n_visits_q3 = as.numeric(quantile(nvis.v, 0.75)),
      n_visits_max = max(nvis.v),
      n_visits_table = table(nvis.v),
      fu_mean = mean(fu.v),
      fu_median = median(fu.v),
      fu_sd = sd(fu.v),
      fu_max = max(fu.v),
      yrs_quantiles = quantile(
        frame.dt$YRS_from_bl, c(0, .25, .5, .75, 1)
      )
    )
  )
}

# (a) A-/CU pooled-sex models, three domains
aneg_params.lst <- lapply(c("MEM", "LAN", "EXF"), function(d) {
  extract_params.fn(
    aneg.res$results$CVR_mimic[[d]]$model$fit,
    sprintf("A-/CU pooled-sex, %s", d), ANEG_REL
  )
})
names(aneg_params.lst) <- paste0("aneg_", c("MEM", "LAN", "EXF"))

# (b) A+ sex-stratified models, 2 sexes x 3 domains
apos_params.lst <- list()
for (sx in c("Female", "Male")) {
  for (d in c("MEM", "LAN", "EXF")) {
    key <- sprintf("apos_%s_%s", sx, d)
    apos_params.lst[[key]] <- extract_params.fn(
      apos.res$stratified$CVR_mimic[[sx]][[d]]$model$fit,
      sprintf("A+ %s, %s", sx, d), APOS_REL
    )
  }
}
all_params.lst <- c(aneg_params.lst, apos_params.lst)

params_table.dt <- rbindlist(lapply(all_params.lst, function(p) {
  data.table(
    Model = p$label,
    N_subj = p$n_subj, N_obs = p$n_obs,
    Beta3 = p$beta3, SE3 = p$se3, p3 = p$p3,
    Var_int = p$var_int, Var_slope = p$var_slp,
    Cov_int_slope = p$cov_int_slp, Var_resid = p$var_res,
    CVR_mean = p$pred$cvr_mean, CVR_sd = p$pred$cvr_sd,
    HVR_mean = p$pred$hvr_mean, HVR_sd = p$pred$hvr_sd,
    r_CVR_HVR = p$pred$cvr_hvr_cor,
    Visits_median = p$visits$n_visits_median,
    FU_mean_yrs = p$visits$fu_mean,
    FU_max_yrs = p$visits$fu_max
  )
}))
log_info("Extracted parameters from %d fitted models",
         nrow(params_table.dt))
print(params_table.dt[, .(Model, N_subj, N_obs, Beta3, SE3, p3)])

# --- Anchors ---
ANCHORS <- list(
  primary = all_params.lst$aneg_MEM,
  secondary = all_params.lst$apos_Female_EXF
)

# --- Verify anchors against manuscript-reported values ---
check_anchor.fn <- function(p, reported, name) {
  flags.v <- character(0)
  if (round(p$beta3, 4) != reported$beta) {
    flags.v <- c(flags.v, sprintf(
      "beta extracted %.4f vs reported %.4f",
      p$beta3, reported$beta
    ))
  }
  if (!is.null(reported$se) && round(p$se3, 4) != reported$se) {
    flags.v <- c(flags.v, sprintf(
      "SE extracted %.4f vs reported %.4f", p$se3, reported$se
    ))
  }
  if (round(p$p3, 3) != reported$p) {
    flags.v <- c(flags.v, sprintf(
      "p extracted %.3f vs reported %.3f", p$p3, reported$p
    ))
  }
  if (length(flags.v) == 0) {
    log_info("%s anchor matches manuscript: beta=%.4f, SE=%.4f, p=%.3f",
             name, p$beta3, p$se3, p$p3)
  } else {
    log_warn("%s anchor DISCREPANCY: %s", name,
             paste(flags.v, collapse = "; "))
  }
  flags.v
}
anchor_flags.lst <- list(
  primary = check_anchor.fn(
    ANCHORS$primary, REPORTED_PRIMARY, "PRIMARY"
  ),
  secondary = check_anchor.fn(
    ANCHORS$secondary, REPORTED_SECONDARY, "SECONDARY"
  )
)

# ============================================================
# 1b. Sex-interaction refits (anchor for the 4-way contrast)
# ============================================================
# The saved models never estimated a SEX x CVR x HVR x YRS term
# (A-/CU: SEX as covariate only; A+: sex-stratified). To anchor
# the contrast on data rather than assumption, the manuscript
# model is refit here with SEX crossed into both interaction
# blocks, on the exact analysis frames stored in the saved fits.
log_section("Sex-interaction refits for the 4-way contrast")

SEX_FORMULA_TEMPLATE <- paste0(
  "%s ~ YRS_from_bl * CVR_mimic * HVR_z * SEX + ",
  "I(YRS_from_bl^2) * CVR_mimic * HVR_z * SEX + ",
  "Age_bl + EDUC + APOE4 + (YRS_from_bl | PTID)"
)

fit_sex_refit.fn <- function(dat.dt, response, desc) {
  log_info("Refitting with SEX interactions: %s", desc)
  dat.dt <- copy(dat.dt)
  drop.v <- grep("^I\\(", names(dat.dt), value = TRUE)
  if (length(drop.v) > 0) dat.dt[, (drop.v) := NULL]
  dat.dt[, SEX := factor(SEX, levels = c("Female", "Male"))]
  fit <- suppressMessages(suppressWarnings(lmerTest::lmer(
    as.formula(sprintf(SEX_FORMULA_TEMPLATE, response)),
    data = dat.dt, REML = TRUE, control = LME_CONTROL
  )))
  if (isSingular(fit)) {
    log_warn("  Sex-interaction refit is singular: %s", desc)
  }
  ct <- summary(fit)$coefficients
  log_info("  %s: delta = %.4f (SE %.4f), p = %.4f", TERM_4WAY,
           ct[TERM_4WAY, "Estimate"], ct[TERM_4WAY, "Std. Error"],
           ct[TERM_4WAY, "Pr(>|t|)"])
  fit
}

# A-/CU Memory: frame already contains SEX
aneg_sex.fit <- fit_sex_refit.fn(
  as.data.table(aneg.res$results$CVR_mimic$MEM$model$fit@frame),
  "MEM", "A-/CU pooled, Memory"
)
# A+ EXF: stack the female and male analysis frames
apos_sex.dt <- rbind(
  as.data.table(
    apos.res$stratified$CVR_mimic$Female$EXF$model$fit@frame
  )[, SEX := "Female"],
  as.data.table(
    apos.res$stratified$CVR_mimic$Male$EXF$model$fit@frame
  )[, SEX := "Male"]
)
apos_sex.fit <- fit_sex_refit.fn(
  apos_sex.dt, "EXF", "A+ both sexes, Executive Function"
)

CONTRASTS <- list(
  primary = extract_params.fn(
    aneg_sex.fit, "A-/CU pooled, MEM, SEX-interaction refit",
    ANEG_REL, term = TERM_4WAY
  ),
  upper = extract_params.fn(
    apos_sex.fit, "A+ both sexes, EXF, SEX-interaction refit",
    APOS_REL, term = TERM_4WAY
  )
)
# Sex-specific 3-way effects implied by each refit
implied_sex_effects.fn <- function(p) {
  b <- p$fixef
  c(women = unname(b[TERM_3WAY]),
    men = unname(b[TERM_3WAY] + b[TERM_4WAY]))
}
for (cn in names(CONTRASTS)) {
  e.v <- implied_sex_effects.fn(CONTRASTS[[cn]])
  log_info("%s refit: 3-way women %.4f, men %.4f, delta %.4f",
           cn, e.v["women"], e.v["men"], CONTRASTS[[cn]]$beta3)
}

# ============================================================
# 2. Simulation design
# ============================================================
log_section("Simulation design")

retention.fn <- function(n_waves, dropout) {
  if (!dropout) return(rep(1, n_waves))
  w.v <- seq_len(n_waves)
  ret.v <- c(1, RETENTION_YEAR1 *
    (1 - DROPOUT_RATE)^(pmax(w.v[-1] - 2, 0)))
  ret.v[seq_len(n_waves)]
}
for (dn in names(DESIGNS)) {
  ret.v <- retention.fn(DESIGNS[[dn]]$n_waves, DESIGNS[[dn]]$dropout)
  DESIGNS[[dn]]$retention <- ret.v
  log_info("%s: retention by wave = %s", DESIGNS[[dn]]$label,
           paste(sprintf("%.3f", ret.v), collapse = " "))
}
log_info("Per-wave attrition after year 1: %.2f%%",
         100 * DROPOUT_RATE)

# Draw each subject's last observed wave from the retention curve
draw_last_wave.fn <- function(n, ret.v) {
  n_waves <- length(ret.v)
  p_last.v <- ret.v - c(ret.v[-1], 0)
  sample.int(n_waves, n, replace = TRUE, prob = p_last.v)
}

sim_data.fn <- function(anchor, design, n_w, n_m, beta_target,
                        term = TERM_3WAY, pred_override = NULL) {
  # Simulate one PREVENT-AD-like dataset from the anchor model.
  # `term` is the coefficient set to beta_target; pred_override
  # is a named list replacing entries of anchor$pred.
  n <- n_w + n_m
  pr <- anchor$pred
  if (!is.null(pred_override)) pr <- modifyList(pr, pred_override)
  # Person-level predictors
  sex.fct <- factor(
    c(rep("Female", n_w), rep("Male", n_m)),
    levels = c("Female", "Male")
  )
  z.mat <- MASS::mvrnorm(
    n, mu = c(0, 0),
    Sigma = matrix(c(1, pr$cvr_hvr_cor, pr$cvr_hvr_cor, 1), 2)
  )
  cvr.v <- pr$cvr_mean + pr$cvr_sd * z.mat[, 1]
  hvr.v <- pr$hvr_mean + pr$hvr_sd * z.mat[, 2]
  age.v <- rnorm(n, pr$age_mean, pr$age_sd)
  educ.v <- rnorm(n, pr$educ_mean, pr$educ_sd)
  apoe.v <- rbinom(n, 1, pr$apoe_rate)
  if (anchor$apoe_logical) apoe.v <- as.logical(apoe.v)
  # Random effects
  re.mat <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = anchor$g.mat)
  # Visit structure
  last.v <- draw_last_wave.fn(n, design$retention)
  id.v <- rep(seq_len(n), last.v)
  yrs.v <- unlist(lapply(last.v, function(k) seq_len(k) - 1))
  dat.dt <- data.table(
    PTID = factor(id.v),
    YRS_from_bl = as.numeric(yrs.v),
    CVR_mimic = cvr.v[id.v],
    HVR_z = hvr.v[id.v],
    Age_bl = age.v[id.v],
    EDUC = educ.v[id.v],
    APOE4 = apoe.v[id.v],
    SEX = sex.fct[id.v]
  )
  # Fixed part via the anchor model's own fixed-effect formula
  x.mat <- model.matrix(anchor$fixed_formula, data = dat.dt)
  beta.v <- anchor$fixef
  stopifnot(identical(colnames(x.mat), names(beta.v)))
  beta.v[term] <- beta_target
  y.v <- as.vector(x.mat %*% beta.v) +
    re.mat[id.v, 1] + re.mat[id.v, 2] * dat.dt$YRS_from_bl +
    rnorm(nrow(dat.dt), 0, sqrt(anchor$var_res))
  dat.dt[, Y := y.v]
  dat.dt
}

fit_formula.fn <- function(pooled, sex_interaction = FALSE) {
  if (sex_interaction) {
    return(as.formula(sprintf(SEX_FORMULA_TEMPLATE, "Y")))
  }
  cov_str <- if (pooled) {
    "Age_bl + EDUC + APOE4 + SEX"
  } else {
    "Age_bl + EDUC + APOE4"
  }
  as.formula(paste0(
    "Y ~ YRS_from_bl * CVR_mimic * HVR_z + ",
    "I(YRS_from_bl^2) * CVR_mimic * HVR_z + ",
    cov_str, " + (YRS_from_bl | PTID)"
  ))
}

fit_rep.fn <- function(dat.dt, formula, term) {
  # Same estimator and test as the manuscript: REML, bobyqa,
  # Satterthwaite p-value for the tested term.
  out.lst <- list(
    est = NA_real_, se = NA_real_, p = NA_real_,
    singular = NA, error = FALSE
  )
  tryCatch({
    fit <- suppressMessages(suppressWarnings(
      lmerTest::lmer(
        formula, data = dat.dt, REML = TRUE,
        control = LME_CONTROL
      )
    ))
    coefs.mat <- suppressMessages(summary(fit)$coefficients)
    out.lst$est <- coefs.mat[term, "Estimate"]
    out.lst$se <- coefs.mat[term, "Std. Error"]
    out.lst$p <- coefs.mat[term, "Pr(>|t|)"]
    out.lst$singular <- isSingular(fit)
  }, error = function(e) {
    out.lst$error <<- TRUE
  })
  out.lst
}

# ============================================================
# 3. Power engine
# ============================================================
condition_counter <- 0L

run_condition.fn <- function(anchor_name, design_name, sample_name,
                             beta_target, n_reps, tag,
                             scenario = "adni",
                             pred_override = NULL,
                             contrast = FALSE) {
  # Run n_reps replicates for one condition. contrast = TRUE
  # switches to the sex-by-moderation (4-way) test: anchor from
  # CONTRASTS, SEX-interaction fit formula, pooled sample only.
  # Per-replicate seeds make results independent of core count.
  condition_counter <<- condition_counter + 1L
  cond_id <- condition_counter
  smp <- SAMPLES[[sample_name]]
  design <- DESIGNS[[design_name]]
  if (contrast) {
    stopifnot(sample_name == "pooled")
    anchor <- CONTRASTS[[anchor_name]]
    term <- TERM_4WAY
    formula <- fit_formula.fn(pooled = TRUE, sex_interaction = TRUE)
  } else {
    anchor <- ANCHORS[[anchor_name]]
    term <- TERM_3WAY
    formula <- fit_formula.fn(pooled = smp$n_w > 0 && smp$n_m > 0)
  }

  reps.lst <- mclapply(seq_len(n_reps), function(r) {
    set.seed(SEED + cond_id * 100000L + r)
    dat.dt <- sim_data.fn(
      anchor, design, smp$n_w, smp$n_m, beta_target,
      term = term, pred_override = pred_override
    )
    res.lst <- fit_rep.fn(dat.dt, formula, term)
    data.table(
      rep = r, est = res.lst$est, se = res.lst$se,
      p = res.lst$p, singular = res.lst$singular,
      error = res.lst$error
    )
  }, mc.cores = N_CORES)
  reps.dt <- rbindlist(reps.lst)
  reps.dt[, `:=`(
    tag = tag, test = if (contrast) "sex_contrast" else "moderation",
    anchor = anchor_name, design = design_name,
    sample = sample_name, scenario = scenario, beta3 = beta_target
  )]

  ok.dt <- reps.dt[!is.na(p)]
  pw <- mean(ok.dt$p < ALPHA)
  summary.dt <- data.table(
    tag = tag, test = if (contrast) "sex_contrast" else "moderation",
    anchor = anchor_name, design = design_name,
    sample = sample_name, scenario = scenario,
    n_subj = smp$n_w + smp$n_m, beta3 = beta_target,
    n_reps = n_reps, n_valid = nrow(ok.dt),
    power = pw,
    power_mc_se = sqrt(pw * (1 - pw) / nrow(ok.dt)),
    mean_est = mean(ok.dt$est),
    mean_se = mean(ok.dt$se),
    singular_rate = mean(ok.dt$singular),
    error_rate = mean(reps.dt$error)
  )
  log_info(
    "[%s] %s | %s | %s | %s | beta=%.4f -> power %.3f (sing %.2f)",
    tag, anchor_name, design_name, sample_name, scenario,
    beta_target, pw, summary.dt$singular_rate
  )
  list(summary = summary.dt, reps = reps.dt)
}

interp_mde.fn <- function(abs_beta.v, power.v, target) {
  # Smallest |beta| reaching target power, by linear
  # interpolation on a monotone (cummax) power curve.
  o <- order(abs_beta.v)
  x.v <- abs_beta.v[o]
  y.v <- cummax(power.v[o])
  if (max(y.v) < target) return(NA_real_)
  i <- which(y.v >= target)[1]
  if (i == 1) return(x.v[1])
  x.v[i - 1] + (target - y.v[i - 1]) /
    (y.v[i] - y.v[i - 1]) * (x.v[i] - x.v[i - 1])
}

# ============================================================
# 4. Power at anchor effect sizes (+ both CI bounds of primary)
# ============================================================
log_section("Power at anchor effect sizes")

beta_specs.lst <- list(
  list(anchor = "primary", tag = "anchor",
       beta3 = ANCHORS$primary$beta3),
  list(anchor = "primary", tag = "ci_lower",
       beta3 = ANCHORS$primary$beta3_ci[["lower"]]),
  list(anchor = "primary", tag = "ci_upper",
       beta3 = ANCHORS$primary$beta3_ci[["upper"]]),
  list(anchor = "secondary", tag = "anchor",
       beta3 = ANCHORS$secondary$beta3)
)

main_runs.lst <- list()
for (spec in beta_specs.lst) {
  for (dn in names(DESIGNS)) {
    for (sn in names(SAMPLES)) {
      key <- paste(spec$anchor, spec$tag, dn, sn, sep = "|")
      main_runs.lst[[key]] <- run_condition.fn(
        spec$anchor, dn, sn, spec$beta3, N_REPS, spec$tag
      )
    }
  }
}
power_main.dt <- rbindlist(lapply(main_runs.lst, `[[`, "summary"))
reps_main.dt <- rbindlist(lapply(main_runs.lst, `[[`, "reps"))

# ============================================================
# 5. Minimum detectable effect (grid scan, 3-way)
# ============================================================
log_section("Minimum detectable effect: grid scan")

mde_runs.lst <- list()
for (an in MDE_ANCHORS) {
  for (dn in names(DESIGNS)) {
    for (sn in names(SAMPLES)) {
      for (b in MDE_GRID) {
        key <- paste(an, dn, sn, b, sep = "|")
        mde_runs.lst[[key]] <- run_condition.fn(
          an, dn, sn, -b, N_REPS, "mde_grid"
        )
      }
    }
  }
}
power_grid.dt <- rbindlist(lapply(mde_runs.lst, `[[`, "summary"))
reps_grid.dt <- rbindlist(lapply(mde_runs.lst, `[[`, "reps"))
power_grid.dt[, abs_beta := abs(beta3)]

mde.dt <- power_grid.dt[
  , .(
    mde_80 = interp_mde.fn(abs_beta, power, TARGET_POWER),
    type1_at_zero = power[abs_beta == 0],
    max_power_on_grid = max(power)
  ),
  by = .(anchor, design, sample, n_subj)
]
mde.dt[, mde_vs_anchor := mde_80 /
  sapply(anchor, function(a) abs(ANCHORS[[a]]$beta3))]

# ============================================================
# 5b. Predictor-spread sensitivity (primary anchor)
# ============================================================
log_section("Predictor-spread sensitivity")

spread_runs.lst <- list()
for (scn in names(SPREAD_SCENARIOS)) {
  for (dn in SPREAD_DESIGNS) {
    for (sn in names(SAMPLES)) {
      key <- paste(scn, dn, sn, sep = "|")
      spread_runs.lst[[key]] <- run_condition.fn(
        "primary", dn, sn, ANCHORS$primary$beta3, N_REPS,
        "spread", scenario = scn,
        pred_override = SPREAD_SCENARIOS[[scn]]$override
      )
    }
  }
}
power_spread.dt <- rbindlist(lapply(spread_runs.lst, `[[`, "summary"))
reps_spread.dt <- rbindlist(lapply(spread_runs.lst, `[[`, "reps"))
# Add the ADNI-spread reference rows from the main run
power_spread.dt <- rbind(
  power_main.dt[
    anchor == "primary" & tag == "anchor" & design %in% SPREAD_DESIGNS
  ],
  power_spread.dt
)

# ============================================================
# 5c. Sex-by-moderation contrast (4-way), pooled sample
# ============================================================
log_section("Sex-by-moderation contrast (4-way)")

sex_runs.lst <- list()
for (cn in names(CONTRASTS)) {
  for (dn in names(DESIGNS)) {
    key <- paste(cn, "anchor", dn, sep = "|")
    sex_runs.lst[[key]] <- run_condition.fn(
      cn, dn, "pooled", CONTRASTS[[cn]]$beta3, N_REPS,
      "sex_anchor", contrast = TRUE
    )
  }
}
for (dn in names(DESIGNS)) {
  for (b in MDE_GRID_SEX) {
    key <- paste("primary", "grid", dn, b, sep = "|")
    sex_runs.lst[[key]] <- run_condition.fn(
      "primary", dn, "pooled", sign(CONTRASTS$primary$beta3) * b,
      N_REPS, "sex_grid", contrast = TRUE
    )
  }
}
power_sex.dt <- rbindlist(lapply(sex_runs.lst, `[[`, "summary"))
reps_sex.dt <- rbindlist(lapply(sex_runs.lst, `[[`, "reps"))
power_sex.dt[, abs_beta := abs(beta3)]

mde_sex.dt <- power_sex.dt[
  tag == "sex_grid",
  .(
    mde_80 = interp_mde.fn(abs_beta, power, TARGET_POWER),
    type1_at_zero = power[abs_beta == 0],
    max_power_on_grid = max(power)
  ),
  by = .(anchor, design, sample, n_subj)
]
mde_sex.dt[, mde_vs_anchor := mde_80 / abs(CONTRASTS$primary$beta3)]

# ============================================================
# 6. Assemble results and report
# ============================================================
log_section("Assembling outputs")

label_design.fn <- function(x) {
  sapply(x, function(d) DESIGNS[[d]]$label)
}
label_sample.fn <- function(x) {
  sapply(x, function(s) SAMPLES[[s]]$label)
}
label_anchor.fn <- function(x) {
  sapply(x, function(a) ANCHORS[[a]]$label)
}
label_contrast.fn <- function(x) {
  sapply(x, function(a) CONTRASTS[[a]]$label)
}
label_scenario.fn <- function(x) {
  sapply(x, function(s) {
    if (s == "adni") {
      sprintf("ADNI spread (HVR_z SD %.2f, CVR SD %.2f)",
              ANCHORS$primary$pred$hvr_sd,
              ANCHORS$primary$pred$cvr_sd)
    } else {
      SPREAD_SCENARIOS[[s]]$label
    }
  })
}
pct.fn <- function(x) sprintf("%.1f%%", 100 * x)

# Main table: power by design x anchor x sample
main_table.dt <- power_main.dt[tag == "anchor", .(
  Anchor = label_anchor.fn(anchor),
  Design = label_design.fn(design),
  Sample = label_sample.fn(sample),
  N = n_subj,
  Beta3 = sprintf("%.4f", beta3),
  Power = pct.fn(power),
  MC_SE = pct.fn(power_mc_se),
  Singular = pct.fn(singular_rate),
  Valid = sprintf("%d/%d", n_valid, n_reps)
)]

# CI-bound table (primary anchor, both bounds)
ci_table.dt <- dcast(
  power_main.dt[anchor == "primary"],
  design + sample + n_subj ~ tag, value.var = "power"
)
setcolorder(ci_table.dt, c(
  "design", "sample", "n_subj", "ci_lower", "anchor", "ci_upper"
))
ci_table_fmt.dt <- ci_table.dt[, .(
  Design = label_design.fn(design),
  Sample = label_sample.fn(sample),
  N = n_subj,
  Power_at_CI_lower = sprintf(
    "%s (beta %.4f)", pct.fn(ci_lower),
    ANCHORS$primary$beta3_ci[["lower"]]
  ),
  Power_at_anchor = sprintf(
    "%s (beta %.4f)", pct.fn(anchor), ANCHORS$primary$beta3
  ),
  Power_at_CI_upper = sprintf(
    "%s (beta %.4f)", pct.fn(ci_upper),
    ANCHORS$primary$beta3_ci[["upper"]]
  )
)]

fmt_mde.fn <- function(m, grid.v) {
  ifelse(is.na(m), sprintf("> %.3f", max(grid.v)), sprintf("%.4f", m))
}
mde_table.dt <- mde.dt[, .(
  Anchor = label_anchor.fn(anchor),
  Design = label_design.fn(design),
  Sample = label_sample.fn(sample),
  N = n_subj,
  MDE_80 = fmt_mde.fn(mde_80, MDE_GRID),
  MDE_over_anchor = ifelse(
    is.na(mde_vs_anchor), "NA", sprintf("%.2fx", mde_vs_anchor)
  ),
  Type1_at_beta0 = pct.fn(type1_at_zero)
)]

spread_table.dt <- power_spread.dt[, .(
  Scenario = label_scenario.fn(scenario),
  Design = label_design.fn(design),
  Sample = label_sample.fn(sample),
  N = n_subj,
  Power = pct.fn(power),
  MC_SE = pct.fn(power_mc_se),
  Singular = pct.fn(singular_rate)
)]

sex_table.dt <- power_sex.dt[tag == "sex_anchor", .(
  Contrast_anchor = label_contrast.fn(anchor),
  Design = label_design.fn(design),
  N = n_subj,
  Delta = sprintf("%.4f", beta3),
  Power = pct.fn(power),
  MC_SE = pct.fn(power_mc_se),
  Singular = pct.fn(singular_rate),
  Valid = sprintf("%d/%d", n_valid, n_reps)
)]
sex_mde_table.dt <- mde_sex.dt[, .(
  Design = label_design.fn(design),
  N = n_subj,
  MDE_80 = fmt_mde.fn(mde_80, MDE_GRID_SEX),
  MDE_over_observed_delta = ifelse(
    is.na(mde_vs_anchor), "NA", sprintf("%.2fx", mde_vs_anchor)
  ),
  Type1_at_delta0 = pct.fn(type1_at_zero)
)]

# --- Headline numbers ---
get_power.fn <- function(an, dn, sn, tg = "anchor") {
  power_main.dt[
    anchor == an & design == dn & sample == sn & tag == tg, power
  ]
}
matched <- DESIGNS[[MATCHED_DESIGN]]
mde_headline <- mde.dt[
  anchor == "primary" & design == MATCHED_DESIGN & sample == "pooled",
  mde_80
]
waves_power.dt <- power_main.dt[
  anchor == "primary" & tag == "anchor" & sample == "pooled" &
    design %in% WAVE_DESIGNS
]
waves_power.dt[, n_waves := sapply(design, function(d) {
  DESIGNS[[d]]$n_waves
})]
setorder(waves_power.dt, n_waves)
waves_mde.dt <- mde.dt[
  anchor == "primary" & sample == "pooled" & design %in% WAVE_DESIGNS
]
waves_mde.dt[, n_waves := sapply(design, function(d) {
  DESIGNS[[d]]$n_waves
})]
setorder(waves_mde.dt, n_waves)

headline_mde <- sprintf(
  paste0(
    "Under the design matching PREVENT-AD's follow-up (%s; years ",
    "0-%d; retention %.2f at year 1 and %.2f at year 10, %.2f at ",
    "the last wave), the pooled sample of N = %d has 80%% power ",
    "for a CVR_mimic x HVR_z x YRS effect of |beta| >= %s, i.e. ",
    "%.2f times the observed A-/CU effect (|beta| = %.4f)."
  ),
  matched$label, matched$n_waves - 1, RETENTION_YEAR1,
  RETENTION_YEAR10, matched$retention[matched$n_waves],
  N_WOMEN + N_MEN,
  fmt_mde.fn(mde_headline, MDE_GRID),
  mde_headline / abs(ANCHORS$primary$beta3),
  abs(ANCHORS$primary$beta3)
)
headline_waves <- paste0(
  "Pooled power at the observed A-/CU effect by follow-up length: ",
  paste(sprintf("%d waves = %.0f%%", waves_power.dt$n_waves,
                100 * waves_power.dt$power), collapse = ", "),
  "; MDE at 80%: ",
  paste(sprintf("%d waves = %s", waves_mde.dt$n_waves,
                fmt_mde.fn(waves_mde.dt$mde_80, MDE_GRID)),
        collapse = ", "),
  "."
)

p_pooled <- get_power.fn("primary", MATCHED_DESIGN, "pooled")
p_women <- get_power.fn("primary", MATCHED_DESIGN, "women")
p_men <- get_power.fn("primary", MATCHED_DESIGN, "men")
final_sentence <- sprintf(
  paste0(
    "Simulation gives %.0f%% power to detect the observed ",
    "moderation in %d women and %d men."
  ),
  100 * p_pooled, N_WOMEN, N_MEN
)
sex_sentence <- sprintf(
  paste0(
    "Analysed separately under the same design, power is %.0f%% ",
    "in the %d women and %.0f%% in the %d men."
  ),
  100 * p_women, N_WOMEN, 100 * p_men, N_MEN
)
ci_sentence <- sprintf(
  paste0(
    "At the two 95%% CI bounds of the anchor (beta = %.4f and ",
    "%.4f), pooled power under the matched design is %.0f%% and ",
    "%.0f%%."
  ),
  ANCHORS$primary$beta3_ci[["lower"]],
  ANCHORS$primary$beta3_ci[["upper"]],
  100 * get_power.fn("primary", MATCHED_DESIGN, "pooled", "ci_lower"),
  100 * get_power.fn("primary", MATCHED_DESIGN, "pooled", "ci_upper")
)
sex_pw <- power_sex.dt[
  tag == "sex_anchor" & anchor == "primary" & design == MATCHED_DESIGN,
  power
]
sex_mde_matched <- mde_sex.dt[design == MATCHED_DESIGN, mde_80]
sex_e.v <- implied_sex_effects.fn(CONTRASTS$primary)
contrast_sentence <- sprintf(
  paste0(
    "Sex-by-moderation contrast (SEX x CVR_mimic x HVR_z x YRS, ",
    "pooled N = %d, matched design): the A-/CU refit gives delta ",
    "= %.4f (SE %.4f, p = %.3f; 3-way %.4f in women vs %.4f in ",
    "men). Power to detect that delta is %.0f%%; the minimum ",
    "detectable delta at 80%% power is %s."
  ),
  N_WOMEN + N_MEN, CONTRASTS$primary$beta3, CONTRASTS$primary$se3,
  CONTRASTS$primary$p3, sex_e.v["women"], sex_e.v["men"],
  100 * sex_pw, fmt_mde.fn(sex_mde_matched, MDE_GRID_SEX)
)
p2_pooled <- get_power.fn("secondary", MATCHED_DESIGN, "pooled")
p2_women <- get_power.fn("secondary", MATCHED_DESIGN, "women")
p2_men <- get_power.fn("secondary", MATCHED_DESIGN, "men")
sex_upper_pw <- power_sex.dt[
  tag == "sex_anchor" & anchor == "upper" & design == MATCHED_DESIGN,
  power
]
secondary_sentence <- sprintf(
  paste0(
    "Upper bound: under the A+ female Executive Function effect ",
    "(beta = %.4f, amyloid-positive participants at a later ",
    "disease stage), matched-design power is %.0f%% pooled, ",
    "%.0f%% in women and %.0f%% in men; the A+ sex contrast ",
    "(delta = %.4f) has %.0f%% power."
  ),
  ANCHORS$secondary$beta3, 100 * p2_pooled, 100 * p2_women,
  100 * p2_men, CONTRASTS$upper$beta3, 100 * sex_upper_pw
)

results.lst <- list(
  meta = list(
    script = "16_power_preventad.R",
    timestamp = Sys.time(),
    seed = SEED,
    n_reps = N_REPS,
    n_cores = N_CORES,
    alpha = ALPHA,
    test = paste(
      "lmerTest Satterthwaite t-test (REML, bobyqa) on",
      TERM_3WAY, "or", TERM_4WAY,
      "- same estimator/test as manuscript scripts 08/10b"
    ),
    simr_installed = requireNamespace("simr", quietly = TRUE),
    lme4_version = as.character(packageVersion("lme4")),
    lmerTest_version = as.character(packageVersion("lmerTest")),
    r_version = R.version.string
  ),
  sources = list(aneg = ANEG_REL, apos = APOS_REL),
  design = list(
    n_women = N_WOMEN, n_men = N_MEN,
    retention_year1 = RETENTION_YEAR1,
    retention_year10 = RETENTION_YEAR10,
    dropout_rate_per_wave = DROPOUT_RATE,
    matched_design = MATCHED_DESIGN,
    designs = DESIGNS, samples = SAMPLES,
    reference = paste(
      "Villeneuve et al. 2025, Alzheimers Dement,",
      "doi 10.1002/alz.70653"
    )
  ),
  anchors = ANCHORS,
  contrasts = CONTRASTS,
  sex_refits = list(
    aneg_coef_table = summary(aneg_sex.fit)$coefficients,
    apos_coef_table = summary(apos_sex.fit)$coefficients,
    formula_template = SEX_FORMULA_TEMPLATE
  ),
  anchor_flags = anchor_flags.lst,
  reported_values = list(
    primary = REPORTED_PRIMARY, secondary = REPORTED_SECONDARY
  ),
  model_params = all_params.lst,
  params_table = params_table.dt,
  power_main = power_main.dt,
  power_grid = power_grid.dt,
  power_spread = power_spread.dt,
  power_sex = power_sex.dt,
  spread_scenarios = SPREAD_SCENARIOS,
  mde = mde.dt,
  mde_sex = mde_sex.dt,
  mde_grid = MDE_GRID,
  mde_grid_sex = MDE_GRID_SEX,
  reps_main = reps_main.dt,
  reps_grid = reps_grid.dt,
  reps_spread = reps_spread.dt,
  reps_sex = reps_sex.dt,
  tables = list(
    main = main_table.dt, ci = ci_table_fmt.dt,
    mde = mde_table.dt, spread = spread_table.dt,
    sex = sex_table.dt, sex_mde = sex_mde_table.dt
  ),
  sentences = list(
    headline_mde = headline_mde, headline_waves = headline_waves,
    final = final_sentence, by_sex = sex_sentence,
    ci_bounds = ci_sentence, sex_contrast = contrast_sentence,
    upper_bound = secondary_sentence
  )
)
ensure_directory(dirname(OUT_RDS))
write_rds_safe(results.lst, OUT_RDS, "PREVENT-AD power results")

# ---------------------------------------------------------
# Markdown report
# ---------------------------------------------------------
md_table.fn <- function(dt) {
  cols.v <- names(dt)
  hdr <- paste0("| ", paste(cols.v, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(cols.v)),
                           collapse = "|"), "|")
  rows.v <- apply(dt, 1, function(r) {
    paste0("| ", paste(r, collapse = " | "), " |")
  })
  paste(c(hdr, sep, rows.v), collapse = "\n")
}

fmt_params.fn <- function(p) {
  pr <- p$pred
  vs <- p$visits
  paste0(
    "- **", p$label, "** (`", p$source_path, "`)\n",
    sprintf("  - Formula: `%s`\n", p$formula),
    sprintf("  - N = %d subjects, %d observations; REML = %s\n",
            p$n_subj, p$n_obs, p$is_reml),
    sprintf(paste0("  - Tested term (%s): beta = %.5f, SE = %.5f, ",
                   "p = %.4f; Wald 95%% CI [%.4f, %.4f]\n"),
            p$term, p$beta3, p$se3, p$p3,
            p$beta3_ci[["lower"]], p$beta3_ci[["upper"]]),
    sprintf(paste0("  - Random effects: var(int) = %.4f, ",
                   "var(slope) = %.5f, cov = %.5f ",
                   "(r = %.3f); residual var = %.4f\n"),
            p$var_int, p$var_slp, p$cov_int_slp,
            p$cov_int_slp / sqrt(p$var_int * p$var_slp),
            p$var_res),
    sprintf(paste0("  - Baseline predictors: CVR_mimic M = %.3f ",
                   "(SD %.3f); HVR_z M = %.3f (SD %.3f); ",
                   "r(CVR, HVR) = %.3f\n"),
            pr$cvr_mean, pr$cvr_sd, pr$hvr_mean, pr$hvr_sd,
            pr$cvr_hvr_cor),
    sprintf(paste0("  - Covariates: Age_bl M = %.1f (SD %.1f); ",
                   "EDUC M = %.1f (SD %.1f); APOE4 rate = %.2f",
                   "%s\n"),
            pr$age_mean, pr$age_sd, pr$educ_mean, pr$educ_sd,
            pr$apoe_rate,
            if (is.na(pr$prop_female)) "" else
              sprintf("; %% female = %.2f", pr$prop_female)),
    sprintf(paste0("  - Visits: median %g per subject (IQR %g-%g, ",
                   "max %d); follow-up mean %.1f y (median %.1f, ",
                   "SD %.1f, max %.1f); YRS quartiles %s\n"),
            vs$n_visits_median, vs$n_visits_q1, vs$n_visits_q3,
            vs$n_visits_max, vs$fu_mean, vs$fu_median, vs$fu_sd,
            vs$fu_max,
            paste(sprintf("%.2f", vs$yrs_quantiles),
                  collapse = " / "))
  )
}

fixef_table.fn <- function(p) {
  ct <- p$coef_table
  data.table(
    Term = rownames(ct),
    Estimate = sprintf("%.5f", ct[, "Estimate"]),
    SE = sprintf("%.5f", ct[, "Std. Error"]),
    p = sprintf("%.4f", ct[, "Pr(>|t|)"])
  )
}

flag_text.fn <- function(flags.v) {
  if (length(flags.v) == 0) {
    "matches the manuscript-reported values (to rounding)."
  } else {
    paste0("DISCREPANCY: ", paste(flags.v, collapse = "; "), ".")
  }
}

other_params.dt <- params_table.dt[, .(
  Model,
  N_subj, N_obs,
  Beta3 = sprintf("%.4f", Beta3),
  SE3 = sprintf("%.4f", SE3),
  p3 = sprintf("%.4f", p3),
  Var_int = sprintf("%.4f", Var_int),
  Var_slope = sprintf("%.5f", Var_slope),
  Cov = sprintf("%.5f", Cov_int_slope),
  Var_resid = sprintf("%.4f", Var_resid),
  CVR_sd = sprintf("%.3f", CVR_sd),
  HVR_mean = sprintf("%.3f", HVR_mean),
  HVR_sd = sprintf("%.3f", HVR_sd),
  r = sprintf("%.3f", r_CVR_HVR),
  Visits_med = Visits_median,
  FU_mean = sprintf("%.1f", FU_mean_yrs)
)]

retention_lines.v <- sapply(names(DESIGNS), function(d) {
  sprintf("  - %s: %s", DESIGNS[[d]]$label, paste(
    sprintf("%.2f", DESIGNS[[d]]$retention), collapse = " / "
  ))
})

md.v <- c(
  "# Monte Carlo power analysis: CVR_mimic x HVR_z x YRS in a",
  "# PREVENT-AD-like cohort",
  "",
  sprintf("Generated by `R/scripts/16_power_preventad.R` on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf(paste0("Seed %d; %d replicates per condition; alpha = ",
                 "%.2f; %s."),
          SEED, N_REPS, ALPHA, results.lst$meta$test),
  sprintf(paste0("`simr` installed in renv: %s (simulation ",
                 "hand-rolled with lme4 %s / lmerTest %s)."),
          results.lst$meta$simr_installed,
          results.lst$meta$lme4_version,
          results.lst$meta$lmerTest_version),
  "",
  "## Headline",
  "",
  paste0("**", headline_mde, "**"),
  "",
  paste0("**", headline_waves, "**"),
  "",
  final_sentence,
  "",
  sex_sentence,
  "",
  ci_sentence,
  "",
  contrast_sentence,
  "",
  secondary_sentence,
  "",
  paste0("The scientifically relevant anchor is the A-/CU model: ",
         "PREVENT-AD is younger, cognitively unimpaired and ",
         "earlier in the disease timeline, so the amyloid-negative ",
         "moderation is the effect it would test. The A+ anchor ",
         "is retained only as an upper bound."),
  "",
  "## 1. Parameter sources",
  "",
  "Fitted `lmerModLmerTest` objects (not just coefficient tables)",
  "were available for all nine models, so all fixed effects,",
  "random-effect (co)variances and residual variances were read",
  "directly from the fits. Baseline predictor distributions and",
  "visit structure were read from each fit's model frame.",
  "",
  sprintf("- A-/CU pooled-sex models (script 10b): `%s`", ANEG_REL),
  sprintf("- A+ sex-stratified models (script 08): `%s`", APOS_REL),
  "",
  "### Anchor models (used for data generation)",
  "",
  fmt_params.fn(ANCHORS$primary),
  sprintf("  - Verification vs manuscript (beta %.4f, SE %.4f, p %.3f): %s",
          REPORTED_PRIMARY$beta, REPORTED_PRIMARY$se,
          REPORTED_PRIMARY$p, flag_text.fn(anchor_flags.lst$primary)),
  "",
  fmt_params.fn(ANCHORS$secondary),
  sprintf("  - Verification vs manuscript (beta %.4f, p %.3f): %s",
          REPORTED_SECONDARY$beta, REPORTED_SECONDARY$p,
          flag_text.fn(anchor_flags.lst$secondary)),
  "",
  "### Sex-interaction refits (anchor for the 4-way contrast)",
  "",
  paste0("No saved model estimated a SEX x CVR_mimic x HVR_z x YRS ",
         "term (the A-/CU model used SEX as a covariate; the A+ ",
         "models were sex-stratified). The manuscript model was ",
         "therefore refit in this script with SEX crossed into ",
         "both interaction blocks, on the exact analysis frames ",
         "stored in the saved fits. These refits are new estimates, ",
         "not manuscript results."),
  "",
  sprintf("Formula: `%s`", sprintf(SEX_FORMULA_TEMPLATE, "Cog")),
  "",
  fmt_params.fn(CONTRASTS$primary),
  sprintf("  - Implied 3-way effect: women %.4f, men %.4f",
          sex_e.v["women"], sex_e.v["men"]),
  "",
  fmt_params.fn(CONTRASTS$upper),
  sprintf("  - Implied 3-way effect: women %.4f, men %.4f",
          implied_sex_effects.fn(CONTRASTS$upper)["women"],
          implied_sex_effects.fn(CONTRASTS$upper)["men"]),
  "",
  "#### Full fixed effects, primary anchor (A-/CU, Memory)",
  "",
  md_table.fn(fixef_table.fn(ANCHORS$primary)),
  "",
  "#### Full fixed effects, upper-bound anchor (A+ females, EXF)",
  "",
  md_table.fn(fixef_table.fn(ANCHORS$secondary)),
  "",
  "#### Full fixed effects, A-/CU SEX-interaction refit",
  "",
  md_table.fn(fixef_table.fn(CONTRASTS$primary)),
  "",
  "### All nine CVR_mimic 3-way models (for reference)",
  "",
  md_table.fn(other_params.dt),
  "",
  "## 2. Simulation design (PREVENT-AD)",
  "",
  sprintf(paste0("- Target cohort: N = %d (%d women, %d men), ",
                 "annual visits (Villeneuve et al. 2025, ",
                 "Alzheimers Dement, doi 10.1002/alz.70653)."),
          N_WOMEN + N_MEN, N_WOMEN, N_MEN),
  sprintf(paste0("- Dropout: monotone. Calibrated to the cohort ",
                 "paper: %.0f%% with at least one follow-up ",
                 "(retention %.2f at year 1) and %.0f%% still ",
                 "followed beyond 10 years (retention %.2f at year ",
                 "10), interpolated geometrically (%.2f%% attrition ",
                 "per wave) and extended at the same rate beyond ",
                 "year 10. Wave w is at year w - 1."),
          100 * RETENTION_YEAR1, RETENTION_YEAR1,
          100 * RETENTION_YEAR10, RETENTION_YEAR10,
          100 * DROPOUT_RATE),
  "- Retention by wave:",
  retention_lines.v,
  sprintf(paste0("- Matched design: %s (years 0-%d) is taken as ",
                 "the design matching PREVENT-AD's real follow-up ",
                 "(79%% followed beyond 10 years). Power for the ",
                 "other designs is reported alongside."),
          DESIGNS[[MATCHED_DESIGN]]$label,
          DESIGNS[[MATCHED_DESIGN]]$n_waves - 1),
  paste0("- Predictors: CVR_mimic and HVR_z drawn per subject from ",
         "a bivariate normal with the anchor model's ADNI baseline ",
         "means, SDs and correlation; Age_bl, EDUC ~ normal and ",
         "APOE4 ~ Bernoulli with ADNI moments. HVR_z is treated as ",
         "a person-level baseline predictor."),
  paste0("- Outcome: all fixed effects from the anchor fit (tested ",
         "term replaced by the target value), random intercept and ",
         "YRS slope from the fitted G matrix, Gaussian residual ",
         "with the fitted variance."),
  paste0("- Fitted model (3-way): `Y ~ YRS*CVR_mimic*HVR_z + ",
         "I(YRS^2)*CVR_mimic*HVR_z + Age_bl + EDUC + APOE4 ",
         "[+ SEX for pooled] + (YRS | PTID)`, REML, bobyqa, ",
         "Satterthwaite p for YRS:CVR_mimic:HVR_z (as in scripts ",
         "08 and 10b)."),
  paste0("- Fitted model (sex contrast): the same with `* SEX` on ",
         "both interaction blocks; Satterthwaite p for ",
         "YRS:CVR_mimic:HVR_z:SEXMale. Pooled sample only. Data ",
         "generated from the SEX-interaction refit (all fixed ",
         "effects including sex-specific 3-way terms), with the ",
         "4-way term set to the observed delta or a grid value."),
  paste0("- For the upper-bound (female-only) anchor, the generating ",
         "model has no SEX term; men are generated from the same ",
         "parameters (sex effect = 0) and SEX is still fitted as a ",
         "covariate in the pooled sample."),
  paste0("- Singular fits are retained (their Satterthwaite p is ",
         "still defined); the singular rate is reported per cell."),
  "",
  "## 3. Power by design x anchor x sample (3-way moderation)",
  "",
  md_table.fn(main_table.dt),
  "",
  "### Primary anchor: power at both Wald 95% CI bounds of beta",
  "",
  sprintf(paste0("The primary anchor beta = %.4f has SE = %.4f ",
                 "(N = %d subjects in the fitted frame), so power ",
                 "is reported at both Wald 95%% CI bounds, ",
                 "%.4f and %.4f."),
          ANCHORS$primary$beta3, ANCHORS$primary$se3,
          ANCHORS$primary$n_subj,
          ANCHORS$primary$beta3_ci[["lower"]],
          ANCHORS$primary$beta3_ci[["upper"]]),
  "",
  md_table.fn(ci_table_fmt.dt),
  "",
  "## 4. Minimum detectable effect at 80% power (3-way)",
  "",
  sprintf(paste0("|beta| grid: %s. MDE is linear interpolation on ",
                 "the (monotonised) power curve; ratio is MDE / ",
                 "|anchor beta|. Type I error is the rejection rate ",
                 "at beta = 0."),
          paste(MDE_GRID, collapse = ", ")),
  "",
  md_table.fn(mde_table.dt),
  "",
  "## 4b. Predictor-spread sensitivity (primary anchor)",
  "",
  paste0("PREVENT-AD sits nearer the UK Biobank norm than ADNI, so ",
         "HVR_z may be less dispersed; CVR_mimic may be ",
         "re-standardised on PREVENT-AD. Only the SDs matter for ",
         "the 3-way term (its coefficient and SE are invariant to ",
         "centring). Power at the ADNI-observed beta under ",
         "alternative predictor SDs:"),
  "",
  md_table.fn(spread_table.dt),
  "",
  "## 5. Sex-by-moderation contrast (4-way), pooled N = 387",
  "",
  "### Power at the observed delta",
  "",
  md_table.fn(sex_table.dt),
  "",
  "### Minimum detectable delta at 80% power (A-/CU refit)",
  "",
  sprintf("|delta| grid: %s.", paste(MDE_GRID_SEX, collapse = ", ")),
  "",
  md_table.fn(sex_mde_table.dt),
  "",
  "## 6. Full power curves",
  "",
  "### 3-way moderation (primary anchor)",
  "",
  md_table.fn(power_grid.dt[, .(
    Design = label_design.fn(design),
    Sample = label_sample.fn(sample),
    abs_beta = sprintf("%.4f", abs_beta),
    Power = sprintf("%.3f", power),
    Singular = sprintf("%.2f", singular_rate)
  )]),
  "",
  "### Sex contrast (A-/CU refit, pooled)",
  "",
  md_table.fn(power_sex.dt[tag == "sex_grid", .(
    Design = label_design.fn(design),
    abs_delta = sprintf("%.4f", abs_beta),
    Power = sprintf("%.3f", power),
    Singular = sprintf("%.2f", singular_rate)
  )]),
  "",
  "## 7. Caveats and discrepancies",
  "",
  sprintf(paste0("- **Predictor variance.** CVR_mimic was z-scored ",
                 "on the ADNI analysis sample (baseline SD %.2f in ",
                 "the A-/CU frame, %.2f in the A+ female frame) and ",
                 "ADNI HVR_z is shifted and dispersed relative to ",
                 "the UKB norm (A-/CU: M %.2f, SD %.2f; A+ females: ",
                 "M %.2f, SD %.2f). FRS has a ceiling in ADNI. ",
                 "PREVENT-AD (younger, cognitively unimpaired, ",
                 "family history) will have its own CVR and HVR ",
                 "spread; because the 3-way term scales with ",
                 "SD(CVR) x SD(HVR), power moves with those SDs ",
                 "(Section 4b)."),
          ANCHORS$primary$pred$cvr_sd,
          ANCHORS$secondary$pred$cvr_sd,
          ANCHORS$primary$pred$hvr_mean,
          ANCHORS$primary$pred$hvr_sd,
          ANCHORS$secondary$pred$hvr_mean,
          ANCHORS$secondary$pred$hvr_sd),
  sprintf(paste0("- **Follow-up length.** 8-12 annual waves is more ",
                 "follow-up than the ADNI models had. A-/CU Memory ",
                 "model: median %g visits per subject (IQR %g-%g), ",
                 "mean follow-up %.1f y (median %.1f, max %.1f). ",
                 "A+ female EXF model: median %g visits (IQR ",
                 "%g-%g), mean follow-up %.1f y (median %.1f, max ",
                 "%.1f). ADNI visits are irregular (6-12 month ",
                 "spacing early, then annual); the quadratic time ",
                 "terms are extrapolated when applied to 7-11 ",
                 "years of annual data."),
          ANCHORS$primary$visits$n_visits_median,
          ANCHORS$primary$visits$n_visits_q1,
          ANCHORS$primary$visits$n_visits_q3,
          ANCHORS$primary$visits$fu_mean,
          ANCHORS$primary$visits$fu_median,
          ANCHORS$primary$visits$fu_max,
          ANCHORS$secondary$visits$n_visits_median,
          ANCHORS$secondary$visits$n_visits_q1,
          ANCHORS$secondary$visits$n_visits_q3,
          ANCHORS$secondary$visits$fu_mean,
          ANCHORS$secondary$visits$fu_median,
          ANCHORS$secondary$visits$fu_max),
  sprintf(paste0("- **Anchor precision.** The primary anchor comes ",
                 "from a pooled-sex A-/CU model with %d subjects in ",
                 "the fitted frame (%d in the cohort; %d rows ",
                 "dropped for missing covariates). Its SE (%.4f) is ",
                 "large relative to beta (%.4f), so power at the ",
                 "two Wald CI bounds spans a wide range (Section 3). ",
                 "Note: the request described the A-/CU anchor as ",
                 "N = 300 with analytic N = 260; the LME fit used ",
                 "%d subjects. N = 260 is the analytic N of the ",
                 "A-/CU LGCM (supplementary), not of this LME."),
          ANCHORS$primary$n_subj,
          aneg.res$sample_info$n_subjects,
          aneg.res$sample_info$n_observations - ANCHORS$primary$n_obs,
          ANCHORS$primary$se3, ANCHORS$primary$beta3,
          ANCHORS$primary$n_subj),
  paste0("- **Winner's curse.** Both anchors are the largest of ",
         "several tested effects (3 domains; 6 sex x domain cells) ",
         "and were selected because they were significant, so the ",
         "true effect is likely smaller than the anchor; the ",
         "CI-bound and MDE results are the more conservative ",
         "reference."),
  paste0("- **Earlier disease timeline.** The anchor was estimated ",
         "in A-/CU participants in their early 70s. If the ",
         "moderation emerges with age or cumulative vascular ",
         "exposure, the effect over any fixed window will be ",
         "smaller in a cohort that starts a decade earlier; the ",
         "MDE grid shows how much smaller an effect can be and ",
         "still be detected. Mean cognitive slopes and slope ",
         "variance are also likely flatter, which lowers signal ",
         "and noise together; the net effect is unknown."),
  sprintf(paste0("- **Retention beyond year 10.** The 12-wave design ",
                 "extends the fitted attrition one year past the ",
                 "last known point (retention %.2f at wave 12). ",
                 "The 5-wave design uses the same curve truncated ",
                 "(retention %.2f at wave 5)."),
          DESIGNS$twelve$retention[12],
          DESIGNS$conservative$retention[5]),
  sprintf(paste0("- **Sex contrast anchor.** The observed delta ",
                 "(%.4f, SE %.4f) comes from a post hoc refit with ",
                 "%d subjects and is very imprecise; it is a ",
                 "data-derived placeholder for the contrast's ",
                 "magnitude, not an established effect. The MDE for ",
                 "delta is the more informative number."),
          CONTRASTS$primary$beta3, CONTRASTS$primary$se3,
          CONTRASTS$primary$n_subj),
  paste0("- **Upper-bound anchor.** The A+ female EXF model was fit ",
         "in amyloid-positive participants at a later disease ",
         "stage; it is reported only as an upper bound and is not ",
         "the effect PREVENT-AD would test."),
  paste0("- **Time-varying HVR.** In the ADNI models HVR_z was the ",
         "most recent MRI at each cognitive visit (time-varying); ",
         "here it is fixed at baseline. Within-person SD of HVR_z ",
         "in ADNI is small relative to between-person SD, so this ",
         "is a minor, slightly conservative simplification."),
  paste0("- **Normality.** Predictors are simulated as (bivariate) ",
         "normal; ADNI HVR_z is left-skewed. Covariate effects ",
         "(Age_bl, EDUC, APOE4) do not affect power for the tested ",
         "terms beyond degrees of freedom."),
  if (length(unlist(anchor_flags.lst)) > 0) {
    paste0("- **Anchor value discrepancies:** ",
           paste(unlist(anchor_flags.lst), collapse = "; "), ".")
  } else {
    paste0("- Extracted anchor values match the manuscript-reported ",
           "values to the reported rounding; no discrepancies.")
  },
  ""
)
writeLines(md.v, OUT_MD)
log_info("Wrote %s", OUT_MD)
log_info("Wrote %s", OUT_RDS)

# ---------------------------------------------------------
# Console summary
# ---------------------------------------------------------
cat("\n==== POWER BY DESIGN x ANCHOR x SAMPLE (3-way) ====\n")
print(main_table.dt, nrows = 100)
cat("\n==== PRIMARY ANCHOR: POWER AT BOTH CI BOUNDS ====\n")
print(ci_table_fmt.dt, nrows = 100)
cat("\n==== MINIMUM DETECTABLE EFFECT (80% power, 3-way) ====\n")
print(mde_table.dt, nrows = 100)
cat("\n==== PREDICTOR-SPREAD SENSITIVITY ====\n")
print(spread_table.dt, nrows = 100)
cat("\n==== SEX-BY-MODERATION CONTRAST (4-way) ====\n")
print(sex_table.dt, nrows = 100)
print(sex_mde_table.dt, nrows = 100)
cat("\n")
cat(headline_mde, "\n")
cat(headline_waves, "\n")
cat(final_sentence, "\n")
cat(sex_sentence, "\n")
cat(ci_sentence, "\n")
cat(contrast_sentence, "\n")
cat(secondary_sentence, "\n\n")

log_script_end("16_power_preventad.R", success = TRUE)
