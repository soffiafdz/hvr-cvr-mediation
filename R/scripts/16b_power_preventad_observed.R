#!/usr/bin/env Rscript

# =============================================================================
# 16b_power_preventad_observed.R
# =============================================================================
# Variant of 16_power_preventad.R (unchanged) with a visit design
# calibrated to PREVENT-AD's recorded follow-up: median 8.0 years
# (SD 3.1) and 3,162 cognitive assessments for 387 participants.
# The 12-wave design in script 16 implies ~4,045 assessments.
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

log_script_start("16b_power_preventad_observed.R")
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
OUT_RDS <- here("outputs/power_preventad_observed.rds")
OUT_MD <- here("outputs/power_preventad_observed.md")

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
# 2. Observed-calibrated visit design
# ============================================================
log_section("Observed-calibrated visit design")
# Verified PREVENT-AD figures (grants/FACTS.md): 387 participants,
# 3,162 RBANS assessments, follow-up median 8.0 y (SD 3.1),
# 95% with >= 1 follow-up, follow-up up to 15 years.
N_ASSESS_OBS <- 3162
FU_MEDIAN <- 8.0
FU_SD <- 3.1
P_NO_FU <- 0.05
FU_MAX <- 14L
draw_fu.fn <- function(n) {
  # Follow-up length in whole years (last visit year)
  fu.v <- pmin(pmax(round(rnorm(n, FU_MEDIAN, FU_SD)), 1L), FU_MAX)
  fu.v[runif(n) < P_NO_FU] <- 0L
  as.integer(fu.v)
}
# Calibrate the share of scheduled follow-up visits that were
# completed so that expected assessments per person match
# 3,162 / 387; baseline is always present.
set.seed(SEED)
fu_cal.v <- draw_fu.fn(1e6)
KEEP_P <- (N_ASSESS_OBS / (N_WOMEN + N_MEN) - 1) / mean(fu_cal.v)
stopifnot(KEEP_P > 0, KEEP_P <= 1)
log_info("Follow-up years: median %.1f, mean %.2f, SD %.2f",
         median(fu_cal.v), mean(fu_cal.v), sd(fu_cal.v))
log_info("Share of follow-up visits completed: %.3f", KEEP_P)
draw_visits.fn <- function(n) {
  fu.v <- draw_fu.fn(n)
  lapply(fu.v, function(f) {
    if (f == 0L) return(0)
    yrs.v <- seq_len(f)
    c(0, yrs.v[runif(f) < KEEP_P])
  })
}
set.seed(SEED + 1L)
chk.lst <- lapply(1:200, function(i) sum(lengths(draw_visits.fn(387L))))
ASSESS_SIM <- mean(unlist(chk.lst))
log_info("Simulated assessments per dataset: mean %.0f (target %d)",
         ASSESS_SIM, N_ASSESS_OBS)
DESIGNS <- list(observed = list(
  label = "Observed-calibrated (median 8 y follow-up, 3,162 assessments)"
))

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
  vis.lst <- draw_visits.fn(n)
  id.v <- rep(seq_len(n), lengths(vis.lst))
  yrs.v <- unlist(vis.lst)
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
# 4. Runs
# ============================================================
log_section("Power, observed-calibrated design")
runs.lst <- list()
for (sn in names(SAMPLES)) {
  runs.lst[[paste("anchor", sn)]] <- run_condition.fn(
    "primary", "observed", sn, ANCHORS$primary$beta3, N_REPS, "anchor")
}
runs.lst[["ci_lower pooled"]] <- run_condition.fn(
  "primary", "observed", "pooled",
  ANCHORS$primary$beta3_ci[["lower"]], N_REPS, "ci_lower")
runs.lst[["ci_upper pooled"]] <- run_condition.fn(
  "primary", "observed", "pooled",
  ANCHORS$primary$beta3_ci[["upper"]], N_REPS, "ci_upper")
runs.lst[["tight_hvr pooled"]] <- run_condition.fn(
  "primary", "observed", "pooled", ANCHORS$primary$beta3, N_REPS,
  "spread", scenario = "tight_hvr",
  pred_override = SPREAD_SCENARIOS$tight_hvr$override)
runs.lst[["sex contrast"]] <- run_condition.fn(
  "primary", "observed", "pooled", CONTRASTS$primary$beta3, N_REPS,
  "sex_anchor", contrast = TRUE)
for (b in MDE_GRID) {
  runs.lst[[paste("mde", b)]] <- run_condition.fn(
    "primary", "observed", "pooled", -b, N_REPS, "mde_grid")
}
power.dt <- rbindlist(lapply(runs.lst, `[[`, "summary"))
grid.dt <- power.dt[tag == "mde_grid"][, abs_beta := abs(beta3)]
mde_80 <- interp_mde.fn(grid.dt$abs_beta, grid.dt$power, TARGET_POWER)

results.lst <- list(
  design = list(n_assess_target = N_ASSESS_OBS,
                n_assess_sim = ASSESS_SIM, keep_p = KEEP_P,
                fu_median = FU_MEDIAN, fu_sd = FU_SD,
                p_no_fu = P_NO_FU, fu_max = FU_MAX),
  power = power.dt, mde_80 = mde_80,
  anchor_beta = ANCHORS$primary$beta3,
  seed = SEED, n_reps = N_REPS, timestamp = Sys.time()
)
write_rds_safe(results.lst, OUT_RDS, "PREVENT-AD power, observed design")

f.fn <- function(x) sprintf("%.0f%%", 100 * x)
pw.fn <- function(t, s, sc = "adni") {
  power.dt[tag == t & sample == s & scenario == sc, power]
}
md.v <- c(
  "# Power, observed-calibrated PREVENT-AD design",
  "",
  sprintf(paste0(
    "Generated by `R/scripts/16b_power_preventad_observed.R` on %s; ",
    "%d replicates per condition; seed %d."
  ), format(Sys.time(), "%Y-%m-%d %H:%M"), N_REPS, SEED),
  "",
  paste0(
    "Same anchors, generating model, fitted model and test as script 16 ",
    "(A-/CU pooled Memory, CVR_mimic x HVR_z x YRS). ",
    "Only the visit structure differs:"
  ),
  sprintf(paste0(
    "- follow-up years ~ round(Normal(%.1f, %.1f)) clamped to 1-%d, ",
    "with %.0f%% having no follow-up;"
  ), FU_MEDIAN, FU_SD, FU_MAX, 100 * P_NO_FU),
  sprintf(paste0(
    "- each scheduled annual follow-up visit completed with probability ",
    "%.3f, calibrated so assessments match %d ",
    "(simulated mean %.0f per dataset of 387)."
  ), KEEP_P, N_ASSESS_OBS, ASSESS_SIM),
  "",
  "| Condition | Power |",
  "|---|---|",
  sprintf("| Pooled (n = 387), anchor beta %.4f | %s |",
          ANCHORS$primary$beta3, f.fn(pw.fn("anchor", "pooled"))),
  sprintf("| Women only (n = 279) | %s |", f.fn(pw.fn("anchor", "women"))),
  sprintf("| Men only (n = 108) | %s |", f.fn(pw.fn("anchor", "men"))),
  sprintf("| Pooled, beta at CI bound nearer zero | %s |",
          f.fn(power.dt[tag == "ci_upper", power])),
  sprintf("| Pooled, beta at CI bound farther from zero | %s |",
          f.fn(power.dt[tag == "ci_lower", power])),
  sprintf("| Pooled, HVR_z SD 0.8 and CVR SD 1.0 | %s |",
          f.fn(power.dt[tag == "spread", power])),
  sprintf("| Sex x moderation contrast (pooled) | %s |",
          f.fn(power.dt[tag == "sex_anchor", power])),
  "",
  sprintf(paste0(
    "Minimum detectable |beta| at 80%% power (pooled): %.4f ",
    "(%.2f x the anchor). Type I error at beta = 0: %s."
  ), mde_80, mde_80 / abs(ANCHORS$primary$beta3),
  f.fn(grid.dt[abs_beta == 0, power])),
  "",
  "## All conditions",
  "",
  "| tag | sample | scenario | beta | power | MC SE | singular | errors |",
  "|---|---|---|---|---|---|---|---|",
  power.dt[, sprintf(
    "| %s | %s | %s | %.4f | %.3f | %.3f | %.2f | %.3f |",
    tag, sample, scenario, beta3, power, power_mc_se, singular_rate,
    error_rate
  )]
)
writeLines(md.v, OUT_MD)
log_info("Wrote %s", OUT_MD)
log_script_end("16b_power_preventad_observed.R", success = TRUE)
