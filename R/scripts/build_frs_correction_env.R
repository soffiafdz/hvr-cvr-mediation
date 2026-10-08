#!/usr/bin/env Rscript
# =============================================================================
# build_frs_correction_env.R - Numbers and figures for the FRS correction deck
# =============================================================================
# Gathers every number and figure that reports-src/frs_correction_slides.qmd
# shows about the ADSP-PHC Framingham sex swap: the checks on the released
# variable, the change in the FRS distribution, and the results before and
# after the corrected re-run. The slides do no analysis of their own.
#
# "Before" is the snapshot taken ahead of the re-run (config
# data.archive.pre_frs_fix; published ADSP-PHC FRS). "After" is the current
# output (FRS recomputed with each participant's own sex, R/utils/frs.R).
# See docs/frs_sex_swap_audit.md.
#
# INPUTS:
#   - ADSP-PHC CVRF release (config data.raw.adsp_phc_cvrf)
#   - ADNIMERGE package (PTGENDER)
#   - data/derivatives/lme/cohort_hvr.rds
#   - <archive>/data_derivatives/lme/cvr_mimic_scores.rds (published FRS)
#   - outputs/manuscript_env.rds and <archive>/outputs/manuscript_env.rds
#   - mediation and LME bootstrap checkpoint folders (current and archived)
# OUTPUTS:
#   - outputs/frs_correction_env.rds
#   - outputs/figures/frs_correction/*.png (slide figures, 16:9)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/frs.R"))
source(here("R/utils/plotting.R"))

log_script_start("build_frs_correction_env.R")
config <- load_config()
validate_config(config)
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "frs_correction_env", default = FALSE
)
output.path <- get_data_path("outputs", "frs_correction_env")
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE; skipping")
  log_script_end("build_frs_correction_env.R", success = TRUE)
  quit(status = 0)
}

# --- Configuration ---
ARCHIVE <- get_data_path("archive", "pre_frs_fix")
FIG_DIR <- get_data_path("outputs", "frs_correction_figures")
OUT_DIR <- get_data_path("outputs", "base_dir")
FRS_COL <- get_parameter("frs", "column")
SEXES <- c("Male", "Female")  # PHC_Sex 1 = Male, 2 = Female
DOMAINS <- unlist(get_parameter("cognitive_domains"))
DOMAIN_LABELS <- unlist(get_parameter("domain_labels"))[DOMAINS]
CAP <- 30  # FRS clamp (percent), as in CVrisk::ascvd_10y_frs_simple
AGE_MAX <- 74  # stock CVrisk returns NA above this age
COLORS <- manuscript_colors()
ensure_directory(FIG_DIR)

# -----------------------------------------------------------------------------
log_section("Part 1: Released variable against both sex equations")
# -----------------------------------------------------------------------------
cvrf.dt <- read_csv_safe(
  get_data_path("raw", "adsp_phc_cvrf"),
  select = c("RID", "PTID", FRS_COL, FRS_INPUT_COLS),
  description = "ADSP-PHC CVRF"
)
stopifnot(all(cvrf.dt$PHC_Sex %in% c(1, 2, NA)))
own_sex.v <- SEXES[cvrf.dt$PHC_Sex]
other_sex.v <- rev(SEXES)[cvrf.dt$PHC_Sex]
frs_with.fn <- function(sex.v) {
  compute_frs.fn(
    sex.v, cvrf.dt$PHC_Age_CardiovascularRisk, cvrf.dt$PHC_BMI,
    cvrf.dt$PHC_SBP, cvrf.dt$PHC_BP_Med, cvrf.dt$PHC_Smoker,
    cvrf.dt$PHC_Diabetes
  )
}
cvrf.dt[, `:=`(
  sex = own_sex.v,
  published = get(FRS_COL),
  own_eq = frs_with.fn(own_sex.v),
  other_eq = frs_with.fn(other_sex.v)
)]
scored.dt <- cvrf.dt[!is.na(published)]
by_sex.dt <- scored.dt[, .(
  n = .N,
  max_diff_own = max(abs(own_eq - published)),
  max_diff_other = max(abs(other_eq - published)),
  n_match_other = sum(abs(other_eq - published) <= 0.01)
), by = sex]
detection.lst <- list(
  n_rows = nrow(cvrf.dt),
  n_scored = nrow(scored.dt),
  n_match_other = scored.dt[abs(other_eq - published) <= 0.01, .N],
  max_diff_other = scored.dt[, max(abs(other_eq - published))],
  by_sex = by_sex.dt,
  n_over74 = scored.dt[PHC_Age_CardiovascularRisk > AGE_MAX, .N],
  oldest = scored.dt[, max(PHC_Age_CardiovascularRisk)],
  n_30_74 = scored.dt[PHC_Age_CardiovascularRisk %between% c(30, 74), .N],
  published_range = range(scored.dt$published)
)
stopifnot(detection.lst$n_match_other == detection.lst$n_scored)
log_info(
  "Opposite-sex equation reproduces %d of %d published values",
  detection.lst$n_match_other, detection.lst$n_scored
)

# Worked example: the first participant in the release (RID removed)
#' Framingham risk (%) before rounding and the 1-30 clamp
#' @param row.dt One CVRF row
#' @param sex_lab "Male" or "Female"
#' @return Numeric risk in percent
unclamped_risk.fn <- function(row.dt, sex_lab) {
  k.dt <- FRS_COEF[SEX == sex_lab]
  sum_x <- k.dt$ln_age * log(row.dt$PHC_Age_CardiovascularRisk) +
    k.dt$ln_bmi * log(row.dt$PHC_BMI) +
    ifelse(
      row.dt$PHC_BP_Med == 1, k.dt$ln_sbp_treated, k.dt$ln_sbp_untreated
    ) * log(row.dt$PHC_SBP) +
    k.dt$smoker * row.dt$PHC_Smoker + k.dt$diabetes * row.dt$PHC_Diabetes
  100 * (1 - k.dt$baseline_survival^exp(sum_x - k.dt$group_mean))
}
ex.dt <- cvrf.dt[RID == REMOVED_ID]
stopifnot(nrow(ex.dt) == 1, ex.dt$sex == "Male")
example.lst <- list(
  rid = ex.dt$RID,
  sex = ex.dt$sex,
  age = ex.dt$PHC_Age_CardiovascularRisk,
  bmi = ex.dt$PHC_BMI,
  sbp = ex.dt$PHC_SBP,
  bp_med = ex.dt$PHC_BP_Med,
  smoker = ex.dt$PHC_Smoker,
  diabetes = ex.dt$PHC_Diabetes,
  published = ex.dt$published,
  female_eq = ex.dt$other_eq,
  male_eq = ex.dt$own_eq,
  male_eq_unclamped = unclamped_risk.fn(ex.dt, "Male")
)
stopifnot(abs(example.lst$female_eq - example.lst$published) <= 0.01)

coef.dt <- FRS_COEF[, .(
  SEX, ln_age, ln_bmi, ln_sbp_untreated, ln_sbp_treated, smoker, diabetes,
  baseline_survival
)]

# -----------------------------------------------------------------------------
log_section("Part 2: Sex labels against ADNIMERGE")
# -----------------------------------------------------------------------------
am.dt <- unique(as.data.table(ADNIMERGE::adnimerge)[
  , .(RID, PTID, PTGENDER = as.character(PTGENDER))
])
phc_labels.dt <- merge(
  unique(cvrf.dt[, .(RID, sex)]), unique(am.dt[, .(RID, PTGENDER)]),
  by = "RID"
)[, .N, by = .(sex, PTGENDER)][order(sex, PTGENDER)]

cohort.dt <- as.data.table(read_rds_safe(
  get_data_path("derivatives", "lme_cohort_hvr"), "LME cohort"
))
bl.dt <- cohort.dt[order(PTID, YRS_from_bl)][, .SD[1], by = PTID]
pipeline_labels.dt <- merge(
  bl.dt[, .(PTID, SEX)], unique(am.dt[, .(PTID, PTGENDER)]), by = "PTID"
)[, .N, by = .(SEX, PTGENDER)][order(SEX, PTGENDER)]
sexcheck.lst <- list(
  phc = phc_labels.dt,
  phc_discordant = phc_labels.dt[sex != PTGENDER, sum(N)],
  pipeline = pipeline_labels.dt,
  pipeline_discordant = pipeline_labels.dt[SEX != PTGENDER, sum(N)],
  phc_total = phc_labels.dt[, sum(N)]
)

# -----------------------------------------------------------------------------
log_section("Part 3: FRS distribution in the analysis cohort")
# -----------------------------------------------------------------------------
old_sc.dt <- as.data.table(read_rds_safe(
  path_join(ARCHIVE, "data_derivatives/lme/cvr_mimic_scores.rds"),
  "Pre-correction CVR scores (published FRS)"
))
base.dt <- merge(
  bl.dt[, .(RID, SEX, AGE, FRS_new = FRS_raw, CVR = CVR_mimic)],
  old_sc.dt[, .(RID, FRS_old = FRS)], by = "RID"
)
check.dt <- merge(
  base.dt[, .(RID, FRS_new, FRS_old)],
  cvrf.dt[, .(RID, own_eq, published)], by = "RID"
)
stopifnot(
  isTRUE(all.equal(check.dt$FRS_new, check.dt$own_eq)),
  isTRUE(all.equal(check.dt$FRS_old, check.dt$published))
)
validate_num_range(base.dt, "FRS_new", 1, CAP, na_allowed = FALSE)
validate_num_range(base.dt, "FRS_old", 1, CAP, na_allowed = FALSE)

# Age at the FRS assessment, for the share scored beyond stock CVrisk's range
age_frs.dt <- merge(
  base.dt[, .(RID)], cvrf.dt[, .(RID, age = PHC_Age_CardiovascularRisk)],
  by = "RID"
)
stopifnot(nrow(age_frs.dt) == nrow(base.dt), !anyNA(age_frs.dt$age))
cohort_over74.lst <- list(
  n = age_frs.dt[age > AGE_MAX, .N],
  total = nrow(age_frs.dt),
  pct = 100 * age_frs.dt[, mean(age > AGE_MAX)],
  age_max = AGE_MAX
)

pct_cap.fn <- function(v) 100 * mean(v == CAP)
dist.dt <- rbind(
  base.dt[, .(
    group = "All", n = .N,
    cap_old = pct_cap.fn(FRS_old), cap_new = pct_cap.fn(FRS_new),
    mean_old = mean(FRS_old), mean_new = mean(FRS_new),
    rho_old_new = cor(FRS_old, FRS_new, method = "spearman")
  )],
  base.dt[, .(
    n = .N,
    cap_old = pct_cap.fn(FRS_old), cap_new = pct_cap.fn(FRS_new),
    mean_old = mean(FRS_old), mean_new = mean(FRS_new),
    rho_old_new = cor(FRS_old, FRS_new, method = "spearman")
  ), by = .(group = SEX)]
)

# Agreement with CVR_mimic: pooled, within sex (sex means removed), by sex
base.dt[, `:=`(
  CVR_w = CVR - mean(CVR),
  FRS_old_w = FRS_old - mean(FRS_old),
  FRS_new_w = FRS_new - mean(FRS_new)
), by = SEX]
agree.dt <- rbind(
  base.dt[, .(
    group = "Pooled", r_old = cor(CVR, FRS_old), r_new = cor(CVR, FRS_new)
  )],
  base.dt[, .(
    group = "Within sex",
    r_old = cor(CVR_w, FRS_old_w), r_new = cor(CVR_w, FRS_new_w)
  )],
  base.dt[, .(r_old = cor(CVR, FRS_old), r_new = cor(CVR, FRS_new)),
          by = .(group = SEX)]
)
cvr_by_sex.dt <- base.dt[, .(cvr_mean = mean(CVR)), by = SEX]

# -----------------------------------------------------------------------------
log_section("Part 4: Results before and after (manuscript environments)")
# -----------------------------------------------------------------------------
old_env.lst <- read_rds_safe(
  path_join(ARCHIVE, "outputs/manuscript_env.rds"), "Pre-correction env"
)
new_env.lst <- read_rds_safe(
  path_join(OUT_DIR, "manuscript_env.rds"), "Current env"
)
pair.fn <- function(get.fn) {
  c(old = get.fn(old_env.lst), new = get.fn(new_env.lst))
}
headline.lst <- list(
  n_subjects = pair.fn(function(e) e$n_subjects),
  frs_ceiling_all = pair.fn(function(e) e$frs_ceiling_all),
  frs_ceiling_male = pair.fn(function(e) e$frs_ceiling_male),
  frs_ceiling_female = pair.fn(function(e) e$frs_ceiling_female),
  frs_age_r = pair.fn(function(e) e$frs_age_r),
  cvr_age_cor = pair.fn(function(e) e$cvr_age_cor),
  cvr_frs_cor = pair.fn(function(e) e$cvr_frs_cor),
  r2_age = pair.fn(function(e) e$frs_var_decomp$r2_age),
  r2_cvrf = pair.fn(function(e) e$frs_var_decomp$r2_cvrf)
)
stopifnot(
  isTRUE(all.equal(
    unname(headline.lst$frs_ceiling_all["new"]),
    dist.dt[group == "All", cap_new], tolerance = 1e-3
  ))
)

COUNT_SPECS.lst <- list(
  "FRS a-path" = c("abs_stats", "frs_a_sig", 6),
  "CVR_mimic a-path" = c("abs_stats", "cvr_a_sig", 6),
  "FRS b-path" = c("abs_stats", "frs_b_sig", 6),
  "CVR_mimic b-path" = c("abs_stats", "cvr_b_sig", 6),
  "FRS indirect effect" = c("abs_stats", "frs_ind_sig", 6),
  "CVR_mimic indirect effect" = c("abs_stats", "cvr_ind_sig", 6),
  "FRS LME interaction (FDR)" = c("fdr_summary.lst", "frs_n_fdr", 6),
  "CVR_mimic LME interaction (FDR)" = c("fdr_summary.lst", "cvr_n_fdr", 6),
  "Coupling, men" = c("abs_stats", "cpl_m_sig", 3),
  "Coupling, women" = c("abs_stats", "cpl_f_sig", 3)
)
counts.dt <- rbindlist(lapply(names(COUNT_SPECS.lst), function(lab) {
  s.v <- COUNT_SPECS.lst[[lab]]
  data.table(
    measure = lab,
    old = old_env.lst[[s.v[1]]][[s.v[2]]],
    new = new_env.lst[[s.v[1]]][[s.v[2]]],
    total = as.integer(s.v[3])
  )
}))

# LGCM mediation paths per model
grid.dt <- CJ(sex = SEXES, domain = DOMAINS, sorted = FALSE)
path.fn <- function(env.lst, predictor, path, version) {
  res.lst <- env.lst$lgcm_results$results
  grid.dt[, {
    p.lst <- require_output(
      res.lst[[paste(sex, domain, predictor, sep = "_")]][[path]],
      paste(sex, domain, predictor, path)
    )
    .(predictor = predictor, version = version, est = p.lst$est,
      lower = p.lst$boot_ci_lower, upper = p.lst$boot_ci_upper)
  }, by = .(sex, domain)]
}
apath.dt <- rbind(
  path.fn(old_env.lst, "FRS", "a_path", "Released"),
  path.fn(new_env.lst, "FRS", "a_path", "Corrected"),
  path.fn(new_env.lst, "CVR_mimic", "a_path", "Corrected")
)
apath.dt[, sig := lower > 0 | upper < 0]
apath_by_sex.dt <- apath.dt[predictor == "FRS", .(
  n_sig = sum(sig), total = .N, mean_est = mean(est),
  mean_width = mean(upper - lower)
), by = .(sex, version)]
width_ratio.v <- apath_by_sex.dt[, .(
  ratio = mean_width[version == "Corrected"] /
    mean_width[version == "Released"]
), by = sex][, setNames(ratio, sex)]
indirect.dt <- rbind(
  path.fn(old_env.lst, "FRS", "indirect", "Released"),
  path.fn(new_env.lst, "FRS", "indirect", "Corrected")
)
indirect.dt[, sig := lower > 0 | upper < 0]

# LME three-way interactions (YRS x CVR x HVR_z), primary stratified models
lme.fn <- function(env.lst, version) {
  as.data.table(env.lst$lme_summary.dt)[
    Analysis == "Primary_Stratified" & CVR_Measure %in% c("FRS", "CVR_mimic"),
    .(predictor = CVR_Measure, sex = Sex, domain = Domain, version,
      beta = Beta, se = SE, p = p_value, p_fdr, sig_fdr = Significant_FDR)
  ]
}
lme.dt <- rbind(lme.fn(old_env.lst, "Released"),
                lme.fn(new_env.lst, "Corrected"))
lme_n.v <- unique(as.data.table(new_env.lst$lme_summary.dt)[
  Analysis == "Primary_Stratified" & CVR_Measure == "FRS",
  .(Sex, N_subjects)
])[, setNames(N_subjects, Sex)]
stopifnot(length(lme_n.v) == 2, setequal(names(lme_n.v), SEXES))
beta_ratio.fn <- function(d.dt) {
  d.dt[predictor == "FRS", mean(abs(beta))] /
    d.dt[predictor == "CVR_mimic", mean(abs(beta))]
}
beta_ratio.v <- c(
  old = beta_ratio.fn(lme.dt[version == "Released"]),
  new = beta_ratio.fn(lme.dt[version == "Corrected"])
)
comparison.fn <- function(env.lst) {
  sum(vapply(env.lst$lme_hvr_z.res$comparison, function(x) {
    isTRUE(x$significant)
  }, logical(1)))
}
lme_diff.v <- c(old = comparison.fn(old_env.lst),
                new = comparison.fn(new_env.lst))

# Amyloid-negative, cognitively unimpaired controls
aneg.fn <- function(env.lst, version) {
  as.data.table(env.lst$aneg_lme.res$summary_table)[
    , .(predictor = CVR_Measure, domain = Domain, version, beta = Beta,
        p = p_value, p_fdr, sig_fdr = Significant_FDR)
  ]
}
aneg.lst <- list(
  n = pair.fn(function(e) e$aneg_lme_summary.lst$n),
  frs_fdr = pair.fn(function(e) e$aneg_lme_summary.lst$frs_fdr),
  cvr_fdr = pair.fn(function(e) e$aneg_lme_summary.lst$cvr_fdr),
  lgcm_a_sig = pair.fn(function(e) e$aneg_lgcm_summary.lst$a_sig),
  lgcm_a_total = pair.fn(function(e) e$aneg_lgcm_summary.lst$a_total),
  lme = rbind(aneg.fn(old_env.lst, "Released"),
              aneg.fn(new_env.lst, "Corrected"))
)

# -----------------------------------------------------------------------------
log_section("Part 5: Re-run bookkeeping")
# -----------------------------------------------------------------------------
chk_count.fn <- function(dir, pattern = "\\.rds$") {
  length(list.files(dir, pattern = pattern))
}
MED_CHK <- "models/results/lgcm/mediation_checkpoints"
LME_CHK <- "models/results/lme/lme_boot_checkpoints"
rerun.lst <- list(
  med_total = chk_count.fn(here(MED_CHK)),
  med_frs = chk_count.fn(here(MED_CHK), "_FRS\\.rds$"),
  med_frs_archived = chk_count.fn(
    path_join(ARCHIVE, "models_results/lgcm/mediation_checkpoints"),
    "_FRS\\.rds$"
  ),
  lme_archived = chk_count.fn(
    path_join(ARCHIVE, "models_results/lme/lme_boot_checkpoints")
  ),
  n_lgcm = sapply(new_env.lst$med_n.lst, identity)
)

# -----------------------------------------------------------------------------
log_section("Part 6: Slide figures")
# -----------------------------------------------------------------------------
VERSION_LEVELS <- c("Released", "Corrected")
SERIES_COLORS <- c(
  "FRS, released" = "#8C8C8C",
  "FRS, corrected" = COLORS$frs_pred,
  "CVR_mimic" = COLORS$cvr_pred
)
slide_theme <- theme_publication(base_size = 20) +
  theme(legend.position = "top", legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 20))

cap.dt <- melt(
  dist.dt[group %in% SEXES, .(group, Released = cap_old,
                              Corrected = cap_new)],
  id.vars = "group", variable.name = "version", value.name = "pct"
)
cap.dt[, `:=`(
  version = factor(version, levels = VERSION_LEVELS),
  group = factor(group, levels = SEXES, labels = c("Men", "Women"))
)]
ceiling.plt <- ggplot(cap.dt, aes(version, pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4,
            size = 7) +
  scale_fill_manual(values = c(Men = COLORS$male, Women = COLORS$female)) +
  scale_x_discrete(labels = c(Released = "Released score",
                              Corrected = "Corrected score")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = sprintf("Participants at the %d%% cap", CAP)) +
  slide_theme

forest_data.fn <- function(d.dt) {
  d.dt[, `:=`(
    series = fifelse(
      predictor == "FRS", paste0("FRS, ", tolower(version)), "CVR_mimic"
    ),
    domain_lab = factor(DOMAIN_LABELS[domain], levels = rev(DOMAIN_LABELS)),
    sex_lab = factor(sex, levels = SEXES, labels = c("Men", "Women"))
  )]
  d.dt[, series := factor(series, levels = names(SERIES_COLORS))]
  # Significant estimates are filled with the series colour, others white
  d.dt[, point_fill := fifelse(sig, as.character(series), "Not significant")]
  d.dt
}
forest.fn <- function(d.dt, x_lab) {
  dodge <- position_dodge(width = 0.6, orientation = "y", reverse = TRUE)
  ggplot(d.dt, aes(est, domain_lab, colour = series, group = series)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.6) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0, linewidth = 1.2, position = dodge) +
    geom_point(aes(fill = point_fill), shape = 21, size = 4.5, stroke = 1.6,
               position = dodge) +
    facet_wrap(~sex_lab, nrow = 1) +
    scale_colour_manual(values = SERIES_COLORS, drop = FALSE) +
    scale_fill_manual(values = c(SERIES_COLORS, "Not significant" = "white"),
                      guide = "none") +
    guides(colour = guide_legend(
      override.aes = list(fill = unname(SERIES_COLORS))
    )) +
    scale_x_continuous(
      labels = scales::label_number(style_negative = "minus")
    ) +
    labs(x = x_lab, y = NULL) +
    slide_theme
}
apath_plot.dt <- forest_data.fn(copy(apath.dt))
apath.plt <- forest.fn(apath_plot.dt, "a-path: risk \u2192 HVR slope")
lme_plot.dt <- forest_data.fn(copy(lme.dt)[, `:=`(
  est = beta, lower = beta - 1.96 * se, upper = beta + 1.96 * se,
  sig = sig_fdr
)])
lme_plot.dt <- lme_plot.dt[!(predictor == "CVR_mimic" &
                               version == "Released")]
lme.plt <- forest.fn(
  lme_plot.dt, "Years \u00d7 risk \u00d7 HVR interaction (\u03b2)"
)

figures.lst <- list(
  ceiling = path_join(FIG_DIR, "frs_ceiling_by_sex.png"),
  apath = path_join(FIG_DIR, "frs_apath_forest.png"),
  lme = path_join(FIG_DIR, "frs_lme_forest.png")
)
ggsave(figures.lst$ceiling, ceiling.plt, width = 13.33, height = 6.5,
       dpi = 150)
# Shorter, to leave room for a note under it on its slide
ggsave(figures.lst$apath, apath.plt, width = 13.33, height = 5.2, dpi = 150)
ggsave(figures.lst$lme, lme.plt, width = 13.33, height = 6.5, dpi = 150)
log_info("Saved slide figures to %s", FIG_DIR)

# -----------------------------------------------------------------------------
log_section("Part 7: Save")
# -----------------------------------------------------------------------------
frs_correction_env.lst <- list(
  frs_column = FRS_COL,
  cap = CAP,
  detection = detection.lst,
  example = example.lst,
  coefficients = coef.dt,
  cohort_over74 = cohort_over74.lst,
  sexcheck = sexcheck.lst,
  distribution = dist.dt,
  agreement = agree.dt,
  cvr_by_sex = cvr_by_sex.dt,
  headline = headline.lst,
  counts = counts.dt,
  apath = apath.dt,
  apath_by_sex = apath_by_sex.dt,
  apath_width_ratio = width_ratio.v,
  indirect = indirect.dt,
  lme = lme.dt,
  lme_n = lme_n.v,
  beta_ratio = beta_ratio.v,
  lme_diff_sig = lme_diff.v,
  aneg = aneg.lst,
  rerun = rerun.lst,
  domain_labels = DOMAIN_LABELS,
  figures = figures.lst,
  timestamp = Sys.time()
)
write_rds_safe(frs_correction_env.lst, output.path, "FRS correction env")
log_script_end("build_frs_correction_env.R", success = TRUE)
