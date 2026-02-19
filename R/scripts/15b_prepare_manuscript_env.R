# =========================================================
# 15b: Prepare Manuscript Environment
# =========================================================
# Loads all pipeline outputs, computes summary statistics
# and adapter objects, and saves a single
# outputs/manuscript_env.rds for the QMD to load.
#
# Usage: Rscript R/scripts/15b_prepare_manuscript_env.R
# =========================================================

library(here)
library(data.table)

source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/tables.R"))

config <- load_config()

# ---------------------------------------------------------
# Helper: load or stop
# ---------------------------------------------------------
load_if_exists <- function(path, name = "data") {
  if (!file.exists(path)) {
    stop(sprintf(
      "%s not found: %s\nRun pipeline first.",
      name, path
    ))
  }
  readRDS(path)
}

# ---------------------------------------------------------
# 1. Directory paths
# ---------------------------------------------------------
derivatives_dir <- here("data/derivatives/lme")
lme_results_dir <- here("models/results/lme")
lgcm_results_dir <- here("models/results/lgcm")

# ---------------------------------------------------------
# 2. Load all RDS files
# ---------------------------------------------------------
cohort.dt <- load_if_exists(
  file.path(derivatives_dir, "cohort_hvr.rds"),
  "Cohort"
)
cvr_scores <- load_if_exists(
  file.path(
    derivatives_dir, "cvr_mimic_scores.rds"
  ),
  "CVR scores"
)
lme_hvr_z.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "lme_hvr_z_results.rds"
  ),
  "LME HVR z"
)
lme_hc_z.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "lme_hc_z_results.rds"
  ),
  "LME HC z"
)
frs_edt <- load_if_exists(
  file.path(
    lme_results_dir,
    "lme_edt_hvr_z_frs_results.rds"
  ),
  "FRS EDT"
)
cvr_edt <- load_if_exists(
  file.path(
    lme_results_dir,
    "lme_edt_hvr_z_cvr_results.rds"
  ),
  "CVR EDT"
)
attrition_analysis <- load_if_exists(
  file.path(
    lme_results_dir,
    "attrition_analysis.rds"
  ),
  "Attrition"
)
lme_raw_sens.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "lme_raw_sensitivity_results.rds"
  ),
  "Raw sensitivity"
)
edt_sens.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "edt_early_middle_sensitivity.rds"
  ),
  "EDT sensitivity"
)
frs_problem.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "frs_problem_analysis.rds"
  ),
  "FRS problem"
)
cvr_model.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "cvr_mimic_model.rds"
  ),
  "CVR MIMIC model"
)
cvr_mi.res <- load_if_exists(
  file.path(
    lme_results_dir,
    "cvr_measurement_invariance.rds"
  ),
  "CVR MI"
)
lgcm_results <- load_if_exists(
  file.path(
    lgcm_results_dir,
    "lgcm_mediation_results.rds"
  ),
  "LGCM mediation"
)
lgcm_parallel <- load_if_exists(
  file.path(
    lgcm_results_dir,
    "lgcm_parallel_results.rds"
  ),
  "LGCM parallel"
)
edt_linearity <- load_if_exists(
  file.path(
    lgcm_results_dir,
    "edt_linearity_test.rds"
  ),
  "EDT linearity"
)
simulation.res <- load_if_exists(
  file.path(
    lgcm_results_dir,
    "simulation_results.rds"
  ),
  "Simulation"
)
inclusion_counts <- load_if_exists(
  file.path(
    derivatives_dir, "inclusion_counts.rds"
  ),
  "Inclusion counts"
)
nv_data <- load_if_exists(
  here(
    "data/derivatives",
    "normative_transfer_validation.rds"
  ),
  "Normative transfer validation"
)

# ---------------------------------------------------------
# 3. Adapter views (per-CVR)
# ---------------------------------------------------------
frs_hvr <- extract_cvr_view.fn(
  lme_hvr_z.res, "FRS"
)
cvr_hvr <- extract_cvr_view.fn(
  lme_hvr_z.res, "CVR_mimic"
)
frs_hc <- extract_cvr_view.fn(
  lme_hc_z.res, "FRS"
)
cvr_hc <- extract_cvr_view.fn(
  lme_hc_z.res, "CVR_mimic"
)

# ---------------------------------------------------------
# 4. LGCM mediation adapters (per-sex)
# ---------------------------------------------------------
stopifnot(
  !is.null(lgcm_results),
  !is.null(lgcm_results$results)
)
lgcm_med_male <- NULL
lgcm_med_female <- NULL
for (sx in c("Male", "Female")) {
  sx_res <- list()
  for (k in names(lgcm_results$results)) {
    parts <- strsplit(k, "_")[[1]]
    if (parts[1] != sx) next
    dom <- parts[2]
    pred <- paste(
      parts[3:length(parts)],
      collapse = "_"
    )
    short <- if (
      pred == "CVR_mimic"
    ) "CVR" else pred
    new_k <- paste(dom, short, sep = "_")
    sx_res[[new_k]] <- remap_med.fn(
      lgcm_results$results[[k]]
    )
  }
  obj <- list(results = sx_res)
  if (sx == "Male") {
    lgcm_med_male <- obj
  } else {
    lgcm_med_female <- obj
  }
}

# ---------------------------------------------------------
# 5. Cohort scalars
# ---------------------------------------------------------
stopifnot(!is.null(cohort.dt))
baseline.dt <- cohort.dt[, .SD[1], by = PTID]
n_subjects <- nrow(baseline.dt)
n_obs <- nrow(cohort.dt)
age_mean <- mean(baseline.dt$AGE, na.rm = TRUE)
age_sd <- sd(baseline.dt$AGE, na.rm = TRUE)
fu_range <- range(
  cohort.dt$YRS_from_bl, na.rm = TRUE
)
pct_male <- 100 * sum(
  baseline.dt$SEX == "Male"
) / n_subjects

# Raw FRS percentage from source CVRF file
frs_raw_col <- get_parameter("frs", "column")
cvrf_raw.dt <- fread(
  get_data_path("raw", "adsp_phc_cvrf"),
  select = c("PTID", frs_raw_col)
)
setnames(cvrf_raw.dt, frs_raw_col, "FRS_pct")
bl_frs <- merge(
  baseline.dt[, .(PTID, SEX)],
  cvrf_raw.dt, by = "PTID"
)
frs_ceiling_all <- 100 * mean(
  bl_frs$FRS_pct == 30, na.rm = TRUE
)
frs_ceiling_female <- 100 * mean(
  bl_frs[SEX == "Female"]$FRS_pct == 30,
  na.rm = TRUE
)
frs_ceiling_male <- 100 * mean(
  bl_frs[SEX == "Male"]$FRS_pct == 30,
  na.rm = TRUE
)
frs_above20_all <- 100 * mean(
  bl_frs$FRS_pct > 20, na.rm = TRUE
)

# Diagnostic composition
dx_counts <- baseline.dt[, .N, by = DX]
n_cu <- dx_counts[DX == "CU", N]
n_mci <- dx_counts[DX == "MCI", N]
n_ad <- dx_counts[DX == "AD", N]

# Visit counts
visits_per_subj <- cohort.dt[, .N, by = PTID]
median_visits <- median(visits_per_subj$N)
range_visits <- range(visits_per_subj$N)

# Age range
age_min <- min(
  baseline.dt$AGE, na.rm = TRUE
)
age_max <- max(
  baseline.dt$AGE, na.rm = TRUE
)

# Education
educ_mean <- mean(
  baseline.dt$EDUC, na.rm = TRUE
)
educ_sd <- sd(
  baseline.dt$EDUC, na.rm = TRUE
)

# APOE4
pct_apoe4 <- 100 * mean(
  baseline.dt$APOE4 > 0, na.rm = TRUE
)

# Follow-up in years
fu_mean <- mean(
  cohort.dt[, max(YRS_from_bl), by = PTID]$V1
)
fu_max <- max(
  cohort.dt$YRS_from_bl, na.rm = TRUE
)

# FRS mean/SD (raw percentage)
frs_mean_raw <- mean(
  bl_frs$FRS_pct, na.rm = TRUE
)
frs_sd_raw <- sd(
  bl_frs$FRS_pct, na.rm = TRUE
)

# CVR merge and correlations
stopifnot(!is.null(cvr_scores))
if (!"CVR_mimic" %in% names(cohort.dt)) {
  cvr_merge_cols <- intersect(
    c("RID", "PTID"), names(cvr_scores)
  )
  merge_key <- cvr_merge_cols[1]
  cohort.dt <- merge(
    cohort.dt,
    cvr_scores[
      ,
      c(merge_key, "CVR_mimic"),
      with = FALSE
    ][!duplicated(get(merge_key))],
    by = merge_key, all.x = TRUE
  )
}
merged_cvr <- merge(
  cohort.dt[
    , .(RID, AGE, FRS)
  ][, .(AGE = AGE[1], FRS = FRS[1]),
    by = RID],
  cvr_scores[, .(RID, CVR_mimic)],
  by = "RID"
)
cvr_age_cor <- cor(
  merged_cvr$CVR_mimic, merged_cvr$AGE,
  use = "complete.obs"
)
cvr_frs_cor <- cor(
  merged_cvr$CVR_mimic, merged_cvr$FRS,
  use = "complete.obs"
)
cvr_mean <- mean(
  cvr_scores$CVR_mimic, na.rm = TRUE
)
cvr_sd <- sd(
  cvr_scores$CVR_mimic, na.rm = TRUE
)

# ---------------------------------------------------------
# 6. MIMIC / MI inline stats
# ---------------------------------------------------------
mimic_fit <- cvr_model.res$fit_indices
mimic_gamma <- cvr_model.res$age_regression$gamma[1]
mimic_r2 <- mimic_gamma^2

mi_verdict <- cvr_mi.res$highest_final
stopifnot(
  !is.null(cvr_mi.res$results),
  !is.null(cvr_mi.res$results$configural)
)
mi_cfg_fit <- cvr_mi.res$results$configural
mi_results <- cvr_mi.res$results

frs_age_r <-
  frs_problem.res$age_correlation$pearson
frs_age_rho <-
  frs_problem.res$age_correlation$spearman
frs_partial_r <-
  frs_problem.res$age_correlation$partial
frs_var_decomp <-
  frs_problem.res$variance_decomposition

# Simulation convergence
stopifnot(!is.null(simulation.res$summary))
sim_conv <- simulation.res$summary

# Common sample N for standardization
stopifnot(
  !is.null(lgcm_results$sample_info)
)
common_n <-
  lgcm_results$sample_info$Male$n +
  lgcm_results$sample_info$Female$n

# ---------------------------------------------------------
# 7. FDR summary from LME results
# ---------------------------------------------------------
stopifnot(
  !is.null(lme_hvr_z.res$summary_table)
)
lme_summary.dt <- as.data.table(
  lme_hvr_z.res$summary_table
)
primary <- lme_summary.dt[
  Analysis == "Primary_Stratified"
]
fdr_summary.lst <- list(
  frs_n_sig = sum(
    primary[CVR_Measure == "FRS"]$Significant
  ),
  frs_n_fdr = sum(
    primary[
      CVR_Measure == "FRS"
    ]$Significant_FDR
  ),
  cvr_n_sig = sum(
    primary[
      CVR_Measure == "CVR_mimic"
    ]$Significant
  ),
  cvr_n_fdr = sum(
    primary[
      CVR_Measure == "CVR_mimic"
    ]$Significant_FDR
  ),
  frs_all_survive = all(
    primary[
      CVR_Measure == "FRS" &
        Significant == TRUE
    ]$Significant_FDR
  ),
  cvr_all_survive = all(
    primary[
      CVR_Measure == "CVR_mimic" &
        Significant == TRUE
    ]$Significant_FDR
  )
)

# HC z-score summary
stopifnot(
  !is.null(lme_hc_z.res$summary_table)
)
hc_summary.dt <- as.data.table(
  lme_hc_z.res$summary_table
)

# ---------------------------------------------------------
# 8. LGCM sample sizes / mediation N
# ---------------------------------------------------------
lgcm_n.lst <- lgcm_results$sample_info
r1_m <- lgcm_results$results[["Male_MEM_FRS"]]
r1_f <- lgcm_results$results[["Female_MEM_FRS"]]
stopifnot(!is.null(r1_m), !is.null(r1_f))
med_n.lst <- list(
  Male = r1_m$n, Female = r1_f$n
)

# ---------------------------------------------------------
# 9. Abstract statistics
# ---------------------------------------------------------
DOMAIN_LABELS <- c(
  MEM = "Memory", LAN = "Language",
  EXF = "Executive Function"
)

r <- lgcm_results$results
sexes <- c("Male", "Female")
doms <- c("MEM", "LAN", "EXF")

boot_sig.fn <- function(lo, hi) {
  !is.null(lo) && !is.null(hi) &&
    !is.na(lo) && !is.na(hi) &&
    (lo > 0 || hi < 0)
}
get_boot_sig <- function(
    sx, d, pred, path) {
  k <- paste(sx, d, pred, sep = "_")
  v <- r[[k]]
  stopifnot(!is.null(v))
  p <- if (path == "a") v$a_path
       else if (path == "b") v$b_path
       else stop("Unknown path: ", path)
  boot_sig.fn(
    p$boot_ci_lower, p$boot_ci_upper
  )
}

frs_a_sig_count <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      get_boot_sig(sx, d, "FRS", "a")
    })
  }
))
cvr_a_sig_count <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      get_boot_sig(
        sx, d, "CVR_mimic", "a"
      )
    })
  }
))
frs_b_sig_count <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      get_boot_sig(sx, d, "FRS", "b")
    })
  }
))
cvr_b_sig_count <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      get_boot_sig(
        sx, d, "CVR_mimic", "b"
      )
    })
  }
))
frs_ind_sig <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      k <- paste(sx, d, "FRS", sep = "_")
      boot_sig.fn(
        r[[k]]$indirect$boot_ci_lower,
        r[[k]]$indirect$boot_ci_upper
      )
    })
  }
))
cvr_ind_sig <- sum(sapply(
  sexes, function(sx) {
    sapply(doms, function(d) {
      k <- paste(
        sx, d, "CVR_mimic", sep = "_"
      )
      boot_sig.fn(
        r[[k]]$indirect$boot_ci_lower,
        r[[k]]$indirect$boot_ci_upper
      )
    })
  }
))

# Coupling significance (parallel process)
cpl_m_sig <- 0L
cpl_f_sig <- 0L
for (d in doms) {
  mk <- paste("Male", d, sep = "_")
  fk <- paste("Female", d, sep = "_")
  lr_m <- lgcm_parallel$linear_results[[mk]]
  lr_f <- lgcm_parallel$linear_results[[fk]]
  if (!is.null(lr_m$coupling) &&
      lr_m$coupling$p < 0.05) {
    cpl_m_sig <- cpl_m_sig + 1L
  }
  if (!is.null(lr_f$coupling) &&
      lr_f$coupling$p < 0.05) {
    cpl_f_sig <- cpl_f_sig + 1L
  }
}
abs_stats <- list(
  frs_a_sig = frs_a_sig_count,
  cvr_a_sig = cvr_a_sig_count,
  frs_b_sig = frs_b_sig_count,
  cvr_b_sig = cvr_b_sig_count,
  frs_ind_sig = frs_ind_sig,
  cvr_ind_sig = cvr_ind_sig,
  cpl_m_sig = cpl_m_sig,
  cpl_f_sig = cpl_f_sig,
  cpl_total_sig = cpl_m_sig + cpl_f_sig
)

# LME abstract stats
abs_lme <- {
  p <- lme_summary.dt[
    Analysis == "Primary_Stratified"
  ]
  frs_m <- p[
    CVR_Measure == "FRS" & Sex == "Male"
  ]
  frs_f <- p[
    CVR_Measure == "FRS" & Sex == "Female"
  ]
  cvr_m <- p[
    CVR_Measure == "CVR_mimic" &
      Sex == "Male"
  ]
  cvr_f <- p[
    CVR_Measure == "CVR_mimic" &
      Sex == "Female"
  ]
  list(
    frs_m_sig = sum(frs_m$Significant),
    frs_f_sig = sum(frs_f$Significant),
    cvr_m_sig = sum(cvr_m$Significant),
    cvr_f_sig = sum(cvr_f$Significant),
    cvr_f_sig_doms = paste(
      DOMAIN_LABELS[
        cvr_f[Significant == TRUE]$Domain
      ],
      collapse = ", "
    )
  )
}

# ---------------------------------------------------------
# 10. LGCM comparison summary
# ---------------------------------------------------------
lgcm_cmp <- lgcm_results$comparisons
m_da <- sapply(
  c("Male_MEM", "Male_LAN", "Male_EXF"),
  function(k) lgcm_cmp[[k]]$delta_a_est
)
f_da <- sapply(
  c("Female_MEM", "Female_LAN",
    "Female_EXF"),
  function(k) lgcm_cmp[[k]]$delta_a_est
)
m_da_mean <- mean(m_da)
f_da_mean <- mean(f_da)
da_ratio_pct <- round(
  (f_da_mean / m_da_mean - 1) * 100
)

# ---------------------------------------------------------
# 11. CONSORT counts
# ---------------------------------------------------------
cog_src <- fread(
  get_data_path("raw", "adsp_phc_cognition"),
  select = c("PTID", "PHC_Diagnosis")
)
db_ptids.v <- unique(cog_src$PTID)

amy_src <- fread(
  get_data_path("raw", "adsp_phc_amyloid"),
  select = "PTID"
)
csf_file <- get_data_path(
  "raw", "adsp_phc_biomarker"
)
csf_ptids.v <- if (file.exists(csf_file)) {
  unique(
    fread(csf_file, select = "PTID")$PTID
  )
} else {
  character(0)
}
amy_tested.v <- intersect(
  union(
    unique(amy_src$PTID), csf_ptids.v
  ),
  db_ptids.v
)

apos_ptids.v <- unique(cohort.dt$PTID)
n_with_amy <- length(amy_tested.v)

amy_pet <- fread(
  get_data_path("raw", "adsp_phc_amyloid"),
  select = c("PTID", "PHC_AMYLOID_STATUS")
)
amy_pet_pos.v <- unique(
  amy_pet[PHC_AMYLOID_STATUS == 1]$PTID
)
csf_apos.v <- if (file.exists(csf_file)) {
  csf_all <- fread(
    csf_file,
    select = c("PTID", "AT_class")
  )
  unique(
    csf_all[grepl("^A\\+", AT_class)]$PTID
  )
} else {
  character(0)
}
all_apos.v <- union(
  amy_pet_pos.v, csf_apos.v
)
all_apos_in_db.v <- intersect(
  all_apos.v, db_ptids.v
)
aneg_in_db.v <- setdiff(
  amy_tested.v, all_apos_in_db.v
)

zs_src <- as.data.table(load_if_exists(
  here("data/derivatives/adni_z-scores.rds"),
  "Z-scores (CONSORT)"
))
mri_ptids.v <- unique(zs_src$PTID)
aneg_mri.v <- intersect(
  aneg_in_db.v, mri_ptids.v
)

bl_dx <- cog_src[, .SD[1], by = PTID]
cu_ptids.v <- bl_dx[
  PHC_Diagnosis == 1
]$PTID
aneg_cu_mri.v <- intersect(
  aneg_mri.v, cu_ptids.v
)

n_age_comp <-
  nv_data$comparable_subsample$HVR$n
ukb_range <-
  nv_data$comparable_subsample$ukb_age_range
ukb_lbl <- sprintf(
  "%d\u2013%d",
  round(ukb_range[1]),
  round(ukb_range[2])
)
n_lgcm_total <-
  lgcm_n.lst$Male$n + lgcm_n.lst$Female$n

consort_counts <- list(
  n_database = inclusion_counts$n_database,
  n_with_amyloid = n_with_amy,
  n_amyloid_pos = length(all_apos_in_db.v),
  n_with_mri = inclusion_counts$n_with_mri,
  n_lme = inclusion_counts$n_analysis,
  n_lgcm = n_lgcm_total,
  n_aneg_mri = length(aneg_mri.v),
  n_aneg_cu_mri = length(aneg_cu_mri.v),
  n_age_comparable = n_age_comp,
  ukb_age_range = ukb_lbl
)

# ---------------------------------------------------------
# 12. Demographics inputs
# ---------------------------------------------------------
gamlss_path <- get_data_path(
  "external", "gamlss_models"
)
g <- load_if_exists(gamlss_path, "GAMLSS")
ukb_crs <- as.data.table(g$DATA$CRS)

ukb_subj <- ukb_crs[
  SPLIT == "train",
  .SD[1],
  by = EID
]
ukb_m <- ukb_subj[SEX == "Male"]
ukb_f <- ukb_subj[SEX == "Female"]

norm_stats <- list(
  male = list(
    n = nrow(ukb_m),
    age_mean = mean(ukb_m$AGE),
    age_sd = sd(ukb_m$AGE),
    educ_mean = mean(
      ukb_m$EDUC_num, na.rm = TRUE
    ),
    educ_sd = sd(
      ukb_m$EDUC_num, na.rm = TRUE
    )
  ),
  female = list(
    n = nrow(ukb_f),
    age_mean = mean(ukb_f$AGE),
    age_sd = sd(ukb_f$AGE),
    educ_mean = mean(
      ukb_f$EDUC_num, na.rm = TRUE
    ),
    educ_sd = sd(
      ukb_f$EDUC_num, na.rm = TRUE
    )
  )
)

vol_stats.fn <- function(dt) {
  list(
    hc_mean = mean(dt[ROI == "HC"]$VAL),
    hc_sd = sd(dt[ROI == "HC"]$VAL),
    hvr_mean = mean(dt[ROI == "HVR"]$VAL),
    hvr_sd = sd(dt[ROI == "HVR"]$VAL)
  )
}
vol_crs <- ukb_crs[
  ADJ == "NON" & SIDE == "LR" &
    SPLIT == "train" &
    ROI %in% c("HC", "HVR")
]
ukb_raw_stats <- list(
  male = vol_stats.fn(
    vol_crs[SEX == "Male"]
  ),
  female = vol_stats.fn(
    vol_crs[SEX == "Female"]
  )
)

# Build cohort_for_tbl.dt (demographics table)
adni_zs <- as.data.table(load_if_exists(
  here("data/derivatives/adni_z-scores.rds"),
  "Z-scores"
))
bl_ptids <- baseline.dt$PTID
hc_bl <- adni_zs[
  ROI == "HC" & ADJ == "NON" & SIDE == "LR" &
    PTID %in% bl_ptids
][, .SD[1], by = PTID]
hvr_bl <- adni_zs[
  ROI == "HVR" & ADJ == "NON" &
    SIDE == "LR" &
    PTID %in% bl_ptids
][, .SD[1], by = PTID]

cohort_bl <- merge(
  baseline.dt,
  hc_bl[, .(PTID, HC_vol = VAL)],
  by = "PTID", all.x = TRUE
)
cohort_bl <- merge(
  cohort_bl,
  hvr_bl[, .(PTID, HVR_vol = VAL)],
  by = "PTID", all.x = TRUE
)
frs_col <- get_parameter("frs", "column")
cvrf_file <- get_data_path(
  "raw", "adsp_phc_cvrf"
)
cvrf.dt <- fread(
  cvrf_file,
  select = c("PTID", frs_col)
)
setnames(cvrf.dt, frs_col, "FRS_pct")

cohort_for_tbl.dt <- copy(cohort.dt)
cohort_for_tbl.dt <- merge(
  cohort_for_tbl.dt,
  cohort_bl[, .(PTID, HC_vol, HVR_vol)],
  by = "PTID", all.x = TRUE
)
cohort_for_tbl.dt <- merge(
  cohort_for_tbl.dt,
  cvrf.dt, by = "PTID", all.x = TRUE
)

# ---------------------------------------------------------
# 13. Normative validation sub-objects
# ---------------------------------------------------------
ukb_demo <-
  nv_data$comparable_subsample$demographics$ukb_normative
adni_demo <-
  nv_data$comparable_subsample$demographics$adni_comparable
comp_hvr <- nv_data$comparable_subsample$HVR
comp_hc <- nv_data$comparable_subsample$HC

# ---------------------------------------------------------
# 14. Assemble and save
# ---------------------------------------------------------
env.lst <- list(
  # Full result objects (needed by narrative)
  lgcm_parallel = lgcm_parallel,
  lgcm_results = lgcm_results,
  lme_hvr_z.res = lme_hvr_z.res,
  lme_hc_z.res = lme_hc_z.res,
  simulation.res = simulation.res,
  attrition_analysis = attrition_analysis,
  lme_raw_sens.res = lme_raw_sens.res,
  frs_edt = frs_edt,
  cvr_edt = cvr_edt,
  edt_sens.res = edt_sens.res,
  edt_linearity = edt_linearity,
  frs_problem.res = frs_problem.res,
  cvr_model.res = cvr_model.res,
  cvr_mi.res = cvr_mi.res,
  inclusion_counts = inclusion_counts,
  nv_data = nv_data,

  # Adapter views
  frs_hvr = frs_hvr,
  cvr_hvr = cvr_hvr,
  frs_hc = frs_hc,
  cvr_hc = cvr_hc,

  # LGCM adapters
  lgcm_med_male = lgcm_med_male,
  lgcm_med_female = lgcm_med_female,

  # Summary tables
  lme_summary.dt = lme_summary.dt,
  hc_summary.dt = hc_summary.dt,
  fdr_summary.lst = fdr_summary.lst,

  # Cohort scalars
  n_subjects = n_subjects,
  n_obs = n_obs,
  age_mean = age_mean,
  age_sd = age_sd,
  fu_range = fu_range,
  pct_male = pct_male,
  frs_ceiling_all = frs_ceiling_all,
  frs_ceiling_female = frs_ceiling_female,
  frs_ceiling_male = frs_ceiling_male,
  frs_above20_all = frs_above20_all,
  n_cu = n_cu,
  n_mci = n_mci,
  n_ad = n_ad,
  median_visits = median_visits,
  range_visits = range_visits,
  age_min = age_min,
  age_max = age_max,
  educ_mean = educ_mean,
  educ_sd = educ_sd,
  pct_apoe4 = pct_apoe4,
  fu_mean = fu_mean,
  fu_max = fu_max,
  frs_mean_raw = frs_mean_raw,
  frs_sd_raw = frs_sd_raw,
  cvr_age_cor = cvr_age_cor,
  cvr_frs_cor = cvr_frs_cor,
  cvr_mean = cvr_mean,
  cvr_sd = cvr_sd,

  # MIMIC stats
  mimic_fit = mimic_fit,
  mimic_gamma = mimic_gamma,
  mimic_r2 = mimic_r2,
  mi_verdict = mi_verdict,
  mi_cfg_fit = mi_cfg_fit,
  mi_results = mi_results,
  frs_age_r = frs_age_r,
  frs_age_rho = frs_age_rho,
  frs_partial_r = frs_partial_r,
  frs_var_decomp = frs_var_decomp,
  sim_conv = sim_conv,
  common_n = common_n,

  # Sample info
  lgcm_n.lst = lgcm_n.lst,
  med_n.lst = med_n.lst,

  # Abstract stats
  abs_stats = abs_stats,
  abs_lme = abs_lme,

  # Comparison stats
  m_da_mean = m_da_mean,
  f_da_mean = f_da_mean,
  da_ratio_pct = da_ratio_pct,

  # CONSORT counts
  consort_counts = consort_counts,
  n_lgcm_total = n_lgcm_total,

  # Demographics inputs
  norm_stats = norm_stats,
  ukb_raw_stats = ukb_raw_stats,
  cohort_for_tbl.dt = cohort_for_tbl.dt,

  # Normative validation
  ukb_demo = ukb_demo,
  adni_demo = adni_demo,
  comp_hvr = comp_hvr,
  comp_hc = comp_hc,

  # Figure data
  cohort.dt = cohort.dt,
  cvr_scores = cvr_scores,
  baseline.dt = baseline.dt
)

write_rds_safe(
  env.lst,
  here("outputs/manuscript_env.rds"),
  "Manuscript environment"
)

cat(sprintf(
  "Saved %d objects to manuscript_env.rds\n",
  length(env.lst)
))
