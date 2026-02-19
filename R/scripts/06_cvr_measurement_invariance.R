#!/usr/bin/env Rscript

# ============================================================
# 06_cvr_measurement_invariance.R
# ============================================================
# PURPOSE:
#   Test measurement invariance of the CVR latent
#   factor across sex. Required before interpreting
#   sex-stratified LME and LGCM results.
#
# APPROACH:
#   Multi-group CFA on the measurement model only
#   (CVR =~ indicators). The Age→CVR path is
#   structural, not measurement — MI tests whether
#   indicators relate to CVR the same way in males
#   and females.
#
# ESTIMATOR:
#   MLR with FIML. The 5-indicator model has one
#   binary indicator (HTN), treated as continuous
#   — an accepted approximation when prevalence is
#   not extreme (~46%). This allows FIML for missing
#   lab values.
#
# INVARIANCE SEQUENCE:
#   Configural: Same structure, all parameters free
#   Metric:     Factor loadings constrained equal
#   Scalar:     Loadings + intercepts constrained
#   Strict:     All above + residual variances
#
# EVALUATION:
#   1. Scaled chi-square difference (lavTestLRT)
#   2. Delta fit: |dCFI| <= 0.01, dRMSEA <= 0.015
#
# INPUTS:
#   - data/raw/adsp_phc/ADSP_PHC_CVRF_*.csv
#   - data/raw/All_Subjects_LABDATA_*.csv
#   - models/results/lme/cvr_mimic_specification.rds
#
# OUTPUTS:
#   - models/results/lme/cvr_measurement_invariance.rds
#
# REFERENCES:
#   Vandenberg & Lance (2000). Measurement invariance.
#   Cheung & Rensvold (2002). MI fit indices.
#   Putnick & Bornstein (2016). MI conventions.
# ============================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(lavaan)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))

log_script_start("06_cvr_measurement_invariance.R")
config <- load_config()
validate_config(config)


# ============================================================
# PART 1: PREPARE DATA
# ============================================================
log_section("Part 1: Prepare Data")

# Load specification from script 05
spec_path <- get_data_path(
  "models", "cvr_mimic_specification"
)
check_files_exist(spec_path)
spec.lst <- read_rds_safe(
  spec_path, "MIMIC specification"
)

# Load and merge data (same as script 05)
cvrf_path <- get_data_path("raw", "adsp_phc_cvrf")
cvrf.dt <- fread(cvrf_path)
cvrf.dt <- cvrf.dt[order(RID, EXAMDATE)]
cvrf.dt <- cvrf.dt[, .SD[1], by = RID]

lab_path <- get_data_path("raw", "adni_labdata")
lab.dt <- fread(lab_path)
for (v in c("RCT11", "RCT20", "RCT392")) {
  lab.dt[, (v) := suppressWarnings(
    as.numeric(get(v))
  )]
}
lab.dt <- lab.dt[order(RID, EXAMDATE)]
lab_bl.dt <- lab.dt[, .(
  GLUCOSE = RCT11[1],
  CHOLESTEROL = RCT20[1],
  CREATININE = RCT392[1]
), by = RID]

analysis.dt <- merge(
  cvrf.dt, lab_bl.dt, by = "RID", all.x = TRUE
)

# Prepare model data with sex
mi.df <- data.frame(
  SEX = ifelse(
    analysis.dt$PHC_Sex == 1, "Male", "Female"
  ),
  SBP_z = as.numeric(scale(analysis.dt$PHC_SBP)),
  HTN = as.numeric(analysis.dt$PHC_Hypertension),
  GLUCOSE_z = as.numeric(
    scale(analysis.dt$GLUCOSE)
  ),
  CHOL_z = as.numeric(
    scale(analysis.dt$CHOLESTEROL)
  ),
  CREAT_z = as.numeric(
    scale(analysis.dt$CREATININE)
  )
)

n_males <- sum(mi.df$SEX == "Male", na.rm = TRUE)
n_females <- sum(
  mi.df$SEX == "Female", na.rm = TRUE
)
log_info(
  "Sample: Males=%d, Females=%d",
  n_males, n_females
)

# CFA measurement model (no Age regression)
# Same residual covariances as the MIMIC
CFA_SYNTAX <- "
  CVR =~ SBP_z + HTN + GLUCOSE_z +
         CHOL_z + CREAT_z
  SBP_z ~~ CHOL_z
  SBP_z ~~ CREAT_z
"

# ============================================================
# PART 2: FIT INVARIANCE SEQUENCE
# ============================================================
log_section("Part 2: Invariance Testing")

MI_LEVELS <- c(
  "configural", "metric", "scalar", "strict"
)
CONSTRAINTS <- list(
  configural = NULL,
  metric = "loadings",
  scalar = c("loadings", "intercepts"),
  strict = c("loadings", "intercepts", "residuals")
)

mi_fits.lst <- list()
mi_results.lst <- list()

for (level in MI_LEVELS) {
  log_info("Fitting %s...", level)
  fit <- tryCatch(
    cfa(
      CFA_SYNTAX,
      data = mi.df,
      group = "SEX",
      estimator = "MLR",
      missing = "fiml",
      group.equal = CONSTRAINTS[[level]]
    ),
    error = function(e) {
      log_warn("  Error: %s", e$message)
      NULL
    }
  )

  converged <- !is.null(fit) && tryCatch(
    lavInspect(fit, "converged"),
    error = function(e) FALSE
  )

  if (converged) {
    fi <- fitMeasures(fit, c(
      "chisq.scaled", "df.scaled",
      "pvalue.scaled",
      "cfi.robust", "tli.robust",
      "rmsea.robust", "srmr"
    ))
    mi_fits.lst[[level]] <- fit
    mi_results.lst[[level]] <- list(
      level = level, converged = TRUE,
      chisq = fi["chisq.scaled"],
      df = fi["df.scaled"],
      pvalue = fi["pvalue.scaled"],
      CFI = fi["cfi.robust"],
      TLI = fi["tli.robust"],
      RMSEA = fi["rmsea.robust"],
      SRMR = fi["srmr"]
    )
    log_info(
      "  CFI=%.4f TLI=%.4f RMSEA=%.4f SRMR=%.4f",
      fi["cfi.robust"], fi["tli.robust"],
      fi["rmsea.robust"], fi["srmr"]
    )
  } else {
    mi_results.lst[[level]] <- list(
      level = level, converged = FALSE
    )
    log_warn("  Did not converge")
  }
}

# ============================================================
# PART 3: COMPARISONS
# ============================================================
log_section("Part 3: Comparisons")

mi_pairs <- list(
  c("configural", "metric"),
  c("metric", "scalar"),
  c("scalar", "strict")
)

# --- Scaled chi-square difference tests ---
log_info("Scaled DIFFTEST:")
mi_comp.lst <- list()
highest_lrt <- "configural"

for (pair in mi_pairs) {
  base_r <- mi_results.lst[[pair[1]]]
  test_r <- mi_results.lst[[pair[2]]]
  if (!base_r$converged || !test_r$converged) {
    log_warn(
      "  %s -> %s: skipped", pair[1], pair[2]
    )
    next
  }
  diff <- tryCatch(
    lavTestLRT(
      mi_fits.lst[[pair[1]]],
      mi_fits.lst[[pair[2]]]
    ),
    error = function(e) NULL
  )
  if (!is.null(diff) && nrow(diff) == 2) {
    d_chi <- diff[2, "Chisq diff"]
    d_df <- diff[2, "Df diff"]
    p <- diff[2, "Pr(>Chisq)"]
    dec <- ifelse(p > 0.05, "HOLD", "REJECT")
    key <- paste(pair, collapse = "_vs_")
    mi_comp.lst[[key]] <- list(
      base = pair[1], test = pair[2],
      d_chi2 = d_chi, d_df = d_df,
      p = p, decision = dec
    )
    log_info(
      "  %s->%s: dChi=%.2f ddf=%.0f p=%.4f %s",
      pair[1], pair[2], d_chi, d_df, p, dec
    )
    if (dec == "HOLD") highest_lrt <- pair[2]
    if (dec == "REJECT") break
  }
}

# --- Delta fit indices ---
log_info("")
log_info("Delta Fit Indices:")
delta_comp.lst <- list()
highest_delta <- "configural"

for (pair in mi_pairs) {
  base_r <- mi_results.lst[[pair[1]]]
  test_r <- mi_results.lst[[pair[2]]]
  if (!base_r$converged || !test_r$converged) next

  d_cfi <- base_r$CFI - test_r$CFI
  d_rmsea <- test_r$RMSEA - base_r$RMSEA
  cfi_ok <- abs(d_cfi) <= 0.01
  rmsea_ok <- abs(d_rmsea) <= 0.015
  dec <- if (cfi_ok && rmsea_ok) "HOLD" else
    "REJECT"

  key <- paste(pair, collapse = "_vs_")
  delta_comp.lst[[key]] <- list(
    base = pair[1], test = pair[2],
    d_cfi = d_cfi, d_rmsea = d_rmsea,
    cfi_ok = cfi_ok, rmsea_ok = rmsea_ok,
    decision = dec
  )
  log_info(
    "  %s->%s: dCFI=%.4f dRMSEA=%.4f %s",
    pair[1], pair[2], d_cfi, d_rmsea, dec
  )
  if (dec == "HOLD") highest_delta <- pair[2]
  if (dec == "REJECT") break
}

# --- Final verdict ---
final_level <- ifelse(
  match(highest_delta, MI_LEVELS) >=
    match(highest_lrt, MI_LEVELS),
  highest_delta, highest_lrt
)

log_info("")
log_info("VERDICT: %s invariance", toupper(final_level))
if (final_level %in% c("scalar", "strict")) {
  log_info("  Scores comparable across sexes")
  log_info("  Sex-stratified analyses VALID")
} else if (final_level == "metric") {
  log_info("  Loadings equal; within-sex valid")
  log_warn("  Do not compare means across sexes")
} else {
  log_warn("  Only configural — interpret cautiously")
}

# ============================================================
# PART 4: SAVE
# ============================================================
log_section("Part 4: Save")

mi_output.lst <- list(
  model = "CFA measurement model for CVR",
  estimator = "MLR + FIML",
  indicators = spec.lst$indicator_z_names,
  sample = list(
    males = n_males, females = n_females
  ),
  results = mi_results.lst,
  comparisons_lrt = mi_comp.lst,
  comparisons_delta = delta_comp.lst,
  highest_lrt = highest_lrt,
  highest_delta = highest_delta,
  highest_final = final_level,
  timestamp = Sys.time()
)

mi_path <- get_data_path(
  "models", "cvr_measurement_invariance"
)
write_rds_safe(
  mi_output.lst, mi_path,
  "CVR measurement invariance"
)

log_info("")
log_info(paste(rep("=", 70), collapse = ""))
log_info("06_cvr_measurement_invariance.R COMPLETE")
log_info(paste(rep("=", 70), collapse = ""))
