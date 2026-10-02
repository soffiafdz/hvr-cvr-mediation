
# =============================================================================
# frs.R
# =============================================================================
# Framingham 10-year general cardiovascular risk, office-based (BMI) model.
#
# The ADSP-PHC variable PHC_ASCVD_10y_FRS_Simple_Ageover30 was computed with
# the opposite sex's equation for every participant
# (docs/frs_sex_swap_audit.md). This file recomputes the score from the
# ADSP-PHC inputs with the correct sex, using the coefficients, rounding and
# 1-30 clamp of CVrisk::ascvd_10y_frs_simple (CVrisk 1.2.1). As in the
# ADSP-PHC release, ages above 74 are extrapolated rather than set to NA.
#
# Reference:
#   D'Agostino RB, Vasan RS, Pencina MJ, et al. General cardiovascular risk
#   profile for use in primary care: the Framingham Heart Study.
#   Circulation. 2008;117(6). doi:10.1161/CIRCULATIONAHA.107.699579
# =============================================================================

library(data.table)

# Sex-specific coefficients (CVrisk frs_simple_coef)
FRS_COEF <- data.table(
  SEX = c("Male", "Female"),
  ln_age = c(3.11296, 2.72107),
  ln_bmi = c(0.79277, 0.51125),
  ln_sbp_untreated = c(1.85508, 2.81291),
  ln_sbp_treated = c(1.92672, 2.88267),
  smoker = c(0.70953, 0.61868),
  diabetes = c(0.53160, 0.77763),
  group_mean = c(23.9388, 26.0145),
  baseline_survival = c(0.88431, 0.94833)
)

# ADSP-PHC input columns (PHC_Sex: 1 = Male, 2 = Female)
FRS_INPUT_COLS <- c(
  "PHC_Sex", "PHC_Age_CardiovascularRisk", "PHC_BMI",
  "PHC_SBP", "PHC_BP_Med", "PHC_Smoker", "PHC_Diabetes"
)

#' Compute the Framingham 10-year risk (%)
#'
#' @param sex.v "Male" or "Female"
#' @param age.v Age in years (not truncated at 74)
#' @param bmi.v Body-mass index (kg/m^2)
#' @param sbp.v Systolic blood pressure (mmHg)
#' @param bp_med.v Antihypertensive treatment (0/1)
#' @param smoker.v Current smoker (0/1)
#' @param diabetes.v Diabetes (0/1)
#' @return Risk in percent, rounded to 2 decimals and clamped to [1, 30]
compute_frs.fn <- function(
  sex.v, age.v, bmi.v, sbp.v, bp_med.v, smoker.v, diabetes.v
) {
  if (!all(sex.v %in% c(FRS_COEF$SEX, NA))) {
    stop("sex.v must be 'Male' or 'Female'", call. = FALSE)
  }
  coef.dt <- FRS_COEF[match(sex.v, FRS_COEF$SEX)]
  ln_sbp.v <- ifelse(
    bp_med.v == 1, coef.dt$ln_sbp_treated, coef.dt$ln_sbp_untreated
  )
  sum.v <- coef.dt$ln_age * log(age.v) +
    coef.dt$ln_bmi * log(bmi.v) +
    ln_sbp.v * log(sbp.v) +
    coef.dt$smoker * smoker.v +
    coef.dt$diabetes * diabetes.v
  risk.v <- round(
    100 * (1 - coef.dt$baseline_survival^exp(sum.v - coef.dt$group_mean)),
    2
  )
  pmin(pmax(risk.v, 1), 30)
}

#' FRS for each row of an ADSP-PHC CVRF table
#'
#' "corrected" recomputes the score with the participant's own sex and
#' stops unless the opposite-sex equation reproduces the published column,
#' which guards against both coding errors here and a changed release.
#' "published" returns the released column unchanged.
#'
#' @param cvrf.dt data.table with FRS_INPUT_COLS and the published column
#' @param frs_source "corrected" or "published" (default from config)
#' @return Numeric vector of FRS (%), aligned with the rows of cvrf.dt
frs_from_cvrf.fn <- function(
  cvrf.dt, frs_source = get_parameter("frs", "source")
) {
  frs_source <- match.arg(frs_source, c("corrected", "published"))
  published.v <- cvrf.dt[[get_parameter("frs", "column")]]
  if (frs_source == "published") return(published.v)

  missing.v <- setdiff(FRS_INPUT_COLS, names(cvrf.dt))
  if (length(missing.v) > 0) {
    stop(
      "CVRF table lacks FRS inputs: ",
      paste(missing.v, collapse = ", "), call. = FALSE
    )
  }
  if (!all(cvrf.dt$PHC_Sex %in% c(1, 2, NA))) {
    stop("PHC_Sex must be coded 1 (Male) or 2 (Female)", call. = FALSE)
  }
  sex.v <- c("Male", "Female")[cvrf.dt$PHC_Sex]
  other.v <- c("Female", "Male")[cvrf.dt$PHC_Sex]
  frs.fn <- function(s.v) {
    compute_frs.fn(
      s.v, cvrf.dt$PHC_Age_CardiovascularRisk, cvrf.dt$PHC_BMI,
      cvrf.dt$PHC_SBP, cvrf.dt$PHC_BP_Med, cvrf.dt$PHC_Smoker,
      cvrf.dt$PHC_Diabetes
    )
  }
  corrected.v <- frs.fn(sex.v)
  swapped.v <- frs.fn(other.v)

  if (!identical(is.na(corrected.v), is.na(published.v))) {
    stop(
      "Corrected and published FRS differ in missingness",
      call. = FALSE
    )
  }
  ok.v <- !is.na(published.v)
  max_diff <- max(abs(swapped.v[ok.v] - published.v[ok.v]))
  if (!isTRUE(max_diff < 0.005)) {
    stop(sprintf(
      paste0(
        "Opposite-sex FRS no longer reproduces the published",
        " column (max |diff| = %.3f); recheck the ADSP-PHC release"
      ),
      max_diff
    ), call. = FALSE)
  }
  corrected.v
}

#' Load FRS per participant from the ADSP-PHC CVRF file
#'
#' @param cvrf_file Path to the CVRF CSV (default from config)
#' @param frs_source "corrected" or "published" (default from config)
#' @return data.table with RID, PTID, FRS
load_frs.fn <- function(
  cvrf_file = get_data_path("raw", "adsp_phc_cvrf"),
  frs_source = get_parameter("frs", "source")
) {
  cvrf.dt <- fread(
    cvrf_file,
    select = c(
      "RID", "PTID", get_parameter("frs", "column"), FRS_INPUT_COLS
    )
  )
  frs.v <- frs_from_cvrf.fn(cvrf.dt, frs_source)
  cvrf.dt[, .(RID, PTID, FRS = frs.v)]
}
