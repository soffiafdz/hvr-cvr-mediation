# Verify the ADSP-PHC Framingham sex swap and recompute the 30% cap.
# Run from the hvr_cvr_mediation repo root:
#   tar -xzf data/external/ADNIMERGE_0.0.1.tar.gz -C /tmp \
#     ADNIMERGE/data/adnimerge.rdata
#   Rscript docs/frs_sex_swap_check.R /tmp/ADNIMERGE/data/adnimerge.rdata
# See docs/frs_sex_swap_audit.md.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
adnimerge_rdata <- args[1]

# D'Agostino 2008 office-based (BMI) model, as in CVrisk::frs_simple_coef
COEF.dt <- data.table(
  sex_lab = c("Male", "Female"),
  ln_age = c(3.11296, 2.72107), ln_bmi = c(0.79277, 0.51125),
  ln_sbp_u = c(1.85508, 2.81291), ln_sbp_t = c(1.92672, 2.88267),
  smoker = c(0.70953, 0.61868), diabetes = c(0.53160, 0.77763),
  grp_mean = c(23.9388, 26.0145), s0 = c(0.88431, 0.94833)
)
# Same arithmetic, rounding and 1-30 clamp as CVrisk::ascvd_10y_frs_simple
# (but without its age > 74 -> NA rule, which ADSP-PHC evidently dropped)
frs.fn <- function(sex.v, age, bmi, sbp, bp_med, smk, dm) {
  idx.v <- match(sex.v, COEF.dt$sex_lab)
  k.dt <- COEF.dt[idx.v]
  x.v <- k.dt$ln_age * log(age) + k.dt$ln_bmi * log(bmi) +
    ifelse(bp_med == 1, k.dt$ln_sbp_t, k.dt$ln_sbp_u) * log(sbp) +
    k.dt$smoker * smk + k.dt$diabetes * dm
  r.v <- round((1 - k.dt$s0^exp(x.v - k.dt$grp_mean)) * 100, 2)
  pmin(pmax(r.v, 1), 30)
}

cv.dt <- fread("data/raw/adsp_phc/ADSP_PHC_CVRF_11Dec2025.csv")
cv.dt[, sex_lab := ifelse(PHC_Sex == 1, "Male", "Female")]
cv.dt[, other_lab := ifelse(PHC_Sex == 1, "Female", "Male")]
cv.dt[, FRS_correct := frs.fn(sex_lab, PHC_Age_CardiovascularRisk,
  PHC_BMI, PHC_SBP, PHC_BP_Med, PHC_Smoker, PHC_Diabetes)]
cv.dt[, FRS_other := frs.fn(other_lab, PHC_Age_CardiovascularRisk,
  PHC_BMI, PHC_SBP, PHC_BP_Med, PHC_Smoker, PHC_Diabetes)]
PHC <- "PHC_ASCVD_10y_FRS_Simple_Ageover30"
ok.dt <- cv.dt[!is.na(get(PHC)) & !is.na(FRS_correct)]
cat("Rows with PHC FRS and complete inputs:", nrow(ok.dt), "\n")
for (sx in c("Male", "Female")) {
  s.dt <- ok.dt[sex_lab == sx]
  cat(sprintf(
    paste0(
      "%s (PHC_Sex=%s, n=%d): max|PHC - own-sex eq| = %.3f;",
      " max|PHC - other-sex eq| = %.3f;",
      " share matching other-sex eq (<=0.01) = %.3f\n"
    ),
    sx, ifelse(sx == "Male", 1, 2), nrow(s.dt),
    max(abs(s.dt[[PHC]] - s.dt$FRS_correct)),
    max(abs(s.dt[[PHC]] - s.dt$FRS_other)),
    mean(abs(s.dt[[PHC]] - s.dt$FRS_other) <= 0.01)))
}

# Sex label check against ADNIMERGE
load(adnimerge_rdata)
am.dt <- unique(as.data.table(adnimerge)[, .(RID, PTGENDER)])
x.dt <- merge(unique(cv.dt[, .(RID, PHC_Sex)]), am.dt, by = "RID")
cat("\nPHC_Sex vs ADNIMERGE PTGENDER:\n")
print(table(PHC_Sex = x.dt$PHC_Sex, PTGENDER = x.dt$PTGENDER))

# Analysis cohort (as in script 15: baseline row per PTID of cohort_hvr)
co.dt <- as.data.table(readRDS("data/derivatives/lme/cohort_hvr.rds"))
bl.dt <- co.dt[, .SD[1], by = PTID]
y.dt <- merge(bl.dt[, .(PTID, SEX)],
  unique(as.data.table(adnimerge)[, .(PTID, PTGENDER)]), by = "PTID")
cat("\nPipeline SEX vs ADNIMERGE PTGENDER (cohort baseline):\n")
print(table(SEX = y.dt$SEX, PTGENDER = y.dt$PTGENDER))
m.dt <- merge(bl.dt[, .(PTID, SEX)],
  cv.dt[, .(PTID, PHC = get(PHC), FRS_correct)], by = "PTID")
pct.fn <- function(v) sprintf("%.1f", 100 * mean(v == 30, na.rm = TRUE))
cat(sprintf("\nAnalysis cohort baseline: n=%d (men %d, women %d)\n",
  nrow(m.dt), m.dt[SEX == "Male", .N], m.dt[SEX == "Female", .N]))
cat("  % at cap, PHC as published: all", pct.fn(m.dt$PHC),
  "| men", pct.fn(m.dt[SEX == "Male"]$PHC),
  "| women", pct.fn(m.dt[SEX == "Female"]$PHC), "\n")
cat("  % at cap, correct-sex eq:   all", pct.fn(m.dt$FRS_correct),
  "| men", pct.fn(m.dt[SEX == "Male"]$FRS_correct),
  "| women", pct.fn(m.dt[SEX == "Female"]$FRS_correct), "\n")
cat(sprintf("  Mean FRS, PHC as published: men %.1f, women %.1f\n",
  mean(m.dt[SEX == "Male"]$PHC, na.rm = TRUE),
  mean(m.dt[SEX == "Female"]$PHC, na.rm = TRUE)))
cat(sprintf("  Mean FRS, correct-sex eq:   men %.1f, women %.1f\n",
  mean(m.dt[SEX == "Male"]$FRS_correct, na.rm = TRUE),
  mean(m.dt[SEX == "Female"]$FRS_correct, na.rm = TRUE)))
cat(sprintf("  Spearman r(PHC, correct) within men %.3f, within women %.3f\n",
  cor(m.dt[SEX == "Male"]$PHC, m.dt[SEX == "Male"]$FRS_correct,
    method = "spearman", use = "complete.obs"),
  cor(m.dt[SEX == "Female"]$PHC, m.dt[SEX == "Female"]$FRS_correct,
    method = "spearman", use = "complete.obs")))
cat("\nRows aged > 74 with a PHC value (stock CVrisk returns NA):",
  cv.dt[PHC_Age_CardiovascularRisk > 74 & !is.na(get(PHC)), .N], "\n")
