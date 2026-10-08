# Compare the published and corrected FRS against CVR_mimic.
# Run from the hvr_cvr_mediation repo root after the corrected re-run:
#   Rscript docs/frs_cvr_comparison.R
# "Old" results come from archive/2026-10-02_pre_frs_fix/ (published
# ADSP-PHC FRS); "new" results from the current outputs (corrected FRS).
# See docs/frs_sex_swap_audit.md.
suppressMessages(library(data.table))
ARCHIVE <- "archive/2026-10-02_pre_frs_fix"
fmt.fn <- function(x) sprintf("%.3f", x)

# --- 1. Score agreement, analysis cohort baseline ---
# The archived cohort's FRS_raw holds z-scores (script 07 was once re-run
# without 03), so take the old FRS percent from the archived score table.
old_sc.dt <- as.data.table(readRDS(file.path(
  ARCHIVE, "data_derivatives/lme/cvr_mimic_scores.rds"
)))
new_sc.dt <- as.data.table(readRDS(
  "data/derivatives/lme/cvr_mimic_scores.rds"
))
stopifnot(
  identical(old_sc.dt$RID, new_sc.dt$RID),
  isTRUE(all.equal(old_sc.dt$CVR_mimic, new_sc.dt$CVR_mimic))
)
coh.dt <- as.data.table(readRDS("data/derivatives/lme/cohort_hvr.rds"))
coh.dt <- coh.dt[order(PTID, YRS_from_bl)][, .SD[1], by = PTID]
coh.dt <- coh.dt[, .(
  RID, SEX, AGE, FRS_new = FRS_raw, CVR = CVR_mimic
)]
coh.dt <- merge(coh.dt, old_sc.dt[, .(RID, FRS_old = FRS)], by = "RID")
chk.dt <- merge(
  coh.dt[, .(RID, FRS_new)], new_sc.dt[, .(RID, FRS)], by = "RID"
)
stopifnot(isTRUE(all.equal(chk.dt$FRS_new, chk.dt$FRS)))
cat("Analysis cohort baseline n =", nrow(coh.dt), "\n\n")

agree.fn <- function(d.dt, lab) {
  res.lst <- lapply(c("FRS_old", "FRS_new"), function(v) {
    resid.v <- resid(lm(d.dt[[v]] ~ d.dt$AGE))
    c(
      r_cvr = cor(d.dt$CVR, d.dt[[v]]),
      rho_cvr = cor(d.dt$CVR, d.dt[[v]], method = "spearman"),
      r_age = cor(d.dt$AGE, d.dt[[v]]),
      r_cvr_agefree = cor(d.dt$CVR, resid.v)
    )
  })
  pair.fn <- function(k) {
    paste(fmt.fn(res.lst[[1]][k]), "->", fmt.fn(res.lst[[2]][k]))
  }
  cat(sprintf(
    paste0(
      "%-7s | r(CVR,FRS) %s | rho %s | r(age,FRS) %s",
      " | r(CVR, age-free FRS) %s | r(age,CVR) %s\n"
    ),
    lab, pair.fn("r_cvr"), pair.fn("rho_cvr"), pair.fn("r_age"),
    pair.fn("r_cvr_agefree"), fmt.fn(cor(d.dt$AGE, d.dt$CVR))
  ))
}
cat("== 1. Score agreement (old -> corrected)\n")
agree.fn(coh.dt, "All")
agree.fn(coh.dt[SEX == "Male"], "Male")
agree.fn(coh.dt[SEX == "Female"], "Female")

# Between- vs within-sex: remove each sex's mean from both scores
sex_means.dt <- coh.dt[, .(
  n = .N, CVR_mean = round(mean(CVR), 3),
  FRS_old_mean = round(mean(FRS_old), 1),
  FRS_new_mean = round(mean(FRS_new), 1)
), by = SEX]
print(sex_means.dt)
coh.dt[, `:=`(
  CVR_w = CVR - mean(CVR),
  FRS_old_w = FRS_old - mean(FRS_old),
  FRS_new_w = FRS_new - mean(FRS_new)
), by = SEX]
cat(sprintf(
  "Within-sex pooled r(CVR,FRS): old %.3f -> new %.3f\n",
  cor(coh.dt$CVR_w, coh.dt$FRS_old_w), cor(coh.dt$CVR_w, coh.dt$FRS_new_w)
))

# --- 2. Script 04: how much of FRS is age ---
cat("\n== 2. FRS variance decomposition (script 04; old -> corrected)\n")
old_frs.res <- readRDS(file.path(
  ARCHIVE, "models_results/lme/frs_problem_analysis.rds"
))
new_frs.res <- readRDS("models/results/lme/frs_problem_analysis.rds")
for (k in c("r2_age", "r2_cvrf", "r2_full")) {
  cat(sprintf(
    "  %-8s %s -> %s\n", k,
    fmt.fn(old_frs.res$variance_decomposition[[k]]),
    fmt.fn(new_frs.res$variance_decomposition[[k]])
  ))
}

# --- 3. LME 3-way betas and paired bootstrap dB (FRS - CVR_mimic) ---
cat("\n== 3. LME YRS x CVR x HVR_z (old -> corrected)\n")
old_lme.res <- readRDS(file.path(
  ARCHIVE, "models_results/lme/lme_hvr_z_results.rds"
))
new_lme.res <- readRDS("models/results/lme/lme_hvr_z_results.rds")
beta_tab.fn <- function(lme.res) {
  s.dt <- lme.res$summary_table[Analysis == "Primary_Stratified"]
  dcast(
    s.dt[CVR_Measure %in% c("FRS", "CVR_mimic")],
    Sex + Domain ~ CVR_Measure, value.var = "Beta"
  )
}
old_tab.dt <- beta_tab.fn(old_lme.res)
new_tab.dt <- beta_tab.fn(new_lme.res)
ci.fn <- function(cmp.lst) {
  sprintf(
    "%7.4f [%7.4f, %7.4f]",
    cmp.lst$delta_beta, cmp.lst$ci_lower, cmp.lst$ci_upper
  )
}
for (i in seq_len(nrow(new_tab.dt))) {
  key <- paste(new_tab.dt$Sex[i], new_tab.dt$Domain[i], sep = "_")
  old_cmp.lst <- old_lme.res$comparison[[key]]
  new_cmp.lst <- new_lme.res$comparison[[key]]
  cat(sprintf(
    paste0(
      "  %-10s FRS b %7.4f -> %7.4f | CVR b %7.4f",
      " | dB %s -> %s | sig %s -> %s\n"
    ),
    key, old_tab.dt$FRS[i], new_tab.dt$FRS[i], new_tab.dt$CVR_mimic[i],
    ci.fn(old_cmp.lst), ci.fn(new_cmp.lst),
    old_cmp.lst$significant, new_cmp.lst$significant
  ))
}
ratio.fn <- function(tab.dt) {
  mean(abs(tab.dt$FRS)) / mean(abs(tab.dt$CVR_mimic))
}
cat(sprintf(
  "  beta_ratio (mean|FRS b| / mean|CVR b|): %.1f -> %.1f\n",
  ratio.fn(old_tab.dt), ratio.fn(new_tab.dt)
))

# --- 4. LGCM mediation a-paths, FRS models finished so far ---
cat("\n== 4. LGCM mediation, FRS models finished so far\n")
CHK_DIR <- "models/results/lgcm/mediation_checkpoints"
path.fn <- function(a.lst) {
  sprintf(
    "%.4f [%.4f, %.4f]",
    a.lst$est, a.lst$boot_ci_lower, a.lst$boot_ci_upper
  )
}
for (f in list.files(CHK_DIR, pattern = "_FRS\\.rds$")) {
  key <- sub("_FRS\\.rds$", "", f)
  new_med.res <- readRDS(file.path(CHK_DIR, f))
  old_med.res <- readRDS(file.path(
    ARCHIVE, "models_results/lgcm/mediation_checkpoints", f
  ))
  cvr_med.res <- readRDS(file.path(CHK_DIR, paste0(key, "_CVR_mimic.rds")))
  cat(sprintf(
    "  %s a-path: FRS old %s | FRS new %s | CVR %s\n", key,
    path.fn(old_med.res$a_path), path.fn(new_med.res$a_path),
    path.fn(cvr_med.res$a_path)
  ))
}
