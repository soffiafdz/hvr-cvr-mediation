#!/usr/bin/env Rscript
# =========================================================
# Export manuscript figures for journal submission
# Generates standalone TIFF/PDF files in
#   outputs/submission_ad/figures/
# =========================================================

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/plotting.R"))
source(here("R/utils/consort_diagram.R"))

config <- load_config()

# ---------------------------------------------------------
# Load pre-computed manuscript environment
# ---------------------------------------------------------
env.lst <- readRDS(here("outputs/manuscript_env.rds"))
list2env(env.lst, envir = environment())
rm(env.lst)

out_dir <- here("outputs/submission_ad/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Exporting submission figures ===\n")

# ---------------------------------------------------------
# 1. fig-consort (PDF via standalone TikZ)
# ---------------------------------------------------------
cat("  [1/6] fig-consort ...\n")

cc <- consort_counts
tikz_code <- generate_consort_tikz(
  n_database = cc$n_database,
  n_with_amyloid = cc$n_with_amyloid,
  n_amyloid_pos = cc$n_amyloid_pos,
  n_with_mri = cc$n_with_mri,
  n_lme = cc$n_lme,
  n_lgcm = cc$n_lgcm,
  n_aneg_mri = cc$n_aneg_mri,
  n_aneg_cu_mri = cc$n_aneg_cu_mri,
  n_age_comparable = cc$n_age_comparable,
  ukb_age_range = cc$ukb_age_range
)

# Write standalone LaTeX document
consort_tex <- file.path(out_dir, "fig-consort.tex")
consort_pdf <- file.path(out_dir, "fig-consort.pdf")

tex_lines <- c(
  "\\documentclass[tikz,border=5pt]{standalone}",
  "\\usepackage[T1]{fontenc}",
  "\\usetikzlibrary{",
  "  shapes.geometric, arrows.meta, positioning",
  "}",
  "\\begin{document}",
  tikz_code,
  "\\end{document}"
)
writeLines(tex_lines, consort_tex)

# Compile with xelatex
system2(
  "xelatex",
  args = c(
    "-interaction=nonstopmode",
    paste0("-output-directory=", out_dir),
    consort_tex
  ),
  stdout = FALSE, stderr = FALSE
)

# Clean auxiliary files
aux_exts <- c(".aux", ".log")
for (ext in aux_exts) {
  f <- file.path(out_dir, paste0("fig-consort", ext))
  if (file.exists(f)) file.remove(f)
}

if (file.exists(consort_pdf)) {
  cat("    -> fig-consort.pdf OK\n")
} else {
  warning("fig-consort.pdf was not generated!")
}

# ---------------------------------------------------------
# 2. fig-cvr-validation (TIFF 500 DPI, 6.5 x 4.5)
# ---------------------------------------------------------
cat("  [2/6] fig-cvr-validation ...\n")

stopifnot(!is.null(cvr_scores), !is.null(cohort.dt))

frs_raw_col <- get_parameter("frs", "column")
cvrf_raw.dt <- fread(
  get_data_path("raw", "adsp_phc_cvrf"),
  select = c("PTID", frs_raw_col)
)
setnames(cvrf_raw.dt, frs_raw_col, "FRS_pct")

baseline_merged.dt <- cohort.dt[, .(
  AGE = AGE[1], SEX = SEX[1]
), by = .(RID, PTID)]
baseline_merged.dt <- merge(
  baseline_merged.dt,
  cvrf_raw.dt, by = "PTID"
)
merged_plot.dt <- merge(
  baseline_merged.dt,
  cvr_scores[, .(RID, CVR_mimic)],
  by = "RID"
)

p_frs_sex <- plot_frs_ceiling_by_sex(merged_plot.dt)
p_frs_age <- plot_frs_validation(
  merged_plot.dt, frs_age_r = frs_age_r
)
p_cvr_age <- plot_cvr_validation(merged_plot.dt)

p_cvr_val <- p_frs_sex /
  (p_frs_age + p_cvr_age) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(out_dir, "fig-cvr-validation.tiff"),
  plot = p_cvr_val,
  width = 6.5, height = 4.5, units = "in",
  dpi = 500, device = "tiff", compression = "lzw"
)
cat("    -> fig-cvr-validation.tiff OK\n")

# ---------------------------------------------------------
# 3. fig-lme-forest (TIFF 500 DPI, 6.5 x 2.5)
# ---------------------------------------------------------
cat("  [3/6] fig-lme-forest ...\n")

p_forest <- plot_lme_forest_comparison(frs_hvr, cvr_hvr)

ggsave(
  file.path(out_dir, "fig-lme-forest.tiff"),
  plot = p_forest,
  width = 6.5, height = 2.5, units = "in",
  dpi = 500, device = "tiff", compression = "lzw"
)
cat("    -> fig-lme-forest.tiff OK\n")

# ---------------------------------------------------------
# 4. fig-coupling (TIFF 500 DPI, 6.5 x 7.5)
# ---------------------------------------------------------
cat("  [4/6] fig-coupling ...\n")

lgcm_m_input <- read_rds_safe(here(
  "data/derivatives/lgcm/lgcm_input_male.rds"
))
lgcm_f_input <- read_rds_safe(here(
  "data/derivatives/lgcm/lgcm_input_female.rds"
))
stopifnot(
  !is.null(lgcm_m_input),
  !is.null(lgcm_f_input),
  !is.null(lgcm_parallel)
)

p_coupling <- plot_coupling_figure(
  lgcm_parallel,
  lgcm_m_input, lgcm_f_input
)

ggsave(
  file.path(out_dir, "fig-coupling.tiff"),
  plot = p_coupling,
  width = 6.5, height = 7.5, units = "in",
  dpi = 500, device = "tiff", compression = "lzw"
)
cat("    -> fig-coupling.tiff OK\n")

# ---------------------------------------------------------
# 5. fig-mediation (TIFF 500 DPI, 6.5 x 8.0)
#    Original uses dev: png — try TIFF first, fall back
#    to PNG -> TIFF conversion if cairo issues arise.
# ---------------------------------------------------------
cat("  [5/6] fig-mediation ...\n")

panels <- build_combined_med_panels(
  lgcm_med_male, lgcm_med_female
)
p_med <- wrap_plots(panels, ncol = 2) +
  plot_annotation(tag_levels = "A")

med_tiff <- file.path(
  out_dir, "fig-mediation.tiff"
)
med_ok <- tryCatch({
  ggsave(
    med_tiff,
    plot = p_med,
    width = 6.5, height = 8.0, units = "in",
    dpi = 500, device = "tiff", compression = "lzw"
  )
  TRUE
}, error = function(e) {
  message("  TIFF failed: ", e$message)
  FALSE
})

if (!med_ok) {
  cat("    Falling back to PNG -> TIFF ...\n")
  med_png <- file.path(
    out_dir, "fig-mediation.png"
  )
  ggsave(
    med_png,
    plot = p_med,
    width = 6.5, height = 8.0, units = "in",
    dpi = 500, device = "png"
  )
  # Convert PNG to TIFF via ImageMagick
  system2(
    "convert",
    args = c(med_png, "-compress", "LZW", med_tiff),
    stdout = FALSE, stderr = FALSE
  )
  if (file.exists(med_tiff)) file.remove(med_png)
}
cat("    -> fig-mediation.tiff OK\n")

# ---------------------------------------------------------
# 6. fig-normative-validation (TIFF 500 DPI, 6.5 x 6.0)
#    (Supplementary figure A1)
# ---------------------------------------------------------
cat("  [6/6] fig-normative-validation ...\n")

stopifnot(!is.null(nv_data))
ukb_d <- nv_data$comparable_subsample$demographics$
  ukb_normative
adni_d <- nv_data$comparable_subsample$demographics$
  adni_comparable

p_age <- plot_age_distribution_comparison(
  ukb_age_mean = ukb_d$age_mean,
  ukb_age_sd = ukb_d$age_sd,
  ukb_age_min = ukb_d$age_min,
  ukb_age_max = ukb_d$age_max,
  ukb_n = ukb_d$n_total,
  adni_age_mean = adni_d$age_mean,
  adni_age_sd = adni_d$age_sd,
  adni_n = adni_d$n_total
)

p_zforest <- plot_zscore_validation_forest(nv_data)

hvr_z.v <- comp_hvr$zscores
hc_z.v <- comp_hc$zscores
stopifnot(!is.null(hvr_z.v), !is.null(hc_z.v))
n_max <- max(length(hvr_z.v), length(hc_z.v))
zscore.dt <- data.table(
  PTID = seq_len(n_max),
  HVR_z = hvr_z.v[seq_len(n_max)],
  HC_z = hc_z.v[seq_len(n_max)]
)

p_dist <- plot_zscore_distribution(
  zscore.dt, sample_label = "Comparable"
)

p_normval <- (p_age | p_zforest) / p_dist +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(
    out_dir, "fig-normative-validation.tiff"
  ),
  plot = p_normval,
  width = 6.5, height = 6.0, units = "in",
  dpi = 500, device = "tiff", compression = "lzw"
)
cat("    -> fig-normative-validation.tiff OK\n")

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------
exported <- list.files(out_dir)
cat("\n=== Export complete ===\n")
cat(sprintf(
  "  %d files in %s:\n", length(exported), out_dir
))
for (f in exported) cat(sprintf("    %s\n", f))
