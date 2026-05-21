#!/usr/bin/env Rscript

# Diff new (heteroscedastic HVR) vs archived (homoscedastic HVR)
# mediation results. Read-only against both.
suppressPackageStartupMessages({
  library(data.table)
})

OLD <- readRDS(here::here(
  "archive/2026-05-19_initial_submission/models_results/lgcm",
  "lgcm_mediation_results.rds"
))
NEW <- readRDS(here::here(
  "models/results/lgcm/lgcm_mediation_results.rds"
))

get_path <- function(res, slot, field) {
  x <- res[[slot]]
  if (is.null(x)) return(NA_real_)
  v <- x[[field]]
  if (is.null(v)) return(NA_real_)
  as.numeric(v)
}

build_row <- function(key, old_res, new_res) {
  data.table(
    cell = key,
    a_old = get_path(old_res, "a_path", "est"),
    a_new = get_path(new_res, "a_path", "est"),
    b_old = get_path(old_res, "b_path", "est"),
    b_new = get_path(new_res, "b_path", "est"),
    cprime_old = get_path(old_res, "cprime", "est"),
    cprime_new = get_path(new_res, "cprime", "est"),
    ind_old = get_path(old_res, "indirect", "est"),
    ind_new = get_path(new_res, "indirect", "est"),
    ind_lo_old = get_path(old_res, "indirect", "boot_ci_lower"),
    ind_lo_new = get_path(new_res, "indirect", "boot_ci_lower"),
    ind_hi_old = get_path(old_res, "indirect", "boot_ci_upper"),
    ind_hi_new = get_path(new_res, "indirect", "boot_ci_upper")
  )
}

keys <- intersect(names(OLD$results), names(NEW$results))
rows <- lapply(keys, function(k) build_row(k, OLD$results[[k]], NEW$results[[k]]))
diff.dt <- rbindlist(rows)

# Flag sign or significance flips on indirect path
diff.dt[, sig_old := !is.na(ind_lo_old) & !is.na(ind_hi_old) &
          (ind_lo_old > 0 | ind_hi_old < 0)]
diff.dt[, sig_new := !is.na(ind_lo_new) & !is.na(ind_hi_new) &
          (ind_lo_new > 0 | ind_hi_new < 0)]
diff.dt[, flip := sig_old != sig_new]

# Print key view
fmt <- function(x) ifelse(is.na(x), "      NA", sprintf("%8.4f", x))

cat("\n=== Path estimate diff (old=homoscedastic, new=heteroscedastic) ===\n")
cat(sprintf("%-28s %s\n", "cell",
            "a_old    a_new      b_old    b_new      c'_old   c'_new    ind_old  ind_new"))
for (i in seq_len(nrow(diff.dt))) {
  r <- diff.dt[i]
  cat(sprintf("%-28s %s %s %s %s %s %s %s %s\n", r$cell,
              fmt(r$a_old), fmt(r$a_new),
              fmt(r$b_old), fmt(r$b_new),
              fmt(r$cprime_old), fmt(r$cprime_new),
              fmt(r$ind_old), fmt(r$ind_new)))
}

cat("\n=== Indirect-effect bootstrap CI shifts ===\n")
cat(sprintf("%-28s %-25s %-25s %s\n", "cell", "old CI", "new CI",
            "sig flip"))
for (i in seq_len(nrow(diff.dt))) {
  r <- diff.dt[i]
  old_ci <- sprintf("[%s, %s]", fmt(r$ind_lo_old), fmt(r$ind_hi_old))
  new_ci <- sprintf("[%s, %s]", fmt(r$ind_lo_new), fmt(r$ind_hi_new))
  flip_str <- if (is.na(r$flip)) "?" else if (r$flip) "YES" else "no"
  cat(sprintf("%-28s %-25s %-25s %s (old=%s, new=%s)\n",
              r$cell, old_ci, new_ci, flip_str,
              r$sig_old, r$sig_new))
}

n_flips <- sum(diff.dt$flip, na.rm = TRUE)
cat(sprintf("\nSummary: %d/%d cells changed significance of indirect effect.\n",
            n_flips, nrow(diff.dt)))

# Sanity: Male_MEM_CVR_mimic should roughly match the LRT heteroscedastic fit
lrt <- readRDS(here::here(
  "docs/openmx_audit/lrt_Male_MEM_CVR_mimic.rds"
))
new_row <- diff.dt[cell == "Male_MEM_CVR_mimic"]
cat("\n=== Sanity check: Male_MEM_CVR_mimic ===\n")
cat(sprintf("  Standalone LRT (heteroscedastic): b = %.4f\n",
            lrt$heteroscedastic$b$Estimate))
cat(sprintf("  Production new run:               b = %.4f\n",
            new_row$b_new))
cat(sprintf("  Standalone LRT (heteroscedastic): a = %.6f\n",
            lrt$heteroscedastic$a$Estimate))
cat(sprintf("  Production new run:               a = %.6f\n",
            new_row$a_new))

# Save diff
saveRDS(diff.dt,
        here::here("docs/openmx_audit/diff_old_vs_new.rds"))
fwrite(diff.dt,
       here::here("docs/openmx_audit/diff_old_vs_new.csv"))
cat("\nSaved: docs/openmx_audit/diff_old_vs_new.{rds,csv}\n")
