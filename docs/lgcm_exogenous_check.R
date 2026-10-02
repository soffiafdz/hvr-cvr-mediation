# Refit script-13 mediation models with the current fixed predictor moments
# (mean 0, variance 1, no covariances) and with estimated moments, and
# compare the a, b and c' paths. See docs/lgcm_exogenous_predictors_audit.md.
# Run from the hvr_cvr_mediation repo root:
#   Rscript docs/lgcm_exogenous_check.R
# It evaluates script 13 up to Part 4 (data preparation and model
# builders), so it writes a log file to logs/ but no results.
suppressMessages(library(data.table))
src.v <- readLines("R/scripts/13_lgcm_mediation.R")
part4 <- grep(
  'log_section("Part 4: Fitting Mediation Models")', src.v, fixed = TRUE
)
skip <- grep(
  "if (!FORCE_REGENERATE && file.exists(output.path))", src.v,
  fixed = TRUE
)
stopifnot(length(part4) == 1, length(skip) == 1)
# Never take script 13's "outputs exist, skipping" exit here
src.v[skip] <- "if (FALSE) {"
eval(parse(text = src.v[seq_len(part4 - 1)]))
mxOption(key = "Number of Threads", value = 4)

free_exo.fn <- function(model, preds.v, dat) {
  s.mat <- cov(dat[, preds.v], use = "pairwise.complete.obs")
  mxModel(
    model,
    mxPath(
      from = preds.v, arrows = 2, connect = "unique.pairs",
      free = TRUE, values = s.mat[lower.tri(s.mat, diag = TRUE)]
    ),
    mxPath(
      from = "one", to = preds.v, free = TRUE,
      values = colMeans(dat[, preds.v], na.rm = TRUE)
    )
  )
}
grab.fn <- function(fit) {
  p.df <- summary(fit)$parameters
  p.df <- p.df[match(c("a", "b", "cprime"), p.df$name), ]
  c(
    setNames(p.df$Estimate, paste0(p.df$name, "_est")),
    setNames(p.df$Std.Error, paste0(p.df$name, "_se")),
    m2ll = fit$output$Minus2LogLikelihood,
    status = fit$output$status$code
  )
}

CASES.lst <- list(
  c("Male", "MEM", "FRS"), c("Male", "MEM", "CVR_mimic"),
  c("Female", "MEM", "CVR_mimic")
)
for (cs.v in CASES.lst) {
  sample_name <- cs.v[1]
  domain <- cs.v[2]
  pred_label <- cs.v[3]
  data.dt <- samples.lst[[sample_name]]
  cog_vars <- paste0("PHC_", domain, "_", VALID_TIMEPOINTS)
  pred_var <- PREDICTORS[[pred_label]]$var
  preds.v <- c(pred_var, "AGE_c", "EDUC_c", if (HAS_APOE) "APOE4_num")
  keep_cols <- unique(c(
    HVR_COLS, cog_vars, EDT_COLS, EDT_SQ_COLS, preds.v
  ))
  model_data <- as.data.frame(data.dt[, ..keep_cols])
  cat(sprintf(
    "\n== %s %s %s (n = %d; missing predictor values: %d)\n",
    sample_name, domain, pred_label, nrow(model_data),
    sum(is.na(model_data[, preds.v]))
  ))
  base.mod <- build_mediation.fn(
    cog_vars, HVR_COLS, predictor = pred_var,
    name = "m", data = data.dt, skip_age_cov = FALSE
  )
  base.mod <- mxModel(
    base.mod, mxData(observed = model_data, type = "raw")
  )
  set.seed(42)
  fixed.fit <- mxTryHard(
    base.mod, extraTries = 50, silent = TRUE, bestInitsOutput = FALSE
  )
  set.seed(42)
  free.fit <- mxTryHard(
    free_exo.fn(base.mod, preds.v, model_data),
    extraTries = 50, silent = TRUE, bestInitsOutput = FALSE
  )
  res.mat <- rbind(fixed = grab.fn(fixed.fit), free = grab.fn(free.fit))
  print(signif(res.mat, 7))
  cat(
    "max |diff| a/b/c' estimates:",
    signif(max(abs(res.mat[1, 1:3] - res.mat[2, 1:3])), 3),
    "| SEs:", signif(max(abs(res.mat[1, 4:6] - res.mat[2, 4:6])), 3), "\n"
  )
  chk_path <- file.path(
    checkpoint_dir,
    paste0(paste(sample_name, domain, pred_label, sep = "_"), ".rds")
  )
  if (file.exists(chk_path)) {
    chk.res <- readRDS(chk_path)
    cat(
      "pipeline checkpoint a:", chk.res$a_path$est,
      "se", chk.res$a_path$se, "| b:", chk.res$b_path$est,
      "se", chk.res$b_path$se, "\n"
    )
  }
}
