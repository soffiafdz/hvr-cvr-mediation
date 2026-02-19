# =============================================================================
# Plotting Utilities
# =============================================================================
# Common plotting functions, themes, and color palettes
# =============================================================================

# ----- Color Palettes -----

#' Get standard color palette
#'
#' Returns consistent color palettes for use across all figures.
#' Use manuscript_colors() for the full named list.
#'
#' @param type Color palette type ("default", "sex", "roi", "hvr", "measure")
#' @return Named vector of colors
get_palette <- function(type = "default") {
  palettes <- list(
    default = c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    ),
    sex = c(
      Female = "#8B0000",      # Dark red
      Male = "#191970",        # Midnight blue
      Combined = "#404040"     # Dark gray
    ),
    sex_md = c(
      Female = "<span style='color: #8B0000;'>Female</span>",
      Male = "<span style='color: #191970;'>Male</span>"
    ),
    roi = c(
      HC = "#0072B2",
      LV = "#D55E00",
      HVR = "#009E73",
      ICC = "#999999"
    ),
    hvr = c(
      `Preserved HVR` = "#228B22",  # Forest green
      `Low HVR` = "#B22222"         # Firebrick
    ),
    measure = c(
      `Z-scored` = "#6A0DAD",       # Purple
      `Raw` = "#D2691E"             # Orange
    )
  )

  palettes[[type]]
}

# ----- Plot Themes -----

#' Standard ggplot theme for publications
#' @param base_size Base font size
#' @param use_markdown Whether to use ggtext markdown
#' @return ggplot2 theme object
theme_publication <- function(base_size = 10, use_markdown = FALSE) {
  base_theme <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      plot.caption = ggplot2::element_text(size = base_size - 3),
      legend.position = "bottom"
    )

  if (use_markdown) {
    base_theme <- base_theme +
      ggplot2::theme(
        axis.text.x = ggtext::element_markdown(),
        axis.text.y = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(),
        plot.title = ggtext::element_markdown()
      )
  }

  base_theme
}

# ----- Plot Helpers -----

# Note: Use ensure_directory() from data_io.R instead
# Kept for backward compatibility
#' @export
ensure_plot_dir <- function(plot_dir) {
  ensure_directory(plot_dir)
}

#' Save plot with standard settings
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
save_plot <- function(plot, filename, width = 7, height = 7, dpi = 600) {
  ensure_plot_dir(dirname(filename))

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  invisible(filename)
}

#' Create comparison plot with error bars
#' @param data Data frame
#' @param x X variable
#' @param y Y variable
#' @param group Grouping variable
#' @param ci_lower Lower CI bound
#' @param ci_upper Upper CI bound
#' @param sig Significance indicator
#' @param colors Color palette
#' @return ggplot object
plot_comparison <- function(
    data, x, y, group, ci_lower, ci_upper, sig = NULL, colors = NULL) {
  if (is.null(colors)) {
    colors <- get_palette("default")
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x]], y = .data[[y]],
      colour = .data[[group]], group = .data[[group]]
    )
  )

  if (!is.null(sig)) {
    p <- p + ggplot2::aes(alpha = .data[[sig]])
  }

  p <- p +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[ci_lower]], ymax = .data[[ci_upper]]),
      position = ggplot2::position_dodge(width = 1),
      width = 0
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data[[group]]),
      position = ggplot2::position_dodge(width = 1),
      size = 0.75
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      alpha = 0.5,
      colour = get_palette("default")[1]
    ) +
    ggplot2::scale_colour_manual(values = colors) +
    theme_publication()

  if (!is.null(sig)) {
    p <- p +
      ggplot2::scale_alpha_manual(
        values = c("non-sig" = 0.5, "sig" = 1),
        guide = "none"
      )
  }

  p + ggplot2::coord_flip()
}

#' Create GAMLSS centile plot
#' @param obs_data Observed data
#' @param cent_data Centile predictions
#' @param age_var Age variable name
#' @param value_var Value variable name
#' @param sex_var Sex variable name
#' @param centiles Vector of centile columns
#' @param linetypes Named vector of linetypes
#' @return ggplot object
plot_gamlss_centiles <- function(
    obs_data, cent_data,
    age_var = "AGE", value_var = "VAL_scl",
    sex_var = "Sex", centiles = c("p10", "p25", "p50", "p75", "p90"),
    linetypes = NULL) {
  if (is.null(linetypes)) {
    linetypes <- c(
      "p10" = "dashed",
      "p25" = "dotdash",
      "p50" = "solid",
      "p75" = "dotdash",
      "p90" = "dashed"
    )
  }

  sex_colors <- get_palette("sex")

  # Melt centile data
  cent_long <- data.table::melt(
    cent_data,
    measure.vars = centiles,
    variable.name = "Centiles",
    value.name = "VAL_scl"
  )

  ggplot2::ggplot() +
    theme_publication(use_markdown = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[age_var]],
        y = .data[[value_var]],
        colour = .data[[sex_var]]
      ),
      data = obs_data,
      alpha = 0.05,
      size = 0.3,
      shape = 21
    ) +
    ggplot2::scale_colour_manual(values = sex_colors) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        x = .data[[age_var]],
        y = VAL_scl,
        linetype = Centiles
      ),
      data = cent_long[cent_long[[sex_var]] == "Male", ],
      se = FALSE,
      colour = sex_colors["Male"],
      method = "loess",
      span = 0.35,
      linewidth = 0.6
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        x = .data[[age_var]],
        y = VAL_scl,
        linetype = Centiles
      ),
      data = cent_long[cent_long[[sex_var]] == "Female", ],
      se = FALSE,
      colour = sex_colors["Female"],
      method = "loess",
      span = 0.35,
      linewidth = 0.6
    ) +
    ggplot2::scale_linetype_manual(values = linetypes)
}

#' Wrap text to specified width
#' @param text Text to wrap
#' @param width Maximum line width
#' @param indent Indentation for wrapped lines
#' @return Wrapped text
wrap_text <- function(text, width = 100, indent = 0) {
  stringr::str_wrap(text, width = width, exdent = indent)
}

# ----- GAMLSS Plotting Functions -----

#' Create GAMLSS centile curve plot
#' @param norm_table Data table with normative centiles
#' @param obs_data Optional observed data to overlay
#' @param centiles Vector of centiles to plot
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @return ggplot object
plot_gamlss_curves <- function(
    norm_table,
    obs_data = NULL,
    centiles = c("p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5"),
    title = "Normative Centiles",
    subtitle = NULL,
    x_label = "Age (years)",
    y_label = "Volume (cc)") {
  # Filter available centiles
  centiles <- centiles[centiles %in% names(norm_table)]

  # Melt normative data
  cent_long.dt <- data.table::melt(
    norm_table,
    id.vars = "AGE",
    measure.vars = centiles,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define line types
  linetypes <- c(
    "p2.5" = "dotted", "p5" = "dotted",
    "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid",
    "p75" = "dotdash", "p90" = "dashed",
    "p95" = "dotted", "p97.5" = "dotted"
  )

  # Get sex from data if available
  sex_col <- if ("SEX" %in% names(norm_table)) {
    unique(norm_table$SEX)[1]
  } else {
    "Male"
  }
  line_color <- get_palette("sex")[sex_col]

  # Base plot
  p <- ggplot2::ggplot()

  # Add observed data if provided
  if (!is.null(obs_data) && nrow(obs_data) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = AGE, y = VAL),
        alpha = 0.1, size = 0.5,
        color = line_color
      )
  }

  # Add centile curves
  p <- p +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        group = Centile, linetype = Centile
      ),
      color = line_color,
      linewidth = 0.8
    ) +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      linetype = "Centile"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )

  p
}

#' Create GAMLSS validation plot (centile calibration)
#' @param validation_data Data table with expected and observed centiles
#' @param title Plot title
#' @return ggplot object
plot_gamlss_validation <- function(
    validation_data,
    title = "Centile Calibration") {
  ggplot2::ggplot(
    validation_data,
    ggplot2::aes(x = CENT_exp, y = CENT_obs)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = manuscript_colors()$calibration) +
    ggplot2::geom_line(color = manuscript_colors()$calibration, linewidth = 1) +
    ggplot2::labs(
      title = title,
      subtitle = "Perfect calibration shown by dashed line",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create longitudinal stability plot
#' @param centile_data Data table with T1 and T2 centiles
#' @param title Plot title
#' @return ggplot object
plot_longitudinal_stability <- function(
    centile_data,
    title = "Longitudinal Centile Stability") {
  ggplot2::ggplot(
    centile_data,
    ggplot2::aes(x = CENT_t1, y = CENT_t2)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      alpha = 0.4, size = 2, color = manuscript_colors()$calibration
    ) +
    ggplot2::geom_smooth(
      method = "lm", se = TRUE,
      color = manuscript_colors()$reference,
      fill = manuscript_colors()$reference,
      alpha = 0.2
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Stable centiles should fall along the diagonal",
      x = "Centile at Time 1 (%)",
      y = "Centile at Time 2 (%)",
      caption = sprintf(
        "r = %.3f | Mean change = %.2f",
        cor(centile_data$CENT_t1, centile_data$CENT_t2, use = "complete.obs"),
        mean(centile_data$CENT_t2 - centile_data$CENT_t1, na.rm = TRUE)
      )
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted centile plot by adjustment method
#' @param norm_tables_list Named list of norm tables by adjustment
#' @param obs_data_list Named list of obs data by adjustment
#' @param roi_label ROI label for title
#' @param sex Sex label for title
#' @param sex_color Color for sex-specific plotting
#' @param y_label Y-axis label
#' @return ggplot object
#' @export
plot_gamlss_faceted_by_adjustment <- function(
    norm_tables_list,
    obs_data_list,
    roi_label,
    sex,
    sex_color,
    y_label = "Volume (cc)") {
  # Combine data across adjustments
  norm_combined.dt <- rbindlist(norm_tables_list, fill = TRUE)
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Get centile columns
  cent_cols <- grep("^p[0-9.]+$", names(norm_combined.dt), value = TRUE)
  cent_cols <- cent_cols[cent_cols %in% c(
    "p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5"
  )]

  # Melt for plotting
  cent_long.dt <- data.table::melt(
    norm_combined.dt,
    id.vars = c("AGE", "ADJ"),
    measure.vars = cent_cols,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define aesthetics
  linetypes <- c(
    "p2.5" = "dotted", "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid", "p75" = "dotdash", "p90" = "dashed",
    "p97.5" = "dotted"
  )

  # Create faceted plot
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = VAL),
      alpha = 0.05, size = 0.2, color = sex_color
    ) +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        linetype = Centile, group = Centile
      ),
      color = sex_color, linewidth = 0.6
    ) +
    ggplot2::facet_wrap(~ ADJ, ncol = 2, scales = "free_y") +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x)
    ) +
    ggplot2::labs(
      title = sprintf("Normative Centiles: %s (%s)", roi_label, sex),
      subtitle = "Faceted by adjustment method",
      x = "Age (years)",
      y = y_label,
      linetype = "Centile"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted sex comparison plot
#' @param norm_tables_list Named list (by adj_sex) of norm tables
#' @param obs_data_list Named list (by adj_sex) of obs data
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @param y_label Y-axis label
#' @return ggplot object
#' @export
plot_gamlss_sexcomp_faceted <- function(
    norm_tables_list,
    obs_data_list,
    roi_label,
    sex_colors,
    y_label = "Volume (cc)") {
  # Combine all data
  norm_combined.dt <- rbindlist(norm_tables_list, fill = TRUE)
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Get centile columns (use fewer for clarity)
  cent_cols <- grep("^p[0-9.]+$", names(norm_combined.dt), value = TRUE)
  cent_cols <- cent_cols[cent_cols %in% c("p25", "p50", "p75")]

  # Melt for plotting
  cent_long.dt <- data.table::melt(
    norm_combined.dt,
    id.vars = c("AGE", "ADJ", "SEX"),
    measure.vars = cent_cols,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define aesthetics
  linetypes <- c("p25" = "dashed", "p50" = "solid", "p75" = "dashed")

  # Create faceted comparison plot
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = VAL, color = SEX),
      alpha = 0.03, size = 0.2
    ) +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        color = SEX, linetype = Centile,
        group = interaction(SEX, Centile)
      ),
      linewidth = 0.7
    ) +
    ggplot2::facet_grid(SEX ~ ADJ, scales = "free_y") +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x),
      name = "Centile"
    ) +
    ggplot2::labs(
      title = sprintf("%s Normative Centiles", roi_label),
      subtitle = "Rows: Sex | Columns: Adjustment method",
      x = "Age (years)",
      y = y_label
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2)
    )
}

#' Create faceted validation plot for both sexes
#' @param valid_data_list Named list of validation data (by adj_sex)
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_validation_sexcomp_faceted <- function(
    valid_data_list, roi_label, sex_colors) {
  valid_combined.dt <- rbindlist(valid_data_list)

  ggplot2::ggplot(
    valid_combined.dt,
    ggplot2::aes(x = CENT_exp, y = CENT_obs, color = SEX)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(SEX ~ ADJ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::labs(
      title = sprintf("Centile Calibration: %s", roi_label),
      subtitle = "Rows: Sex | Columns: Adjustment method",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted validation plot (single sex)
#' @param valid_data_list Named list of validation data by adjustment
#' @param roi_label ROI label for title
#' @param sex Sex label for title
#' @return ggplot object
#' @export
plot_validation_faceted <- function(valid_data_list, roi_label, sex) {
  valid_combined.dt <- rbindlist(valid_data_list)

  ggplot2::ggplot(
    valid_combined.dt,
    ggplot2::aes(x = CENT_exp, y = CENT_obs)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.7, color = manuscript_colors()$calibration) +
    ggplot2::geom_line(color = manuscript_colors()$calibration, linewidth = 0.8) +
    ggplot2::facet_wrap(~ ADJ, ncol = 2) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(
      title = sprintf("Centile Calibration: %s (%s)", roi_label, sex),
      subtitle = "Faceted by adjustment method",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted stability plot with summary metrics
#' @param stability_data_list Named list of stability summary data (by adj_sex)
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_stability_faceted <- function(
    stability_data_list, roi_label, sex_colors) {
  stab_combined.dt <- rbindlist(stability_data_list)

  # Reshape to long format for faceting by metric
  stab_long.dt <- data.table::melt(
    stab_combined.dt,
    id.vars = c("SEX", "ADJ"),
    measure.vars = c("CORR", "MEAN_CHG", "SD_CHG", "CROSS"),
    variable.name = "METRIC",
    value.name = "VALUE"
  )

  # Create labels for metrics
  metric_labels <- c(
    CORR = "Correlation (r)",
    MEAN_CHG = "Mean Change",
    SD_CHG = "SD Change",
    CROSS = "Threshold Crossing"
  )

  stab_long.dt[, METRIC_LABEL := metric_labels[METRIC]]

  ggplot2::ggplot(
    stab_long.dt,
    ggplot2::aes(x = ADJ, y = VALUE, fill = SEX)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggplot2::facet_wrap(~ METRIC_LABEL, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_manual(values = sex_colors, name = "Sex") +
    ggplot2::labs(
      title = sprintf("Longitudinal Stability: %s", roi_label),
      subtitle = "Test-retest metrics across adjustment methods",
      x = "Adjustment Method",
      y = "Value"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# ----- Manuscript-Specific Colors -----

#' Color palettes for manuscript figures
#'
#' Consistent color palette used throughout all manuscript figures.
#' IMPORTANT: Use these colors consistently - do not redefine in individual
#' plotting functions.
#'
#' @return Named list of color values
#' @export
manuscript_colors <- function() {
  # Paul Tol Bright + Muted (colorblind-safe)
  # Each pair maximises warm/cool contrast
  list(
    female = "#8B0000",        # Dark red
    male = "#191970",          # Midnight blue
    ukb = "#CCBB44",           # Tol yellow
    adni = "#66CCEE",          # Tol cyan
    hvr_measure = "#AA3377",   # Tol purple
    hc_measure = "#228833",    # Tol green
    cvr_high = "#CC6677",      # Tol muted rose
    cvr_low = "#44AA99",       # Tol seafoam
    frs_pred = "#EE7733",      # Tol bright orange
    cvr_pred = "#4477AA",      # Tol blue
    traj_low = "#117733",      # Forest green
    traj_high = "#882255",     # Wine

    # Pipeline-only (script 15)
    combined = "gray40",
    hvr_preserved = "#228833",
    hvr_low = "#CC6677",
    frs_hist = "#88CCEE",
    ceiling = "gray30",
    age_color = "#882255",
    cvrf_color = "#999933",
    summary_color = "gray50",
    method_se = "#332288",
    method_free = "#DDCC77"
  )
}

# =============================================================================
# LME / LGCM Plotting Functions (Added for manuscript generation)
# =============================================================================

#' Create forest plot for LME effect sizes
#'
#' @param results_dt Data.table with columns: Domain, Sex, Beta, SE, Significant
#' @return ggplot object
plot_forest_effects <- function(results_dt) {
  if (is.null(results_dt) || nrow(results_dt) == 0) return(NULL)

  dt <- copy(results_dt)
  colors <- manuscript_colors()

  # Map domain codes to full names
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  if ("Domain" %in% names(dt) && any(dt$Domain %in% names(domain_map))) {
    dt[, Domain := domain_map[Domain]]
  }

  # Add CIs
  dt[, CI_lower := Beta - 1.96 * SE]
  dt[, CI_upper := Beta + 1.96 * SE]

  # Create plot (no title - use Quarto fig-cap)
  p <- ggplot(dt, aes(x = Beta, y = Domain, color = Sex)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                  width = 0.2,
                  position = position_dodge(width = 0.4)) +
    geom_point(aes(shape = Significant),
               size = 3,
               position = position_dodge(width = 0.4)) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
      name = ""
    ) +
    scale_color_manual(
      values = c("Male" = colors$male,
                 "Female" = colors$female,
                 "Combined" = colors$combined),
      name = "Sex"
    ) +
    labs(
      x = expression("Standardized " * beta * " (95% CI)"),
      y = NULL
    ) +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom",
          legend.box = "horizontal")

  p
}

#' Create mediation path diagram
#'
#' @param mediation_result Single mediation result
#' @param domain_label Optional title (e.g. "Memory")
#' @param predictor_label Node label (default "FRS")
#' @return ggplot object
# -- old plot_mediation_diagram commented out --
# -- replaced by improved version below (was
# -- plot_mediation_improved); delete after
# -- verifying render
# plot_mediation_diagram_old <- function(
#     mediation_result, domain_label = NULL,
#     predictor_label = "FRS") {
#   ... (171 lines removed) ...
# }

#' Create interaction plot for FRS x HVR effects
#'
#' Shows cognitive trajectories by HVR group, faceted by FRS level.
#' Groups are defined at BASELINE and held constant throughout follow-up.
#' The key message: At low FRS, high HVR strongly protects cognition (large gap
#' between High/Low HVR). At high FRS, this protection
#' is weakened (smaller gap).
#'
#' @param cohort_dt Data.table with cognitive and HVR data
#' @param outcome_var Name of cognitive outcome variable (default: "MEM")
#' @return ggplot object
plot_interaction_frs_hvr <- function(cohort_dt, outcome_var = "MEM") {
  if (is.null(cohort_dt) || nrow(cohort_dt) == 0) return(NULL)

  dt <- copy(cohort_dt)
  colors <- manuscript_colors()

  # Get baseline for median splits - BASELINE ONLY
  if ("YRS_from_bl" %in% names(dt)) {
    baseline <- dt[YRS_from_bl == 0]
  } else if ("VISIT" %in% names(dt)) {
    baseline <- dt[VISIT == 1]
  } else {
    baseline <- dt[!duplicated(PTID)]
  }

  hvr_median <- median(baseline$HVR_z, na.rm = TRUE)
  frs_median <- median(baseline$FRS, na.rm = TRUE)

  # Assign groups at BASELINE and merge back to all timepoints
  baseline_groups <- baseline[, .(
    PTID,
    HVR_group = factor(
      fifelse(HVR_z > hvr_median, "Preserved HVR", "Low HVR"),
      levels = c("Preserved HVR", "Low HVR")
    ),
    FRS_group = factor(
      fifelse(FRS > frs_median, "High CV Risk", "Low CV Risk"),
      levels = c("Low CV Risk", "High CV Risk")
    )
  )]

  # Merge baseline groups to all timepoints
  dt <- merge(dt, baseline_groups, by = "PTID", all.x = TRUE)

  # Determine time variable
  time_var <- if ("YRS_from_bl" %in% names(dt)) "YRS_from_bl" else "YRS"
  if (!time_var %in% names(dt)) {
    warning("Time variable not found in data")
    return(NULL)
  }

  # Aggregate trajectories by year bins
  traj_summary <- dt[!is.na(get(outcome_var)) & !is.na(HVR_group), .(
    Mean_Cog = mean(get(outcome_var), na.rm = TRUE),
    SE_Cog = sd(get(outcome_var), na.rm = TRUE) / sqrt(.N),
    N = .N
  ), by = .(HVR_group, FRS_group, YRS_bin = round(get(time_var)))]

  traj_summary <- traj_summary[N >= 10]

  if (nrow(traj_summary) == 0) {
    warning("No data after filtering for minimum cell size")
    return(NULL)
  }

  # CVR = color, HVR = linetype (consistent)
  frs_colors <- c(
    "Low CV Risk" = colors$cvr_low,
    "High CV Risk" = colors$cvr_high
  )
  hvr_ltypes <- c(
    "Preserved HVR" = "solid",
    "Low HVR" = "dashed"
  )
  # Point shapes: filled for Preserved, open for Low
  hvr_shapes <- c(
    "Preserved HVR" = 16,
    "Low HVR" = 21
  )

  p <- ggplot(
    traj_summary,
    aes(
      x = YRS_bin, y = Mean_Cog,
      color = FRS_group,
      linetype = HVR_group
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Mean_Cog - SE_Cog,
        ymax = Mean_Cog + SE_Cog,
        fill = FRS_group
      ),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(
      aes(shape = HVR_group),
      size = 2, fill = "white", stroke = 0.8
    ) +
    scale_color_manual(
      values = frs_colors,
      name = "FRS (Baseline)"
    ) +
    scale_fill_manual(
      values = frs_colors, guide = "none"
    ) +
    scale_linetype_manual(
      values = hvr_ltypes,
      name = "HVR (Baseline)"
    ) +
    scale_shape_manual(
      values = hvr_shapes,
      name = "HVR (Baseline)"
    ) +
    labs(
      x = "Years from Baseline",
      y = switch(
        outcome_var,
        MEM = "Memory Composite",
        LAN = "Language Composite",
        EXF = "Executive Composite",
        paste(outcome_var, "Score")
      )
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(1.5, "cm")
    )

  p
}

#' Plot HVR group difference over time
#'
#' Creates a plot showing the difference in cognitive scores between
#' Preserved HVR and Low HVR groups over time, stratified by CV risk.
#'
#' @param cohort_dt Data.table with cognitive and HVR data
#' @param outcome_var Name of cognitive outcome variable (default: "MEM")
#' @return ggplot object
#' @export
plot_hvr_difference <- function(cohort_dt, outcome_var = "MEM") {
  if (is.null(cohort_dt) || nrow(cohort_dt) == 0) return(NULL)

  dt <- copy(cohort_dt)
  colors <- manuscript_colors()

  # Get baseline for median splits
  if ("YRS_from_bl" %in% names(dt)) {
    baseline <- dt[YRS_from_bl == 0]
  } else if ("VISIT" %in% names(dt)) {
    baseline <- dt[VISIT == 1]
  } else {
    baseline <- dt[!duplicated(PTID)]
  }

  hvr_median <- median(baseline$HVR_z, na.rm = TRUE)
  frs_median <- median(baseline$FRS, na.rm = TRUE)

  # Assign groups at baseline
  baseline_groups <- baseline[, .(
    PTID,
    HVR_group = fifelse(HVR_z > hvr_median, "Preserved", "Low"),
    FRS_group = factor(
      fifelse(FRS > frs_median, "High CV Risk", "Low CV Risk"),
      levels = c("Low CV Risk", "High CV Risk")
    )
  )]

  dt <- merge(dt, baseline_groups, by = "PTID", all.x = TRUE)

  # Determine time variable
  time_var <- if ("YRS_from_bl" %in% names(dt)) "YRS_from_bl" else "YRS"

  # Aggregate by HVR group, FRS group, and year
  traj_summary <- dt[!is.na(get(outcome_var)) & !is.na(HVR_group), .(
    Mean_Cog = mean(get(outcome_var), na.rm = TRUE),
    SE_Cog = sd(get(outcome_var), na.rm = TRUE) / sqrt(.N),
    N = .N
  ), by = .(HVR_group, FRS_group, YRS_bin = round(get(time_var)))]

  traj_summary <- traj_summary[N >= 5]

  # Calculate difference (Preserved - Low) at each time point by FRS group
  diff_dt <- dcast(traj_summary, FRS_group + YRS_bin ~ HVR_group,
                   value.var = c("Mean_Cog", "SE_Cog", "N"))

  # Filter to time points with both groups
  diff_dt <- diff_dt[!is.na(Mean_Cog_Preserved) & !is.na(Mean_Cog_Low)]

  diff_dt[, `:=`(
    Difference = Mean_Cog_Preserved - Mean_Cog_Low,
    SE_Diff = sqrt(SE_Cog_Preserved^2 + SE_Cog_Low^2)
  )]

  if (nrow(diff_dt) == 0) {
    warning("No data for difference calculation")
    return(NULL)
  }

  # Create difference plot
  p <- ggplot(diff_dt, aes(x = YRS_bin, y = Difference, color = FRS_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = Difference - SE_Diff, ymax = Difference + SE_Diff,
                    fill = FRS_group), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(
      "Low CV Risk" = colors$cvr_low,
      "High CV Risk" = colors$cvr_high
    ), name = "CV Risk Level") +
    scale_fill_manual(values = c(
      "Low CV Risk" = colors$cvr_low,
      "High CV Risk" = colors$cvr_high
    ), name = "CV Risk Level") +
    labs(
      x = "Years from Baseline",
      y = paste("Difference (Preserved - Low HVR)")
    ) +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom")

  p
}

#' Trajectory panel: one sex, all domains, with diff
#'
#' Creates a compact trajectory figure for one sex.
#' Left column: 4 trajectory lines (2 HVR x 2 CVR).
#' Right column: HVR difference (Preserved - Low) per
#' CVR level. Row facets for cognitive domain.
#'
#' @param sex_dt Data.table filtered to one sex
#' @param cvr_col Name of CVR column ("FRS" or
#'   "CVR_mimic")
#' @param cvr_label Label for CVR grouping in legend
#' @param min_n Minimum cell size (default 10)
#' @return patchwork object
#' @export
plot_lme_trajectory_panel <- function(
    sex_dt,
    cvr_col = "FRS",
    cvr_label = "CV Risk",
    min_n = 10
) {
  if (is.null(sex_dt) || nrow(sex_dt) == 0) {
    return(NULL)
  }

  dt <- copy(sex_dt)
  colors <- manuscript_colors()
  domains.v <- c("MEM", "LAN", "EXF")
  domain_labels <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )

  # -- baseline median splits --
  if ("YRS_from_bl" %in% names(dt)) {
    baseline <- dt[YRS_from_bl == 0]
  } else if ("VISIT" %in% names(dt)) {
    baseline <- dt[VISIT == 1]
  } else {
    baseline <- dt[!duplicated(PTID)]
  }

  hvr_med <- median(baseline$HVR_z, na.rm = TRUE)
  cvr_med <- median(
    baseline[[cvr_col]], na.rm = TRUE
  )

  baseline_grp <- baseline[, .(
    PTID,
    HVR_group = factor(
      fifelse(
        HVR_z > hvr_med,
        "Preserved HVR", "Low HVR"
      ),
      levels = c("Preserved HVR", "Low HVR")
    ),
    CVR_group = factor(
      fifelse(
        get(cvr_col) >= cvr_med,
        paste("High", cvr_label),
        paste("Low", cvr_label)
      ),
      levels = c(
        paste("Low", cvr_label),
        paste("High", cvr_label)
      )
    )
  )]

  dt <- merge(
    dt, baseline_grp, by = "PTID", all.x = TRUE
  )

  time_var <- if ("YRS_from_bl" %in% names(dt)) {
    "YRS_from_bl"
  } else {
    "YRS"
  }

  # -- pivot to long on domain --
  long.dt <- rbindlist(lapply(domains.v, function(d) {
    out <- dt[
      !is.na(get(d)) & !is.na(HVR_group),
      .(
        Mean = mean(get(d), na.rm = TRUE),
        SE = sd(get(d), na.rm = TRUE) / sqrt(.N),
        N = .N
      ),
      by = .(HVR_group, CVR_group,
             YRS_bin = round(get(time_var)))
    ]
    out[, Domain := d]
    out
  }))

  long.dt <- long.dt[N >= min_n]
  long.dt[, Domain := factor(
    domain_labels[Domain],
    levels = domain_labels
  )]

  if (nrow(long.dt) == 0) return(NULL)

  # -- compute HVR difference per domain/CVR/year --
  diff.dt <- dcast(
    long.dt,
    CVR_group + YRS_bin + Domain ~ HVR_group,
    value.var = c("Mean", "SE", "N")
  )
  setnames(
    diff.dt,
    c("Mean_Preserved HVR", "Mean_Low HVR",
      "SE_Preserved HVR", "SE_Low HVR"),
    c("Mean_P", "Mean_L", "SE_P", "SE_L"),
    skip_absent = TRUE
  )
  diff.dt <- diff.dt[
    !is.na(Mean_P) & !is.na(Mean_L)
  ]
  diff.dt[, `:=`(
    Diff = Mean_P - Mean_L,
    SE_Diff = sqrt(SE_P^2 + SE_L^2)
  )]

  hvr_cols <- c(
    "Preserved HVR" = colors$hvr_preserved,
    "Low HVR" = colors$hvr_low
  )
  cvr_cols <- c(colors$cvr_low, colors$cvr_high)
  names(cvr_cols) <- levels(long.dt$CVR_group)

  # -- left panel: trajectories --
  p_traj <- ggplot(
    long.dt,
    aes(
      x = YRS_bin, y = Mean,
      color = HVR_group, linetype = CVR_group
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Mean - SE, ymax = Mean + SE,
        fill = HVR_group, group = interaction(
          HVR_group, CVR_group
        )
      ),
      alpha = 0.12, color = NA
    ) +
    geom_line(
      aes(group = interaction(
        HVR_group, CVR_group
      )),
      linewidth = 0.9
    ) +
    geom_point(size = 1.5) +
    facet_wrap(
      ~ Domain, ncol = 1, scales = "free_y"
    ) +
    scale_color_manual(
      values = hvr_cols, name = "HVR (Baseline)"
    ) +
    scale_fill_manual(
      values = hvr_cols, guide = "none"
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed"),
      name = cvr_label
    ) +
    labs(x = "Years from Baseline", y = "Z-Score") +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = element_text(
        size = 9, face = "bold"
      )
    )

  # -- right panel: difference --
  p_diff <- ggplot(
    diff.dt,
    aes(
      x = YRS_bin, y = Diff,
      color = CVR_group, fill = CVR_group
    )
  ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray30", linewidth = 0.6
    ) +
    geom_ribbon(
      aes(
        ymin = Diff - SE_Diff,
        ymax = Diff + SE_Diff
      ),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    facet_wrap(
      ~ Domain, ncol = 1, scales = "free_y"
    ) +
    scale_color_manual(
      values = cvr_cols, name = cvr_label
    ) +
    scale_fill_manual(
      values = cvr_cols, guide = "none"
    ) +
    labs(
      x = "Years from Baseline",
      y = "Preserved \u2212 Low HVR"
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(
        size = 9, face = "bold"
      )
    )

  p_traj + p_diff +
    plot_layout(
      widths = c(3, 2), guides = "collect"
    ) &
    theme(legend.position = "bottom")
}

#' Domain-level trajectory figure (redesign)
#'
#' Single domain, both sexes, both CVR sources.
#' 2x2 layout: rows = Sex, cols = panel type.
#' Each panel faceted by CVR source (FRS / MIMIC).
#'
#' @param cohort_dt Full cohort data.table
#' @param domain One of "MEM", "LAN", "EXF"
#' @param min_n Minimum cell size to plot
#' @return patchwork object or NULL
#' @export
plot_domain_trajectory <- function(
    cohort_dt, domain, min_n = 5) {
  if (is.null(cohort_dt) ||
      nrow(cohort_dt) == 0) {
    return(NULL)
  }

  colors <- manuscript_colors()
  domain_labels <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )

  time_var <- if (
    "YRS_from_bl" %in% names(cohort_dt)
  ) "YRS_from_bl" else "YRS"

  # -- helper: per-sex/CVR source data --
  build_panel.fn <- function(
      sex, cvr_col, cvr_label) {
    dt <- copy(cohort_dt[SEX == sex])
    if (nrow(dt) == 0 ||
        !cvr_col %in% names(dt)) {
      return(NULL)
    }

    # Baseline median splits
    if ("YRS_from_bl" %in% names(dt)) {
      bl <- dt[YRS_from_bl == 0]
    } else if ("VISIT" %in% names(dt)) {
      bl <- dt[VISIT == 1]
    } else {
      bl <- dt[!duplicated(PTID)]
    }

    hvr_med <- median(
      bl$HVR_z, na.rm = TRUE
    )
    cvr_vals <- bl[[cvr_col]]
    if (all(is.na(cvr_vals))) return(NULL)
    cvr_med <- median(
      cvr_vals, na.rm = TRUE
    )

    grp <- bl[
      !is.na(HVR_z) & !is.na(get(cvr_col)),
      .(
        PTID,
        HVR_group = factor(fifelse(
          HVR_z > hvr_med,
          "Preserved HVR", "Low HVR"
        ), levels = c(
          "Preserved HVR", "Low HVR"
        )),
        CVR_group = factor(fifelse(
          get(cvr_col) >= cvr_med,
          "High CV Risk", "Low CV Risk"
        ), levels = c(
          "Low CV Risk", "High CV Risk"
        ))
      )
    ]

    dt <- merge(
      dt, grp, by = "PTID", all.x = TRUE
    )

    # Aggregate
    agg <- dt[
      !is.na(get(domain)) &
        !is.na(HVR_group) &
        !is.na(CVR_group),
      .(
        Mean = mean(
          get(domain), na.rm = TRUE
        ),
        SE = sd(
          get(domain), na.rm = TRUE
        ) / sqrt(.N),
        N = .N
      ),
      by = .(
        HVR_group, CVR_group,
        YRS_bin = round(get(time_var))
      )
    ][N >= min_n]

    if (nrow(agg) == 0) return(NULL)

    agg[, Sex := sex]
    agg[, CVR_source := cvr_label]

    # HVR difference
    wide <- dcast(
      agg,
      CVR_group + YRS_bin + Sex +
        CVR_source ~ HVR_group,
      value.var = c("Mean", "SE", "N")
    )
    setnames(
      wide,
      c(
        "Mean_Preserved HVR",
        "Mean_Low HVR",
        "SE_Preserved HVR",
        "SE_Low HVR"
      ),
      c("Mean_P", "Mean_L",
        "SE_P", "SE_L"),
      skip_absent = TRUE
    )
    wide <- wide[
      !is.na(Mean_P) & !is.na(Mean_L)
    ]
    if (nrow(wide) > 0) {
      wide[, `:=`(
        Diff = Mean_P - Mean_L,
        SE_Diff = sqrt(SE_P^2 + SE_L^2)
      )]
    }

    list(traj = agg, diff = wide)
  }

  # Build all 4 combos (2 sexes x 2 CVR)
  panels <- list()
  cvr_specs <- list(
    c("FRS", "FRS"),
    c("CVR_mimic", "CVR[MIMIC]")
  )
  for (sx in c("Male", "Female")) {
    for (cv in cvr_specs) {
      r <- build_panel.fn(
        sx, cv[1], cv[2]
      )
      if (!is.null(r)) {
        panels <- c(panels, list(r))
      }
    }
  }
  if (length(panels) == 0) return(NULL)

  traj.dt <- rbindlist(
    lapply(panels, `[[`, "traj")
  )
  diff.dt <- rbindlist(
    lapply(panels, `[[`, "diff"),
    fill = TRUE
  )

  # Factor ordering
  sex_lvl <- c("Male", "Female")
  sex_lbl <- c("Males", "Females")
  cvr_lvl <- c("FRS", "CVR[MIMIC]")
  traj.dt[, Sex := factor(
    Sex, levels = sex_lvl,
    labels = sex_lbl
  )]
  diff.dt[, Sex := factor(
    Sex, levels = sex_lvl,
    labels = sex_lbl
  )]
  traj.dt[, CVR_source := factor(
    CVR_source, levels = cvr_lvl
  )]
  diff.dt[, CVR_source := factor(
    CVR_source, levels = cvr_lvl
  )]

  cvr_cols <- c(
    "Low CV Risk" = colors$traj_low,
    "High CV Risk" = colors$traj_high
  )

  base_sz <- 10
  strip_sz <- 9

  # -- helper: trajectory panel --
  make_traj.fn <- function(data, title) {
    ggplot(
      data,
      aes(
        x = YRS_bin, y = Mean,
        color = CVR_group,
        linetype = HVR_group,
        shape = HVR_group,
        group = interaction(
          HVR_group, CVR_group
        )
      )
    ) +
      geom_ribbon(
        aes(
          ymin = Mean - SE,
          ymax = Mean + SE,
          fill = CVR_group
        ),
        alpha = 0.10, color = NA,
        show.legend = FALSE
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(
        size = 1.5, fill = "white"
      ) +
      facet_grid(
        CVR_source ~ .,
        scales = "free_y",
        labeller = label_parsed
      ) +
      scale_color_manual(
        values = cvr_cols,
        name = "CV Risk"
      ) +
      scale_fill_manual(
        values = cvr_cols,
        guide = "none"
      ) +
      scale_linetype_manual(
        values = c(
          "Preserved HVR" = "solid",
          "Low HVR" = "dashed"
        ),
        name = "HVR"
      ) +
      scale_shape_manual(
        values = c(
          "Preserved HVR" = 16,
          "Low HVR" = 21
        ),
        name = "HVR"
      ) +
      labs(
        x = "Years",
        y = paste(
          domain_labels[domain],
          "Z-Score"
        ),
        title = title
      ) +
      theme_publication(
        base_size = base_sz
      ) +
      theme(
        strip.text = element_text(
          size = strip_sz,
          face = "bold"
        ),
        plot.title = element_text(
          size = 10, face = "bold",
          hjust = 0.5
        )
      )
  }

  # -- helper: difference panel --
  make_diff.fn <- function(data, title) {
    ggplot(
      data,
      aes(
        x = YRS_bin, y = Diff,
        color = CVR_group,
        fill = CVR_group
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "gray50"
      ) +
      geom_ribbon(
        aes(
          ymin = Diff - SE_Diff,
          ymax = Diff + SE_Diff
        ),
        alpha = 0.15, color = NA
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(
        size = 1.5, fill = "white"
      ) +
      facet_grid(
        CVR_source ~ .,
        scales = "free_y",
        labeller = label_parsed
      ) +
      scale_color_manual(
        values = cvr_cols,
        name = "CV Risk"
      ) +
      scale_fill_manual(
        values = cvr_cols,
        guide = "none"
      ) +
      labs(
        x = "Years",
        y = "Cog Difference\n(Preserved \u2212 Low HVR)",
        title = title
      ) +
      theme_publication(
        base_size = base_sz
      ) +
      theme(
        strip.text = element_text(
          size = strip_sz,
          face = "bold"
        ),
        plot.title = element_text(
          size = 10, face = "bold",
          hjust = 0.5
        )
      )
  }

  # Build 4 panels (2x2)
  m_traj <- traj.dt[Sex == "Males"]
  f_traj <- traj.dt[Sex == "Females"]
  m_diff <- diff.dt[Sex == "Males"]
  f_diff <- diff.dt[Sex == "Females"]

  p_traj_m <- if (nrow(m_traj) > 0) {
    make_traj.fn(m_traj, "Males")
  } else {
    plot_spacer()
  }
  p_diff_m <- if (nrow(m_diff) > 0) {
    make_diff.fn(m_diff, "")
  } else {
    plot_spacer()
  }
  p_traj_f <- if (nrow(f_traj) > 0) {
    make_traj.fn(f_traj, "Females")
  } else {
    plot_spacer()
  }
  p_diff_f <- if (nrow(f_diff) > 0) {
    make_diff.fn(f_diff, "")
  } else {
    plot_spacer()
  }

  # 2x2: rows=sex, cols=panel type
  (p_traj_m | p_diff_m) /
    (p_traj_f | p_diff_f) +
    plot_layout(
      widths = c(3, 2),
      guides = "collect"
    ) +
    plot_annotation(
      tag_levels = "A"
    ) &
    theme(legend.position = "bottom")
}

#' Create sensitivity comparison figure
#'
#' Compares effects from two analyses (e.g., z-scored vs raw HVR) using a
#' paired forest plot that clearly shows concordance between methods.
#'
#' @param results1_dt First results data.table with Domain, Sex, Beta, SE, p
#' @param results2_dt Second results data.table with Domain, Sex, Beta, SE, p
#' @param label1 Label for first analysis
#' @param label2 Label for second analysis
#' @return ggplot object
plot_sensitivity_comparison <- function(results1_dt, results2_dt,
                                         label1 = "Z-scored HVR",
                                         label2 = "Raw HVR") {

  if (is.null(results1_dt) || is.null(results2_dt)) return(NULL)

  colors <- manuscript_colors()

  # Prepare data - stack the two analyses
  dt1 <- as.data.table(copy(results1_dt))
  dt2 <- as.data.table(copy(results2_dt))

  # Determine which columns to keep (only those present in both)
  common_cols <- intersect(names(dt1), names(dt2))
  keep_cols <- intersect(c("Domain", "Sex", "Beta", "SE"), common_cols)

  if (length(keep_cols) < 3) {
    warning("Insufficient columns for sensitivity comparison")
    return(NULL)
  }

  dt1 <- dt1[, .SD, .SDcols = keep_cols]
  dt2 <- dt2[, .SD, .SDcols = keep_cols]

  dt1[, Analysis := label1]
  dt2[, Analysis := label2]

  comp_dt <- rbind(dt1, dt2)

  if (nrow(comp_dt) == 0) return(NULL)

  # Map domain codes to full names
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  has_domain <- "Domain" %in% names(comp_dt)
  if (has_domain &&
      any(comp_dt$Domain %in% names(domain_map))) {
    comp_dt[, Domain := domain_map[Domain]]
  }

  # Add CIs
  comp_dt[, CI_lower := Beta - 1.96 * SE]
  comp_dt[, CI_upper := Beta + 1.96 * SE]

  # Create Domain-Sex label for y-axis
  if ("Sex" %in% names(comp_dt)) {
    comp_dt[, Label := paste(Domain, Sex, sep = " - ")]
  } else {
    comp_dt[, Label := Domain]
  }

  # Order by domain then sex
  comp_dt[, Label := factor(Label, levels = unique(Label))]

  # Colors for the two analyses
  analysis_colors <- c(colors$adni, colors$ukb)
  names(analysis_colors) <- c(label1, label2)

  # Paired forest plot with facets by Domain
  p <- ggplot(comp_dt, aes(x = Beta, y = Analysis, color = Analysis)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                  width = 0.3,
                  linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = analysis_colors, name = "Brain Measure") +
    facet_wrap(~Domain, ncol = 3) +
    labs(
      x = expression("Effect Size (" * beta * ") with 95% CI"),
      y = NULL
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold")
    )

  p
}

# =============================================================================
# Trajectory and Distribution Plots (for manuscript figures)
# =============================================================================

#' Plot Predicted Cognitive Trajectories by FRS x HVR Stratification
#'
#' Creates visualization of model-predicted cognitive trajectories stratified
#' by High/Low FRS and High/Low HVR combinations.
#' Groups are defined at baseline.
#'
#' @param model.fit lme model object from nlme (or list containing $model)
#' @param data.dt Analysis cohort data.table with YRS_from_bl, FRS, HVR_z
#' @param domain Character: cognitive domain label ("MEM", "LAN", or "EXF")
#' @param sex Character: sex label ("Male", "Female", or NULL for combined)
#' @return ggplot object
#' @export
plot_trajectory_by_frs_hvr <- function(model.fit, data.dt, domain, sex = NULL) {

  if (is.null(model.fit) || is.null(data.dt)) return(NULL)

  colors <- manuscript_colors()

  # Handle case where model is wrapped in a list
  if (is.list(model.fit) && "model" %in% names(model.fit)) {
    model.fit <- model.fit$model
  }
  if (is.null(model.fit)) return(NULL)

  # Ensure data.table
  dt <- as.data.table(data.dt)

  # Compute baseline age if not present (Age_bl = AGE at YRS_from_bl == 0)
  if (!"Age_bl" %in% names(dt) && "AGE" %in% names(dt)) {
    baseline_ages <- dt[YRS_from_bl == 0, .(Age_bl = AGE[1]), by = PTID]
    dt <- merge(dt, baseline_ages, by = "PTID", all.x = TRUE)
  }

  # Get max follow-up
  max_yrs <- max(dt$YRS_from_bl, na.rm = TRUE)
  yrs_seq.v <- seq(0, max_yrs, length.out = 50)

  # Define High/Low based on median split at BASELINE
  frs_med <- median(dt$FRS, na.rm = TRUE)
  frs_sd <- sd(dt$FRS, na.rm = TRUE)
  hvr_med <- median(dt$HVR_z, na.rm = TRUE)

  # Get mean baseline age
  age_mean <- mean(dt$Age_bl, na.rm = TRUE)

  # Create prediction grid
  pred_grid.dt <- data.table(
    YRS_from_bl = rep(yrs_seq.v, 4),
    FRS = rep(c(frs_med - frs_sd, frs_med - frs_sd,
                frs_med + frs_sd, frs_med + frs_sd), each = length(yrs_seq.v)),
    HVR_z = rep(c(hvr_med + 1, hvr_med - 1,
                  hvr_med + 1, hvr_med - 1), each = length(yrs_seq.v)),
    Age_bl = age_mean,
    EDUC = as.integer(mean(dt$EDUC, na.rm = TRUE)),
    APOE4 = FALSE
  )

  # Create group labels
  pred_grid.dt[, FRS_group := fifelse(FRS < frs_med, "Low FRS", "High FRS")]
  pred_grid.dt[, HVR_group := fifelse(
    HVR_z > hvr_med, "Preserved HVR", "Low HVR"
  )]
  pred_grid.dt[, Group := paste(FRS_group, HVR_group, sep = " + ")]

  # Get predictions (fixed effects only)
  tryCatch({
    is_lmer <- inherits(model.fit, "lmerMod") ||
      inherits(model.fit, "lmerModLmerTest")
    if (is_lmer) {
      pred_grid.dt[, predicted := predict(
        model.fit, newdata = pred_grid.dt,
        re.form = NA, allow.new.levels = TRUE
      )]
    } else {
      pred_grid.dt[, predicted := predict(
        model.fit,
        newdata = pred_grid.dt, level = 0
      )]
    }
  }, error = function(e) {
    warning("Could not generate predictions: ", e$message)
    return(NULL)
  })

  if (!"predicted" %in% names(pred_grid.dt)) {
    warning("Predictions not generated")
    return(NULL)
  }

  colors <- manuscript_colors()
  frs_colors <- c(
    "Low FRS" = colors$cvr_low,
    "High FRS" = colors$cvr_high
  )
  hvr_ltypes <- c(
    "Preserved HVR" = "solid",
    "Low HVR" = "dashed"
  )

  hvr_shapes <- c(
    "Preserved HVR" = 16,
    "Low HVR" = 21
  )

  p <- ggplot(
    pred_grid.dt,
    aes(
      x = YRS_from_bl, y = predicted,
      color = FRS_group,
      linetype = HVR_group
    )
  ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = frs_colors,
      name = "FRS (Baseline)"
    ) +
    scale_linetype_manual(
      values = hvr_ltypes,
      name = "HVR (Baseline)"
    ) +
    scale_shape_manual(
      values = hvr_shapes,
      name = "HVR (Baseline)"
    ) +
    labs(
      x = "Years from Baseline",
      y = switch(
        domain,
        MEM = "Memory Composite",
        LAN = "Language Composite",
        EXF = "Executive Composite",
        paste(domain, "Composite")
      )
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.key.width = unit(1.5, "cm")
    )

  p
}

#' Plot z-score distributions for HVR and HC
#'
#' Faceted histogram+density comparing observed ADNI
#' z-scores against expected N(0,1).
#'
#' @param data.dt data.table with HVR_z and HC_z columns
#' @param sample_label Label prefix for mean annotation
#' @return ggplot object
#' @export
plot_zscore_distribution <- function(
    data.dt, sample_label = "Comparable") {

  if (is.null(data.dt)) return(NULL)
  dt <- as.data.table(data.dt)

  colors <- manuscript_colors()
  measure_cols <- c(
    HVR = colors$hvr_measure,
    HC = colors$hc_measure
  )
  measure_labels <- c(
    HVR = "Hippocampal-to-Ventricle Ratio",
    HC = "Hippocampal Volume"
  )

  # Melt to long format
  long.dt <- melt(
    dt[, .(PTID, HVR_z, HC_z)],
    id.vars = "PTID",
    variable.name = "Measure",
    value.name = "z_score"
  )
  long.dt[, Measure := gsub("_z$", "", Measure)]
  long.dt[, Measure := factor(
    Measure, levels = c("HVR", "HC")
  )]

  # Per-facet means for annotations
  means.dt <- long.dt[, .(
    mean_z = mean(z_score, na.rm = TRUE)
  ), by = Measure]

  # Symmetric x-axis
  z_range <- range(long.dt$z_score, na.rm = TRUE)
  x_lim <- max(abs(z_range)) + 0.5

  p <- ggplot(long.dt, aes(x = z_score)) +
    geom_histogram(
      aes(
        y = after_stat(density),
        fill = Measure
      ),
      bins = 40, alpha = 0.7, color = "white"
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      geom = "area",
      fill = colors$ukb, alpha = 0.1
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      color = colors$ukb,
      linewidth = 0.5, linetype = "solid"
    ) +
    geom_vline(
      data = means.dt,
      aes(xintercept = mean_z, color = Measure),
      linewidth = 0.6, linetype = "dashed"
    ) +
    geom_vline(
      xintercept = 0, color = colors$ukb,
      linewidth = 0.6, linetype = "dashed"
    ) +
    geom_text(
      data = means.dt,
      aes(
        x = mean_z - 0.15, y = 0.42,
        label = sprintf("%.2f", mean_z),
        color = Measure
      ),
      hjust = 1, size = 2.5
    ) +
    annotate(
      "text", x = 0.15, y = 0.42,
      label = "Expected: 0",
      color = colors$ukb,
      hjust = 0, size = 2.5
    ) +
    scale_fill_manual(
      values = measure_cols, guide = "none"
    ) +
    scale_color_manual(
      values = measure_cols, guide = "none"
    ) +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    facet_wrap(
      ~Measure,
      labeller = labeller(
        Measure = measure_labels
      )
    ) +
    labs(
      x = "Z-Score (vs UK Biobank)", y = "Density"
    ) +
    theme_publication(base_size = 10)

  p
}

#' Plot age distribution comparison: UKB vs ADNI
#'
#' Density plot comparing UKB normative and ADNI comparable
#' subsample age distributions. Both simulated from summary
#' statistics (individual UKB data not distributable).
#'
#' @param ukb_age_mean UKB mean age
#' @param ukb_age_sd UKB age SD
#' @param ukb_age_min UKB min age for truncation
#' @param ukb_age_max UKB max age for truncation
#' @param ukb_n UKB N for simulation
#' @param adni_age_mean ADNI comparable subsample mean age
#' @param adni_age_sd ADNI comparable subsample age SD
#' @param adni_n ADNI comparable subsample N
#' @return ggplot object
plot_age_distribution_comparison <- function(
    ukb_age_mean,
    ukb_age_sd,
    ukb_age_min,
    ukb_age_max,
    ukb_n,
    adni_age_mean,
    adni_age_sd,
    adni_n) {

  set.seed(42)

  # Simulate UKB (truncated to actual range)
  sim_n <- min(ukb_n, 5000)
  ukb_sim <- rnorm(
    sim_n * 2, mean = ukb_age_mean, sd = ukb_age_sd
  )
  ukb_sim <- ukb_sim[
    ukb_sim >= ukb_age_min & ukb_sim <= ukb_age_max
  ]
  ukb_sim <- ukb_sim[seq_len(
    min(sim_n, length(ukb_sim))
  )]

  # ADNI comparable subsample (from summary stats)
  adni_sim <- rnorm(
    adni_n, mean = adni_age_mean, sd = adni_age_sd
  )
  adni_sim <- adni_sim[
    adni_sim >= ukb_age_min & adni_sim <= ukb_age_max
  ]
  adni_mean <- adni_age_mean
  adni_label <- "ADNI Comparable"

  combined <- rbind(
    data.table(Age = ukb_sim, Sample = "UKB Normative"),
    data.table(Age = adni_sim, Sample = adni_label)
  )
  combined[, Sample := factor(
    Sample, levels = c("UKB Normative", adni_label)
  )]

  colors <- manuscript_colors()

  ggplot(
    combined,
    aes(x = Age, fill = Sample, color = Sample)
  ) +
    geom_density(alpha = 0.4, linewidth = 0.6) +
    scale_fill_manual(values = setNames(
      c(colors$ukb, colors$adni),
      c("UKB Normative", adni_label)
    )) +
    scale_color_manual(values = setNames(
      c(colors$ukb, colors$adni),
      c("UKB Normative", adni_label)
    )) +
    geom_vline(
      xintercept = ukb_age_mean,
      color = colors$ukb,
      linetype = "dashed", linewidth = 0.8
    ) +
    geom_vline(
      xintercept = adni_mean,
      color = colors$adni,
      linetype = "dashed", linewidth = 0.8
    ) +
    annotate(
      "text", x = ukb_age_mean - 1, y = Inf,
      label = sprintf("%.1f", ukb_age_mean),
      color = colors$ukb,
      hjust = 1, vjust = 1.5, size = 2.5
    ) +
    annotate(
      "text", x = adni_mean + 1, y = Inf,
      label = sprintf("%.1f", adni_mean),
      color = colors$adni,
      hjust = 0, vjust = 1.5, size = 2.5
    ) +
    scale_y_continuous(expand = expansion(
      mult = c(0.02, 0.12)
    )) +
    labs(
      x = "Age (years)", y = "Density",
      fill = "Sample", color = "Sample"
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.size = unit(0.5, "lines")
    )
}

# =============================================================================
# CVR Comparison Plotting Functions (FRS vs CVR_mimic)
# =============================================================================

#' Plot CVR_mimic validation: correlation with age
#'
#' Shows that CVR_mimic is orthogonal to age (r ≈ 0) by construction,
#' contrasting with FRS which is strongly age-correlated.
#'
#' @param data.dt Data.table with CVR_mimic and AGE columns
#' @param frs_data.dt Optional: separate data.table with FRS for comparison
#' @return ggplot object
#' @export
plot_cvr_validation <- function(data.dt) {
  if (is.null(data.dt)) return(NULL)

  dt <- as.data.table(data.dt)
  colors <- manuscript_colors()

  if ("YRS_from_bl" %in% names(dt)) {
    baseline.dt <- dt[
      YRS_from_bl == 0 | !duplicated(PTID)
    ]
  } else {
    baseline.dt <- dt[!duplicated(PTID)]
  }

  cvr_age_cor <- cor(
    baseline.dt$CVR_mimic, baseline.dt$AGE,
    use = "complete"
  )

  ggplot(
    baseline.dt, aes(x = AGE, y = CVR_mimic)
  ) +
    geom_point(
      shape = 1, alpha = 0.3, size = 1,
      color = colors$cvr_pred
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = colors$cvr_pred,
      fill = colors$cvr_pred,
      alpha = 0.15, linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray30", linewidth = 0.6
    ) +
    annotate(
      "label",
      x = max(baseline.dt$AGE, na.rm = TRUE) - 1,
      y = min(
        baseline.dt$CVR_mimic, na.rm = TRUE
      ) * 0.95,
      label = sprintf("r = %.3f", cvr_age_cor),
      size = 2.5, fontface = "bold", hjust = 1,
      fill = alpha("white", 0.6),
      label.size = 0.3
    ) +
    labs(
      x = "Age (years)",
      y = expression(CVR[mimic] ~ "(age-adjusted)")
    ) +
    theme_publication(base_size = 10)
}

#' Plot FRS vs Age scatter
#'
#' @param data.dt Data with FRS_pct and AGE columns
#' @return ggplot object
plot_frs_validation <- function(data.dt) {
  if (is.null(data.dt)) return(NULL)

  dt <- as.data.table(data.dt)
  colors <- manuscript_colors()

  if ("YRS_from_bl" %in% names(dt)) {
    baseline.dt <- dt[
      YRS_from_bl == 0 | !duplicated(PTID)
    ]
  } else {
    baseline.dt <- dt[!duplicated(PTID)]
  }

  frs_age_cor <- cor(
    baseline.dt$FRS_pct, baseline.dt$AGE,
    use = "complete"
  )

  ggplot(
    baseline.dt, aes(x = AGE, y = FRS_pct)
  ) +
    geom_point(
      shape = 1, alpha = 0.3, size = 1,
      color = colors$frs_pred
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = colors$frs_pred,
      fill = colors$frs_pred,
      alpha = 0.15, linewidth = 0.8
    ) +
    annotate(
      "label",
      x = max(baseline.dt$AGE, na.rm = TRUE) - 1,
      y = min(
        baseline.dt$FRS_pct, na.rm = TRUE
      ) + 1,
      label = sprintf("r = %.3f", frs_age_cor),
      size = 2.5, fontface = "bold", hjust = 1,
      fill = alpha("white", 0.6),
      label.size = 0.3
    ) +
    labs(
      x = "Age (years)",
      y = "Framingham Risk Score (%)"
    ) +
    theme_publication(base_size = 10)
}

#' Plot FRS vs CVR_mimic comparison forest plot
#'
#' Side-by-side comparison of effect sizes from FRS and CVR_mimic analyses.
#' Highlights the key finding: FRS shows significant a-path, CVR_mimic does not.
#'
#' @param frs_results Data.table with FRS results (Domain, Beta, SE, p)
#' @param cvr_results Data.table with CVR_mimic results (Domain, Beta, SE, p)
#' @param effect_type Label for effect type
#'   (e.g., "a-path", "3-way interaction")
#' @return ggplot object
#' @export
plot_frs_cvr_comparison <- function(frs_results, cvr_results,
                                     effect_type = "Effect") {
  if (is.null(frs_results) || is.null(cvr_results)) return(NULL)

  colors <- manuscript_colors()

  # Prepare data
  frs_dt <- as.data.table(copy(frs_results))
  cvr_dt <- as.data.table(copy(cvr_results))

  frs_dt[, Predictor := "FRS"]
  cvr_dt[, Predictor := "CVR_mimic"]

  # Ensure consistent columns
  keep_cols <- c("Domain", "Beta", "SE", "Predictor")
  if ("p" %in% names(frs_dt)) {
    frs_dt[, Significant := p < 0.05]
    cvr_dt[, Significant := p < 0.05]
    keep_cols <- c(keep_cols, "Significant")
  }

  comp_dt <- rbind(
    frs_dt[, .SD, .SDcols = keep_cols],
    cvr_dt[, .SD, .SDcols = keep_cols]
  )

  # Map domain codes to full names
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  if (any(comp_dt$Domain %in% names(domain_map))) {
    comp_dt[, Domain := domain_map[Domain]]
  }

  # Add CIs
  comp_dt[, CI_lower := Beta - 1.96 * SE]
  comp_dt[, CI_upper := Beta + 1.96 * SE]

  pred_colors <- c(
    "FRS" = colors$frs_pred,
    "CVR_mimic" = colors$cvr_pred
  )

  p <- ggplot(
    comp_dt,
    aes(x = Beta, y = Domain, color = Predictor)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                  width = 0.2,
                  position = position_dodge(width = 0.5),
                  linewidth = 0.8) +
    geom_point(size = 3, position = position_dodge(width = 0.5))

  # Add significance shapes if available
  if ("Significant" %in% names(comp_dt)) {
    p <- p +
      geom_point(aes(shape = Significant),
                 size = 3,
                 position = position_dodge(width = 0.5)) +
      scale_shape_manual(
        values = c("TRUE" = 16, "FALSE" = 1),
        labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
        name = ""
      )
  }

  p <- p +
    scale_color_manual(values = pred_colors, name = "Predictor") +
    labs(
      x = expression("Effect Size (" * beta * ") with 95% CI"),
      y = NULL
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  p
}

#' Create side-by-side mediation comparison panel
#'
#' Two-panel figure comparing FRS and CVR_mimic mediation diagrams.
#'
#' @param frs_result FRS mediation result
#' @param cvr_result CVR_mimic mediation result
#' @param domain Domain label (e.g., "Memory")
#' @return patchwork combined plot
#' @export
plot_mediation_comparison <- function(frs_result, cvr_result, domain = NULL) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package required for mediation comparison")
    return(NULL)
  }

  p_frs <- plot_mediation_diagram(frs_result,
                                   domain_label = "FRS",
                                   predictor_label = "FRS")
  p_cvr <- plot_mediation_diagram(cvr_result,
                                   domain_label = "CVR_mimic",
                                   predictor_label = "CVR")

  if (is.null(p_frs) || is.null(p_cvr)) return(NULL)

  combined <- p_frs + p_cvr +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
      title = if (!is.null(domain)) domain else NULL,
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )

  combined
}

#' Create FRS ceiling effect visualization
#'
#' Shows the FRS distribution in elderly cohorts with ceiling annotation.
#'
#' @param data.dt Data.table with FRS column
#' @return ggplot object
#' @export
plot_frs_ceiling <- function(data.dt) {
  if (is.null(data.dt)) return(NULL)
  colors <- manuscript_colors()
  dt <- as.data.table(data.dt)

  # Use raw FRS if available (FRS may be z-scored)
  frs_col <- if ("FRS_raw" %in% names(dt)) {
    "FRS_raw"
  } else {
    "FRS"
  }

  # Get baseline only
  if ("YRS_from_bl" %in% names(dt)) {
    baseline.dt <- dt[
      YRS_from_bl == 0 | !duplicated(PTID)
    ]
  } else {
    baseline.dt <- dt[!duplicated(PTID)]
  }

  ceiling_pct <- mean(
    baseline.dt[[frs_col]] == 30, na.rm = TRUE
  ) * 100

  p <- ggplot(
    baseline.dt, aes(x = .data[[frs_col]])
  ) +
    geom_histogram(
      bins = 30, fill = colors$frs_hist,
      alpha = 0.7, color = "white"
    ) +
    geom_vline(
      xintercept = 30, linetype = "dashed",
      color = colors$ceiling, linewidth = 1
    ) +
    annotate(
      "text", x = 28, y = Inf,
      label = sprintf(
        "%.0f%% at ceiling (30)", ceiling_pct
      ),
      color = colors$ceiling, hjust = 1,
      vjust = 1.5, size = 3
    ) +
    labs(
      x = "Framingham Risk Score (%)",
      y = "Count"
    ) +
    theme_publication(base_size = 10)

  p
}

#' FRS ceiling histogram by sex
#'
#' Shows FRS distributions overlaid for males and
#' females with sex colors and ceiling line.
#'
#' @param data.dt Data.table with FRS and SEX columns
#' @return ggplot object
#' @export
plot_frs_ceiling_by_sex <- function(data.dt) {
  if (is.null(data.dt)) return(NULL)

  dt <- as.data.table(data.dt)
  colors <- manuscript_colors()

  if ("YRS_from_bl" %in% names(dt)) {
    bl.dt <- dt[
      YRS_from_bl == 0 | !duplicated(PTID)
    ]
  } else {
    bl.dt <- dt[!duplicated(PTID)]
  }

  stopifnot("FRS_pct" %in% names(bl.dt))
  bl.dt[, SEX := factor(
    SEX,
    levels = c("Female", "Male"),
    labels = c("Females", "Males")
  )]

  # Ceiling pct per sex
  ceil.dt <- bl.dt[, .(
    ceil_pct = mean(
      FRS_pct == 30, na.rm = TRUE
    ) * 100
  ), by = SEX]
  ceil.dt[, label := sprintf(
    "%.0f%% at ceiling", ceil_pct
  )]

  ggplot(bl.dt, aes(x = FRS_pct, fill = SEX)) +
    geom_histogram(
      bins = 30, alpha = 0.7, color = "white"
    ) +
    geom_vline(
      xintercept = 30, linetype = "dashed",
      color = "gray30", linewidth = 1
    ) +
    geom_text(
      data = ceil.dt,
      aes(x = 28, y = Inf, label = label),
      color = "black", hjust = 1, vjust = 2.5,
      size = 2.7, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c(
        Males = colors$male,
        Females = colors$female
      )
    ) +
    facet_wrap(~SEX, scales = "free_y") +
    labs(
      x = "Framingham Risk Score (%)",
      y = "Count"
    ) +
    theme_publication(base_size = 10) +
    theme(legend.position = "none")
}

#' Create interaction plot for CVR x HVR effects
#'
#' CVR_mimic version of interaction plot, showing cognitive trajectories
#' by HVR group, faceted by CVR_mimic level.
#'
#' @param cohort_dt Data.table with cognitive, HVR, and CVR_mimic data
#' @param outcome_var Name of cognitive outcome variable (default: "MEM")
#' @return ggplot object
#' @export
plot_interaction_cvr_hvr <- function(cohort_dt, outcome_var = "MEM") {
  if (is.null(cohort_dt) || nrow(cohort_dt) == 0) return(NULL)

  dt <- copy(cohort_dt)
  colors <- manuscript_colors()

  # Get baseline for median splits
  if ("YRS_from_bl" %in% names(dt)) {
    baseline <- dt[YRS_from_bl == 0]
  } else if ("VISIT" %in% names(dt)) {
    baseline <- dt[VISIT == 1]
  } else {
    baseline <- dt[!duplicated(PTID)]
  }

  hvr_median <- median(baseline$HVR_z, na.rm = TRUE)
  cvr_median <- median(baseline$CVR_mimic, na.rm = TRUE)

  # Assign groups at baseline
  baseline_groups <- baseline[, .(
    PTID,
    HVR_group = factor(
      fifelse(HVR_z > hvr_median, "Preserved HVR", "Low HVR"),
      levels = c("Preserved HVR", "Low HVR")
    ),
    CVR_group = factor(
      fifelse(CVR_mimic > cvr_median, "High CVR", "Low CVR"),
      levels = c("Low CVR", "High CVR")
    )
  )]

  dt <- merge(dt, baseline_groups, by = "PTID", all.x = TRUE)

  # Determine time variable
  time_var <- if ("YRS_from_bl" %in% names(dt)) "YRS_from_bl" else "YRS"
  if (!time_var %in% names(dt)) {
    warning("Time variable not found in data")
    return(NULL)
  }

  # Aggregate trajectories
  traj_summary <- dt[!is.na(get(outcome_var)) & !is.na(HVR_group), .(
    Mean_Cog = mean(get(outcome_var), na.rm = TRUE),
    SE_Cog = sd(get(outcome_var), na.rm = TRUE) / sqrt(.N),
    N = .N
  ), by = .(HVR_group, CVR_group, YRS_bin = round(get(time_var)))]

  traj_summary <- traj_summary[N >= 10]

  if (nrow(traj_summary) == 0) {
    warning("No data after filtering for minimum cell size")
    return(NULL)
  }

  # CVR = color, HVR = linetype (consistent)
  cvr_colors <- c(
    "Low CVR" = colors$cvr_low,
    "High CVR" = colors$cvr_high
  )
  hvr_ltypes <- c(
    "Preserved HVR" = "solid",
    "Low HVR" = "dashed"
  )

  p <- ggplot(
    traj_summary,
    aes(
      x = YRS_bin, y = Mean_Cog,
      color = CVR_group,
      linetype = HVR_group
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Mean_Cog - SE_Cog,
        ymax = Mean_Cog + SE_Cog,
        fill = CVR_group
      ),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(
      aes(shape = HVR_group),
      size = 2, fill = "white", stroke = 0.8
    ) +
    scale_color_manual(
      values = cvr_colors,
      name = "CVR (Baseline)"
    ) +
    scale_fill_manual(
      values = cvr_colors, guide = "none"
    ) +
    scale_linetype_manual(
      values = hvr_ltypes,
      name = "HVR (Baseline)"
    ) +
    scale_shape_manual(
      values = c(
        "Preserved HVR" = 16,
        "Low HVR" = 21
      ),
      name = "HVR (Baseline)"
    ) +
    labs(
      x = "Years from Baseline",
      y = switch(
        outcome_var,
        MEM = "Memory Composite",
        LAN = "Language Composite",
        EXF = "Executive Composite",
        paste(outcome_var, "Score")
      )
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(1.5, "cm")
    )

  p
}

#' Create LME forest comparison plot (FRS vs CVR)
#'
#' Builds a two-panel (Males/Females) forest plot comparing
#' 3-way interaction effect sizes from FRS and CVR LME analyses.
#' Returns a patchwork composition with panel tags A, B.
#'
#' @param frs_res LME results list for FRS predictor
#' @param cvr_res LME results list for CVR predictor
#' @return patchwork object
#' @export
plot_lme_forest_comparison <- function(frs_res, cvr_res) {
  if (is.null(frs_res) || is.null(cvr_res)) {
    return(NULL)
  }
  colors <- manuscript_colors()

  extract_stratified.fn <- function(res, pred) {
    if (is.null(res$stratified)) {
      return(data.table())
    }
    lapply(c("Male", "Female"), function(sx) {
      lapply(c("MEM", "LAN", "EXF"), function(d) {
        r <- res$stratified[[sx]][[d]]
        if (!is.null(r) && !is.null(r$interaction)) {
          data.table(
            Domain = d, Sex = sx,
            Beta = r$interaction$beta,
            SE = r$interaction$se,
            p = r$interaction$p_value,
            Predictor = pred
          )
        }
      }) |> rbindlist()
    }) |> rbindlist()
  }

  frs_forest.dt <- extract_stratified.fn(
    frs_res, "FRS"
  )
  cvr_forest.dt <- extract_stratified.fn(
    cvr_res, "CVR_MIMIC"
  )

  if (nrow(frs_forest.dt) == 0 ||
      nrow(cvr_forest.dt) == 0) {
    return(NULL)
  }

  forest.dt <- rbind(frs_forest.dt, cvr_forest.dt)
  forest.dt[, Significant := p < 0.05]

  domain_map <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )
  forest.dt[, Domain := domain_map[Domain]]
  forest.dt[, CI_lower := Beta - 1.96 * SE]
  forest.dt[, CI_upper := Beta + 1.96 * SE]

  # Facet label with "s" suffix
  forest.dt[, Sex_label := paste0(Sex, "s")]
  forest.dt[, Sex_label := factor(
    Sex_label, levels = c("Females", "Males")
  )]

  pred_labels <- c(
    "FRS" = "FRS",
    "CVR_MIMIC" = expression(CVR[mimic])
  )
  pred_colors <- c(
    "FRS" = colors$frs_pred,
    "CVR_MIMIC" = colors$cvr_pred
  )

  forest.dt[, Predictor := factor(
    Predictor, levels = names(pred_colors)
  )]
  forest.dt[, Domain := factor(
    Domain,
    levels = rev(c(
      "Memory", "Language",
      "Executive Function"
    ))
  )]

  forest.dt[, point_fill := fifelse(
    Significant,
    pred_colors[as.character(Predictor)],
    "white"
  )]

  pd <- position_dodge(width = 0.5)

  ggplot(
    forest.dt,
    aes(
      x = Beta, y = Domain,
      group = Predictor,
      color = Predictor
    )
  ) +
    geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "gray50"
    ) +
    geom_errorbarh(
      aes(xmin = CI_lower, xmax = CI_upper),
      height = 0, position = pd,
      linewidth = 0.5
    ) +
    geom_point(
      aes(fill = point_fill),
      shape = 21, size = 2.5, stroke = 0.6,
      position = pd
    ) +
    scale_fill_identity() +
    scale_color_manual(
      values = pred_colors,
      labels = pred_labels,
      name = "Predictor",
      drop = FALSE
    ) +
    facet_wrap(~ Sex_label, ncol = 2) +
    labs(
      x = expression(beta * " (95% CI)"),
      y = NULL
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9)
    )
}

# =========================================================
# New / Improved Plot Functions (Phase 1b)
# =========================================================

#' Marginal effects of brain on cognitive slope
#' at different CVR levels
#'
#' Replaces the noisy spaghetti interaction plots
#' with clean model-predicted marginal effects.
#'
#' @param model.fit lmerMod or lmerModLmerTest
#' @param data.dt Analysis cohort data.table
#' @param domain Cognitive domain (MEM/LAN/EXF)
#' @param sex Sex label for title
#' @param cvr_var CVR column name (FRS or CVR_mimic)
#' @param brain_var Brain column name (HVR_z)
#' @return ggplot object
plot_marginal_effects <- function(
    model.fit, data.dt, domain, sex,
    cvr_var = "FRS", brain_var = "HVR_z") {
  if (is.null(model.fit) ||
      is.null(data.dt)) return(NULL)

  # Unwrap model if nested
  if (is.list(model.fit) &&
      "model" %in% names(model.fit)) {
    model.fit <- model.fit$model
  }
  if (is.list(model.fit) &&
      "fit" %in% names(model.fit)) {
    model.fit <- model.fit$fit
  }
  if (is.null(model.fit)) return(NULL)

  colors <- manuscript_colors()
  dt <- as.data.table(data.dt)

  # Compute baseline age if needed
  if (!"Age_bl" %in% names(dt) &&
      "AGE" %in% names(dt)) {
    bl_ages.dt <- dt[
      YRS_from_bl == 0,
      .(Age_bl = AGE[1]),
      by = PTID
    ]
    dt <- merge(
      dt, bl_ages.dt, by = "PTID", all.x = TRUE
    )
  }

  # CVR tertile cutoffs from baseline
  bl.dt <- dt[, .SD[1], by = PTID]
  cvr_vals.v <- bl.dt[[cvr_var]]
  q_cuts.v <- quantile(
    cvr_vals.v, c(1 / 3, 2 / 3), na.rm = TRUE
  )

  max_yrs <- max(
    dt$YRS_from_bl, na.rm = TRUE
  )
  yrs_seq.v <- seq(0, max_yrs, length.out = 40)

  # Prediction grid: 3 CVR levels
  cvr_levels.v <- c(
    Low = median(
      cvr_vals.v[cvr_vals.v <= q_cuts.v[1]],
      na.rm = TRUE
    ),
    Medium = median(cvr_vals.v, na.rm = TRUE),
    High = median(
      cvr_vals.v[cvr_vals.v >= q_cuts.v[2]],
      na.rm = TRUE
    )
  )

  brain_mean <- mean(
    bl.dt[[brain_var]], na.rm = TRUE
  )
  age_mean <- mean(bl.dt$Age_bl, na.rm = TRUE)
  educ_mean <- as.integer(
    mean(bl.dt$EDUC, na.rm = TRUE)
  )

  grids.lst <- list()
  for (lvl in names(cvr_levels.v)) {
    g.dt <- data.table(
      YRS_from_bl = yrs_seq.v,
      Age_bl = age_mean,
      EDUC = educ_mean,
      APOE4 = FALSE
    )
    g.dt[[cvr_var]] <- cvr_levels.v[[lvl]]
    g.dt[[brain_var]] <- brain_mean
    g.dt[, CVR_Level := lvl]
    grids.lst[[lvl]] <- g.dt
  }
  pred.dt <- rbindlist(grids.lst)

  tryCatch({
    pred.dt[, predicted := predict(
      model.fit,
      newdata = pred.dt,
      re.form = NA,
      allow.new.levels = TRUE
    )]
  }, error = function(e) {
    return(NULL)
  })

  if (!"predicted" %in% names(pred.dt)) {
    return(NULL)
  }

  pred.dt[, CVR_Level := factor(
    CVR_Level,
    levels = c("Low", "Medium", "High")
  )]

  DOMAIN_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )
  dom_label <- if (
    domain %in% names(DOMAIN_LABELS)
  ) DOMAIN_LABELS[domain] else domain

  cvr_label <- if (
    cvr_var == "CVR_mimic"
  ) "CVR" else cvr_var

  pal.v <- c(
    Low = "#1A9850",
    Medium = "gray50",
    High = "#D73027"
  )

  ggplot(
    pred.dt,
    aes(
      x = YRS_from_bl,
      y = predicted,
      color = CVR_Level
    )
  ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = pal.v,
      name = paste(cvr_label, "Level")
    ) +
    labs(
      x = "Years from Baseline",
      y = paste(dom_label, "Composite"),
      title = paste(sex, "\u2013", dom_label)
    ) +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom")
}

#' Coupling scatter: HVR slope vs Cognitive slope
#'
#' Scatter of LGCM-estimated individual slopes
#' colored by CVR tertile.
#'
#' @param lgcm_parallel.res lgcm_parallel_results.rds
#' @param sex "Male" or "Female"
#' @param domain "MEM", "LAN", or "EXF"
#' @return ggplot object
plot_coupling_scatter <- function(
    lgcm_parallel.res, sex, domain) {
  if (is.null(lgcm_parallel.res)) return(NULL)
  key <- paste(sex, domain, sep = "_")
  res <- lgcm_parallel.res$linear_results[[key]]
  if (is.null(res)) return(NULL)

  colors <- manuscript_colors()

  DOMAIN_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )
  dom_label <- if (
    domain %in% names(DOMAIN_LABELS)
  ) DOMAIN_LABELS[domain] else domain

  # Extract params summary for annotation
  cpl <- res$coupling

  # Build label
  cov_label <- sprintf(
    "Coupling cov = %.4f\np = %s",
    cpl$cov,
    format.pval(cpl$p, digits = 2)
  )

  # We cannot extract individual scores from
  # OpenMx summary, so use the summary stats
  # and plot the group-level result as annotated
  pe.dt <- as.data.table(res$params)

  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = paste(
        sex, "\u2013", dom_label,
        "\n\n", cov_label
      ),
      size = 3.5, hjust = 0.5
    ) +
    labs(
      title = paste(
        sex, "\u2013", dom_label,
        "Slope Coupling"
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 10
      )
    )
}

#' SEM-style path diagram for parallel process
#'
#' Draws a ggplot-based path diagram showing
#' parallel process model with coupling.
#'
#' @param lgcm_parallel.res lgcm_parallel_results.rds
#' @param sex "Male" or "Female"
#' @param domain "MEM", "LAN", or "EXF"
#' @return ggplot object
plot_coupling_path_diagram <- function(
    lgcm_parallel.res, sex, domain) {
  if (is.null(lgcm_parallel.res)) return(NULL)
  key <- paste(sex, domain, sep = "_")
  res <- lgcm_parallel.res$linear_results[[key]]
  if (is.null(res)) return(NULL)

  DOMAIN_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )
  dom_label <- if (
    domain %in% names(DOMAIN_LABELS)
  ) DOMAIN_LABELS[domain] else domain

  cpl <- res$coupling
  mod <- res$moderation

  # Node positions
  nodes.dt <- data.table(
    label = c(
      "HVR\nIntercept", "HVR\nSlope",
      paste0(dom_label, "\nIntercept"),
      paste0(dom_label, "\nSlope"),
      "FRS"
    ),
    x = c(0, 0, 2, 2, 1),
    y = c(1, 0, 1, 0, -0.6)
  )

  # Significance formatting
  sig_fmt.fn <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else ""
  }

  cpl_label <- sprintf(
    "cov = %.4f%s",
    cpl$cov, sig_fmt.fn(cpl$p)
  )
  mod_label <- sprintf(
    "\u03b2 = %.4f%s",
    mod$est, sig_fmt.fn(mod$p)
  )

  cpl_color <- if (
    !is.na(cpl$p) && cpl$p < 0.05
  ) "firebrick" else "gray60"

  mod_color <- if (
    !is.na(mod$p) && mod$p < 0.05
  ) "steelblue" else "gray60"

  p <- ggplot() +
    # Coupling arrow (hvr_s ~~ cog_s)
    geom_segment(
      aes(x = 0.15, y = 0, xend = 1.85, yend = 0),
      arrow = arrow(
        length = unit(0.1, "cm"),
        ends = "both"
      ),
      linewidth = 1.2,
      color = cpl_color
    ) +
    annotate(
      "text", x = 1, y = 0.12,
      label = cpl_label,
      size = 2.7, color = cpl_color,
      fontface = "bold"
    ) +
    # FRS x HVR -> cog_s
    geom_segment(
      aes(
        x = 1, y = -0.45,
        xend = 1.85, yend = -0.05
      ),
      arrow = arrow(length = unit(0.1, "cm")),
      linewidth = 0.8,
      color = mod_color, linetype = "dashed"
    ) +
    annotate(
      "text", x = 1.6, y = -0.35,
      label = mod_label,
      size = 2.5, color = mod_color
    ) +
    # Nodes
    geom_point(
      data = nodes.dt,
      aes(x = x, y = y),
      size = 16, shape = 21,
      fill = "white", color = "black",
      stroke = 1
    ) +
    geom_text(
      data = nodes.dt,
      aes(x = x, y = y, label = label),
      size = 2.3, fontface = "bold"
    ) +
    coord_cartesian(
      xlim = c(-0.5, 2.5),
      ylim = c(-0.8, 1.3)
    ) +
    labs(
      title = paste(
        sex, "\u2013", dom_label
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 10
      )
    )

  p
}

#' Z-score validation forest plot (HVR + HC)
#'
#' Combined forest plot showing both HVR and HC mean
#' z-scores (+/- 95% CI) by age bin and overall.
#'
#' @param nv.res normative_transfer_validation.rds
#' @return ggplot object
plot_zscore_validation_forest <- function(nv.res) {
  if (is.null(nv.res)) return(NULL)

  cs <- nv.res$comparable_subsample
  if (is.null(cs)) return(NULL)

  colors <- manuscript_colors()
  measure_cols <- c(
    HVR = colors$hvr_measure,
    HC = colors$hc_measure
  )

  build_rows.fn <- function(
    strat.dt, overall, measure
  ) {
    rows.lst <- list()
    if (!is.null(strat.dt) && nrow(strat.dt) > 0) {
      strat.dt[,
        bin_label := as.character(age_bin)
      ]
      for (i in seq_len(nrow(strat.dt))) {
        rows.lst[[i]] <- data.table(
          Group = strat.dt$bin_label[i],
          Mean = strat.dt$mean_z[i],
          CI_lower = strat.dt$ci_lower[i],
          CI_upper = strat.dt$ci_upper[i],
          N = strat.dt$n[i],
          Type = "Age Bin",
          Measure = measure
        )
      }
    }
    if (!is.null(overall)) {
      rows.lst[["overall"]] <- data.table(
        Group = "Overall",
        Mean = overall$mean,
        CI_lower = overall$ci_95[["lower"]],
        CI_upper = overall$ci_95[["upper"]],
        N = overall$n,
        Type = "Overall",
        Measure = measure
      )
    }
    rows.lst
  }

  rows.lst <- c(
    build_rows.fn(
      as.data.table(cs$stratified_hvr),
      cs$HVR, "HVR"
    ),
    build_rows.fn(
      as.data.table(cs$stratified_hc),
      cs$HC, "HC"
    )
  )

  if (length(rows.lst) == 0) return(NULL)
  forest.dt <- rbindlist(rows.lst)

  # Order: age bins then overall (shared y-axis)
  grp_levels <- unique(forest.dt$Group)
  overall_idx <- which(grp_levels == "Overall")
  if (length(overall_idx) > 0) {
    grp_levels <- c(
      setdiff(grp_levels, "Overall"), "Overall"
    )
  }
  forest.dt[, Group := factor(
    Group, levels = rev(grp_levels)
  )]
  forest.dt[, Measure := factor(
    Measure, levels = c("HVR", "HC")
  )]

  # Dodge offset for two measures
  pd <- position_dodge(width = 0.5)

  ggplot(
    forest.dt,
    aes(
      x = Mean, y = Group,
      xmin = CI_lower, xmax = CI_upper,
      color = Measure,
      linetype = Type
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed", color = "gray50"
    ) +
    geom_errorbarh(
      aes(xmin = CI_lower, xmax = CI_upper),
      height = 0, position = pd
    ) +
    # Age-bin points: white fill, drawn over bars
    geom_point(
      data = forest.dt[Type == "Age Bin"],
      aes(fill = Measure),
      shape = 21, size = 2,
      position = pd, stroke = 0.6,
      fill = "white"
    ) +
    # Overall points: solid fill
    geom_point(
      data = forest.dt[Type == "Overall"],
      aes(fill = Measure),
      shape = 21, size = 2.5,
      position = pd, stroke = 0.6
    ) +
    scale_color_manual(
      values = measure_cols,
      labels = c(HVR = "HVR", HC = "HC")
    ) +
    scale_fill_manual(
      values = measure_cols,
      labels = c(HVR = "HVR", HC = "HC")
    ) +
    scale_linetype_manual(
      values = c(
        "Age Bin" = "solid",
        "Overall" = "solid"
      ),
      guide = "none"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 16, size = 3,
          linetype = "blank"
        )
      ),
      fill = "none"
    ) +
    labs(x = "Z-score", y = NULL) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
}

#' FRS variance decomposition bar chart
#'
#' Grouped bar: Age-unique, CVRF-unique, Shared R2.
#'
#' @param frs_problem.res frs_problem_analysis.rds
#' @return ggplot object
plot_frs_variance_decomposition <- function(
    frs_problem.res) {
  if (is.null(frs_problem.res) ||
      is.null(
        frs_problem.res$variance_decomposition
      )) {
    return(NULL)
  }

  vd <- frs_problem.res$variance_decomposition
  colors <- manuscript_colors()

  # Decompose shared variance (handle naming variants)
  r2_age <- vd$r2_age %||%
    vd$r2_age_only %||%
    vd$r_squared_age_only
  r2_cvrf <- vd$r2_cvrf %||%
    vd$r2_cvrf_only %||%
    vd$r_squared_cvrf_only
  r2_full <- vd$r2_full %||%
    vd$r_squared_full

  shared <- r2_age + r2_cvrf - r2_full
  age_unique <- r2_age - shared
  cvrf_unique <- r2_cvrf - shared

  # Clamp negatives
  shared <- max(shared, 0)
  age_unique <- max(age_unique, 0)
  cvrf_unique <- max(cvrf_unique, 0)

  bar.dt <- data.table(
    Component = factor(
      c(
        "Age (unique)",
        "CV Risk Factors (unique)",
        "Shared"
      ),
      levels = c(
        "Age (unique)",
        "CV Risk Factors (unique)",
        "Shared"
      )
    ),
    R2 = c(age_unique, cvrf_unique, shared)
  )

  ggplot(
    bar.dt, aes(x = Component, y = R2, fill = Component)
  ) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = sprintf("%.1f%%", R2 * 100)),
      vjust = -0.5, size = 3
    ) +
    scale_fill_manual(
      values = c(
        "Age (unique)" = colors$age_color,
        "CV Risk Factors (unique)" =
          colors$cvrf_color,
        "Shared" = colors$summary_color
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = NULL,
      y = expression(R^2),
      title = paste(
        "FRS Variance Decomposition:",
        "Age vs CV Risk Factors"
      )
    ) +
    theme_publication(base_size = 10)
}

#' Simulation comparison plot (bias + RMSE)
#'
#' Two-panel: (A) Bias by N, (B) RMSE by N.
#'
#' @param simulation.res simulation_results.rds
#' @return ggplot object
plot_simulation_comparison <- function(
    simulation.res) {
  if (is.null(simulation.res) ||
      is.null(simulation.res$summary)) {
    return(NULL)
  }

  sim.dt <- as.data.table(simulation.res$summary)
  colors <- manuscript_colors()

  # Reshape for plotting
  bias.dt <- melt(
    sim.dt,
    id.vars = "n",
    measure.vars = c("bias_se", "bias_free"),
    variable.name = "Method",
    value.name = "Bias"
  )
  bias.dt[, Method := fifelse(
    Method == "bias_se",
    "SE-Constrained", "Free Loading"
  )]

  rmse.dt <- melt(
    sim.dt,
    id.vars = "n",
    measure.vars = c("rmse_se", "rmse_free"),
    variable.name = "Method",
    value.name = "RMSE"
  )
  rmse.dt[, Method := fifelse(
    Method == "rmse_se",
    "SE-Constrained", "Free Loading"
  )]

  pal.v <- c(
    "SE-Constrained" = colors$method_se,
    "Free Loading" = colors$method_free
  )

  p_bias <- ggplot(
    bias.dt,
    aes(
      x = factor(n), y = Bias,
      fill = Method
    )
  ) +
    geom_col(
      position = position_dodge(0.7),
      width = 0.6
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed", color = "gray50"
    ) +
    scale_fill_manual(values = pal.v) +
    labs(x = "N", y = "Bias") +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom")

  p_rmse <- ggplot(
    rmse.dt,
    aes(
      x = factor(n), y = RMSE,
      fill = Method
    )
  ) +
    geom_col(
      position = position_dodge(0.7),
      width = 0.6
    ) +
    scale_fill_manual(values = pal.v) +
    labs(x = "N", y = "RMSE") +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom")

  (p_bias | p_rmse) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}

#' Improved mediation path diagram
#'
#' Mediation path diagram with significance
#' color coding and bootstrap CI.
#'
#' @param mediation_result Mediation result list
#' @param domain_label Domain label
#' @param predictor_label CVR predictor label
#' @return ggplot object
plot_mediation_diagram <- function(
    mediation_result, domain_label = NULL,
    predictor_label = "FRS") {
  if (is.null(mediation_result) ||
      !isTRUE(mediation_result$converged)) {
    return(NULL)
  }

  r <- mediation_result

  # Significance color helper
  path_color.fn <- function(p) {
    if (is.na(p)) return("gray40")
    if (p < 0.05) "firebrick" else "gray40"
  }

  sig_fmt.fn <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else ""
  }

  # Format predictor label for display
  disp_label <- if (
    grepl("CVR", predictor_label)
  ) "CVR\nMIMIC" else predictor_label

  # Nodes
  nodes.dt <- data.table(
    name = c(
      disp_label,
      "HVR\nSlope",
      "Cog\nSlope"
    ),
    x = c(0, 1, 2),
    y = c(0, 0.35, 0)
  )

  a_col <- path_color.fn(r$a$p)
  b_col <- path_color.fn(r$b$p)
  c_col <- path_color.fn(r$cprime$p)

  a_label <- sprintf(
    "a = %.4f%s", r$a$est, sig_fmt.fn(r$a$p)
  )
  b_label <- sprintf(
    "b = %.3f%s", r$b$est, sig_fmt.fn(r$b$p)
  )
  has_cp_boot <- !is.null(r$cprime$boot_ci_lower) &&
    !is.na(r$cprime$boot_ci_lower)
  cp_label <- if (has_cp_boot) {
    sprintf(
      "c' = %.4f%s\n[%.4f, %.4f]",
      r$cprime$est, sig_fmt.fn(r$cprime$p),
      r$cprime$boot_ci_lower,
      r$cprime$boot_ci_upper
    )
  } else {
    sprintf(
      "c' = %.4f%s",
      r$cprime$est, sig_fmt.fn(r$cprime$p)
    )
  }

  # Indirect effect with bootstrap CI
  has_boot <- !is.null(r$indirect$boot_ci_lower) &&
    !is.na(r$indirect$boot_ci_lower)

  if (has_boot) {
    boot_sig <- if (
      isTRUE(r$indirect$boot_significant)
    ) "*" else ""
    ind_label <- sprintf(
      "ab = %.5f\n[%.5f, %.5f]%s",
      r$indirect$est,
      r$indirect$boot_ci_lower,
      r$indirect$boot_ci_upper,
      boot_sig
    )
    ind_col <- if (
      isTRUE(r$indirect$boot_significant)
    ) "firebrick" else "black"
  } else {
    ind_label <- sprintf(
      "ab = %.5f\n(Sobel p = %.3f)",
      r$indirect$est,
      r$indirect$sobel_p
    )
    ind_col <- if (
      !is.na(r$indirect$sobel_p) &&
        r$indirect$sobel_p < 0.05
    ) "firebrick" else "black"
  }

  # Proportion mediated
  prop_med <- NULL
  pm <- r$prop_mediated %||%
    r$proportion_mediated
  if (!is.null(pm) && !is.na(pm)) {
    pm_pct <- min(pm * 100, 100)
    prop_med <- sprintf(
      "Mediated: %.1f%%", pm_pct
    )
  }

  p <- ggplot() +
    # Path a
    geom_segment(
      aes(
        x = 0.34, y = 0.12,
        xend = 0.66, yend = 0.23
      ),
      arrow = arrow(length = unit(0.15, "cm")),
      linewidth = 0.4, color = a_col
    ) +
    annotate(
      "text", x = -0.45, y = 0.22,
      label = a_label, size = 1.8,
      fontface = "bold", color = a_col,
      hjust = 0
    ) +
    # Path b
    geom_segment(
      aes(
        x = 1.34, y = 0.23,
        xend = 1.66, yend = 0.12
      ),
      arrow = arrow(
        length = unit(0.15, "cm")
      ),
      linewidth = 0.4, color = b_col
    ) +
    annotate(
      "text", x = 2.45, y = 0.22,
      label = b_label, size = 1.8,
      fontface = "bold", color = b_col,
      hjust = 1
    ) +
    # Path c' (direct)
    geom_segment(
      aes(
        x = 0.49, y = -0.04,
        xend = 1.51, yend = -0.04
      ),
      arrow = arrow(
        length = unit(0.15, "cm")
      ),
      linewidth = 0.4,
      color = c_col
    ) +
    annotate(
      "text", x = 0.96, y = -0.12,
      label = cp_label, size = 1.8,
      fontface = "bold", color = c_col
    ) +
    # Indirect effect box
    annotate(
      "label", x = 1, y = 0.56,
      label = ind_label, size = 1.8,
      fill = "white", color = ind_col,
      fontface = "bold",
      label.padding = unit(0.3, "lines")
    ) +
    # Nodes
    geom_point(
      data = nodes.dt,
      aes(x = x, y = y),
      size = 13, shape = 21,
      fill = "white", color = "black",
      stroke = 0.8
    ) +
    geom_text(
      data = nodes.dt,
      aes(x = x, y = y, label = name),
      size = 2.2, fontface = "bold"
    ) +
    coord_cartesian(
      xlim = c(-0.8, 2.8),
      ylim = c(-0.18, 0.72)
    ) +
    labs(title = domain_label) +
    theme_void() +
    theme(
      plot.title = element_text(
        face = "bold", size = 9, hjust = 0.5
      )
    )

  # Add proportion mediated
  if (!is.null(prop_med)) {
    p <- p + annotate(
      "text", x = 1, y = 0.67,
      label = prop_med, size = 2,
      fontface = "italic", color = "black"
    )
  }

  p
}

# ----------------------------------------------------------
# Combined mediation panel (both sexes in one panel)
# ----------------------------------------------------------

#' Single mediation path panel with male + female
#'
#' @param male_result Male mediation list (from
#'   make_med_list)
#' @param female_result Female mediation list
#' @param domain_label Panel title
#' @param predictor_label "FRS" or "CVR MIMIC"
#' @return ggplot object
plot_mediation_panel.fn <- function(
    male_result, female_result,
    domain_label = NULL,
    predictor_label = "FRS") {

  if (is.null(male_result) &&
      is.null(female_result)) {
    return(NULL)
  }

  # -- helpers --
  # Significance via bootstrap CI excluding zero
  ci_sig.fn <- function(lo, hi) {
    !is.null(lo) && !is.null(hi) &&
      !is.na(lo) && !is.na(hi) &&
      (lo > 0 || hi < 0)
  }

  sig_star.fn <- function(lo, hi) {
    if (ci_sig.fn(lo, hi)) "*" else ""
  }

  is_sig.fn <- function(lo, hi) {
    ci_sig.fn(lo, hi)
  }

  # "a: 0.032*"
  fmt_est.fn <- function(
      sx, path, est, lo, hi,
      dig = 4L) {
    sprintf(
      "%s: %s%s", path,
      formatC(
        est, format = "f", digits = dig
      ),
      sig_star.fn(lo, hi)
    )
  }

  # "c': 0.032 [0.01, 0.05]*"
  fmt_ci.fn <- function(
      sx, path, est,
      lo, hi, dig = 4L) {
    if (is.null(lo) || is.na(lo)) {
      return(sprintf(
        "%s: %s", path,
        formatC(
          est, format = "f", digits = dig
        )
      ))
    }
    d.fn <- function(x) {
      formatC(
        x, format = "f", digits = dig
      )
    }
    sprintf(
      "%s: %s [%s, %s]%s",
      path,
      d.fn(est), d.fn(lo), d.fn(hi),
      sig_star.fn(lo, hi)
    )
  }

  # -- colours & styling --
  COL_M <- "midnightblue"
  COL_F <- "darkred"
  A_SIG <- 1.0
  A_NS <- 0.55
  SZ <- 1.5

  get_alpha.fn <- function(lo, hi) {
    if (is_sig.fn(lo, hi)) A_SIG else A_NS
  }

  get_face.fn <- function(lo, hi) {
    if (is_sig.fn(lo, hi)) {
      "bold"
    } else {
      "plain"
    }
  }

  # -- predictor display --
  pred_disp <- if (
    grepl("CVR", predictor_label)
  ) "CVR\nMIMIC" else predictor_label

  # -- nodes (wider triangle) --
  nodes.dt <- data.table(
    name = c(
      pred_disp, "HVR\nSlope",
      "Cog\nSlope"
    ),
    x = c(0, 1.5, 3),
    y = c(0, 0.50, 0)
  )

  m <- male_result
  f <- female_result

  # Arrow alpha (sig if either sex)
  a_arr <- if (
    is_sig.fn(
      m$a$boot_ci_lower,
      m$a$boot_ci_upper
    ) || is_sig.fn(
      f$a$boot_ci_lower,
      f$a$boot_ci_upper
    )
  ) 0.9 else 0.3
  b_arr <- if (
    is_sig.fn(
      m$b$boot_ci_lower,
      m$b$boot_ci_upper
    ) || is_sig.fn(
      f$b$boot_ci_lower,
      f$b$boot_ci_upper
    )
  ) 0.9 else 0.3
  cp_arr <- if (
    is_sig.fn(
      m$cprime$boot_ci_lower,
      m$cprime$boot_ci_upper
    ) || is_sig.fn(
      f$cprime$boot_ci_lower,
      f$cprime$boot_ci_upper
    )
  ) 0.9 else 0.3

  # -- text labels (all use boot CI) --
  # Path a (est + CI + star)
  m_a <- fmt_ci.fn(
    "M", "a", m$a$est,
    m$a$boot_ci_lower,
    m$a$boot_ci_upper
  )
  f_a <- fmt_ci.fn(
    "F", "a", f$a$est,
    f$a$boot_ci_lower,
    f$a$boot_ci_upper
  )

  # Path b (est + CI + star)
  m_b <- fmt_ci.fn(
    "M", "b", m$b$est,
    m$b$boot_ci_lower,
    m$b$boot_ci_upper,
    dig = 3L
  )
  f_b <- fmt_ci.fn(
    "F", "b", f$b$est,
    f$b$boot_ci_lower,
    f$b$boot_ci_upper,
    dig = 3L
  )

  # Path c' (est + CI + star)
  m_cp <- fmt_ci.fn(
    "M", "c'", m$cprime$est,
    m$cprime$boot_ci_lower,
    m$cprime$boot_ci_upper
  )
  f_cp <- fmt_ci.fn(
    "F", "c'", f$cprime$est,
    f$cprime$boot_ci_lower,
    f$cprime$boot_ci_upper
  )

  # Indirect ab (est + CI + star via CI)
  m_ab <- fmt_ci.fn(
    "M", "ab", m$indirect$est,
    m$indirect$boot_ci_lower,
    m$indirect$boot_ci_upper,
    dig = 5L
  )
  f_ab <- fmt_ci.fn(
    "F", "ab", f$indirect$est,
    f$indirect$boot_ci_lower,
    f$indirect$boot_ci_upper,
    dig = 5L
  )

  # -- build ggplot --
  p <- ggplot() +
    # Arrow: a (pred -> HVR)
    geom_segment(
      aes(
        x = 0.35, y = 0.12,
        xend = 1.15, yend = 0.38
      ),
      arrow = arrow(
        length = unit(0.12, "cm")
      ),
      linewidth = 0.5,
      color = alpha("gray30", a_arr)
    ) +
    # Arrow: b (HVR -> cog)
    geom_segment(
      aes(
        x = 1.85, y = 0.38,
        xend = 2.65, yend = 0.12
      ),
      arrow = arrow(
        length = unit(0.12, "cm")
      ),
      linewidth = 0.5,
      color = alpha("gray30", b_arr)
    ) +
    # Arrow: c' (pred -> cog, pushed down)
    geom_segment(
      aes(
        x = 0.45, y = -0.12,
        xend = 2.55, yend = -0.12
      ),
      arrow = arrow(
        length = unit(0.12, "cm")
      ),
      linewidth = 0.5,
      color = alpha("gray30", cp_arr)
    ) +
    # a text (upper-left, pushed out)
    annotate(
      "text", x = -0.45, y = 0.40,
      label = m_a, size = SZ,
      fontface = get_face.fn(
        m$a$boot_ci_lower,
        m$a$boot_ci_upper
      ),
      color = alpha(
        COL_M, get_alpha.fn(
          m$a$boot_ci_lower,
          m$a$boot_ci_upper
        )
      ),
      hjust = 0
    ) +
    annotate(
      "text", x = -0.45, y = 0.35,
      label = f_a, size = SZ,
      fontface = get_face.fn(
        f$a$boot_ci_lower,
        f$a$boot_ci_upper
      ),
      color = alpha(
        COL_F, get_alpha.fn(
          f$a$boot_ci_lower,
          f$a$boot_ci_upper
        )
      ),
      hjust = 0
    ) +
    # b text (upper-right, pushed out)
    annotate(
      "text", x = 3.45, y = 0.40,
      label = m_b, size = SZ,
      fontface = get_face.fn(
        m$b$boot_ci_lower,
        m$b$boot_ci_upper
      ),
      color = alpha(
        COL_M, get_alpha.fn(
          m$b$boot_ci_lower,
          m$b$boot_ci_upper
        )
      ),
      hjust = 1
    ) +
    annotate(
      "text", x = 3.45, y = 0.35,
      label = f_b, size = SZ,
      fontface = get_face.fn(
        f$b$boot_ci_lower,
        f$b$boot_ci_upper
      ),
      color = alpha(
        COL_F, get_alpha.fn(
          f$b$boot_ci_lower,
          f$b$boot_ci_upper
        )
      ),
      hjust = 1
    ) +
    # c' text (below arrow, nudged up)
    annotate(
      "text", x = 1.5, y = -0.21,
      label = m_cp, size = SZ,
      fontface = get_face.fn(
        m$cprime$boot_ci_lower,
        m$cprime$boot_ci_upper
      ),
      color = alpha(
        COL_M, get_alpha.fn(
          m$cprime$boot_ci_lower,
          m$cprime$boot_ci_upper
        )
      ),
      hjust = 0.5
    ) +
    annotate(
      "text", x = 1.5, y = -0.26,
      label = f_cp, size = SZ,
      fontface = get_face.fn(
        f$cprime$boot_ci_lower,
        f$cprime$boot_ci_upper
      ),
      color = alpha(
        COL_F, get_alpha.fn(
          f$cprime$boot_ci_lower,
          f$cprime$boot_ci_upper
        )
      ),
      hjust = 0.5
    ) +
    # ab text (above HVR, nudged down)
    annotate(
      "text", x = 1.5, y = 0.73,
      label = m_ab, size = SZ,
      fontface = get_face.fn(
        m$indirect$boot_ci_lower,
        m$indirect$boot_ci_upper
      ),
      color = alpha(
        COL_M, get_alpha.fn(
          m$indirect$boot_ci_lower,
          m$indirect$boot_ci_upper
        )
      ),
      hjust = 0.5
    ) +
    annotate(
      "text", x = 1.5, y = 0.68,
      label = f_ab, size = SZ,
      fontface = get_face.fn(
        f$indirect$boot_ci_lower,
        f$indirect$boot_ci_upper
      ),
      color = alpha(
        COL_F, get_alpha.fn(
          f$indirect$boot_ci_lower,
          f$indirect$boot_ci_upper
        )
      ),
      hjust = 0.5
    ) +
    # Nodes
    geom_point(
      data = nodes.dt,
      aes(x = x, y = y),
      size = 12, shape = 21,
      fill = "white", color = "black",
      stroke = 0.8
    ) +
    geom_text(
      data = nodes.dt,
      aes(
        x = x, y = y, label = name
      ),
      size = 2.0, fontface = "bold"
    ) +
    coord_cartesian(
      xlim = c(-1.4, 4.4),
      ylim = c(-0.38, 0.82)
    ) +
    labs(title = domain_label) +
    theme_void() +
    theme(
      plot.title = element_text(
        face = "bold", size = 9,
        hjust = 0.5
      )
    )

  p
}

# ----------------------------------------------------------
# Coupling Figure: Trajectories + Scatter
# ----------------------------------------------------------

#' Melt wide LGCM input to long spaghetti data
#'
#' @param inp Wide-format LGCM input
#' @param domain Domain code (MEM, LAN, EXF)
#' @param sex Sex label
#' @return data.table with PTID, Sex, EDT, HVR_z,
#'   COG_z columns
melt_lgcm_long.fn <- function(
    inp, domain, sex) {
  edt_c <- paste0("EDT_T", 0:5)
  hvr_c <- paste0("HVR_Z_T", 0:5)
  cog_c <- paste0(
    "PHC_", domain, "_T", 0:5
  )
  n <- nrow(inp)
  long.dt <- data.table(
    PTID = rep(inp$PTID, each = 6L),
    Sex = sex,
    visit = rep(0:5, times = n),
    EDT = as.numeric(
      unlist(inp[, ..edt_c])
    ),
    HVR_z = as.numeric(
      unlist(inp[, ..hvr_c])
    ),
    COG_z = as.numeric(
      unlist(inp[, ..cog_c])
    )
  )
  long.dt[!is.na(EDT)]
}

#' Plot coupling figure (2-panel stacked)
#'
#' Two-panel figure showing HVR-cognitive
#' coupling from LGCM parallel process models.
#'
#' Panel A: Trajectories faceted Measure (rows)
#'   x Sex (cols). HVR shown once + 3 cognitive
#'   domains. Spaghetti + bold mean curves.
#' Panel B: Latent slope scatter faceted by
#'   domain (cols), both sexes overlaid, with
#'   model covariance annotations.
#'
#' @param parallel_results From script 12
#' @param lgcm_m_input Male wide-format data
#' @param lgcm_f_input Female wide-format data
#' @param pct_subsample Fraction of subjects
#'   for spaghetti (default 0.20 = 20%)
#' @return patchwork object
plot_coupling_figure <- function(
    parallel_results,
    lgcm_m_input,
    lgcm_f_input,
    pct_subsample = 0.20) {

  # -- constants --------------------------------
  DOMAINS <- c("MEM", "LAN", "EXF")
  DOM_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )
  DOM_FCT <- as.character(DOM_LABELS)
  # 4-level row factor for trajectory panel
  ROW_LEVS <- c("HVR", DOM_FCT)
  # Sex: remap to plural for display
  SEX_FROM <- c("Female", "Male")
  SEX_DISP <- c("Females", "Males")
  sex_cols <- get_palette("sex")
  names(sex_cols)[1:2] <- SEX_DISP
  lr <- parallel_results$linear_results
  edt_mean <-
    parallel_results$config$edt_mean
  BASE_SZ <- 10
  STRIP_SZ <- 9

  set.seed(42)

  # Shared strip theme (Figs 5-7 style)
  strip_thm <- ggplot2::theme(
    strip.text = ggplot2::element_text(
      size = STRIP_SZ, face = "bold"
    )
  )

  # ============================================
  # DATA PREPARATION
  # ============================================

  # Melt MEM domain for HVR + shared EDT
  m_ref <- melt_lgcm_long.fn(
    lgcm_m_input, "MEM", "Male"
  )
  f_ref <- melt_lgcm_long.fn(
    lgcm_f_input, "MEM", "Female"
  )
  edt_rng <- range(
    c(m_ref$EDT, f_ref$EDT),
    na.rm = TRUE
  )
  edt_seq <- seq(
    edt_rng[1], edt_rng[2],
    length.out = 100
  )

  # Subsample IDs (shared across domains)
  m_ids <- unique(m_ref$PTID)
  f_ids <- unique(f_ref$PTID)
  n_m <- ceiling(
    length(m_ids) * pct_subsample
  )
  n_f <- ceiling(
    length(f_ids) * pct_subsample
  )
  sub_m.v <- sample(m_ids, n_m)
  sub_f.v <- sample(f_ids, n_f)

  # HVR spaghetti (from MEM model data)
  hvr_spag.dt <- rbindlist(list(
    m_ref[
      PTID %in% sub_m.v,
      .(PTID, Sex, EDT,
        value = HVR_z,
        Measure = "HVR")
    ],
    f_ref[
      PTID %in% sub_f.v,
      .(PTID, Sex, EDT,
        value = HVR_z,
        Measure = "HVR")
    ]
  ))

  # HVR mean (one curve from MEM model)
  # Clip quadratic curves past turning point
  # to avoid spurious upswing at tails
  clip_quad_curve.fn <- function(
    edt, vals, mu_s, q_coef, edt_mu
  ) {
    if (q_coef >= 0 || mu_s >= 0) return(vals)
    # Turning pt: d/dt = mu_s + 2*q*(t-mu) = 0
    edt_peak <- -mu_s / (2 * q_coef) + edt_mu
    ifelse(edt < edt_peak, NA_real_, vals)
  }
  hvr_mean.lst <- list()
  for (sx in SEX_FROM) {
    key <- paste0(sx, "_MEM")
    r <- lr[[key]]
    fs <- r$factor_scores
    ec <- edt_seq - edt_mean
    hq <- if (!is.null(
      r$quad_means
    )) r$quad_means$hvr_q else 0
    raw_vals <- mean(fs$hvr_i) +
      mean(fs$hvr_s) * ec + hq * ec^2
    hvr_mean.lst[[sx]] <- data.table(
      Sex = sx,
      Measure = "HVR",
      EDT = edt_seq,
      value = clip_quad_curve.fn(
        edt_seq, raw_vals,
        mean(fs$hvr_s), hq, edt_mean
      )
    )
  }

  # Per-domain: COG spaghetti, COG means,
  # scatter data, covariance labels
  cog_spag.lst <- list()
  cog_mean.lst <- list()
  scat.lst <- list()
  cov.lst <- list()

  for (di in seq_along(DOMAINS)) {
    d <- DOMAINS[di]
    d_lab <- DOM_LABELS[d]
    m_long <- melt_lgcm_long.fn(
      lgcm_m_input, d, "Male"
    )
    f_long <- melt_lgcm_long.fn(
      lgcm_f_input, d, "Female"
    )

    # COG spaghetti
    cog_spag.lst[[d]] <- rbindlist(list(
      m_long[
        PTID %in% sub_m.v,
        .(PTID, Sex, EDT,
          value = COG_z,
          Measure = d_lab)
      ],
      f_long[
        PTID %in% sub_f.v,
        .(PTID, Sex, EDT,
          value = COG_z,
          Measure = d_lab)
      ]
    ))

    for (sx in SEX_FROM) {
      key <- paste0(sx, "_", d)
      r <- lr[[key]]
      fs <- r$factor_scores
      ec <- edt_seq - edt_mean
      cq <- if (!is.null(
        r$quad_means
      )) r$quad_means$cog_q else 0

      raw_cog <- mean(fs$cog_i) +
        mean(fs$cog_s) * ec + cq * ec^2
      cog_mean.lst[[key]] <-
        data.table(
          Sex = sx,
          Measure = d_lab,
          EDT = edt_seq,
          value = clip_quad_curve.fn(
            edt_seq, raw_cog,
            mean(fs$cog_s), cq,
            edt_mean
          )
        )

      # Scatter data
      fs_cp <- copy(fs)
      fs_cp[, Sex := sx]
      fs_cp[, Domain := d_lab]
      scat.lst[[key]] <- fs_cp[
        , .(PTID, hvr_s, cog_s,
            Sex, Domain)
      ]

      # Model covariance
      cpl <- r$coupling
      sig <- ""
      if (cpl$p < 0.001) {
        sig <- "***"
      } else if (cpl$p < 0.01) {
        sig <- "**"
      } else if (cpl$p < 0.05) {
        sig <- "*"
      }
      sx_pre <- substr(sx, 1, 1)
      r_val <- cpl$cor
      if (is.null(r_val)) {
        pe <- r$params
        v_h <- pe$Estimate[
          pe$name == "var_hvr_s"
        ]
        v_c <- pe$Estimate[
          pe$name == "var_cog_s"
        ]
        r_val <- cpl$cov / sqrt(v_h * v_c)
      }
      cov.lst[[key]] <- data.table(
        Sex = sx,
        Domain = d_lab,
        label = sprintf(
          "%s: r = %.3f%s",
          sx_pre, r_val, sig
        )
      )
    }
  }

  # -- combine ---------------------------------
  all_spag.dt <- rbindlist(
    c(list(hvr_spag.dt),
      cog_spag.lst)
  )
  all_mean.dt <- rbindlist(
    c(hvr_mean.lst, cog_mean.lst)
  )
  all_scat.dt <- rbindlist(scat.lst)
  cov_lbl.dt <- rbindlist(cov.lst)

  # Set factors (remap Sex to plural)
  for (dt in list(
    all_spag.dt, all_mean.dt
  )) {
    dt[, Sex := factor(
      Sex, levels = SEX_FROM,
      labels = SEX_DISP
    )]
    dt[, Measure := factor(
      Measure, levels = ROW_LEVS
    )]
  }
  for (dt in list(
    all_scat.dt, cov_lbl.dt
  )) {
    dt[, Sex := factor(
      Sex, levels = SEX_FROM,
      labels = SEX_DISP
    )]
    dt[, Domain := factor(
      Domain, levels = DOM_FCT
    )]
  }

  # Negate slopes: positive = decline
  all_scat.dt[, hvr_s := -hvr_s]
  all_scat.dt[, cog_s := -cog_s]

  # Trim HVR outliers in spaghetti
  HVR_CLIP <- 4
  all_spag.dt[
    Measure == "HVR" &
      value < -HVR_CLIP,
    value := NA_real_
  ]

  # Covariance label positions
  # Stacked vertically, top-left corner
  # x = cog_s, y = hvr_s
  lbl_stats.dt <- all_scat.dt[, .(
    y_max = max(
      hvr_s, na.rm = TRUE
    ),
    x_min = min(
      cog_s, na.rm = TRUE
    ),
    y_gap = diff(range(
      hvr_s, na.rm = TRUE
    )) * 0.06
  ), by = Domain]
  cov_lbl.dt <- merge(
    cov_lbl.dt, lbl_stats.dt,
    by = "Domain"
  )
  cov_lbl.dt[, x_pos := x_min]
  cov_lbl.dt[
    Sex == "Females",
    y_pos := y_max
  ]
  cov_lbl.dt[
    Sex == "Males",
    y_pos := y_max - y_gap
  ]

  # ============================================
  # PANEL A: Trajectories (2x2, sexes merged)
  # ============================================
  p_traj <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = all_spag.dt,
      ggplot2::aes(
        x = EDT, y = value,
        group = PTID,
        colour = Sex
      ),
      alpha = 0.12,
      linewidth = 0.25
    ) +
    ggplot2::geom_line(
      data = all_mean.dt,
      ggplot2::aes(
        x = EDT, y = value,
        colour = Sex
      ),
      linewidth = 1.0
    ) +
    ggplot2::scale_colour_manual(
      values = sex_cols,
      guide = "none"
    ) +
    ggplot2::facet_wrap(
      ~ Measure, nrow = 2
    ) +
    ggplot2::labs(
      x = "EDT (years)",
      y = "z-score"
    ) +
    theme_publication(
      base_size = BASE_SZ
    ) +
    strip_thm

  # ============================================
  # PANEL B: Coupling scatter (domain cols)
  # ============================================
  p_scat <- ggplot2::ggplot(
    all_scat.dt,
    ggplot2::aes(
      x = cog_s, y = hvr_s,
      colour = Sex
    )
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      colour = "grey75",
      linetype = "solid",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      shape = 1, alpha = 0.4,
      size = 1.2, stroke = 0.3
    ) +
    ggplot2::geom_smooth(
      method = "lm", se = FALSE,
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::geom_label(
      data = cov_lbl.dt,
      ggplot2::aes(
        x = x_pos, y = y_pos,
        label = label,
        colour = Sex
      ),
      fill = scales::alpha(
        "white", 0.7
      ),
      label.size = 0,
      label.padding = ggplot2::unit(
        0.15, "lines"
      ),
      size = 2.5,
      hjust = 0, vjust = 1,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = sex_cols,
      guide = "none"
    ) +
    ggplot2::facet_wrap(
      ~ Domain, nrow = 1,
      scales = "free_x"
    ) +
    ggplot2::labs(
      x = "Latent cognitive decline",
      y = "Latent HVR decline"
    ) +
    theme_publication(
      base_size = BASE_SZ
    ) +
    strip_thm

  # ============================================
  # COMBINE: A / B with patchwork
  # ============================================
  (p_traj / p_scat) +
    patchwork::plot_layout(
      heights = c(3, 2)
    ) +
    patchwork::plot_annotation(
      tag_levels = "A"
    )
}

# -----------------------------------------------------------
# Mediation figure helpers
# -----------------------------------------------------------

#' Convert remapped mediation result to panel input
#'
#' @param r  Single result with a/b/cprime/indirect
#' @return List for plot_mediation_panel.fn
#' @export
make_med_list <- function(r) {
  list(
    converged = TRUE,
    a = list(
      est = r$a$est, p = r$a$p,
      boot_ci_lower = r$a$boot_ci_lower,
      boot_ci_upper = r$a$boot_ci_upper
    ),
    b = list(
      est = r$b$est, p = r$b$p,
      boot_ci_lower = r$b$boot_ci_lower,
      boot_ci_upper = r$b$boot_ci_upper
    ),
    cprime = list(
      est = r$cprime$est,
      p = r$cprime$p,
      boot_ci_lower = r$cprime$boot_ci_lower,
      boot_ci_upper = r$cprime$boot_ci_upper
    ),
    indirect = list(
      est = r$indirect$est,
      sobel_p = r$indirect$sobel_p,
      boot_ci_lower =
        r$indirect$boot_ci_lower,
      boot_ci_upper =
        r$indirect$boot_ci_upper,
      boot_significant =
        r$indirect$boot_significant
    ),
    prop_mediated = r$prop_mediated
  )
}

#' Build 3x2 mediation panel grid
#'
#' @param male_res  Per-sex adapter ($results)
#' @param female_res  Per-sex adapter
#' @param domain_labels  Named vec c(MEM=,LAN=,EXF=)
#' @return List of ggplot panels
#' @export
build_combined_med_panels <- function(
    male_res, female_res,
    domain_labels = c(
      MEM = "Memory",
      LAN = "Language",
      EXF = "Executive Function"
    )) {
  panels <- list()
  pred_cfg <- list(
    FRS = list(key = "FRS", lab = "FRS"),
    CVR = list(
      key = "CVR", lab = "CVR MIMIC"
    )
  )
  for (d in c("MEM", "LAN", "EXF")) {
    for (pc in pred_cfg) {
      k <- paste0(d, "_", pc$key)
      m_r <- male_res$results[[k]]
      f_r <- female_res$results[[k]]
      if (is.null(m_r) &&
          is.null(f_r)) next
      panels[[length(panels) + 1]] <-
        plot_mediation_panel.fn(
          make_med_list(m_r),
          make_med_list(f_r),
          domain_label = paste(
            domain_labels[d], "-",
            pc$lab
          ),
          predictor_label = pc$lab
        )
    }
  }
  panels[!sapply(panels, is.null)]
}
