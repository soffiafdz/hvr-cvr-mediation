# =============================================================================
# Table Utilities
# =============================================================================
# Common table formatting and generation functions using gt
# =============================================================================

# ----- Helper Functions -----

#' Wrap LaTeX table fragment in standalone document
#' @param tex_fragment_path Path to .tex file with table code
#' @param output_path Path for standalone .tex document (optional)
#' @param title Document title (optional, if NULL table title is used)
#' @return Path to standalone document
wrap_latex_table <- function(
    tex_fragment_path, output_path = NULL, title = NULL) {
  if (is.null(output_path)) {
    output_path <- sub("\\.tex$", "_standalone.tex", tex_fragment_path)
  }

  # Read the table fragment
  table_code <- readLines(tex_fragment_path)

  # Create standalone document - minimal header, table has its own title
  doc_header <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{geometry}",
    "\\geometry{a4paper, margin=0.75in}",
    "\\usepackage{caption}",
    "\\begin{document}",
    "\\pagestyle{empty}" # No page numbers
  )

  # Only add section title if explicitly provided
  # Otherwise the table's own title/subtitle will be used
  if (!is.null(title)) {
    doc_header <- c(doc_header, "", paste0("\\section*{", title, "}"), "")
  }

  doc_footer <- c("", "\\end{document}")

  # Combine
  standalone_doc <- c(doc_header, "", table_code, doc_footer)

  # Write
  writeLines(standalone_doc, output_path)

  return(output_path)
}

#' Compile multiple LaTeX table fragments into organized document
#' @param table_dir Directory containing .tex table fragments
#' @param output_path Path for compiled document
#' @param pattern File pattern to match (default "normative_*.tex")
#' @param title Document title
#' @return Path to compiled document
compile_normative_tables <- function(
    table_dir,
    output_path = NULL,
    pattern = "normative_.*\\.tex$",
    title = "Normative Centile Tables") {
  if (is.null(output_path)) {
    output_path <- file.path(table_dir, "all_normative_tables.tex")
  }

  # Find all matching table files
  table_files <- list.files(table_dir, pattern = pattern, full.names = TRUE)

  if (length(table_files) == 0) {
    stop("No table files found matching pattern: ", pattern, call. = FALSE)
  }

  # Parse filenames to organize by ROI → ADJ → SEX → SIDE
  # Format: normative_ROI_ADJ_SIDE_sex.tex
  file_info <- lapply(basename(table_files), \(fname) {
    parts <- strsplit(sub("\\.tex$", "", fname), "_")[[1]]
    list(
      file = fname,
      roi = parts[2],
      adj = parts[3],
      side = parts[4],
      sex = parts[5]
    )
  })

  # Sort by ROI, ADJ, SEX, SIDE
  file_info <- file_info[order(
    sapply(file_info, `[[`, "roi"),
    sapply(file_info, `[[`, "adj"),
    sapply(file_info, `[[`, "sex"),
    # sapply(file_info, `[[`, "side")
    factor(sapply(file_info, `[[`, "side"), levels = c("L", "R", "LR"))
  )]

  # Create document header with TOC
  doc_lines <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{geometry}",
    "\\geometry{a4paper, margin=1in}",
    "\\usepackage{caption}",
    "\\usepackage{titlesec}",
    "% Make section headers larger",
    paste0(
      "\\titleformat{\\section}{\\Large\\bfseries}",
      "{\\thesection}{1em}{}"
    ),
    paste0(
      "\\titleformat{\\subsection}{\\large\\bfseries}",
      "{\\thesubsection}{1em}{}"
    ),
    "\\begin{document}",
    "\\pagestyle{plain}",
    "",
    sprintf("\\begin{center}\\Large\\textbf{%s}\\end{center}", title),
    "\\vspace{1em}",
    "",
    "\\tableofcontents",
    "\\clearpage",
    ""
  )

  # Get labels from config
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()

  # Add tables organized by sections
  current_roi <- NULL
  current_adj <- NULL

  for (info in file_info) {
    # ROI section
    if (is.null(current_roi) || current_roi != info$roi) {
      roi_label <- if (info$roi %in% names(roi_labels)) {
        roi_labels[info$roi]
      } else {
        info$roi
      }
      doc_lines <- c(
        doc_lines,
        "",
        sprintf("\\section{%s}", roi_label),
        ""
      )
      current_roi <- info$roi
      current_adj <- NULL
    }

    # Adjustment subsection
    if (is.null(current_adj) || current_adj != info$adj) {
      adj_label <- if (info$adj %in% names(adj_labels)) {
        adj_labels[info$adj]
      } else {
        info$adj
      }
      doc_lines <- c(
        doc_lines,
        sprintf("\\subsection{%s}", adj_label),
        ""
      )
      current_adj <- info$adj
    }

    # Include table - sex and side are in table title
    table_path <- file.path(table_dir, info$file)
    doc_lines <- c(
      doc_lines,
      readLines(table_path),
      "\\clearpage",
      ""
    )
  }

  # Document footer
  doc_lines <- c(doc_lines, "\\end{document}")

  # Write
  writeLines(doc_lines, output_path)

  return(output_path)
}

#' Bin ages into intervals
#' @param ages Vector of ages
#' @param bin_width Width of age bins in years
#' @return Character vector of age bin labels
bin_ages <- function(ages, bin_width = 5) {
  breaks <- seq(
    floor(min(ages) / bin_width) * bin_width,
    ceiling(max(ages) / bin_width) * bin_width,
    by = bin_width
  )
  bins <- cut(ages, breaks = breaks, include.lowest = TRUE, right = FALSE)
  # Format as "45-50", "50-55", etc.
  levels(bins) <- sapply(seq_along(breaks[-length(breaks)]), \(i) {
    sprintf("%d-%d", breaks[i], breaks[i + 1])
  })
  return(bins)
}

# ----- Table Styling -----

#' Apply manuscript-style formatting to gt tables
#'
#' Consistent styling for publication-ready tables.
#' Adapted from headsize project.
#'
#' @param gt_tbl A gt table object
#' @param table_type Type: "main" (larger), "supplementary" (smaller), "compact"
#' @return Formatted gt table
#' @export
style_manuscript_table <- function(gt_tbl, table_type = "main") {
  # Font sizes based on table type (in pt for PDF compatibility)
  if (table_type == "main") {
    title_size <- 12
    subtitle_size <- 10
    body_size <- 9
    footnote_size <- 8
    source_size <- 8
  } else if (table_type == "supplementary") {
    title_size <- 11
    subtitle_size <- 9
    body_size <- 8
    footnote_size <- 7
    source_size <- 7
  } else if (table_type == "compact") {
    title_size <- 10
    subtitle_size <- 8
    body_size <- 7
    footnote_size <- 6
    source_size <- 6
  } else {
    stop("table_type must be 'main', 'supplementary', or 'compact'")
  }

  # Note: In LaTeX, cols_align() sets the tabular column spec which affects both

  # headers AND data. To get centered headers with right-aligned data, we must
  # use centered columns and then style data cells specifically.
  gt_tbl |>
    gt::tab_options(
      table.font.size = gt::px(body_size),
      quarto.disable_processing = FALSE,
      heading.title.font.size = gt::px(title_size),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(subtitle_size),
      column_labels.font.size = gt::px(body_size),
      column_labels.font.weight = "bold",
      footnotes.font.size = gt::px(footnote_size),
      footnotes.multiline = FALSE,
      source_notes.font.size = gt::px(source_size),
      row_group.font.size = gt::px(body_size),
      row_group.font.weight = "bold",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      heading.border.bottom.style = "solid"
    ) |>
    # Center all columns (affects both headers and data in LaTeX)
    gt::cols_align(align = "center", columns = gt::everything()) |>
    # Left-align first column (row labels)
    gt::cols_align(align = "left", columns = 1) |>
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = gt::cells_column_spanners()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    )
}

#' Apply significance styling to table
#' @param gt_table gt table object
#' @param sig_col Column indicating significance
#' @param style_cols Columns to style
#' @return Styled gt table
apply_significance_style <- function(
    gt_table, sig_col = "SIGN", style_cols = NULL) {
  if (is.null(style_cols)) {
    style_cols <- gt::everything()
  }

  gt_table |>
    gt::tab_style(
      style = gt::cell_text(style = "italic"),
      locations = gt::cells_body(
        columns = style_cols,
        rows = .data[[sig_col]] == FALSE
      )
    ) |>
    gt::cols_hide(columns = gt::contains(sig_col))
}

#' Format p-values consistently
#' @param pval P-value
#' @param threshold Threshold for "<" notation
#' @param stars If TRUE, return significance stars instead
#' @return Formatted p-value string
format_pval <- function(pval, threshold = 0.001, stars = FALSE) {
  if (stars) {
    # Return significance stars
    ifelse(pval < 0.001, "***",
      ifelse(pval < 0.01, "**",
        ifelse(pval < 0.05, "*", "")
      )
    )
  } else {
    # Return formatted p-value
    ifelse(pval < threshold,
      sprintf("<%s", threshold),
      sprintf("%.3f", pval)
    )
  }
}

#' Format p-value for display (vectorized, NA-safe)
#' Consolidated from reports-src/_common.R for single source of truth
#' @param p P-value (scalar or vector)
#' @return Formatted p-value string(s)
#' @export
format_p <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.0001) return("<0.0001")
    if (x < 0.001) return("<0.001")
    if (x < 0.01) return(sprintf("%.3f", x))
    return(sprintf("%.2f", x))
  })
}

#' Alias for format_p
#' @param p P-value
#' @return Formatted p-value string
#' @export
fmt_p <- function(p) format_p(p)

#' Get significance stars (vectorized, NA-safe)
#' Consolidated from reports-src/_common.R for single source of truth
#' @param p P-value (scalar or vector)
#' @return Significance stars string(s)
#' @export
sig_stars <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("")
    if (x < 0.001) return("***")
    if (x < 0.01) return("**")
    if (x < 0.05) return("*")
    return("")
  })
}

#' Format effect size with CI
#' @param estimate Point estimate
#' @param ci_lower Lower CI bound
#' @param ci_upper Upper CI bound
#' @param digits Number of decimal places
#' @return Formatted string
format_effect_ci <- function(estimate, ci_lower, ci_upper, digits = 2) {
  fmt <- sprintf("%%.%df [%%.%df, %%.%df]", digits, digits, digits)
  sprintf(fmt, estimate, ci_lower, ci_upper)
}

# ----- Abbreviations and Labels -----

#' Get family/distribution abbreviations and labels
#' @return Named list with abbrev and label vectors
get_family_info <- function() {
  list(
    abbrev = c(
      NO = "Norm",
      L_NO = "L-Norm",
      BCCG = "BCCG",
      BE = "Beta"
    ),
    label = c(
      NO = "Normal (Gaussian)",
      L_NO = "Logit-Normal",
      BCCG = "Box-Cox Cole & Green",
      BE = "Beta"
    )
  )
}

#' Get adjustment method abbreviations (lowercase)
#' @return Named vector of abbreviations
get_adjustment_abbrev <- function() {
  c(
    NON = "Unadj.",
    PRP = "Prop.",
    STX = "Stereo.",
    RES = "Resid."
  )
}

#' Get side labels and abbreviations
#' @return Named list with abbrev and label vectors
get_side_info <- function() {
  list(
    abbrev = c(L = "L", R = "R", LR = "Bilat."),
    label = c(L = "Left", R = "Right", LR = "Bilateral")
  )
}

#' Get interpretation symbols for validation quality
#' @return Named vector mapping interpretation to symbol
get_interpretation_symbols <- function() {
  c(
    Excellent = "***",
    Good = "**",
    Acceptable = "*",
    Poor = "\u2020" # dagger symbol
  )
}

#' Create abbreviation legend for source notes
#' @param type One of "adjustment", "family", "side", "interpretation"
#' @return Character string with abbreviation definitions
get_abbrev_legend <- function(
    type = c("adjustment", "family", "side", "interpretation"),
    filter_values = NULL) {
  type <- match.arg(type)

  switch(type,
    adjustment = paste(
      "Adjustment methods:",
      "Unadj. = Unadjusted (raw volumes);",
      "Prop. = Proportions (volume/TIV);",
      "Stereo. = Stereotaxic (atlas-normalized);",
      "Resid. = Residuals (TIV regressed out)"
    ),
    family = {
      all_fams <- c(
        "Norm" = "Norm = Normal (Gaussian)",
        "L-Norm" = "L-Norm = Logit-Normal",
        "BCCG" = "BCCG = Box-Cox Cole & Green",
        "Beta" = "Beta = Beta distribution"
      )
      if (!is.null(filter_values)) {
        all_fams <- all_fams[names(all_fams) %in% filter_values]
      }
      if (length(all_fams) == 0) {
        "Distributions: None"
      } else {
        paste(
          "Distributions:", paste(all_fams, collapse = "; ")
        )
      }
    },
    side = paste(
      "L = Left hemisphere; R = Right hemisphere;",
      "Bilat. = Bilateral (average)"
    ),
    interpretation = paste(
      "*** = Excellent (MAE < 5); ** = Good (5-10);",
      "* = Acceptable (10-15); \u2020 = Poor (> 15)"
    )
  )
}

# ----- LaTeX Post-processing -----

#' Fix LaTeX source notes formatting
#' @param tex_path Path to .tex file
#' @param table_type Type of table: "normative" or "summary"
#' @return Invisibly returns the path
fix_latex_source_notes <- function(
    tex_path, table_type = "normative", bilateral = FALSE) {
  if (!file.exists(tex_path)) {
    return(invisible(tex_path))
  }

  lines <- readLines(tex_path, warn = FALSE)

  # Title/subtitle sizes - consistent across all tables
  title_size <- "16"
  title_leading <- "20"
  subtitle_size <- "14"
  subtitle_leading <- "18"

  # Define table body and source notes font sizes based on table type
  if (table_type == "normative" || table_type == "diagnostic") {
    # Normative and validation diagnostic tables: same formatting
    table_size <- "11"
    table_leading <- "14"
    source_notes_size <- "9"
  } else if (table_type == "diagnostic_stability") {
    # Stability diagnostic tables: slightly smaller
    table_size <- "9"
    table_leading <- "11.5"
    source_notes_size <- "8"
  } else if (table_type == "diagnostic_summary") {
    # Model summary diagnostic tables: smallest (wide tables)
    table_size <- "7.5"
    table_leading <- "9.5"
    source_notes_size <- "7"
  } else if (table_type == "summary") {
    # Summary tables: smaller fonts to fit wide tables
    if (bilateral) {
      table_size <- "9.5"
      table_leading <- "10"
      source_notes_size <- "6"
    } else {
      table_size <- "6"
      table_leading <- "7.5"
      source_notes_size <- "6"
    }
  } else {
    stop(
      "table_type must be 'normative', 'summary', ",
      "'diagnostic', 'diagnostic_stability', ",
      "or 'diagnostic_summary'", call. = FALSE
    )
  }

  # Fix table font size
  fontsize_line <- grep("^\\\\fontsize\\{", lines)[1]
  if (!is.na(fontsize_line)) {
    lines[fontsize_line] <- sprintf(
      "\\fontsize{%spt}{%spt}\\selectfont",
      table_size, table_leading
    )
  }

  # Fix title/subtitle font sizes in caption
  for (i in seq_along(lines)) {
    if (grepl("\\\\caption\\*\\{", lines[i])) {
      # Next few lines should contain title and subtitle
      # Title line
      if (
        i + 1 <= length(lines) &&
          grepl("\\{\\\\fontsize\\{", lines[i + 1])
      ) {
        lines[i + 1] <- gsub(
          "\\{\\\\fontsize\\{[0-9.]+\\}\\{[0-9.]+\\}",
          sprintf("{\\\\fontsize{%s}{%s}", title_size, title_leading),
          lines[i + 1]
        )
      }
      # Subtitle line
      if (
        i + 2 <= length(lines) &&
          grepl("\\{\\\\fontsize\\{", lines[i + 2])
      ) {
        lines[i + 2] <- gsub(
          "\\{\\\\fontsize\\{[0-9.]+\\}\\{[0-9.]+\\}",
          sprintf(
            "{\\\\fontsize{%s}{%s}",
            subtitle_size, subtitle_leading
          ),
          lines[i + 2]
        )
      }
      # Add vertical space between title and subtitle
      # Replace \\ at end of title line with \\[6pt]
      if (i + 1 <= length(lines)) {
        lines[i + 1] <- gsub(
          "\\\\\\\\\\s*$",
          "\\\\\\\\[6pt]",
          lines[i + 1]
        )
      }
      break
    }
  }

  # Add longtable centering if not present
  ltpost_line <- grep("\\\\setlength\\{\\\\LTpost\\}", lines)
  if (length(ltpost_line) > 0) {
    # Check if LTleft/LTright already exist
    has_ltleft <- any(grepl("\\\\setlength\\{\\\\LTleft\\}", lines))
    if (!has_ltleft) {
      # Insert centering settings after LTpost
      lines <- c(
        lines[1:ltpost_line[1]],
        "\\setlength{\\LTleft}{0pt plus 1fill}",
        "\\setlength{\\LTright}{0pt plus 1fill}",
        lines[(ltpost_line[1] + 1):length(lines)]
      )
    }
  }

  # Find minipage or center with source notes
  minipage_start <- grep("\\\\begin\\{minipage\\}", lines)
  center_start <- grep("\\\\begin\\{center\\}", lines)

  # Process minipage if found
  if (length(minipage_start) > 0) {
    for (start_idx in minipage_start) {
      # Find end of minipage
      end_idx <- start_idx
      while (
        end_idx <= length(lines) &&
          !grepl("\\\\end\\{minipage\\}", lines[end_idx])
      ) {
        end_idx <- end_idx + 1
      }

      if (end_idx <= length(lines)) {
        # Extract source note lines (between minipage tags)
        note_lines <- lines[(start_idx + 1):(end_idx - 1)]

        # Remove gt's formatting commands
        note_lines <- note_lines[
          !grepl("^\\\\centering", note_lines) &
            !grepl("^\\\\fontsize", note_lines)
        ]

        # Clean up trailing \\ but keep separate lines if multiple
        note_lines <- gsub("\\\\\\\\$", "", note_lines)

        # Remove empty or whitespace-only lines
        note_lines <- note_lines[trimws(note_lines) != ""]

        # Only consolidate if there are many short fragments (old behavior)
        # Keep separate if they're already organized as intended
        if (length(note_lines) > 3) {
          # Many fragments - consolidate
          combined_note <- paste(note_lines, collapse = " ")
          note_lines <- combined_note
        } else if (length(note_lines) > 1) {
          # 2-3 lines: add LaTeX line breaks between them
          note_lines <- paste(note_lines, collapse = " \\\\\n")
        }

        # Replace minipage - remove width spec and use centering
        # Add negative vspace to reduce gap between table and notes
        lines <- c(
          lines[1:(start_idx - 1)],
          "\\vspace{-12pt}",
          "\\begin{center}",
          sprintf(
            "\\fontsize{%s}{%s}\\selectfont", source_notes_size,
            as.character(as.numeric(source_notes_size) * 1.2)
          ),
          note_lines,
          "\\end{center}",
          lines[(end_idx + 1):length(lines)]
        )
      }
    }
  } else if (length(center_start) > 0) {
    # If center already exists, add vspace before it if missing
    for (start_idx in center_start) {
      # Check if vspace already exists before this center
      has_vspace <- start_idx > 1 &&
        grepl("\\\\vspace", lines[start_idx - 1])

      if (!has_vspace) {
        # Insert vspace before center
        lines <- c(
          lines[1:(start_idx - 1)],
          "\\vspace{-12pt}",
          lines[start_idx:length(lines)]
        )
      }
    }
  }

  # Fix escaped dollar signs (gt converts $ to \$)
  lines <- gsub("\\\\\\$", "$", lines)

  # Fix \times that gt escaped as \textbackslash{}times
  lines <- gsub(
    "\\\\textbackslash\\{\\}times", "\\\\times", lines
  )

  # Fix superscripts: \textasciicircum{}\{-3\} -> ^{-3}
  lines <- gsub(
    "\\\\textasciicircum\\{\\}\\\\\\{(-?[0-9]+)\\\\\\}",
    "^{\\1}",
    lines
  )

  # Fix percent signs - multiple patterns
  lines <- gsub("\\\\textbackslash\\{\\}\\\\%", "\\\\%", lines)
  lines <- gsub("\\\\textbackslash\\{\\}%", "\\\\%", lines)
  lines <- gsub("textbackslash\\{\\}%", "\\\\%", lines)

  # Fix Unicode symbols
  lines <- gsub("\u2020", "\\\\dag", lines)
  lines <- gsub("\u00d7", "$\\\\times$", lines)
  lines <- gsub("\u207b\u00b3", "$^{-3}$", lines)

  writeLines(lines, tex_path)

  invisible(tex_path)
}

# ----- Table Saving -----

# Note: Use ensure_directory() from data_io.R instead
# Kept for backward compatibility
#' @export
ensure_table_dir <- function(table_dir) {
  ensure_directory(table_dir)
}

#' Save table in multiple formats
#' @param gt_table gt table object
#' @param filename Base filename (without extension)
#' @param output_dir Output directory
#' @param formats Vector of formats ("html", "tex", "rtf")
#' @param fix_latex Whether to apply LaTeX post-processing fixes
#'   (default TRUE)
#' @param create_standalone Whether to create standalone tex documents
#'   (default TRUE)
save_table <- function(
    gt_table, filename, output_dir, formats = c("html", "tex"),
    fix_latex = TRUE, create_standalone = TRUE) {
  ensure_table_dir(output_dir)

  for (fmt in formats) {
    fpath <- file.path(output_dir, paste0(filename, ".", fmt))
    gt::gtsave(gt_table, fpath)

    # Apply LaTeX fixes and create standalone document
    if (fmt == "tex") {
      if (fix_latex) {
        fix_latex_source_notes(fpath)
      }
      if (create_standalone) {
        wrap_latex_table(fpath)
      }
    }
  }

  invisible(filename)
}

#' Save table to both HTML and LaTeX
#' @param gt_table gt table object
#' @param filename Base filename
#' @param html_dir HTML output directory
#' @param tex_dir LaTeX output directory
save_table_dual <- function(gt_table, filename, html_dir, tex_dir) {
  ensure_table_dir(html_dir)
  ensure_table_dir(tex_dir)

  gt::gtsave(gt_table, file.path(html_dir, paste0(filename, ".html")))
  gt::gtsave(gt_table, file.path(tex_dir, paste0(filename, ".tex")))

  invisible(filename)
}

# ----- Demographic Table Helpers -----

#' Create variable labels for demographics
#' @param use_markdown Whether to use LaTeX/markdown formatting
#' @return Named vector of labels
get_demographic_labels <- function(use_markdown = TRUE) {
  if (use_markdown) {
    c(
      AGE          = "$\\text{Age (years)}$",
      TIME_diff    = "$\\text{Time (months)}$",
      EDUC         = "$\\text{Education (years)}$",
      EDUC_NA      = "$\\text{Missing}$",
      ICC          = "$\\text{Head-size (TIV)}$",
      HC           = "$\\text{Hippocampus (cc)}$",
      VC           = "$\\text{Lat. ventricles (cc)}$",
      MEM_EXEC_z   = "$\\text{Memory \\& exec. func.}$",
      MEM_EXEC_NA  = "$\\text{Missing}$",
      PRSP_z       = "$\\text{Processing speed}$",
      PRSP_NA      = "$\\text{Missing}$"
    )
  } else {
    c(
      AGE          = "Age (years)",
      TIME_diff    = "Time (months)",
      EDUC         = "Education (years)",
      EDUC_NA      = "Missing",
      ICC          = "Head-size (TIV)",
      HC           = "Hippocampus (cc)",
      VC           = "Lat. ventricles (cc)",
      MEM_EXEC_z   = "Memory & exec. func.",
      MEM_EXEC_NA  = "Missing",
      PRSP_z       = "Processing speed",
      PRSP_NA      = "Missing"
    )
  }
}

#' Format mean and SD
#' @param values Vector of values
#' @param is_icc Whether this is ICC (use comma separator)
#' @param use_cog_format Whether to use cognitive format (3 decimals)
#' @return Formatted string
format_mean_sd <- function(values, is_icc = FALSE, use_cog_format = FALSE) {
  if (use_cog_format) {
    sprintf("%.3f (%.2f)", mean(values, na.rm = TRUE), sd(values, na.rm = TRUE))
  } else if (is_icc) {
    sprintf(
      "%s (%s)",
      format(mean(values, na.rm = TRUE), big.mark = ",", digits = 4),
      format(sd(values, na.rm = TRUE), big.mark = ",", digits = 4)
    )
  } else {
    sprintf("%.1f (%.1f)", mean(values, na.rm = TRUE), sd(values, na.rm = TRUE))
  }
}

#' Calculate t-test and format results
#' @param data Data table
#' @param variable Variable name
#' @param group_var Grouping variable
#' @return List with test statistics
ttest_summary <- function(data, variable, group_var = "SEX") {
  tt <- data[, t.test(get(variable) ~ get(group_var), na.rm = TRUE)]
  list(
    X = variable,
    Tstat = tt$statistic,
    DF = tt$parameter,
    Pval = tt$p.value,
    Pval_fmt = format_pval(tt$p.value)
  )
}

# ----- SEM Table Helpers -----

#' Create path labels for SEM models
#' @param use_markdown Whether to use LaTeX formatting
#' @return Named vector of path labels
get_sem_path_labels <- function(use_markdown = TRUE) {
  if (use_markdown) {
    list(
      tot_SES_COG = "$\\text{SES} \\to \\text{Cog (Total)}$",
      dir_SES_COG = "$\\text{SES} \\to \\text{Cog (Direct)}$",
      SES_ICC_HC_COG = paste0(
        "$\\text{SES} \\to \\text{TIV} \\to ",
        "\\text{HC} \\to \\text{Cog}$"
      ),
      SES_ICC_VC_COG = paste0(
        "$\\text{SES} \\to \\text{TIV} \\to ",
        "\\text{VC} \\to \\text{Cog}$"
      ),
      SES_HC_COG = "$\\text{SES} \\to \\text{HC} \\to \\text{Cog}$",
      SES_VC_COG = "$\\text{SES} \\to \\text{VC} \\to \\text{Cog}$",
      SES_ICC_COG = "$\\text{SES} \\to \\text{TIV} \\to \\text{Cog}$",
      tot_SEX_COG = "$\\text{Sex} \\to \\text{Cog (Total)}$",
      dir_SEX_COG = "$\\text{Sex} \\to \\text{Cog (Direct)}$",
      SEX_ICC_HC_COG = paste0(
        "$\\text{Sex} \\to \\text{TIV} \\to ",
        "\\text{HC} \\to \\text{Cog}$"
      ),
      SEX_ICC_VC_COG = paste0(
        "$\\text{Sex} \\to \\text{TIV} \\to ",
        "\\text{VC} \\to \\text{Cog}$"
      ),
      SEX_HC_COG = "$\\text{Sex} \\to \\text{HC} \\to \\text{Cog}$",
      SEX_VC_COG = "$\\text{Sex} \\to \\text{VC} \\to \\text{Cog}$",
      SEX_ICC_COG = "$\\text{Sex} \\to \\text{TIV} \\to \\text{Cog}$"
    )
  } else {
    list(
      tot_SES_COG     = "SES → Cog (Total)",
      dir_SES_COG     = "SES → Cog (Direct)",
      SES_ICC_HC_COG  = "SES → TIV → HC → Cog",
      SES_ICC_VC_COG  = "SES → TIV → VC → Cog",
      SES_HC_COG      = "SES → HC → Cog",
      SES_VC_COG      = "SES → VC → Cog",
      SES_ICC_COG     = "SES → TIV → Cog",
      tot_SEX_COG     = "Sex → Cog (Total)",
      dir_SEX_COG     = "Sex → Cog (Direct)",
      SEX_ICC_HC_COG  = "Sex → TIV → HC → Cog",
      SEX_ICC_VC_COG  = "Sex → TIV → VC → Cog",
      SEX_HC_COG      = "Sex → HC → Cog",
      SEX_VC_COG      = "Sex → VC → Cog",
      SEX_ICC_COG     = "Sex → TIV → Cog"
    )
  }
}

#' Format standardized coefficient
#' @param value Coefficient value
#' @param threshold Threshold for scientific notation
#' @return Formatted string
format_beta_std <- function(value, threshold = 0.01) {
  ifelse(abs(value) < threshold,
    sprintf("%.0e", value),
    sprintf("%.2f", value)
  )
}

#' Format confidence interval
#' @param lower Lower bound
#' @param upper Upper bound
#' @param threshold Threshold for scientific notation
#' @return Formatted string
format_ci <- function(lower, upper, threshold = 0.01) {
  lower_fmt <- ifelse(abs(lower) < threshold, "%.0e", "%.2f")
  upper_fmt <- ifelse(abs(upper) < threshold, "%.0e", "%.2f")
  sprintf(paste(lower_fmt, upper_fmt, sep = ", "), lower, upper)
}

# ----- Normative Table Helpers -----

#' Get adjustment method labels from config
#' @return Named vector
get_adjustment_labels <- function() {
  tryCatch(
    {
      params <- get_parameter("adjustment_methods")
      setNames(unlist(params), names(params))
    },
    error = function(e) {
      # Fallback if config not available (standalone mode)
      c(
        NON = "Unadjusted",
        PRP = "Proportions",
        STX = "Stereotaxic",
        RES = "Residuals"
      )
    }
  )
}

#' Get ROI labels from config
#' @return Named vector
get_roi_labels <- function() {
  tryCatch(
    {
      params <- get_parameter("rois")
      setNames(unlist(params), names(params))
    },
    error = function(e) {
      # Fallback if config not available (standalone mode)
      c(
        HC  = "Hippocampus",
        LV  = "Lateral Ventricles",
        HVR = "Hippocampus-to-Ventricle ratio",
        ICC = "Intracranial volume"
      )
    }
  )
}

#' Create publication-ready normative table
#' @param norm_data Data table with normative centiles
#' @param title Table title
#' @param subtitle Table subtitle
#' @param roi_type Type of ROI ("HC", "LV", "HVR") for clinical
#'   thresholds.
#' @param adj_type Adjustment type ("NON", "PRP", "STX", "RES")
#' @param family Whether there is a column with Family data
#' @param sample_size Whether there is a column with the subsample size
#' @param add_clinical_labels Whether to add clinical threshold
#'   spanners
#' @param format Output format ("html" or "latex")
#' @return gt table object
create_normative_table <- function(
    norm_data,
    title = "Normative Centiles",
    subtitle = NULL,
    roi_type = "HC",
    adj_type = "NON",
    family = FALSE,
    sample_size = FALSE,
    add_clinical_labels = FALSE,
    format = "html") {
  # Work with a copy to avoid modifying original data

  norm_data <- copy(norm_data)

  # Get centile columns
  cent_cols <- grep("^p[0-9.]+$", names(norm_data), value = TRUE)

  # For Proportions (PRP) on HC/LV, values are very small - multiply by 1000

  is_scaled <- FALSE
  if (adj_type == "PRP" && roi_type %in% c("HC", "LV")) {
    norm_data[, (cent_cols) := lapply(.SD, `*`, 1000), .SDcols = cent_cols]
    is_scaled <- TRUE
  }

  # Create labels - Age without (y), will add note in source notes
  cent_labels <- as.list(
    setNames(
      paste0(as.numeric(sub("p", "", cent_cols)), "th"),
      cent_cols
    )
  )
  cent_labels$AGE <- "Age"

  # Create base table
  gt_table <- norm_data[, c("AGE", cent_cols), with = FALSE] |>
    gt() |>
    tab_header(
      title = title,
      subtitle = subtitle
    ) |>
    cols_label(.list = cent_labels) |>
    fmt_number(columns = all_of(cent_cols), decimals = 2)

  gt_table <- gt_table |>
    cols_align(align = "left", columns = "AGE") |>
    tab_style(
      style = cell_text(align = "left"),
      locations = cells_column_labels(columns = "AGE")
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )

  # Source notes - each on separate line
  source_notes.lst <- list()
  source_notes.lst[[1]] <- character(0)

  if (sample_size && "N" %in% names(norm_data)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("N = %s;", norm_data[, format(sum(N), big.mark = ",")])
    )
  }

  if (family && "FAMILY" %in% names(norm_data)) {
    fam_info <- get_family_info()
    fam_code <- unique(norm_data$FAMILY)
    fam_label <- ifelse(
      fam_code %in% names(fam_info$label),
      fam_info$label[fam_code],
      fam_code
    )
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("%s distribution.", fam_label)
    )
  }

  # Add Age note
  source_notes.lst[[1]] <- c(
    source_notes.lst[[1]], "Age in years."
  )

  # Add scaling note for proportions
  if (is_scaled) {
    if (format == "latex") {
      source_notes.lst[[2]] <- paste0(
        "Values $\\times 10^{-3}$ ",
        "(multiply by 0.001 to obtain proportions)."
      )
    } else {
      source_notes.lst[[2]] <- paste0(
        "Values \u00d710\u207b\u00b3 ",
        "(multiply by 0.001 to obtain proportions)."
      )
    }
  }

  # Add clinical threshold information to source notes
  if (add_clinical_labels) {
    list_n <- length(source_notes.lst) + 1
    if (roi_type %in% c("HC", "HVR")) {
      if (format == "latex") {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: $<$2.5\\% = Severe atrophy (-2SD), ",
          "2.5-5\\% = Borderline (-1.6SD)"
        )
      } else {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: <2.5% = Severe atrophy (-2SD), ",
          "2.5-5% = Borderline (-1.6SD)"
        )
      }
    } else if (roi_type == "LV") {
      if (format == "latex") {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: $>$97.5\\% = ",
          "Abnormal enlargement (+2SD), ",
          "95-97.5\\% = Borderline (+1.6SD)"
        )
      } else {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: >97.5% = ",
          "Abnormal enlargement (+2SD), ",
          "95-97.5% = Borderline (+1.6SD)"
        )
      }
    }
  }

  # Format-specific styling
  if (format == "latex") {
    # Formal LaTeX: no colors, clean borders, centered table
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML: with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_fill(color = "#E8F4F8"),
        locations = cells_body(columns = "AGE")
      ) |>
      tab_options(
        table.font.size = px(12),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(16),
        heading.subtitle.font.size = px(14),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Add clinical Z-score spanners to gt table
#' @param gt_table gt table object
#' @param col_names Column names (excluding SEX and AGE_BIN)
#' @param roi_type Type of ROI ("HC", "LV", or "HVR") to determine thresholds
#' @return Modified gt table with Z-score spanners and clinical footnotes
#' @details Hierarchy:
#'   Level 1 (percentiles) <- Level 2 (clinical) <- Level 3 (side)
add_clinical_spanners <- function(gt_table, col_names, roi_type = "HC") {
  # --- 1. Define Z-Score Map (Universal) ---
  # Maps percentile numeric value -> Z-score label
  # Note: 50th percentile is "0 SD" (Average)
  z_map <- c(
    "1"    = "-2.3 SD",
    "2.5"  = "-2.0 SD",
    "5"    = "-1.6 SD",
    "10"   = "-1.3 SD",
    "25"   = "-0.7 SD",
    "50"   = "0 SD", # Median
    "75"   = "+0.7 SD",
    "90"   = "+1.3 SD",
    "95"   = "+1.6 SD",
    "97.5" = "+2.0 SD",
    "99"   = "+2.3 SD"
  )

  # Detect whether columns have:
  # - "SIDE_" prefix (Summary Table)
  # - just "pXX" (Main Table)
  has_side_prefix <- any(grepl("_[p][0-9]", col_names))

  # Extract numeric values
  if (has_side_prefix) {
    # Format: L_p1, R_p99
    p_vals <- unique(sub(".*_p", "", col_names))
    sides <- unique(sub("_p.*", "", col_names))
  } else {
    # Format: p1, p99 (Main Table)
    p_vals <- unique(sub("^p", "", col_names))
    sides <- "dummy" # No side prefix to span
  }

  # Filter map to only what exists in data
  z_map <- z_map[names(z_map) %in% p_vals]

  # Apply Spanners
  for (pz in names(z_map)) {
    lbl <- z_map[[pz]]

    if (has_side_prefix) {
      for (side in sides) {
        target <- paste0(side, "_p", pz)
        if (target %in% col_names) {
          gt_table <- gt_table |>
            tab_spanner(
              label = lbl,
              columns = all_of(target),
              level = 2,
              id = paste0(side, pz)
            )
        }
      }
    } else {
      # Main table (no side prefix)
      target <- paste0("p", pz)
      if (target %in% col_names) {
        gt_table <- gt_table |>
          tab_spanner(
            label = lbl,
            columns = all_of(target),
            level = 2,
            id = paste0("z", pz)
          )
      }
    }
  }

  return(gt_table)
}

#' Create summary normative table with hierarchical year-bins structure
#' @param norm_data_list Named list of normative tables in "Sex_Side" format
#'   (e.g., "Male_L", "Female_R")
#' @param title Table title
#' @param subtitle Table subtitle (optional, NULL for cleaner tables)
#' @param format Output format ("html" or "latex")
#' @param bin_width Age bin width in years (default 5)
#' @param centiles Vector of centiles to display
#'   (default c(.05, .25, .5, .75, .95))
#' @param roi_type Type of ROI for clinical thresholds ("HC", "LV", "HVR")
#' @param adj_type Adjustment type ("NON", "PRP", "STX", "RES")
#' @param add_clinical_labels Whether to add clinical threshold spanners
#' @return gt table object with hierarchical structure:
#'   - Row groups by Sex (Female/Male) with N counts
#'   - Column spanners by Side (Left/Right hemisphere)
#'   - Optional clinical threshold spanners
#'   - Columns for centiles within each side
#' @details Creates compact tables suitable for portrait pages by using
#'   vertical grouping for sex and horizontal spanners for sides
create_summary_normative_table <- function(
    norm_data_list,
    title = "Normative Centiles Summary",
    subtitle = NULL,
    format = "html",
    bin_width = 5,
    centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
    roi_type = "HC",
    adj_type = "NON",
    family = FALSE,
    sample_size = FALSE,
    add_clinical_labels = FALSE) {
  # Process and combine data with proper hierarchical structure
  # Structure: Sex (row groups) → Side (column spanners) → Centiles (columns)
  # For bilateral: Sex (row groups) → Centiles (columns only)

  # Determine if we need to scale for proportions
  is_scaled <- adj_type == "PRP" && roi_type %in% c("HC", "LV")

  # Detect if this is bilateral data (all sides are "LR")
  all_sides <- unique(sapply(strsplit(names(norm_data_list), "_"), `[`, 2))
  is_bilateral <- length(all_sides) == 1 && all_sides[1] == "LR"

  # Extract family information if available (check first table)
  family_label <- NULL
  if (family && length(norm_data_list) > 0) {
    first_dt <- norm_data_list[[1]]
    if ("FAMILY" %in% names(first_dt)) {
      fam_info <- get_family_info()
      fam_code <- unique(first_dt$FAMILY)
      family_label <- ifelse(
        fam_code %in% names(fam_info$label),
        fam_info$label[fam_code],
        fam_code
      )
    }
  }

  processed_list <- lapply(names(norm_data_list), function(key) {
    dt <- copy(norm_data_list[[key]])

    # Parse key to extract sex and side (format: "Sex_Side")
    parts <- strsplit(key, "_")[[1]]
    sex <- parts[1]
    side <- parts[2]

    # Add age bins
    dt[, AGE_BIN := bin_ages(AGE, bin_width)]

    # Get centile column names
    cent_cols <- paste0("p", centiles * 100)

    # Scale proportions if needed (multiply by 1000)
    if (is_scaled) {
      dt[, (cent_cols) := lapply(.SD, `*`, 1000), .SDcols = cent_cols]
    }

    # Average within bins
    binned <- merge(
      dt[, .(N = sum(N)), "AGE_BIN"],
      dt[, lapply(.SD, mean), AGE_BIN, .SDcols = cent_cols]
    )

    # Add grouping variables for hierarchical structure
    binned[, `:=`(SEX = sex, SIDE = side)]

    binned
  })

  # Combine all data
  combined_dt <- rbindlist(processed_list)
  id_cols <- c("SEX", "AGE_BIN")
  sex_labels <- total_n <- NULL
  # Create sex labels with N counts (pluralized)
  if (sample_size && "N" %in% names(combined_dt)) {
    id_cols <- c(id_cols, "N")
    sex_n <- combined_dt[, .(N = sum(N)), .(SEX, SIDE)]
    if (any(sex_n[, any(uniqueN(N) != 1), SEX]$V1)) {
      stop(
        "Left/Right sample sizes are unequal ",
        "within at least one SEX",
        call. = FALSE
      )
    }
    total_n <- sex_n[, unique(.SD[, -"SIDE"])][, sum(N)]
    sex_labels <- sex_n[, unique(.SD[, -"SIDE"])][, setNames(
      sprintf("%ss (N = %s)", SEX, format(N, big.mark = ",")),
      SEX
    )]
  }

  # Reshape data based on whether it's bilateral or unilateral
  cent_cols <- paste0("p", centiles * 100)

  if (is_bilateral) {
    # Bilateral: Simple structure - Sex | Age | p5 | p25 | p50 | ...
    # No SIDE dimension needed
    display_dt <- copy(combined_dt)
    display_dt[, SIDE := NULL]
    setorder(display_dt, SEX, AGE_BIN)
    all_cols <- cent_cols
  } else {
    # Unilateral: Age | 5th_L | 5th_R | 25th_L | 25th_R | ...
    # Centile spanners on top, L/R spanners below

    # Melt to long
    long_dt <- melt(
      combined_dt,
      measure.vars = cent_cols,
      variable.name = "CENTILE",
      value.name = "VALUE"
    )

    # Create centile-side identifier for columns
    long_dt[, CENT_SIDE := paste(CENTILE, SIDE, sep = "_")]

    # Cast to wide: one row per Sex + Age_Bin (+ N) combination
    display_dt <- dcast(
      long_dt[, -c("SIDE", "CENTILE")],
      ... ~ CENT_SIDE,
      value.var = "VALUE"
    )

    # Sort by sex and age
    setorder(display_dt, SEX, AGE_BIN)

    # Reorder columns to ensure proper grouping
    # Expected order: AGE_BIN, p5_L, p5_R, p25_L, p25_R, ...
    data_cols <- setdiff(names(display_dt), id_cols)

    # Separate by centile and ensure L/R order within each centile
    centiles_in_data <- unique(
      sapply(data_cols, \(x) strsplit(x, "_")[[1]][1])
    )
    # Sort centiles numerically
    cent_values <- as.numeric(sub("p", "", centiles_in_data))
    centiles_in_data <- centiles_in_data[order(cent_values)]

    ordered_cols <- unlist(lapply(centiles_in_data, \(c) {
      cent_cols_side <- grep(paste0("^", c, "_"), data_cols, value = TRUE)
      # Ensure L comes before R
      cent_cols_side[order(sub(".*_", "", cent_cols_side))]
    }))

    # Reorder display_dt columns
    setcolorder(display_dt, c(id_cols, ordered_cols))
    all_cols <- ordered_cols
  }

  # Replace SEX values with labeled versions (with N counts)
  if (!is.null(sex_labels)) {
    display_dt[, SEX := sex_labels[SEX]]
  } else {
    display_dt[, SEX := paste0(SEX, "s")]
  }

  # Create the gt table - use AGE_BIN as row name, hide N if present
  # Remove N column and AGE_BIN column from data columns for display
  if (sample_size && "N" %in% names(display_dt)) {
    # Create combined Age(N) label for row names
    display_dt[, AGE_LABEL := sprintf("%s (%d)", AGE_BIN, N)]
    display_dt[, c("AGE_BIN", "N") := NULL]

    gt_table <- display_dt |>
      gt(groupname_col = "SEX", rowname_col = "AGE_LABEL") |>
      tab_stubhead(label = "Age (N)") |>
      fmt_number(columns = all_of(all_cols), decimals = 2)
  } else {
    gt_table <- display_dt |>
      gt(groupname_col = "SEX", rowname_col = "AGE_BIN") |>
      tab_stubhead(label = "Age") |>
      fmt_number(columns = all_of(all_cols), decimals = 2)
  }

  # Add title if provided
  if (!is.null(title)) {
    gt_table <- gt_table |>
      tab_header(
        title = title,
        subtitle = subtitle
      )
  }

  # Create column labels and spanners based on data type
  if (is_bilateral) {
    # Bilateral: Simple centile labels (1th, 2.5th, 5th, etc.)
    col_labels <- setNames(
      paste0(as.numeric(sub("p", "", all_cols)), "th"),
      all_cols
    )
    gt_table <- gt_table |> cols_label(.list = col_labels)
    # No spanners needed for bilateral
  } else {
    # Unilateral: Show side (L or R) as column headers
    col_labels <- setNames(
      sub("^.*_", "", all_cols),
      all_cols
    )
    gt_table <- gt_table |> cols_label(.list = col_labels)

    # Add centile spanners (groups L/R columns by centile)
    # Extract unique centiles from column names
    centiles_in_data <- unique(
      sapply(all_cols, \(x) strsplit(x, "_")[[1]][1])
    )
    for (centile in centiles_in_data) {
      cent_label <- paste0(as.numeric(sub("p", "", centile)), "th")
      cent_pattern_cols <- grep(
        paste0("^", centile, "_"), all_cols,
        value = TRUE
      )

      gt_table <- gt_table |>
        tab_spanner(
          label = cent_label,
          id = paste0("centile_", centile),
          columns = all_of(cent_pattern_cols),
          level = 1
        )
    }
  }

  # Source notes - each on separate line
  source_notes.lst <- list()
  source_notes.lst[[1]] <- character(0)

  if (!is.null(total_n)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("Total N = %s.", format(total_n, big.mark = ","))
    )
  }

  # Add family distribution info if available
  if (!is.null(family_label)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("%s distribution.", family_label)
    )
  }

  # Side legend - different for bilateral vs unilateral
  if (!is_bilateral) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      ifelse(
        is_scaled,
        "L: Left hemisphere & R: Right hemisphere.",
        "L/R: Left/Right hemisphere."
      )
    )
  }

  if (is_scaled) {
    if (is_bilateral || format == "latex") i <- 2 else i <- 1
    if (i == 2) source_notes.lst[[2]] <- character(0)
    source_notes.lst[[i]] <- c(
      source_notes.lst[[i]],
      sprintf(
        "%s %s averaged within %d-year age bins.",
        ifelse(is_bilateral, "Sum of hemispheres", "Values"),
        ifelse(format == "latex", "$\\times 10^{-3}$", "\u00d710\u207b\u00b3"),
        bin_width
      )
    )
  } else {
    if (format == "latex" && is_bilateral) i <- 2 else i <- 1
    if (i == 2) source_notes.lst[[2]] <- character(0)
    source_notes.lst[[i]] <- c(
      source_notes.lst[[i]],
      sprintf(
        "%s within %d-year age bins.",
        ifelse(is_bilateral, "Sum of hemispheres averaged", "Means"),
        bin_width
      )
    )
  }


  # Add clinical threshold information with proper LaTeX escaping
  i <- length(source_notes.lst) + 1
  source_notes.lst[[i]] <- character(0)
  if (add_clinical_labels) {
    if (roi_type %in% c("HC", "HVR")) {
      if (format == "latex") {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: $<$2.5\\% = Severe atrophy (-2SD),",
          "2.5-5\\% = Borderline (-1.6SD)."
        )
      } else {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: <2.5% = Severe atrophy (-2SD),",
          "2.5-5% = Borderline (-1.6SD)."
        )
      }
    } else if (roi_type == "LV") {
      if (format == "latex") {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: $>$97.5\\% =",
          "Abnormal enlargement (+2SD),",
          "95-97.5\\% = Borderline (+1.6SD)."
        )
      } else {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: >97.5% =",
          "Abnormal enlargement (+2SD),",
          "95-97.5% = Borderline (+1.6SD)."
        )
      }
    }
  }

  # Format-specific styling
  if (format == "latex") {
    # Formal LaTeX style: smaller font to fit wide tables, no colors,
    # clean typography
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML style with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(10),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(13),
        heading.subtitle.font.size = px(10),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Create GAMLSS model validation table
#' @param validation_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   MAE, INTERPRETATION
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
create_validation_table <- function(
    validation_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()
  interp_symbols <- get_interpretation_symbols()

  # Filter for this ROI
  val.dt <- copy(validation_data)[ROI == roi_code]
  if (nrow(val.dt) == 0) return(NULL)

  # Check if all interpretations are the same (no need for indicators)
  use_indicators <- length(unique(val.dt$INTERPRETATION)) > 1

  # Add symbols if needed
  if (use_indicators) {
    val.dt[, MAE_DISPLAY := paste0(MAE, interp_symbols[INTERPRETATION])]
  } else {
    val.dt[, MAE_DISPLAY := as.character(MAE)]
  }

  # Add labels
  val.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  val.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Reshape: rows = SEX + ADJ, columns = SIDE
  val_wide.dt <- dcast(
    val.dt, SEX + ADJ_LABEL ~ SIDE_LABEL,
    value.var = "MAE_DISPLAY"
  )

  # Column order
  side_order <- c("Left", "Right", "Bilateral")
  side_cols <- side_order[side_order %in% names(val_wide.dt)]
  setcolorder(val_wide.dt, c("SEX", "ADJ_LABEL", side_cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Model Validation (MAE)")
  gt_table <- val_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Column labels (no "MAE" label, just the side names)
  col_labels <- setNames(side_cols, side_cols)
  gt_table <- gt_table |> cols_label(.list = col_labels)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid",
        heading.align = "center"
      )
  } else {
    # HTML with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(10),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(13),
        heading.align = "center",
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Create longitudinal stability table
#' @param stability_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   CORR, MEAN_CHG, SD_CHG, CROSS
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
create_stability_table <- function(
    stability_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()

  # Filter for this ROI
  stab.dt <- copy(stability_data)[ROI == roi_code]
  if (nrow(stab.dt) == 0) return(NULL)

  # Add labels
  stab.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  stab.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Reshape: rows = SEX + ADJ, columns = SIDE x METRIC
  stab_long.dt <- melt(
    stab.dt,
    id.vars = c("SEX", "ADJ_LABEL", "SIDE_LABEL"),
    measure.vars = c("CORR", "MEAN_CHG", "SD_CHG", "CROSS"),
    variable.name = "METRIC", value.name = "VALUE"
  )
  stab_long.dt[, COL := paste(SIDE_LABEL, METRIC, sep = "_")]

  stab_wide.dt <- dcast(
    stab_long.dt, SEX + ADJ_LABEL ~ COL,
    value.var = "VALUE"
  )

  # Build column order: Left (r, Mean, SD, Cross), Right (...), Bilateral (...)
  cols <- character(0)
  for (side in c("Left", "Right", "Bilateral")) {
    for (metric in c("CORR", "MEAN_CHG", "SD_CHG", "CROSS")) {
      col <- paste(side, metric, sep = "_")
      if (col %in% names(stab_wide.dt)) cols <- c(cols, col)
    }
  }
  setcolorder(stab_wide.dt, c("SEX", "ADJ_LABEL", cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Longitudinal Stability")
  gt_table <- stab_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Add side spanners
  for (side in c("Left", "Right", "Bilateral")) {
    side_cols <- grep(paste0("^", side, "_"), cols, value = TRUE)
    if (length(side_cols) > 0) {
      gt_table <- gt_table |>
        tab_spanner(
          label = side, columns = all_of(side_cols),
          id = paste0("spanner_", side)
        )
    }
  }

  # Column labels
  labels <- list()
  for (col in cols) {
    labels[[col]] <- fcase(
      grepl("_CORR$", col), "r",
      grepl("_MEAN_CHG$", col), "Mean",
      grepl("_SD_CHG$", col), "SD",
      grepl("_CROSS$", col), "Cross"
    )
  }
  gt_table <- gt_table |> cols_label(.list = labels)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML with colors (smaller font for stability)
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(9),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(12),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Create GAMLSS model summary table
#' @param summary_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   FAMILY, AIC, BIC
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
create_model_summary_table <- function(
    summary_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()
  fam_info <- get_family_info()

  # Filter for this ROI
  summ.dt <- copy(summary_data)[ROI == roi_code]
  if (nrow(summ.dt) == 0) return(NULL)

  # Add labels
  summ.dt[, FAM_ABBREV := as.character(fam_info$abbrev[FAMILY])]
  summ.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  summ.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Convert numeric columns to character for consistent melting
  summ.dt[, AIC := as.character(AIC)]
  summ.dt[, BIC := as.character(BIC)]

  # Reshape: rows = SEX + ADJ, columns = SIDE x METRIC
  summ_long.dt <- melt(
    summ.dt,
    id.vars = c("SEX", "ADJ_LABEL", "SIDE_LABEL"),
    measure.vars = c("FAM_ABBREV", "AIC", "BIC"),
    variable.name = "METRIC", value.name = "VALUE"
  )
  summ_long.dt[, COL := paste(SIDE_LABEL, METRIC, sep = "_")]

  summ_wide.dt <- dcast(
    summ_long.dt, SEX + ADJ_LABEL ~ COL,
    value.var = "VALUE"
  )

  # Build column order: Left (Family, AIC, BIC), Right (...), Bilateral (...)
  cols <- character(0)
  for (side in c("Left", "Right", "Bilateral")) {
    for (metric in c("FAM_ABBREV", "AIC", "BIC")) {
      col <- paste(side, metric, sep = "_")
      if (col %in% names(summ_wide.dt)) cols <- c(cols, col)
    }
  }
  setcolorder(summ_wide.dt, c("SEX", "ADJ_LABEL", cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Model Summary")
  gt_table <- summ_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Add side spanners
  for (side in c("Left", "Right", "Bilateral")) {
    side_cols <- grep(paste0("^", side, "_"), cols, value = TRUE)
    if (length(side_cols) > 0) {
      gt_table <- gt_table |>
        tab_spanner(
          label = side, columns = all_of(side_cols),
          id = paste0("spanner_", side)
        )
    }
  }

  # Column labels
  labels <- list()
  for (col in cols) {
    labels[[col]] <- fcase(
      grepl("_FAM_ABBREV$", col), "Family",
      grepl("_AIC$", col), "AIC",
      grepl("_BIC$", col), "BIC"
    )
  }
  gt_table <- gt_table |> cols_label(.list = labels)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML with colors (smaller font for wide model summary)
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(7.5),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(11),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(2),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

# =============================================================================
# Analysis-Specific Table Functions
# =============================================================================

#' Create LGCM coupling results table
#'
#' Shows slope coupling (covariance between HVR
#' and cognitive slopes) from parallel process LGCM.
#'
#' @param lgcm_results LGCM parallel process results
#' @param title Table title
#' @return gt table object
#' @export
create_lgcm_coupling_table <- function(
    lgcm_results, title = NULL) {
  if (is.null(lgcm_results) ||
      is.null(lgcm_results$summary)) {
    return(NULL)
  }

  domain_map <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )
  domains <- c("MEM", "LAN", "EXF")
  sexes <- c("Female", "Male")

  # Flag degenerate estimates (SE > 100
  # = non-positive-definite Hessian)
  SE_DEGEN <- 100

  rows.lst <- lapply(sexes, function(sx) {
    lapply(domains, function(d) {
      key <- paste0(sx, "_", d)
      lr <- lgcm_results$linear_results[[key]]
      if (is.null(lr) ||
          is.null(lr$coupling)) {
        return(data.table(
          Sex = sx,
          Domain = domain_map[d],
          Beta = NA_real_,
          p = NA_real_
        ))
      }
      se <- lr$coupling$se
      degen <- !is.null(se) && !is.na(se) &&
        se > SE_DEGEN
      data.table(
        Sex = sx,
        Domain = domain_map[d],
        Beta = if (degen) NA_real_ else {
          lr$coupling$cov
        },
        p = if (degen || is.na(lr$coupling$p)) {
          NA_real_
        } else {
          lr$coupling$p
        }
      )
    }) |> rbindlist()
  }) |> rbindlist()

  if (nrow(rows.lst) == 0) return(NULL)

  # Pivot: rows = Domain, cols = Sex (beta + p)
  wide_rows <- lapply(
    domain_map, function(dom) {
      row <- list(Domain_label = dom)
      for (s in sexes) {
        ss <- rows.lst[
          Sex == s & Domain == dom
        ]
        row[[paste0(s, "_b")]] <- ss$Beta[1]
        row[[paste0(s, "_p")]] <- ss$p[1]
      }
      as.data.table(row)
    }
  ) |> rbindlist()

  tbl <- gt(
    wide_rows,
    rowname_col = "Domain_label"
  )

  sex_label.fn <- function(sx) paste0(sx, "s")

  for (s in sexes) {
    bc <- paste0(s, "_b")
    pc <- paste0(s, "_p")
    lab <- sex_label.fn(s)
    tbl <- tbl |>
      sub_missing(
        columns = all_of(c(bc, pc)),
        missing_text = "\u2014"
      ) |>
      fmt_number(
        columns = all_of(bc),
        rows = !is.na(.data[[bc]]),
        decimals = 4
      ) |>
      fmt_number(
        columns = all_of(pc),
        rows = !is.na(.data[[pc]]) &
          .data[[pc]] >= 0.001,
        decimals = 3
      ) |>
      fmt_scientific(
        columns = all_of(pc),
        rows = !is.na(.data[[pc]]) &
          .data[[pc]] < 0.001,
        decimals = 2
      ) |>
      tab_spanner(
        label = lab,
        columns = all_of(c(bc, pc))
      ) |>
      cols_label(
        !!bc := "\u03b2",
        !!pc := "p"
      ) |>
      tab_style(
        style = cell_text(
          weight = "bold"
        ),
        locations = cells_body(
          columns = all_of(c(bc, pc)),
          rows =
            !is.na(.data[[pc]]) &
            .data[[pc]] < 0.05
        )
      )
  }

  if (!is.null(title)) {
    tbl <- tbl |>
      tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' LGCM mediation model fit and LRT table
#'
#' Reports -2LL, AIC, BIC and nested LRT results for
#' each path (a, b, c') per sex/domain/predictor.
#' @param lgcm_results Mediation results list
#' @param title Table title
#' @return gt table object
#' Create mediation model fit table
#'
#' Shows -2LL, AIC, BIC per sex/domain/predictor.
#' Spanners = Sex, row groups = Domain,
#' rows = FRS / CVR MIMIC.
#' @param lgcm_results Mediation results list
#' @param title Table title
#' @return gt table object
create_mediation_fit_table <- function(
    lgcm_results, title = NULL) {
  if (is.null(lgcm_results) ||
      is.null(lgcm_results$results)) {
    return(NULL)
  }
  domains <- c("MEM", "LAN", "EXF")
  domain_map <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )
  sexes <- c("Female", "Male")
  preds <- c("FRS", "CVR_mimic")
  metrics <- c("-2LL", "AIC", "BIC")

  # Build wide: Domain | Metric | F_FRS F_CVR M_FRS M_CVR
  wide.dt <- lapply(
    domain_map, function(dom) {
      d <- names(domain_map)[domain_map == dom]
      lapply(metrics, function(met) {
        row <- list(
          Domain = dom, Metric = met
        )
        for (s in sexes) {
          for (pr in preds) {
            k <- paste(s, d, pr, sep = "_")
            r <- lgcm_results$results[[k]]
            fi <- r$fit_indices
            val <- if (!is.null(fi)) {
              switch(met,
                "-2LL" = fi$minus2LL,
                "AIC" = fi$AIC,
                "BIC" = fi$BIC
              )
            } else NA_real_
            col <- paste0(
              s, "_", gsub(
                "CVR_mimic", "CVR", pr
              )
            )
            row[[col]] <- val
          }
        }
        as.data.table(row)
      }) |> rbindlist()
    }
  ) |> rbindlist()

  if (nrow(wide.dt) == 0) return(NULL)

  tbl <- gt(
    wide.dt,
    rowname_col = "Metric",
    groupname_col = "Domain"
  )

  for (s in sexes) {
    cols <- paste0(s, c("_FRS", "_CVR"))
    tbl <- tbl |>
      fmt_number(
        columns = all_of(cols),
        decimals = 1
      ) |>
      tab_spanner(
        label = paste0(s, "s"),
        columns = all_of(cols)
      ) |>
      cols_label(
        !!cols[1] := "FRS",
        !!cols[2] := md("CVR~mimic~")
      )
  }

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }
  tbl |> style_manuscript_table()
}

#' Create mediation LRT table
#'
#' Shows a/b LRT chi-squared, ab/c' bootstrap CIs
#' per sex/domain/predictor.
#' Rows = path per domain, cols = Sex x {FRS, CVR}.
#' @param lgcm_results Mediation results list
#' @param title Table title
#' @return gt table object
create_mediation_lrt_table <- function(
    lgcm_results, title = NULL) {
  if (is.null(lgcm_results) ||
      is.null(lgcm_results$results)) {
    return(NULL)
  }
  domains <- c("MEM", "LAN", "EXF")
  domain_map <- c(
    MEM = "Memory", LAN = "Language",
    EXF = "Executive Function"
  )
  sexes <- c("Female", "Male")
  preds <- c("FRS", "CVR_mimic")
  path_rows <- c(
    "a-path", "b-path",
    "ab (indirect)", "c'-path"
  )

  # Format estimate [boot CI] cell
  fmt_est_ci.fn <- function(est, lo, hi) {
    if (is.null(est) || is.na(est) ||
        is.null(lo) || is.na(lo) ||
        is.null(hi) || is.na(hi)) {
      return(NA_character_)
    }
    sprintf(
      "%.4f [%.4f, %.4f]", est, lo, hi
    )
  }

  # Boot CI excludes zero?
  ci_sig.fn <- function(lo, hi) {
    !is.null(lo) && !is.null(hi) &&
      !is.na(lo) && !is.na(hi) &&
      (lo > 0 || hi < 0)
  }

  # Build wide: Domain | Path | F_FRS F_CVR M_FRS M_CVR
  wide.dt <- lapply(
    domain_map, function(dom) {
      d <- names(domain_map)[
        domain_map == dom
      ]
      lapply(path_rows, function(pth) {
        row <- list(
          Domain = dom, Path = pth
        )
        for (s in sexes) {
          for (pr in preds) {
            k <- paste(s, d, pr, sep = "_")
            r <- lgcm_results$results[[k]]
            col <- paste0(
              s, "_",
              gsub("CVR_mimic", "CVR", pr)
            )
            val <- if (pth == "a-path") {
              a <- r$a_path
              fmt_est_ci.fn(
                a$est, a$boot_ci_lower,
                a$boot_ci_upper
              )
            } else if (pth == "b-path") {
              b <- r$b_path
              fmt_est_ci.fn(
                b$est, b$boot_ci_lower,
                b$boot_ci_upper
              )
            } else if (
              pth == "ab (indirect)"
            ) {
              ind <- r$indirect
              fmt_est_ci.fn(
                ind$est,
                ind$boot_ci_lower,
                ind$boot_ci_upper
              )
            } else {
              cp <- r$cprime
              fmt_est_ci.fn(
                cp$est,
                cp$boot_ci_lower,
                cp$boot_ci_upper
              )
            }
            row[[col]] <- val
          }
        }
        as.data.table(row)
      }) |> rbindlist()
    }
  ) |> rbindlist()

  tbl <- gt(
    wide.dt,
    rowname_col = "Path",
    groupname_col = "Domain"
  )

  for (s in sexes) {
    cols <- paste0(s, c("_FRS", "_CVR"))
    tbl <- tbl |>
      sub_missing(
        columns = all_of(cols),
        missing_text = "\u2014"
      ) |>
      tab_spanner(
        label = paste0(s, "s"),
        columns = all_of(cols)
      ) |>
      cols_label(
        !!cols[1] := "FRS",
        !!cols[2] := md("CVR~mimic~")
      )
  }

  # Bold significant cells (CI excludes 0)
  for (i in seq_len(nrow(wide.dt))) {
    pth <- wide.dt$Path[i]
    d <- names(domain_map)[
      domain_map == wide.dt$Domain[i]
    ]
    for (s in sexes) {
      for (pr in preds) {
        k <- paste(s, d, pr, sep = "_")
        r <- lgcm_results$results[[k]]
        col <- paste0(
          s, "_",
          gsub("CVR_mimic", "CVR", pr)
        )
        is_sig <- FALSE
        if (pth == "a-path") {
          a <- r$a_path
          is_sig <- ci_sig.fn(
            a$boot_ci_lower,
            a$boot_ci_upper
          )
        } else if (pth == "b-path") {
          b <- r$b_path
          is_sig <- ci_sig.fn(
            b$boot_ci_lower,
            b$boot_ci_upper
          )
        } else if (
          pth == "ab (indirect)"
        ) {
          ind <- r$indirect
          is_sig <- ci_sig.fn(
            ind$boot_ci_lower,
            ind$boot_ci_upper
          )
        } else {
          cp <- r$cprime
          is_sig <- ci_sig.fn(
            cp$boot_ci_lower,
            cp$boot_ci_upper
          )
        }
        if (is_sig) {
          tbl <- tbl |> tab_style(
            style = cell_text(
              weight = "bold"
            ),
            locations = cells_body(
              columns = all_of(col),
              rows = i
            )
          )
        }
      }
    }
  }

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }
  tbl |> style_manuscript_table()
}

#' Create mediation results table
#'
#' @param mediation_results List with mediation analysis results
#' @param title Table title
#' @return gt table object
#' @export
create_mediation_table <- function(mediation_results, title = NULL) {
  if (is.null(mediation_results) || is.null(mediation_results$results)) {
    return(NULL)
  }

  # Extract results into data.table
  res.lst <- lapply(names(mediation_results$results), function(nm) {
    r <- mediation_results$results[[nm]]
    if (is.null(r) || !isTRUE(r$converged)) return(NULL)

    # Parse name: "Sample_Domain_Predictor"
    # or legacy "Sample_Domain"
    parts <- strsplit(nm, "_")[[1]]
    sample <- parts[1]
    domain <- if (length(parts) > 1) {
      parts[2]
    } else {
      "Unknown"
    }
    predictor <- if (length(parts) > 2) {
      paste(parts[3:length(parts)],
            collapse = "_")
    } else {
      NA_character_
    }

    # Extract path coefficients
    a_est <- r$a_path$est %||% NA
    a_p <- r$a_path$p %||% NA
    b_est <- r$b_path$est %||% NA
    b_p <- r$b_path$p %||% NA
    ind_est <- r$indirect$est %||% NA
    ci_lo <- r$indirect$boot_ci_lower %||% NA
    ci_hi <- r$indirect$boot_ci_upper %||% NA
    prop <- r$proportion_mediated %||%
      r$prop_mediated %||% NA

    data.table(
      Sample = sample,
      Domain = domain,
      Predictor = predictor,
      Path_a = a_est,
      Path_a_p = a_p,
      Path_b = b_est,
      Path_b_p = b_p,
      Indirect = ind_est,
      CI_lower = ci_lo,
      CI_upper = ci_hi,
      Proportion = prop
    )
  })

  summ <- rbindlist(res.lst[!sapply(res.lst, is.null)])
  if (nrow(summ) == 0) return(NULL)

  # Map domain codes to labels
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  summ[, Domain := fifelse(Domain %in% names(domain_map),
                            domain_map[Domain], Domain)]

  # Convert Male/Female to Males/Females
  summ[, Sample := fifelse(Sample == "Male", "Males",
                    fifelse(Sample == "Female", "Females",
                    fifelse(Sample == "Full", "Full Sample", Sample)))]

  # Format proportion: cap at 100%, show "Full" if >= 100%
  summ[, Prop_display := fifelse(Proportion >= 1.0, "Full",
                                  sprintf("%.0f%%", Proportion * 100))]

  # Hide Predictor if all NA (legacy format)
  has_pred <- "Predictor" %in% names(summ) &&
    any(!is.na(summ$Predictor))

  tbl <- summ |>
    gt(
      groupname_col = "Sample",
      rowname_col = "Domain"
    ) |>
    fmt_number(
      columns = c(
        Path_a, Path_b, Indirect,
        CI_lower, CI_upper
      ),
      decimals = 4
    ) |>
    fmt_scientific(
      columns = c(Path_a_p, Path_b_p),
      decimals = 2
    ) |>
    cols_merge(
      columns = c(CI_lower, CI_upper),
      pattern = "[{1}, {2}]"
    ) |>
    cols_hide(columns = "Proportion") |>
    cols_label(
      Path_a = "\u03b2",
      Path_a_p = "p",
      Path_b = "\u03b2",
      Path_b_p = "p",
      Indirect = "a \u00d7 b",
      CI_lower = "95% CI",
      Prop_display = "% Mediated"
    ) |>
    tab_spanner(
      label = "Path a",
      columns = c(Path_a, Path_a_p)
    ) |>
    tab_spanner(
      label = "Path b",
      columns = c(Path_b, Path_b_p)
    ) |>
    tab_spanner(
      label = "Indirect Effect",
      columns = c(
        Indirect, CI_lower, Prop_display
      )
    )
  if (!has_pred) {
    tbl <- tbl |>
      cols_hide(columns = "Predictor")
  }

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

# =============================================================================
# FRS vs CVR_mimic Comparison Tables
# =============================================================================

#' Create LME 3-way interaction results table
#'
#' Formats LME results for manuscript tables, showing sex-stratified
#' YRS × CVR × Brain interaction effects.
#'
#' @param lme_results List with LME results (from 06a/06b scripts)
#' @param title Table title
#' @param subtitle Table subtitle (optional)
#' @return gt table object
#' @export
create_lme_results_table <- function(lme_results, title = NULL,
                                     subtitle = NULL) {
  if (is.null(lme_results) || is.null(lme_results$summary_table)) return(NULL)

  summ <- as.data.table(lme_results$summary_table)

  # Filter to primary stratified analyses only
  if ("Analysis" %in% names(summ)) {
    summ <- summ[Analysis == "Primary_Stratified"]
  }

  if (nrow(summ) == 0) return(NULL)

  # Map domain codes to labels
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  if ("Domain" %in% names(summ)) {
    summ[, Domain := fifelse(Domain %in% names(domain_map),
                              domain_map[Domain], Domain)]
  }

  # Get CVR measure label
  cvr_label <- if (!is.null(lme_results$cvr_measure)) {
    lme_results$cvr_measure
  } else {
    "FRS"
  }

  brain_label <- if (!is.null(lme_results$brain_measure)) {
    lme_results$brain_measure
  } else {
    "Brain"
  }

  # Format p-values with significance stars
  summ[, p_display := format_p(p_value)]
  summ[, p_stars := sig_stars(p_value)]

  # Format beta with SE
  summ[, Effect := sprintf("%.4f (%.4f)%s", Beta, SE, p_stars)]

  # Create display table
  display.dt <- dcast(summ, Domain ~ Sex, value.var = "Effect")

  # Create gt table
  tbl <- display.dt |>
    gt(rowname_col = "Domain") |>
    tab_stubhead(label = "Outcome") |>
    cols_label(
      Male = "Males",
      Female = "Females"
    ) |>
    tab_spanner(
      label = sprintf("YRS \u00d7 %s \u00d7 %s", cvr_label, brain_label),
      columns = c(Male, Female)
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title, subtitle = subtitle)
  }

  tbl <- tbl |>
    tab_source_note(
      source_note = paste(
        "Values are \u03b2 (SE) with significance stars.",
        "* p < 0.05, ** p < 0.01, *** p < 0.001"
      )
    )

  tbl |> style_manuscript_table()
}

#' Create side-by-side FRS vs CVR_mimic comparison table
#'
#' Compares 3-way interaction effects between FRS and CVR_mimic
#' to demonstrate age-confounding in FRS.
#'
#' @param frs_results List with FRS LME results
#' @param cvr_results List with CVR_mimic LME results
#' @param title Table title
#' @param subtitle Table subtitle
#' @return gt table object
#' @export
create_comparison_table <- function(frs_results, cvr_results,
                                    title = NULL, subtitle = NULL) {
  if (is.null(frs_results) || is.null(cvr_results)) return(NULL)

  frs_summ <- as.data.table(frs_results$summary_table)
  cvr_summ <- as.data.table(cvr_results$summary_table)

  # Filter to primary stratified analyses
  if ("Analysis" %in% names(frs_summ)) {
    frs_summ <- frs_summ[Analysis == "Primary_Stratified"]
  }
  if ("Analysis" %in% names(cvr_summ)) {
    cvr_summ <- cvr_summ[Analysis == "Primary_Stratified"]
  }

  if (nrow(frs_summ) == 0 || nrow(cvr_summ) == 0) return(NULL)

  # Map domain codes to labels
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  frs_summ[, Domain := fifelse(Domain %in% names(domain_map),
                                domain_map[Domain], Domain)]
  cvr_summ[, Domain := fifelse(Domain %in% names(domain_map),
                                domain_map[Domain], Domain)]

  # Format effects with stars
  frs_summ[, FRS_effect := sprintf("%.4f%s", Beta, sig_stars(p_value))]
  frs_summ[, FRS_p := format_p(p_value)]

  cvr_summ[, CVR_effect := sprintf("%.4f%s", Beta, sig_stars(p_value))]
  cvr_summ[, CVR_p := format_p(p_value)]

  # Merge FRS and CVR results
  combined.dt <- merge(
    frs_summ[, .(Sex, Domain, FRS_effect, FRS_p)],
    cvr_summ[, .(Sex, Domain, CVR_effect, CVR_p)],
    by = c("Sex", "Domain")
  )

  # Convert Male/Female to Males/Females
  combined.dt[, Sex := fifelse(Sex == "Male", "Males",
                        fifelse(Sex == "Female", "Females", Sex))]

  # Create gt table with row groups by sex
  tbl <- combined.dt |>
    gt(groupname_col = "Sex", rowname_col = "Domain") |>
    tab_stubhead(label = "Outcome") |>
    cols_label(
      FRS_effect = "\u03b2",
      FRS_p = "p",
      CVR_effect = "\u03b2",
      CVR_p = "p"
    ) |>
    tab_spanner(label = "FRS", columns = c(FRS_effect, FRS_p)) |>
    tab_spanner(label = "CVR (age-adj.)", columns = c(CVR_effect, CVR_p))

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title, subtitle = subtitle)
  }

  tbl <- tbl |>
    tab_source_note(
      source_note = paste(
        "FRS = Framingham Risk Score;",
        "CVR = cardiovascular risk",
        "from MIMIC model.",
        "* p < 0.05, ** p < 0.01, *** p < 0.001"
      )
    )

  tbl |> style_manuscript_table()
}

#' Create LGCM mediation comparison table
#'
#' Compares mediation effects (a-path, b-path, indirect) between
#' FRS and CVR_mimic predictors.
#'
#' @param frs_mediation List with FRS LGCM mediation results
#' @param cvr_mediation List with CVR_mimic LGCM mediation results
#' @param title Table title
#' @return gt table object
#' @export
create_mediation_comparison_table <- function(frs_mediation, cvr_mediation,
                                               title = NULL) {
  if (is.null(frs_mediation) || is.null(cvr_mediation)) return(NULL)

  # Helper to extract mediation results
  extract_mediation.fn <- function(results, predictor_label) {
    if (is.null(results$results)) return(NULL)

    res.lst <- lapply(names(results$results), function(nm) {
      r <- results$results[[nm]]
      if (is.null(r) || !isTRUE(r$converged)) return(NULL)

      parts <- strsplit(nm, "_")[[1]]
      sample <- parts[1]
      domain <- if (length(parts) > 1) parts[2] else "Unknown"

      data.table(
        Predictor = predictor_label,
        Sample = sample,
        Domain = domain,
        a_est = r$a_path$est %||% NA,
        a_p = r$a_path$p %||% NA,
        b_est = r$b_path$est %||% NA,
        b_p = r$b_path$p %||% NA,
        ind_est = r$indirect$est %||% NA,
        ind_p = r$indirect$p %||% NA
      )
    })

    rbindlist(res.lst[!sapply(res.lst, is.null)])
  }

  frs_dt <- extract_mediation.fn(frs_mediation, "FRS")
  cvr_dt <- extract_mediation.fn(cvr_mediation, "CVR_mimic")

  if (is.null(frs_dt) || is.null(cvr_dt)) return(NULL)
  if (nrow(frs_dt) == 0 || nrow(cvr_dt) == 0) return(NULL)

  combined.dt <- rbind(frs_dt, cvr_dt)

  # Map domain codes
  domain_map <- c(MEM = "Memory", LAN = "Language", EXF = "Executive")
  combined.dt[, Domain := fifelse(Domain %in% names(domain_map),
                                   domain_map[Domain], Domain)]

  # Filter to Full sample or specific sex
  combined.dt <- combined.dt[Sample == "Full"]

  # Format columns
  combined.dt[, a_display := sprintf("%.4f%s", a_est, sig_stars(a_p))]
  combined.dt[, b_display := sprintf("%.4f%s", b_est, sig_stars(b_p))]
  combined.dt[, ind_display := sprintf("%.4f%s", ind_est, sig_stars(ind_p))]

  # Reshape for side-by-side comparison
  display.dt <- dcast(
    combined.dt,
    Domain ~ Predictor,
    value.var = c("a_display", "b_display", "ind_display")
  )

  # Reorder columns
  setcolorder(display.dt, c(
    "Domain",
    "a_display_FRS", "a_display_CVR_mimic",
    "b_display_FRS", "b_display_CVR_mimic",
    "ind_display_FRS", "ind_display_CVR_mimic"
  ))

  tbl <- display.dt |>
    gt(rowname_col = "Domain") |>
    tab_stubhead(label = "Outcome") |>
    cols_label(
      a_display_FRS = "FRS",
      a_display_CVR_mimic = "CVR",
      b_display_FRS = "FRS",
      b_display_CVR_mimic = "CVR",
      ind_display_FRS = "FRS",
      ind_display_CVR_mimic = "CVR"
    ) |>
    tab_spanner(
      label = "a-path (Pred \u2192 HVR)",
      columns = c(a_display_FRS, a_display_CVR_mimic)
    ) |>
    tab_spanner(
      label = "b-path (HVR \u2192 Cog)",
      columns = c(b_display_FRS, b_display_CVR_mimic)
    ) |>
    tab_spanner(
      label = "Indirect (a \u00d7 b)",
      columns = c(ind_display_FRS, ind_display_CVR_mimic)
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl <- tbl |>
    tab_source_note(
      source_note = paste(
        "FRS = Framingham Risk Score; CVR = age-adjusted cardiovascular risk.",
        "* p < 0.05, ** p < 0.01, *** p < 0.001"
      )
    )

  tbl |> style_manuscript_table()
}

# =============================================================================
# CVR Comparison Manuscript Table Functions
# =============================================================================

#' Extract LME 3-way interaction effects from results
#'
#' Extracts sex-stratified Time x CVR x Brain interaction terms
#' from LME results objects. Used by table and figure functions.
#'
#' @param res LME results list (with $stratified element)
#' @param predictor Label for predictor (e.g., "FRS", "CVR")
#' @param domain_labels Named vector mapping codes to labels
#' @return data.table with Sex, Domain, Predictor, Beta, SE, p
#' @export
extract_lme_effects <- function(
    res, predictor,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(res$stratified)) return(data.table())
  lapply(c("Male", "Female"), function(sex) {
    lapply(c("MEM", "LAN", "EXF"), function(dom) {
      r <- res$stratified[[sex]][[dom]]
      if (!is.null(r) && !is.null(r$interaction)) {
        data.table(
          Sex = sex,
          Domain = domain_labels[dom],
          Predictor = predictor,
          Beta = r$interaction$beta,
          SE = r$interaction$se,
          p = r$interaction$p_value,
          N = r$n_subjects
        )
      }
    }) |> rbindlist()
  }) |> rbindlist()
}

#' Create LME comparison table (FRS vs CVR)
#'
#' Builds a gt table comparing 3-way interaction effects
#' between two LME result sets (e.g., FRS vs CVR_mimic).
#' Reusable for both HVR (primary) and HC (sensitivity).
#'
#' @param frs_res LME results list for FRS predictor
#' @param cvr_res LME results list for CVR predictor
#' @param domain_labels Named vector for domain labels
#' @param title Table title
#' @return gt table object
#' @export
create_lme_comparison_table <- function(
    frs_res, cvr_res,
    title = NULL,
    summary_dt = NULL,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(frs_res) || is.null(cvr_res)) {
    return(NULL)
  }

  frs_lme <- extract_lme_effects(
    frs_res, "FRS", domain_labels
  )
  cvr_lme <- extract_lme_effects(
    cvr_res, "CVR", domain_labels
  )

  if (nrow(frs_lme) == 0 ||
      nrow(cvr_lme) == 0) return(NULL)

  # Build wide: Domain | Predictor |
  #   Female_b Female_se Female_p Female_pfdr
  #   Male_b Male_se Male_p Male_pfdr
  both <- rbind(frs_lme, cvr_lme)
  sexes <- c("Female", "Male")
  preds <- c("FRS", "CVR")
  dom_cols <- unique(both$Domain)
  # Reverse domain label lookup
  inv_dom <- setNames(
    names(domain_labels), domain_labels
  )

  fmt_p.fn <- function(x) {
    ifelse(
      x < 0.001, "< 0.001",
      sprintf("%.3f", x)
    )
  }

  # FDR lookup helper
  get_fdr <- function(cvr_meas, sex, dom_lab) {
    if (is.null(summary_dt)) return(NA_real_)
    dom_code <- inv_dom[dom_lab]
    if (is.na(dom_code)) return(NA_real_)
    row <- summary_dt[
      Analysis == "Primary_Stratified" &
      CVR_Measure == cvr_meas &
      Sex == sex & Domain == dom_code
    ]
    if (nrow(row) > 0 &&
        "p_fdr" %in% names(row)) {
      row$p_fdr[1]
    } else NA_real_
  }

  has_fdr <- !is.null(summary_dt) &&
    "p_fdr" %in% names(summary_dt)

  wide_rows <- lapply(dom_cols, function(d) {
    lapply(preds, function(pr) {
      cvr_meas <- if (pr == "FRS") {
        "FRS"
      } else "CVR_mimic"
      row <- list(
        Domain = d, Predictor = pr
      )
      for (s in sexes) {
        dd <- both[
          Sex == s & Domain == d &
            Predictor == pr
        ]
        if (nrow(dd) > 0) {
          row[[paste0(s, "_b")]] <-
            dd$Beta[1]
          row[[paste0(s, "_se")]] <-
            dd$SE[1]
          row[[paste0(s, "_p")]] <-
            dd$p[1]
          if (has_fdr) {
            row[[paste0(s, "_pfdr")]] <-
              get_fdr(cvr_meas, s, d)
          }
        } else {
          row[[paste0(s, "_b")]] <- NA
          row[[paste0(s, "_se")]] <- NA
          row[[paste0(s, "_p")]] <- NA
          if (has_fdr) {
            row[[paste0(s, "_pfdr")]] <- NA
          }
        }
      }
      as.data.table(row)
    }) |> rbindlist()
  }) |> rbindlist()

  # Use md subscript notation in data
  wide_rows[
    Predictor == "CVR",
    Predictor := "CVR~mimic~"
  ]

  tbl <- gt(
    wide_rows,
    rowname_col = "Predictor",
    groupname_col = "Domain"
  ) |>
    fmt_markdown(columns = "Predictor")

  for (s in sexes) {
    bc <- paste0(s, "_b")
    sc <- paste0(s, "_se")
    pc <- paste0(s, "_p")
    cols_s <- c(bc, sc, pc)
    if (has_fdr) {
      fc <- paste0(s, "_pfdr")
      cols_s <- c(cols_s, fc)
    }
    tbl <- tbl |>
      fmt_number(
        columns = all_of(c(bc, sc)),
        decimals = 4
      ) |>
      fmt(
        columns = all_of(pc),
        fns = fmt_p.fn
      )
    if (has_fdr) {
      tbl <- tbl |>
        fmt(
          columns = all_of(fc),
          fns = fmt_p.fn
        )
    }
    tbl <- tbl |>
      tab_spanner(
        label = paste0(s, "s"),
        columns = all_of(cols_s)
      ) |>
      cols_label(
        !!bc := "\u03b2",
        !!sc := "SE",
        !!pc := "p"
      )
    if (has_fdr) {
      tbl <- tbl |>
        cols_label(!!fc := md("p~FDR~"))
    }
    tbl <- tbl |>
      tab_style(
        style = cell_text(
          weight = "bold"
        ),
        locations = cells_body(
          columns = all_of(
            c(bc, sc, pc)
          ),
          rows = .data[[pc]] < 0.05
        )
      )
  }

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' Create LME covariate effects table
#'
#' Extracts covariate effects from LME model results
#' with domain spanners and sex row groups with N.
#'
#' @param lme_res LME results list ($stratified)
#' @param cvr_label Label for CVR predictor column
#' @param title Table title
#' @param domain_labels Named vector for domain labels
#' @return gt table object
#' @export
create_lme_covariate_table <- function(
    lme_res,
    cvr_label = "FRS",
    title = NULL,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(lme_res) ||
      is.null(lme_res$stratified)) {
    return(NULL)
  }

  # Flexible covariate patterns + labels
  COVAR_PATS <- list(
    Age = "^AGE(_c)?$",
    Education = "^EDUC$",
    "APOE e4+" = "^APOE4(1|TRUE)$",
    CVR = "^(FRS|CVR_mimic)(_z)?$",
    `HVR z-score` = "^HVR_z$"
  )
  covar_names.v <- names(COVAR_PATS)
  # Replace generic CVR label
  covar_names.v[
    covar_names.v == "CVR"
  ] <- cvr_label

  rows.lst <- lapply(
    c("Male", "Female"), function(sex) {
      lapply(
        c("MEM", "LAN", "EXF"),
        function(dom) {
          r <- lme_res$stratified[[sex]][[dom]]
          if (is.null(r) || is.null(r$model) ||
              is.null(r$model$coefficients)) {
            return(NULL)
          }
          coefs.mat <- r$model$coefficients
          rn <- rownames(coefs.mat)
          n_obs <- if (
            !is.null(r$model$dims)
          ) r$model$dims$N[1]
          else if (
            !is.null(r$n)
          ) r$n else NA_integer_

          lapply(
            seq_along(COVAR_PATS),
            function(i) {
              idx <- grep(
                COVAR_PATS[[i]], rn
              )
              if (length(idx) == 0) {
                return(NULL)
              }
              data.table(
                Sex = sex,
                N = n_obs,
                Domain = dom,
                Covariate = covar_names.v[i],
                Beta = coefs.mat[
                  idx[1], "Estimate"
                ],
                SE = coefs.mat[
                  idx[1], "Std. Error"
                ],
                p = coefs.mat[
                  idx[1], "Pr(>|t|)"
                ]
              )
            }
          ) |> rbindlist()
        }
      ) |> rbindlist()
    }
  ) |> rbindlist()

  if (nrow(rows.lst) == 0) return(NULL)

  # p = 0 from lme → replace with machine eps
  rows.lst[
    !is.na(p) & p < .Machine$double.eps,
    p := .Machine$double.eps
  ]

  # Sex spanner labels with N
  n_by_sex <- rows.lst[
    , .(N = N[1]), by = Sex
  ]
  sex_label.fn <- function(sx) {
    n <- n_by_sex[Sex == sx, N]
    plural <- paste0(sx, "s")
    if (length(n) == 0 || is.na(n)) {
      plural
    } else {
      sprintf(
        "%s (N = %s)", plural,
        format(n, big.mark = ",")
      )
    }
  }

  # Domain row groups
  rows.lst[
    , Domain_label := domain_labels[Domain]
  ]

  # Pivot: Covariates x sex columns
  sex_codes <- c("Female", "Male")
  dom_codes <- c("MEM", "LAN", "EXF")
  wide_rows <- lapply(
    dom_codes, function(dom) {
      lapply(
        unique(rows.lst$Covariate),
        function(cv) {
          sub <- rows.lst[
            Domain == dom &
              Covariate == cv
          ]
          row <- list(
            Domain_label =
              domain_labels[dom],
            Covariate = cv
          )
          for (s in sex_codes) {
            ss <- sub[Sex == s]
            if (nrow(ss) > 0) {
              row[[paste0(s, "_b")]] <-
                ss$Beta[1]
              row[[paste0(s, "_se")]] <-
                ss$SE[1]
              row[[paste0(s, "_p")]] <-
                ss$p[1]
            } else {
              row[[paste0(s, "_b")]] <- NA
              row[[paste0(s, "_se")]] <- NA
              row[[paste0(s, "_p")]] <- NA
            }
          }
          as.data.table(row)
        }
      ) |> rbindlist()
    }
  ) |> rbindlist()

  tbl <- gt(
    wide_rows,
    rowname_col = "Covariate",
    groupname_col = "Domain_label"
  ) |>
    fmt_markdown(columns = "Covariate")

  for (s in sex_codes) {
    bc <- paste0(s, "_b")
    sc <- paste0(s, "_se")
    pc <- paste0(s, "_p")
    lab <- sex_label.fn(s)
    tbl <- tbl |>
      fmt_number(
        columns = all_of(c(bc, sc)),
        decimals = 4
      ) |>
      fmt_number(
        columns = all_of(pc),
        rows = !is.na(.data[[pc]]) &
          .data[[pc]] >= 0.001,
        decimals = 3
      ) |>
      fmt_scientific(
        columns = all_of(pc),
        rows = !is.na(.data[[pc]]) &
          .data[[pc]] < 0.001,
        decimals = 2
      ) |>
      tab_spanner(
        label = lab,
        columns = all_of(c(bc, sc, pc))
      ) |>
      cols_label(
        !!bc := "\u03b2",
        !!sc := "SE",
        !!pc := "p"
      ) |>
      tab_style(
        style = cell_text(
          weight = "bold"
        ),
        locations = cells_body(
          columns = all_of(
            c(bc, sc, pc)
          ),
          rows =
            !is.na(.data[[pc]]) &
            .data[[pc]] < 0.05
        )
      )
  }

  if (!is.null(title)) {
    tbl <- tbl |>
      tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' Create EDT comparison table (FRS vs CVR)
#'
#' Builds a gt table comparing EDT 3-way moderation effects
#' between FRS and CVR_mimic predictors.
#'
#' @param frs_edt EDT results for FRS predictor
#' @param cvr_edt EDT results for CVR predictor
#' @param domain_labels Named vector for domain labels
#' @param title Table title
#' @return gt table object
#' @export
create_edt_comparison_table <- function(
    frs_edt, cvr_edt,
    title = NULL,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(frs_edt) || is.null(cvr_edt)) return(NULL)

  extract_edt.fn <- function(res, predictor) {
    if (is.null(res$continuous)) return(data.table())
    lapply(c("Male", "Female"), function(sex) {
      lapply(c("MEM", "LAN", "EXF"), function(dom) {
        r <- res$continuous[[sex]][[dom]]$interaction_3way
        if (!is.null(r)) {
          data.table(
            Sex = sex, Domain = domain_labels[dom],
            Beta = r$beta, SE = r$se,
            p = r$p_value, Predictor = predictor
          )
        }
      }) |> rbindlist()
    }) |> rbindlist()
  }

  frs_edt.dt <- extract_edt.fn(frs_edt, "FRS")
  cvr_edt.dt <- extract_edt.fn(cvr_edt, "CVR")

  if (nrow(frs_edt.dt) == 0 || nrow(cvr_edt.dt) == 0) {
    return(NULL)
  }

  fmt_beta.fn <- function(b, se, p) {
    stars <- sig_stars(p)
    sprintf("%.5f (%.5f)%s", b, se, stars)
  }

  edt_wide <- merge(
    frs_edt.dt[, .(
      Sex, Domain,
      FRS_beta = fmt_beta.fn(Beta, SE, p), FRS_p = p
    )],
    cvr_edt.dt[, .(
      Sex, Domain,
      CVR_beta = fmt_beta.fn(Beta, SE, p), CVR_p = p
    )],
    by = c("Sex", "Domain")
  )

  tbl <- edt_wide |>
    gt(rowname_col = "Domain", groupname_col = "Sex") |>
    cols_hide(columns = c(FRS_p, CVR_p)) |>
    tab_spanner(
      label = "\u03b2 (SE)",
      columns = c(FRS_beta, CVR_beta)
    ) |>
    cols_label(
      FRS_beta = "FRS",
      CVR_beta = "CVR"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = FRS_beta, rows = FRS_p < 0.05
      )
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = CVR_beta, rows = CVR_p < 0.05
      )
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' Create attrition analysis table
#'
#' Compares baseline characteristics between completers
#' and dropouts.
#'
#' @param attrition_res Attrition analysis results list
#' @param title Table title
#' @return gt table object
#' @export
create_attrition_table <- function(
    attrition_res, title = NULL) {
  if (is.null(attrition_res) ||
      is.null(attrition_res$comparison_table) ||
      nrow(attrition_res$comparison_table) == 0) {
    return(NULL)
  }

  attr_comp <- copy(attrition_res$comparison_table)

  var_map <- c(
    SEX = "Sex (% Male)",
    AGE = "Age", EDUC = "Education", FRS = "FRS",
    HVR_z = "HVR (z-score)", HC_z = "HC (z-score)",
    MEM = "Memory", LAN = "Language", EXF = "Executive",
    APOE4 = "APOE \u03b54 (%)"
  )
  # Preferred order: Sex first
  VAR_ORDER <- names(var_map)
  if ("Variable" %in% names(attr_comp)) {
    attr_comp[, Variable := fifelse(
      Variable %in% names(var_map),
      var_map[Variable], Variable
    )]
    # Reorder: Sex first, then Age, etc.
    attr_comp[, order_idx := match(
      Variable, var_map, nomatch = 999L
    )]
    setorder(attr_comp, order_idx)
    attr_comp[, order_idx := NULL]
  }

  tbl <- attr_comp |>
    gt(rowname_col = "Variable") |>
    sub_missing(
      columns = everything(), missing_text = "\u2014"
    ) |>
    fmt_number(
      columns = any_of(c(
        "Completers_mean", "Completers_sd",
        "Dropouts_mean", "Dropouts_sd", "Diff"
      )),
      decimals = 2
    ) |>
    fmt_number(
      columns = any_of("t_statistic"), decimals = 2
    ) |>
    fmt_scientific(
      columns = any_of("p_value"), decimals = 2
    ) |>
    cols_hide(columns = any_of("Test")) |>
    cols_label(
      Completers_mean = "Mean/%",
      Completers_sd = "SD",
      Dropouts_mean = "Mean/%",
      Dropouts_sd = "SD",
      Diff = "\u0394",
      t_statistic = "\u03c7\u00b2/t",
      p_value = "p"
    ) |>
    tab_spanner(
      label = "Completers",
      columns = c(Completers_mean, Completers_sd)
    ) |>
    tab_spanner(
      label = "Dropouts",
      columns = c(Dropouts_mean, Dropouts_sd)
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(rows = p_value < 0.05)
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' Create EDT linearity test table
#'
#' Shows linear vs quadratic model comparison for EDT.
#'
#' @param edt_linearity EDT linearity test results list
#' @param domain_labels Named vector for domain labels
#' @param title Table title
#' @return gt table object
#' @export
create_linearity_table <- function(
    edt_linearity,
    title = NULL,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(edt_linearity) ||
      is.null(edt_linearity$results_table)) {
    return(NULL)
  }

  lin_dt <- copy(edt_linearity$results_table)

  if ("Domain" %in% names(lin_dt)) {
    lin_dt[, Domain := fifelse(
      Domain %in% names(domain_labels),
      domain_labels[Domain], Domain
    )]
  }

  cols_to_show <- intersect(
    c("Domain", "N_obs", "AIC_diff",
      "LRT_chisq", "LRT_df", "LRT_p"),
    names(lin_dt)
  )
  lin_dt <- lin_dt[, ..cols_to_show]

  tbl <- lin_dt |>
    gt(rowname_col = "Domain") |>
    fmt_number(
      columns = any_of(c("AIC_diff", "LRT_chisq")),
      decimals = 2
    ) |>
    fmt_integer(columns = any_of(c("N_obs", "LRT_df"))) |>
    fmt(
      columns = any_of("LRT_p"),
      fns = function(x) {
        ifelse(x == 0, "< 2.2e-16",
          ifelse(x < 0.001,
            sprintf("%.2e", x),
            sprintf("%.3f", x)))
      }
    ) |>
    cols_label(
      N_obs = "N", AIC_diff = "\u0394AIC",
      LRT_chisq = "\u03c7\u00b2",
      LRT_df = "df", LRT_p = "p"
    ) |>
    tab_spanner(
      label = "Likelihood Ratio Test",
      columns = c(LRT_chisq, LRT_df, LRT_p)
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}

#' Create LGCM mediation comparison table (gt)
#'
#' Compares FRS vs CVR_mimic mediation paths (a-path, b-path,
#' indirect) across cognitive domains. Designed for the CVR
#' comparison manuscript.
#'
#' @param lgcm_results LGCM CVR comparison results list
#' @param domain_labels Named vector for domain labels
#' @param title Table title
#' @return gt table object
#' @export
create_mediation_comparison_table_gt <- function(
    lgcm_results, title = NULL,
    domain_labels = c(
      MEM = "Memory", LAN = "Language",
      EXF = "Executive Function"
    )) {
  if (is.null(lgcm_results) ||
      is.null(lgcm_results$results)) {
    return(NULL)
  }

  domains <- c("MEM", "LAN", "EXF")

  fmt_est.fn <- function(est, p) {
    stars <- sig_stars(p)
    p_fmt <- if (p < 0.001) {
      "< .001"
    } else {
      sprintf("%.3f", p)
    }
    sprintf("%.4f%s (%s)", est, stars, p_fmt)
  }

  fmt_ci.fn <- function(lo, hi) {
    if (is.null(lo) || is.na(lo) ||
        is.null(hi) || is.na(hi)) {
      return("\u2014")
    }
    lo_fmt <- if (
      abs(lo) < 1e-04 && lo != 0
    ) {
      sprintf("%.2e", lo)
    } else {
      sprintf("%.4f", lo)
    }
    sprintf("%s, %.4f", lo_fmt, hi)
  }

  # Detect key format: per-sex ("MEM_FRS") or

  # full ("Male_MEM_FRS")
  keys.v <- names(lgcm_results$results)
  has_sex_pfx <- any(
    grepl("^(Male|Female)_", keys.v)
  )
  if (has_sex_pfx) {
    sexes.v <- c("Male", "Female")
  } else {
    sexes.v <- ""
  }

  med_rows <- lapply(sexes.v, function(sx) {
    pfx <- if (nchar(sx) > 0) {
      paste0(sx, "_")
    } else {
      ""
    }
    lapply(domains, function(d) {
    frs_k <- paste0(pfx, d, "_FRS")
    cvr_k <- paste0(pfx, d, "_CVR")
    frs <- lgcm_results$results[[frs_k]]
    cvr <- lgcm_results$results[[cvr_k]]
    # Try CVR_mimic suffix (raw format)
    if (is.null(cvr)) {
      cvr_k2 <- paste0(
        pfx, d, "_CVR_mimic"
      )
      cvr <- lgcm_results$results[[cvr_k2]]
    }
    if (is.null(frs) || is.null(cvr)) {
      return(NULL)
    }

    # Remap a_path->a, b_path->b if needed
    if (is.null(frs$a) && !is.null(frs$a_path)) {
      frs$a <- frs$a_path
      frs$b <- frs$b_path
    }
    if (is.null(cvr$a) && !is.null(cvr$a_path)) {
      cvr$a <- cvr$a_path
      cvr$b <- cvr$b_path
    }

    # Handle non-converged results
    frs_ok <- isTRUE(frs$converged)
    cvr_ok <- isTRUE(cvr$converged)
    dash <- "\u2014"

    row <- data.table(
      Domain = domain_labels[d],
      FRS_a_est = if (frs_ok) {
        fmt_est.fn(frs$a$est, frs$a$p)
      } else {
        dash
      },
      FRS_a_p = if (frs_ok) {
        frs$a$p
      } else {
        NA_real_
      },
      CVR_a_est = if (cvr_ok) {
        fmt_est.fn(cvr$a$est, cvr$a$p)
      } else {
        dash
      },
      CVR_a_p = if (cvr_ok) {
        cvr$a$p
      } else {
        NA_real_
      },
      FRS_b_est = if (frs_ok) {
        fmt_est.fn(frs$b$est, frs$b$p)
      } else {
        dash
      },
      FRS_b_p = if (frs_ok) {
        frs$b$p
      } else {
        NA_real_
      },
      CVR_b_est = if (cvr_ok) {
        fmt_est.fn(cvr$b$est, cvr$b$p)
      } else {
        dash
      },
      CVR_b_p = if (cvr_ok) {
        cvr$b$p
      } else {
        NA_real_
      },
      FRS_ind = if (frs_ok) {
        sprintf(
          "%.4f [%s]",
          frs$indirect$est,
          fmt_ci.fn(
            frs$indirect$boot_ci_lower,
            frs$indirect$boot_ci_upper
          )
        )
      } else {
        dash
      },
      FRS_ind_sig = if (frs_ok) {
        frs$indirect$boot_significant
      } else {
        NA
      },
      CVR_ind = if (cvr_ok) {
        sprintf(
          "%.4f [%s]",
          cvr$indirect$est,
          fmt_ci.fn(
            cvr$indirect$boot_ci_lower,
            cvr$indirect$boot_ci_upper
          )
        )
      } else {
        dash
      },
      CVR_ind_sig = if (cvr_ok) {
        cvr$indirect$boot_significant
      } else {
        NA
      }
    )
    if (nchar(sx) > 0) row[, Sex := sx]
    row
    }) |> rbindlist()
  }) |> rbindlist()

  med_dt <- med_rows
  if (is.null(med_dt) || nrow(med_dt) == 0) {
    return(NULL)
  }

  grp_col <- if (
    "Sex" %in% names(med_dt)
  ) "Sex" else NULL
  tbl <- med_dt |>
    gt(
      rowname_col = "Domain",
      groupname_col = grp_col
    ) |>
    cols_hide(columns = c(
      FRS_a_p, CVR_a_p, FRS_b_p, CVR_b_p,
      FRS_ind_sig, CVR_ind_sig
    )) |>
    tab_spanner(
      label = "a-path (CVR \u2192 HVR slope)",
      columns = c(FRS_a_est, CVR_a_est)
    ) |>
    tab_spanner(
      label = "b-path (HVR \u2192 Cog slope)",
      columns = c(FRS_b_est, CVR_b_est)
    ) |>
    tab_spanner(
      label = "Indirect effect [95% CI]",
      columns = c(FRS_ind, CVR_ind)
    ) |>
    cols_label(
      FRS_a_est = "FRS", CVR_a_est = "CVR",
      FRS_b_est = "FRS", CVR_b_est = "CVR",
      FRS_ind = "FRS", CVR_ind = "CVR"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = FRS_ind,
        rows = FRS_ind_sig == TRUE
      )
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = CVR_ind,
        rows = CVR_ind_sig == TRUE
      )
    )

  if (!is.null(title)) {
    tbl <- tbl |> tab_header(title = title)
  }

  tbl |> style_manuscript_table()
}


# ---- Demographics Table (ADNI vs UKB Normative) ----

#' Create demographics comparison table
#'
#' Compares ADNI analysis cohort (overall and by sex) with
#' UKB normative sample. ADNI stats computed from cohort
#' data; UKB stats from normative_transfer_validation.rds.
#'
#' @param cohort.dt data.table with longitudinal cohort data
#' @param normative_stats list from normative validation RDS
#' @return gt table object
create_demographics_table <- function(cohort.dt,
                                      title,
                                      normative_stats = NULL,
                                      ukb_raw = NULL) {
  bl.dt <- cohort.dt[, .SD[1], by = PTID]
  males.dt <- bl.dt[SEX == "Male"]
  females.dt <- bl.dt[SEX == "Female"]
  n_m <- nrow(males.dt)
  n_f <- nrow(females.dt)

  ms <- function(m, s, d = 1) {
    sprintf(
      paste0("%.", d, "f (%.", d, "f)"), m, s
    )
  }

  grp <- function(dt, col, d = 1) {
    ms(mean(dt[[col]], na.rm = TRUE),
       sd(dt[[col]], na.rm = TRUE), d)
  }

  # --- UKB sex-stratified helpers ---------------
  ukb_n <- function(sex) {
    st <- normative_stats[[sex]]
    if (is.null(st)) return("\u2014")
    format(st$n, big.mark = ",")
  }
  ukb_ms <- function(sex, m_f, s_f, d = 1) {
    st <- normative_stats[[sex]]
    if (is.null(st)) return("\u2014")
    ms(st[[m_f]], st[[s_f]], d)
  }
  ukb_vol <- function(sex, roi, d = 2) {
    st <- ukb_raw[[sex]]
    m_f <- paste0(roi, "_mean")
    s_f <- paste0(roi, "_sd")
    if (is.null(st) ||
        is.null(st[[m_f]])) return("\u2014")
    ms(st[[m_f]], st[[s_f]], d)
  }

  dash <- "\u2014"

  dx_levels <- c("CU", "MCI", "AD")
  dx_row <- function(dx, dt, n) {
    ct <- sum(dt$DX == dx, na.rm = TRUE)
    sprintf("%d (%.1f)", ct, 100 * ct / n)
  }

  rows.lst <- list(
    data.table(
      Variable = "N",
      UKB_Males = ukb_n("male"),
      UKB_Females = ukb_n("female"),
      ADNI_Males = as.character(n_m),
      ADNI_Females = as.character(n_f)
    ),
    data.table(
      Variable = "Age, years",
      UKB_Males = ukb_ms(
        "male", "age_mean", "age_sd"
      ),
      UKB_Females = ukb_ms(
        "female", "age_mean", "age_sd"
      ),
      ADNI_Males = grp(males.dt, "AGE"),
      ADNI_Females = grp(
        females.dt, "AGE"
      )
    ),
    data.table(
      Variable = "Education, years",
      UKB_Males = ukb_ms(
        "male", "educ_mean", "educ_sd"
      ),
      UKB_Females = ukb_ms(
        "female", "educ_mean", "educ_sd"
      ),
      ADNI_Males = grp(males.dt, "EDUC"),
      ADNI_Females = grp(
        females.dt, "EDUC"
      )
    ),
    {
      hc_col <- if ("HC_vol" %in% names(bl.dt))
        "HC_vol" else "HC_raw"
      data.table(
        Variable = "HC bilateral, cm\u00B3",
        UKB_Males = ukb_vol("male", "hc"),
        UKB_Females = ukb_vol(
          "female", "hc"
        ),
        ADNI_Males = grp(
          males.dt, hc_col, 2
        ),
        ADNI_Females = grp(
          females.dt, hc_col, 2
        )
      )
    },
    {
      hvr_col <- if ("HVR_vol" %in% names(bl.dt))
        "HVR_vol" else "HVR_raw"
      data.table(
        Variable = "HVR (ratio)",
        UKB_Males = ukb_vol(
          "male", "hvr", 3
        ),
        UKB_Females = ukb_vol(
          "female", "hvr", 3
        ),
        ADNI_Males = grp(
          males.dt, hvr_col, 3
        ),
        ADNI_Females = grp(
          females.dt, hvr_col, 3
        )
      )
    },
    data.table(
      Variable = "HC *z*-score",
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = grp(
        males.dt, "HC_z", 2
      ),
      ADNI_Females = grp(
        females.dt, "HC_z", 2
      )
    ),
    data.table(
      Variable = "HVR *z*-score",
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = grp(
        males.dt, "HVR_z", 2
      ),
      ADNI_Females = grp(
        females.dt, "HVR_z", 2
      )
    ),
    data.table(
      Variable = "FRS, %",
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = grp(
        males.dt, "FRS_pct"
      ),
      ADNI_Females = grp(
        females.dt, "FRS_pct"
      )
    ),
    data.table(
      Variable = "APOE \u03b54+, %",
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = sprintf(
        "%.1f",
        100 * mean(
          males.dt$APOE4 == TRUE,
          na.rm = TRUE
        )
      ),
      ADNI_Females = sprintf(
        "%.1f",
        100 * mean(
          females.dt$APOE4 == TRUE,
          na.rm = TRUE
        )
      )
    ),
    data.table(
      Variable = "Diagnosis, *n* (%)",
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = "",
      ADNI_Females = ""
    )
  )

  for (dx in dx_levels) {
    rows.lst <- c(rows.lst, list(data.table(
      Variable = paste0("\n", dx),
      UKB_Males = dash,
      UKB_Females = dash,
      ADNI_Males = dx_row(
        dx, males.dt, n_m
      ),
      ADNI_Females = dx_row(
        dx, females.dt, n_f
      )
    )))
  }

  tbl.dt <- rbindlist(rows.lst)

  data_cols.v <- c(
    "UKB_Males", "UKB_Females",
    "ADNI_Males", "ADNI_Females"
  )
  tbl <- gt(tbl.dt) |>
    fmt_markdown(columns = "Variable") |>
    cols_label(
      Variable = "Characteristic",
      UKB_Males = "Males",
      UKB_Females = "Females",
      ADNI_Males = "Males",
      ADNI_Females = "Females"
    ) |>
    tab_spanner(
      label = "UKB Normative",
      columns = c(
        "UKB_Males", "UKB_Females"
      )
    ) |>
    tab_spanner(
      label = "ADNI Analysis Cohort",
      columns = c(
        "ADNI_Males", "ADNI_Females"
      )
    ) |>
    cols_align(
      align = "left",
      columns = "Variable"
    ) |>
    cols_align(
      align = "center",
      columns = all_of(data_cols.v)
    ) |>
    tab_header(title = title) |>
    style_manuscript_table()
}

# =========================================================
# New Output Tables (Phase 1a)
# =========================================================

#' MIMIC factor loadings table
#'
#' @param mimic_model.res cvr_mimic_model.rds contents
#' @param title Table title
#' @return gt table
create_mimic_loadings_table <- function(
    mimic_model.res, title) {
  if (is.null(mimic_model.res)) return(NULL)
  ld.dt <- as.data.table(
    mimic_model.res$loadings
  )
  if (nrow(ld.dt) == 0) return(NULL)

  INDICATOR_LABELS <- c(
    SBP_z = "Systolic BP, z-scored",
    HTN = "Hypertension, binary",
    GLUCOSE_z = "Glucose, z-scored",
    CHOL_z = "Total Cholesterol, z-scored",
    CREAT_z = "Creatinine, z-scored"
  )

  ind_names <- ld.dt$indicator
  tbl.dt <- data.table(
    Indicator = fifelse(
      ind_names %in% names(INDICATOR_LABELS),
      INDICATOR_LABELS[ind_names],
      ind_names
    ),
    Lambda  = ld.dt$lambda,
    Std     = ld.dt$std_all,
    SE      = ld.dt$se,
    p       = ld.dt$pvalue
  )
  # lavaan reports p = 0 when < machine eps
  tbl.dt[
    !is.na(p) & p < .Machine$double.eps,
    p := .Machine$double.eps
  ]

  gt(tbl.dt) |>
    fmt_number(
      columns = c("Lambda", "Std", "SE"),
      decimals = 3
    ) |>
    sub_missing(
      columns = "p",
      missing_text = "\u2014"
    ) |>
    fmt_number(
      columns = "p",
      rows = !is.na(p) & p >= 0.001,
      decimals = 3
    ) |>
    fmt_scientific(
      columns = "p",
      rows = !is.na(p) & p < 0.001,
      decimals = 2
    ) |>
    cols_label(
      Indicator = "Indicator",
      Lambda = "\u03bb",
      Std = "Std. \u03bb",
      SE = "SE",
      p = "p"
    ) |>
    tab_header(title = title) |>
    style_manuscript_table()
}

#' Measurement invariance table
#'
#' @param mi.res cvr_measurement_invariance.rds contents
#' @param title Table title
#' @return gt table
create_mi_table <- function(mi.res, title) {
  if (is.null(mi.res)) return(NULL)

  levels.v <- c(
    "configural", "metric", "scalar", "strict"
  )

  # --- Fit section ---
  fit_rows.lst <- list()
  for (lvl in levels.v) {
    inv <- mi.res$results[[lvl]]
    if (is.null(inv)) next
    fit_rows.lst[[lvl]] <- data.table(
      Level    = tools::toTitleCase(lvl),
      Converge = inv$converged,
      ChiSq    = inv$chisq %||% NA_real_,
      df       = inv$df %||% NA_integer_,
      CFI      = inv$CFI %||% NA_real_,
      TLI      = inv$TLI %||% NA_real_,
      RMSEA    = inv$RMSEA %||% NA_real_,
      SRMR     = inv$SRMR %||% NA_real_
    )
  }
  if (length(fit_rows.lst) == 0) return(NULL)
  fit.dt <- rbindlist(fit_rows.lst, fill = TRUE)

  # --- Comparison section ---
  cmp <- mi.res$comparisons_delta
  cmp_rows.lst <- list()
  if (!is.null(cmp) && is.list(cmp) &&
      length(cmp) > 0) {
    for (nm in names(cmp)) {
      entry <- cmp[[nm]]
      if (is.null(entry)) next
      constrained_lvl <- entry$test
      if (is.null(constrained_lvl)) {
        constrained_lvl <- sub(
          ".*vs_?", "", nm
        )
      }
      cmp_rows.lst[[nm]] <- data.table(
        Level = tools::toTitleCase(
          constrained_lvl
        ),
        dCFI = entry$d_cfi %||%
          NA_real_,
        dRMSEA = entry$d_rmsea %||%
          NA_real_,
        Decision = entry$decision %||%
          "\u2014"
      )
    }
  }

  # Merge fit + comparison
  if (length(cmp_rows.lst) > 0) {
    cmp.dt <- rbindlist(
      cmp_rows.lst, fill = TRUE
    )
    merged.dt <- merge(
      fit.dt, cmp.dt,
      by = "Level", all.x = TRUE
    )
    merged.dt[
      is.na(Decision), Decision := "\u2014"
    ]
  } else {
    merged.dt <- fit.dt
    merged.dt[, `:=`(
      dCFI = NA_real_,
      dRMSEA = NA_real_,
      Decision = "\u2014"
    )]
  }

  # Reorder
  lvl_order.v <- c(
    "Configural", "Metric", "Scalar", "Strict"
  )
  merged.dt[, Level := factor(
    Level, levels = lvl_order.v
  )]
  setorder(merged.dt, Level)

  gt(merged.dt) |>
    sub_missing(missing_text = "\u2014") |>
    fmt_number(
      columns = c("ChiSq", "CFI", "TLI"),
      decimals = 3
    ) |>
    fmt_number(
      columns = c("RMSEA", "SRMR"),
      decimals = 3
    ) |>
    fmt_number(
      columns = any_of(c("dCFI", "dRMSEA")),
      decimals = 3
    ) |>
    cols_hide(columns = any_of("Converge")) |>
    cols_label(
      Level = "Level",
      ChiSq = "\u03c7\u00b2",
      df = "df",
      CFI = "CFI",
      TLI = "TLI",
      RMSEA = "RMSEA",
      SRMR = "SRMR",
      dCFI = "\u0394CFI",
      dRMSEA = "\u0394RMSEA",
      Decision = "Decision"
    ) |>
    tab_spanner(
      label = "Model Fit",
      columns = c(
        "ChiSq", "df", "CFI", "TLI",
        "RMSEA", "SRMR"
      )
    ) |>
    tab_spanner(
      label = "Comparison",
      columns = any_of(
        c("dCFI", "dRMSEA", "Decision")
      )
    ) |>
    tab_header(title = title) |>
    style_manuscript_table()
}

#' Attrition logistic regression table
#'
#' @param attrition.res attrition_analysis.rds contents
#' @param title Table title
#' @return gt table
create_attrition_logistic_table <- function(
    attrition.res, title) {
  if (is.null(attrition.res) ||
      is.null(attrition.res$logistic_model)) {
    return(NULL)
  }
  lm_info <- attrition.res$logistic_model
  or.dt <- as.data.table(lm_info$odds_ratios)

  # Normalize CI column names
  ci_lo_col <- intersect(
    c("CI_lower", "2.5 %", "ci_lower"),
    names(or.dt)
  )[1]
  ci_hi_col <- intersect(
    c("CI_upper", "97.5 %", "ci_upper"),
    names(or.dt)
  )[1]

  if (!is.na(ci_lo_col)) {
    setnames(or.dt, ci_lo_col, "CI_lower")
  }
  if (!is.na(ci_hi_col)) {
    setnames(or.dt, ci_hi_col, "CI_upper")
  }

  PRED_LABELS <- c(
    SEXMale = "Sex (Male)",
    AGE = "Age",
    FRS = "FRS",
    HVR_z = "HVR *z*-score",
    HC_z = "HC *z*-score",
    APOE4TRUE = "APOE \u03b54+"
  )
  # Preferred row order
  PRED_ORDER <- names(PRED_LABELS)

  if ("term" %in% names(or.dt)) {
    or.dt[, Predictor := fifelse(
      term %in% names(PRED_LABELS),
      PRED_LABELS[term], term
    )]
  } else {
    or.dt[, Predictor := paste0("V", .I)]
  }

  # Remove intercept
  or.dt <- or.dt[
    !grepl("ntercept", term, ignore.case = TRUE)
  ]

  # Compute p-value from CI (approximate)
  or.dt[, p_value := 2 * pnorm(
    -abs(log(OR)) /
      ((log(CI_upper) - log(CI_lower)) /
        (2 * 1.96))
  )]

  # Reorder rows: Sex first, Age second
  or.dt[, order_idx := match(
    term, PRED_ORDER, nomatch = 999L
  )]
  setorder(or.dt, order_idx)

  fmt_p.fn <- function(x) {
    ifelse(
      x < 0.001, "< 0.001",
      sprintf("%.3f", x)
    )
  }

  display.dt <- or.dt[, .(
    Predictor,
    OR = sprintf("%.2f", OR),
    CI = sprintf(
      "[%.2f, %.2f]", CI_lower, CI_upper
    ),
    p_value
  )]

  gt(display.dt) |>
    fmt_markdown(columns = "Predictor") |>
    fmt(columns = "p_value", fns = fmt_p.fn) |>
    cols_label(
      Predictor = "Predictor",
      OR = "OR",
      CI = "95% CI",
      p_value = "p"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = p_value < 0.05
      )
    ) |>
    tab_header(title = title) |>
    style_manuscript_table()
}

#' EDT early/middle sensitivity table
#'
#' @param edt_sens.res edt_early_middle_sensitivity.rds
#' @param title Table title
#' @return gt table
create_edt_sensitivity_table <- function(
    edt_sens.res, title) {
  if (is.null(edt_sens.res) ||
      is.null(edt_sens.res$results_table)) {
    return(NULL)
  }
  rt.dt <- copy(edt_sens.res$results_table)

  DOMAIN_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )
  if ("Domain" %in% names(rt.dt)) {
    rt.dt[, Domain := fifelse(
      Domain %in% names(DOMAIN_LABELS),
      DOMAIN_LABELS[Domain], Domain
    )]
  }

  BRAIN_LABELS <- c(
    HVR_z = "HVR *z*",
    HC_z = "HC *z*"
  )
  if ("Brain_Measure" %in% names(rt.dt)) {
    rt.dt[, Brain_Measure := fifelse(
      Brain_Measure %in% names(BRAIN_LABELS),
      BRAIN_LABELS[Brain_Measure],
      Brain_Measure
    )]
  }

  disp_cols.v <- intersect(
    c(
      "Brain_Measure", "Domain",
      "Full_Beta", "Full_p",
      "Restricted_Beta", "Restricted_p",
      "Pct_Change", "Consistent"
    ),
    names(rt.dt)
  )
  display.dt <- rt.dt[, ..disp_cols.v]

  tbl <- gt(display.dt) |>
    sub_missing(missing_text = "\u2014") |>
    fmt_markdown(
      columns = any_of("Brain_Measure")
    )

  num_cols.v <- intersect(
    c(
      "Full_Beta", "Restricted_Beta",
      "Pct_Change"
    ),
    names(display.dt)
  )
  if (length(num_cols.v) > 0) {
    tbl <- tbl |>
      fmt_number(
        columns = all_of(num_cols.v),
        decimals = 4
      )
  }

  p_cols.v <- intersect(
    c("Full_p", "Restricted_p"),
    names(display.dt)
  )
  if (length(p_cols.v) > 0) {
    tbl <- tbl |>
      fmt_scientific(
        columns = all_of(p_cols.v),
        decimals = 2
      )
  }

  label_map.lst <- list(
    Brain_Measure = "Brain",
    Domain = "Domain",
    Full_Beta = "\u03b2 (Full)",
    Full_p = "*p* (Full)",
    Restricted_Beta = "\u03b2 (E+M)",
    Restricted_p = "*p* (E+M)",
    Pct_Change = "% Change",
    Consistent = "Consistent"
  )
  label_map.lst <- label_map.lst[
    names(label_map.lst) %in% disp_cols.v
  ]
  tbl <- tbl |>
    cols_label(.list = label_map.lst) |>
    fmt_markdown(
      columns = any_of(
        c("Full_p", "Restricted_p")
      )
    )

  tbl |>
    tab_header(title = title) |>
    style_manuscript_table()
}

#' LGCM parallel process growth table
#'
#' @param lgcm_parallel.res lgcm_parallel_results.rds
#' @param title Table title
#' @return gt table
create_lgcm_growth_table <- function(
    lgcm_parallel.res, title) {
  if (is.null(lgcm_parallel.res) ||
      is.null(lgcm_parallel.res$linear_results)) {
    return(NULL)
  }

  DOMAIN_LABELS <- c(
    MEM = "Memory",
    LAN = "Language",
    EXF = "Executive"
  )

  rows.lst <- list()
  for (key in names(
    lgcm_parallel.res$linear_results
  )) {
    res <- lgcm_parallel.res$linear_results[[key]]
    pe.dt <- as.data.table(res$params)

    get_est.fn <- function(lhs, op, rhs) {
      row <- pe.dt[
        label == paste0(lhs, op, rhs) |
          (name == paste0(lhs, op, rhs))
      ]
      if (nrow(row) == 0) {
        row <- pe.dt[
          grepl(lhs, name) &
            grepl(rhs, name) &
            grepl(op, name, fixed = TRUE)
        ]
      }
      if (nrow(row) > 0) row$Estimate[1]
      else NA_real_
    }

    # Extract means and variances from params
    cog_i_mean <- pe.dt[
      grepl("cog_i", name) &
        grepl("mean", name, ignore.case = TRUE),
      Estimate
    ]
    cog_s_mean <- pe.dt[
      grepl("cog_s", name) &
        grepl("mean", name, ignore.case = TRUE),
      Estimate
    ]
    cog_i_var <- pe.dt[
      grepl("cog_i", name) &
        grepl("var", name, ignore.case = TRUE),
      Estimate
    ]
    cog_s_var <- pe.dt[
      grepl("cog_s", name) &
        grepl("var", name, ignore.case = TRUE),
      Estimate
    ]

    parts <- strsplit(key, "_")[[1]]
    sex_lbl <- parts[1]
    dom_lbl <- if (
      parts[2] %in% names(DOMAIN_LABELS)
    ) DOMAIN_LABELS[parts[2]] else parts[2]

    rows.lst[[key]] <- data.table(
      Sex = sex_lbl,
      Domain = dom_lbl,
      N = res$n,
      Int_Mean = if (length(cog_i_mean) > 0)
        cog_i_mean[1] else NA_real_,
      Int_Var = if (length(cog_i_var) > 0)
        cog_i_var[1] else NA_real_,
      Slope_Mean = if (length(cog_s_mean) > 0)
        cog_s_mean[1] else NA_real_,
      Slope_Var = if (length(cog_s_var) > 0)
        cog_s_var[1] else NA_real_
    )
  }

  if (length(rows.lst) == 0) return(NULL)
  growth.dt <- rbindlist(rows.lst)

  # Sex row groups with N
  growth.dt[, Sex_group := fifelse(
    !is.na(N),
    sprintf(
      "%ss (N = %s)", Sex,
      format(N, big.mark = ",")
    ),
    paste0(Sex, "s")
  )]

  gt(
    growth.dt,
    rowname_col = "Domain",
    groupname_col = "Sex_group"
  ) |>
    fmt_number(
      columns = c(
        "Int_Mean", "Int_Var",
        "Slope_Mean", "Slope_Var"
      ),
      decimals = 3
    ) |>
    cols_hide(columns = any_of(c("Sex", "N"))) |>
    cols_label(
      Int_Mean = "Mean",
      Int_Var = "Variance",
      Slope_Mean = "Mean",
      Slope_Var = "Variance"
    ) |>
    tab_spanner(
      label = "Intercept",
      columns = c("Int_Mean", "Int_Var")
    ) |>
    tab_spanner(
      label = "Slope",
      columns = c("Slope_Mean", "Slope_Var")
    ) |>
    tab_header(title = title) |>
    tab_source_note(
      "Growth parameters are CVR-independent."
    ) |>
    style_manuscript_table()
}

#' Simulation comparison table
#'
#' @param simulation.res simulation_results.rds
#' @param title Table title
#' @return gt table
create_simulation_table <- function(
    simulation.res, title) {
  if (is.null(simulation.res) ||
      is.null(simulation.res$summary)) {
    return(NULL)
  }
  sim.dt <- as.data.table(simulation.res$summary)

  # Build display columns
  disp_cols.v <- intersect(
    c(
      "n",
      "n_converged_se", "n_converged_free",
      "bias_se", "rmse_se", "coverage_se",
      "bias_free", "rmse_free", "coverage_free"
    ),
    names(sim.dt)
  )

  display.dt <- sim.dt[, ..disp_cols.v]

  tbl <- gt(display.dt) |>
    sub_missing(missing_text = "\u2014")

  bias_cols.v <- grep(
    "^bias", names(display.dt), value = TRUE
  )
  rmse_cols.v <- grep(
    "^rmse", names(display.dt), value = TRUE
  )
  cov_cols.v <- grep(
    "^coverage", names(display.dt),
    value = TRUE
  )

  if (length(bias_cols.v) > 0) {
    tbl <- tbl |> fmt_number(
      columns = all_of(bias_cols.v),
      decimals = 4
    )
  }
  if (length(rmse_cols.v) > 0) {
    tbl <- tbl |> fmt_number(
      columns = all_of(rmse_cols.v),
      decimals = 4
    )
  }
  if (length(cov_cols.v) > 0) {
    tbl <- tbl |> fmt_number(
      columns = all_of(cov_cols.v),
      decimals = 1
    )
  }

  label_map.lst <- list(n = "N")
  se_labels.lst <- list(
    n_converged_se = "Converged",
    bias_se = "Bias",
    rmse_se = "RMSE",
    coverage_se = "Coverage (%)"
  )
  free_labels.lst <- list(
    n_converged_free = "Converged",
    bias_free = "Bias",
    rmse_free = "RMSE",
    coverage_free = "Coverage (%)"
  )
  all_labels.lst <- c(
    label_map.lst,
    se_labels.lst[
      names(se_labels.lst) %in% disp_cols.v
    ],
    free_labels.lst[
      names(free_labels.lst) %in% disp_cols.v
    ]
  )
  tbl <- tbl |>
    cols_label(.list = all_labels.lst)

  se_cols.v <- intersect(
    c(
      "n_converged_se",
      "bias_se", "rmse_se", "coverage_se"
    ),
    disp_cols.v
  )
  free_cols.v <- intersect(
    c(
      "n_converged_free",
      "bias_free", "rmse_free", "coverage_free"
    ),
    disp_cols.v
  )
  if (length(se_cols.v) > 0) {
    tbl <- tbl |> tab_spanner(
      label = "SE-Constrained",
      columns = all_of(se_cols.v)
    )
  }
  if (length(free_cols.v) > 0) {
    tbl <- tbl |> tab_spanner(
      label = "Free Loading",
      columns = all_of(free_cols.v)
    )
  }

  tbl |>
    tab_header(title = title) |>
    style_manuscript_table()
}

#' FRS-Age correlation table
#'
#' @param frs_problem.res frs_problem_analysis.rds
#' @param title Table title
#' @return gt table
create_frs_correlation_table <- function(
    frs_problem.res, title) {
  if (is.null(frs_problem.res) ||
      is.null(frs_problem.res$age_correlation)) {
    return(NULL)
  }
  ac <- frs_problem.res$age_correlation

  rows.lst <- list()

  # Overall Pearson
  if (!is.null(ac$pearson)) {
    rows.lst[["overall_p"]] <- data.table(
      Group = "Overall",
      Method = "Pearson",
      r = ac$pearson
    )
  }

  # Overall Spearman
  if (!is.null(ac$spearman)) {
    rows.lst[["overall_s"]] <- data.table(
      Group = "Overall",
      Method = "Spearman",
      r = ac$spearman
    )
  }

  # By sex
  if (!is.null(ac$by_sex)) {
    for (sx in names(ac$by_sex)) {
      sx_info <- ac$by_sex[[sx]]
      if (!is.null(sx_info$pearson)) {
        rows.lst[[paste0(sx, "_p")]] <-
          data.table(
            Group = sx,
            Method = "Pearson",
            r = sx_info$pearson
          )
      }
    }
  }

  # Partial
  if (!is.null(ac$partial)) {
    rows.lst[["partial"]] <- data.table(
      Group = "Partial",
      Method = "Partial *r*",
      r = ac$partial
    )
  }

  if (length(rows.lst) == 0) return(NULL)
  corr.dt <- rbindlist(rows.lst)

  gt(corr.dt,
     groupname_col = "Group") |>
    fmt_number(columns = "r", decimals = 3) |>
    fmt_markdown(columns = "Method") |>
    cols_label(
      Method = "Method",
      r = "r"
    ) |>
    tab_header(title = title) |>
    style_manuscript_table()
}

# -----------------------------------------------------------
# Manuscript inline helpers
# -----------------------------------------------------------

#' Format p-value for inline reporting
#'
#' @param p  Numeric p-value
#' @return Character string, e.g. "p < 0.001"
#' @export
fmt_p_inline <- function(p) {
  if (p < 0.001) "p < 0.001"
  else sprintf("p = %.3f", p)
}

#' Check whether bootstrap CI excludes zero
#'
#' @param r  List with boot_ci_lower, boot_ci_upper
#' @return TRUE if CI excludes zero
#' @export
ci_sig_check.fn <- function(r) {
  !is.null(r$boot_ci_lower) &&
    !is.null(r$boot_ci_upper) &&
    !is.na(r$boot_ci_lower) &&
    !is.na(r$boot_ci_upper) &&
    (r$boot_ci_lower > 0 ||
       r$boot_ci_upper < 0)
}
