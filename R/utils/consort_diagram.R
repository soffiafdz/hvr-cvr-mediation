# ============================================
# consort_diagram.R
# CONSORT Flow Diagram generation (TikZ)
# ============================================

#' Generate CONSORT flow diagram as TikZ code
#'
#' Two-track diagram:
#' Left:  ADSP-PHC -> A+ -> MRI+Cog -> LME -> LGCM
#' Right: A-neg -> MRI -> CU -> Age-comparable
#'
#' @param n_database    Total in ADSP-PHC
#' @param n_with_amyloid Subjects with amyloid data
#' @param n_amyloid_pos  Amyloid-positive
#' @param n_with_mri     A+ with MRI + cognition
#' @param n_lme          LME analysis cohort
#' @param n_lgcm         LGCM subsample (optional)
#' @param n_aneg_mri     A-neg with MRI (optional)
#' @param n_aneg_cu_mri  A-neg CU with MRI (optional)
#' @param n_age_comparable Age-comparable (optional)
#' @param ukb_age_range  Character label, e.g. "46--82"
#' @return TikZ string for inclusion in LaTeX
generate_consort_tikz <- function(
    n_database,
    n_with_amyloid,
    n_amyloid_pos,
    n_with_mri,
    n_lme,
    n_lgcm = NULL,
    n_aneg_mri = NULL,
    n_aneg_cu_mri = NULL,
    n_age_comparable = NULL,
    ukb_age_range = NULL) {

  fmt <- function(x) {
    format(x, big.mark = ",")
  }

  # Derived counts
  n_no_amy <- n_database - n_with_amyloid
  n_aneg <- n_with_amyloid - n_amyloid_pos
  exc_mri_l <- n_amyloid_pos - n_with_mri
  exc_temp <- n_with_mri - n_lme

  has_lgcm <- !is.null(n_lgcm)
  has_right <- !is.null(n_aneg_mri) &&
    !is.null(n_aneg_cu_mri) &&
    !is.null(n_age_comparable)

  if (has_lgcm) {
    exc_lgcm <- n_lme - n_lgcm
  }
  if (has_right) {
    exc_mri_r <- n_aneg - n_aneg_mri
    exc_dx <- n_aneg_mri - n_aneg_cu_mri
    exc_age <- n_aneg_cu_mri - n_age_comparable
    age_lbl <- if (!is.null(ukb_age_range)) {
      paste0("Age ", ukb_age_range)
    } else {
      "Age-Comparable"
    }
  }

  # Helper: build a TikZ node line
  node.fn <- function(style, id, x, y, label) {
    paste0(
      "  \\node[", style, "] (", id,
      ") at (", x, ",", y, ")\n",
      "    {", label, "};\n"
    )
  }

  # ---- Layout geometry ----
  # Full text width = 16.5cm (1in margins)
  # Process boxes: 3.6cm text + ~0.2cm padding = ~3.8cm
  # Excluded boxes: 2.6cm text + ~0.15cm padding = ~2.75cm
  # Total: 2.75 + 0.5 + 3.8 + 0.8 + 3.8 + 0.5 + 2.75
  #      = 14.9cm (centered in 16.5cm)
  lx <- -2.5   # left track center
  rx <- 2.5    # right track center
  cx <- 0      # center (db, with_amy)
  exl <- -3.8  # left exclusion offset from track
  exr <- 3.8   # right exclusion offset from track

  # Row y-positions
  y_db <- 0
  y_amy <- -2.0
  y_split <- -4.5
  y_mri <- -6.5
  y_lme <- -8.5
  y_lgcm <- -10.5

  # No-hyphenation directive for nodes
  nohyph <- paste0(
    "execute at begin node={",
    "\\hyphenpenalty=10000 ",
    "\\exhyphenpenalty=10000}"
  )

  # ---- Preamble ----
  lines.v <- c(
    "\\centering",
    "\\begin{tikzpicture}[",
    paste0(
      "  every node/.style=",
      "{font=\\footnotesize},"
    ),
    "  process/.style={",
    "    rectangle, draw=black!70, thick,",
    paste0(
      "    text width=3.6cm, inner sep=3pt,"
    ),
    paste0(
      "    minimum height=0.8cm, ",
      nohyph, ","
    ),
    paste0(
      "    align=center, ",
      "rounded corners=2pt},"
    ),
    "  final/.style={",
    paste0(
      "    rectangle, draw=black!70, ",
      "very thick,"
    ),
    paste0(
      "    text width=3.6cm, inner sep=3pt,"
    ),
    paste0(
      "    minimum height=0.8cm, ",
      nohyph, ","
    ),
    paste0(
      "    align=center, ",
      "rounded corners=2pt,"
    ),
    "    fill=black!5},",
    "  excluded/.style={",
    "    rectangle, draw=black!50, thick,",
    paste0(
      "    text width=2.6cm, inner sep=2pt,"
    ),
    paste0(
      "    minimum height=0.5cm, ",
      nohyph, ","
    ),
    paste0(
      "    align=center, ",
      "rounded corners=2pt,"
    ),
    paste0(
      "    fill=black!8, text=black!70, ",
      "font=\\scriptsize},"
    ),
    "  arrow/.style={",
    paste0(
      "    -{Stealth[length=5pt]}, ",
      "thick, black!70},"
    ),
    "  dasharrow/.style={",
    paste0(
      "    -{Stealth[length=4pt]}, ",
      "dashed, black!50}"
    ),
    "]"
  )

  tikz <- paste0(
    paste(lines.v, collapse = "\n"), "\n"
  )

  # ---- Row 0: Database ----
  tikz <- paste0(tikz, node.fn(
    "process", "db", cx, y_db,
    paste0(
      "ADSP-PHC Database\\\\[2pt]",
      "\\textbf{N\\,=\\,", fmt(n_database), "}"
    )
  ))

  # ---- Row 1: With amyloid + exclusion ----
  tikz <- paste0(tikz, node.fn(
    "process", "with_amy", cx, y_amy,
    paste0(
      "With Amyloid Data\\\\[2pt]",
      "\\textbf{N\\,=\\,",
      fmt(n_with_amyloid), "}"
    )
  ))
  exc_amy_x <- 4.5
  tikz <- paste0(tikz, node.fn(
    "excluded", "exc_amy", exc_amy_x, y_amy,
    paste0(
      "Excluded:\\\\No amyloid data\\\\",
      "n\\,=\\,", fmt(n_no_amy)
    )
  ))

  # ---- Row 2: A+ / A- split ----
  tikz <- paste0(tikz, node.fn(
    "process", "apos", lx, y_split,
    paste0(
      "Amyloid-Positive\\\\",
      "(PET >25 CL or CSF A+)\\\\[2pt]",
      "\\textbf{N\\,=\\,",
      fmt(n_amyloid_pos), "}"
    )
  ))
  tikz <- paste0(tikz, node.fn(
    "process", "aneg", rx, y_split,
    paste0(
      "Amyloid-Negative\\\\[2pt]",
      "\\textbf{N\\,=\\,", fmt(n_aneg), "}"
    )
  ))

  # ---- Row 3: MRI + exclusions ----
  tikz <- paste0(tikz, node.fn(
    "process", "mri_l", lx, y_mri,
    paste0(
      "Structural MRI +\\\\",
      "Cognitive Data\\\\[2pt]",
      "\\textbf{N\\,=\\,",
      fmt(n_with_mri), "}"
    )
  ))
  tikz <- paste0(tikz, node.fn(
    "excluded", "exc_mri_l",
    lx + exl, y_mri,
    paste0(
      "Excluded:\\\\No structural MRI\\\\",
      "n\\,=\\,", fmt(exc_mri_l)
    )
  ))
  if (has_right) {
    tikz <- paste0(tikz, node.fn(
      "process", "mri_r", rx, y_mri,
      paste0(
        "With Structural MRI\\\\[2pt]",
        "\\textbf{N\\,=\\,",
        fmt(n_aneg_mri), "}"
      )
    ))
    tikz <- paste0(tikz, node.fn(
      "excluded", "exc_mri_r",
      rx + exr, y_mri,
      paste0(
        "Excluded:\\\\No structural MRI\\\\",
        "n\\,=\\,", fmt(exc_mri_r)
      )
    ))
  }

  # ---- Row 4: LME / CU + exclusions ----
  tikz <- paste0(tikz, node.fn(
    "final", "lme", lx, y_lme,
    paste0(
      "\\textbf{LME Analysis Cohort}",
      "\\\\[2pt]",
      "\\textbf{N\\,=\\,", fmt(n_lme), "}"
    )
  ))
  tikz <- paste0(tikz, node.fn(
    "excluded", "exc_temp",
    lx + exl, y_lme,
    paste0(
      "Excluded:\\\\MRI after cognition\\\\",
      "n\\,=\\,", fmt(exc_temp)
    )
  ))
  if (has_right) {
    tikz <- paste0(tikz, node.fn(
      "process", "cu", rx, y_lme,
      paste0(
        "Cognitively Unimpaired\\\\[2pt]",
        "\\textbf{N\\,=\\,",
        fmt(n_aneg_cu_mri), "}"
      )
    ))
    tikz <- paste0(tikz, node.fn(
      "excluded", "exc_dx",
      rx + exr, y_lme,
      paste0(
        "Excluded:\\\\MCI or AD\\\\",
        "n\\,=\\,", fmt(exc_dx)
      )
    ))
  }

  # ---- Row 5: LGCM / Age-comparable ----
  if (has_lgcm) {
    tikz <- paste0(tikz, node.fn(
      "final", "lgcm", lx, y_lgcm,
      paste0(
        "\\textbf{LGCM Subsample}\\\\",
        "($\\geq$ 2 visits + EDT)",
        "\\\\[2pt]",
        "\\textbf{N\\,=\\,",
        fmt(n_lgcm), "}"
      )
    ))
    tikz <- paste0(tikz, node.fn(
      "excluded", "exc_lgcm",
      lx + exl, y_lgcm,
      paste0(
        "Excluded:\\\\",
        "< 2 visits or no EDT\\\\",
        "n\\,=\\,", fmt(exc_lgcm)
      )
    ))
  }
  if (has_right) {
    tikz <- paste0(tikz, node.fn(
      "final", "age_comp", rx, y_lgcm,
      paste0(
        "\\textbf{", age_lbl,
        " Subsample}\\\\[2pt]",
        "\\textbf{N\\,=\\,",
        fmt(n_age_comparable), "}"
      )
    ))
    tikz <- paste0(tikz, node.fn(
      "excluded", "exc_age",
      rx + exr, y_lgcm,
      paste0(
        "Excluded:\\\\",
        "Outside UKB age range\\\\",
        "n\\,=\\,", fmt(exc_age)
      )
    ))
  }

  # ---- Arrows ----
  # Vertical main flow
  arrows.v <- c(
    "  % Main vertical flow",
    "  \\draw[arrow] (db) -- (with_amy);",
    "  \\draw[arrow] (apos) -- (mri_l);",
    "  \\draw[arrow] (mri_l) -- (lme);"
  )

  # T-junction: with_amy splits to apos / aneg
  # Vertical down, then right-angle to each
  mid_y <- (y_amy + y_split) / 2
  arrows.v <- c(
    arrows.v,
    "  % T-junction: amyloid split",
    paste0(
      "  \\path (with_amy.south)",
      " -- ++(0,", mid_y - y_amy, ")",
      " coordinate (split);"
    ),
    paste0(
      "  \\draw[thick, black!70]",
      " (with_amy.south) -- (split);"
    ),
    "  \\draw[arrow] (split) -| (apos.north);",
    "  \\draw[arrow] (split) -| (aneg.north);"
  )

  # Exclusion arrows (all horizontal)
  arrows.v <- c(
    arrows.v,
    "  % Exclusion arrows",
    paste0(
      "  \\draw[dasharrow]",
      " (with_amy.east) -- (exc_amy.west);"
    ),
    paste0(
      "  \\draw[dasharrow]",
      " (mri_l.west) -- (exc_mri_l.east);"
    ),
    paste0(
      "  \\draw[dasharrow]",
      " (lme.west) -- (exc_temp.east);"
    )
  )

  if (has_lgcm) {
    arrows.v <- c(
      arrows.v,
      "  \\draw[arrow] (lme) -- (lgcm);",
      paste0(
        "  \\draw[dasharrow]",
        " (lgcm.west) -- (exc_lgcm.east);"
      )
    )
  }

  if (has_right) {
    arrows.v <- c(
      arrows.v,
      "  % Right-track vertical flow",
      "  \\draw[arrow] (aneg) -- (mri_r);",
      "  \\draw[arrow] (mri_r) -- (cu);",
      "  \\draw[arrow] (cu) -- (age_comp);",
      "  % Right-track exclusions",
      paste0(
        "  \\draw[dasharrow]",
        " (mri_r.east) -- (exc_mri_r.west);"
      ),
      paste0(
        "  \\draw[dasharrow]",
        " (cu.east) -- (exc_dx.west);"
      ),
      paste0(
        "  \\draw[dasharrow]",
        " (age_comp.east) -- (exc_age.west);"
      )
    )
  }

  tikz <- paste0(
    tikz,
    paste(arrows.v, collapse = "\n"), "\n",
    "\\end{tikzpicture}\n"
  )

  return(tikz)
}
