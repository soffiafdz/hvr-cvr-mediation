# ============================================
# consort_diagram.R
# CONSORT Flow Diagram generation (DiagrammeR)
# ============================================

#' Generate CONSORT flow diagram DOT specification
#'
#' Two-track diagram:
#' Left:  ADSP-PHC -> A+ -> MRI+Cog -> LME -> LGCM
#' Right: A-neg -> MRI -> CU -> Age-comparable
#' Left exclusions branch LEFT; right exclusions RIGHT
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
#' @param ukb_age_range  Character label, e.g. "46\u201382"
#' @return DOT string
generate_consort_dot <- function(
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

  fmt <- function(x) format(x, big.mark = ",")

  # Derived counts
  n_no_amy <- n_database - n_with_amyloid
  n_aneg <- n_with_amyloid - n_amyloid_pos
  exc_mri_l <- n_amyloid_pos - n_with_mri
  exc_temp <- n_with_mri - n_lme

  # --- Style constants ---
  wh <- ', fillcolor="white"'
  gr <- paste0(
    ', fillcolor="#E8F4E8"',
    ', color="#2E7D32", penwidth=2'
  )
  ex <- paste0(
    ', fillcolor="#F5F5F5"',
    ', color="#999999"',
    ', width=2.2, height=0.6',
    ', fontsize=10'
  )
  dash_l <- paste0(
    ' [style=dashed, color="#999999"',
    ', constraint=false]'
  )
  dash_r <- dash_l
  invis <- " [style=invis]"

  # ---- Nodes ----
  nodes <- sprintf(paste0(
    '  db [label="ADSP-PHC Database\\n',
    'N = %s"%s]\n',
    '  with_amy [label=',
    '"With Amyloid Data\\n',
    'N = %s"%s]\n',
    '  exc_amy [label=',
    '"Excluded:\\n',
    'No amyloid data\\n',
    'n = %s"%s]\n',
    '  apos [label=',
    '"Amyloid-Positive\\n',
    '(PET >25 CL or CSF A+)\\n',
    'N = %s"%s]\n',
    '  aneg [label=',
    '"Amyloid-Negative\\n',
    'N = %s"%s]\n',
    '  exc_mri_l [label=',
    '"Excluded:\\n',
    'No structural MRI\\n',
    'n = %d"%s]\n',
    '  mri_l [label=',
    '"Structural MRI +\\n',
    'Cognitive Data\\n',
    'N = %d"%s]\n',
    '  exc_temp [label=',
    '"Excluded:\\n',
    'MRI after cognition\\n',
    'n = %d"%s]\n',
    '  lme [label=',
    '"LME Analysis Cohort\\n',
    'N = %d"%s]\n'
  ),
    fmt(n_database), wh,
    fmt(n_with_amyloid), wh,
    fmt(n_no_amy), ex,
    fmt(n_amyloid_pos), wh,
    fmt(n_aneg), wh,
    exc_mri_l, ex,
    n_with_mri, wh,
    exc_temp, ex,
    n_lme, gr
  )

  # ---- Edges (left track) ----
  edges <- paste0(
    '  db -> with_amy\n',
    '  db -> exc_amy', dash_r, '\n',
    '  with_amy -> apos\n',
    '  with_amy -> aneg\n',
    '  apos -> mri_l\n',
    '  apos -> exc_mri_l', dash_l, '\n',
    '  mri_l -> lme\n',
    '  mri_l -> exc_temp', dash_l, '\n'
  )

  # ---- Ranks ----
  ranks <- paste0(
    '  {rank=same; with_amy; exc_amy}\n',
    '  {rank=same; apos; aneg}\n',
    '  {rank=same; exc_mri_l; mri_l}\n',
    '  {rank=same; exc_temp; lme}\n'
  )

  # ---- LGCM (optional) ----
  lgcm_block <- ""
  if (!is.null(n_lgcm)) {
    exc_lgcm <- n_lme - n_lgcm
    lgcm_block <- sprintf(paste0(
      '  exc_lgcm [label=',
      '"Excluded:\\n',
      '< 2 visits or no EDT\\n',
      'n = %d"%s]\n',
      '  lgcm [label=',
      '"LGCM Subsample\\n',
      '(\u2265 2 visits + EDT)\\n',
      'N = %d"%s]\n',
      '  lme -> lgcm\n',
      '  lme -> exc_lgcm', dash_l, '\n',
      '  {rank=same; exc_lgcm; lgcm}\n'
    ),
      exc_lgcm, ex,
      n_lgcm, gr
    )
  }

  # ---- Right track (optional) ----
  right_block <- ""
  has_right <- !is.null(n_aneg_mri) &&
    !is.null(n_aneg_cu_mri) &&
    !is.null(n_age_comparable)

  if (has_right) {
    exc_mri_r <- n_aneg - n_aneg_mri
    exc_dx <- n_aneg_mri - n_aneg_cu_mri
    exc_age <- n_aneg_cu_mri - n_age_comparable
    age_lbl <- if (!is.null(ukb_age_range)) {
      sprintf("Age %s", ukb_age_range)
    } else {
      "Age-Comparable"
    }

    right_block <- sprintf(paste0(
      '  mri_r [label=',
      '"With Structural MRI\\n',
      'N = %d"%s]\n',
      '  exc_mri_r [label=',
      '"Excluded:\\n',
      'No structural MRI\\n',
      'n = %d"%s]\n',
      '  cu [label=',
      '"Cognitively Unimpaired\\n',
      'N = %d"%s]\n',
      '  exc_dx [label=',
      '"Excluded:\\n',
      'MCI or AD\\n',
      'n = %d"%s]\n',
      '  age_comp [label=',
      '"%s\\n',
      'Subsample\\n',
      'N = %d"%s]\n',
      '  exc_age [label=',
      '"Excluded:\\n',
      'Outside UKB age range\\n',
      'n = %d"%s]\n',
      '  aneg -> mri_r\n',
      '  aneg -> exc_mri_r', dash_r, '\n',
      '  mri_r -> cu\n',
      '  mri_r -> exc_dx', dash_r, '\n',
      '  cu -> age_comp\n',
      '  cu -> exc_age', dash_r, '\n',
      '  {rank=same; mri_l; mri_r; exc_mri_r}\n',
      '  {rank=same; lme; cu; exc_dx}\n'
    ),
      n_aneg_mri, wh,
      exc_mri_r, ex,
      n_aneg_cu_mri, wh,
      exc_dx, ex,
      age_lbl,
      n_age_comparable, gr,
      exc_age, ex
    )

    # Full L-R chains for each rank row
    # Row: exc_mri_l | mri_l | mri_r | exc_mri_r
    # Row: exc_temp  | lme   | cu    | exc_dx
    right_block <- paste0(
      right_block,
      '  exc_mri_l -> mri_l', invis, '\n',
      '  mri_l -> mri_r', invis, '\n',
      '  mri_r -> exc_mri_r', invis, '\n',
      '  exc_temp -> lme', invis, '\n',
      '  lme -> cu', invis, '\n',
      '  cu -> exc_dx', invis, '\n'
    )

    # LGCM + age-comparable same rank
    if (!is.null(n_lgcm)) {
      right_block <- paste0(
        right_block,
        '  {rank=same; lgcm; age_comp; exc_age}\n',
        '  exc_lgcm -> lgcm', invis, '\n',
        '  lgcm -> age_comp', invis, '\n',
        '  age_comp -> exc_age', invis, '\n'
      )
    } else {
      right_block <- paste0(
        right_block,
        '  {rank=same; age_comp; exc_age}\n',
        '  age_comp -> exc_age', invis, '\n'
      )
    }
  }

  # ---- Assemble ----
  dot <- paste0(
    'digraph CONSORT {\n',
    '  graph [rankdir=TB',
    ', splines=polyline',
    ', nodesep=0.5',
    ', ranksep=0.7]\n',
    '  node [shape=box',
    ', style="filled"',
    ', fontname="Helvetica"',
    ', fontsize=11',
    ', width=2.5',
    ', height=0.8]\n',
    '  edge [arrowsize=0.8]\n\n',
    nodes, '\n',
    edges, '\n',
    ranks, '\n',
    lgcm_block,
    right_block,
    '}\n'
  )

  return(dot)
}

#' Render CONSORT flow diagram
#'
#' @param ... Arguments passed to generate_consort_dot()
#' @return DiagrammeR grViz object
render_consort_diagram <- function(...) {
  if (!requireNamespace(
    "DiagrammeR", quietly = TRUE
  )) {
    stop(
      "DiagrammeR package required",
      " for CONSORT diagram"
    )
  }
  dot <- generate_consort_dot(...)
  DiagrammeR::grViz(dot)
}
