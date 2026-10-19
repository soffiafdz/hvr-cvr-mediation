# R analysis project conventions

Reference implementation: `/ipl/ipl27/sfernandez/hvr_cvr_mediation` (ADNI and
ADSP-PHC data; linear mixed-effects and OpenMx latent growth curve models; Quarto
manuscript for *Alzheimer's & Dementia*). Start every new R analysis project from
this setup. Where a new project needs to differ, make the change deliberately and
write it down in that project's README.

Sections 1 and 12 are rules. The rest describes how the reference project is built
and what to copy.

---

## 1. Non-negotiables

- **Never fabricate.** No numbers, citations, file contents or interpretations from
  memory. Every number in a manuscript, report or grant traces to an RDS file
  produced by pipeline code. Never hardcode a result in prose or code.
- **Citations:** verify against Crossref, the DOI or the publisher, or ask. Do not
  "recall" a reference.
- **Verify before stating.** Compute a number from the files before quoting it;
  read a function before describing its signature.
- **If uncertain, ask.**

---

## 2. Toolchain

| Layer | What | Details in the reference project |
|---|---|---|
| System environment | mamba/conda env, exported to `environment.yml` | Env `ukb_adni_lgcm`: R 4.5.1 (`r-base`), Quarto 1.7.33, pandoc, compilers, `radian`, `direnv`, plus Python 3.10/Jupyter tooling. Exported with full build pins. |
| Activation | direnv, `.envrc` (gitignored) | `source <mamba>/etc/profile.d/conda.sh`, `eval "$(direnv stdlib)"`, `conda activate ukb_adni_lgcm` |
| R packages | renv: `renv.lock`, `renv/`, `.Rprofile` containing `source("renv/activate.R")` | Repository: Posit Package Manager `cran/latest`. `renv/settings.json`: `snapshot.type: implicit`, `use.cache: false`. Local packages (e.g. the ADNIMERGE tarball) are recorded as `Source: Local`. renv 1.1.5. |
| Paths | `here` plus an empty `.here` file at the project root | All paths are built with `here()` |
| TeX | system `xelatex` (`/usr/bin/xelatex`) | PDF rendering; font TeX Gyre Termes for Times |
| Documents | Quarto (`reports-src/_quarto.yml`) | Render to **PDF only** |

The conda env also installs many R packages (`r-essentials`, `r-tidyverse`), but the
renv library takes precedence once `.Rprofile` activates renv.

Setup (from the README):

```bash
conda env create -f environment.yml
conda activate <env>          # or: direnv allow
Rscript -e 'renv::restore()'
```

Known gaps in the reference project. Fix these in new projects rather than copying
them:
- Every R run prints "The project is out-of-sync": the library does not match
  `renv.lock`. Keep `renv::status()` clean.
- `lintr` is not installed, so the `.lintr` rules are checked by hand
  (`awk 'length>80' file` for line length). Add `lintr` to the lock file.

---

## 3. Repository layout

```
project/
├── .here                     # empty; anchors here()
├── .Rprofile                 # source("renv/activate.R")
├── .lintr                    # naming regex (§6)
├── .envrc                    # (gitignored) direnv: activate the conda env
├── environment.yml           # conda env export
├── renv.lock, renv/
├── README.md                 # overview, data access, setup, pipeline order,
│                             # project rules (copy §1, §6, §8, §11)
├── config/pipeline_config.yaml
├── R/
│   ├── scripts/              # NN_verb_object.R, run in numeric order
│   └── utils/                # sourced helper modules (§5)
├── reports-src/              # .qmd sources, _quarto.yml, references.bib, CSL,
│                             # title.tex, build_docx.py
├── docs/                     # audits, check scripts, conventions (tracked)
├── data/{raw,external,derivatives}/     # gitignored, .gitkeep tracked
├── models/{fits,results,diagnostics}/   # gitignored
├── outputs/{figures,tables,reports}/    # gitignored
│   └── manuscript_env.rds    # cached environment for the manuscript (§8)
├── logs/                     # gitignored; one log file per R process
├── archive/<date>_<label>/   # gitignored frozen snapshots + INDEX.md (§9)
└── submission/               # gitignored assembled journal package
```

`.gitignore` policy (copy it):
- `data/`, `models/`, `outputs/`, `logs/` ignored, each with `!<dir>/.gitkeep`.
- `archive/` and `submission/`.
- R session files, R CMD build/check output, `.Rproj.user/`.
- `renv/library/`, `renv/staging/`, `renv/local/`.
- `_freeze/`, `*_cache/`, `/.quarto/`; `*.tex` except `reports-src/title.tex`.
- `.envrc`, `.Renviron`.
- Root-level `/*.pdf`, `/*.png`, `/*.svg` for posters, screenshots and other
  personal material.

`reports-src/.gitignore`: `/.quarto/`, `*.html`, `*_files/`.

---

## 4. Configuration: one YAML file

Everything a script needs to know lives in `config/pipeline_config.yaml`. Scripts
contain no paths and no analysis parameters.

| Section | Holds | Read with |
|---|---|---|
| `general` | `project_name`, `random_seed` (1618), `log_level`, `n_cores` | `get_config("general", ...)`, `set_seed()` |
| `paths` | directory layout | `get_config("paths", ...)` |
| `data` | every input and output file, grouped `raw`, `external`, `derivatives`, `models`, `outputs` | `get_data_path(group, key)` returns an absolute `here()` path |
| `parameters` | analysis parameters: domains and labels, variable names, thresholds, specification switches | `get_parameter(...)` |
| method sections (`lme`, `lgcm`) | estimator and optimizer settings | `get_config("lgcm", ...)` |
| `output` | table and figure formats, DPI, sizes, journal limits | `get_config("output", ...)` |
| `scripts.force_regenerate` | one flag per script, commented with the script name | `get_script_setting("force_regenerate", key, default = FALSE)` |
| `packages` | required and optional packages | `validate_packages()` |

Rules:
- **Dated raw file names live only in the config**, e.g.
  `adsp_phc_cvrf: data/raw/adsp_phc/ADSP_PHC_CVRF_11Dec2025.csv`.
- **One loader per derived variable.** If several scripts need the same variable,
  write one function (e.g. `load_frs.fn()` in `R/utils/frs.R`) instead of repeating
  `fread(select = ...)`. In the reference project the Framingham score was read in
  7 places, several with the column name hardcoded, and all of them had to be
  found and changed.
- **Alternative specifications are config switches**, not code edits, so earlier
  results stay reproducible. Example: `parameters.frs.source: corrected | published`.
- **Every script with a cache check has its key** under `scripts.force_regenerate`.
- **To re-run a script:** set its flag to `true`, run it, set it back to `false`.
  Never delete outputs to force a re-run; manuscripts and theses depend on the
  files existing.

---

## 5. Helper library (`R/utils/`)

Plain files loaded with `source(here("R/utils/<file>.R"))`, not a package. Each
module sources its own dependencies (e.g. `data_io.R` sources `logging.R`).

**`config.R`** (copy verbatim)
- `load_config(config_file = here("config/pipeline_config.yaml"))`: reads the YAML
  into the global `.config`.
- `get_config(..., default = NULL)`: nested lookup; stops on a missing key unless a
  default is given.
- `get_data_path(...)`: `here()` + `data:` lookup.
- `get_parameter(..., default = NULL)`, `get_script_setting(..., default = NULL)`.
- `get_seed()`, `set_seed()`.
- `validate_config(config = NULL)`: requires the sections `general`, `data`,
  `parameters` and `scripts`, `general.random_seed`, and `data.raw`,
  `data.derivatives` and `data.models`.
- Project-specific: `get_valid_timepoints()`, `get_n_timepoints()`,
  `get_factor_levels(variable_name)` (canonical factor levels, e.g.
  `SEX = c("Female", "Male")`).

**`logging.R`** (copy verbatim)
- Initialises itself on the first message and writes to
  `logs/pipeline_YYYYMMDD_HHMMSS.log`, one file per R process. Messages at INFO and
  above are also printed to the console.
- `log_debug()`, `log_info()`, `log_warn()`, `log_error()` take `(message, ...)`
  in `sprintf` style.
- `log_script_start(name)`, `log_script_end(name, success = TRUE)`,
  `log_section(name)` (prints `>>> name`), `log_time(description, expr)`.
- Warnings printed by packages (e.g. OpenMx status codes) go to the console only,
  not to the log file. Capture stdout and stderr when running long jobs.

**`data_io.R`** (copy verbatim)
- `check_files_exist(files.v, stop_on_missing = TRUE)`.
- `read_rds_safe(file_path, description = NULL)`.
- `write_rds_safe(object, file_path, description = NULL)`: creates the directory and
  logs `Saved: <path>`.
- `read_csv_safe(file_path, ..., description = NULL)`: wraps `fread`.
- `needs_regeneration(output_path, input_paths.v = NULL, force_regenerate = FALSE)`:
  compares modification times.
- `ensure_directory(dir_path)`, `path_join(...)`.
- `require_output(x, label)`: stops if `x` is `NULL`. Use it for anything a
  reported number depends on.

**`validation.R`** (copy verbatim)
- `validate_columns(data.dt, required_cols.v, data_name)`,
  `validate_not_empty(x, data_name)`.
- `validate_num_range(data.dt, col_name, min_val, max_val, na_allowed = TRUE)`,
  `validate_categorical(data.dt, col_name, expected_levels.v, allow_extra = FALSE)`.
- `validate_age(...)`, `validate_unique_ids(data.dt, id_cols.v)`,
  `validate_packages(packages.v)`.

**`tables.R`** (gt; take the generic part)
- `style_manuscript_table(gt_tbl, table_type = "main")`.
- `save_table(gt_table, filename, output_dir, formats = c("html", "tex"),
  fix_latex = TRUE, create_standalone = TRUE)`, `save_table_dual(...)`.
- Formatters: `format_p`, `fmt_p`, `format_pval`, `format_ci`, `format_effect_ci`,
  `sig_stars`.
- One `create_<thing>_table()` per manuscript table.
- Manuscript helpers: `ci_sig_check.fn()` (significant if the bootstrap CI excludes
  0), `n_of_total.fn()` ("4 of 6").

**`plotting.R`** (ggplot2; take the generic part)
- `theme_publication(base_size = 10, use_markdown = FALSE)` (based on
  `theme_classic`).
- `manuscript_colors()`: a named list of Paul Tol colourblind-safe colours (e.g.
  `female`, `male`, one colour per predictor and measure). Every figure takes its
  colours from it.
- `get_palette(type)`.
- `save_plot(plot, filename, width = 7, height = 7, dpi = 600)`.
- One `plot_<thing>()` per figure.

**Domain modules:** one file per method or construct, e.g. `disease_time.R`
(estimated disease time via progmod), `lavaan_helpers.R`, `consort_diagram.R`
(TikZ flow diagram) and `frs.R` (Framingham recomputation with a built-in check).

---

## 6. Code style

- **Naming** (enforced by `.lintr`): `snake_case` plus a type suffix.
  - `.dt` data.table, `.df` data.frame, `.v` vector, `.lst` list, `.fn` function,
    `.fit` fitted model, `.mod` model, `.res` results object, `.tbl` table,
    `.plt` plot, `.mat` matrix, `.fct` factor, `.path` path, `.subdt` subset.
  - UPPERCASE for constants (`COGNITIVE_DOMAINS`, `N_BOOT`).
  - The `.lintr` regex:
    `^([a-z][a-z0-9_]*(\.df|\.dt|\.fct|\.fit|\.fn|\.lst|\.mat|\.mod|\.path|\.plt|\.subdt|\.res|\.tbl|\.v)?|[A-Z][A-Z0-9_]*)$`.
  - Older helpers predate the `.fn` suffix. New code uses it.
- **Lines of 80 characters or fewer.**
- **data.table throughout.**
  - Write `:=` against each table explicitly. Don't loop over `list(a.dt, b.dt)`
    to add columns by reference.
  - Never name a function argument the same as a column used inside `[`:
    `coef.dt[match(sex, coef.dt$sex)]` silently uses the column.
- **Roxygen-style `#'` headers** (params, return) on every helper; `# ---` section
  banners in scripts.
- **Idempotent scripts.** A script that reads a file and writes it back must give
  the same result when run twice. Derive values from the raw columns every time.
- **Fail loudly.** Use `stopifnot()` or `require_output()` for anything a reported
  number depends on. Never fall back silently (`next`, an empty string) in a
  manuscript chunk.

---

## 7. Pipeline scripts

- Named `NN_verb_object.R`; the number is the execution order. Add-ons take a
  letter suffix (`10b_`, `13b_`, `16b_`).
- Run from the project root: `Rscript R/scripts/NN_verb_object.R`.
- The header block documents purpose, model specification, rationale, verified
  references, `INPUTS` and `OUTPUTS`. Keep the file name in the header in sync when
  renumbering (one reference script still calls itself by its old number).
- Results are saved as one named list per script. Include the results, sample
  sizes, methodology notes, the relevant config values and a timestamp.

Template (from `08_lme_hvr_z.R`):

```r
#!/usr/bin/env Rscript
# =============================================================================
# NN_verb_object.R - One-line purpose
# =============================================================================
# What and why; model specification; references (verified).
#
# INPUTS:
#   - data/derivatives/...
# OUTPUTS:
#   - models/results/...
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

log_script_start("NN_verb_object.R")
config <- load_config()
validate_config(config)
set_seed()

# --- Cache check ---
FORCE_REGENERATE <- get_script_setting(
  "force_regenerate", "nn_key", default = FALSE
)
output.path <- get_data_path("models", "nn_output")
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output exists and force_regenerate=FALSE; skipping")
  log_script_end("NN_verb_object.R", success = TRUE)
  quit(status = 0)
}

# --- Configuration ---
DOMAINS <- get_parameter("cognitive_domains")

log_section("Part 1: Load data")
# ...
log_section("Part N: Save results")
write_rds_safe(results.lst, output.path, "NN results")
log_script_end("NN_verb_object.R", success = TRUE)
```

---

## 8. From results to manuscript

- **One environment script.** The last pipeline script
  (`15_prepare_manuscript_env.R`) loads every result file, computes every number
  the text uses and saves a single named list to `outputs/manuscript_env.rds`
  (95 objects in the reference project).
- **Manuscript setup chunk.** It sets the working directory to the project root
  (found by checking for `.here`), sources the utils, calls `load_config()`, then
  runs `list2env(readRDS(here("outputs/manuscript_env.rds")))`.
- **Prose** uses inline R only (`` `r sprintf("%.1f", pct_white)` ``). A new
  number goes into script 15, not into a chunk.
- **After any upstream change**, re-run script 15 before rendering.
- **`reports-src/_quarto.yml`:** `output-dir: ../outputs/reports`,
  `execute-dir: project`, `echo`, `warning` and `message` all false, plus an
  explicit render list.
- **Render to PDF only:**
  `cd reports-src && quarto render <file>.qmd --to pdf`.
- **Word copies for journals:** `python3 reports-src/build_docx.py` converts the
  rendered `.tex`. It handles the TikZ diagram, gt tables and equations.
- **gt to LaTeX:** a custom `knit_print.gt_tbl` in the setup chunk fixes table
  widths and source notes.
- **Captions:** every `tbl-cap` and `fig-cap` starts with a bold title of 3 to 8
  words ending in a full stop, then the details:
  `"**Short title.** Details..."`. The bold part becomes the short caption in lists
  of tables and figures.
- **References:** `references.bib` (Zotero export, abbreviated journal names) and
  a CSL file in `reports-src/`. Verify every entry against Crossref.
- **Prose style:** no em-dashes (use semicolons, colons or a new sentence). Avoid
  stock openers such as "Of note", "It is worth noting" and "Several limitations
  warrant note".

---

## 9. Provenance and snapshots

- **Seeds** come from the config (`set_seed()` at the top of every script).
  Every run writes a log.
- **Snapshot before any re-run that changes results.**
  - Copy `outputs/`, `models/results/` and `data/derivatives/` to
    `archive/<YYYY-MM-DD>_<label>/`.
  - Check the copy with `diff -rq`.
  - Write an `INDEX.md`: purpose, contents, recovery commands, related docs and
    branch.
  - The latest snapshot is the "before" state for any diff.
- **Investigations** go in `docs/<topic>_audit.md` with a companion
  `docs/<topic>_check.R` that reproduces every number in the document when run
  from the project root. Name them as audits and checks.
- **Derived variables from external releases:** recompute them from their inputs
  when possible, and assert agreement with the released version. In the reference
  project `frs.R` stops unless the opposite-sex equation reproduces the released
  column exactly. That check is how a sex swap in the ADSP-PHC Framingham score
  was confirmed.

---

## 10. Long runs, bootstraps and OpenMx

- **Checkpoints don't follow `force_regenerate`.** Bootstrap checkpoints
  (`models/results/<method>/<name>_checkpoints/<key>.rds`) are reused whatever the
  flags say. Their validity check looks only at stored flags (`bootstrap_n`,
  `shared_seed`, `fixed_quad`, `boot_all_ci`). When data or specification changes:
  - either add a flag for it to the validity check;
  - or remove the stale checkpoints after copying them into the snapshot.
- **Runtime reference:** an LGCM mediation model with 1,000 bootstraps takes about
  1.5 hours on 24 cores. Script 13 took 9 to 19 hours, depending on how many
  checkpoints it could reuse.
- **Detached runs.** Run long pipelines with
  `nohup setsid bash runner.sh > run.log 2>&1 &`. Have the runner write
  `START`/`END`/`FAILED`/`ALL_DONE` lines per script to a status file, and stop at
  the first failure. Check `ps` and the logs, not the terminal.
- **OpenMx: exogenous predictors.** When predictors are manifest variables:
  - never fix their means to 0 and variances to 1;
  - fix means, variances and covariances at the sample values (lavaan's
    `fixed.x`), or estimate them;
  - centre every predictor so the latent means describe an average participant.

  With raw age fixed at N(0, 1), −2LL carried a constant of about 3 million, the
  Hessian was not convex, SEs were missing and nested LRTs failed. See
  `docs/lgcm_exogenous_predictors_audit.md`.
- **OpenMx: multiple starts.** `mxTryHard` stops at the first acceptable solution,
  which can be a local optimum. Fit from several starts (default values and the
  estimates of a nested model) and keep the lowest −2LL. In the reference project
  a single start missed the optimum by 54.5 −2LL units.
- **OpenMx: check every fit.** Look at the status code (5 means "Hessian not
  convex"), missing SEs and whether nested fits succeeded. To compare two
  specifications of the same model, subtract the predictors' own log-likelihood
  and compare the conditional −2LL.
- **OpenMx: factor scores.** `mxFactorScores` levels can move with the way
  predictor means are parameterised. Draw population mean curves from
  model-implied means, not from averaged factor scores. (This is still under
  investigation in the reference project.)

---

## 11. Git

- One-line commit messages. No AI attribution or co-author lines.
- Stage files by name; never `git add -A`. One logical change per commit.
  Interactive `git add -p` is unavailable to agents; to split one file's changes,
  stage an intermediate version with
  `git hash-object -w` + `git update-index --cacheinfo`.
- Never `git checkout` or switch branches unless asked. To change another branch,
  use a temporary `git worktree`.
- Branches: `fix/<topic>` for corrections. Findings with consequences outside the
  repository (a submitted paper, a preprint, a data provider) stay on an unpushed
  local branch until their disclosure is decided.
- Commit and push only when asked.
- No "handoff", "Claude" or similar in file names or document text that others
  will see.

---

## 12. Working agreements

- **A question is not a request to edit.** Answer it, propose the change, wait for
  a go-ahead.
- **Confirm first** before anything outward-facing (push, publish) or
  destructive.
- **Deletions.** Agent auto mode blocks `rm` on outputs. Give the exact command and
  let the user run it (prefix `!` in Claude Code).
- **Report faithfully,** with numbers, including failures and checks that did not
  pass. Correct earlier statements explicitly when new evidence contradicts them.
- **Results feed other work.** Grant proposals (`/ipl/ipl27/sfernandez/grants`,
  with `FACTS.md`), posters and preprints reuse these results. Flag any claim
  there that a change affects.

---

## 13. Lessons from the reference project

| Problem | What prevents it |
|---|---|
| A released derived variable was wrong: ADSP-PHC's Framingham score used the opposite sex's equation | Recompute from the inputs; assert agreement with the release |
| The same variable was read in 7 places, several with the column name hardcoded | One loader function, with the column name in the config |
| A script overwrote its own input and wasn't safe to run twice (`FRS_raw` became a z-score) | Idempotent scripts; derive from raw columns every run |
| Bootstrap checkpoints outlived a data change | Data and specification flags in the checkpoint validity check |
| Predictor moments fixed at 0 and 1 in OpenMx | Sample moments plus centring (§10) |
| A single-start fit hit a local optimum | Multiple starts; compare −2LL |
| A manuscript chunk silently printed nothing for months | `stopifnot()` on every element a chunk reads |
| A "calibrated" constant was hardcoded in a simulation script (`slope_mean = -0.19`) | Load it from the results file it came from |
| The config said 62, the code hardcoded 62, the paper said 74 (Framingham age range) | One source of truth (the config), citing the source |
| A unit bug: the ceiling was checked at 0.30 on a percent scale | `validate_num_range()` on inputs |
| Script headers kept old numbers after renumbering | Update headers when renaming |

---

## 14. Starting a new project

1. Create the conda env (adapt `environment.yml`) and `.envrc`; run `renv::init()`;
   put `lintr` in the lock.
2. Copy `R/utils/config.R`, `logging.R`, `data_io.R` and `validation.R` unchanged.
   Start `tables.R` and `plotting.R` from the theme, palette, formatting and save
   helpers only.
3. Write `config/pipeline_config.yaml` with the sections `general`, `paths`,
   `data`, `parameters`, `output`, `scripts.force_regenerate` and `packages`.
4. Add `.here`, `.Rprofile`, `.lintr` and `.gitignore` (copy them), plus `data/`,
   `models/`, `outputs/` and `logs/` with `.gitkeep`.
5. Put §1, §6, §8 (captions and prose) and §11 of this document in the README.
6. Build the numbered scripts from the template in §7. The last one writes
   `outputs/<document>_env.rds`, which the `.qmd` loads.
7. Set up `reports-src/` with `_quarto.yml`, `references.bib`, a CSL file and
   `title.tex`. Render to PDF.
