# OpenMx Codebase Audit

Cross-reference of the OpenMx code in `R/scripts/12_lgcm_parallel.R`,
`R/scripts/13_lgcm_mediation.R`, `R/scripts/14_lgcm_simulation.R`, and
`R/utils/openmx_helpers.R` against authoritative OpenMx documentation.

## Critical — likely affects results

### 1. HVR residuals silently constrained equal across timepoints (script 13)

**Files:** `R/scripts/13_lgcm_mediation.R:670` (`build_mediation.fn`) and
`R/scripts/13_lgcm_mediation.R:906` (`build_mimic_mediation.fn`).

```r
hvr_resid <- mxPath(
  from = hvr_vars, arrows = 2,
  free = TRUE, values = 0.5,
  labels = "hvr_resid"    # singular → recycled across all 6 timepoints
)
```

OpenMx treats shared labels as a hard equality constraint. With `from =
hvr_vars` (6 variables) and a single recycled label, **all 6 HVR residual
variances are forced equal across T0–T5** — a strong homoscedasticity
assumption that was not stated in the model description.

Script 12 specifies them as time-specific
(`R/scripts/12_lgcm_parallel.R:269, 717, 956`):

```r
labels = paste0("hvr_resid_", VALID_TIMEPOINTS)
```

This inconsistency means mediation results are estimated under a different
measurement model than the parallel-process results. Fix: use
`paste0("hvr_resid_", VALID_TIMEPOINTS)` in both
`build_mediation.fn` and `build_mimic_mediation.fn`.

### 2. `extrapolate_edt_hybrid()` does not guarantee NA-free EDT columns

**File:** `R/utils/disease_time.R:508-544`

Subjects with **zero** valid EDT observations fall through both branches
(`length(valid_idx.v) >= 2` and `== 1`) and leave NAs in place. OpenMx
documentation is explicit: *"NA value for a definition variable is Not
Yet Implemented"* — the fit errors out, it does not silently drop those
rows. The comment at `13_lgcm_mediation.R:1448-1449`:

```r
# OpenMx FIML handles all missingness —
# no complete.cases() filtering needed.
```

is misleading. FIML handles NAs in **endogenous** variables, not in
definition variables. EDT_T0…EDT_T5 are loaded as `data.<var>` labels
and must be non-NA for every row OpenMx sees.

Fix: add an explicit `stop()` (or warning + filter) after
`extrapolate_edt_hybrid()` when any row remains all-NA.

### 3. Helper `extract_openmx_fit_indices.fn` calls `mxRefModels()` on ITVS models

**File:** `R/utils/openmx_helpers.R:619`

The helper computes CFI/TLI/RMSEA via `mxRefModels(model.fit, run = TRUE)`.
Scripts 12 and 13 correctly avoid this and have an in-source comment
explaining why (`12_lgcm_parallel.R:617-621`):

> "ITVS models (definition variables for EDT) cannot use mxRefModels: it
> ignores individually-varying loadings, making CFI/TLI/RMSEA/chi-sq
> meaningless."

The OpenMx maintainers confirm this on the forums: *"The saturated model
will ignore your definition variables… mxRefModels() is behaving as it
normally does, which is not appropriate for the kind of growth models
in which one is interested."*

The helper is currently not sourced by the pipeline scripts in production
(12/13/14 define their own local `extract_fit_indices.fn`), but it sits
in `R/utils/openmx_helpers.R` and will mislead any future user who
imports it.

Fix: either delete `extract_openmx_fit_indices.fn` (and the other unused
helpers — `create_bivariate_lgcm_openmx`, `fit_openmx_model`,
`compare_lavaan_openmx`, `assess_concordance`), or guard the function
with an explicit ITVS check that errors instead of returning meaningless
indices.

## Medium — inconsistencies and small issues

### 4. Optimizer not pinned in script 13

- `12_lgcm_parallel.R:115`: `mxOption(NULL, "Default optimizer", "SLSQP")`
- `14_lgcm_simulation.R:64`: `mxOption(NULL, "Default optimizer", "SLSQP")`
- `13_lgcm_mediation.R:148`: only sets `Number of Threads`, leaves
  optimizer at OpenMx's default (CSOLNP or NPSOL depending on build)

The mediation results are therefore optimized with a different algorithm
than the parallel-process and simulation results.

Fix: add `mxOption(NULL, "Default optimizer", "SLSQP")` near line 148 of
script 13.

### 5. `Number of Threads = detectCores()` during a sequential bootstrap loop

**File:** `R/scripts/13_lgcm_mediation.R:148-151`

`bootstrap_indirect.fn` (lines 1174+) is a sequential `for` loop, but
each `mxRun` is told to use all cores. Per-fit thread overhead dominates
for small models. Script 12 explicitly caps at 2 threads, which is the
better choice.

### 6. Status-code acceptance is inconsistent across scripts

- Script 13 (mediation): accepts `code %in% c(0, 1)` —
  `13_lgcm_mediation.R:1536`
- Script 12 (linear/fixed-quad main loop): requires `code == 0`
  exactly — `12_lgcm_parallel.R:1248`
- Script 12 (quadratic attempt): accepts `code %in% c(0, 1)` —
  `12_lgcm_parallel.R:442`
- Script 14 (simulation): accepts only `code == 0` —
  `14_lgcm_simulation.R:359, 379`

`mxTryHard()` documents code 1 ("OK with warnings") as acceptable when
`greenOK=TRUE`. The asymmetry means script 12's main fixed-quadratic
analysis may discard fits that 13/14 would keep.

Fix: pick one rule and apply consistently.

## Minor — code-hygiene and documented workarounds

### 7. `1e-10` variance on fixed-quadratic factors

**Files:** `12_lgcm_parallel.R:874-878`, `13_lgcm_mediation.R:796-800,
1098-1102`

Folklore, not a documented OpenMx recommendation. Works because
`mxFactorScores` and the Cholesky/inverse operations need a strictly
positive variance. Alternatives:

- Drop the quadratic latent factor entirely and add a single fixed
  effect on the manifests via a path with a free coefficient on
  `"one"`-style means.
- Leave the quadratic variance free with a small lower bound via
  `lbound`.

Flag in the methods section that this is an ad-hoc epsilon, not
canonical practice.

### 8. Manual case-resampling bootstrap instead of `mxBootstrap()`

**File:** `R/scripts/13_lgcm_mediation.R:1152-1212`

The manual loop is correct in spirit but loses (a) `mxBootstrapEval` /
`mxBootstrapStdizeRAMpaths` machinery for downstream derived quantities,
(b) automatic filtering of non-converged replicates, (c) the parallel
infrastructure built into `mxBootstrap()`. Since `mxAlgebra(a*b, name =
"indirect")` already exists in the model, `mxBootstrap(fit) +
mxBootstrapEval("indirect", boot)` would give the same answer with
less code and the convergence filter baked in.

### 9. Dead / unsynced helper code

`R/utils/openmx_helpers.R` contains `create_bivariate_lgcm_openmx()`,
`fit_openmx_model()`, `compare_lavaan_openmx()`, `assess_concordance()`,
`extract_openmx_fit_indices.fn()` — none are sourced or called by the
pipeline scripts (12/13/14 each define their own local helpers).
Either delete or update them; in their current form they propagate the
`mxRefModels`-on-ITVS error (issue 3).

## What's correct

- Definition-variable syntax (`labels = paste0("data.", edt_cols)` with
  `free = FALSE`) matches the official OpenMx specification.
- `mxAlgebra(a * b, name = "indirect")` correctly resolves to
  free-parameter labels in a RAM model.
- `coef(fit)["a"]` extraction in the bootstrap loop is the canonical
  short form (OpenMx aliases `coef` to `omxGetParameters`).
- The decision to omit CFI/TLI/RMSEA for ITVS models
  (`12_lgcm_parallel.R:617-621` / `13_lgcm_mediation.R:477-498`) is
  correct per OpenMx documentation.

## Sources

- OpenMx Definition Variables documentation:
  <https://openmx.ssri.psu.edu/docs/OpenMx/2.5.1/DefinitionMeans_Path.html>
- `mxBootstrap` CRAN reference:
  <https://search.r-project.org/CRAN/refmans/OpenMx/html/mxBootstrap.html>
- `mxTryHard` reference:
  <https://www.rdocumentation.org/packages/OpenMx/versions/2.21.11/topics/mxTryHard>
- `mxRefModels` / `omxSaturatedModel` reference:
  <https://openmx.ssri.psu.edu/docs/OpenMx/2.6.7/_static/Rdoc/omxSaturatedModel.html>
- OpenMx forum on NAs in definition variables:
  <https://openmx.ssri.psu.edu/forums/opensem-forums/behavioral-genetics-models/using-fiml-impute-missing-covariate-data>
- OpenMx forum on `mxRefModels` with ITVS (Boker / Estabrook):
  <https://openmx.ssri.psu.edu/node/4456>
- `omxGetParameters` wiki:
  <https://openmx.ssri.psu.edu/wiki/omxgetparameters>
- `mxFactorScores` reference:
  <https://rdrr.io/cran/OpenMx/man/mxFactorScores.html>
