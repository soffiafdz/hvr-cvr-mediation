# Framingham Risk Score sex swap in the ADSP-PHC data

Found 1 October 2026 during a review of a grant application, then verified
independently three ways. No pipeline script, config, output or manuscript file
has been changed yet. The correction will be made on the branch `fix/frs-sex-swap`.

---

## 1. Summary

The variable `PHC_ASCVD_10y_FRS_Simple_Ageover30` in the ADSP Phenotype
Harmonization Consortium (ADSP-PHC) release
`data/raw/adsp_phc/ADSP_PHC_CVRF_11Dec2025.csv` was computed with the **opposite
sex's Framingham equation for every participant**. Men were scored with the women's
equation and women with the men's.

The sex labels in the same file (`PHC_Sex`) are correct. The pipeline reads the
variable unmodified. So the error is **upstream, in the ADSP-PHC derivation**, not in
this repository.

Consequence: the manuscript's ceiling result is reversed, and every Framingham-based
result in the manuscript was computed with the wrong-sex score.

| Analysis cohort baseline (n = 923; 513 men, 410 women) | All | Men | Women |
|---|---|---|---|
| % at the 30% cap, PHC variable as published (manuscript) | 48.9 | 28.1 | 74.9 |
| % at the 30% cap, correct-sex equation | **54.8** | **80.1** | **23.2** |
| Mean FRS (%), PHC as published | | 21.1 | 28.1 |
| Mean FRS (%), correct-sex equation | | 28.7 | 20.3 |

The corrected direction is what Framingham is designed to give: men carry more
predicted risk than women of the same age and risk-factor profile.

Within each sex, the published and correct scores correlate only moderately
(Spearman 0.688 in men, 0.745 in women). Sex-stratified Framingham results can
therefore change materially, not just shift in scale.

---

## 2. Where the swap occurred

**In the ADSP-PHC computation of the Framingham score, at the step that passes
participant sex to the risk equation.** The evidence rules out every other place.

| Candidate location | Verdict | Evidence |
|---|---|---|
| This repository's pipeline | **Not here** | `config/pipeline_config.yaml` (`frs: column:`) names the PHC variable. Scripts 03, 07, 10b and 15 read it directly from the raw CSV, with no recoding. The pipeline's `SEX` matches ADNIMERGE `PTGENDER` (table in §3). Applying the published variable to the pipeline's own sex labels reproduces the manuscript numbers exactly (48.9 / 28.1 / 74.9). |
| `PHC_Sex` in the ADSP-PHC file | **Correct** | Against ADNIMERGE `PTGENDER`: code 1 is 1,269 men and 2 women; code 2 is 1,154 women and 0 men. The data dictionary also says 1 = Male, 2 = Female. |
| The CVrisk R package | **Correct** | CRAN source `R/ascvd_10y_frs_simple.R` and coefficient table `data/frs_simple_coef.rda` label male and female correctly: male `ln_age` 3.11296, baseline survival 0.88431; female `ln_age` 2.72107, baseline survival 0.94833. The function merges coefficients on the `gender` string and re-sorts by input order, so it does not scramble rows. |
| **ADSP-PHC's derivation** | **Here** | The stored value equals CVrisk's arithmetic (same rounding to 2 decimals, same clamp to 1 to 30) applied with the **other** sex's coefficients, in 3,485 of 3,485 rows, with a maximum difference of 0.000. |

**Narrowing it further.** The ADSP-PHC methods document
(`data/raw/adsp_phc/ADNI_ADSP_PHC_CVRF_METHODS_20250217.pdf`, Step B1) says the
score was computed with `CVrisk::ascvd_10y_frs_simple` from "the participant's
gender (male or female)". Sex came from `PTDEMOG_23Feb2024.csv`, where ADNI codes
`PTGENDER` 1 = Male, 2 = Female.

However, they did **not** run stock CVrisk unchanged. Stock CVrisk sets age above 74
to `NA` and returns `NA`, yet 1,887 rows aged over 74 have a value; the oldest is
98.9. The variable name's suffix `_Ageover30` points the same way. So ADSP-PHC used a
modified copy or a re-implementation of the function.

The swap therefore sits either in their mapping from `PTGENDER` (or a harmonized sex
code) to the strings `"male"`/`"female"`, or in a coefficient table inside their
modified copy. **Their code is not available here, so the exact line cannot be
identified.**

The release history in the methods PDF names the intermediate file as
`ADNI_ phc_VRF_2024Oct_n3501_frs_pc1_fixed.csv` (downloaded 1 November 2024, V2,
3,500 participants). Whether V1 (Dec 2023, 1,577 participants) had the same error is
unknown.

---

## 3. Evidence and how to reproduce

Run from the repo root:

```bash
tar -xzf data/external/ADNIMERGE_0.0.1.tar.gz -C /tmp ADNIMERGE/data/adnimerge.rdata
Rscript docs/frs_sex_swap_check.R /tmp/ADNIMERGE/data/adnimerge.rdata
```

Output on 1 October 2026:

```
Rows with PHC FRS and complete inputs: 3485
Male (PHC_Sex=1, n=1796): max|PHC - own-sex eq| = 20.850; max|PHC - other-sex eq| = 0.000; share matching other-sex eq (<=0.01) = 1.000
Female (PHC_Sex=2, n=1689): max|PHC - own-sex eq| = 19.880; max|PHC - other-sex eq| = 0.000; share matching other-sex eq (<=0.01) = 1.000

PHC_Sex vs ADNIMERGE PTGENDER:
       PTGENDER
PHC_Sex Female Male
      1      2 1269
      2   1154    0

Pipeline SEX vs ADNIMERGE PTGENDER (cohort baseline):
        PTGENDER
SEX      Female Male
  Female    410    0
  Male        2  511

Analysis cohort baseline: n=923 (men 513, women 410)
  % at cap, PHC as published: all 48.9 | men 28.1 | women 74.9
  % at cap, correct-sex eq:   all 54.8 | men 80.1 | women 23.2
  Mean FRS, PHC as published: men 21.1, women 28.1
  Mean FRS, correct-sex eq:   men 28.7, women 20.3
  Spearman r(PHC, correct) within men 0.688, within women 0.745

Rows aged > 74 with a PHC value (stock CVrisk returns NA): 1887
```

The same two participants disagree between ADNI sources in both tables: coded male
in PHC and in the pipeline, but `Female` in ADNIMERGE. This is unrelated to the swap.

**Hand check**, first row of the CSV: recomputing it by hand with the coefficients
below, the opposite-sex equation reproduces the published value to the second
decimal; the own-sex equation does not.

**Pitfall when re-implementing in R/data.table.** Do not name a function argument
the same as a column of the coefficient table. `cf[match(g, cf$g)]` silently uses the
column `g`, not the argument, and returns garbage. An earlier check in this
investigation hit exactly this. `docs/frs_sex_swap_check.R` avoids it (`sex.v` vs
`sex_lab`).

**The coefficients.** Taken from CVrisk's `frs_simple_coef`
(D'Agostino et al. 2008, *Circulation* 117:743-753, office-based BMI model):

| | ln age | ln BMI | ln SBP untreated | ln SBP treated | smoker | diabetes | group mean | S0(10) |
|---|---|---|---|---|---|---|---|---|
| Male | 3.11296 | 0.79277 | 1.85508 | 1.92672 | 0.70953 | 0.53160 | 23.9388 | 0.88431 |
| Female | 2.72107 | 0.51125 | 2.81291 | 2.88267 | 0.61868 | 0.77763 | 26.0145 | 0.94833 |

Risk = round((1 − S0^exp(Σβx − group mean)) × 100, 2), clamped to [1, 30].

---

## 4. What in this repository is affected

### Code that reads the published variable

| Where | What it does with FRS |
|---|---|
| `config/pipeline_config.yaml`, `parameters: frs: column:` (about line 177) | names `PHC_ASCVD_10y_FRS_Simple_Ageover30` |
| `R/scripts/03_prepare_analysis_data.R` (about lines 268-280) | merges FRS into `data/derivatives/lme/cohort_hvr.rds` |
| `R/scripts/04_frs_problem_analysis.R`, Part 1 | ceiling, age-FRS correlations, variance decomposition (R² age vs risk factors) |
| `R/scripts/07_cvr_score_extraction.R` (line 148; lines 260-299; lines 380-384) | convergent validity r(η, FRS), r(ζ, FRS) overall and by sex; z-scores FRS into the cohort |
| `R/scripts/10b_lme_amyloid_negative.R` (lines 214-260) | reads FRS raw from the CVRF CSV and z-scores it (A−/CU models) |
| `R/scripts/15_prepare_manuscript_env.R` (lines 245-268) | `frs_ceiling_all/_female/_male`, `frs_above20_all` |
| Downstream consumers of the cohort FRS | `08_lme_hvr_z.R`, `09_lme_hc_z.R`, `10_lme_sensitivity.R`, `11_prepare_lgcm_data.R`, `12_lgcm_parallel.R`, `13_lgcm_mediation.R`, `13b_lgcm_amyloid_negative.R`, `export_submission_figures.R`, `R/utils/plotting.R` (`plot_frs_ceiling_by_sex`), `R/utils/tables.R` |

### Results that must be re-run

- **The ceiling by sex.** Manuscript inline values at about lines 1570-1576 of
  `reports-src/manuscript_ad_submission.qmd`, the figure "Cardiovascular risk
  characterization" (panel A, FRS distribution by sex), and any text that interprets
  which sex is at the ceiling.
- **The FRS LGCM a-path** ("FRS predicted hippocampal atrophy in 6 of 6 models"), the
  b-paths and indirect effects for FRS.
- **The FRS LME three-way interactions** (years × FRS × HVR), and their FDR counts.
- **A−/CU FRS models** (the "FRS none" comparison in script 10b / 13b).
- **FRS convergent validity** with η and ζ (script 07), including the by-sex values.
- **Variance decomposition** of FRS into age and risk factors (script 04 Part 1c;
  manuscript inline `frs_var_decomp`).
- The EDT (disease-time) FRS models (`lme_edt_hvr_z_frs`) and the quadratic attempts
  (`lgcm_frs_quadratic_attempts`).

### Not affected

- **CVR_mimic (the MIMIC measure).** FRS is not an indicator. Indicators are SBP,
  hypertension, glucose, cholesterol and creatinine (scripts 05 to 07). Only FRS
  *comparisons* with it change.
- **The measurement invariance tests** (script 06), which do not use FRS.
- **The PREVENT-AD power simulation** (`16_power_preventad.R`). Its anchors are
  CVR_mimic terms; FRS appears only in a descriptive sentence (line 1433).
- **HVR measures and cognition.**

### Frozen snapshot

`archive/2026-05-19_initial_submission/` (tag `initial-submission`, commit
`0a1595a`) is the as-submitted state for the *Alzheimer's & Dementia* initial
submission. It contains the swapped results. Use it as the "before" for a diff.

---

## 5. Pre-existing issues that interact with the fix

1. **`PHC_Diabetes` is itself suspect.** Script 04, Part 2, already documents that
   `PHC_Diabetes`, `PHC_Heart` and `PHC_Stroke` are nearly collinear (phi 0.91 to
   0.96), with implausible prevalences (about 28-29%). It also documents
   HTN–Diabetes phi = 0.004 and Age–Diabetes r = −0.17. FRS uses `PHC_Diabetes` as a
   direct input. A corrected FRS that still uses `PHC_Diabetes` fixes the sex error
   but not this one. Script 04, Part 3 verifies against raw ADNI files; consider
   deriving diabetes from raw ADNI (MEDHIST / RECMHIST / medications) for the
   recomputed score, and report both.
2. **Script 04's ceiling check uses the wrong unit.** It tests
   `abs(FRS - 0.30) < 0.001`, but the PHC values are in percent (4.21 to 30.00). The
   saved `models/results/lme/frs_problem_analysis.rds` therefore reports
   0 participants at the ceiling. The manuscript does not use this output: it takes
   the ceiling from script 15, which correctly tests `== 30`. Fix while rerunning.
3. **The FRS validity age range is inconsistent.** `config` (`valid_age_max: 62`) and
   the comment in script 03 (about line 282: "validated for ages 30-62 (D'Agostino
   et al., 2008)") disagree with the CVrisk documentation ("between 30 and 74") and
   with the applicant's grant materials (30 to 74). Check against D'Agostino 2008
   before relying on either.
4. **Ages over 74.** ADSP-PHC extrapolated the equation beyond 74 (§2). Decide whether
   the recomputed score does the same (to stay comparable to PHC) or follows stock
   CVrisk (`NA` above 74). In ADNI that would remove most of the cohort, so state the
   choice explicitly.

---

## 6. Suggested fix plan

Nothing below has been done yet. Corrected numbers must come from the re-run
outputs (RDS files), not be hardcoded.

1. **Compute a corrected FRS inside the pipeline.** Use PHC inputs (age, BMI, SBP,
   BP medication, smoking, diabetes) with the correct sex: a small utility in
   `R/utils/`, or add CVrisk to `renv` and call it with an explicit
   `gender = ifelse(PHC_Sex == 1, "male", "female")`.
   - Do not overwrite the published variable. Add a new column (e.g. `FRS_corrected`)
     and a config switch, so the published-variable results stay reproducible for
     the correction notice.
   - Add a unit test reproducing the hand check, and an assertion that the
     recomputed score differs from the PHC variable as in §3.
2. **Fix the script 04 unit bug** (§5.2) and the age-range inconsistency (§5.3).
3. **Re-run in order:** 03 → 04 → 07 → 08 → 09 → 10 → 10b → 11 → 12 → 13 → 13b →
   15 → figures and tables.
   - Set the relevant `force_regenerate` flags in config. Script 04 skips itself if
     its outputs exist.
   - The LGCM bootstraps are slow.
4. **Diff against** `archive/2026-05-19_initial_submission/`. Report which FRS
   conclusions survive: the 6/6 a-path, the A−/CU FRS null, and the FRS LME
   interactions.
5. **Rewrite the manuscript's ceiling and sex framing** from the new outputs. The
   argument that Framingham saturates in an older cohort survives (54.8% at the cap),
   but now it saturates mostly in **men**.
6. **Sensitivity analysis:** the corrected FRS with diabetes from raw ADNI (§5.1).

---

## 7. Implications outside this repository

- **The *Alzheimer's & Dementia* submission** (initial submission 19 May 2026)
  reports the swapped ceiling and FRS results.
- **The medRxiv preprint** doi:10.64898/2026.06.01.26354601 (PubMed PMID 42282182)
  reports 74.9% women / 28.1% men and the FRS model results.
- **The AAIC 2026 materials** (`AAIC-2026.pdf`). They include "reduced FRS
  variability (76% at ceiling)", and the poster title "Cardiovascular Risk Decouples
  Hippocampal Structure From Cognitive Function In Amyloid-positive Adults" rests on
  FRS effects.
- **The grant applications** (separate repo `/ipl/ipl27/sfernandez/grants`). The CIHR
  CPRA 2027 package (submitted, cannot be changed) and the FRQS and AARFA drafts use
  the "75% of women vs 28% of men" framing. `grants/FACTS.md` now flags it `[!]`.
- **Co-authors** (Villeneuve, Collins) need to know before any correction.
- **ADSP-PHC and other users.** The methods document says the same derivation was
  applied to ADNI, EFIGA, HABS-HD, NACC, ROS, MAP, MARS, TARCC, WHICAP and WRAP. Other
  cohorts' files have **not** been checked; HABS-HD should be checked before any use
  in the AARFA application. The methods PDF lists Dr. Timothy Hohman (Vanderbilt) as
  the contact. Whether and when to report is not yet decided.

---

## 8. Open questions

1. Is the swap also present in the other ADSP-PHC cohorts (start with HABS-HD, which
   is already in hand)? Run the §3 comparison on each cohort's CVRF file.
2. Did the V1 release (Dec 2023) have it, and does any newer ADSP-PHC release fix it?
3. Is `PHC_Diabetes` usable at all (§5.1)?
4. After the re-run: does the FRS a-path to hippocampal atrophy still hold in each
   sex, and does the A−/CU contrast with CVR_mimic still hold?
