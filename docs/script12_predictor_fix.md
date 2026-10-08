---
title: "Script 12 predictor specification fix"
date: "2026-10-06"
author: "Sofia Fernandez-Lozano"
---

## What script 12 models

Script 12 fits six parallel-process latent growth curve models (2 sexes × 3 cognitive domains). The result the manuscript reports from it is the coupling between hippocampal and cognitive decline.

Each model has two trajectories over estimated disease time (EDT):

- **HVR trajectory:** a person-level intercept (`hvr_i`), a person-level slope (`hvr_s`) and a population-level quadratic term (the "fixed quadratic").
- **Cognitive trajectory:** the same three parts (`cog_i`, `cog_s`, quadratic).

Two kinds of observed variable feed the growth factors:

- **Exogenous predictors:** baseline age, education, APOE ε4 and FRS. "Exogenous" means nothing in the model predicts them.
- **Baseline HVR (`HVR_Z_T0`)**, which also predicts the cognitive slope. It is a measured outcome as well (the first HVR time point), so it is not exogenous.

The **coupling** is the covariance between `hvr_s` and `cog_s` after the predictors are accounted for. It answers whether people whose hippocampus declines faster also decline faster cognitively.

## The issue

The predictors' means and variances were fixed at values the data do not have. This did not bias the estimates, but it made the models numerically fragile and made the stored latent means meaningless.

### How the predictors were specified

OpenMx fits these models by full-information maximum likelihood on all observed variables, predictors included. Every variable listed as observed in a RAM model gets a distribution, so OpenMx needs a mean and a variance for each predictor. The code fixed them as follows:

| Predictor | Model assumed | Men, actual | Women, actual |
| --- | --- | --- | --- |
| Age (raw years) | mean 0, variance 1 | mean 74.8, variance 46.5 | mean 72.6, variance 51.5 |
| Education (raw years) | mean 0, variance 1 | mean 16.5, variance 7.2 | mean 15.6, variance 7.4 |
| APOE ε4 (0/1/2) | mean 0, variance 1 | mean 0.60, variance 0.24 | mean 0.62, variance 0.24 |
| Correlations between predictors | all 0 | FRS–age r = 0.48 | FRS–age r = 0.41 |

### Why the estimates were not biased

The likelihood splits into two parts: the probability of the predictors, p(x), and the probability of the outcomes given the predictors, p(y | x).

Every predictor mean and variance was fixed, so p(x) had nothing to estimate. It was the same number whatever values the other parameters took. Maximising the total therefore gives exactly the estimates of p(y | x) alone, the standard "fixed-x" regression. Adding a constant to a score moves the total but not where its maximum is.

So the paths and the coupling were always proper partial effects, adjusted for the other predictors. The refits confirmed this: every model reaches the same optimum (next section).

### Where it did damage: numerical fragility

Two features made the optimisation fragile.

1. **A huge constant.** Under a standard normal, each person's age adds roughly age² to −2LL, about 75² ≈ 5,600 per person. The predictor term came to about 3.0 million in men and 2.2 million in women. The informative part of the likelihood was about 3,000 to 4,500, so more than 99.8% of the number being optimised was a constant.
2. **Intercepts defined at age 0 with no education.** With raw predictors, each latent mean is the value for a 0-year-old with 0 years of education, about 75 years from any participant. The intercept and the age path then become almost interchangeable: nudge the age path and the intercept must move about 75 times as much to keep predictions in range. The likelihood surface becomes a long, narrow valley, which optimisers and Hessians handle badly.

The 2 October run showed the symptoms:

- **OpenMx status code 5** ("the Hessian at the solution does not appear to be convex", Mx status RED), 8 times, on the fixed-quadratic models. The reported coupling p-values are Wald tests, which rely on that Hessian.
- **Missing SEs:** 10 of 36 parameters in Male Memory, 2 in Male Executive Function.
- **Failed nested tests:** the coupling LRT returned a result in only 2 of 6 models; the linear-vs-quadratic comparison in only 2 of 6.
- **A difficult fit:** Male Memory ended with "Retry limit reached" (3 of 6 attempts errored).

These runs cannot separate how much came from the constant and how much from the age-0 parameterisation, because the fix removed both together.

### The latent means were extrapolations

Because the intercepts referred to age 0, the stored slope means did too. Male Memory's HVR slope mean was **+0.157**, a rising HVR for an imaginary newborn. For a participant of average age, education and APOE ε4 status, the same model gives about **−0.21**, the decline one would expect. The manuscript's growth sentence was meant to report these means; it only escaped publishing them because a bug made that chunk print nothing.

## The fix, and why each part works

Three changes to `R/scripts/12_lgcm_parallel.R` remove the fragility while leaving the model itself unchanged; all six models reach the same optimum as before.

### 1. Predictor moments fixed at their sample values (`exo_moments.fn`)

The predictors' means, variances and covariances are now fixed at what the data show, which is lavaan's `fixed.x` behaviour.

- They are still fixed, so p(x) is still constant and the estimates are still the fixed-x estimates.
- The constant drops from millions to a size that fits the data: −2LL is now about 10,000 (Male Memory 10,018).
- The predictor correlations are no longer forced to 0, so the joint model describes the data honestly.
- The model-implied predictor moments were checked to equal the sample moments exactly.

### 2. Every predictor centred within each sample

Centring is a reparameterisation: it changes what the intercept means, not the model.

- **Paths are unchanged.** In the refit they moved by at most about 10⁻⁵.
- **Latent means become sample averages.** They now describe a participant with the sample's average age, education, APOE ε4 status and FRS.
- **The narrow valley disappears.** With the intercept defined in the middle of the data, the intercept and the predictor paths no longer stand in for each other, so the Hessian is well conditioned.

### 3. Two starting points (`fit_fixedquad.fn`)

`mxTryHard` stops at the first acceptable solution, even a local optimum. With the new specification, Male Executive Function's default start settled about **54.5 units worse** than the best solution.

The fixed-quadratic model is now fitted twice and the lower −2LL is kept:

- from the default starting values;
- from the estimates of its linear version (quadratic means fixed at 0), a simpler nested model whose solution sits close to the full model's.

Both starts agree to two decimals in 5 of 6 models. In Male Executive Function the linear start wins (11,495.14 vs 11,549.60). The old specification happened to find the better solution from its default start; two starts make that independent of luck.

### How we know it is the same model

The old and new −2LL cannot be compared directly, because each contains a different constant p(x). That constant can be computed straight from the data, so it was subtracted from both. What remains is the conditional −2LL, the part that depends on the parameters:

| Model | Old | New | Difference |
| --- | --- | --- | --- |
| Male Memory | 2946.28 | 2946.28 | 0.00 |
| Male Language | 3940.61 | 3940.61 | 0.00 |
| Male Executive Function | 4423.18 | 4423.18 | 0.00 |
| Female Memory | 2559.42 | 2559.23 | −0.19 |
| Female Language | 2908.61 | 2908.82 | +0.21 |
| Female Executive Function | 3356.86 | 3356.84 | −0.02 |

The women's ±0.2 comes from one woman missing APOE ε4. Under FIML her row is averaged over the assumed APOE ε4 distribution, which used to be the wrong N(0, 1) and is now the sample's. Her row cannot be split exactly, so the subtraction is approximate for the women. The same correction moves the women's APOE ε4 → intercept paths by about 0.01, and nothing else.

The final run, with all four predictors centred, gives the same values.

## Factor scores, and why the figure now uses implied means

The coupling figure's population curves now come from model-implied latent means, because averaged factor scores proved unreliable at the level of a few hundredths. No reported estimate depends on factor scores.

**What they are.** Factor scores are each person's estimated `hvr_i`, `hvr_s`, `cog_i` and `cog_s`, given that person's data. Only the coupling figure uses them: panel A's population mean curves were the *average* factor score at each point, and panel B's scatter plots each person's two slopes.

**What happened.** With only age and education centred, the average factor scores for `hvr_i`, `hvr_s` and `cog_i` shifted by exactly the APOE ε4 term. For Male Memory `hvr_i` that is −0.138 × 0.602 = −0.083: the average came out as if every participant were an APOE ε4 non-carrier. This was established from the numbers, not traced inside `mxFactorScores`.

**Why centring all four predictors fixes those three.** Once every predictor has mean 0, there is no mean contribution left to drop. In a Male Memory test, `hvr_i`, `hvr_s` and `cog_i` then matched the old factor scores to within 10⁻⁵.

**Why `cog_s` is different.** `cog_s` is the only latent with a non-exogenous observed variable pointing at it: baseline HVR, which averages −2.27 and cannot be centred away because it is also an outcome. The level of its averaged factor scores depends on the parameterisation, and none of the versions is exact. Individual positions barely move (r = 0.98 between old and new).

| Source | Mean `cog_s`, Male Memory |
| --- | --- |
| Model-implied, average participant | −0.1351 |
| Old factor scores (raw predictors) | −0.1386 |
| Factor scores, all predictors centred | −0.1273 |

**The change (option b).** Script 12 now stores `implied_means` for each model, computed directly from the fitted parameters: each latent's mean parameter plus, for every observed predictor, path × that predictor's sample mean. The centred predictors contribute 0; baseline HVR contributes its path × −2.27. The figure draws its population curves from these values, so they are exact by construction. The scatter still shows factor scores, since individual positions are what it displays, and its printed correlations come from the model's covariance estimates, which are unchanged.

## What it means

No reported result changes; the coupling inference becomes sturdier, and the growth parameters become reportable.

**Reported results are unchanged.** Coupling covariances, SEs, p-values and correlations are identical. The coupling is significant in all three male domains (Memory p = 0.0014, Language 0.0082, Executive Function 0.0012) and, among women, only in Executive Function (p = 0.0034). The mediation results (script 13) were checked separately and are unaffected there too: paths matched to within about 10⁻⁵ (`docs/lgcm_exogenous_check.R`).

**The coupling inference is sturdier.** The fits are clean, with no "Hessian not convex" warnings and no missing SEs. A likelihood-ratio test of the coupling now exists for all six models and agrees with the Wald tests. Wald tests depend on a well-behaved Hessian and LRTs do not, so reporting both is a real robustness gain, the kind reviewers ask for.

| Model | LRT χ² | p |
| --- | --- | --- |
| Male Memory | 10.77 | 0.0010 |
| Male Language | 7.26 | 0.0071 |
| Male Executive Function | 10.95 | 0.0009 |
| Female Memory | 2.13 | 0.14 |
| Female Language | 1.63 | 0.20 |
| Female Executive Function | 9.03 | 0.0027 |

**Growth parameters become reportable.** Slope means now describe an average participant. HVR slopes are −0.21 per EDT year in men and −0.17 in women, and cognitive slopes −0.11 to −0.16 across domains. The manuscript's growth sentence can be repaired to report them, and its claim ("negative across all domains") will then be true.

**The figure's mean curves become exact model curves.** Expect a slight visual shift, mainly in the cognitive curves.

**The fixed-quadratic justification stands.** The fully random quadratic models still fail under the corrected specification, so the manuscript's reason for the fixed-quadratic model was not caused by this issue.

**Still outstanding:**

- **Scripts 13 and 13b** use the same moment specification. It is harmless there because age was already centred (the constant is about 25,000), but they should match at the next full re-run, with a checkpoint flag added.
- **The manuscript's fit table** (from script 13) still prints the inflated absolute −2LL, AIC and BIC. It should be dropped or reworked.

## Status

Done and committed on 6 October on branch `fix/frs-sex-swap` (`0b7de98` script 12, `12189b8` figure). The all-predictor re-run of scripts 12 and 15 finished at 15:45, and `outputs/manuscript_env.rds` now matches the committed code.

The re-run passed every check:

- [x] All six models match the old optimum on the conditional −2LL (men 0.00; women within ±0.2).
- [x] Coupling estimates, SEs and p-values unchanged; no missing SEs; no status-5 warnings; coupling LRT for all six models.
- [x] `hvr_i`, `hvr_s` and `cog_i` factor scores match the old ones (r = 1.0000).
- [x] `implied_means` stored for every model; they equal the factor-score means for `hvr_i`, `hvr_s` and `cog_i`, and Male Memory `cog_s` = −0.1351.
- [x] Coupling figure drawn and inspected: HVR and Memory curves unchanged; Language and Executive Function curves shift slightly; scatter correlations unchanged.

Male Executive Function again needed the second start (default 11,921.44; linear start 11,495.14). Pre-fix outputs remain in `archive/2026-10-05_pre_exo_fix_script12/`.
