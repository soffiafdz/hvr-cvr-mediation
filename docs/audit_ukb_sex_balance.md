# UKB Normative Sample — Sex-Balance Audit

**Audit date:** 2026-05-14
**Auditor:** S. Fernandez-Lozano
**Reproducible:** `R/scripts/audit_ukb_sex_balance.R` (saves
`outputs/audit/ukb_sex_balance.rds`)

## Question

The UKB normative reference used to z-score ADNI HVR/HC values
(`data/external/gamlss.rds`) is **56.2% female** (15,412 F vs 12,022 M
total; 12,327 F vs 9,615 M in the train split reported in
@tbl-demographics). Why?

## Verified numbers

| Stage | Female | Male | Total | % F |
|---|---:|---:|---:|---:|
| 1. Raw imaging cohort (ses-2 or ses-3 MRI) | 38,861 | 35,976 | 74,837 | 51.9 |
| 2. DARQ QC failures (DARQ < 0.25) | 3,252 | 5,372 | 8,624 | 37.7 |
| 3. AssemblyNet QC failures (grade ≠ A) | 163 | 114 | 277 | 58.8 |
| 4. ICD-10 F-code (any) | 1,007 | 899 | 1,906 | 52.8 |
| 5. ICD-10 G-code (any) | 1,233 | 982 | 2,215 | 55.7 |
| 6. ICD-10 Q0-code (any) | 0 | 0 | 0 | — |
| 7. ICD-10 primary exclusion (F\|G\|Q0) | 2,240 | 1,881 | 4,121 | 54.4 |
| 8. Final GAMLSS sample (all splits) | 15,412 | 12,022 | 27,434 | 56.2 |
| 9. Final GAMLSS sample (train only) | 12,327 | 9,615 | 21,942 | 56.2 |

QC failure rates within sex:

| Sex | N (scans w/ QC) | DARQ fail % | AssemblyNet fail % |
|---|---:|---:|---:|
| Female | 25,002 | 13.01 | 0.65 |
| Male | 22,394 | 23.99 | 0.51 |

## Mechanism

The 56:44 imbalance has two stacked drivers:

1. **The UKB imaging cohort is already 52% female** (38,861 F vs 35,976 M).
   UK Biobank's imaging extension enrolled more females than males at
   source.
2. **DARQ QC fails at ~1.85× the rate in males** (24.0% vs 13.0%), removing
   ~2,100 more males than females during the QC stage. This is the
   dominant cause of the post-QC shift.

ICD-10 exclusions remove slightly more females than males overall (2,240 F
vs 1,881 M = +359 more F excluded), so they nudge the ratio *back* toward
parity, but the differential is small compared to DARQ's.

## What is NOT verified by this audit

These are plausible inferences that the data here cannot confirm:

- **The mechanism of DARQ failures in males.** DARQ is a deep-learning QC
  score; whether it is tracking motion specifically, intensity outliers,
  or some other artifact is not visible from these data. Citing "male
  motion" without checking the DARQ paper would be speculation.
- **The mechanism of the F-code sex pattern.** Females have slightly more
  F-code (psychiatric) flags in the imaging cohort, but whether this
  reflects underlying prevalence, healthcare-contact differences, or
  diagnostic patterns is not derivable from these counts.

## Conclusion

The sex imbalance is real and reproducible, originating in the upstream
UKB normative pipeline. **It does not bias z-scores in this manuscript**
because GAMLSS normative models are fit **separately by sex** — each sex's
reference curves are built only from its own data. The imbalance only
affects relative statistical support per sex (slightly less precision for
male tails of the normative distribution).

## Caveats on the audit

- Two upstream filters were not reproduced in this script:
  - **Outlier removal** (|Z| > 3 on volume residuals from
    `lm(VAL ~ SEX + AGE + ICC)`). The regression is sex-adjusted so the
    cut is approximately symmetric within sex.
  - **Manually-curated bad-EID list** loaded via `exclude_ids_file` in
    `05_adjust_headsize.R`. Small in absolute numbers.
- Numbers will shift slightly if UKB derivatives are updated.

## How to reproduce

```bash
# fst is not in this project's renv; the upstream UKB project has it.
cd /data/ipl/ipl27/sfernandez/headsize_sexeffects_ukb
Rscript /ipl/ipl27/sfernandez/hvr_cvr_mediation/R/scripts/\
audit_ukb_sex_balance.R
```

Output: `outputs/audit/ukb_sex_balance.rds` and a console summary.

## Upstream transfer

See `docs/upstream_ukb_audit_handoff.md` for the plan to move this audit
into the upstream UKB-paper repository, where it is a more natural fit.
