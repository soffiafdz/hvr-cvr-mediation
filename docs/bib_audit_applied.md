# Bibliography Audit — Applied Corrections

Source: `/ipl/ipl27/sfernandez/thesis_phd/BIB_AUDIT.md` (Crossref-verified
audit of the thesis bib, which is the union of paper1–paper4 bibs).
Date applied: 2026-05-25.

## Summary

21 corrections from the upstream audit were evaluated against this repo's
`reports-src/references.bib`. **10 corrections were applicable and applied;
11 entries did not exist in this repo and were skipped (no invention).**

## Applied corrections

For each entry, the `author = {...}` field was replaced wholesale with the
Crossref-verified value. Entry keys preserved unchanged.

| DOI | Entry key | Nature of fix |
|---|---|---|
| 10.7326/0003-4819-157-6-201209180-00002 | `koller_2012` (around line 583) | Replaced `Koller and Leening and Wolbers and Steyerberg and Et., Al` with full 11-author list |
| 10.1037/0894-4105.21.4.412 | `kramer_longitudinal_2007` | Replaced `Kramer and Mungas and Reed and Wetzel and Et., Al` with full 8-author list |
| 10.1136/jech.57.8.634 | `marrugat_adaptation_2003` | Replaced single-token `Marrugat` with full 12-author list |
| 10.1037/1082-989x.5.1.23 | `mehta_putting_2000` | Added first names: `Mehta, Paras D.` and `West, Stephen G.` |
| 10.1007/s11336-014-9435-8 | `neale_openmx_2016` | Replaced `Neale and Hunter and Pritikin and Zahery and Et., Al` with full 10-author list |
| 10.12968/bjca.2010.5.5.47882 | `ouldred_vascular_2010` | Added first names: `Ouldred, Emma` and `Bryant, Catherine` |
| 10.1007/BF02296192 | `satorra_scaled_1994` | Added first names: `Satorra, Albert` and `Bentler, Peter M.` |
| 10.1212/01.wnl.0000341271.90478.8e | `viswanathan_vascular_2009` | Added first names: `Viswanathan, Anand`, `Rocca, Walter A.`, `Tzourio, Christophe` |
| 10.1016/S0140-6736(24)01296-0 | `livingston_dementia_2024` | Added 6 missing authors (Shirai, Singh-Manoux, Schneider, Walsh, Yao, Sommerlad); Mukadam moved from position 21 to 27 |
| 10.1016/j.neuroimage.2019.116108 | `schoemaker_hippocampal_ventricle_2019` | Added `J.` middle initial: `Lupien, Sonia J.` |

## Skipped (entry not present in this repo)

The following DOIs from the upstream audit have no matching entry in
`references.bib` — skipped per the source-of-truth rule (no invention).

| # | DOI | Paper |
|---|---|---|
| 1 | 10.1037/0882-7974.8.2.156 | Morse 1993 |
| 2 | 10.3389/fninf.2011.00037 | Das et al. 2011 (LORIS) |
| 12 | 10.3389/fnagi.2017.00445 | Todd et al. 2017 |
| 13 | 10.1002/hbm.26407 | Morrison et al. 2023 |
| 14 | 10.3389/fnins.2016.00558 | Backhausen et al. 2016 |
| 15 | 10.1002/hbm.23894 | Dadar et al. 2018 |
| 16 | 10.3389/fninf.2016.00052 | Pizarro et al. 2016 |
| 17 | 10.1002/hbm.24811 | Goubran et al. 2020 |
| 18 | 10.1212/WNL.49.3.786 | Jack et al. 1997 (also searched by title; not found) |
| 19 | 10.1212/WNL.55.4.484 | Jack et al. 2000 (also searched by title; not found) |
| 21 | 10.1002/hbm.25320 | Dima et al. 2022 (ENIGMA Lifespan) |

## Sanity checks (post-change)

```text
grep -c "and others" reports-src/references.bib  → 0  (baseline: 0)
grep -c "Et\..*Al"   reports-src/references.bib  → 0  (baseline: 3)
```

Three `Et., Al` artefacts (Neale, Kramer, Koller) were eliminated by the
corresponding full-author replacements.

## Review state

Changes left unstaged for user review per the upstream task instruction.
Run `git diff reports-src/references.bib` to inspect the 10 author-line
swaps (10 insertions, 10 deletions, no other field touched).
