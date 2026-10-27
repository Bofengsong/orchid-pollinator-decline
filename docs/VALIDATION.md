# Validation record

Checked on 27 October 2026 with R 4.5.1. Package versions are in `outputs/sessionInfo.txt`.

- Evidence: 595 observations from 170 publications.
- Primary Meta: 251 effects from 71 publications; all effects are finite and all variances are positive.
- The single-sheet Excel workbook matches `Final_data.csv` in every cell (595 rows, 70 columns).
- `Analysis.R` and all 15 figure scripts ran successfully from a working directory outside the repository.
- Recomputed model summaries, sensitivity analyses, bootstrap risk and spatial priorities match the manuscript results to numerical tolerance 1e-10, excluding runtime metadata.
- All 15 generated figures match the manuscript figure images pixel for pixel. Figure 5 is the retained ecological illustration.

The principal model gives lnRR = -0.8789369, equivalent to a 58.5% decline. The order moderator test gives Qm = 6.354331 and p = 0.09558778. The publication-level Egger test gives z = 0.09055293 and p = 0.9278478.

Warnings were reviewed: sparse method subgroups have random-effect components fixed at zero when only one level is present; the all-scaled annualisation sensitivity has a large variance ratio; the map-centroid warning concerns constant attributes during geometry conversion. These warnings do not change the reproduced numerical results.

Figure S1 uses the confirmed screening totals: the original path retains 136 publications, and the supplementary path retains 34 (129 full texts assessed minus 95 excluded). Their combined evidence comprises 170 publications and 595 observations; 71 publications contribute 251 primary effects.
