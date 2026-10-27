# Orchid pollinator decline: data and R code

Data and code for **Where have all the orchid pollinators gone? A global synthesis of orchid pollinator declines and implications for conservation**.

The evidence set contains **595 observations from 170 source publications**. The primary meta-analysis uses **251 effect sizes from 71 source publications**. Twelve seasonal or within-cycle comparisons remain in the evidence set and are excluded from the primary meta-analysis. The pooled estimate is lnRR = −0.87894, equivalent to a 58.5% reduction relative to the study-specific reference or comparator.

## Files

- `data/Final_data.csv`: the single analysis table. `meta_included = 1` identifies the primary Meta subset.
- `data/Supporting_data.rds`: processed spatial inputs, comparison intervals and supporting inputs used by the analyses.
- `data/source_data/Orchid_Pollinator_Decline_Source_Data.xlsx`: the same data in a single-sheet Excel workbook.
- `code/Analysis.R`: statistical analyses and spatial prioritisation.
- `code/Figure1.R`–`Figure4.R` and `code/FigureS1.R`–`FigureS11.R`: one script for each statistical figure.
- `outputs/`: results and the figures used in the manuscript. Figure 5 is an explanatory illustration and has no statistical plotting script. The S1 PDF and PNG show the same figure.

## Run

Download the complete repository, keeping `code`, `data` and `outputs` together. In RStudio, open `code/Analysis.R` and click **Source**. Then open and source the required figure script. 图形会保存到 outputs，同时可在 RStudio 的 Plots 中预览。

From the repository root:

```r
source("code/Analysis.R", encoding = "UTF-8")
source("code/Figure2.R", encoding = "UTF-8")
```

Or use `Rscript code/Analysis.R` from a terminal. Scripts locate the repository from their own location; no local OneDrive paths are needed. If data change, run `Analysis.R` again before plotting. Saved results allow the figures to be drawn without refitting models.

R package requirements and concise data definitions are in `docs/`. Raw source publications and private review documents are not part of this repository.
