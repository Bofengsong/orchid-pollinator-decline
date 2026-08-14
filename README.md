# Orchid pollinator decline: analysis and source data

Private peer-review repository for the manuscript:

> *Where have all the orchids gone? A global synthesis of orchid pollinator declines and implications for conservation*

This repository contains the analysis-ready evidence set, meta-analysis data,
driver-level long data, source-data workbook, and R scripts used for the main
statistical and spatial analyses. Manuscript files, reviewer correspondence,
downloaded publications, and historical working files are intentionally
excluded.

## Repository structure

- `code/`: R scripts for the meta-analysis, Figures 1-4, and the PRISMA diagram.
- `data/current/`: canonical analysis tables.
- `data/source_data/`: publication-level and row-level source-data workbook.
- `data/spatial_inputs/`: processed spatial inputs required by the supplied scripts.
- `data/audit/`: locked PRISMA count table.
- `docs/`: data definitions, run order, and source notes.
- `outputs/`: generated files (not versioned; created when scripts are run).

## Locked dataset counts

- 164 unique publications and 729 pollinator species observation records in
  the full evidence set.
- 57 unique publications and 257 study-level lnRR effect sizes in the
  multilevel meta-analysis.
- 610 driver-effect rows in the driver-level long table.

The unit of analysis is not a publication count. A publication may contribute
multiple pollinator species observation records and, where quantitative
requirements are met, multiple effect sizes. Publication identity is included
as a random effect in the multilevel models.

## Quick start

Install the R packages listed in `docs/R_PACKAGES.md`, then run scripts from the
repository root in the order given in `docs/RUN_ORDER.md`. All scripts resolve
paths relative to their own location; no private local paths are required.

## Data and code availability

This repository is private during editorial assessment and peer review. Data
and code will be deposited in a public archival repository upon acceptance.

