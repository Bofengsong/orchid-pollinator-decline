# Data definitions

## `data/current/MainData.xlsx`

Full evidence set. Each row is a pollinator species observation record. The
table contains 729 records from 164 unique source publications.

## `data/current/MetaData.csv`

Meta-analysis table. Each row is a study-level log response ratio (`lnRR`),
with its sampling-variance proxy and hierarchical identifiers. The table
contains 257 effect sizes from 57 publications.

Key model fields include:

- `reference_no`: source-publication identifier and publication-level random effect.
- `country`: country-level random effect.
- `pollinator_family`: taxonomic-family random effect.
- `yi`: effect size, defined as `ln(1 + decline_value)` in the locked workflow.
- `vi`: variance proxy, defined as `1 / sample_size` in the locked workflow.

## `data/current/MainMetaDataLong.csv`

Driver-level long table used for driver analyses and plotting. One effect size
may appear under more than one reported driver, producing 610 driver-effect
rows. This table must not be treated as 610 independent effect sizes.

## `data/source_data/Orchid_Pollinator_Decline_Source_Data.xlsx`

Reviewer-facing workbook containing:

- publication-level source information;
- all 729 evidence records;
- all 257 meta-analysis records;
- all 610 driver-effect rows;
- extraction notes and source-location fields where available.

