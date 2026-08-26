# Validation record

Validated on 14 August 2026 with R 4.5.1.

## Dataset counts

- Evidence set: 729 rows; 164 unique `reference_no` values.
- Meta-analysis: 257 rows; 57 unique `reference_no` values.
- Driver-level long table: 610 rows.
- Source workbook: 164 publication rows, 729 evidence rows, 257 meta rows,
  and 610 driver rows.

## Script checks

All seven scripts completed successfully from the repository layout and wrote
their expected files under `outputs/`.

The principal model reproduced:

- pooled lnRR = -1.3708978;
- 95% CI = -1.9526333 to -0.7891622;
- k = 257; publications = 57;
- regional moderator: Qm = 11.2067, df = 5, p = 0.04743.



