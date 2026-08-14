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

## SHA-256 checksums for canonical tables

```text
bbaf7f5296178a88067bfb2df81f6a5f7de3a1c46688b8455f894dbbc9046295  data/current/MainData.xlsx
a1f7b15bb8ee9499f995f647a8678377cc2f8c78d34f675324a0b5055ec5625d  data/current/MetaData.csv
0a6abe12e8759aa31a13e632d6d84eb28264f384e6690c86ff3da2e17a35d8bc  data/current/MainMetaDataLong.csv
1adbf44b87b855d65e53357cda292630655922070edd3c176d3915ed7486c0c2  data/source_data/Orchid_Pollinator_Decline_Source_Data.xlsx
```

