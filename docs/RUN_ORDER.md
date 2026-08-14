# Run order

Run from a clean R session. Each script writes to a subdirectory of `outputs/`.

1. `Rscript code/01_run_meta_analysis.R`
2. `Rscript code/02_make_figure2.R`
3. `Rscript code/03_validate_predicted_risk_confidence.R`
4. `Rscript code/04_make_priority_map.R`
5. `Rscript code/05_make_wcvp_lineage_maps.R`
6. `Rscript code/06_make_prisma_figure.R`
7. `Rscript code/07_make_figure1_base_maps.R`

Notes:

- Script 03 validates and stages the processed 1-degree risk layers used by the
  priority-map workflow. Regeneration from the very large original climate,
  pesticide and habitat rasters is outside this compact review package; source
  products are documented in `docs/SPATIAL_SOURCE_NOTES.md`.
- Script 07 generates the data-driven Figure 1 base maps. The small pollinator
  pictograms in the submitted figure are graphical legend elements and do not
  alter the mapped data.

