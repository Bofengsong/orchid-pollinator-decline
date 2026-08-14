# Figure 3 Kingsley mask update and WCVP masked previews

Generated on: 2026-05-27

## What changed

1. The previous strict orchid-presence Figure 3 was retained as the starting point.
2. Additional cells annotated with crosses in Kingsley's iPad review were matched to the nearest currently coloured 2-degree grid cells and removed as NA.
3. The resulting Figure 3 keep/NA layer was used as a global mask for all WCVP subfamily and selected-genus preview maps.
4. For WCVP panels, a cell is coloured only when both conditions are met: the focal subfamily/genus has WCVP native richness > 0, and the updated Figure 3 global mask retains that cell.

## Important caveat

The Kingsley annotations were supplied as screenshot marks rather than exact GIS coordinates. Therefore, the removals were implemented by nearest-coloured-cell matching. The exact matched cells are listed in `kingsley_manual_annotation_matched_cells.csv` and `kingsley_manual_removed_cells.csv` for manual checking.

## Main outputs

- `Figure3_confidence_priority_STRICT_Kingsley20260526.png`
- `Figure3_confidence_priority_STRICT_Kingsley20260526_with_removed_marks.png`
- `Figure3_STRICT_before_vs_after_Kingsley20260526.png`
- `WCVP_subfamily_richness_5panel_Figure3KingsleyMask.png`
- `WCVP_selected_genus_richness_12panel_Figure3KingsleyMask.png`
- `WCVP_subfamily_priority_preview_5panel_Figure3KingsleyMask.png`
- `WCVP_selected_genus_priority_preview_12panel_Figure3KingsleyMask.png`
