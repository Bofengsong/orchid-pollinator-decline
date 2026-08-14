# Final WCVP Subfamily and Selected-Genus Maps

Generated on: 2026-05-27

## Data logic

These maps use Kew/WCVP accepted native/extant/non-doubtful Orchidaceae species distributions aggregated to TDWG botanical-country units, converted to the 2-degree Figure 3 grid and masked by the latest strict Figure 3 orchid-presence mask.

The maps are exploratory lineage-specific previews. They do not replace the main Figure 3 hotspot layer because WCVP provides richness/distribution, not lineage-specific threatened proportions.

## Included groups

Subfamilies: Epidendroideae, Orchidoideae, Cypripedioideae, Vanilloideae, Apostasioideae.

Additional tribe-level panel: Diurideae, defined from APWeb code `Odiu` in the genus crosswalk. This is included as an exploratory lineage-focused panel within Orchidoideae, not as a sixth subfamily.

Selected genera: Cypripedium, Paphiopedilum, Phragmipedium, Cymbidium, Orchis, Corybas, Platanthera, Vanilla, Caladenia, Disa, Bulbophyllum, Dendrobium.

## Outputs

- `WCVP_subfamily_plus_Diurideae_priority_6panel_ROBINSON_3x2.png`
- `WCVP_subfamily_plus_Diurideae_relative_richness_6panel_ROBINSON_3x2.png`
- `WCVP_selected_genus_priority_12panel_ROBINSON_3x4.png`
- `WCVP_selected_genus_relative_richness_12panel_ROBINSON_3x4.png`
- `WCVP_subfamily_genus_final_masked_2deg.csv`
- `WCVP_final_maps_summary.csv`
- `build_final_wcvp_maps.R`

## Interpretation

Priority maps reuse the Figure 3 logic, replacing the overall orchid hotspot component with focal subfamily/genus WCVP richness. Richness maps show within-taxon relative richness after the same final Figure 3 global mask.
