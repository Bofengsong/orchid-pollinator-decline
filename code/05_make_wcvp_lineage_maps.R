suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(patchwork)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source_dir <- file.path(root, "data/spatial_inputs/wcvp_lineage")
outdir <- file.path(root, "outputs/wcvp_lineage")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

subfamilies <- c("Epidendroideae", "Orchidoideae", "Cypripedioideae",
                 "Vanilloideae", "Apostasioideae")
tribes <- c("Diurideae")
selected_genera <- c("Cypripedium", "Paphiopedilum", "Phragmipedium",
                     "Cymbidium", "Orchis", "Corybas",
                     "Platanthera", "Vanilla", "Caladenia",
                     "Disa", "Bulbophyllum", "Dendrobium")

priority_levels <- c("Very Low", "Low", "Medium", "High", "Very High")
cols_priority <- c(
  "Very Low" = "#bdbdbd",
  "Low" = "#cbddf0",
  "Medium" = "#9ecae1",
  "High" = "#fc9272",
  "Very High" = "#de2d26"
)
cols_richness <- c("#f7fbff", "#deebf7", "#9ecae1", "#fdae6b", "#de2d26")

make_cell_polygons <- function(dat) {
  polys <- lapply(seq_len(nrow(dat)), function(i) {
    x <- dat$lon_bin[i]
    y <- dat$lat_bin[i]
    st_polygon(list(matrix(c(
      x - 1, y - 1,
      x + 1, y - 1,
      x + 1, y + 1,
      x - 1, y + 1,
      x - 1, y - 1
    ), ncol = 2, byrow = TRUE)))
  })
  st_sf(dat, geometry = st_sfc(polys, crs = 4326))
}

make_robinson_frame <- function() {
  frame_coords <- rbind(
    cbind(seq(-180, 180, by = 1), -85),
    cbind(180, seq(-85, 85, by = 1)),
    cbind(seq(180, -180, by = -1), 85),
    cbind(-180, seq(85, -85, by = -1)),
    c(-180, -85)
  )
  st_sfc(st_polygon(list(frame_coords)), crs = 4326)
}

plate_crs <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"

all_grid <- read_csv(
  file.path(source_dir, "wcvp_subfamily_genus_priority_preview_2deg_Figure3KingsleyMasked.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    group = as.character(group),
    group_type = as.character(group_type),
    Priority_final = factor(Priority_final_chr, levels = priority_levels)
  )

write_csv(all_grid, file.path(outdir, "WCVP_subfamily_genus_final_masked_2deg.csv"))

summary_table <- all_grid %>%
  filter(group %in% c(subfamilies, tribes, selected_genera)) %>%
  group_by(group_type, group) %>%
  summarise(
    retained_cells = sum(final_keep, na.rm = TRUE),
    very_high_cells = sum(Priority_final_chr == "Very High", na.rm = TRUE),
    high_cells = sum(Priority_final_chr == "High", na.rm = TRUE),
    max_richness_in_retained_cells = suppressWarnings(max(richness_masked, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    max_richness_in_retained_cells = ifelse(is.infinite(max_richness_in_retained_cells),
                                            NA_real_,
                                            max_richness_in_retained_cells)
  )
write_csv(summary_table, file.path(outdir, "WCVP_final_maps_summary.csv"))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(plate_crs) %>%
  st_make_valid()
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_transform(plate_crs) %>%
  st_make_valid()
graticule <- st_graticule(
  lat = seq(-60, 60, by = 30),
  lon = seq(-180, 180, by = 60),
  crs = st_crs(4326)
) %>%
  st_transform(plate_crs)
frame <- make_robinson_frame() %>%
  st_transform(plate_crs)

bb <- st_bbox(frame)
label_x <- bb["xmin"] + 0.09 * (bb["xmax"] - bb["xmin"])
label_y <- bb["ymax"] - 0.30 * (bb["ymax"] - bb["ymin"])

base_panel <- function() {
  list(
    geom_sf(data = frame, fill = "white", color = "black", linewidth = 0.35),
    geom_sf(data = world, fill = "#f1f1f1", color = "grey65", linewidth = 0.08),
    geom_sf(data = graticule, color = "grey67", linewidth = 0.16, linetype = "dotted")
  )
}

panel_theme <- theme_void(base_family = "Times", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5,
                              margin = margin(b = 1)),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(0, 1, 1, 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

make_priority_panel <- function(dat_sf, taxon, panel_letter = NULL) {
  taxon_sf <- dat_sf %>% filter(group == taxon)
  label <- if (is.null(panel_letter)) "" else paste0("(", panel_letter, ")")
  title_face <- if (taxon %in% selected_genera) "bold.italic" else "bold"

  ggplot() +
    base_panel() +
    geom_sf(data = taxon_sf, aes(fill = Priority_final), color = NA) +
    geom_sf(data = coast, fill = NA, color = "grey25", linewidth = 0.11) +
    geom_sf(data = frame, fill = NA, color = "black", linewidth = 0.35) +
    annotate("text", x = label_x, y = label_y, label = label,
             hjust = 0, vjust = 1, size = 3.1, family = "Times") +
    scale_fill_manual(values = cols_priority, drop = FALSE,
                      name = "Priority class", na.value = "transparent") +
    coord_sf(crs = plate_crs, datum = NA, expand = FALSE) +
    labs(title = taxon) +
    panel_theme +
    theme(plot.title = element_text(face = title_face, size = 9.5, hjust = 0.5,
                                    margin = margin(b = 1)))
}

make_richness_panel <- function(dat_sf, taxon, panel_letter = NULL) {
  taxon_sf <- dat_sf %>% filter(group == taxon)
  label <- if (is.null(panel_letter)) "" else paste0("(", panel_letter, ")")
  title_face <- if (taxon %in% selected_genera) "bold.italic" else "bold"

  ggplot() +
    base_panel() +
    geom_sf(data = taxon_sf, aes(fill = richness_relative_masked), color = NA) +
    geom_sf(data = coast, fill = NA, color = "grey25", linewidth = 0.11) +
    geom_sf(data = frame, fill = NA, color = "black", linewidth = 0.35) +
    annotate("text", x = label_x, y = label_y, label = label,
             hjust = 0, vjust = 1, size = 3.1, family = "Times") +
    scale_fill_gradientn(colors = cols_richness, limits = c(0, 1),
                         name = "Relative richness",
                         breaks = c(0, 0.5, 1),
                         labels = c("0", "0.5", "1.0")) +
    coord_sf(crs = plate_crs, datum = NA, expand = FALSE) +
    labs(title = taxon) +
    panel_theme +
    theme(plot.title = element_text(face = title_face, size = 9.5, hjust = 0.5,
                                    margin = margin(b = 1)))
}

priority_sf <- all_grid %>%
  filter(group %in% c(subfamilies, tribes, selected_genera), !is.na(Priority_final_chr)) %>%
  mutate(
    group = factor(group, levels = c(subfamilies, tribes, selected_genera)),
    Priority_final = factor(Priority_final_chr, levels = priority_levels)
  ) %>%
  make_cell_polygons() %>%
  st_transform(plate_crs)

richness_sf <- all_grid %>%
  filter(group %in% c(subfamilies, tribes, selected_genera),
         !is.na(richness_relative_masked),
         richness_relative_masked > 0,
         final_keep) %>%
  mutate(group = factor(group, levels = c(subfamilies, tribes, selected_genera))) %>%
  make_cell_polygons() %>%
  st_transform(plate_crs)

make_subfamily_tribe_layout <- function(dat_sf, type = c("priority", "richness")) {
  type <- match.arg(type)
  panel_fun <- if (type == "priority") make_priority_panel else make_richness_panel

  p1 <- panel_fun(dat_sf, "Epidendroideae", "a")
  p2 <- panel_fun(dat_sf, "Orchidoideae", "b")
  p3 <- panel_fun(dat_sf, "Cypripedioideae", "c")
  p4 <- panel_fun(dat_sf, "Vanilloideae", "d")
  p5 <- panel_fun(dat_sf, "Apostasioideae", "e")
  p6 <- panel_fun(dat_sf, "Diurideae", "f")

  (p1 + p2 + p3) / (p4 + p5 + p6) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.65, "cm"),
      plot.margin = margin(1, 2, 1, 2)
    )
}

make_genus_layout <- function(dat_sf, type = c("priority", "richness")) {
  type <- match.arg(type)
  panel_fun <- if (type == "priority") make_priority_panel else make_richness_panel
  panels <- Map(panel_fun, MoreArgs = list(dat_sf = dat_sf),
                taxon = selected_genera,
                panel_letter = letters[seq_along(selected_genera)])
  wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(ifelse(type == "priority", 0.65, 1.2), "cm"),
      plot.margin = margin(1, 2, 1, 2)
    )
}

subfamily_priority <- make_subfamily_tribe_layout(priority_sf, "priority")
subfamily_richness <- make_subfamily_tribe_layout(richness_sf, "richness")
genus_priority <- make_genus_layout(priority_sf, "priority")
genus_richness <- make_genus_layout(richness_sf, "richness")

ggsave(file.path(outdir, "WCVP_subfamily_plus_Diurideae_priority_6panel_ROBINSON_3x2.png"),
       subfamily_priority, width = 12.2, height = 5.7, dpi = 300, bg = "white")
ggsave(file.path(outdir, "WCVP_subfamily_plus_Diurideae_relative_richness_6panel_ROBINSON_3x2.png"),
       subfamily_richness, width = 12.2, height = 5.7, dpi = 300, bg = "white")
ggsave(file.path(outdir, "WCVP_selected_genus_priority_12panel_ROBINSON_3x4.png"),
       genus_priority, width = 12.2, height = 9.0, dpi = 300, bg = "white")
ggsave(file.path(outdir, "WCVP_selected_genus_relative_richness_12panel_ROBINSON_3x4.png"),
       genus_richness, width = 12.2, height = 9.0, dpi = 300, bg = "white")

readme <- c(
  "# Final WCVP Subfamily and Selected-Genus Maps",
  "",
  paste0("Generated on: ", Sys.Date()),
  "",
  "## Data logic",
  "",
  "These maps use Kew/WCVP accepted native/extant/non-doubtful Orchidaceae species distributions aggregated to TDWG botanical-country units, converted to the 2-degree Figure 3 grid and masked by the latest strict Figure 3 orchid-presence mask.",
  "",
  "The maps are exploratory lineage-specific previews. They do not replace the main Figure 3 hotspot layer because WCVP provides richness/distribution, not lineage-specific threatened proportions.",
  "",
  "## Included groups",
  "",
  "Subfamilies: Epidendroideae, Orchidoideae, Cypripedioideae, Vanilloideae, Apostasioideae.",
  "",
  "Additional tribe-level panel: Diurideae, defined from APWeb code `Odiu` in the genus crosswalk. This is included as an exploratory lineage-focused panel within Orchidoideae, not as a sixth subfamily.",
  "",
  paste0("Selected genera: ", paste(selected_genera, collapse = ", "), "."),
  "",
  "## Outputs",
  "",
  "- `WCVP_subfamily_plus_Diurideae_priority_6panel_ROBINSON_3x2.png`",
  "- `WCVP_subfamily_plus_Diurideae_relative_richness_6panel_ROBINSON_3x2.png`",
  "- `WCVP_selected_genus_priority_12panel_ROBINSON_3x4.png`",
  "- `WCVP_selected_genus_relative_richness_12panel_ROBINSON_3x4.png`",
  "- `WCVP_subfamily_genus_final_masked_2deg.csv`",
  "- `WCVP_final_maps_summary.csv`",
  "- `build_final_wcvp_maps.R`",
  "",
  "## Interpretation",
  "",
  "Priority maps reuse the Figure 3 logic, replacing the overall orchid hotspot component with focal subfamily/genus WCVP richness. Richness maps show within-taxon relative richness after the same final Figure 3 global mask."
)
writeLines(readme, file.path(outdir, "README_WCVP_final_maps.md"))

cat("Final WCVP map package written to:\n", outdir, "\n", sep = "")
