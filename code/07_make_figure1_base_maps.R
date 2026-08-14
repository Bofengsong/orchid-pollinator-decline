suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(rnaturalearth)
  library(metafor)
  library(cowplot)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
outdir <- file.path(root, "outputs/figure1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

data_dir <- file.path(root, "data/current")
meta_results <- file.path(root, "outputs", "meta_analysis")

main_path <- file.path(data_dir, "MainData.xlsx")
meta_path <- file.path(data_dir, "MetaData.csv")
hotspot_rds <- file.path(root, "data/spatial_inputs/global_orchid_hotspot", "Orchids_summary_sf.rds")

stopifnot(file.exists(main_path), file.exists(meta_path), file.exists(hotspot_rds))

theme_map <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.margin = margin(0, 8, 0, 8)
    )
}

cols_order <- c(
  "Hymenoptera" = "#B2D7B0",
  "Diptera"     = "#F7B7D2",
  "Coleoptera"  = "#8CC7DD",
  "Lepidoptera" = "#F6D993"
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  filter(name != "Antarctica") |>
  st_transform(4326)

base_map_layers <- list(
  geom_sf(data = world, fill = "grey97", color = "grey78", linewidth = 0.18),
  coord_sf(xlim = c(-180, 180), ylim = c(-48, 85), expand = FALSE)
)

read_main_data <- function() {
  sheet_use <- if ("MainData" %in% excel_sheets(main_path)) "MainData" else excel_sheets(main_path)[1]
  read_excel(main_path, sheet = sheet_use) |>
    mutate(
      latitude = suppressWarnings(as.numeric(latitude)),
      longitude = suppressWarnings(as.numeric(longitude)),
      pollinator_order = str_to_title(str_squish(as.character(pollinator_order))),
      reference_no = as.character(reference_no)
    ) |>
    filter(is.finite(latitude), is.finite(longitude), pollinator_order %in% names(cols_order))
}

read_meta_data <- function() {
  read_csv(meta_path, show_col_types = FALSE) |>
    mutate(
      latitude = suppressWarnings(as.numeric(latitude)),
      longitude = suppressWarnings(as.numeric(longitude)),
      pollinator_order = str_to_title(str_squish(as.character(pollinator_order))),
      reference_no = as.character(reference_no),
      country = as.character(country),
      pollinator_family = as.character(pollinator_family),
      yi = suppressWarnings(as.numeric(yi)),
      vi = suppressWarnings(as.numeric(vi))
    ) |>
    filter(is.finite(latitude), is.finite(longitude), pollinator_order %in% names(cols_order))
}

grid_order_counts <- function(dat, count_label = "n") {
  dat |>
    mutate(
      lon = floor(longitude / 2) * 2 + 1,
      lat = floor(latitude / 2) * 2 + 1
    ) |>
    count(pollinator_order, lon, lat, name = count_label) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
}

make_bubble_panel <- function(points, orders, panel_label, size_title, breaks = c(1, 5, 10, 20)) {
  ggplot() +
    base_map_layers +
    annotate("text", x = -179, y = 83, label = panel_label, hjust = 0, vjust = 1,
             fontface = "plain", size = 5.1, color = "black") +
    geom_sf(
      data = points |> filter(pollinator_order %in% orders),
      aes(size = n, fill = pollinator_order),
      shape = 21, color = "grey15", alpha = 0.82, stroke = 0.45
    ) +
    scale_fill_manual(
      values = cols_order, limits = names(cols_order), drop = FALSE,
      name = "Pollinator order"
    ) +
    scale_size_area(
      max_size = 12, breaks = breaks, limits = c(1, max(max(points$n, na.rm = TRUE), max(breaks))),
      name = size_title
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 4), order = 2),
      size = guide_legend(order = 1)
    ) +
    theme_map(11)
}

make_bubble_legend <- function(size_title, breaks = c(1, 5, 10, 20)) {
  size_df <- tibble(
    x = c(4.15, 4.95, 5.95, 7.15),
    y = 1.38,
    n = breaks,
    lab = as.character(breaks)
  )
  fill_df <- tibble(
    x = c(4.15, 6.0, 7.35, 9.0),
    y = 0.66,
    pollinator_order = factor(names(cols_order), levels = names(cols_order))
  )

  ggplot() +
    annotate("text", x = 0.35, y = 1.38, label = size_title, hjust = 0, size = 4.1) +
    geom_point(data = size_df, aes(x = x, y = y, size = n),
               shape = 21, fill = "white", color = "grey20", stroke = 0.45) +
    geom_text(data = size_df, aes(x = x + 0.32, y = y, label = lab),
              hjust = 0, vjust = 0.5, size = 3.8) +
    annotate("text", x = 0.35, y = 0.66, label = "Pollinator order", hjust = 0, size = 4.1) +
    geom_point(data = fill_df, aes(x = x, y = y, fill = pollinator_order),
               shape = 21, size = 4.7, color = "grey20", stroke = 0.45) +
    geom_text(data = fill_df, aes(x = x + 0.28, y = y, label = pollinator_order),
              hjust = 0, vjust = 0.5, size = 3.8) +
    scale_size_area(max_size = 11, limits = c(1, max(breaks)), guide = "none") +
    scale_fill_manual(values = cols_order, guide = "none") +
    coord_cartesian(xlim = c(0, 11.6), ylim = c(0.35, 1.68), clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

save_plot <- function(plot, name, width, height, dpi = 300) {
  ggsave(file.path(outdir, paste0(name, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

# -------------------------------------------------------------------------
# Figure 1: corrected evidence-set bubble map
# -------------------------------------------------------------------------
main_dat <- read_main_data()
main_pts <- grid_order_counts(main_dat)

fig1_a <- make_bubble_panel(
  main_pts, c("Hymenoptera", "Diptera"),
  "(a)",
  "Number of species observations",
  breaks = c(1, 5, 10, 20)
)
fig1_b <- make_bubble_panel(
  main_pts, c("Coleoptera", "Lepidoptera"),
  "(b)",
  "Number of species observations",
  breaks = c(1, 5, 10, 20)
)

figure1 <- wrap_plots(
  fig1_a, fig1_b, make_bubble_legend("Number of species observations"),
  ncol = 1, heights = c(1, 1, 0.22)
)
save_plot(figure1, "Figure 1", width = 10, height = 8.4)

# -------------------------------------------------------------------------
# Figure S2: corrected meta-analysis bubble map, matched to Figure 1 style
# -------------------------------------------------------------------------
meta_dat <- read_meta_data()
meta_pts <- grid_order_counts(meta_dat)

figS2_a <- make_bubble_panel(
  meta_pts, c("Hymenoptera", "Diptera"),
  "(a)",
  "Meta-analysis observations",
  breaks = c(1, 5, 10, 20)
)
figS2_b <- make_bubble_panel(
  meta_pts, c("Coleoptera", "Lepidoptera"),
  "(b)",
  "Meta-analysis observations",
  breaks = c(1, 5, 10, 20)
)

figureS2 <- wrap_plots(
  figS2_a, figS2_b, make_bubble_legend("Meta-analysis observations"),
  ncol = 1, heights = c(1, 1, 0.22)
)
save_plot(figureS2, "Figure S2", width = 10, height = 8.4)

# -------------------------------------------------------------------------
# Figure S3: corrected study-method subgroup estimates
# -------------------------------------------------------------------------
method <- read_csv(file.path(meta_results, "Part2B_subgroup_study_method_corrected.csv"), show_col_types = FALSE) |>
  mutate(
    level = factor(level, levels = rev(c("Historical", "Historical + Field", "Field"))),
    p_lab = ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p))),
    label = paste0("k = ", k, ", ", p_lab)
  )

x_min <- floor(min(method$ci_lb, na.rm = TRUE) * 2) / 2
x_max <- ceiling(max(method$ci_ub, na.rm = TRUE) * 2) / 2

figureS3 <- ggplot(method, aes(x = yi, y = level)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.6) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), orientation = "y",
                width = 0, linewidth = 1.15, color = "grey35") +
  geom_point(shape = 21, size = 5.5, fill = "#9ecae1", color = "black", stroke = 0.65) +
  geom_text(aes(label = label), x = x_max + 0.35, hjust = 0, size = 4.1) +
  coord_cartesian(xlim = c(x_min - 0.1, x_max + 1.45), clip = "off") +
  scale_x_continuous(breaks = seq(floor(x_min), ceiling(x_max), by = 1)) +
  labs(x = "Effect size (lnRR)", y = NULL, title = "Study-method subgroup meta-analysis") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", size = 15),
    plot.margin = margin(10, 60, 10, 10)
  )
save_plot(figureS3, "Figure S3", width = 8.8, height = 4.8)

# -------------------------------------------------------------------------
# Figure S9: corrected hotspot-evidence mismatch map
# -------------------------------------------------------------------------
bbox <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -60, ymax = 85), crs = 4326))
grid2 <- st_sf(geometry = st_make_grid(bbox, cellsize = c(2, 2), square = TRUE)) |>
  mutate(
    index = row_number(),
    centroid = st_centroid(geometry),
    lon_bin = st_coordinates(centroid)[, 1],
    lat_bin = st_coordinates(centroid)[, 2]
  ) |>
  st_set_crs(4326) |>
  select(index, lon_bin, lat_bin, geometry)

kew <- readRDS(hotspot_rds) |>
  st_make_valid() |>
  st_transform(4326) |>
  mutate(ml_index = log1p(as.numeric(total)) * as.numeric(pc_threat_pred))

grid_cent <- st_as_sf(
  grid2 |> st_drop_geometry() |> select(index, lon_bin, lat_bin),
  coords = c("lon_bin", "lat_bin"), crs = 4326
)

hot_join <- st_join(grid_cent, kew, join = st_within, left = TRUE) |>
  st_drop_geometry()

hot_grid <- grid2 |>
  st_drop_geometry() |>
  left_join(hot_join |> select(index, ml_index), by = "index") |>
  mutate(hot_ml = ifelse(is.na(ml_index), 0, ml_index)) |>
  select(index, hot_ml)

main_sf <- st_as_sf(main_dat, coords = c("longitude", "latitude"), crs = 4326)
main_grid <- st_join(main_sf, grid2, join = st_intersects) |>
  st_drop_geometry()

evidence_grid <- main_grid |>
  distinct(reference_no, index, .keep_all = TRUE) |>
  count(index, name = "evidence")

land <- world |>
  st_make_valid() |>
  st_union()

mismatch <- grid2 |>
  left_join(hot_grid, by = "index") |>
  left_join(evidence_grid, by = "index") |>
  mutate(
    hot_ml = coalesce(hot_ml, 0),
    evidence = coalesce(evidence, 0)
  ) |>
  filter(!(hot_ml == 0 & evidence == 0)) |>
  mutate(
    z_hotspot = as.numeric(scale(log1p(hot_ml))),
    z_evidence = as.numeric(scale(log1p(evidence))),
    M = z_hotspot - z_evidence,
    mismatch_class = cut(
      M,
      breaks = c(-Inf, -2, -1, 1, 2, Inf),
      labels = c(
        "Severely over-studied",
        "Moderately over-studied",
        "Balanced",
        "Moderately under-studied",
        "Severely under-studied"
      ),
      include.lowest = TRUE
    )
  )
mismatch <- mismatch[lengths(st_intersects(mismatch, land)) > 0, ]

mismatch_cols <- c(
  "Severely over-studied" = "#2166ac",
  "Moderately over-studied" = "#67a9cf",
  "Balanced" = "#f7f7f7",
  "Moderately under-studied" = "#f4a582",
  "Severely under-studied" = "#b2182b"
)

figureS9 <- ggplot() +
  geom_sf(data = world, fill = "grey97", color = "grey82", linewidth = 0.15) +
  geom_sf(data = mismatch, aes(fill = mismatch_class), color = NA) +
  scale_fill_manual(values = mismatch_cols, drop = FALSE, name = "Mismatch class", na.value = "transparent") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  labs(title = "Spatial mismatch between orchid hotspots and research effort") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
save_plot(figureS9, "Figure S9", width = 12, height = 6)

# -------------------------------------------------------------------------
# Figure S12 and S13: corrected publication-bias and influence diagnostics
# -------------------------------------------------------------------------
clean_continent <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  case_when(
    str_detect(x, "^asia$") ~ "Asia",
    str_detect(x, "^europe$") ~ "Europe",
    str_detect(x, "^north\\s*america$") ~ "North America",
    str_detect(x, "^south\\s*america$") ~ "South America",
    str_detect(x, "oceania|austral") ~ "Oceania",
    TRUE ~ str_to_title(x)
  )
}

meta_model <- read_csv(meta_path, show_col_types = FALSE) |>
  mutate(
    reference_no = as.character(reference_no),
    country = as.character(country),
    pollinator_family = as.character(pollinator_family),
    continent = clean_continent(continent),
    yi = suppressWarnings(as.numeric(yi)),
    vi = suppressWarnings(as.numeric(vi))
  ) |>
  filter(is.finite(yi), is.finite(vi), vi > 0,
         !is.na(reference_no), !is.na(country), !is.na(pollinator_family)) |>
  mutate(
    reference_no = factor(reference_no),
    country = factor(country),
    pollinator_family = factor(pollinator_family)
  )

rand_list <- list(~1 | reference_no, ~1 | country, ~1 | pollinator_family)
m_main <- rma.mv(yi, vi, random = rand_list, data = meta_model, method = "REML")

agg <- meta_model |>
  as_tibble() |>
  group_by(reference_no) |>
  group_modify(~{
    d <- .x
    if (nrow(d) == 1) {
      tibble(yi = d$yi[1], vi = d$vi[1])
    } else {
      fit_i <- rma(yi, vi, data = d, method = "REML")
      tibble(yi = as.numeric(fit_i$b), vi = as.numeric(fit_i$se)^2)
    }
  }) |>
  ungroup()

m_study <- rma(yi, vi, data = agg, method = "REML")
tf_study <- trimfill(m_study)

png(file.path(outdir, "Figure S12.png"), width = 3600, height = 1200, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.1, 3.2, 1.2), oma = c(0, 0, 1.2, 0))
funnel(m_main, xlab = "Effect size (lnRR)", ylab = "Standard Error",
       main = "a  Multilevel model", refline = 0)
funnel(m_study, main = "b  Study-level model", xlab = "lnRR", ylab = "SE")
funnel(tf_study, main = "c  Trim-and-fill", xlab = "lnRR", ylab = "SE")
par(mfrow = c(1, 1))
dev.off()
pdf(file.path(outdir, "Figure S12.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4.2, 4.1, 3.2, 1.2), oma = c(0, 0, 1.2, 0))
funnel(m_main, xlab = "Effect size (lnRR)", ylab = "Standard Error",
       main = "a  Multilevel model", refline = 0)
funnel(m_study, main = "b  Study-level model", xlab = "lnRR", ylab = "SE")
funnel(tf_study, main = "c  Trim-and-fill", xlab = "lnRR", ylab = "SE")
par(mfrow = c(1, 1))
dev.off()

inf <- influence(m_study)
png(file.path(outdir, "Figure S13.png"), width = 1600, height = 1200, res = 200)
plot(inf)
dev.off()
pdf(file.path(outdir, "Figure S13.pdf"), width = 8, height = 6)
plot(inf)
dev.off()

# Keep a small audit table next to the figures.
figure_checks <- tibble(
  item = c("Evidence rows", "Evidence publications", "Meta rows", "Meta publications",
           "Figure S3 Field k", "Figure S9 evidence publications"),
  value = c(nrow(main_dat), n_distinct(main_dat$reference_no), nrow(meta_dat),
            n_distinct(meta_dat$reference_no),
            method$k[method$level == "Field"], n_distinct(main_dat$reference_no))
)
write_csv(figure_checks, file.path(outdir, "updated_figure_source_checks.csv"))

cat("Updated figures saved to:\n", outdir, "\n", sep = "")
