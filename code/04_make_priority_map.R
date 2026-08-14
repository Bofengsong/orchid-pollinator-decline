suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(readxl)
  library(readr)
  library(janitor)
  library(stringr)
  library(patchwork)
})

# ============================================================
# Figure 3 test:
# Replace the environmental-risk layer with:
#   A) predicted risk
#   B) confidence of high-risk classification
#   C) confidence-weighted risk = predicted risk * confidence
#
# Ecological recommendation:
#   C is the most defensible replacement because it prioritises areas
#   that are both high-risk and consistently classified as high-risk.
# ============================================================

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
outdir <- file.path(root, "outputs/priority_maps")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

f_kew_rds <- file.path(root, "data/spatial_inputs/global_orchid_hotspot/Orchids_summary_sf.rds")
f_main <- file.path(root, "data/current/MainData.xlsx")
f_pred_risk <- file.path(root, "data/spatial_inputs/risk_layers",
                         "Predicted_orchid_pollinator_risk_1deg.tif")
f_conf <- file.path(root, "data/spatial_inputs/risk_layers",
                    "Confidence_high_risk_1deg.tif")

stopifnot(file.exists(f_kew_rds), file.exists(f_main), file.exists(f_pred_risk), file.exists(f_conf))

# ------------------------------------------------------------
# 1) 2° grid
# ------------------------------------------------------------
bbox <- st_as_sfc(
  st_bbox(c(xmin = -180, xmax = 180, ymin = -60, ymax = 85), crs = 4326)
)

grid2 <- st_sf(geometry = st_make_grid(bbox, cellsize = c(2, 2), square = TRUE)) %>%
  mutate(
    index = row_number(),
    centroid = st_centroid(geometry),
    lon_bin = st_coordinates(centroid)[, 1],
    lat_bin = st_coordinates(centroid)[, 2]
  ) %>%
  st_set_crs(4326) %>%
  select(index, lon_bin, lat_bin, geometry)

# ------------------------------------------------------------
# 2) Kew/WCVP orchid hotspot layer
# ------------------------------------------------------------
dat_kew <- readRDS(f_kew_rds) %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  mutate(ml_index = log1p(as.numeric(total)) * as.numeric(pc_threat_pred))

grid_cent <- st_as_sf(
  grid2 %>% st_drop_geometry() %>% select(index, lon_bin, lat_bin),
  coords = c("lon_bin", "lat_bin"),
  crs = 4326
)

grid_join <- st_join(grid_cent, dat_kew, join = st_within, left = TRUE) %>%
  st_drop_geometry()

hot_grid_ml <- grid2 %>%
  left_join(
    grid_join %>% select(index, total, pc_threat_pred, ml_index),
    by = "index"
  ) %>%
  st_drop_geometry() %>%
  transmute(
    index,
    lon_bin,
    lat_bin,
    hot_ml = ifelse(is.na(ml_index), 0, ml_index)
  )

# ------------------------------------------------------------
# 3) Evidence layer from corrected evidence set
# ------------------------------------------------------------
sheets <- excel_sheets(f_main)
sheet_use <- if ("MainData" %in% sheets) "MainData" else sheets[1]

dat_all <- read_excel(f_main, sheet = sheet_use) %>%
  clean_names() %>%
  mutate(
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    trend = str_to_lower(str_squish(as.character(trend)))
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

pts_sf <- st_as_sf(dat_all, coords = c("longitude", "latitude"), crs = 4326)
pts_on_grid <- st_join(pts_sf, grid2, join = st_intersects) %>% st_drop_geometry()

evidence_grid <- pts_on_grid %>%
  distinct(reference_no, index, .keep_all = TRUE) %>%
  group_by(index) %>%
  summarise(evidence = n(), .groups = "drop")

# ------------------------------------------------------------
# 4) Research mismatch layer
# ------------------------------------------------------------
mismatch_map <- function(hot_grid, hot_col, evi_grid = evidence_grid,
                         sigma_breaks = c(-2, -1, 1, 2),
                         land_mask = TRUE) {
  combo <- grid2 %>%
    left_join(hot_grid %>% select(index, !!sym(hot_col)), by = "index") %>%
    left_join(evi_grid, by = "index") %>%
    mutate(
      hot = coalesce(.data[[hot_col]], 0),
      evidence = coalesce(evidence, 0)
    ) %>%
    filter(!(hot == 0 & evidence == 0)) %>%
    mutate(
      z_hot = scale(log1p(hot))[, 1],
      z_evi = scale(log1p(evidence))[, 1],
      M = z_hot - z_evi
    )
  
  if (land_mask) {
    world_poly <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
      st_transform(4326) %>%
      st_make_valid() %>%
      st_union()
    combo <- combo[lengths(st_intersects(combo, world_poly)) > 0, ]
  }
  
  combo %>%
    mutate(
      bin = cut(
        M,
        breaks = c(-Inf, sigma_breaks[1], sigma_breaks[2], sigma_breaks[3], sigma_breaks[4], Inf),
        labels = c("<= -2σ (over-studied)", "-2σ~ -1σ", "-1σ~ +1σ",
                   "+1σ~ +2σ", ">= +2σ (under-studied)"),
        include.lowest = TRUE,
        right = TRUE,
        ordered_result = TRUE
      )
    )
}

gap_df <- mismatch_map(hot_grid_ml, "hot_ml") %>%
  st_drop_geometry() %>%
  select(index, M, bin) %>%
  mutate(
    gap3 = case_when(
      bin %in% c(">= +2σ (under-studied)", "+1σ~ +2σ") ~ "Under",
      bin %in% c("<= -2σ (over-studied)", "-2σ~ -1σ") ~ "Over",
      TRUE ~ "Balanced"
    )
  )

# ------------------------------------------------------------
# 5) Risk variants aggregated to 2°
# ------------------------------------------------------------
r_pred <- rast(f_pred_risk)
r_conf <- rast(f_conf)
r_weighted <- r_pred * r_conf
names(r_pred) <- "predicted_risk"
names(r_conf) <- "confidence_high_risk"
names(r_weighted) <- "confidence_weighted_risk"

extract_mean <- function(r, name) {
  x <- terra::extract(r, vect(grid2), fun = mean, na.rm = TRUE)
  names(x)[1:2] <- c("index", name)
  x[[name]][is.na(x[[name]])] <- 0
  x
}

risk_df <- extract_mean(r_pred, "predicted_risk") %>%
  left_join(extract_mean(r_conf, "confidence_high_risk"), by = "index") %>%
  left_join(extract_mean(r_weighted, "confidence_weighted_risk"), by = "index")

write_csv(risk_df, file.path(outdir, "risk_variants_2deg.csv"))

# ------------------------------------------------------------
# 6) Build priority map using selected risk column
# ------------------------------------------------------------
build_priority <- function(risk_col) {
  hot_vals <- hot_grid_ml$hot_ml[hot_grid_ml$hot_ml > 0]
  hot_q <- quantile(hot_vals, probs = c(0.5, 0.75), na.rm = TRUE)
  names(hot_q) <- c("med", "high")
  
  risk_vals <- risk_df[[risk_col]][risk_df[[risk_col]] > 0]
  risk_q <- quantile(risk_vals, probs = c(1/3, 2/3), na.rm = TRUE)
  names(risk_q) <- c("low_high", "med_high")
  
  grid2 %>%
    left_join(hot_grid_ml, by = "index") %>%
    left_join(evidence_grid, by = "index") %>%
    left_join(gap_df %>% select(index, gap3), by = "index") %>%
    left_join(risk_df %>% select(index, all_of(risk_col)), by = "index") %>%
    mutate(
      hot = coalesce(hot_ml, 0),
      evid = coalesce(evidence, 0),
      risk = coalesce(.data[[risk_col]], 0),
      B_class = case_when(
        hot >= hot_q["high"] ~ "High",
        hot >= hot_q["med"] ~ "Medium",
        hot > 0 ~ "Low",
        TRUE ~ "None"
      ),
      R_class = case_when(
        risk >= risk_q["med_high"] ~ "High",
        risk >= risk_q["low_high"] ~ "Medium",
        risk > 0 ~ "Low",
        TRUE ~ "None"
      ),
      gap3 = recode(gap3, .missing = "Under"),
      BasePriority = case_when(
        B_class == "High" & R_class == "High" ~ "Very High",
        B_class == "High" & R_class == "Medium" ~ "High",
        B_class == "High" & R_class == "Low" ~ "Medium",
        B_class == "High" & R_class == "None" ~ "Low",
        B_class == "Medium" & R_class == "High" ~ "High",
        B_class == "Medium" & R_class == "Medium" ~ "Medium",
        B_class == "Medium" & R_class == "Low" ~ "Low",
        B_class == "Medium" & R_class == "None" ~ "Low",
        B_class == "Low" & R_class == "High" ~ "Medium",
        B_class == "Low" & R_class == "Medium" ~ "Low",
        B_class == "Low" & R_class == "Low" ~ "Very Low",
        B_class == "Low" & R_class == "None" ~ "Very Low",
        B_class == "None" & R_class %in% c("High", "Medium") ~ "Low",
        B_class == "None" & R_class %in% c("Low", "None") ~ "Very Low",
        TRUE ~ "Very Low"
      ),
      BaseRank = case_when(
        BasePriority == "Very Low" ~ 1L,
        BasePriority == "Low" ~ 2L,
        BasePriority == "Medium" ~ 3L,
        BasePriority == "High" ~ 4L,
        BasePriority == "Very High" ~ 5L,
        TRUE ~ 1L
      ),
      AdjRank = case_when(
        gap3 == "Under" ~ pmin(BaseRank + 1L, 5L),
        gap3 == "Over" ~ pmax(BaseRank - 1L, 1L),
        TRUE ~ BaseRank
      ),
      AdjRank = case_when(
        gap3 == "Under" & risk > 0 & AdjRank < 2L ~ 2L,
        TRUE ~ AdjRank
      ),
      Priority_chr = case_when(
        AdjRank == 1L ~ "Very Low",
        AdjRank == 2L ~ "Low",
        AdjRank == 3L ~ "Medium",
        AdjRank == 4L ~ "High",
        AdjRank == 5L ~ "Very High",
        TRUE ~ "Very Low"
      ),
      Priority_chr = ifelse(hot == 0 & risk == 0 & evid == 0, NA_character_, Priority_chr),
      Priority = factor(Priority_chr, levels = c("Very Low", "Low", "Medium", "High", "Very High")),
      risk_layer = risk_col
    )
}

priority_pred <- build_priority("predicted_risk")
priority_conf <- build_priority("confidence_high_risk")
priority_weighted <- build_priority("confidence_weighted_risk")

priority_all <- bind_rows(
  priority_pred %>% st_drop_geometry(),
  priority_conf %>% st_drop_geometry(),
  priority_weighted %>% st_drop_geometry()
)
write_csv(priority_all, file.path(outdir, "priority_classes_all_variants_2deg.csv"))

# ------------------------------------------------------------
# 7) Plot
# ------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326)

cols_priority <- c(
  "Very Low" = "#bdbdbd",
  "Low" = "#cbddf0",
  "Medium" = "#9ecae1",
  "High" = "#fc9272",
  "Very High" = "#de2d26"
)

plot_priority <- function(sf_obj, title) {
  ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
    geom_sf(data = sf_obj, aes(fill = Priority), color = NA) +
    scale_fill_manual(values = cols_priority, drop = FALSE, name = "Priority class",
                      na.value = "transparent") +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

p1 <- plot_priority(priority_pred, "A. Priority using predicted environmental risk")
p2 <- plot_priority(priority_conf, "B. Priority using confidence as risk layer")
p3 <- plot_priority(priority_weighted, "C. Priority using confidence-weighted risk")

ggsave(file.path(outdir, "Priority_predicted_risk_2deg.png"),
       p1, width = 12, height = 6, dpi = 300)
ggsave(file.path(outdir, "Priority_confidence_as_risk_2deg.png"),
       p2, width = 12, height = 6, dpi = 300)
ggsave(file.path(outdir, "Priority_confidence_weighted_risk_2deg.png"),
       p3, width = 12, height = 6, dpi = 300)

p_combined <- p1 / p2 / p3 + plot_layout(guides = "collect")
ggsave(file.path(outdir, "Priority_risk_variant_comparison_2deg.png"),
       p_combined, width = 12, height = 15, dpi = 300)

cat("\nSaved priority risk-variant maps to:\n", outdir, "\n", sep = "")
