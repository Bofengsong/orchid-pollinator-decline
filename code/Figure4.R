# Figure 4：五个兰科亚科的保护优先级。
local({
# 在 RStudio 打开本脚本，点击 Source；也可以在仓库根目录用 Rscript 运行。
paths <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
paths <- c(rev(paths), sub("^--file=", "", args))
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  paths <- c(paths, tryCatch(rstudioapi::getSourceEditorContext()$path,
                            error = function(e) ""))
}
paths <- paths[!is.na(paths) & nzchar(paths)]
dirs <- dirname(normalizePath(gsub("~+~", " ", paths, fixed = TRUE), mustWork = FALSE))
candidates <- unique(c(dirs, dirname(dirs), getwd(), dirname(getwd()),
                       getOption("chapter1.project_folder")))
has_data <- function(p) all(file.exists(file.path(p, "data",
                                      c("Final_data.csv", "Supporting_data.rds"))))
found <- candidates[vapply(candidates, has_data, logical(1))]
if (!length(found) && interactive()) {
  message("请选择仓库 data 文件夹中的 Final_data.csv。")
  selected <- file.choose()
  candidate <- dirname(dirname(normalizePath(selected)))
  if (has_data(candidate)) found <- candidate
}
if (!length(found)) stop("未找到 data/Final_data.csv 和 data/Supporting_data.rds，请完整下载仓库。")
root <- normalizePath(found[1])
data_dir <- file.path(root, "data")
output_dir <- file.path(root, "outputs")
dir.create(output_dir, showWarnings = FALSE)
options(chapter1.project_folder = root)
old_wd <- setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)
# Mac 使用自带 Quartz；避免未安装 XQuartz 时的 Cairo 报错。
if (Sys.info()[["sysname"]] == "Darwin") options(bitmapType = "quartz")
# 字体布局在临时设备中完成，避免额外生成 Rplots.pdf。
grDevices::pdf(file = NULL)
layout_device <- grDevices::dev.cur()
on.exit({
    if (!is.na(layout_device) && layout_device %in% grDevices::dev.list()) grDevices::dev.off(layout_device)
}, add = TRUE)
support <- readRDS(file.path(data_dir, "Supporting_data.rds"))
if (!file.exists("Results.rds")) stop("请先运行 Analysis.R。")
results <- readRDS("Results.rds")
if (!identical(unname(tools::md5sum(file.path(data_dir, "Final_data.csv"))), results$input_md5)) {
  stop("Final_data.csv 已变化，请先重新运行 Analysis.R，再画图。")
}
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(sf)
    library(patchwork)
})

subfamilies <- c(
  "Epidendroideae",
  "Orchidoideae",
  "Cypripedioideae",
  "Vanilloideae",
  "Apostasioideae"
)

priority_levels <- c("Very Low", "Low", "Medium", "High", "Very High")
cols_priority <- c(
  "Very Low" = "#9a99a0",
  "Low" = "#b7d3e6",
  "Medium" = "#f6c49a",
  "High" = "#fe8c3a",
  "Very High" = "#d44c02"
)

plate_crs <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"

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

all_grid <- results$priority$lineages %>%
mutate(
  Priority_final_chr = ifelse(is.finite(rank_masked), priority_levels[rank_masked], NA_character_),

  group = as.character(group),
  group_type = as.character(group_type),
  Priority_final = factor(Priority_final_chr, levels = priority_levels)
)

priority_sf <- all_grid %>%
filter(
  group %in% subfamilies,
  group_type == "Subfamily",
  !is.na(Priority_final_chr)
) %>%
mutate(
  group = factor(group, levels = subfamilies),
  Priority_final = factor(Priority_final_chr, levels = priority_levels)
) %>%
make_cell_polygons() %>%
st_transform(plate_crs)

world <- support$world %>%
st_transform(plate_crs) %>%
st_make_valid()
coast <- sf::st_boundary(sf::st_union(support$world)) %>%
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
  plot.title = element_text(
  face = "bold", size = 10, hjust = 0.5,
  margin = margin(b = 0)
),
  legend.position = "bottom",
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8),
  plot.margin = margin(0.2, 0.2, 0.2, 0.2),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)

make_priority_panel <- function(dat_sf, taxon, panel_letter) {
  taxon_sf <- dat_sf %>% filter(group == taxon)

  ggplot() +
  base_panel() +
  geom_sf(data = taxon_sf, aes(fill = Priority_final), color = NA) +
  geom_sf(data = frame, fill = NA, color = "black", linewidth = 0.35) +
  annotate(
    "text",
    x = label_x, y = label_y, label = paste0("(", panel_letter, ")"),
    hjust = 0, vjust = 1, size = 3.1, family = "Times"
  ) +
  scale_fill_manual(
    values = cols_priority,
    drop = FALSE,
    name = "Priority class",
    na.value = "transparent"
  ) +
  coord_sf(crs = plate_crs, datum = NA, expand = FALSE) +
  labs(title = taxon) +
  panel_theme
}

icon_files <- c(Epidendroideae="Epidendroideae_Cattleya_labiata_icon.png",Orchidoideae="Orchidoideae_Orchis_mascula_icon.png",
  Cypripedioideae="Cypripedioideae_Cypripedium_calceolus_icon.png",Vanilloideae="Vanilloideae_Vanilla_planifolia_icon.png",
  Apostasioideae="Apostasioideae_Apostasia_wallichii_icon.png")
p1 <- make_priority_panel(priority_sf, "Epidendroideae", "a")
p2 <- make_priority_panel(priority_sf, "Orchidoideae", "b")
p3 <- make_priority_panel(priority_sf, "Cypripedioideae", "c")
p4 <- make_priority_panel(priority_sf, "Vanilloideae", "d")
p5 <- make_priority_panel(priority_sf, "Apostasioideae", "e")

subfamily_design <- c(
  patchwork::area(t = 1, l = 4, b = 1, r = 7),
  patchwork::area(t = 3, l = 2, b = 3, r = 5),
  patchwork::area(t = 3, l = 6, b = 3, r = 9),
  patchwork::area(t = 5, l = 2, b = 5, r = 5),
  patchwork::area(t = 5, l = 6, b = 5, r = 9)
)

subfamily_5panel <- (p1 + p2 + p3 + p4 + p5) +
plot_layout(
  design = subfamily_design,
  guides = "collect",
  heights = c(1, 0.001, 1, 0.001, 1),
  widths = rep(1, 10)
) &
theme(
  legend.position = "bottom",
  legend.key.height = unit(0.35, "cm"),
  legend.key.width = unit(0.65, "cm"),
  legend.box.margin = margin(t = -2),
  plot.margin = margin(4, 0.5, 4, 0.5)
)

# Add the original pictograms to the title grobs; keep the five-map vector layout intact.
subfamily_5panel <- patchwork::patchworkGrob(subfamily_5panel)
for(i in seq_along(subfamilies)) {
  idx<-which(subfamily_5panel$layout$name==paste0("title-",i))
  if(length(idx)!=1) stop("Missing title grob for panel",i)
  taxon<-subfamilies[i]
  im<-support$icons[[icon_files[[taxon]]]]
  alpha<-if(length(dim(im))==3&&dim(im)[3]==4) im[,,4] else matrix(1,nrow(im),ncol(im))
  rr<-which(rowSums(alpha)>0)
  cc<-which(colSums(alpha)>0)
  im<-im[min(rr):max(rr),min(cc):max(cc),,drop=FALSE]
  tt<-grid::textGrob(taxon,x=.5,y=.5,gp=grid::gpar(fontfamily="Times",fontsize=10,fontface="bold"))
  icon<-grid::rasterGrob(im,x=grid::unit(.5,"npc")-grid::grobWidth(tt)/2-grid::unit(.15,"in"),y=.5,width=grid::unit(.24,
    "in"),height=grid::unit(.24,"in"))
  subfamily_5panel$grobs[[idx]]<-grid::grobTree(tt,icon)
}

ggsave(
  "Figure4.png",
  subfamily_5panel,
  width = 10.0,
  height = 7.4,
  dpi = 300,
  bg = "white"
)


# 保存后回到可见的图形设备，在 Plots / Quartz 中预览。
if (layout_device %in% grDevices::dev.list()) grDevices::dev.off(layout_device)
layout_device <- NA_integer_
if (interactive()) {
  visible_devices <- c("RStudioGD", "quartz", "windows", "X11", "X11cairo")
  preview_ready <- TRUE
  if (!names(grDevices::dev.cur()) %in% visible_devices) {
    preview_ready <- tryCatch({
        if (Sys.info()[["sysname"]] == "Darwin" && Sys.getenv("RSTUDIO") != "1") {
          suppressWarnings(grDevices::quartz())
      } else {
          grDevices::dev.new()
      }
        TRUE
    }, error = function(e) {
        message("PNG已保存。当前R会话没有可用图形窗口；在RStudio或R.app中运行可预览。")
        FALSE
    })
  }
  if (preview_ready) {
    # 预览按整张原图缩放，避免小Plots面板把森林图文字挤在一起。
    preview_image <- png::readPNG("Figure4.png", native = TRUE)
    image_ratio <- nrow(preview_image) / ncol(preview_image)
    preview_scale <- min(1, 1600 / max(dim(preview_image)))
    if (preview_scale < 1) {
      rows <- round(seq(1, nrow(preview_image), length.out = round(nrow(preview_image) * preview_scale)))
      cols <- round(seq(1, ncol(preview_image), length.out = round(ncol(preview_image) * preview_scale)))
      # nativeRaster按行储存像素，不能直接用普通矩阵的行列下标缩放。
      pixels <- rep((rows - 1L) * ncol(preview_image), each = length(cols)) + rep(cols, times = length(rows))
      preview_image <- as.integer(preview_image)[pixels]
      dim(preview_image) <- c(length(rows), length(cols))
      class(preview_image) <- "nativeRaster"
      attr(preview_image, "channels") <- 4L
    }
    preview_plot <- ggplot2::ggplot() +
    ggplot2::annotation_raster(preview_image, 0, 1, 0, 1, interpolate = TRUE) +
    ggplot2::coord_fixed(ratio = image_ratio, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    print(preview_plot)
  }
}
})
