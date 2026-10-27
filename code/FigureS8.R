# FigureS8：沿用已固定的环境及兰花分布资料。
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
suppressPackageStartupMessages({library(terra)
    library(sf)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)})
world <- support$world
save_all <- function(p,name,w,h) ggsave(paste0(name,'.png'),p,width=w,height=h,dpi=300,bg='white')
bbox<-st_as_sfc(st_bbox(c(xmin=-180,xmax=180,ymin=-60,ymax=85),crs=4326))
g<-st_sf(geometry=st_make_grid(bbox,cellsize=c(2,2),square=TRUE))
g$index<-seq_len(nrow(g))
xy<-st_coordinates(st_centroid(g))
g$lon_bin<-xy[,1]
g$lat_bin<-xy[,2]
kew<-support$kew |> st_make_valid() |> st_transform(4326) |> mutate(ml_index=log1p(as.numeric(total))*as.numeric(pc_threat_pred))
cent<-st_as_sf(st_drop_geometry(g),coords=c('lon_bin','lat_bin'),crs=4326)
j<-st_join(cent,kew,join=st_within,left=TRUE) |> st_drop_geometry() |> select(index,ml_index)
g<-left_join(g,j,by='index') |> mutate(hot_ml=ifelse(is.na(ml_index),0,ml_index))
br<-quantile(g$hot_ml[g$hot_ml>0],c(.5,.75,.9),na.rm=TRUE)
g$hotspot_class<-cut(g$hot_ml,c(-Inf,0,br,Inf),labels=c('No hotspot signal','Low','Medium','High','Very high'),
  include.lowest=TRUE)
p8<-ggplot()+geom_sf(data=world,fill='grey96',color='grey80',linewidth=.15)+geom_sf(data=g,aes(fill=hotspot_class),
  color=NA)+
scale_fill_manual(values=c('No hotspot signal'='#f0f0f0','Low'='#d9f0d3','Medium'='#addd8e','High'='#31a354',
  'Very high'='#006d2c'),drop=FALSE,name='Hotspot class')+
coord_sf(xlim=c(-180,180),ylim=c(-60,85),expand=FALSE)+labs(title='Global orchid hotspot layer (2° grid)')+theme_void(base_size=12)+theme(plot.title=element_text(hjust=.5,
  face='bold'),legend.position='bottom',plot.background=element_rect(fill='white',color=NA))
save_all(p8,'FigureS8',12,6)

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
    preview_image <- png::readPNG("FigureS8.png", native = TRUE)
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
