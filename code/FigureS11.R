# FigureS11：沿用已固定的环境及兰花分布资料。
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
selected_genera<-c('Cypripedium','Paphiopedilum','Phragmipedium','Cymbidium','Orchis','Corybas','Platanthera',
  'Vanilla','Caladenia','Disa','Bulbophyllum','Dendrobium')
d<-support$lineage |> filter(group %in% selected_genera,final_keep,!is.na(richness_relative_masked),richness_relative_masked>0)
poly<-lapply(seq_len(nrow(d)),function(i){x<-d$lon_bin[i]
    y<-d$lat_bin[i]
    st_polygon(list(matrix(c(x-1,y-1,x+1,y-1,x+1,y+1,x-1,y+1,x-1,y-1),ncol=2,byrow=TRUE)))})
crs<-'+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs'
ds<-st_sf(d,geometry=st_sfc(poly,crs=4326)) |> st_transform(crs)
coast<-st_boundary(st_union(world)) |> st_transform(crs) |> st_make_valid()
worldr<-st_transform(world,crs) |> st_make_valid()
fc<-rbind(cbind(seq(-180,180,1),-85),cbind(180,seq(-85,85,1)),cbind(seq(180,-180,-1),85),cbind(-180,seq(85,
  -85,-1)),c(-180,-85))
frame<-st_sfc(st_polygon(list(fc)),crs=4326) |> st_transform(crs)
grat<-st_graticule(lat=seq(-60,60,30),lon=seq(-180,180,60),crs=st_crs(4326)) |> st_transform(crs)
bb<-st_bbox(frame)
lx<-bb['xmin']+.09*(bb['xmax']-bb['xmin'])
ly<-bb['ymax']-.30*(bb['ymax']-bb['ymin'])
ps<-lapply(seq_along(selected_genera),function(i){
    ggplot()+geom_sf(data=frame,fill='white',color='black',linewidth=.35)+geom_sf(data=worldr,fill='#f1f1f1',
    color='grey65',linewidth=.08)+geom_sf(data=grat,color='grey67',linewidth=.16,linetype='dotted')+
    geom_sf(data=filter(ds,group==selected_genera[i]),aes(fill=richness_relative_masked),color=NA)+geom_sf(data=coast,
    color='grey25',linewidth=.11)+geom_sf(data=frame,fill=NA,color='black',linewidth=.35)+
    annotate('text',x=lx,y=ly,label=paste0('(',letters[i],')'),hjust=0,vjust=1,size=3.1,family='Times')+
    scale_fill_gradientn(colors=c('#f7fbff','#deebf7','#9ecae1','#fdae6b','#de2d26'),limits=c(0,1),name='Relative richness',
    breaks=c(0,.5,1),labels=c('0','0.5','1.0'))+
    coord_sf(crs=crs,datum=NA,expand=FALSE)+labs(title=selected_genera[i])+theme_void(base_family='Times',
    base_size=9)+theme(plot.title=element_text(face='bold.italic',size=9.5,hjust=.5,margin=margin(b=1)),legend.position='bottom',
    plot.margin=margin(1,2,1,2),plot.background=element_rect(fill='white',color=NA))
})
p11<-wrap_plots(ps,ncol=3,guides='collect') & theme(legend.position='bottom',legend.key.width=grid::unit(1.2,
  'cm'))
save_all(p11,'FigureS11',12.2,9)

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
    preview_image <- png::readPNG("FigureS11.png", native = TRUE)
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
