# Figure1：全球记录分布，沿用正文图样式。
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
suppressPackageStartupMessages({library(readr)
    library(dplyr)
    library(sf)
    library(ggplot2)
    library(patchwork)
    library(grid)})
world<-support$world
cols<-c(Hymenoptera='#B2D7B0',Diptera='#F7B7D2',Coleoptera='#8CC7DD',Lepidoptera='#F6D993')
legend_grobs<-list(textGrob('Pollinator order',x=.02,y=.5,just='left',gp=gpar(fontsize=10)))
for(i in seq_along(cols)) {
  taxon<-names(cols)[i]
  x<-.25+(i-1)*.19
  im<-support$icons[[paste0(tolower(taxon),'_pictogram.png')]]
  legend_grobs<-c(legend_grobs,list(circleGrob(x=x,y=.5,r=unit(.115,'in'),gp=gpar(fill=cols[[i]],col='#777777')),
    rasterGrob(im,x=x,y=.5,width=unit(.18,'in'),height=unit(.18,'in')),textGrob(taxon,x=x+.025,y=.5,just='left',
    gp=gpar(fontsize=9))))
}
leg<-do.call(grobTree,legend_grobs)
counts<-list()
bubbles<-function(raw,out) {
  d<-raw%>%mutate(latitude=as.numeric(latitude),longitude=as.numeric(longitude))%>%filter(is.finite(latitude),
    is.finite(longitude),between(latitude,-90,90),between(longitude,-180,180))
  counts[[out]]<<-data.frame(figure=out,total_records=nrow(raw),coordinate_records=nrow(d),sources_with_coordinates=n_distinct(d$reference_no),
    unmapped_records=nrow(raw)-nrow(d))
  d<-d%>%mutate(lat=floor(latitude/2)*2+1,lon=floor(longitude/2)*2+1)%>%count(pollinator_order,lat,lon,
    name='n')
  stopifnot(sum(d$n)==counts[[out]]$coordinate_records)
  breaks<-c(1,5,10,20)
  breaks<-breaks[breaks<=max(d$n)]
  p<-lapply(list(c('Hymenoptera','Diptera'),c('Coleoptera','Lepidoptera')),function(g)ggplot()+geom_sf(data=world,
    fill='#f5f5f5',color='grey75',linewidth=.13)+geom_point(data=filter(d,pollinator_order%in%g),aes(lon,
    lat,size=n,fill=pollinator_order),shape=21,alpha=.85,show.legend=TRUE)+scale_fill_manual(values=cols,
    limits=names(cols),drop=FALSE,name='Pollinator order',guide='none')+scale_size_area(max_size=12,limits=c(0,
    max(d$n)),breaks=breaks,name=if(out=='FigureS3')'Number of effect sizes' else 'Number of species observations')+coord_sf(ylim=c(-60,
    85))+theme_void()+theme(legend.position='bottom'))
  q<-((p[[1]]/p[[2]])+plot_layout(guides='collect')+plot_annotation(tag_levels='a')) & theme(legend.position='bottom',
    legend.box='vertical')
  q<-(wrap_elements(full=patchworkGrob(q))/wrap_elements(full=leg))+plot_layout(heights=c(10,.5))
  ggsave(paste0(out,'.png'),q,width=11,height=10.5,dpi=300,bg='white')
}
bubbles(results$evidence,"Figure1")

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
    preview_image <- png::readPNG("Figure1.png", native = TRUE)
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
