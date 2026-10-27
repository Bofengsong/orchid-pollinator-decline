# Figure S1：PRISMA-EcoEvo 筛选与纳入流程。
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
# 字号按投稿页面设计；检索路径横向排列，阶段从上到下。
# 历史计数使用现有可追溯记录；未知项留空，最终数量直接读取 Final_data.csv。
dat <- read.csv(file.path(data_dir, "Final_data.csv"), check.names = FALSE, fileEncoding = 'UTF-8-BOM')
meta <- dat[dat$meta_included == 1, ]
n_evidence <- nrow(dat)
n_evidence_sources <- length(unique(dat$reference_no))
n_effects <- nrow(meta)
n_meta_sources <- length(unique(meta$reference_no))
# 筛选计数由作者确认；最终纳入数量从主数据表核对。
stopifnot(n_evidence == 595, n_evidence_sources == 170,
          n_effects == 251, n_meta_sources == 71,
          9786 - 1065 == 8721, 8721 - 8523 == 198,
          198 - 17 == 181, 181 - 45 == 136,
          944 - 273 == 671, 671 - 538 == 133,
          133 - 4 == 129, 129 - 95 == 34,
          136 + 34 == 170, 506 + 89 == 595, 170 - 71 == 99,
          40 + 3 + 2 == 45, 7 + 92 == 99)
if (!exists('language')) language <- 'en'
zh <- identical(language, 'zh')
L <- function(en, cn) if (zh) cn else en
library(grid)
page_w <- 7.2
page_h <- 9.0
ink <- '#162C3A'
line_col <- '#526574'
stage_fill <- '#DCE9F1'
head_fill <- '#ECD28E'
font <- if (zh) 'ReviewChinese' else 'Helvetica'
g_x <- function(x) unit(x / page_w, 'npc')
g_y <- function(y) unit(1 - y / page_h, 'npc')
label <- function(t, x, y, size = 9.4, bold = FALSE, rot = 0) {
  grid.text(t, x = g_x(x), y = g_y(y), rot = rot,
    gp = gpar(fontfamily = font, fontsize = size,
      fontface = if (bold) 'bold' else 'plain', col = ink, lineheight = 1.2))
}
box <- function(x, y, w, h, lines, sizes = 9.4, bold = FALSE,
                fill = 'white', border = line_col, lwd = .85) {
  grid.rect(x = g_x(x+w/2), y = g_y(y+h/2),
    width = unit(w/page_w, 'npc'), height = unit(h/page_h, 'npc'),
    gp = gpar(fill = fill, col = border, lwd = lwd))
  sizes <- rep(sizes, length.out = length(lines))
  bold <- rep(bold, length.out = length(lines))
  step <- if (length(lines) == 1) 0 else min(.19, (h-.15)/(length(lines)-1))
  yy <- y+h/2+(seq_along(lines)-(length(lines)+1)/2)*step
  for (j in seq_along(lines)) {
    tg <- textGrob(lines[j], gp = gpar(fontfamily = font, fontsize = sizes[j],
      fontface = if (bold[j]) 'bold' else 'plain'))
    actual <- convertWidth(grobWidth(tg), 'inches', valueOnly = TRUE)
    if (actual > w-.10) stop(paste('Text exceeds box:', lines[j], actual, w))
    label(lines[j], x+w/2, yy[j], sizes[j], bold[j])
  }
}
phase <- function(y,h,en,cn,fill=stage_fill) {
  grid.rect(x=g_x(.255), y=g_y(y+h/2), width=unit(.31/page_w,'npc'),
    height=unit(h/page_h,'npc'), gp=gpar(fill=fill,col=line_col,lwd=.8))
  label(L(en,paste(strsplit(cn,'')[[1]],collapse='\n')), .255,y+h/2,9.8,TRUE,
    rot=if (zh) 0 else 90)
}
link <- function(xs,ys,end=TRUE) {
  grid.lines(x=g_x(xs),y=g_y(ys),
    arrow=if(end) arrow(angle=25,length=unit(.05,'inches'),type='closed') else NULL,
    gp=gpar(col=line_col,fill=line_col,lwd=.95,lineend='butt',linejoin='mitre'))
}
draw_flow <- function() {
  grid.newpage()
  box(.56,.12,3.14,.46,L('Primary database search','原始数据库检索'),11,TRUE,head_fill)
  box(3.93,.12,3.14,.46,
    L(c('Supplementary neutral-term search','and citation tracing'),c('补充中性词检索','与引文追溯')),10.3,TRUE,head_fill)
  phase(.78,1.26,'Identification','识别')
  phase(2.20,1.46,'Screening','筛选')
  phase(3.92,1.46,'Eligibility','资格评估')
  phase(5.62,3.24,'Included','纳入','#C6E0DA')

  # 左右两路采用相同网格；箭头只经过留白。
  for (x in c(.56,3.93)) {
    mid <- x+.95
    link(c(x+1.90,x+2.06),c(1.41,1.41))
    link(c(mid,mid),c(2.04,2.20))
    link(c(x+1.90,x+2.06),c(2.48,2.48))
    link(c(mid,mid),c(2.76,3.08))
    link(c(x+1.90,x+2.06),c(3.37,3.37))
    link(c(mid,mid),c(3.66,4.13))
    link(c(x+1.90,x+2.06),c(4.45,4.45))
    link(c(mid,mid),c(4.77,5.62))
  }
  link(c(2.13,2.13,5.50,5.50),c(6.40,6.57,6.57,6.40),FALSE)
  link(c(3.815,3.815),c(6.57,6.78))
  link(c(2.13,2.13),c(7.43,8.30))
  link(c(2.13,4.00),c(7.88,7.88))

  box(.56,.78,1.90,1.26,
    L(c('Web of Science','9,786 records','Hymenoptera 3,371','Diptera 1,862',
        'Coleoptera 2,181','Lepidoptera 2,372'),
      c('Web of Science','9,786 条检索记录','膜翅目 3,371','双翅目 1,862','鞘翅目 2,181','鳞翅目 2,372')),
    c(9.6,10.5,9.3,9.3,9.3,9.3),c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
  box(2.62,1.08,1.08,.66,L(c('Duplicates','removed','n = 1,065'),c('去除重复','记录','n = 1,065')),
    c(9.3,9.3,10),c(FALSE,FALSE,TRUE))
  box(3.93,.78,1.90,1.26,
    L(c('Supplementary search','and citation tracing','944 records in total','Search: 13 Sep 2026','Cut-off: 2 Oct 2025'),
      c('补充检索与引文追溯','合计 944 条检索记录','检索：2026-09-13','发表截止：2025-10-02')),
    if(zh)c(9.6,10.5,8.5,8.5) else c(9.6,9.6,10.5,8.5,8.5),
    if(zh)c(TRUE,TRUE,FALSE,FALSE) else c(TRUE,TRUE,TRUE,FALSE,FALSE))
  box(5.99,1.08,1.08,.66,L(c('Duplicates','removed','n = 273'),c('去除重复','记录','n = 273')),
    c(9.3,9.3,10),c(FALSE,FALSE,TRUE))
  for(x in c(.56,3.93)) {
    primary <- x == .56
    box(x,2.20,1.90,.56,L(c('Records screened',if(primary)'n = 8,721' else 'n = 671'),
      c('题名与摘要筛选',if(primary)'n = 8,721' else 'n = 671')),c(9.6,10.5),c(FALSE,TRUE))
    box(x+2.06,2.20,1.08,.56,L(c('Title / abstract','exclusions',if(primary)'n = 8,523' else 'n = 538'),
      c('题名／摘要','排除记录',if(primary)'n = 8,523' else 'n = 538')),c(8.6,9.0,10),c(FALSE,FALSE,TRUE))
    box(x,3.08,1.90,.58,L(c('Reports sought for retrieval',if(primary)'n = 198' else 'n = 133'),c('拟获取全文',if(primary)'n = 198' else 'n = 133')),
      c(9,10.5),c(FALSE,TRUE))
    box(x+2.06,3.08,1.08,.58,L(c('Not retrieved',if(primary)'n = 17' else 'n = 4'),c('未取得全文',if(primary)'n = 17' else 'n = 4')),
      c(9.2,10),c(FALSE,TRUE))
    box(x,4.13,1.90,.64,L(c('Full texts assessed',if(primary)'n = 181' else 'n = 129'),
      c('全文资格评估',if(primary)'n = 181' else 'n = 129')),c(9.6,10.5),c(FALSE,TRUE))
  }
  box(2.62,3.92,1.08,1.46,
    L(c('Reports excluded','n = 45','Orchid link: 40','Response: 3','Mapping: 2'),
      c('全文排除','n = 45','兰花关联：40','响应指标：3','记录对应：2')),
    c(8.5,10.2,8.8,8.8,8.8),c(FALSE,TRUE,FALSE,FALSE,FALSE))
  box(5.99,3.92,1.08,1.46,L(c('Reports excluded','n = 95'),
    c('全文排除','n = 95')),c(8.5,10.2),c(FALSE,TRUE))
  box(.56,5.62,3.14,.78,L(c('Included from primary search','136 publications','506 observations'),
    c('原始检索纳入','136 篇文献','506 条观察记录')),c(10.3,11.2,10.2),c(TRUE,TRUE,FALSE))
  box(3.93,5.62,3.14,.78,L(c('Included from supplementary sources','34 publications; 89 observations','30 direct + 4 citation-traced'),
    c('补充来源纳入','34 篇文献；89 条观察记录','直接获取 30 篇＋引文追溯 4 篇')),c(9.9,10.6,9.6),c(TRUE,TRUE,FALSE))
  box(.56,6.78,6.51,.65,L(c('Included in the systematic evidence set',
    paste(n_evidence_sources,'publications;',n_evidence,'observations')),
    c('系统证据集纳入',paste(n_evidence_sources,'篇文献；',n_evidence,'条观察记录'))),
    c(11.0,12.0),TRUE,'white','#315D58',1.2)
  box(4.00,7.51,3.07,.74,L(c('No effects in the primary meta-analysis','99 publications',
    'Seasonal / cycle: 7; other sources: 92'),c('未贡献主 Meta 效应量','99 篇文献','季节／周期比较：7；其他来源：92')),
    c(9.3,10.7,9.1),c(FALSE,TRUE,FALSE))
  box(.56,8.30,6.51,.56,L(c('Included in the primary meta-analysis',
    paste(n_meta_sources,'publications;',n_effects,'effect sizes')),
    c('主 Meta 分析纳入',paste(n_meta_sources,'篇文献；',n_effects,'个效应量'))),
    c(11.0,12.0),TRUE,'white','#315D58',1.2)
}

# 正式图只保留一种版式。PNG用于Word，PDF用于矢量排版。
png('FigureS1.png', width = page_w, height = page_h, units = 'in', res = 600)
draw_flow()
dev.off()
grDevices::pdf('FigureS1.pdf', width = page_w, height = page_h, family = 'Helvetica', useDingbats = FALSE)
draw_flow()
dev.off()

# 保存后回到可见的图形设备，在 Plots / Quartz 中预览。
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
    preview_image <- png::readPNG("FigureS1.png", native = TRUE)
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
