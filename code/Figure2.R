# Figure 2：森林图，k为效应量数，n为文献数；地区面板不显示Africa。
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
if (!file.exists("Results.rds")) stop("请先运行 Analysis.R。")
results <- readRDS("Results.rds")
if (!identical(unname(tools::md5sum(file.path(data_dir, "Final_data.csv"))), results$input_md5)) {
  stop("Final_data.csv 已变化，请先重新运行 Analysis.R，再画图。")
}
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(metafor)
    library(scales)
    library(tibble)
    library(grid)
})

to_pct <- function(x) 100 * (exp(x) - 1)

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

# Broad country-level biogeographic grouping used for the regional moderator.
# Records from Costa Rica and Panama in Central America, the Dominican Republic
# in the Caribbean, and all South American records are combined operationally
# as Central-South America. Mexico is retained in North America because a
# whole-country assignment cannot resolve its Nearctic-Neotropical transition;
# no Mexican record entered the meta-analysis, so this choice does not affect
# the regional estimates reported here.
assign_biogeographic_region <- function(country, continent) {
  country <- str_squish(as.character(country))
  out <- as.character(continent)
  out[country %in% c("Costa Rica", "Dominican Republic", "Panama",
    "Brazil", "Chile", "Ecuador")] <- "Central-South America"
  out[country %in% c("Canada", "Canade", "USA", "USA, Canada",
    "Canada/USA", "Mexico")] <- "North America"
  out
}

clean_order <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  case_when(
    x == "coleoptera" ~ "Coleoptera",
    x == "diptera" ~ "Diptera",
    x == "hymenoptera" ~ "Hymenoptera",
    x == "lepidoptera" ~ "Lepidoptera",
    TRUE ~ str_to_title(x)
  )
}

fit_level_mv <- function(df_level,
  rand_list = list(~1 | reference_no, ~1 | country, ~1 | pollinator_family)) {
  if (nrow(df_level) < 2) return(NULL)
  df_level <- df_level %>%
  mutate(yi = as.numeric(yi), vi = as.numeric(vi)) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)
  if (nrow(df_level) < 2) return(NULL)

  fit <- try(rma.mv(yi = yi, V = vi, random = rand_list, data = df_level, method = "REML"), silent = TRUE)
  if (inherits(fit, "try-error") || any(!is.finite(fit$b))) {
    fit <- try(rma.mv(yi, vi, random = rand_list, data = df_level, method = "REML", control=list(optimizer="optim",
      optmethod="BFGS")), silent = TRUE)
    if (inherits(fit, "try-error") || any(!is.finite(fit$b))) stop("Subgroup failed without simplifying random effects")
  }

  beta <- as.numeric(fit$b)
  se <- as.numeric(fit$se)
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  ci_lb <- as.numeric(fit$ci.lb)
  ci_ub <- as.numeric(fit$ci.ub)

  pi_lb <- pi_ub <- NA_real_
  pi_obj <- try(if (inherits(fit, "rma.mv")) predict(fit, pi.type = "t") else predict(fit), silent = TRUE)
  if (!inherits(pi_obj, "try-error") && all(c("pi.lb", "pi.ub") %in% names(pi_obj))) {
    pi_lb <- as.numeric(pi_obj$pi.lb)
    pi_ub <- as.numeric(pi_obj$pi.ub)
  }
  if (is.na(pi_lb) || is.na(pi_ub)) {
    tau2_total <- if (inherits(fit, "rma.mv")) sum(as.numeric(fit$sigma2)) else as.numeric(fit$tau2)
    hw <- 1.96 * sqrt(se^2 + tau2_total)
    pi_lb <- beta - hw
    pi_ub <- beta + hw
  }

  if (nrow(df_level) <= 2) { pi_lb <- pi_ub <- NA_real_ }
  tibble(k = nrow(df_level), yi = beta, se = se, zval = zval, pval = pval,
    ci_lb = ci_lb, ci_ub = ci_ub, pi_lb = pi_lb, pi_ub = pi_ub)
}

group_density <- function(x) {
  if (length(unique(x)) < 2) return(rep(1, length(x)))
  d <- density(x)
  approx(d$x, d$y, xout = x, rule = 2)$y
}

sig_star <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

fmt_pct_text <- function(est, lb, ub, pi_lb, pi_ub, k) {
  fmt <- function(v) ifelse(abs(v)>9999,formatC(v,format="e",digits=1),sprintf("%.1f",v))
  paste0(fmt(abs(to_pct(est))), "% ",ifelse(est<0,"reduction","increase"),"\n",
    "CI: [",fmt(to_pct(lb)),"; ",fmt(to_pct(ub)),"]%\n",
    if(k<=2) "PI: not estimated (k <= 2)" else paste0("PI: [",fmt(to_pct(pi_lb)),"; ",fmt(to_pct(pi_ub)),
    "]%"))
}

meta <- results$meta %>%
mutate(
  reference_no = as.character(reference_no),
  country = ifelse(is.na(country) & reference_no == "D419", "Regional_D419_Southeast_Europe", as.character(country)),

  pollinator_family = as.character(pollinator_family),
  pollinator_order = clean_order(pollinator_order),
  continent = clean_continent(continent),
  biogeographic_region = assign_biogeographic_region(country, continent),
  yi = suppressWarnings(as.numeric(yi)),
  vi = suppressWarnings(as.numeric(vi))
) %>%
filter(is.finite(yi), is.finite(vi), vi > 0,
  !is.na(reference_no), !is.na(country), !is.na(pollinator_family))

long <- results$driver_data %>%
mutate(
  reference_no = as.character(reference_no),
  country = ifelse(is.na(country) & reference_no == "D419", "Regional_D419_Southeast_Europe", as.character(country)),

  pollinator_family = as.character(pollinator_family),
  pollinator_order = clean_order(pollinator_order),
  continent = clean_continent(continent),
  biogeographic_region = assign_biogeographic_region(country, continent),
  driver = str_to_title(str_squish(as.character(driver))),
  yi = suppressWarnings(as.numeric(yi)),
  vi = suppressWarnings(as.numeric(vi))
) %>%
filter(is.finite(yi), is.finite(vi), vi > 0,
  !is.na(reference_no), !is.na(country), !is.na(pollinator_family), !is.na(driver))

rand_list <- list(~1 | reference_no, ~1 | country, ~1 | pollinator_family)

driver_levels <- c("Multiple", "Climate", "Habitat", "Invasives", "Pesticides", "Urbanization")
order_levels <- c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera")
# Display selection requested by the author; underlying analysis data are unchanged.
region_levels <- c("Asia", "Europe", "North America", "Oceania", "Central-South America")

build_block <- function(raw_df, block_name, group_var, levels_keep, fill_col) {
  raw_df <- raw_df %>% filter(.data[[group_var]] %in% levels_keep)
  pts <- raw_df %>%
  mutate(level = .data[[group_var]]) %>%
  group_by(level) %>%
  mutate(
    dens = group_density(yi),
    alpha_pt = scales::rescale(dens, to = c(0.18, 0.72)),
    weight_raw = scales::rescale(1 / vi, to = c(2.0, 4.0))
  ) %>%
  ungroup() %>%
  mutate(block = block_name, fill_col = fill_col) %>%
  select(block, level, yi, vi, alpha_pt, weight_raw, fill_col, reference_no, country, pollinator_family)

  summ <- map_dfr(levels_keep, function(lv) {
      dsub <- raw_df %>% filter(.data[[group_var]] == lv)
      fit <- fit_level_mv(dsub, rand_list = rand_list)
      if (is.null(fit)) return(NULL)
      fit %>%
      mutate(block = block_name, level = lv, fill_col = fill_col,
      p_star = sig_star(pval), publications = n_distinct(dsub$reference_no),
      label_k = paste0("k = ", k, "\nn = ", publications))
  })

  unestimable <- map_dfr(levels_keep, function(lv) {
      dsub <- raw_df %>% filter(.data[[group_var]] == lv)
      if (nrow(dsub) >= 2) return(NULL)
      tibble(
      block = block_name,
      level = lv,
      k = nrow(dsub),
      fill_col = fill_col,
      label_left = if (nrow(dsub) == 1) {
        "Not estimated\n(single effect size)"
    } else {
        "Not estimated\n(no effect sizes)"
    },
      label_k = paste0("k = ", nrow(dsub))
    )
  })

  list(pts = pts, summ = summ, unestimable = unestimable)
}

block_driver <- build_block(long %>% filter(driver %in% driver_levels), "Drivers", "driver", driver_levels,
  "#fabdb7")
block_order <- build_block(meta %>% filter(pollinator_order %in% order_levels), "Pollinator Orders", "pollinator_order",
  order_levels, "#bcd689")
block_cont <- build_block(
  meta %>% filter(biogeographic_region %in% region_levels),
  "Biogeographic Regions", "biogeographic_region", region_levels, "#81dee0"
)

pts_all <- bind_rows(block_driver$pts, block_order$pts, block_cont$pts)
summ_all <- bind_rows(block_driver$summ, block_order$summ, block_cont$summ)
unestimable_all <- bind_rows(
  block_driver$unestimable,
  block_order$unestimable,
  block_cont$unestimable
)
if (ncol(unestimable_all) == 0) {
  unestimable_all <- tibble(
    block = character(), level = character(), k = integer(),
    fill_col = character(), label_left = character(), label_k = character()
  )
}

layout_tbl <- tibble(
  block = c(rep("Drivers", length(driver_levels)),
  rep("Pollinator Orders", length(order_levels)),
  rep("Biogeographic Regions", length(region_levels))),
  level = c(driver_levels, order_levels, region_levels)
) %>%
mutate(
  display_level = level,
  row_index = row_number(),
  group_gap = case_when(
  block == "Pollinator Orders" ~ 1.05,
  block == "Biogeographic Regions" ~ 2.00,
  TRUE ~ 0
),
  y = 24.8 - (row_index - 1) * (20.10 / (n() - 1)) - group_gap
)

pts_all <- pts_all %>% left_join(layout_tbl, by = c("block", "level"))
summ_all <- summ_all %>%
left_join(layout_tbl, by = c("block", "level")) %>%
mutate(
  label_left = pmap_chr(
  list(yi, ci_lb, ci_ub, pi_lb, pi_ub, k),
  ~ fmt_pct_text(..1, ..2, ..3, ..4, ..5, ..6)
)
)

unestimable_all <- unestimable_all %>%
left_join(layout_tbl, by = c("block", "level"))

block_label_tbl <- tibble(
  block = c("Drivers", "Pollinator Orders", "Biogeographic Regions"),
  y_mid = c(mean(layout_tbl$y[layout_tbl$block == "Drivers"]),
  mean(layout_tbl$y[layout_tbl$block == "Pollinator Orders"]),
  mean(layout_tbl$y[layout_tbl$block == "Biogeographic Regions"])),
  y_min = c(min(layout_tbl$y[layout_tbl$block == "Drivers"]),
  min(layout_tbl$y[layout_tbl$block == "Pollinator Orders"]),
  min(layout_tbl$y[layout_tbl$block == "Biogeographic Regions"])),
  y_max = c(max(layout_tbl$y[layout_tbl$block == "Drivers"]),
  max(layout_tbl$y[layout_tbl$block == "Pollinator Orders"]),
  max(layout_tbl$y[layout_tbl$block == "Biogeographic Regions"]))
)

# -------------------------------------------------------------------------
# Layout coordinate system. lnRR values are mapped into a dedicated plot
# column so uncertainty intervals cannot overlap the left text column.
# -------------------------------------------------------------------------
x_ln_min <- min(-6.2, floor(min(c(pts_all$yi,summ_all$pi_lb,summ_all$ci_lb),na.rm=TRUE)))
x_ln_max <- max(3.2, ceiling(max(c(pts_all$yi,summ_all$pi_ub,summ_all$ci_ub),na.rm=TRUE)))
x_plot_left <- 8.10
x_plot_right <- 13.65
x_plot_w <- x_plot_right - x_plot_left

x_map <- function(x) x_plot_left + (pmin(pmax(x, x_ln_min), x_ln_max) - x_ln_min) / (x_ln_max - x_ln_min) * x_plot_w
x_map_interval <- function(x) {
  raw_x <- x_plot_left + (x - x_ln_min) / (x_ln_max - x_ln_min) * x_plot_w
  pmin(pmax(raw_x, x_plot_left), 14.35)
}

pts_all <- pts_all %>%
mutate(x = x_map(yi))

summ_all <- summ_all %>%
mutate(
  x = x_map(yi),
  x_ci_lb = x_map_interval(ci_lb),
  x_ci_ub = x_map_interval(ci_ub),
  x_pi_lb = x_map_interval(pi_lb),
  x_pi_ub = x_map_interval(pi_ub),
  sig_y = y + 0.50
)

set.seed(42)

summ_all$stat_colour <- summ_all$k_colour <- summ_all$star_colour <- "black"
p <- ggplot() +
annotate("rect", xmin = 3.55, xmax = 16.70, ymin = 1.30, ymax = 25.75,
  fill = NA, color = "black", linewidth = 0.9) +
annotate("segment", x = x_map(-5), xend = x_map(-5), y = 1.30, yend = 25.75,
  color = "grey78", linewidth = 0.65) +
annotate("segment", x = x_map(0), xend = x_map(0), y = 1.30, yend = 25.75,
  color = "black", linewidth = 0.75, linetype = "dashed") +
geom_segment(
  data = summ_all,
  aes(x = x_pi_lb, xend = x_pi_ub, y = y, yend = y),
  linewidth = 0.46, color = "grey55", lineend = "round"
) +
geom_segment(
  data = summ_all,
  aes(x = x_ci_lb, xend = x_ci_ub, y = y, yend = y),
  linewidth = 1.30, color = "black", lineend = "round"
) +
geom_point(
  data = pts_all,
  aes(x = x, y = y, size = weight_raw, alpha = alpha_pt, fill = fill_col),
  shape = 21, stroke = 0.25, colour = "black",
  position = position_jitter(width = 0, height = 0.12)
) +
geom_point(
  data = summ_all,
  aes(x = x, y = y, size = k, fill = fill_col),
  shape = 21, stroke = 0.85, colour = "black"
) +
geom_text(
  data = summ_all,
  aes(x = x, y = sig_y, label = p_star),
  family = "sans", fontface = "bold", size = 4.05, color = summ_all$star_colour
) +
geom_text(
  data = summ_all,
  aes(x = 5.25, y = y + 0.06, label = label_left),
  family = "sans", size = 3.15, lineheight = 0.88, hjust = 0.5, vjust = 0.5, color = summ_all$stat_colour
) +
geom_text(
  data = summ_all,
  aes(x = 15.15, y = y, label = label_k),
  family = "sans", fontface = "italic", size = 4.35, hjust = 0, vjust = 0.5, color = summ_all$k_colour
) +
geom_text(
  data = unestimable_all,
  aes(x = 5.25, y = y + 0.06, label = label_left),
  family = "sans", size = 3.15, lineheight = 0.88, hjust = 0.5, vjust = 0.5,
  color = "grey35"
) +
geom_text(
  data = unestimable_all,
  aes(x = 15.15, y = y, label = label_k),
  family = "sans", fontface = "italic", size = 4.35, hjust = 0, vjust = 0.5,
  color = "black"
) +
geom_text(
  data = layout_tbl,
  aes(x = 0.56, y = y, label = display_level),
  family = "sans", size = 4.9, hjust = 0, vjust = 0.5
) +
geom_text(
  data = block_label_tbl,
  aes(x = 0.18, y = y_mid, label = block),
  family = "sans", fontface = "bold", size = 6.8, angle = 90, color = "black"
) +
geom_segment(
  data = block_label_tbl,
  aes(x = 0.48, xend = 0.48, y = y_min - 0.52, yend = y_max + 0.52),
  linewidth = 1.15, color = "black"
) +
annotate("segment", x = 3.55, xend = 16.70, y = 1.30, yend = 1.30,
  color = "black", linewidth = 0.9) +
annotate("text", x = x_map(-5), y = 0.82, label = "-5",
  family = "sans", size = 6.2, color = "black") +
annotate("text", x = x_map(0), y = 0.82, label = "0",
  family = "sans", size = 6.2, color = "black") +
annotate("text", x = x_map(5), y = 0.82, label = "5", family="sans", size=6.2) +
annotate("text", x = mean(c(x_plot_left, x_plot_right)), y = 0.28,
  label = "Effect Size (lnRR)", family = "sans", fontface = "bold",
  size = 6.5, color = "black") +
scale_size(range = c(2.4, 9.2), guide = "none") +
scale_alpha_identity() +
scale_fill_identity() +
coord_cartesian(xlim = c(0, 17.0), ylim = c(0.05, 26.05), clip = "off") +
theme_void() +
theme(
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  legend.position = "none",
  plot.margin = margin(10, 12, 18, 10)
)


ggsave("Figure2.png", p, width = 12.8, height = 11.7, dpi = 600, bg = "white")

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
    preview_image <- png::readPNG("Figure2.png", native = TRUE)
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
