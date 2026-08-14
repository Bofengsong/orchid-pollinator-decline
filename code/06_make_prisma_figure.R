script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
audit_dir <- file.path(package_root, "data/audit")
out_dir <- file.path(package_root, "outputs/prisma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
audit <- read.csv(file.path(audit_dir, "FigureS1_PRISMA_counts_audit_2026-06-08.csv"), check.names = FALSE)
get_n <- function(item) audit$n[match(item, audit$`PRISMA item`)]

counts <- list(
  hym = get_n("Records identified: Hymenoptera"),
  dip = get_n("Records identified: Diptera"),
  col = get_n("Records identified: Coleoptera"),
  lep = get_n("Records identified: Lepidoptera"),
  total = get_n("Records identified: total"),
  dup = get_n("Duplicate records removed"),
  screened = get_n("Records screened"),
  title_abs = get_n("Records excluded at title/abstract screening"),
  sought = get_n("Reports sought for retrieval"),
  not_retrieved = get_n("Reports not retrieved"),
  fulltext = get_n("Full-text reports assessed for eligibility"),
  full_excl = get_n("Full-text reports excluded after orchid-pollinator eligibility check"),
  evidence_np = get_n("Studies included in global evidence base"),
  evidence_no = get_n("Species observations in global evidence base"),
  meta_excl = get_n("Reports excluded from quantitative meta-analysis"),
  meta_np = get_n("Studies included in quantitative meta-analysis"),
  meta_k = get_n("Independent effect sizes in quantitative meta-analysis")
)

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

png(file.path(out_dir, "FigureS1_PRISMA_EcoEvo_flow_updated_2026-06-08.png"),
    width = 4200, height = 2800, res = 300, bg = "white")

library(grid)
grid.newpage()
pushViewport(viewport(xscale = c(0, 1), yscale = c(0, 1)))

box_gp <- gpar(fill = "white", col = "black", lwd = 1.4)
head_gp <- gpar(fill = "#FFC000", col = "#B38A00", lwd = 1.4)
side_gp <- gpar(fill = "#9DC3E6", col = "black", lwd = 1.2)
txt_gp <- gpar(fontsize = 10.5, col = "black", fontfamily = "Helvetica")
small_gp <- gpar(fontsize = 9.5, col = "black", fontfamily = "Helvetica")
head_txt_gp <- gpar(fontsize = 11, fontface = "bold", fontfamily = "Helvetica")
side_txt_gp <- gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica")

draw_box <- function(x, y, w, h, label, gp = box_gp, text_gp = txt_gp) {
  grid.roundrect(x, y, w, h, r = unit(0, "npc"), gp = gp)
  grid.text(label, x, y, gp = text_gp, just = "center")
}

arrow_down <- function(x, y1, y2) {
  grid.lines(c(x, x), c(y1, y2), arrow = arrow(type = "closed", length = unit(0.10, "inches")), gp = gpar(lwd = 1.1))
}

arrow_right <- function(x1, x2, y) {
  grid.lines(c(x1, x2), c(y, y), arrow = arrow(type = "closed", length = unit(0.10, "inches")), gp = gpar(lwd = 1.1))
}

grid.text("PRISMA-EcoEvo screening flow diagram for the orchid pollinator evidence synthesis",
          x = 0.04, y = 0.955, just = "left", gp = gpar(fontsize = 13, fontface = "bold"))

draw_box(0.37, 0.885, 0.50, 0.043, "Identification of studies for the systematic evidence base", head_gp, head_txt_gp)
draw_box(0.77, 0.885, 0.27, 0.043, "Eligibility assessment for quantitative meta-analysis", head_gp, head_txt_gp)

grid.roundrect(0.065, 0.74, 0.035, 0.19, r = unit(0.015, "npc"), gp = side_gp)
grid.text("Identification", 0.065, 0.74, rot = 90, gp = side_txt_gp)
grid.roundrect(0.065, 0.45, 0.035, 0.35, r = unit(0.015, "npc"), gp = side_gp)
grid.text("Screening", 0.065, 0.45, rot = 90, gp = side_txt_gp)
grid.roundrect(0.065, 0.17, 0.035, 0.13, r = unit(0.015, "npc"), gp = side_gp)
grid.text("Included", 0.065, 0.17, rot = 90, gp = side_txt_gp)

draw_box(0.23, 0.73, 0.23, 0.15,
         paste0("Web of Science:\n",
                "- Hymenoptera (n = ", fmt(counts$hym), ")\n",
                "- Diptera (n = ", fmt(counts$dip), ")\n",
                "- Coleoptera (n = ", fmt(counts$col), ")\n",
                "- Lepidoptera (n = ", fmt(counts$lep), ")\n",
                "Total (n = ", fmt(counts$total), ")"),
         text_gp = small_gp)
draw_box(0.47, 0.73, 0.22, 0.12,
         paste0("Records removed before screening:\n",
                "- Duplicate records removed\n",
                "(n = ", fmt(counts$dup), ")"),
         text_gp = small_gp)
arrow_right(0.345, 0.36, 0.73)

draw_box(0.23, 0.57, 0.23, 0.075,
         paste0("Records screened\n(n = ", fmt(counts$screened), ")"))
draw_box(0.47, 0.57, 0.22, 0.095,
         paste0("Records excluded:\n- Title and abstract screening\n(n = ", fmt(counts$title_abs), ")"),
         text_gp = small_gp)
arrow_down(0.23, 0.655, 0.61)
arrow_right(0.345, 0.36, 0.57)

draw_box(0.23, 0.44, 0.23, 0.075,
         paste0("Reports sought for retrieval\n(n = ", fmt(counts$sought), ")"))
draw_box(0.47, 0.44, 0.22, 0.075,
         paste0("Reports not retrieved\n(n = ", fmt(counts$not_retrieved), ")"))
arrow_down(0.23, 0.532, 0.48)
arrow_right(0.345, 0.36, 0.44)

draw_box(0.23, 0.305, 0.23, 0.085,
         paste0("Full-text reports assessed\nfor eligibility\n(n = ", fmt(counts$fulltext), ")"),
         text_gp = small_gp)
draw_box(0.47, 0.305, 0.22, 0.11,
         paste0("Reports excluded:\n",
                "- Not retained after\n",
                "  orchid-pollinator eligibility check\n",
                "(n = ", fmt(counts$full_excl), ")"),
         text_gp = small_gp)
arrow_down(0.23, 0.402, 0.352)
arrow_right(0.345, 0.36, 0.305)

draw_box(0.23, 0.15, 0.23, 0.105,
         paste0("Studies included in the global\n",
                "evidence base\n",
                "(np = ", fmt(counts$evidence_np), "; no = ", fmt(counts$evidence_no), ")"),
         text_gp = small_gp)
arrow_down(0.23, 0.26, 0.205)

draw_box(0.76, 0.305, 0.24, 0.105,
         paste0("Reports excluded:\n",
                "- Missing decline value\n",
                "  and/or sample size for lnRR\n",
                "(n = ", fmt(counts$meta_excl), ")"),
         text_gp = small_gp)
draw_box(0.76, 0.15, 0.24, 0.105,
         paste0("Studies included in quantitative\n",
                "meta-analysis\n",
                "(np = ", fmt(counts$meta_np), "; k = ", fmt(counts$meta_k), ")"),
         text_gp = small_gp)
arrow_right(0.58, 0.64, 0.305)
arrow_down(0.76, 0.252, 0.205)

grid.text("Note: np = number of publications; no = species observations; k = independent lnRR effect sizes.",
          x = 0.04, y = 0.045, just = "left", gp = gpar(fontsize = 9))

dev.off()
