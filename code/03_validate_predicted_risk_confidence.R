suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(terra)
})

# Package-safe risk-layer step.
# The original full raw-raster workflow requires large external files that are
# intentionally not copied into clean_package. This script validates the current
# generated 1-degree risk/confidence rasters and copies them into results.

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

input_dir <- file.path(package_root, "data/spatial_inputs/risk_layers")
outdir <- file.path(package_root, "outputs/spatial_risk")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "Predicted_orchid_pollinator_risk_1deg.tif",
  "Confidence_high_risk_1deg.tif",
  "Prediction_uncertainty_sd_1deg.tif",
  "PredictedRisk_1deg_grid_values.csv",
  "RiskBootstrap_weight_summary.csv",
  "RiskBootstrap_weights_all_iterations.csv",
  "RiskBootstrap_effect_size_exposure_table.csv"
)

missing <- required[!file.exists(file.path(input_dir, required))]
if (length(missing) > 0) {
  stop(
    "Missing packaged risk-layer input(s):\n",
    paste(file.path(input_dir, missing), collapse = "\n"),
    "\n\nThe full raw-raster regeneration requires the external sources listed in ",
    "data/raw_source_manifest/RAW_SOURCES_NOT_COPIED.md."
  )
}

file.copy(file.path(input_dir, required), file.path(outdir, required), overwrite = TRUE)

risk <- rast(file.path(outdir, "Predicted_orchid_pollinator_risk_1deg.tif"))
conf <- rast(file.path(outdir, "Confidence_high_risk_1deg.tif"))
uncert <- rast(file.path(outdir, "Prediction_uncertainty_sd_1deg.tif"))

summary_tbl <- tibble(
  layer = c("predicted_risk", "confidence_high_risk", "prediction_uncertainty_sd"),
  file = c(
    "Predicted_orchid_pollinator_risk_1deg.tif",
    "Confidence_high_risk_1deg.tif",
    "Prediction_uncertainty_sd_1deg.tif"
  ),
  ncell = c(ncell(risk), ncell(conf), ncell(uncert)),
  n_non_na = c(global(!is.na(risk), "sum", na.rm = TRUE)[1, 1],
               global(!is.na(conf), "sum", na.rm = TRUE)[1, 1],
               global(!is.na(uncert), "sum", na.rm = TRUE)[1, 1]),
  min = c(global(risk, "min", na.rm = TRUE)[1, 1],
          global(conf, "min", na.rm = TRUE)[1, 1],
          global(uncert, "min", na.rm = TRUE)[1, 1]),
  max = c(global(risk, "max", na.rm = TRUE)[1, 1],
          global(conf, "max", na.rm = TRUE)[1, 1],
          global(uncert, "max", na.rm = TRUE)[1, 1])
)

write_csv(summary_tbl, file.path(outdir, "packaged_risk_layer_validation_summary.csv"))

cat("Validated packaged risk/confidence layers in:\n", outdir, "\n", sep = "")
