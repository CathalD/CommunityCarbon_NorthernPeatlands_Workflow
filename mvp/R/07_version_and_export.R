# =============================================================================
# 07_version_and_export.R
#
# Snapshot this run into a numbered, dated version folder with a short
# metadata record, build an optional hexagon reporting layer, and export
# the core points as GIS-ready files.
#
# This is what turns "every field season improves the map" from a claim into
# something you can check: outputs/versions/carbon_map_v1/, v2/, ... each
# with a metadata.csv showing what changed.
#
# INPUT   everything in mvp/outputs/current/
# OUTPUT  mvp/outputs/versions/carbon_map_vN/   full snapshot + metadata.csv
#         mvp/outputs/current/cores.gpkg
#         mvp/outputs/current/hex_carbon_layer.gpkg
#
# REQUIRES: sf, terra
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  # Works three ways: Rscript (finds --file=), source() from ANY working
  # directory (finds the sourced file's own path), and run_all.R (sets
  # MVP_R_DIR). Only falls back to getwd() if all three fail.
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

library(sf)
library(terra)
library(dplyr)

msg("07  version, archive & export")

req <- c("carbon_posterior_mean.tif", "carbon_posterior_sd.tif",
        "carbon_difference_from_prior.tif", "cores_clean.geojson",
        "validation_metrics.csv", "prediction_residuals.csv")
missing <- req[!file.exists(file.path(CFG$dir_current, req))]
if (length(missing)) {
  stop("Missing from outputs/current/: ", paste(missing, collapse = ", "),
      " -- run steps 1-6 first.")
}

# ---- 1. optional hexagon reporting layer -----------------------------------

posterior_mean <- rast(file.path(CFG$dir_current, "carbon_posterior_mean.tif"))
posterior_sd   <- rast(file.path(CFG$dir_current, "carbon_posterior_sd.tif"))
cores <- st_read(file.path(CFG$dir_current, "cores_clean.geojson"), quiet = TRUE)

aoi_extent <- st_as_sfc(st_bbox(posterior_mean), crs = CFG$crs_geographic)
hex_cellsize_m <- CFG$gee$scale_m * 40   # a few dozen pixels per hexagon
hex_geom <- st_make_grid(st_transform(aoi_extent, CFG$crs_equal_area),
                         cellsize = hex_cellsize_m, square = FALSE)
hex <- st_transform(st_sf(hex_id = seq_along(hex_geom), geometry = hex_geom),
                    CFG$crs_geographic)

hex_mean <- terra::extract(posterior_mean, vect(hex), fun = mean, na.rm = TRUE)
hex_sd   <- terra::extract(posterior_sd, vect(hex), fun = mean, na.rm = TRUE)
hex$mean_carbon_kgm2 <- hex_mean[, 2]
hex$mean_uncertainty_kgm2 <- hex_sd[, 2]
hex$n_cores_within <- lengths(st_intersects(hex, cores))
hex <- hex[!is.na(hex$mean_carbon_kgm2), ]

st_write(hex, file.path(CFG$dir_current, "hex_carbon_layer.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
msg("wrote hex_carbon_layer.gpkg  (", nrow(hex), " hexagons)")

# ---- 2. export core points as GIS-ready GeoPackage -------------------------

st_write(cores, file.path(CFG$dir_current, "cores.gpkg"), delete_dsn = TRUE, quiet = TRUE)
msg("wrote cores.gpkg")

# ---- 3. metadata for this version -------------------------------------------

metrics <- read.csv(file.path(CFG$dir_current, "validation_metrics.csv"), stringsAsFactors = FALSE)
residuals <- read.csv(file.path(CFG$dir_current, "prediction_residuals.csv"), stringsAsFactors = FALSE)

prev_dir <- latest_version_dir(CFG$dir_versions)
prev_uncertainty <- if (!is.null(prev_dir) &&
                        file.exists(file.path(prev_dir, "carbon_posterior_sd.tif"))) {
  global(rast(file.path(prev_dir, "carbon_posterior_sd.tif")), "mean", na.rm = TRUE)[1, 1]
} else NA_real_
this_uncertainty <- global(posterior_sd, "mean", na.rm = TRUE)[1, 1]

area_updated_km2 <- {
  diff_r <- rast(file.path(CFG$dir_current, "carbon_difference_from_prior.tif"))
  cell_km2 <- prod(res(diff_r)) * 111 * 111  # rough deg->km conversion at this latitude
  meaningfully_changed <- abs(diff_r) > 0.1  # kg/m2, an arbitrary "visible change" threshold
  global(meaningfully_changed, "sum", na.rm = TRUE)[1, 1] * cell_km2
}

version_num <- next_version_number(CFG$dir_versions)
version_dir <- file.path(CFG$dir_versions, sprintf("carbon_map_v%d", version_num))
ensure_dir(version_dir)

metadata <- data.frame(
  field = c("version", "date", "prior_source", "n_cores_used",
           "loco_rmse_kgm2", "loco_r2", "beats_null",
           "mean_uncertainty_kgm2", "prev_mean_uncertainty_kgm2",
           "area_meaningfully_updated_km2"),
  value = c(version_num, as.character(Sys.Date()), CFG$bayes$prior_source,
           nrow(residuals),
           metrics$value[metrics$metric == "loco_rmse_kgm2"],
           metrics$value[metrics$metric == "loco_r2"],
           metrics$value[metrics$metric == "beats_null"],
           round(this_uncertainty, 3), round(prev_uncertainty, 3),
           round(area_updated_km2, 2))
)
write.csv(metadata, file.path(version_dir, "metadata.csv"), row.names = FALSE)

# ---- 4. archive everything in outputs/current/ into this version folder ----

current_files <- list.files(CFG$dir_current, full.names = TRUE)
file.copy(current_files, version_dir, overwrite = TRUE, recursive = TRUE)

msg("archived this run to ", version_dir)
print(metadata)

msg("07 complete -- carbon_map_v", version_num, " is the new current best estimate")
