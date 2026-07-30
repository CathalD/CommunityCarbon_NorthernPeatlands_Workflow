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

# Grid over the AOI POLYGON, not the raster's bounding box. The polygon is
# ~5,400 km2; its bounding box is nearly four times that, and every extra
# hexagon costs a zonal statistic. Fall back to the bbox only if step 02's AOI
# file is missing.
aoi_path <- file.path(CFG$dir_current, "aoi.geojson")
aoi_poly <- if (file.exists(aoi_path)) {
  st_union(st_geometry(st_read(aoi_path, quiet = TRUE)))
} else {
  msg("no aoi.geojson; falling back to the raster bounding box (more hexagons)")
  st_as_sfc(st_bbox(posterior_mean), crs = CFG$crs_geographic)
}

#' Build one hexagon reporting layer at a given cell size, with the posterior
#' mean and uncertainty averaged inside each cell.
#'
#' RASTERIZE THEN zonal(), NOT extract() PER POLYGON. terra::extract() with a
#' summary function walks the raster once per polygon, so at 500 m over this AOI
#' it is tens of thousands of passes and effectively hangs. rasterize() writes
#' the hexagon id into the raster grid once, and zonal() then computes every
#' hexagon's mean in a single pass.
build_hex <- function(cellsize_m) {
  aoi_ea <- st_transform(aoi_poly, CFG$crs_equal_area)
  g <- st_make_grid(aoi_ea, cellsize = cellsize_m, square = FALSE)
  h <- st_sf(hex_id = seq_along(g), geometry = g)
  # Drop hexagons that never touch the study area before doing any raster work.
  h <- h[lengths(st_intersects(h, aoi_ea)) > 0, ]
  h <- st_transform(h, CFG$crs_geographic)
  msg("    ", cellsize_m, " m: ", nrow(h), " candidate hexagons over the AOI")

  hz <- terra::rasterize(vect(h), posterior_mean, field = "hex_id")
  zm <- terra::zonal(posterior_mean, hz, fun = "mean", na.rm = TRUE)
  zs <- terra::zonal(posterior_sd,   hz, fun = "mean", na.rm = TRUE)
  names(zm) <- c("hex_id", "mean_carbon_kgm2")
  names(zs) <- c("hex_id", "mean_uncertainty_kgm2")

  h <- merge(h, zm, by = "hex_id", all.x = TRUE)
  h <- merge(h, zs, by = "hex_id", all.x = TRUE)
  h$n_cores_within <- lengths(st_intersects(h, cores))
  h$hex_size_m     <- cellsize_m
  h[!is.na(h$mean_carbon_kgm2), ]
}

hex_path <- file.path(CFG$dir_current, "hex_carbon_layer.gpkg")
if (file.exists(hex_path)) unlink(hex_path)

hex_counts <- integer(0)
for (sz in CFG$hex_sizes_m) {
  h <- build_hex(sz)
  # An empty layer means the posterior raster carries no data anywhere -- don't
  # write a blank file and call it a deliverable.
  if (!nrow(h)) {
    stop("Hexagon layer at ", sz, " m is empty: carbon_posterior_mean.tif has ",
        "no valid cells. Re-check step 6's coverage line -- if it reported NaN ",
        "or 0%, the fusion produced nothing and this archive would be blank.")
  }
  lyr <- sprintf("hex_%dm", sz)
  st_write(h, hex_path, layer = lyr, append = file.exists(hex_path), quiet = TRUE)
  hex_counts[lyr] <- nrow(h)
  msg("  ", lyr, ": ", nrow(h), " hexagons")
}
msg("wrote hex_carbon_layer.gpkg  (", length(CFG$hex_sizes_m), " scales: ",
   paste0(CFG$hex_sizes_m, "m", collapse = ", "), ")")

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
  # cellSize() gives true per-cell ground area. A flat degrees->km factor is
  # wrong here: at 56N one degree of longitude spans ~62 km, not 111, so a
  # squared-111 conversion overstates the area by roughly 1.8x.
  changed <- abs(diff_r) > 0.1   # kg/m2 -- the "visibly different" threshold
  global(changed * cellSize(diff_r, unit = "km"), "sum", na.rm = TRUE)[1, 1]
}

version_num <- next_version_number(CFG$dir_versions)
version_dir <- file.path(CFG$dir_versions, sprintf("carbon_map_v%d", version_num))
ensure_dir(version_dir)

metadata <- data.frame(
  field = c("version", "date", "prior_source", "n_cores_used",
           "loco_rmse_kgm2", "loco_r2_variance_explained", "beats_null",
           "mean_uncertainty_kgm2", "prev_mean_uncertainty_kgm2",
           "area_meaningfully_updated_km2"),
  value = c(version_num, as.character(Sys.Date()), CFG$bayes$prior_source,
           nrow(residuals),
           metrics$value[metrics$metric == "loco_rmse_kgm2"],
           metrics$value[metrics$metric == "loco_r2_variance_explained"],
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
