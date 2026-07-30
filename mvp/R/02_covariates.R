# =============================================================================
# 02_covariates.R
#
# Build a small AOI around the cores, pull a covariate stack + the published
# prior maps from Earth Engine, and download everything as local GeoTIFFs.
#
# REQUIRES an authenticated Earth Engine session:
#   install.packages("rgee"); rgee::ee_install(); rgee::ee_Initialize(project = "your-project")
# Run once interactively before this script if you haven't already.
#
# INPUT   mvp/outputs/current/cores_clean.geojson   (from 01)
# OUTPUT  mvp/outputs/current/predictors.tif          multi-band covariate stack
#         mvp/outputs/current/prior_li2025.tif         full-column prior, kg C/m2
#         mvp/outputs/current/prior_sothe_0_30.tif     0-30cm prior, kg C/m2
#         mvp/outputs/current/aoi.geojson
#
# REQUIRES: sf, rgee, terra
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
library(rgee)
library(dplyr)   # for the %>% pipe below

msg("02  covariates + priors (Earth Engine)")

cores_path <- file.path(CFG$dir_current, "cores_clean.geojson")
if (!file.exists(cores_path)) {
  stop("Missing ", cores_path, " -- run 01_clean_and_stocks.R first.")
}
cores <- st_read(cores_path, quiet = TRUE)

ee_Initialize(project = CFG$gee$project)
msg("Earth Engine initialised, project: ", CFG$gee$project)

# ---- 1. AOI: a simple buffered hull, not a coast-oriented rectangle -------

aoi_hull <- cores %>%
  st_transform(CFG$crs_equal_area) %>%
  st_union() %>%
  st_convex_hull() %>%
  st_buffer(CFG$gee$aoi_buffer_km * 1000) %>%
  st_transform(CFG$crs_geographic)

st_write(st_sf(geometry = aoi_hull), file.path(CFG$dir_current, "aoi.geojson"),
         delete_dsn = TRUE, quiet = TRUE)
msg("AOI: convex hull of cores + ", CFG$gee$aoi_buffer_km, " km buffer")

aoi_ee <- sf_as_ee(st_sf(geometry = aoi_hull))$geometry()

# ---- 2. covariate stack ----------------------------------------------------

season_filter <- function(coll) {
  coll$
    filter(ee$Filter$calendarRange(CFG$gee$season_start_doy,
                                   CFG$gee$season_end_doy, "day_of_year"))$
    filter(ee$Filter$calendarRange(min(CFG$gee$season_years),
                                   max(CFG$gee$season_years), "year"))
}

# GLO30 is an ImageCollection of tiles, so mosaic it into one image first.
elevation <- ee$ImageCollection(CFG$gee$asset_dem)$select("DEM")$mosaic()$
  rename("elevation")
slope     <- ee$Terrain$slope(elevation)$rename("slope")

s1 <- ee$ImageCollection(CFG$gee$asset_s1)$
  filterBounds(aoi_ee)$
  filter(ee$Filter$eq("instrumentMode", "IW"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))
s1 <- season_filter(s1)
s1_med <- s1$select(c("VV", "VH"))$median()$rename(c("s1_vv", "s1_vh"))

s2 <- ee$ImageCollection(CFG$gee$asset_s2)$
  filterBounds(aoi_ee)$
  filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", CFG$gee$s2_cloud_pct))
s2 <- season_filter(s2)
s2_med <- s2$median()
ndvi <- s2_med$normalizedDifference(c("B8", "B4"))$rename("s2_ndvi")

water <- ee$Image(CFG$gee$asset_jrc_water)$select("occurrence")$
  rename("water_occurrence")$unmask(0)

# Earth Engine refuses to export a multi-band image whose bands have
# different data types, and these do: the DEM and slope come back Float32
# while median() and normalizedDifference() produce Float64, and JRC
# occurrence is uint8. Casting the whole stack to Float32 with $toFloat()
# makes them uniform. Do this on ANY multi-band export you add later.
predictors <- elevation$addBands(slope)$addBands(s1_med)$
  addBands(ndvi)$addBands(water)$clip(aoi_ee)$toFloat()

msg("covariate stack: ", paste(predictors$bandNames()$getInfo(), collapse = ", "))

#' Download an ee$Image to a local GeoTIFF.
#'
#' rgee renamed ee_as_raster -> ee_as_rast (terra-based) and will drop the old
#' name in 1.2.0, so prefer the new one and fall back for older installs.
download_ee_image <- function(img, dsn) {
  fn <- if (exists("ee_as_rast")) ee_as_rast else ee_as_raster
  fn(img, region = aoi_ee, dsn = dsn, scale = CFG$gee$scale_m,
     crs = CFG$gee$export_crs, via = "drive", quiet = FALSE)
}

pred_path <- file.path(CFG$dir_current, "predictors.tif")
download_ee_image(predictors, pred_path)
msg("wrote ", pred_path)

# ---- 3. published priors ---------------------------------------------------
# Single-band, so no type-mixing possible, but cast anyway so every raster
# this pipeline writes has the same dtype.

li_path <- file.path(CFG$dir_current, "prior_li2025.tif")
download_ee_image(ee$Image(CFG$gee$asset_li2025)$clip(aoi_ee)$toFloat(), li_path)
msg("wrote ", li_path, "  (Li et al. 2025, full peat column, kg C/m2)")

sothe_path <- file.path(CFG$dir_current, "prior_sothe_0_30.tif")
download_ee_image(ee$Image(CFG$gee$asset_sothe_sc_0_30)$clip(aoi_ee)$toFloat(), sothe_path)
msg("wrote ", sothe_path, "  (Sothe et al. 2022, 0-30cm, kg C/m2)")

msg("02 complete")
