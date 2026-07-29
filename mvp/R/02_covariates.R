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
  .args <- commandArgs(trailingOnly = FALSE)
  .this_file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  .this_dir <- if (length(.this_file)) dirname(.this_file[[1]]) else getwd()
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

elevation <- ee$Image(CFG$gee$asset_dem)$select(0)$rename("elevation")
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

predictors <- elevation$addBands(slope)$addBands(s1_med)$
  addBands(ndvi)$addBands(water)$clip(aoi_ee)

msg("covariate stack: ", paste(predictors$bandNames()$getInfo(), collapse = ", "))

pred_path <- file.path(CFG$dir_current, "predictors.tif")
ee_as_raster(predictors, region = aoi_ee, dsn = pred_path,
            scale = CFG$gee$scale_m, crs = CFG$gee$export_crs,
            via = "drive", quiet = FALSE)
msg("wrote ", pred_path)

# ---- 3. published priors ---------------------------------------------------

li_path <- file.path(CFG$dir_current, "prior_li2025.tif")
ee_as_raster(ee$Image(CFG$gee$asset_li2025)$clip(aoi_ee), region = aoi_ee,
            dsn = li_path, scale = CFG$gee$scale_m, crs = CFG$gee$export_crs,
            via = "drive", quiet = FALSE)
msg("wrote ", li_path, "  (Li et al. 2025, full peat column, kg C/m2)")

sothe_path <- file.path(CFG$dir_current, "prior_sothe_0_30.tif")
ee_as_raster(ee$Image(CFG$gee$asset_sothe_sc_0_30)$clip(aoi_ee), region = aoi_ee,
            dsn = sothe_path, scale = CFG$gee$scale_m, crs = CFG$gee$export_crs,
            via = "drive", quiet = FALSE)
msg("wrote ", sothe_path, "  (Sothe et al. 2022, 0-30cm, kg C/m2)")

msg("02 complete")
