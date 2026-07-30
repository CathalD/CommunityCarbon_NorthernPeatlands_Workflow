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

# ---- 1. AOI ----------------------------------------------------------------
# Prefer the hand-drawn coast-following polygon in mvp/data/aoi.geojson. It
# follows the shoreline trend, so it contains the sampled transect without
# spending most of its area over Hudson Bay. Fall back to a buffered convex
# hull of the cores only if that file is absent.

if (file.exists(CFG$file_aoi)) {
  aoi_sf <- st_read(CFG$file_aoi, quiet = TRUE) %>% st_transform(CFG$crs_geographic)
  aoi_geom <- st_union(st_geometry(aoi_sf))
  msg("AOI: supplied polygon from ", basename(CFG$file_aoi))
} else {
  aoi_geom <- cores %>%
    st_transform(CFG$crs_equal_area) %>%
    st_union() %>% st_convex_hull() %>%
    st_buffer(CFG$gee$aoi_buffer_km * 1000) %>%
    st_transform(CFG$crs_geographic)
  msg("AOI: no ", basename(CFG$file_aoi), " found, fell back to core hull + ",
     CFG$gee$aoi_buffer_km, " km buffer")
}

# Every core must be inside the AOI, or its covariates come back NA and the
# core silently drops out of the model at step 4.
outside <- cores$core_id[!apply(st_intersects(cores, aoi_geom, sparse = FALSE), 1, any)]
if (length(outside)) {
  stop("These cores fall OUTSIDE the AOI: ", paste(outside, collapse = ", "),
      ". Their covariates would all be NA. Widen mvp/data/aoi.geojson.")
}
area_km2 <- as.numeric(st_area(st_transform(aoi_geom, CFG$crs_equal_area))) / 1e6
msg(sprintf("AOI area %.0f km2; all %d cores inside", area_km2, nrow(cores)))

st_write(st_sf(geometry = aoi_geom), file.path(CFG$dir_current, "aoi.geojson"),
         delete_dsn = TRUE, quiet = TRUE)

aoi_ee <- sf_as_ee(st_sf(geometry = aoi_geom))$geometry()

# ---- 2. covariate stack ----------------------------------------------------

season_filter <- function(coll) {
  coll$
    filter(ee$Filter$calendarRange(CFG$gee$season_start_doy,
                                   CFG$gee$season_end_doy, "day_of_year"))$
    filter(ee$Filter$calendarRange(min(CFG$gee$season_years),
                                   max(CFG$gee$season_years), "year"))
}

# --- terrain ---------------------------------------------------------------
# Copernicus GLO30 (an ImageCollection of tiles) as the base, ArcticDEM
# mosaicked ON TOP so it wins where it exists and Copernicus fills its gaps.
# ArcticDEM alone leaves a masked band here, which is what cost the first run
# both elevation and slope.
cop <- ee$ImageCollection(CFG$gee$asset_dem_copernicus)$select("DEM")$
  mosaic()$rename("elevation")
arctic <- ee$Image(CFG$gee$asset_dem_arctic)$select("elevation")
dem_native <- ee$ImageCollection(list(cop, arctic))$mosaic()$rename("elevation")

# PIN THE WORKING RESOLUTION BEFORE ANY TERRAIN OPERATION. This line is
# load-bearing, not tidying. Earth Engine runs slope and focal operations at
# the image's native scale unless told otherwise. ArcticDEM is 2 m, so a 2 km
# focal radius on the native grid asks for a kernel a thousand pixels across
# and Earth Engine refuses the request outright ("Reprojection output too
# large"). Reprojecting to the analysis scale first makes the kernels the size
# their names claim, and costs nothing scientifically here: this landscape has
# metres of relief over tens of kilometres, so 2 m detail is not the signal.
dem <- dem_native$resample("bilinear")$
  reproject(crs = CFG$gee$export_crs, scale = CFG$gee$scale_m)

slope <- ee$Terrain$slope(dem)$rename("slope_deg")

# Topographic position: elevation minus the local mean, at two scales. On a
# plain this flat, POSITION carries the peat-vs-mineral signal that absolute
# elevation does not.
tpi <- function(img, radius_m, name) {
  k <- ee$Kernel$circle(radius = radius_m, units = "meters")
  img$subtract(img$focalMean(kernel = k))$rename(name)
}
tpi_small <- tpi(dem, CFG$gee$tpi_radii_m[["small"]], "tpi_300m")
tpi_large <- tpi(dem, CFG$gee$tpi_radii_m[["large"]], "tpi_2km")

# --- distance from the coast -----------------------------------------------
# Isostatic rebound makes distance inland a proxy for time since the sea left,
# and therefore for how long peat has had to accumulate. Built from the JRC
# water layer so it needs no extra asset.
gsw <- ee$Image(CFG$gee$asset_jrc_water)
sea <- gsw$select("occurrence")$gt(90)$selfMask()   # permanent water
dist_coast <- sea$fastDistanceTransform(4096)$sqrt()$
  multiply(ee$Image$pixelArea()$sqrt())$
  rename("dist_coast_m")

# --- Sentinel-1 ------------------------------------------------------------
# Require BOTH polarisations: filtering on VV alone can admit scenes with no
# VH band, and median() over a partly-absent band gives a masked result.
s1 <- ee$ImageCollection(CFG$gee$asset_s1)$
  filterBounds(aoi_ee)$
  filter(ee$Filter$eq("instrumentMode", "IW"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))
s1 <- season_filter(s1)
s1_med <- s1$select(c("VV", "VH"))$median()$rename(c("s1_vv", "s1_vh"))

s2 <- ee$ImageCollection(CFG$gee$asset_s2)$
  filterBounds(aoi_ee)$
  filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", CFG$gee$s2_cloud_pct))
s2 <- season_filter(s2)
s2_med <- s2$median()
ndvi <- s2_med$normalizedDifference(c("B8", "B4"))$rename("s2_ndvi")

water <- gsw$select("occurrence")$unmask(0)$rename("water_occurrence")

# Earth Engine refuses to export a multi-band image whose bands have
# different data types, and these do: the DEM and slope come back Float32
# while median() and normalizedDifference() produce Float64, and JRC
# occurrence is uint8. Casting the whole stack to Float32 with $toFloat()
# makes them uniform. Do this on ANY multi-band export you add later.
predictors <- dem$
  addBands(slope)$
  addBands(tpi_small)$
  addBands(tpi_large)$
  addBands(dist_coast)$
  addBands(s1_med)$
  addBands(ndvi)$
  addBands(water)$
  clip(aoi_ee)$toFloat()

band_names <- predictors$bandNames()$getInfo()
msg(length(band_names), "-band covariate stack: ", paste(band_names, collapse = ", "))

# Sample the stack AT THE CORES before spending minutes on the export. A band
# that comes back NA here has no coverage over this AOI, and finding that out
# now is far cheaper than finding it out in step 3.
check <- tryCatch(
  ee_extract(predictors, sf_as_ee(cores["core_id"]), scale = CFG$gee$scale_m,
            fun = ee$Reducer$first(), sf = FALSE),
  error = function(e) { msg("pre-export check skipped: ", conditionMessage(e)); NULL })

if (!is.null(check)) {
  band_cols <- intersect(band_names, names(check))
  na_bands <- band_cols[vapply(check[band_cols], function(v) all(is.na(v)), logical(1))]
  if (length(na_bands)) {
    msg("WARNING: these bands are NA at EVERY core: ", paste(na_bands, collapse = ", "))
    msg("  They have no coverage over this AOI. Fix the asset in config.R before ",
       "continuing -- step 4 would drop them and you would be modelling on less ",
       "than you think.")
  } else {
    msg("pre-export check: all ", length(band_cols),
       " bands have values at all ", nrow(cores), " cores")
  }
}

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
