# =============================================================================
# 03_gee_covariates.R
#
# Build a covariate stack for the study area in Google Earth Engine and extract
# it at the eight core locations.
#
# REQUIRES rgee and an authenticated Earth Engine account. This is the first
# script in the workflow with an external dependency; everything upstream runs
# on base R alone.
#
# ---------------------------------------------------------------------------
# COVARIATE SELECTION, AND THE REASONING BEHIND IT
#
# The Hudson Bay Lowlands are the second-largest peatland complex on earth: a
# nearly flat plain of marine clay left by the Tyrrell Sea, still rebounding
# isostatically, with peat accumulating on it since emergence. That history
# dictates which covariates can plausibly carry information about soil carbon
# here, and which are decoration.
#
#   1. RELATIVE ELEVATION / MICROTOPOGRAPHY               [strong]
#      On a plain with metres of relief over tens of kilometres, centimetres of
#      relative elevation decide whether ground is saturated peat or drained
#      mineral soil. This is the dominant local control. ArcticDEM is used in
#      preference to SRTM or NASADEM because it is built for high latitudes and
#      resolves the low ridges and flarks that matter. Multi-scale topographic
#      position is more useful here than absolute elevation.
#
#   2. DISTANCE FROM THE HUDSON BAY COAST                 [strong]
#      Isostatic rebound makes distance inland a proxy for time since marine
#      emergence, and therefore for how long peat has had to accumulate. This
#      is the clearest large-scale gradient in the region and it is nearly free
#      to compute. It is a chronosequence proxy, not a causal variable.
#
#   3. SENTINEL-1 BACKSCATTER, GROWING SEASON             [strong]
#      C-band VV/VH responds to surface soil moisture and to standing water
#      beneath low vegetation. It separates saturated peat surfaces from
#      drained mineral ground in a way optical data cannot, and it works
#      through the persistent cloud of the Hudson Bay coast.
#
#   4. SENTINEL-2 GROWING-SEASON COMPOSITE + INDICES      [moderate]
#      Vegetation community tracks substrate closely in this landscape: sedge
#      fen, Sphagnum bog and black spruce upland are optically distinct. NDVI,
#      NDWI and a wetness index are the useful summaries. Treated as a proxy
#      for the surface, not for what lies beneath it.
#
#   5. JRC GLOBAL SURFACE WATER OCCURRENCE                [moderate]
#      Frequency of standing water maps the pool-and-ridge structure of the
#      peat complexes directly.
#
#   6. ESA WORLDCOVER                                     [used for strata]
#      Not a model input. Used to define the mapping strata for the
#      design-based product in 05, which is the role a categorical layer can
#      honestly play with eight cores.
#
# DELIBERATELY EXCLUDED
#
#   TERRACLIMATE / ERA5 CLIMATE NORMALS
#      Across a study area roughly 35 km on a side, gridded climate normals are
#      very nearly constant. They would enter a model as a near-constant, which
#      cannot explain within-area variation but can absorb degrees of freedom
#      and inflate apparent fit. Excluded on those grounds, not because climate
#      does not matter for peat: it matters enormously, at scales larger than
#      this map.
#
#   ALPHAEARTH ANNUAL EMBEDDINGS AS MODEL PREDICTORS
#      The 64-band embedding is a genuinely powerful representation, and with
#      eight observations it is exactly the wrong tool: 64 predictors on 8
#      points fits any target perfectly and generalises to nothing. It is
#      therefore NOT offered to the model. It is extracted and used instead for
#      an unsupervised task where n = 8 is not the binding constraint --
#      delineating landscape strata, and measuring how far each core sits from
#      the rest of the study area in embedding space, which is a defensible way
#      to quantify where a map is extrapolating. See 06.
#
# ---------------------------------------------------------------------------
# INPUT   data/derived/01_cores_qc.csv
# OUTPUT  data/derived/03_covariates_at_cores.csv
#         outputs/gee/03_covariate_dictionary.csv
#         plus GEE export tasks for the raster stack (to Google Drive)
# =============================================================================

## --- bootstrap ---------------------------------------------------------------
.root <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  d <- if (length(f)) dirname(f[[1]]) else getwd()
  while (!file.exists(file.path(d, "config.R"))) {
    p <- dirname(d)
    if (identical(p, d)) stop("config.R not found above '", getwd(), "'", call. = FALSE)
    d <- p
  }
  normalizePath(d)
})
source(file.path(.root, "config.R"))

log_step("03  EARTH ENGINE COVARIATES")

# =============================================================================
# 0. Dependency gate
# =============================================================================

if (!requireNamespace("rgee", quietly = TRUE)) {
  fail_loudly(
    "rgee is not installed",
    "Scripts 03 and 04 are the only parts of this workflow that need Earth",
    "Engine. Everything upstream (01, 02) and the design-based products (05)",
    "run on base R alone and do not require this.",
    remedy = paste0(
      'install.packages("rgee"); rgee::ee_install(); rgee::ee_Initialize(',
      'user = "<you>", drive = TRUE)   then re-run this script.')
  )
}
suppressPackageStartupMessages(library(rgee))

log_info("initialising Earth Engine for project: ", CFG$gee$project)
ee_Initialize(project = CFG$gee$project, drive = TRUE)
log_ok("Earth Engine initialised")

cores <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                          "scripts/01_ingest_qc.R")

# Re-assert the geolocation invariant here. This script is the first point at
# which a sign error would start silently sampling the wrong continent, so the
# check is repeated rather than trusted from upstream.
if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the Earth Engine step",
              "Cores would be sampled in the eastern hemisphere.",
              remedy = "Re-run 01_ingest_qc.R and inspect the geolocation gate.")
}
log_ok("geolocation invariant re-checked at the Earth Engine boundary")

# =============================================================================
# 1. Geometry
# =============================================================================

bb  <- CFG$aoi_bbox
aoi <- ee$Geometry$Rectangle(c(bb[["xmin"]], bb[["ymin"]],
                               bb[["xmax"]], bb[["ymax"]]))
log_info(sprintf("AOI: %.3f..%.3f E, %.3f..%.3f N",
                 bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]]))

pts <- ee$FeatureCollection(lapply(seq_len(nrow(cores)), function(i) {
  ee$Feature(
    ee$Geometry$Point(c(cores$longitude[i], cores$latitude[i])),
    list(core_id = cores$core_id[i], campaign = cores$campaign[i])
  )
}))

# =============================================================================
# 2. Terrain
# =============================================================================

log_step("03a  TERRAIN")

#' ArcticDEM where it exists, Copernicus GLO-30 as the fallback. ArcticDEM has
#' gaps; silently returning an empty band would be worse than saying so.
build_terrain <- function(aoi) {
  arctic <- ee$Image(CFG$gee$asset_dem_arctic)$select("elevation")
  cop    <- ee$ImageCollection(CFG$gee$asset_dem_copernicus)$
              select("DEM")$mosaic()$rename("elevation")
  # Prefer ArcticDEM, fill holes from Copernicus.
  dem <- ee$ImageCollection(list(cop, arctic))$mosaic()$rename("elevation")

  # Multi-scale topographic position. In a landscape this flat the SCALE at
  # which position is measured determines what it means: at 300 m it captures
  # ridge-and-flark microtopography, at 2 km it captures the broad peat plateaus.
  tpi <- function(img, radius_m, name) {
    k <- ee$Kernel$circle(radius = radius_m, units = "meters")
    img$subtract(img$focalMean(kernel = k))$rename(name)
  }

  slope <- ee$Terrain$slope(dem)$rename("slope_deg")

  dem$
    addBands(tpi(dem, 300,  "tpi_300m"))$
    addBands(tpi(dem, 2000, "tpi_2km"))$
    addBands(slope)
}

terrain <- build_terrain(aoi)
log_ok("terrain: elevation, tpi_300m, tpi_2km, slope_deg")

# --- distance from the coast -------------------------------------------------
# The rebound chronosequence proxy. Built from the JRC water layer rather than
# a coastline vector so it needs no extra asset.
log_info("building distance-from-coast (isostatic rebound proxy)")
gsw <- ee$Image(CFG$gee$asset_jrc_water)
# Permanent water: occurrence above 90% of observations.
sea <- gsw$select("occurrence")$gt(90)$selfMask()
dist_coast <- sea$fastDistanceTransform(4096)$sqrt()$
  multiply(ee$Image$pixelArea()$sqrt())$
  rename("dist_coast_m")$
  clip(aoi)

water_occurrence <- gsw$select("occurrence")$unmask(0)$rename("water_occurrence")
log_ok("hydrology: dist_coast_m, water_occurrence")

# =============================================================================
# 3. Sentinel-1
# =============================================================================

log_step("03b  SENTINEL-1")

build_s1 <- function(aoi) {
  s1 <- ee$ImageCollection(CFG$gee$asset_s1)$
    filterBounds(aoi)$
    filter(ee$Filter$eq("instrumentMode", "IW"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
    filter(ee$Filter$eq("orbitProperties_pass", "DESCENDING"))$
    filter(ee$Filter$calendarRange(CFG$gee$season_start_doy,
                                   CFG$gee$season_end_doy, "day_of_year"))$
    filter(ee$Filter$calendarRange(min(CFG$gee$season_years),
                                   max(CFG$gee$season_years), "year"))$
    select(c("VV", "VH"))

  med <- s1$median()
  # The VV/VH ratio separates surface-scattering wet ground from
  # volume-scattering vegetation better than either band alone.
  med$addBands(med$select("VV")$subtract(med$select("VH"))$rename("VV_VH_diff"))$
    rename(c("s1_vv", "s1_vh", "s1_vv_vh_diff"))$
    clip(aoi)
}

s1 <- build_s1(aoi)
log_ok("sentinel-1: s1_vv, s1_vh, s1_vv_vh_diff (descending, JUL-AUG median)")

# =============================================================================
# 4. Sentinel-2
# =============================================================================

log_step("03c  SENTINEL-2")

build_s2 <- function(aoi) {
  mask_s2 <- function(img) {
    scl <- img$select("SCL")
    # Drop saturated, shadow, cloud (medium+high), cirrus and snow.
    bad <- scl$eq(1)$Or(scl$eq(3))$Or(scl$eq(8))$Or(scl$eq(9))$
      Or(scl$eq(10))$Or(scl$eq(11))
    img$updateMask(bad$Not())$divide(10000)
  }

  s2 <- ee$ImageCollection(CFG$gee$asset_s2)$
    filterBounds(aoi)$
    filter(ee$Filter$calendarRange(CFG$gee$season_start_doy,
                                   CFG$gee$season_end_doy, "day_of_year"))$
    filter(ee$Filter$calendarRange(min(CFG$gee$season_years),
                                   max(CFG$gee$season_years), "year"))$
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", CFG$gee$s2_cloud_pct))$
    map(mask_s2)

  m <- s2$median()

  ndvi <- m$normalizedDifference(c("B8", "B4"))$rename("ndvi")
  ndwi <- m$normalizedDifference(c("B3", "B8"))$rename("ndwi")   # McFeeters
  ndmi <- m$normalizedDifference(c("B8", "B11"))$rename("ndmi")  # moisture
  # Tasseled-cap wetness (Sentinel-2 coefficients, Shi & Xu 2019). Responds to
  # soil and canopy moisture together, which is the relevant quantity here.
  tcw <- m$expression(
    "0.1509*B + 0.1973*G + 0.3279*R + 0.3406*N - 0.7112*S1 - 0.4572*S2",
    list(B = m$select("B2"), G = m$select("B3"), R = m$select("B4"),
         N = m$select("B8"), S1 = m$select("B11"), S2 = m$select("B12"))
  )$rename("tc_wetness")

  m$select(c("B2", "B3", "B4", "B8", "B11", "B12"))$
    rename(c("s2_blue", "s2_green", "s2_red", "s2_nir", "s2_swir1", "s2_swir2"))$
    addBands(ndvi)$addBands(ndwi)$addBands(ndmi)$addBands(tcw)$
    clip(aoi)
}

s2 <- build_s2(aoi)
log_ok("sentinel-2: 6 bands + ndvi, ndwi, ndmi, tc_wetness (JUL-AUG median)")

# =============================================================================
# 5. Land cover (strata) and AlphaEarth (NOT a model predictor)
# =============================================================================

log_step("03d  LAND COVER AND EMBEDDINGS")

worldcover <- ee$ImageCollection(CFG$gee$asset_worldcover)$first()$
  select("Map")$rename("worldcover")$clip(aoi)
log_ok("worldcover: used to define mapping strata in 05, not as a predictor")

alphaearth <- ee$ImageCollection(CFG$gee$asset_alphaearth)$
  filterBounds(aoi)$
  filter(ee$Filter$calendarRange(max(CFG$gee$season_years),
                                 max(CFG$gee$season_years), "year"))$
  mosaic()$clip(aoi)
log_ok("alphaearth: 64-band embedding extracted for stratification and ")
log_info("            extrapolation diagnostics only. With n=8 it is NOT ")
log_info("            offered to any predictive model. See the header.")

# =============================================================================
# 6. Assemble, extract, export
# =============================================================================

log_step("03e  EXTRACTION AT CORE LOCATIONS")

# Predictor stack: what the model in 06 may legitimately see.
predictors <- terrain$
  addBands(dist_coast)$
  addBands(water_occurrence)$
  addBands(s1)$
  addBands(s2)$
  clip(aoi)

# Full stack adds the categorical and embedding layers, which are extracted but
# ring-fenced from the model.
stack <- predictors$addBands(worldcover)$addBands(alphaearth)

log_info("sampling ", nrow(cores), " core locations at ",
         CFG$gee$export_scale_m, " m")
samp <- stack$sampleRegions(
  collection = pts,
  scale      = CFG$gee$export_scale_m,
  geometries = FALSE,
  tileScale  = 4
)

cov_at_cores <- ee_as_sf(samp, maxFeatures = 10000)
cov_at_cores <- sf::st_drop_geometry(cov_at_cores)

if (nrow(cov_at_cores) != nrow(cores)) {
  log_warn("extraction returned ", nrow(cov_at_cores), " rows for ",
           nrow(cores), " cores. A missing row means the core fell on a ",
           "masked pixel in at least one layer; sampleRegions drops such ",
           "features entirely rather than returning NA.")
  log_warn("missing: ",
           paste(setdiff(cores$core_id, cov_at_cores$core_id), collapse = ", "))
}

cov_at_cores <- merge(
  cores[, c("core_id", "campaign", "latitude", "longitude")],
  cov_at_cores, by = c("core_id", "campaign"), all.x = TRUE)

write_csv_logged(cov_at_cores,
                 file.path(CFG$dir_derived, "03_covariates_at_cores.csv"),
                 "covariates sampled at the 8 core locations")

# --- dictionary --------------------------------------------------------------
dict <- data.frame(
  covariate = c("elevation", "tpi_300m", "tpi_2km", "slope_deg",
                "dist_coast_m", "water_occurrence",
                "s1_vv", "s1_vh", "s1_vv_vh_diff",
                "s2_blue", "s2_green", "s2_red", "s2_nir", "s2_swir1",
                "s2_swir2", "ndvi", "ndwi", "ndmi", "tc_wetness",
                "worldcover", "A00..A63"),
  role = c(rep("predictor", 19), "stratum", "diagnostic"),
  rationale = c(
    "Absolute height. Weak here on its own; the plain is nearly flat.",
    "Topographic position at 300 m: ridge-and-flark microtopography, the scale that decides saturation.",
    "Topographic position at 2 km: broad peat plateaus versus drained margins.",
    "Slope. Near zero almost everywhere; retained to identify the few drained margins.",
    "Distance from permanent water/coast. Proxy for time since marine emergence under isostatic rebound, hence for peat accumulation time.",
    "Frequency of standing water. Maps pool-and-ridge structure of peat complexes.",
    "C-band VV backscatter. Responds to surface soil moisture; penetrates cloud.",
    "C-band VH backscatter. Volume scattering from vegetation structure.",
    "VV minus VH. Separates surface-scattering wet ground from volume-scattering canopy.",
    rep("Sentinel-2 growing-season median reflectance. Vegetation community as a substrate proxy.", 6),
    "Greenness. Distinguishes productive upland forest from sparse bog.",
    "Open-water index.",
    "Moisture index from the NIR/SWIR contrast.",
    "Tasseled-cap wetness. Combined soil and canopy moisture.",
    "ESA WorldCover class. Used to define mapping strata in 05, NOT as a model predictor.",
    "AlphaEarth 64-band annual embedding. Extracted for unsupervised stratification and for measuring extrapolation distance. NEVER a predictor: 64 dimensions on 8 observations cannot be fitted honestly."
  ),
  stringsAsFactors = FALSE
)
write_csv_logged(dict, file.path(CFG$dir_gee, "03_covariate_dictionary.csv"),
                 "what each covariate is and why it is here")

# =============================================================================
# 7. Raster export for mapping
# =============================================================================

log_step("03f  RASTER EXPORT")

task_pred <- ee_image_to_drive(
  image       = predictors$toFloat(),
  description = "ccnp_predictors_30m",
  folder      = "CCNP_SOC",
  region      = aoi,
  scale       = CFG$gee$export_scale_m,
  crs         = CFG$gee$export_crs,
  maxPixels   = 1e13
)
task_pred$start()
log_ok("started export: ccnp_predictors_30m -> Drive/CCNP_SOC")

task_lc <- ee_image_to_drive(
  image       = worldcover,
  description = "ccnp_worldcover_30m",
  folder      = "CCNP_SOC",
  region      = aoi,
  scale       = CFG$gee$export_scale_m,
  crs         = CFG$gee$export_crs,
  maxPixels   = 1e13
)
task_lc$start()
log_ok("started export: ccnp_worldcover_30m -> Drive/CCNP_SOC")

log_info("monitor with rgee::ee_monitoring(), then place the GeoTIFFs in ",
         CFG$dir_gee)
log_ok("03 complete")
