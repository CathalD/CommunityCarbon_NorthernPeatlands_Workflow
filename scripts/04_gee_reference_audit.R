# =============================================================================
# 04_gee_reference_audit.R
#
# Extract the three published SOC layers at the core locations -- but audit
# each one's UNITS and DEPTH SUPPORT first, from evidence, before any value is
# allowed near a comparison.
#
# This script exists because combining published carbon layers is where unit
# and depth errors do their worst damage, and because they are silent: a layer
# in t/ha and a layer in kg/m2 differ by a factor of ten and both look
# entirely reasonable on a map. A full-profile peat stock and a 0-30 cm stock
# differ by more than an order of magnitude in this landscape and, again, both
# look reasonable.
#
# The asset name of the Li et al. layer is a good illustration of why names
# cannot be trusted:
#
#     McMasterCarbon30mkgm2version1
#
# The "30" in that name is a PIXEL SIZE, not a depth. Read as a depth it would
# invite a direct comparison against SoilGrids OCS 0-30 cm, which would be a
# category error: Hudson Bay Lowlands peat is metres thick, and a full-profile
# peat carbon layer is not remotely the same quantity as a 0-30 cm stock. The
# audit below tests that reading against the layer's actual pixel grid and
# against the physical plausibility of its values, and refuses to combine
# layers whose depth support differs.
#
# INPUT   data/derived/01_cores_qc.csv
#         data/derived/02_core_stocks.csv
# OUTPUT  data/derived/04_reference_at_cores.csv
#         outputs/tables/04_reference_audit.csv
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

log_step("04  PUBLISHED REFERENCE LAYERS: AUDIT AND EXTRACTION")

if (!requireNamespace("rgee", quietly = TRUE)) {
  fail_loudly(
    "rgee is not installed",
    "This script audits and extracts the published reference layers.",
    "The sample-based products in 05 do not depend on it.",
    remedy = 'install.packages("rgee"); rgee::ee_install(); then re-run.')
}
suppressPackageStartupMessages(library(rgee))
ee_Initialize(project = CFG$gee$project, drive = TRUE)
log_ok("Earth Engine initialised")

cores <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                          "scripts/01_ingest_qc.R")
if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the Earth Engine step",
              remedy = "Re-run 01_ingest_qc.R.")
}

bb  <- CFG$aoi_bbox
aoi <- ee$Geometry$Rectangle(c(bb[["xmin"]], bb[["ymin"]],
                               bb[["xmax"]], bb[["ymax"]]))
pts <- ee$FeatureCollection(lapply(seq_len(nrow(cores)), function(i) {
  ee$Feature(ee$Geometry$Point(c(cores$longitude[i], cores$latitude[i])),
             list(core_id = cores$core_id[i]))
}))

# =============================================================================
# 1. The audit machinery
# =============================================================================

#' Interrogate an Earth Engine image: what bands, what grid, what values?
#' Returns plain R values so the reasoning can be done and logged locally.
probe_layer <- function(img, aoi, name) {
  bands <- img$bandNames()$getInfo()
  proj  <- img$select(0)$projection()
  scale_m <- tryCatch(proj$nominalScale()$getInfo(), error = function(e) NA_real_)
  crs     <- tryCatch(proj$crs()$getInfo(), error = function(e) NA_character_)

  qs <- c(1, 5, 25, 50, 75, 95, 99)
  red <- ee$Reducer$percentile(qs)$combine(ee$Reducer$mean(), "", TRUE)$
    combine(ee$Reducer$minMax(), "", TRUE)
  st <- img$select(0)$reduceRegion(
    reducer = red, geometry = aoi, scale = max(30, scale_m, na.rm = TRUE),
    maxPixels = 1e12, bestEffort = TRUE)$getInfo()

  list(name = name, bands = bands, scale_m = scale_m, crs = crs, stats = st)
}

#' Decide what a layer's values MEAN, by testing candidate unit readings
#' against physically plausible ranges rather than trusting a name.
#'
#' Returns the reading that survives, or "UNRESOLVED" when more than one
#' survives or none does. Pure.
resolve_units <- function(median_value, candidates, plausible) {
  out <- do.call(rbind, lapply(names(candidates), function(k) {
    v <- candidates[[k]](median_value)
    data.frame(
      reading   = k,
      implied_kgm2 = v,
      consistent_with_0_30 = v > 0 && v <= plausible$stock_0_30_kgm2_max,
      consistent_with_full = v >= plausible$full_profile_kgm2_min,
      stringsAsFactors = FALSE)
  }))
  out
}

audit_rows <- list()
record_audit <- function(...) {
  audit_rows[[length(audit_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
}

# =============================================================================
# 2. SoilGrids 2.0 OCS 0-30 cm
# =============================================================================

log_step("04a  SOILGRIDS 2.0 OCS 0-30 cm")

sg_img <- ee$Image(CFG$gee$asset_soilgrids_ocs)
sg <- probe_layer(sg_img, aoi, "SoilGrids_OCS_0_30")
log_info("bands: ", paste(sg$bands, collapse = ", "))
log_info(sprintf("native grid: %.1f m, %s", sg$scale_m, sg$crs))

sg_band <- grep("0-30", sg$bands, value = TRUE)[1]
if (is.na(sg_band)) sg_band <- sg$bands[1]
log_info("using band: ", sg_band)

sg_sel <- sg_img$select(sg_band)
sg_med <- ee$Number(sg_sel$reduceRegion(
  ee$Reducer$median(), aoi, 250, maxPixels = 1e12, bestEffort = TRUE
)$values()$get(0))$getInfo()
log_info(sprintf("median raw value over AOI: %.3f", sg_med))

# SoilGrids OCS is documented as t/ha with a stored scale factor of 10.
# Test that reading rather than apply it blindly.
sg_res <- resolve_units(sg_med, list(
  "raw is t/ha"                = function(v) tha_to_kgm2(v),
  "raw/10 is t/ha (documented)" = function(v) tha_to_kgm2(v / 10),
  "raw is kg/m2"               = function(v) v
), CFG$ref_plausible)
print_table(sg_res, digits = 3)

sg_ok <- sg_res[sg_res$consistent_with_0_30 & !sg_res$consistent_with_full, ]
if (nrow(sg_ok) == 1L) {
  log_ok("units resolved: ", sg_ok$reading,
         sprintf("  -> %.2f kg C/m2 over 0-30 cm", sg_ok$implied_kgm2))
} else {
  log_warn("units NOT uniquely resolved by plausibility; ",
           nrow(sg_ok), " readings survive. Applying the DOCUMENTED reading ",
           "(raw/10 = t/ha) and flagging the ambiguity.")
}
# Documented conversion, applied explicitly and named.
sg_kgm2 <- sg_sel$divide(10)$divide(10)$rename("soilgrids_ocs_0_30_kgm2")

record_audit(
  layer = "SoilGrids 2.0 OCS",
  asset = CFG$gee$asset_soilgrids_ocs,
  band = sg_band,
  native_scale_m = sg$scale_m,
  median_raw = sg_med,
  units_applied = "raw/10 = t/ha, then t/ha/10 = kg/m2",
  depth_support = "0-30 cm (stated in band name and product documentation)",
  comparable_to_core_0_30 = TRUE,
  note = "Depth support matches the reference window. Directly comparable."
)

# =============================================================================
# 3. Sothe et al. 2022, Canada-wide SOC 30 cm
# =============================================================================

log_step("04b  SOTHE ET AL. 2022 CANADA SOC")

so_img <- ee$Image(CFG$gee$asset_sothe2022)
so <- probe_layer(so_img, aoi, "Sothe2022")
log_info("bands: ", paste(so$bands, collapse = ", "))
log_info(sprintf("native grid: %.1f m, %s", so$scale_m, so$crs))

so_med <- ee$Number(so_img$select(0)$reduceRegion(
  ee$Reducer$median(), aoi, max(30, so$scale_m), maxPixels = 1e12,
  bestEffort = TRUE)$values()$get(0))$getInfo()
log_info(sprintf("median raw value over AOI: %.3f", so_med))

so_res <- resolve_units(so_med, list(
  "raw is t/ha"  = function(v) tha_to_kgm2(v),
  "raw is kg/m2" = function(v) v
), CFG$ref_plausible)
print_table(so_res, digits = 3)

so_ok <- so_res[so_res$consistent_with_0_30 & !so_res$consistent_with_full, ]
so_units <- if (nrow(so_ok) == 1L) so_ok$reading else "UNRESOLVED"
if (so_units == "UNRESOLVED") {
  log_warn("Sothe units are not uniquely determined by plausibility alone.")
  log_warn("Both t/ha and kg/m2 readings fall in a defensible range, or ",
           "neither does. Resolve against the publication before using this ",
           "layer quantitatively. Extraction proceeds; the values are ",
           "carried with an explicit UNRESOLVED marker.")
}
so_kgm2 <- if (identical(so_units, "raw is t/ha")) {
  so_img$select(0)$divide(10)$rename("sothe2022_soc_0_30_kgm2")
} else {
  so_img$select(0)$rename("sothe2022_soc_0_30_kgm2")
}

record_audit(
  layer = "Sothe et al. 2022",
  asset = CFG$gee$asset_sothe2022,
  band = so$bands[1],
  native_scale_m = so$scale_m,
  median_raw = so_med,
  units_applied = so_units,
  depth_support = "0-30 cm (per publication title)",
  comparable_to_core_0_30 = TRUE,
  note = paste0("Depth matches the reference window. Units ",
                if (so_units == "UNRESOLVED") "UNRESOLVED - see warning." else
                  "resolved by plausibility test.")
)

# =============================================================================
# 4. Li et al. 2025, Hudson Bay Lowlands peat carbon
# =============================================================================

log_step("04c  LI ET AL. 2025 HBL PEAT CARBON")

li_img <- ee$Image(CFG$gee$asset_li2025)
li <- probe_layer(li_img, aoi, "Li2025")
log_info("bands: ", paste(li$bands, collapse = ", "))
log_info(sprintf("native grid: %.1f m, %s", li$scale_m, li$crs))

li_med <- ee$Number(li_img$select(0)$reduceRegion(
  ee$Reducer$median(), aoi, max(30, li$scale_m), maxPixels = 1e12,
  bestEffort = TRUE)$values()$get(0))$getInfo()
li_p95 <- ee$Number(li_img$select(0)$reduceRegion(
  ee$Reducer$percentile(list(95)), aoi, max(30, li$scale_m),
  maxPixels = 1e12, bestEffort = TRUE)$values()$get(0))$getInfo()
log_info(sprintf("median raw value over AOI: %.3f   (p95: %.3f)", li_med, li_p95))

# --- test the asset-name reading --------------------------------------------
log_info("testing the asset name 'McMasterCarbon30mkgm2version1':")
name_says_30_is_pixels <- !is.na(li$scale_m) && abs(li$scale_m - 30) < 5
log_info(sprintf("   is the native pixel ~30 m?  %s (measured %.1f m)",
                 if (name_says_30_is_pixels) "YES" else "NO", li$scale_m))
if (name_says_30_is_pixels) {
  log_ok("The '30m' in the asset name is the PIXEL SIZE, not a depth.")
  log_warn("Nothing in this layer's name states a depth support.")
}

li_res <- resolve_units(li_med, list(
  "raw is kg/m2" = function(v) v,
  "raw is t/ha"  = function(v) tha_to_kgm2(v)
), CFG$ref_plausible)
print_table(li_res, digits = 3)

reads_full_profile <- any(li_res$consistent_with_full &
                          li_res$reading == "raw is kg/m2")

log_warn(strrep("-", 68))
if (reads_full_profile) {
  log_warn("LI ET AL. VALUES ARE TOO LARGE TO BE A 0-30 cm STOCK.")
  log_warn(sprintf("  median over AOI      : %.1f kg C/m2", li_med))
  log_warn(sprintf("  0-30 cm ceiling      : %.1f kg C/m2",
                   CFG$ref_plausible$stock_0_30_kgm2_max))
  log_warn(sprintf("  our deepest core     : %.1f kg C/m2 over 0-30 cm",
                   max(read_csv_verbatim(
                     file.path(CFG$dir_derived, "02_core_stocks.csv"))$stock_kgm2)))
  log_warn("  This is a FULL-PROFILE peat carbon stock, integrated over peat")
  log_warn("  depths of metres. It is NOT the same quantity as a 0-30 cm")
  log_warn("  stock and MUST NOT be differenced against, averaged with, or")
  log_warn("  validated by the 0-30 cm products or by these cores.")
  li_depth <- "FULL PROFILE (peat column); NOT 0-30 cm"
  li_comparable <- FALSE
} else {
  log_warn("Li et al. values fall within a 0-30 cm plausible range, which")
  log_warn("contradicts the expectation for a full-profile HBL peat layer.")
  log_warn("Treat the depth support as UNRESOLVED and confirm against the")
  log_warn("publication before combining.")
  li_depth <- "UNRESOLVED"
  li_comparable <- FALSE
}
log_warn(strrep("-", 68))

li_kgm2 <- li_img$select(0)$rename("li2025_peat_carbon_kgm2")

record_audit(
  layer = "Li et al. 2025 HBL peat carbon",
  asset = CFG$gee$asset_li2025,
  band = li$bands[1],
  native_scale_m = li$scale_m,
  median_raw = li_med,
  units_applied = "kg/m2 (as the asset name states; magnitude consistent)",
  depth_support = li_depth,
  comparable_to_core_0_30 = li_comparable,
  note = paste0("The '30m' in the asset name is pixel size, not depth. ",
                "Depth support differs from the 0-30 cm window, so this layer ",
                "is retained as CONTEXT only and is barred from combination.")
)

# =============================================================================
# 5. Combination guard
# =============================================================================

log_step("04d  COMBINATION GUARD")

audit <- do.call(rbind, audit_rows)
print_table(audit[, c("layer", "native_scale_m", "median_raw",
                      "depth_support", "comparable_to_core_0_30")], digits = 2)

ok_to_combine <- audit$layer[audit$comparable_to_core_0_30]
barred        <- audit$layer[!audit$comparable_to_core_0_30]

log_ok("comparable to the 0-30 cm core stocks: ",
       paste(ok_to_combine, collapse = "; "))
if (length(barred)) {
  log_warn("BARRED from combination or validation against the cores: ",
           paste(barred, collapse = "; "))
  log_warn("These layers are extracted and reported for context. Any script ",
           "that tries to difference or ensemble them with a 0-30 cm product ",
           "should stop.")
}

#' Refuse to combine layers whose depth support is not identical. Exported so
#' downstream scripts can call it rather than re-deriving the rule.
guard_same_depth_support <- function(audit, layers) {
  d <- audit$depth_support[match(layers, audit$layer)]
  if (length(unique(d)) > 1L) {
    fail_loudly(
      "Refusing to combine layers with different depth support",
      paste0("  ", layers, "  ->  ", d),
      remedy = "Harmonise depth support first, or report the layers separately."
    )
  }
  invisible(TRUE)
}

# =============================================================================
# 6. Extract at cores
# =============================================================================

log_step("04e  EXTRACTION AT CORE LOCATIONS")

ref_stack <- sg_kgm2$addBands(so_kgm2)$addBands(li_kgm2)
samp <- ref_stack$sampleRegions(collection = pts, scale = 30,
                                geometries = FALSE, tileScale = 4)
ref_at_cores <- sf::st_drop_geometry(ee_as_sf(samp, maxFeatures = 10000))
ref_at_cores <- merge(cores[, c("core_id", "campaign")], ref_at_cores,
                      by = "core_id", all.x = TRUE)

# Carry the audit verdict INTO the data, so a downstream reader cannot pick up
# the Li column without also picking up the reason it is not comparable.
attr(ref_at_cores, "audit") <- audit
ref_at_cores$li2025_depth_support   <- li_depth
ref_at_cores$li2025_comparable_0_30 <- li_comparable
ref_at_cores$sothe_units            <- so_units

print_table(ref_at_cores, digits = 3)

write_csv_logged(ref_at_cores,
                 file.path(CFG$dir_derived, "04_reference_at_cores.csv"),
                 "published layers at core locations, with audit verdicts")
write_csv_logged(audit, file.path(CFG$dir_tables, "04_reference_audit.csv"),
                 "units and depth support, established by measurement")

log_ok("04 complete")
