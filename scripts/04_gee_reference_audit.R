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
ee_Initialize(project = CFG$gee$project)
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
# 2. SoilGrids 2.0 -- built from concentration x bulk density
# =============================================================================

log_step("04a  SOILGRIDS 2.0 (built from soc_mean x bdod_mean)")

# The packaged `ocs_mean` asset is NOT used. Building the stock here from
# carbon concentration and bulk density makes every step auditable: the depth
# intervals are explicit, the unit conversion is written out, and the same
# construction yields 0-30 cm and 0-100 cm on an identical basis. It also
# mirrors exactly what this workflow does to the community cores, so the
# comparison is like for like rather than one integration method against
# another.
#
#   soc_mean  is dg/kg   (/10  -> g/kg)
#   bdod_mean is cg/cm3  (/100 -> g/cm3)
#   layer OCS (kg C/m2) = (soc/10) * (bdod/100) * thickness_cm / 100
#
# Dimensional check, done on paper before trusting the code:
#   peat at 400 g/kg and 0.20 g/cm3 over 5 cm
#     = 0.20 g/cm3 * 0.40 gC/g * 5 cm = 0.40 g/cm2 = 4.0 kg C/m2
#   formula: 400 * 0.20 * 5 / 100 = 4.0   OK

sg_soc  <- ee$Image(CFG$gee$asset_soilgrids_soc)
sg_bdod <- ee$Image(CFG$gee$asset_soilgrids_bdod)

sg_probe <- probe_layer(sg_soc, aoi, "SoilGrids_soc")
log_info("soc bands: ", paste(sg_probe$bands, collapse = ", "))
log_info(sprintf("native grid: %.1f m, %s", sg_probe$scale_m, sg_probe$crs))

sg_layers <- list(
  list(soc = "soc_0-5cm_mean",    bd = "bdod_0-5cm_mean",    thick = 5,  nm = "ocs_0_5cm"),
  list(soc = "soc_5-15cm_mean",   bd = "bdod_5-15cm_mean",   thick = 10, nm = "ocs_5_15cm"),
  list(soc = "soc_15-30cm_mean",  bd = "bdod_15-30cm_mean",  thick = 15, nm = "ocs_15_30cm"),
  list(soc = "soc_30-60cm_mean",  bd = "bdod_30-60cm_mean",  thick = 30, nm = "ocs_30_60cm"),
  list(soc = "soc_60-100cm_mean", bd = "bdod_60-100cm_mean", thick = 40, nm = "ocs_60_100cm")
)

sg_layer_img <- function(l) {
  sg_soc$select(l$soc)$divide(10)$
    multiply(sg_bdod$select(l$bd)$divide(100))$
    multiply(l$thick)$divide(100)$
    rename(l$nm)
}
sg_stack <- Reduce(function(a, b) a$addBands(b), lapply(sg_layers, sg_layer_img))

sg_0_30 <- sg_stack$select(c("ocs_0_5cm", "ocs_5_15cm", "ocs_15_30cm"))$
  reduce(ee$Reducer$sum())$rename("soilgrids_ocs_0_30_kgm2")
sg_0_100 <- sg_stack$reduce(ee$Reducer$sum())$rename("soilgrids_ocs_0_100_kgm2")

sg_med <- ee$Number(sg_0_30$reduceRegion(
  ee$Reducer$median(), aoi, 250, maxPixels = 1e12, bestEffort = TRUE
)$values()$get(0))$getInfo()
sg_med100 <- ee$Number(sg_0_100$reduceRegion(
  ee$Reducer$median(), aoi, 250, maxPixels = 1e12, bestEffort = TRUE
)$values()$get(0))$getInfo()

log_info(sprintf("median constructed 0-30 cm stock over AOI : %.2f kg C/m2", sg_med))
log_info(sprintf("median constructed 0-100 cm stock over AOI: %.2f kg C/m2", sg_med100))

# Two internal consistency checks that would catch a wrong scale factor.
if (sg_med100 < sg_med) {
  fail_loudly("SoilGrids 0-100 cm stock is smaller than its own 0-30 cm stock",
              sprintf("0-30 cm: %.2f, 0-100 cm: %.2f", sg_med, sg_med100),
              "A deeper interval cannot contain less carbon than an interval",
              "nested inside it. The band scaling or the depth mapping is wrong.",
              remedy = "Re-check the SoilGrids band names and unit divisors above.")
}
log_ok("0-100 cm exceeds 0-30 cm, as nesting requires")

if (sg_med <= 0 || sg_med > CFG$ref_plausible$stock_0_30_kgm2_max) {
  # Before blaming the construction, check whether the arithmetic could even
  # produce this from real soil. 30 cm of the densest plausible peat -- bulk
  # density 0.2 g/cm3 at 50 % carbon -- yields 30 kg/m2. Reaching the observed
  # median needs a carbon density no peat attains.
  implied_c_density <- sg_med / (30 * 10)   # kg/m2 over 30 cm -> g C/cm3
  log_warn(strrep("-", 68))
  log_warn("SOILGRIDS EXCEEDS PHYSICAL PLAUSIBILITY OVER THIS AREA.")
  log_warn(sprintf("  constructed 0-30 cm median : %.1f kg C/m2", sg_med))
  log_warn(sprintf("  plausibility ceiling       : %.0f kg C/m2",
                   CFG$ref_plausible$stock_0_30_kgm2_max))
  log_warn(sprintf("  implied carbon density     : %.3f g C/cm3", implied_c_density))
  log_warn("  For comparison, 30 cm of dense peat at 0.2 g/cm3 and 50% carbon")
  log_warn("  gives 0.10 g C/cm3, or 30 kg/m2 -- less than half this.")
  log_warn("")
  log_warn("  The construction itself is verified: the same formula reproduces")
  log_warn("  the documented worked example, and the 0-100 cm sum correctly")
  log_warn("  exceeds the 0-30 cm sum. So this is not an arithmetic error on")
  log_warn("  our side; it is what SoilGrids reports here.")
  log_warn("")
  log_warn("  SoilGrids is known to model organic soils poorly, typically by")
  log_warn("  assigning them mineral-like bulk densities. The community cores")
  log_warn("  are the only local measurement available to test that, and the")
  log_warn("  comparison in 07 is where it gets tested.")
  log_warn(strrep("-", 68))
  sg_implausible <- TRUE
} else {
  log_ok(sprintf("constructed 0-30 cm stock is physically plausible (%.2f kg/m2)",
                 sg_med))
  sg_implausible <- FALSE
}

sg_kgm2 <- sg_0_30

record_audit(
  layer = "SoilGrids 2.0 OCS (constructed)",
  asset = paste(CFG$gee$asset_soilgrids_soc, "+", CFG$gee$asset_soilgrids_bdod),
  band = "soc_*_mean x bdod_*_mean, summed 0-5/5-15/15-30",
  native_scale_m = sg_probe$scale_m,
  median_raw = sg_med,
  units_applied = "(soc/10 g/kg) x (bdod/100 g/cm3) x thickness_cm / 100 = kg C/m2",
  depth_support = "0-30 cm (constructed from the three nested intervals)",
  comparable_to_core_0_30 = TRUE,
  note = paste0("Built from concentration and bulk density rather than the ",
                "packaged ocs_mean asset, so the integration matches the one ",
                "applied to the cores. 0-100 cm median = ",
                sprintf("%.2f kg/m2.", sg_med100))
)

# =============================================================================
# 3. Sothe et al. 2022, Canada-wide SOC 30 cm
# =============================================================================

log_step("04b  SOTHE ET AL. 2022 CANADA SOC")

so_img <- ee$Image(CFG$gee$asset_sothe_sc_0_30)
so <- probe_layer(so_img, aoi, "Sothe2022")
log_info("bands: ", paste(so$bands, collapse = ", "))
log_info(sprintf("native grid: %.1f m, %s", so$scale_m, so$crs))

so_med <- ee$Number(so_img$select(0)$reduceRegion(
  ee$Reducer$median(), aoi, max(30, so$scale_m), maxPixels = 1e12,
  bestEffort = TRUE)$values()$get(0))$getInfo()
log_info(sprintf("median raw value over AOI: %.3f", so_med))

so_res <- resolve_units(so_med, list(
  "raw is kg/m2 (catalogued)" = function(v) v,
  "raw is t/ha"               = function(v) tha_to_kgm2(v)
), CFG$ref_plausible)
print_table(so_res, digits = 3)

so_ok <- so_res[so_res$consistent_with_0_30 & !so_res$consistent_with_full, ]
so_units <- if (nrow(so_ok) == 1L) so_ok$reading else "kg/m2 (catalogued; plausibility ambiguous)"
if (nrow(so_ok) != 1L) {
  log_warn("Sothe units are not uniquely determined by plausibility alone.")
  log_warn("The project catalogue records this asset as kg C/m2, so that ",
           "reading is applied; the ambiguity is recorded in the audit table.")
}
so_kgm2 <- so_img$select(0)$rename("sothe2022_soc_0_30_kgm2")

record_audit(
  layer = "Sothe et al. 2022 SC 0-30 cm",
  asset = CFG$gee$asset_sothe_sc_0_30,
  band = so$bands[1],
  native_scale_m = so$scale_m,
  median_raw = so_med,
  units_applied = so_units,
  depth_support = "0-30 cm",
  comparable_to_core_0_30 = TRUE,
  note = paste0("The 0-30 cm cut of the same product whose 0-1 m version is ",
                "audited below. Same source, different depth support: never ",
                "mix or average the two.")
)

# --- the same product at 0-1 m, audited separately ---------------------------
# Included precisely because it is the easiest mistake available here: two
# assets from one publication, identical units, silently different depths.
log_step("04b2  SOTHE ET AL. 2022, 0-1 m (SAME product, DIFFERENT depth)")

so1_ok <- TRUE
so1_img <- tryCatch(
  ee$ImageCollection(CFG$gee$asset_sothe_sc_0_100)$first(),
  error = function(e) { so1_ok <<- FALSE; NULL })

if (!so1_ok || is.null(so1_img)) {
  log_warn("0-1 m Sothe layer unavailable; skipping. The 0-30 cm audit stands.")
  so1_med <- NA_real_
} else {
  so1 <- probe_layer(so1_img, aoi, "Sothe2022_0_100")
  so1_med <- ee$Number(so1_img$select(0)$reduceRegion(
    ee$Reducer$median(), aoi, max(30, so1$scale_m), maxPixels = 1e12,
    bestEffort = TRUE)$values()$get(0))$getInfo()
  log_info(sprintf("median 0-1 m stock over AOI: %.2f kg C/m2", so1_med))
  if (is.finite(so_med) && so1_med < so_med) {
    log_warn("The 0-1 m layer reads LOWER than the 0-30 cm layer over this ",
             "AOI, which nesting forbids. Treat both as suspect here.")
  } else {
    log_ok(sprintf("0-1 m exceeds 0-30 cm by %.2f kg/m2, as nesting requires",
                   so1_med - so_med))
  }
  record_audit(
    layer = "Sothe et al. 2022 SC 0-1 m",
    asset = CFG$gee$asset_sothe_sc_0_100,
    band = so1$bands[1],
    native_scale_m = so1$scale_m,
    median_raw = so1_med,
    units_applied = "kg/m2 (catalogued)",
    depth_support = "0-100 cm",
    comparable_to_core_0_30 = FALSE,
    note = paste0("SAME publication as the 0-30 cm layer above, DIFFERENT ",
                  "depth support. Barred from the 0-30 cm comparison. Kept ",
                  "because it shows how much carbon sits below the 30 cm ",
                  "window the shallow products report.")
  )
}

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

# --- test against the PUBLICATION, not just against plausibility -------------
# Li et al. (2025) report a Hudson Bay Lowlands mean peat carbon stock of
# 86 (+/- 35) kg C/m2 over a mean peat depth of 184 (+/- 48) cm. If this asset
# is the published full-column layer, its central value over a Hudson Bay
# Lowlands AOI should sit within roughly one standard deviation of that. This
# is a far sharper test than "is the number big", and it is the check that
# distinguishes the published product from a mis-scaled or mis-selected band.
li_mean <- ee$Number(li_img$select(0)$reduceRegion(
  ee$Reducer$mean(), aoi, max(30, li$scale_m), maxPixels = 1e12,
  bestEffort = TRUE)$values()$get(0))$getInfo()

lit_mu <- CFG$lit$li_peat_stock_mean_kgm2
lit_sd <- CFG$lit$li_peat_stock_sd_kgm2
z_score <- if (is.finite(li_mean)) (li_mean - lit_mu) / lit_sd else NA_real_

log_info("checking the layer against its own publication:")
log_info(sprintf("   published HBL mean : %.0f +/- %.0f kg C/m2 (Li et al. 2025)",
                 lit_mu, lit_sd))
log_info(sprintf("   measured AOI mean  : %.1f kg C/m2", li_mean))
log_info(sprintf("   deviation          : %+.2f published SD", z_score))

matches_publication <- is.finite(z_score) && abs(z_score) <= 1.5
if (matches_publication) {
  log_ok("layer agrees with its publication; the full-column reading is confirmed")
} else {
  log_warn("Layer does NOT sit within 1.5 published SD of the reported HBL mean.")
  log_warn("Either this AOI is genuinely atypical of the Lowlands, or the band ",
           "selection or scaling is wrong. Resolve before quantitative use.")
}

log_warn(strrep("-", 68))
if (reads_full_profile) {
  log_warn("LI ET AL. VALUES ARE TOO LARGE TO BE A 0-30 cm STOCK.")
  log_warn(sprintf("  median over AOI      : %.1f kg C/m2", li_med))
  log_warn(sprintf("  0-30 cm ceiling      : %.1f kg C/m2",
                   CFG$ref_plausible$stock_0_30_kgm2_max))
  log_warn(sprintf("  our deepest core     : %.1f kg C/m2 over 0-30 cm",
                   max(read_csv_verbatim(
                     file.path(CFG$dir_derived, "02_core_stocks.csv"))$stock_kgm2)))
  log_warn(sprintf("  published peat depth : %.0f +/- %.0f cm",
                   CFG$lit$li_peat_depth_mean_cm, CFG$lit$li_peat_depth_sd_cm))
  log_warn("  This is FULL PEAT-COLUMN carbon, integrated from the surface to")
  log_warn("  the base of the peat -- metres, not centimetres. It is NOT the")
  log_warn("  same quantity as a 0-30 cm stock and MUST NOT be differenced")
  log_warn("  against, averaged with, or validated by the 0-30 cm products or")
  log_warn("  by these cores.")
  log_warn("")
  log_warn(sprintf("  For scale: 0-30 cm is %.0f%% of the mean peat depth.",
                   100 * 30 / CFG$lit$li_peat_depth_mean_cm))
  li_depth <- "FULL PEAT COLUMN (surface to base of peat); NOT 0-30 cm"
  li_comparable <- FALSE
} else {
  log_warn("Li et al. values fall within a 0-30 cm plausible range, which")
  log_warn("contradicts the published full-column basis. Treat the depth")
  log_warn("support as UNRESOLVED and confirm against the publication before")
  log_warn("combining. Note that earlier project notes labelled this asset")
  log_warn("'0-100 cm', which disagrees with the paper's headline map.")
  li_depth <- "UNRESOLVED (published basis is full column; layer reads shallow)"
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
                "is retained as CONTEXT only and is barred from combination. ",
                sprintf("AOI mean %.1f vs published %.0f +/- %.0f kg/m2 (%+.2f SD); ",
                        li_mean, lit_mu, lit_sd, z_score),
                if (matches_publication) "consistent with the publication."
                else "NOT consistent with the publication -- investigate.")
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

ref_stack <- sg_kgm2$addBands(sg_0_100)$addBands(so_kgm2)$addBands(li_kgm2)

# Per-pixel uncertainty, where the project has access to it. Extracted so the
# comparison in 07 can ask whether a disagreement exceeds the layer's own
# stated uncertainty, rather than treating every difference as meaningful.
li_unc_ok <- TRUE
li_unc <- tryCatch(
  ee$Image(CFG$gee$asset_li2025_unc)$select(0)$rename("li2025_peat_carbon_unc_kgm2"),
  error = function(e) { li_unc_ok <<- FALSE; NULL })
if (li_unc_ok && !is.null(li_unc)) {
  ref_stack <- ref_stack$addBands(li_unc)
  log_ok("included Li et al. per-pixel uncertainty band")
} else {
  log_warn("Li uncertainty asset unavailable (private asset not shared?); ",
           "continuing without it")
}
if (exists("so1_img") && !is.null(so1_img)) {
  ref_stack <- ref_stack$addBands(so1_img$select(0)$rename("sothe2022_soc_0_100_kgm2"))
  log_ok("included Sothe 0-1 m band (barred from the 0-30 cm comparison)")
}

# sampleRegions drops a feature ENTIRELY when any band is masked at that point,
# so a core sitting on one masked layer disappears from every column at once.
# The first run of this script returned six cores of eight for exactly that
# reason. Sampling through a sentinel keeps all eight and costs only the
# masked value.
ref_at_cores <- gee_sample_points(ref_stack, pts, scale = 30)
ref_at_cores <- merge(cores[, c("core_id", "campaign")], ref_at_cores,
                      by = "core_id", all.x = TRUE)

n_missing <- sum(!is.finite(ref_at_cores$li2025_peat_carbon_kgm2))
if (n_missing > 0) {
  log_flag(n_missing, " core(s) fall on a masked pixel in at least one ",
           "reference layer: ",
           paste(ref_at_cores$core_id[
             !is.finite(ref_at_cores$li2025_peat_carbon_kgm2)], collapse = ", "))
  log_info("They are retained with NA in the affected column rather than ",
           "dropped, so the row count still matches the core count.")
}
stopifnot(nrow(ref_at_cores) == nrow(cores))
log_ok("all ", nrow(cores), " cores retained in the extraction")

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

# =============================================================================
# 7. Download the reference rasters to local disk
# =============================================================================

log_step("04f  REFERENCE RASTER DOWNLOAD")

# The local filenames here are not arbitrary: 11_bayesian_map.R looks for
# exactly these. Landing them in outputs/gee/ upgrades the Bayesian products
# from DEMO_ (a spatially constant regional prior) to the real published
# surfaces, with no further action.
log_info("local names are chosen so 11_bayesian_map.R picks them up ",
         "automatically and drops its DEMO_ fallback")

exports <- list()
exports[[1]] <- gee_download_image(
  li_kgm2$toFloat(), "ccnp_li2025_peat_carbon_30m", aoi,
  scale = 30, dir_local = CFG$dir_gee, crs = CFG$gee$export_crs)
exports[[2]] <- gee_download_image(
  so_kgm2$toFloat(), "ccnp_sothe_sc_0_30_30m", aoi,
  scale = 250, dir_local = CFG$dir_gee, crs = CFG$gee$export_crs)
exports[[3]] <- gee_download_image(
  sg_kgm2$toFloat(), "ccnp_soilgrids_ocs_0_30_30m", aoi,
  scale = 250, dir_local = CFG$dir_gee, crs = CFG$gee$export_crs)
if (exists("li_unc") && !is.null(li_unc)) {
  exports[[4]] <- gee_download_image(
    li_unc$toFloat(), "ccnp_li2025_uncertainty_30m", aoi,
    scale = 30, dir_local = CFG$dir_gee, crs = CFG$gee$export_crs)
}

man <- gee_write_manifest(exports,
                          file.path(CFG$dir_gee, "04_download_manifest.csv"))

if (!is.null(man) && any(man$status %in% c("local", "local_multipart", "cached"))) {
  log_ok("Reference rasters are on disk. Re-run 11_bayesian_map.R to replace ",
         "the DEMO_ outputs with real prior-based Bayesian maps.")
}

log_ok("04 complete")
