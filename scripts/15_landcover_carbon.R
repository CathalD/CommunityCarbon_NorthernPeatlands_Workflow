# =============================================================================
# 15_landcover_carbon.R
#
# Carbon stock per ecosystem type: zonal statistics of a carbon surface within
# land-cover classes, over the coast-oriented AOI from 13.
#
# WHICH LAND COVER, AND WHY IT MATTERS HERE
#
# ESA WorldCover cannot separate treed bog growing on metres of peat from
# upland black spruce on mineral soil -- both read as tree cover. In the Hudson
# Bay Lowlands that is precisely the distinction that decides carbon, so
# WorldCover alone would blur the two ecosystems whose difference this whole
# project is about.
#
# GWL_FCS30 (Zhang et al. 2023) carries a wetland classification and is used as
# the primary ecosystem layer, with WorldCover reported alongside so the two
# can be compared. Where they disagree is itself informative.
#
# WHAT IS BEING AVERAGED
#
# The carbon surface is a PUBLISHED product, not the community cores. Eight
# cores cannot populate per-class means for a dozen classes; 05 already shows
# they cannot support two. So these are the published layers' own estimates,
# summarised by ecosystem, with the cores overlaid as ground truth wherever
# they happen to fall. Class means carry their sample count and spread, and a
# class containing no core is marked as having no local verification at all.
#
# REQUIRES rgee and an authenticated Earth Engine session.
#
# OUTPUT  outputs/tables/15_carbon_by_landcover.csv
#         outputs/tables/15_cores_by_landcover.csv
#         GEE exports: ecosystem carbon raster + class map over the AOI
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

log_step("15  CARBON BY ECOSYSTEM TYPE")

if (!requireNamespace("rgee", quietly = TRUE)) {
  fail_loudly(
    "rgee is not installed",
    "This script needs Earth Engine for the land-cover and carbon rasters.",
    remedy = 'install.packages("rgee"); rgee::ee_install(); then re-run.')
}
suppressPackageStartupMessages(library(rgee))
ee_Initialize(project = CFG$gee$project, drive = TRUE)
log_ok("Earth Engine initialised")

cores <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                          "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")
if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the ecosystem step",
              remedy = "Re-run 01_ingest_qc.R.")
}

# =============================================================================
# 1. AOI -- the oriented polygon from 13, so every product shares one boundary
# =============================================================================

log_step("15a  AOI")

dir_sp <- file.path(CFG$root, "outputs", "spatial")
aoi_src <- tryCatch(read_aoi_ring(dir_sp), error = function(e) NULL)

if (!is.null(aoi_src)) {
  ring <- aoi_src$ring
  aoi  <- ee$Geometry$Polygon(list(ring_to_ee_coords(ring)))
  log_ok("using the coast-oriented AOI from 13 (", nrow(ring),
         " vertices, read from ", aoi_src$source, ")")
  log_info(sprintf("ring spans %.4f..%.4f E, %.4f..%.4f N",
                   min(ring$lon), max(ring$lon), min(ring$lat), max(ring$lat)))
  # Confirm the geometry Earth Engine actually built matches what was read.
  a_km2 <- aoi$area(maxError = 100)$getInfo() / 1e6
  log_info(sprintf("Earth Engine reports the AOI area as %.0f km2", a_km2))
  if (!is.finite(a_km2) || a_km2 < 100 || a_km2 > 50000) {
    fail_loudly("The AOI Earth Engine built is not the AOI that was written",
                sprintf("Area came back as %.1f km2, which is implausible.", a_km2),
                remedy = "Re-run 13_aoi_boundary.R and check aoi_vertices.csv.")
  }
  log_ok("AOI round-trip verified against Earth Engine")
} else {
  bb <- CFG$aoi_bbox
  aoi <- ee$Geometry$Rectangle(c(bb[["xmin"]], bb[["ymin"]],
                                 bb[["xmax"]], bb[["ymax"]]))
  log_warn("No usable AOI file found; falling back to the axis-aligned box. ",
           "Run 13_aoi_boundary.R for the coast-oriented AOI.")
}

pts <- ee$FeatureCollection(lapply(seq_len(nrow(cores)), function(i) {
  ee$Feature(ee$Geometry$Point(c(cores$longitude[i], cores$latitude[i])),
             list(core_id = cores$core_id[i], campaign = cores$campaign[i]))
}))

# =============================================================================
# 2. Ecosystem layers
# =============================================================================

log_step("15b  ECOSYSTEM LAYERS")

gwl <- ee$ImageCollection(CFG$gee$asset_gwl_fcs30)$
  filterBounds(aoi)$mosaic()$rename("gwl_class")$clip(aoi)
wc <- ee$ImageCollection(CFG$gee$asset_worldcover)$first()$
  select("Map")$rename("worldcover")$clip(aoi)

gwl_labels <- data.frame(
  class = c(180L, 181L, 182L, 183L, 184L, 185L, 186L, 187L, 188L),
  label = c("Non-wetland", "Permanent water", "Swamp", "Marsh",
            "Flooded flat", "Saline wetland", "Mangrove", "Salt marsh",
            "Tidal flat"),
  stringsAsFactors = FALSE)
wc_labels <- data.frame(
  class = c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 90L, 95L, 100L),
  label = c("Tree cover", "Shrubland", "Grassland", "Cropland", "Built-up",
            "Bare/sparse", "Snow and ice", "Water", "Herbaceous wetland",
            "Mangroves", "Moss and lichen"),
  stringsAsFactors = FALSE)
log_ok("GWL_FCS30 (primary, wetland-aware) and ESA WorldCover (comparison)")

# =============================================================================
# 3. Carbon surfaces
# =============================================================================

log_step("15c  CARBON SURFACES")

# SoilGrids, rebuilt the same way as in 04 so the integration matches the
# cores'. Kept at 0-30 cm: the depth the cores actually measure.
sg_soc  <- ee$Image(CFG$gee$asset_soilgrids_soc)
sg_bdod <- ee$Image(CFG$gee$asset_soilgrids_bdod)
sg_layer <- function(sb, bb2, th) {
  sg_soc$select(sb)$divide(10)$
    multiply(sg_bdod$select(bb2)$divide(100))$multiply(th)$divide(100)
}
sg_0_30 <- sg_layer("soc_0-5cm_mean", "bdod_0-5cm_mean", 5)$
  add(sg_layer("soc_5-15cm_mean", "bdod_5-15cm_mean", 10))$
  add(sg_layer("soc_15-30cm_mean", "bdod_15-30cm_mean", 15))$
  rename("soilgrids_0_30_kgm2")$clip(aoi)

li <- ee$Image(CFG$gee$asset_li2025)$select(0)$
  rename("li_full_column_kgm2")$clip(aoi)

log_ok("SoilGrids 0-30 cm (constructed) and Li et al. full peat column")
log_warn("These two have DIFFERENT depth support and are summarised in ",
         "separate columns. They are never averaged, differenced or ranked ",
         "against each other -- see the guard in 04.")

# =============================================================================
# 4. Zonal statistics
# =============================================================================

log_step("15d  ZONAL STATISTICS BY ECOSYSTEM")

#' Mean, sd, count and area of a carbon band within each class of a
#' categorical band, returned as a plain data frame.
zonal_by_class <- function(carbon, classes, class_band, scale_m, labels, tag) {
  img <- carbon$addBands(classes)
  red <- ee$Reducer$mean()$
    combine(ee$Reducer$stdDev(), "", TRUE)$
    combine(ee$Reducer$count(), "", TRUE)$
    group(groupField = 1L, groupName = "class")
  res <- img$reduceRegion(reducer = red, geometry = aoi, scale = scale_m,
                          maxPixels = 1e13, bestEffort = TRUE)$getInfo()
  g <- res$groups
  if (is.null(g) || !length(g)) return(NULL)
  out <- do.call(rbind, lapply(g, function(x) data.frame(
    class = as.integer(x$class),
    mean_kgm2 = if (is.null(x$mean)) NA_real_ else x$mean,
    sd_kgm2   = if (is.null(x$stdDev)) NA_real_ else x$stdDev,
    n_pixels  = if (is.null(x$count)) NA_real_ else x$count,
    stringsAsFactors = FALSE)))
  out$label <- labels$label[match(out$class, labels$class)]
  out$label[is.na(out$label)] <- paste0("class ", out$class[is.na(out$label)])
  out$area_km2 <- out$n_pixels * (scale_m / 1000)^2
  out$pct_of_aoi <- 100 * out$area_km2 / sum(out$area_km2, na.rm = TRUE)
  out$carbon_layer <- tag
  out[order(-out$area_km2), ]
}

tabs <- list()
for (spec in list(
  list(c = sg_0_30, tag = "SoilGrids 0-30 cm", depth = "0-30 cm", s = 250),
  list(c = li,      tag = "Li et al. full peat column", depth = "full column", s = 30)
)) {
  for (lc in list(list(img = gwl, band = "gwl_class", lab = gwl_labels, nm = "GWL_FCS30"),
                  list(img = wc,  band = "worldcover", lab = wc_labels,  nm = "WorldCover"))) {
    z <- tryCatch(
      zonal_by_class(spec$c, lc$img, lc$band, spec$s, lc$lab, spec$tag),
      error = function(e) { log_warn("zonal failed for ", spec$tag, " x ",
                                     lc$nm, ": ", conditionMessage(e)); NULL })
    if (!is.null(z)) {
      z$landcover_source <- lc$nm
      z$depth_support <- spec$depth
      tabs[[length(tabs) + 1L]] <- z
      log_info(spec$tag, " x ", lc$nm, ": ", nrow(z), " classes")
    }
  }
}

if (!length(tabs)) {
  fail_loudly("No zonal statistics were produced",
              remedy = "Check Earth Engine asset access and the AOI.")
}
zt <- do.call(rbind, tabs)
rownames(zt) <- NULL

log_info("carbon by ecosystem (GWL_FCS30, 0-30 cm):")
print_table(zt[zt$landcover_source == "GWL_FCS30" &
               zt$depth_support == "0-30 cm",
               c("label", "area_km2", "pct_of_aoi", "mean_kgm2", "sd_kgm2")],
            digits = 2)

# =============================================================================
# 5. Which ecosystems do the cores actually verify?
# =============================================================================

log_step("15e  GROUND VERIFICATION BY ECOSYSTEM")

samp <- gwl$addBands(wc)$sampleRegions(collection = pts, scale = 30,
                                       geometries = FALSE, tileScale = 4)
cl <- sf::st_drop_geometry(ee_as_sf(samp, maxFeatures = 1000))
cl$gwl_label <- gwl_labels$label[match(cl$gwl_class, gwl_labels$class)]
cl$wc_label  <- wc_labels$label[match(cl$worldcover, wc_labels$class)]

r30 <- stocks[stocks$scenario == "as_supplied" &
              stocks$window == "reference_0_30", ]
cl$stock_0_30_kgm2 <- r30$stock_kgm2[match(cl$core_id, r30$core_id)]
cl$is_lower_bound  <- r30$is_lower_bound[match(cl$core_id, r30$core_id)]

log_info("which ecosystem each core actually sits in:")
print_table(cl[, c("core_id", "campaign", "gwl_label", "wc_label",
                   "stock_0_30_kgm2", "is_lower_bound")], digits = 2)

# Disagreement between the two land-cover products at the cores is a direct,
# checkable illustration of why the wetland-aware layer was chosen.
if (nrow(cl)) {
  mism <- cl$gwl_label != cl$wc_label
  log_info(sum(mism, na.rm = TRUE), " of ", nrow(cl),
           " cores are labelled differently by the two products")
  peat_as_tree <- grepl("peat", cl$campaign) & cl$wc_label == "Tree cover"
  if (any(peat_as_tree, na.rm = TRUE)) {
    log_warn("WorldCover calls ", sum(peat_as_tree, na.rm = TRUE),
             " peat core site(s) 'Tree cover', which is exactly the ",
             "confusion GWL_FCS30 was brought in to resolve.")
  }
}

verified <- unique(cl$gwl_label[!is.na(cl$gwl_label)])
gwl_area <- zt[zt$landcover_source == "GWL_FCS30" & zt$depth_support == "0-30 cm", ]
gwl_area$has_core <- gwl_area$label %in% verified
unver <- gwl_area[!gwl_area$has_core, ]
log_warn(sprintf("%.1f%% of the AOI lies in ecosystem classes containing NO core.",
                 sum(unver$pct_of_aoi, na.rm = TRUE)))
log_warn("Those class means come entirely from the published layers and have ",
         "no local verification whatsoever. They are marked accordingly.")

# =============================================================================
# 6. Exports and outputs
# =============================================================================

log_step("15f  OUTPUTS")

write_csv_logged(zt, file.path(CFG$dir_tables, "15_carbon_by_landcover.csv"),
                 "carbon per ecosystem class, by layer and land-cover source")
write_csv_logged(cl, file.path(CFG$dir_tables, "15_cores_by_landcover.csv"),
                 "which ecosystem each community core sits in")
write_csv_logged(gwl_area[, c("label", "area_km2", "pct_of_aoi", "mean_kgm2",
                              "has_core")],
                 file.path(CFG$dir_tables, "15_ecosystem_verification.csv"),
                 "which ecosystems have any ground truth at all")

# Each product goes to Drive AND to local disk. The local copy is what the
# reports embed and what `terra` can open; the Drive copy is what gets shared.
exports <- list(
  gee_export_image(gwl$addBands(wc)$toInt16(), "ccnp_ecosystem_classes_30m",
                   aoi, scale = CFG$gee$export_scale_m,
                   crs = CFG$gee$export_crs, dir_local = CFG$dir_gee),
  gee_export_image(sg_0_30$toFloat(), "ccnp_soilgrids_ocs_0_30_30m",
                   aoi, scale = 250, crs = CFG$gee$export_crs,
                   dir_local = CFG$dir_gee, overwrite = FALSE),
  gee_export_image(li$toFloat(), "ccnp_li2025_peat_carbon_30m",
                   aoi, scale = CFG$gee$export_scale_m,
                   crs = CFG$gee$export_crs, dir_local = CFG$dir_gee,
                   overwrite = FALSE)
)

# A carbon-per-ecosystem raster: every pixel painted with its own class's mean.
# This is the layer people actually want to look at, and it is a CLASS-MEAN
# ASSIGNMENT in exactly the sense 05 defines -- not a per-pixel estimate.
gwl_tbl <- zt[zt$landcover_source == "GWL_FCS30" & zt$depth_support == "0-30 cm", ]
if (nrow(gwl_tbl)) {
  eco_mean <- gwl$remap(as.integer(gwl_tbl$class),
                        as.numeric(gwl_tbl$mean_kgm2))$
    rename("ecosystem_mean_soc_0_30_kgm2")$toFloat()
  exports[[length(exports) + 1L]] <- gee_export_image(
    eco_mean, "ccnp_ecosystem_mean_soc_0_30_kgm2", aoi,
    scale = CFG$gee$export_scale_m, crs = CFG$gee$export_crs,
    dir_local = CFG$dir_gee)
  log_warn("ccnp_ecosystem_mean_soc_0_30_kgm2 is a CLASS-MEAN ASSIGNMENT: ",
           "every pixel carries its ecosystem's average, not its own value. ",
           "Within-class variation is discarded by construction.")
}

gee_write_manifest(exports, file.path(CFG$dir_gee, "15_export_manifest.csv"))
log_info("local GeoTIFFs are in ", CFG$dir_gee,
         "; re-run 11 to build the Bayesian products on the real priors")
log_ok("15 complete")
