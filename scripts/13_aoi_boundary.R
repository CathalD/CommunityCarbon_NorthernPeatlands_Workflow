# =============================================================================
# 13_aoi_boundary.R
#
# Build a coast-following area of interest: a rectangle aligned to the trend of
# the Hudson Bay shoreline at Fort Severn, with a large buffer, rather than a
# north-up bounding box.
#
# WHY ORIENTED
#
# The shoreline here runs north-west to south-east and the community sampling
# follows it. A north-up box either clips the ends of the sampled transect or,
# if enlarged to contain them, spends most of its area on open water to the
# north-east and on ground far inland that was never visited. An oriented
# rectangle contains the same sampling in less area, and its axes mean
# something: distance ALONG the coast, and distance ACROSS it, which in an
# isostatically rebounding landscape is a proxy for time since the sea left.
#
# The orientation is DERIVED from the cores' principal axis rather than drawn
# by eye, and the script refuses to orient at all if the cores are not
# elongated enough for a direction to be meaningful.
#
# OUTPUT  outputs/spatial/aoi_coast_oriented.geojson
#         outputs/spatial/aoi_coast_oriented_mask.tif
#         outputs/spatial/aoi_for_earthengine.js
#         outputs/tables/13_aoi_definition.csv
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

log_step("13  COAST-FOLLOWING AREA OF INTEREST")

dir_sp <- file.path(CFG$root, "outputs", "spatial")
if (!dir.exists(dir_sp)) dir.create(dir_sp, recursive = TRUE)

cores <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                          "scripts/01_ingest_qc.R")
if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the AOI step",
              "The AOI would be drawn in the wrong hemisphere.",
              remedy = "Re-run 01_ingest_qc.R.")
}

# =============================================================================
# 1. Orientation
# =============================================================================

log_step("13a  ORIENTATION")

pa <- principal_axis_bearing(cores$longitude, cores$latitude)
log_info(sprintf("cores' principal axis: %.1f degrees from north", pa$bearing_deg))
log_info(sprintf("variance explained by that axis: %.1f%%", 100 * pa$var_explained))
log_info(sprintf("core spread: %.1f km along, %.1f km across",
                 pa$extent_along_km, pa$extent_across_km))

cfg_bearing <- CFG$aoi_oriented$bearing_deg
use_oriented <- TRUE

if (!is.null(cfg_bearing)) {
  bearing <- cfg_bearing
  log_info(sprintf("using the bearing set in config: %.1f degrees", bearing))
} else if (pa$var_explained < CFG$aoi_oriented$min_var_explained) {
  use_oriented <- FALSE
  bearing <- 0
  log_warn(sprintf(paste0("The cores explain only %.0f%% of their variance on ",
                          "one axis, below the %.0f%% threshold."),
                   100 * pa$var_explained,
                   100 * CFG$aoi_oriented$min_var_explained))
  log_warn("They are closer to a cluster than a line, so there is no ",
           "defensible direction to orient to. Falling back to an ",
           "axis-aligned box.")
} else {
  bearing <- pa$bearing_deg
  log_ok(sprintf("orientation derived from the sampling: %.1f degrees", bearing))
  log_info("For a transect that follows a shoreline this recovers the ")
  log_info("shoreline's local trend, which is why it is not drawn by eye.")
}

# Sanity: the Hudson Bay coast at Fort Severn trends roughly NW-SE, so a
# derived bearing far from that deserves a second look before it is trusted.
if (use_oriented) {
  expected <- 135
  off <- min(abs(bearing - expected), 180 - abs(bearing - expected))
  if (off > 45) {
    log_warn(sprintf(paste0("Derived bearing (%.0f deg) is %.0f deg from the ",
                            "NW-SE trend expected for this coast. Check the ",
                            "core layout before relying on the orientation."),
                     bearing, off))
  } else {
    log_ok(sprintf("consistent with the NW-SE coastal trend (%.0f deg off)", off))
  }
}

# =============================================================================
# 2. Rectangle
# =============================================================================

log_step("13b  RECTANGLE")

ab <- CFG$aoi_oriented$along_buffer_km
cb <- CFG$aoi_oriented$across_buffer_km
log_info(sprintf("buffer: %.0f km beyond the cores at each end, %.0f km each side",
                 ab, cb))

rect <- oriented_rectangle(cores$longitude, cores$latitude,
                           if (use_oriented) bearing else 0, ab, cb)

log_ok(sprintf("AOI: %.1f km long x %.1f km wide = %.0f km2",
               rect$length_km, rect$width_km, rect$area_km2))
log_info(sprintf("axis-aligned envelope: %.3f..%.3f E, %.3f..%.3f N",
                 rect$bbox[["xmin"]], rect$bbox[["xmax"]],
                 rect$bbox[["ymin"]], rect$bbox[["ymax"]]))

# The oriented rectangle should be meaningfully smaller than its own envelope,
# otherwise the rotation bought nothing.
env_km2 <- haversine_km(rect$bbox[["xmin"]], mean(cores$latitude),
                        rect$bbox[["xmax"]], mean(cores$latitude)) *
           haversine_km(0, rect$bbox[["ymin"]], 0, rect$bbox[["ymax"]])
log_info(sprintf("envelope area: %.0f km2 -- the oriented AOI is %.0f%% of it",
                 env_km2, 100 * rect$area_km2 / env_km2))

stopifnot(all(point_in_oriented_rect(cores$longitude, cores$latitude, rect)))
log_ok("all ", nrow(cores), " cores fall inside the AOI, as required")

# Where does each core sit in coast-relative coordinates?
oc <- oriented_coordinates(cores$longitude, cores$latitude, rect)
oc <- cbind(data.frame(core_id = cores$core_id, campaign = cores$campaign,
                       stringsAsFactors = FALSE), round(oc, 2))
log_info("cores in coast-relative coordinates (km along / across the trend):")
print_table(oc[order(oc$along_km), ], digits = 2)
write_csv_logged(oc, file.path(CFG$dir_tables, "13_cores_coast_relative.csv"),
                 "core positions along and across the coastal trend")

# =============================================================================
# 3. Write vector
# =============================================================================

log_step("13c  OUTPUTS")

write_geojson_polygon(rect$lon, rect$lat,
                      file.path(dir_sp, "aoi_coast_oriented.geojson"),
                      props = list(
                        name = "Fort Severn coast-oriented AOI",
                        bearing_deg = round(rect$bearing_deg, 2),
                        length_km = round(rect$length_km, 2),
                        width_km = round(rect$width_km, 2),
                        area_km2 = round(rect$area_km2, 1),
                        along_buffer_km = ab,
                        across_buffer_km = cb,
                        orientation_source = if (use_oriented)
                          "derived from core principal axis" else
                          "axis-aligned fallback",
                        crs = "EPSG:4326"))
log_ok("wrote aoi_coast_oriented.geojson")

# Earth Engine literal, so the AOI used here and in GEE are the same shape.
snippet <- c(
  "// Fort Severn coast-oriented AOI",
  sprintf("// bearing %.2f deg from north; %.1f x %.1f km; %.0f km2",
          rect$bearing_deg, rect$length_km, rect$width_km, rect$area_km2),
  sprintf("// generated by 13_aoi_boundary.R on %s", format(Sys.Date())),
  ee_polygon_snippet(rect, "aoi"),
  "Map.addLayer(aoi, {color: 'red'}, 'Fort Severn AOI');",
  "Map.centerObject(aoi, 9);"
)
writeLines(snippet, file.path(dir_sp, "aoi_for_earthengine.js"))
log_ok("wrote aoi_for_earthengine.js (paste straight into the Code Editor)")

# =============================================================================
# 4. Raster mask
# =============================================================================

log_step("13d  RASTER MASK")

# Grid the ENVELOPE, then burn the oriented rectangle into it, so the mask can
# be stacked with the other EPSG:4326 products.
gm <- make_grid(rect$bbox, CFG$spatial$cell_m)
LON <- matrix(rep(gm$x, each = gm$ny), nrow = gm$ny)
LAT <- matrix(rep(gm$y, times = gm$nx), nrow = gm$ny)
inside <- point_in_oriented_rect(as.vector(LON), as.vector(LAT), rect)
mask <- matrix(as.integer(inside), nrow = gm$ny, ncol = gm$nx)

log_info(sprintf("mask grid: %d x %d at ~%d m", gm$nx, gm$ny, CFG$spatial$cell_m))
log_info(sprintf("inside the oriented AOI: %.1f%% of the envelope",
                 100 * mean(inside)))

write_geotiff(mask, file.path(dir_sp, "aoi_coast_oriented_mask.tif"),
              xmin = gm$xmin, ymax = gm$ymax, xres = gm$xres, yres = gm$yres,
              epsg = CFG$spatial$epsg_out, datatype = "Int32", nodata = -9999L)
log_ok("wrote aoi_coast_oriented_mask.tif (1 = inside, 0 = outside)")

defn <- data.frame(
  property = c("bearing_deg", "orientation_source", "length_km", "width_km",
               "area_km2", "along_buffer_km", "across_buffer_km",
               "envelope_area_km2", "pct_of_envelope", "n_cores_inside",
               "var_explained_axis1", "crs"),
  value = c(round(rect$bearing_deg, 2),
            if (use_oriented) "derived from core principal axis" else "axis-aligned fallback",
            round(rect$length_km, 2), round(rect$width_km, 2),
            round(rect$area_km2, 1), ab, cb, round(env_km2, 1),
            round(100 * rect$area_km2 / env_km2, 1), nrow(cores),
            round(pa$var_explained, 4), "EPSG:4326"),
  stringsAsFactors = FALSE)
write_csv_logged(defn, file.path(CFG$dir_tables, "13_aoi_definition.csv"),
                 "how the AOI was defined")

# =============================================================================
# 5. Figure
# =============================================================================

if (capabilities("png")) {
  f <- file.path(CFG$dir_figures, "13_aoi_coast_oriented.png")
  grDevices::png(f, width = 1500, height = 1200, res = 170)
  graphics::par(mar = c(5.6, 4.4, 3.2, 1.2), mgp = c(2.6, 0.7, 0))
  graphics::plot(rect$lon, rect$lat, type = "n",
                 asp = 1 / cos(mean(cores$latitude) * pi / 180),
                 xlab = "Longitude", ylab = "Latitude",
                 main = "Coast-following area of interest")
  graphics::rect(CFG$aoi_bbox[["xmin"]], CFG$aoi_bbox[["ymin"]],
                 CFG$aoi_bbox[["xmax"]], CFG$aoi_bbox[["ymax"]],
                 border = "grey65", lty = 3, lwd = 1.4)
  graphics::polygon(rect$lon, rect$lat, border = "#d7301f", lwd = 2.4,
                    col = grDevices::adjustcolor("#fdd49e", 0.30))
  graphics::points(cores$longitude, cores$latitude, pch = 21,
                   bg = ifelse(grepl("peat", cores$campaign), "#1b9e77", "#7570b3"),
                   col = "white", cex = 1.6, lwd = 1.6)
  graphics::points(CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]],
                   pch = 8, cex = 1.8, lwd = 2)
  graphics::legend("topright", bty = "n", cex = 0.72,
                   legend = c("coast-oriented AOI", "previous axis-aligned box",
                              "peat cores", "mineral cores", "Fort Severn"),
                   pch = c(22, 22, 21, 21, 8),
                   pt.bg = c(grDevices::adjustcolor("#fdd49e", 0.6), NA,
                             "#1b9e77", "#7570b3", NA),
                   col = c("#d7301f", "grey65", "white", "white", "black"),
                   pt.cex = c(1.7, 1.7, 1.2, 1.2, 1.2))
  graphics::mtext(sprintf(paste0("Bearing %.0f deg from north, derived from the ",
                                 "cores' principal axis. %.0f x %.0f km, %.0f km2."),
                          rect$bearing_deg, rect$length_km, rect$width_km,
                          rect$area_km2),
                  side = 1, line = 4.2, cex = 0.68, adj = 0, col = "grey30")
  grDevices::dev.off()
  log_ok("wrote ", basename(f))
}

log_ok("13 complete")
