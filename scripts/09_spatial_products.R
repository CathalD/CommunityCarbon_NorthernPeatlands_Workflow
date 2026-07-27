# =============================================================================
# 09_spatial_products.R
#
# Spatial outputs that can be produced from the cores alone, with no Earth
# Engine session and no GDAL. Written with the base-R GeoTIFF writer in
# R/geotiff.R, so this script runs anywhere R runs.
#
# WHAT IS PRODUCED, AND WHAT EACH ONE IS
#
#   cores.geojson                      OBSERVATION. The eight cores with their
#                                      stocks, flags and depth support.
#   aoi.geojson                        The study area rectangle.
#
#   dist_to_nearest_core_m.tif         GEOMETRY. Distance from every pixel to
#                                      the nearest core. This is the single
#                                      most honest spatial product here: it
#                                      shows how far any map would be reaching.
#
#   nearest_core_index.tif             DIAGNOSTIC. Thiessen polygons as a
#                                      raster: which core is closest. Shows
#                                      which core any interpolator would lean
#                                      on at each pixel.
#
#   DIAGNOSTIC_soc_nearest_core_kgm2.tif
#                                      DIAGNOSTIC, NOT A CARBON MAP. The
#                                      nearest core's stock painted across its
#                                      Thiessen cell. Included because it makes
#                                      the sampling gap visible, not because
#                                      the value at an unsampled pixel is
#                                      believable.
#
#   PRODUCT1b_soc_within_credible_radius_kgm2.tif
#                                      SAMPLE-BASED, MASKED. The same surface,
#                                      but NoData wherever no core is within
#                                      the credible radius. This is the only
#                                      version of a core-derived surface this
#                                      workflow is willing to call a product,
#                                      and the fraction of the study area it
#                                      covers is reported.
#
# Note what is deliberately absent: no kriged surface, no IDW with a fitted
# structure. 05 shows why. Nearest-neighbour assignment is used here precisely
# because it assumes NO spatial covariance structure -- it claims only "this is
# the closest measurement", which is a statement the data can support.
#
# INPUT   data/derived/01_cores_qc.csv, 02_core_stocks.csv
# OUTPUT  outputs/spatial/*.tif, *.geojson
#         outputs/tables/09_spatial_coverage.csv
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

log_step("09  SPATIAL PRODUCTS (no Earth Engine required)")

dir_sp <- file.path(CFG$root, "outputs", "spatial")
if (!dir.exists(dir_sp)) dir.create(dir_sp, recursive = TRUE)

cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

# Re-assert the geolocation invariant. Everything below writes georeferenced
# files; a sign error here would produce confidently wrong maps.
if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the spatial-output step",
              "Every raster and GeoJSON written here would be in the wrong hemisphere.",
              remedy = "Re-run 01_ingest_qc.R and inspect the geolocation gate.")
}
log_ok("geolocation invariant re-checked before writing georeferenced output")

# =============================================================================
# 1. Vector outputs
# =============================================================================

log_step("09a  VECTOR")

cs <- stocks[stocks$scenario == "as_supplied" & stocks$window == "common_support", ]
r30 <- stocks[stocks$scenario == "as_supplied" & stocks$window == "reference_0_30", ]

pts <- cores[, c("core_id", "campaign", "year", "latitude", "longitude",
                 "n_segments", "core_length_cm", "bd_min_gcm3", "bd_max_gcm3",
                 "om_max_pct", "n_flagged_segs", "dist_from_fort_severn_km")]
pts$stock_common_kgm2   <- cs$stock_kgm2[match(pts$core_id, cs$core_id)]
pts$common_window_cm    <- cs$window_hi_cm[match(pts$core_id, cs$core_id)]
pts$stock_0_30_kgm2     <- r30$stock_kgm2[match(pts$core_id, r30$core_id)]
pts$stock_0_30_is_lower_bound <- r30$is_lower_bound[match(pts$core_id, r30$core_id)]
pts$depth_covered_0_30_cm <- r30$covered_cm[match(pts$core_id, r30$core_id)]

write_geojson_points(pts, file.path(dir_sp, "cores.geojson"))
log_ok("wrote cores.geojson (", nrow(pts), " points, ",
       ncol(pts) - 2, " attributes)")

write_geojson_bbox(CFG$aoi_bbox, file.path(dir_sp, "aoi.geojson"),
                   props = list(name = "Fort Severn SOC study area",
                                crs = "EPSG:4326"))
log_ok("wrote aoi.geojson")

# =============================================================================
# 2. Grid
# =============================================================================

log_step("09b  GRID")

g <- make_grid(CFG$aoi_bbox, CFG$spatial$cell_m)
log_info(sprintf("grid: %d x %d cells at ~%d m (%.5f x %.5f deg)",
                 g$nx, g$ny, CFG$spatial$cell_m, g$xres, g$yres))
log_info(sprintf("extent: %.4f..%.4f E, %.4f..%.4f N",
                 CFG$aoi_bbox[["xmin"]], CFG$aoi_bbox[["xmax"]],
                 CFG$aoi_bbox[["ymin"]], CFG$aoi_bbox[["ymax"]]))
log_info(sprintf("approximate ground extent: %.1f x %.1f km",
                 haversine_km(CFG$aoi_bbox[["xmin"]], g$lat_mid,
                              CFG$aoi_bbox[["xmax"]], g$lat_mid),
                 haversine_km(0, CFG$aoi_bbox[["ymin"]], 0, CFG$aoi_bbox[["ymax"]])))

# Cell-centre coordinates as full matrices, rows north to south.
LON <- matrix(rep(g$x, each = g$ny), nrow = g$ny, ncol = g$nx)
LAT <- matrix(rep(g$y, times = g$nx), nrow = g$ny, ncol = g$nx)

# =============================================================================
# 3. Distance to nearest core, and which core that is
# =============================================================================

log_step("09c  DISTANCE AND NEAREST-CORE FIELDS")

n_core <- nrow(cores)
best_d <- matrix(Inf, g$ny, g$nx)
best_i <- matrix(NA_integer_, g$ny, g$nx)

for (i in seq_len(n_core)) {
  d <- haversine_km(LON, LAT, cores$longitude[i], cores$latitude[i])
  hit <- d < best_d
  best_d[hit] <- d[hit]
  best_i[hit] <- i
}
log_ok(sprintf("computed %s cell-to-core distances",
               format(g$ny * g$nx * n_core, big.mark = ",")))
log_info(sprintf("distance to nearest core: %.2f - %.2f km (median %.2f)",
                 min(best_d), max(best_d), median(best_d)))

write_geotiff(best_d * 1000, file.path(dir_sp, "dist_to_nearest_core_m.tif"),
              xmin = g$xmin, ymax = g$ymax, xres = g$xres, yres = g$yres,
              epsg = CFG$spatial$epsg_out, datatype = "Float32")
log_ok("wrote dist_to_nearest_core_m.tif")

write_geotiff(best_i, file.path(dir_sp, "nearest_core_index.tif"),
              xmin = g$xmin, ymax = g$ymax, xres = g$xres, yres = g$yres,
              epsg = CFG$spatial$epsg_out, datatype = "Int32")
legend <- data.frame(index = seq_len(n_core), core_id = cores$core_id,
                     campaign = cores$campaign, stringsAsFactors = FALSE)
write_csv_logged(legend, file.path(dir_sp, "nearest_core_index_legend.csv"),
                 "value -> core id for nearest_core_index.tif")
log_ok("wrote nearest_core_index.tif (Thiessen assignment)")

# =============================================================================
# 4. Nearest-core carbon surface, and the credible-radius mask
# =============================================================================

log_step("09d  NEAREST-CORE CARBON SURFACE")

stock_by_index <- cs$stock_kgm2[match(cores$core_id, cs$core_id)]
soc <- matrix(stock_by_index[best_i], g$ny, g$nx)

write_geotiff(soc, file.path(dir_sp, "DIAGNOSTIC_soc_nearest_core_kgm2.tif"),
              xmin = g$xmin, ymax = g$ymax, xres = g$xres, yres = g$yres,
              epsg = CFG$spatial$epsg_out, datatype = "Float32")
log_ok("wrote DIAGNOSTIC_soc_nearest_core_kgm2.tif")
log_warn("That file is a DIAGNOSTIC, not a carbon map. It assigns each pixel ",
         "the stock of whichever core happens to be closest, with no claim ",
         "that the value holds there.")

# The credible radius: how far a single core can reasonably speak for. Taken
# from the data as the median nearest-neighbour spacing between cores, which is
# the scale at which this sampling design actually resolves anything.
nn <- nn_distance_summary(cores$longitude, cores$latitude, cores$core_id)
radius_km <- median(nn$nn_dist_km)
log_info(sprintf("credible radius = median core-to-core spacing = %.2f km",
                 radius_km))
log_info("Beyond that distance no core has a neighbour to corroborate it, so ")
log_info("there is no evidence the value carries further.")

soc_masked <- soc
soc_masked[best_d > radius_km] <- NA_real_
frac <- mean(!is.na(soc_masked))

write_geotiff(soc_masked,
              file.path(dir_sp, "PRODUCT1b_soc_within_credible_radius_kgm2.tif"),
              xmin = g$xmin, ymax = g$ymax, xres = g$xres, yres = g$yres,
              epsg = CFG$spatial$epsg_out, datatype = "Float32")
log_ok("wrote PRODUCT1b_soc_within_credible_radius_kgm2.tif")

log_warn(strrep("-", 68))
log_warn(sprintf("THE EIGHT CORES SPEAK FOR %.1f%% OF THE STUDY AREA.", 100 * frac))
log_warn(sprintf("  Within %.1f km of a core: %.1f%% of the AOI.",
                 radius_km, 100 * frac))
log_warn(sprintf("  The remaining %.1f%% is NoData in the masked product,",
                 100 * (1 - frac)))
log_warn("  because nothing in this dataset constrains it.")
log_warn("  That single number is the most useful spatial result here, and it")
log_warn("  is the argument for the next field season rather than for a")
log_warn("  better interpolator.")
log_warn(strrep("-", 68))

# =============================================================================
# 5. Coverage summary
# =============================================================================

log_step("09e  COVERAGE SUMMARY")

cell_km2 <- (CFG$spatial$cell_m / 1000)^2
bands <- c(1, 2, radius_km, 5, 10, 20)
bands <- sort(unique(round(bands, 2)))
cov <- do.call(rbind, lapply(bands, function(r) {
  data.frame(within_km = r,
             pct_of_aoi = 100 * mean(best_d <= r),
             area_km2   = sum(best_d <= r) * cell_km2,
             is_credible_radius = isTRUE(all.equal(r, round(radius_km, 2))),
             stringsAsFactors = FALSE)
}))
print_table(cov, digits = 2)

aoi_km2 <- g$nx * g$ny * cell_km2
log_info(sprintf("study area: %.0f km2 at %d m resolution (%d x %d cells)",
                 aoi_km2, CFG$spatial$cell_m, g$nx, g$ny))

summary_tbl <- rbind(
  data.frame(metric = "aoi_area_km2", value = round(aoi_km2, 1)),
  data.frame(metric = "cell_size_m", value = CFG$spatial$cell_m),
  data.frame(metric = "n_cores", value = n_core),
  data.frame(metric = "credible_radius_km", value = round(radius_km, 3)),
  data.frame(metric = "pct_aoi_within_credible_radius", value = round(100 * frac, 2)),
  data.frame(metric = "area_km2_within_credible_radius",
             value = round(sum(best_d <= radius_km) * cell_km2, 1)),
  data.frame(metric = "max_dist_to_core_km", value = round(max(best_d), 2))
)
write_csv_logged(summary_tbl,
                 file.path(CFG$dir_tables, "09_spatial_coverage.csv"),
                 "how much of the study area the cores can speak for")
write_csv_logged(cov, file.path(CFG$dir_tables, "09_coverage_by_distance.csv"))

# =============================================================================
# 6. Quick-look figures
# =============================================================================

log_step("09f  FIGURES")

if (!capabilities("png")) {
  log_warn("this R build has no PNG device; skipping figures (rasters are ",
           "written regardless)")
} else {
  # Sequential palette, colour-blind safe, dark = more carbon. Built by hand so
  # the workflow keeps its no-package promise.
  ramp_seq <- grDevices::colorRampPalette(
    c("#fff7ec", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"))
  ramp_dist <- grDevices::colorRampPalette(
    c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))

  #' Draw one raster with core locations over it. Rows run north to south, so
  #' the matrix is flipped for image(), which expects y to increase upward.
  # Cores sitting within 2 km of a neighbour get their labels alternated above
  # and below the point, so the tight PM cluster stays readable.
  dmat <- pairwise_km(cores$longitude, cores$latitude); diag(dmat) <- Inf
  crowded <- apply(dmat, 1, min) < 2
  lab_pos <- rep(3L, nrow(cores))
  lab_pos[crowded] <- c(3L, 1L)[seq_len(sum(crowded)) %% 2L + 1L]

  draw_map <- function(m, main, ramp, legend_title, path, log_scale = FALSE) {
    grDevices::png(path, width = 1400, height = 1400, res = 170)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(5.8, 4.4, 3.2, 1.2), mgp = c(2.6, 0.7, 0))

    v <- if (log_scale) log10(pmax(m, 1)) else m
    zl <- range(v, na.rm = TRUE)
    graphics::image(x = g$x, y = rev(g$y), z = t(v[nrow(v):1, ]),
                    col = ramp(200), zlim = zl, asp = 1 / cos(g$lat_mid * pi / 180),
                    xlab = "Longitude", ylab = "Latitude", main = main,
                    useRaster = TRUE)
    graphics::points(cores$longitude, cores$latitude, pch = 21,
                     bg = ifelse(grepl("peat", cores$campaign), "#1b9e77", "#7570b3"),
                     col = "white", cex = 1.5, lwd = 1.6)
    graphics::text(cores$longitude, cores$latitude, cores$core_id,
                   pos = lab_pos, offset = 0.45, cex = 0.55, col = "grey15")
    graphics::legend("topright", legend = c("PM peat (2024)", "FS mineral (2025)"),
                     pt.bg = c("#1b9e77", "#7570b3"), pch = 21, col = "white",
                     bg = "white", cex = 0.7, pt.cex = 1.2)
    graphics::mtext(legend_title, side = 1, line = 4.3, cex = 0.66, adj = 0,
                    col = "grey30")
    graphics::box()
    invisible(path)
  }

  f1 <- file.path(CFG$dir_figures, "09_dist_to_nearest_core.png")
  draw_map(best_d, "Distance to nearest community core",
           ramp_dist,
           sprintf(paste0("km. Darker = further from any measurement. ",
                          "Max %.1f km. Status: GEOMETRY (exact)."), max(best_d)),
           f1)
  log_ok("wrote ", basename(f1))

  f2 <- file.path(CFG$dir_figures, "09_soc_within_credible_radius.png")
  draw_map(soc_masked,
           sprintf("SOC 0-%.1f cm, within %.1f km of a core",
                   cs$window_hi_cm[1], radius_km),
           ramp_seq,
           sprintf(paste0("kg C/m2. Blank = unconstrained by these cores ",
                          "(%.1f%% of the area). Status: SAMPLE-BASED, MASKED."),
                   100 * (1 - frac)),
           f2)
  log_ok("wrote ", basename(f2))

  f3 <- file.path(CFG$dir_figures, "09_thiessen_diagnostic.png")
  draw_map(soc, "Nearest-core assignment (DIAGNOSTIC, not a carbon map)",
           ramp_seq,
           paste0("kg C/m2. Each pixel takes the closest core's value. Shown ",
                  "to expose the sampling gap, not to estimate carbon."),
           f3)
  log_ok("wrote ", basename(f3))

  # Core stocks by campaign and depth window: the numbers behind the maps.
  f4 <- file.path(CFG$dir_figures, "09_core_stocks.png")
  grDevices::png(f4, width = 1500, height = 900, res = 170)
  graphics::par(mar = c(6.5, 4.6, 3.0, 1.0), mgp = c(2.8, 0.7, 0))
  o   <- order(cores$campaign, -cs$stock_kgm2[match(cores$core_id, cs$core_id)])
  ids <- cores$core_id[o]
  v_c <- cs$stock_kgm2[match(ids, cs$core_id)]
  v_3 <- r30$stock_kgm2[match(ids, r30$core_id)]
  lb  <- r30$is_lower_bound[match(ids, r30$core_id)]
  bp <- graphics::barplot(rbind(v_c, v_3), beside = TRUE, names.arg = ids,
                          las = 2, col = c("#4292c6", "#f16913"),
                          border = NA, ylim = c(0, max(v_3, v_c) * 1.25),
                          ylab = "kg C / m2",
                          main = "Core carbon stocks by depth window")
  graphics::legend("topright", bty = "n", cex = 0.75, border = NA,
                   fill = c("#4292c6", "#f16913"),
                   legend = c(sprintf("0-%.1f cm (common support)", cs$window_hi_cm[1]),
                              "0-30 cm (reference)"))
  if (any(lb)) {
    graphics::text(bp[2, lb], v_3[lb], "LOWER\nBOUND", pos = 3,
                   cex = 0.52, col = "#a63603", font = 2)
  }
  graphics::mtext(paste0("Orange bars marked 'lower bound' are cores that ",
                         "stop short of 30 cm; their true stock is higher."),
                  side = 1, line = 5.1, cex = 0.62, adj = 0, col = "grey30")
  grDevices::dev.off()
  log_ok("wrote ", basename(f4))
}

# =============================================================================
# 7. Manifest
# =============================================================================

log_step("09f  MANIFEST")

manifest <- data.frame(
  file = c("cores.geojson", "aoi.geojson", "dist_to_nearest_core_m.tif",
           "nearest_core_index.tif", "DIAGNOSTIC_soc_nearest_core_kgm2.tif",
           "PRODUCT1b_soc_within_credible_radius_kgm2.tif"),
  type = c("vector point", "vector polygon", "raster Float32", "raster Int32",
           "raster Float32", "raster Float32"),
  units = c("mixed", "n/a", "metres", "core index (see legend)",
            "kg C/m2", "kg C/m2"),
  scientific_status = c(
    "OBSERVATION -- measured cores",
    "n/a -- study area outline",
    "GEOMETRY -- exact, no inference",
    "DIAGNOSTIC -- Thiessen assignment",
    "DIAGNOSTIC -- NOT a carbon map; nearest core's value, unmasked",
    "SAMPLE-BASED, MASKED -- nearest core's value within the credible radius only"),
  depth_support = c("see attributes", "n/a", "n/a", "n/a",
                    sprintf("0-%.1f cm (common support)", cs$window_hi_cm[1]),
                    sprintf("0-%.1f cm (common support)", cs$window_hi_cm[1])),
  stringsAsFactors = FALSE
)
print_table(manifest)
write_csv_logged(manifest, file.path(dir_sp, "MANIFEST.csv"),
                 "what each spatial file is and what it may be used for")

log_ok("09 complete -- ", length(list.files(dir_sp)), " files in outputs/spatial/")
