# =============================================================================
# 17_package_for_community.R   --   STEP 1 of the data-sharing plan
#
# Package the community cores as GIS point layers for a First Nation governing
# council: what was collected, where, with the quality-control record attached
# and nothing hidden.
#
# OUTPUTS
#   deliverables/step1_cores/fort_severn_cores.gpkg
#       layers: segments_raw (22 rows), cores_qc (8 rows)
#   deliverables/step1_cores/fort_severn_cores_shp.zip
#       ESRI Shapefile, zipped, ready to upload to Earth Engine
#   deliverables/step1_cores/*.csv
#       plain tables with EPSG:3978 coordinates, readable anywhere
#   deliverables/step1_cores/FIELD_DICTIONARY.csv
#       every column, its meaning, and its shapefile-truncated name
#   deliverables/step1_cores/README.md
#       what the layers are, what the flags mean, what is assumed
#
# CRS: EPSG:3978, NAD83 / Canada Atlas Lambert, throughout.
#
# TWO THINGS THIS SCRIPT REFUSES TO DO SILENTLY
#
#   1. Reproject without checking. Eight core locations are the foundation of
#      every downstream product. The transform is done by sf/PROJ AND by the
#      independent base-R implementation in R/proj.R, and the run stops if the
#      two disagree by more than a centimetre. A projection error would move
#      the cores plausibly -- a few hundred metres, right province, nothing to
#      notice.
#
#   2. Truncate field names by accident. Shapefiles cap column names at ten
#      characters and will silently mangle or collide them. Short names are
#      assigned explicitly here and the mapping is shipped as a dictionary.
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

log_step("17  STEP 1 -- PACKAGE THE CORES FOR THE COMMUNITY")

out_dir <- file.path(CFG$root, "deliverables", "step1_cores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seg    <- require_artifact(file.path(CFG$dir_derived, "01_segments_qc.csv"),
                           "scripts/01_ingest_qc.R")
cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

pd_path <- file.path(CFG$dir_tables, "11_peat_depth_at_cores.csv")
peatd <- if (file.exists(pd_path)) read_csv_verbatim(pd_path) else NULL

# =============================================================================
# 1. Geolocation gate, once more at the delivery boundary
# =============================================================================

log_step("17a  GEOLOCATION")

if (any(seg$longitude > 0) || any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the community deliverable",
              "These files are the ones that leave the project.",
              remedy = "Re-run 01_ingest_qc.R and inspect the geolocation gate.")
}
log_ok("all longitudes are western, as the study area requires")

p3978 <- lcc_params(3978)
log_info("target CRS: EPSG:", p3978$epsg, " -- ", p3978$name)
log_info("NOTE: EPSG:3978 (NAD83) and EPSG:3979 (NAD83(CSRS)) share identical ")
log_info("      projection parameters and differ only in datum realisation, ")
log_info("      under a metre in Canada. The plan specifies 3978; that is used ")
log_info("      here, while CFG$crs_equal_area (3979) is left as-is for the ")
log_info("      Earth Engine side. They are not interchangeable in metadata.")

#' Project with sf if available, always also with base R, and compare.
project_checked <- function(lon, lat, id) {
  base <- lonlat_to_lcc(lon, lat, p3978)
  if (!requireNamespace("sf", quietly = TRUE)) {
    log_warn("sf is not installed; using the base-R projection unchecked. ",
             "Install sf to enable the independent cross-check.")
    return(list(x = base$x, y = base$y, checked = FALSE, worst_m = NA_real_))
  }
  pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                      coords = c("lon", "lat"), crs = 4326)
  tr  <- sf::st_transform(pts, 3978)
  co  <- sf::st_coordinates(tr)
  cmp <- compare_projections(base$x, base$y, co[, 1], co[, 2], id)
  worst <- max(cmp$offset_m)
  if (worst > 0.01) {
    print_table(cmp[order(-cmp$offset_m), ], digits = 4)
    fail_loudly(
      "sf/PROJ and the base-R projection disagree",
      sprintf("Worst offset: %.3f m (tolerance 0.01 m)", worst),
      "One of the two is wrong, and the cores are the foundation of every",
      "downstream product.",
      remedy = "Check the EPSG code and the PROJ installation before shipping.")
  }
  log_ok(sprintf("sf/PROJ and base-R agree to %.4f m across %d points",
                 worst, length(lon)))
  list(x = co[, 1], y = co[, 2], checked = TRUE, worst_m = worst)
}

pc <- project_checked(cores$longitude, cores$latitude, cores$core_id)
ps <- project_checked(seg$longitude, seg$latitude, seg$sample_id)

cores$x_3978 <- pc$x; cores$y_3978 <- pc$y
seg$x_3978   <- ps$x; seg$y_3978   <- ps$y

# =============================================================================
# 2. Build the delivery tables
# =============================================================================

log_step("17b  LAYERS")

# --- core level: attach the stocks the council will actually ask about -------
cs  <- stocks[stocks$scenario == "as_supplied" & stocks$window == "common_support", ]
r30 <- stocks[stocks$scenario == "as_supplied" & stocks$window == "reference_0_30", ]

cores_out <- data.frame(
  core_id        = cores$core_id,
  campaign       = cores$campaign,
  soil_type      = ifelse(grepl("peat", cores$campaign), "organic/peat",
                          "mineral/forest"),
  year           = cores$year,
  latitude       = cores$latitude,
  longitude      = cores$longitude,
  x_3978         = cores$x_3978,
  y_3978         = cores$y_3978,
  n_segments     = cores$n_segments,
  core_cm        = round(cores$core_length_cm, 2),
  bd_min         = round(cores$bd_min_gcm3, 4),
  bd_max         = round(cores$bd_max_gcm3, 4),
  om_max_pct     = round(cores$om_max_pct, 2),
  stringsAsFactors = FALSE)

cores_out$soc_shal   <- round(cs$stock_kgm2[match(cores_out$core_id, cs$core_id)], 3)
cores_out$shal_cm    <- round(cs$window_hi_cm[match(cores_out$core_id, cs$core_id)], 2)
cores_out$soc_30cm   <- round(r30$stock_kgm2[match(cores_out$core_id, r30$core_id)], 3)
cores_out$soc30_lb   <- as.integer(r30$is_lower_bound[match(cores_out$core_id, r30$core_id)])
cores_out$depth30_cm <- round(r30$covered_cm[match(cores_out$core_id, r30$core_id)], 2)

if (!is.null(peatd)) {
  m <- match(cores_out$core_id, peatd$core_id)
  cores_out$peat_cm  <- round(peatd$peat_depth_cm[m], 2)
  cores_out$peat_sta <- peatd$status[m]
}
cores_out$n_flagged <- cores$n_flagged_segs
cores_out$fs_dist_km <- round(cores$dist_from_fort_severn_km, 3)

# --- segment level: the raw record, every flag preserved ---------------------
seg_out <- data.frame(
  sample_id   = seg$sample_id,
  core_id     = seg$core_id,
  campaign    = seg$campaign,
  year        = seg$year,
  seg_index   = seg$seg_index,
  latitude    = seg$latitude,
  longitude   = seg$longitude,
  x_3978      = seg$x_3978,
  y_3978      = seg$y_3978,
  top_cm      = round(seg$depth_top_cm, 3),
  bot_cm      = round(seg$depth_bottom_cm, 3),
  thick_cm    = round(seg$thickness_cm, 3),
  bd_gcm3     = round(seg$bulk_density_gcm3, 6),
  om_pct      = round(seg$om_pct, 3),
  soc_pct     = round(seg$soc_pct, 3),
  cdens_gcm3  = round(seg$carbon_density_gcm3, 8),
  ocd_supp    = round(seg$ocd_supplied, 8),
  ocd_ok      = as.integer(seg$ocd_reconciles),
  lon_fixed   = as.integer(seg$lon_repaired),
  qc_flags    = seg$qc_flags,
  qc_n        = seg$qc_n_flags,
  stringsAsFactors = FALSE)

log_ok("segments_raw: ", nrow(seg_out), " rows x ", ncol(seg_out), " fields")
log_ok("cores_qc:     ", nrow(cores_out), " rows x ", ncol(cores_out), " fields")

# Every field name is already <= 10 characters, so the shapefile writer has
# nothing to truncate. Verified rather than assumed.
long_s <- names(seg_out)[nchar(names(seg_out)) > 10]
long_c <- names(cores_out)[nchar(names(cores_out)) > 10]
if (length(long_s) || length(long_c)) {
  fail_loudly("Field names exceed the 10-character shapefile limit",
              paste0("segments: ", paste(long_s, collapse = ", ")),
              paste0("cores: ", paste(long_c, collapse = ", ")),
              remedy = "Shorten them here, deliberately, rather than letting the shapefile writer mangle them.")
}
log_ok("all field names are 10 characters or fewer; nothing will be truncated")

# =============================================================================
# 3. Field dictionary
# =============================================================================

log_step("17c  FIELD DICTIONARY")

dict <- rbind(
  data.frame(layer = "cores_qc", field = names(cores_out), stringsAsFactors = FALSE),
  data.frame(layer = "segments_raw", field = names(seg_out), stringsAsFactors = FALSE))
meanings <- c(
  core_id = "Core identifier as recorded in the field",
  sample_id = "Segment identifier; trailing integer gives order down the core",
  campaign = "Field campaign: PM_2024_peat or FS_2025_mineral",
  soil_type = "Plain-language soil type",
  year = "Year of collection",
  seg_index = "Position of this segment down the core (1 = surface)",
  latitude = "Latitude, decimal degrees, EPSG:4326",
  longitude = "Longitude, decimal degrees, EPSG:4326 (NEGATIVE = west)",
  x_3978 = "Easting, metres, EPSG:3978 NAD83 / Canada Atlas Lambert",
  y_3978 = "Northing, metres, EPSG:3978 NAD83 / Canada Atlas Lambert",
  n_segments = "Number of segments cut from this core",
  core_cm = "Total depth reached by this core, cm",
  top_cm = "Depth to the top of this segment, cm below the surface",
  bot_cm = "Depth to the bottom of this segment, cm below the surface",
  thick_cm = "Thickness of this segment, cm (the raw `Depth` column)",
  bd_gcm3 = "Dry bulk density, g/cm3",
  bd_min = "Minimum bulk density in this core, g/cm3",
  bd_max = "Maximum bulk density in this core, g/cm3",
  om_pct = "Organic matter, percent by dry mass",
  om_max_pct = "Maximum organic matter in this core, percent",
  soc_pct = "Organic carbon, percent by dry mass (see README on conversion)",
  cdens_gcm3 = "Carbon density, g C per cm3 of soil = bd_gcm3 x soc_pct/100",
  ocd_supp = "The supplied 'Organic Carbon Density' value, unaltered",
  ocd_ok = "1 if the supplied value reconciles with our recomputation",
  lon_fixed = "1 if this row's longitude sign was corrected on ingest",
  qc_flags = "Quality-control flags, semicolon separated (see README)",
  qc_n = "Number of quality-control flags on this row",
  n_flagged = "Number of flagged segments in this core",
  soc_shal = "Soil carbon over the shallow window every core reached, kg C/m2",
  shal_cm = "Depth of that shallow window, cm",
  soc_30cm = "Soil carbon over 0-30 cm, kg C/m2",
  soc30_lb = "1 if soc_30cm is a LOWER BOUND (core stopped before 30 cm)",
  depth30_cm = "Depth actually covered within the 0-30 cm window, cm",
  peat_cm = "Depth to the peat/mineral contact, cm",
  peat_sta = "How peat depth was determined: contact_observed, still_peat, no_peat",
  fs_dist_km = "Distance from Fort Severn, km")
dict$meaning <- meanings[dict$field]
dict$meaning[is.na(dict$meaning)] <- "(see project documentation)"
dict$units <- ifelse(grepl("_cm$|_cm ", dict$field), "cm",
              ifelse(grepl("gcm3", dict$field), "g/cm3",
              ifelse(grepl("^soc_shal|^soc_30cm", dict$field), "kg C/m2",
              ifelse(grepl("pct", dict$field), "percent",
              ifelse(grepl("x_3978|y_3978", dict$field), "metres",
              ifelse(grepl("km$", dict$field), "km", ""))))))
write_csv_logged(dict, file.path(out_dir, "FIELD_DICTIONARY.csv"),
                 "every field, its meaning and its units")

# =============================================================================
# 4. Write the spatial formats
# =============================================================================

log_step("17d  SPATIAL FORMATS")

written <- list()
have_sf <- requireNamespace("sf", quietly = TRUE)

if (!have_sf) {
  log_warn(strrep("-", 68))
  log_warn("`sf` IS NOT INSTALLED, so no GeoPackage or Shapefile is written.")
  log_warn("  install.packages('sf')  then re-run this script.")
  log_warn("  The CSV and GeoJSON below carry the same data and coordinates,")
  log_warn("  so nothing is lost, but the plan asks for a GeoPackage.")
  log_warn(strrep("-", 68))
} else {
  to_sf <- function(df) sf::st_as_sf(df, coords = c("x_3978", "y_3978"),
                                     crs = 3978, remove = FALSE)
  sf_cores <- to_sf(cores_out)
  sf_seg   <- to_sf(seg_out)

  gpkg <- file.path(out_dir, "fort_severn_cores.gpkg")
  if (file.exists(gpkg)) unlink(gpkg)
  sf::st_write(sf_cores, gpkg, layer = "cores_qc", quiet = TRUE)
  sf::st_write(sf_seg,   gpkg, layer = "segments_raw", quiet = TRUE,
               append = TRUE)
  log_ok("wrote fort_severn_cores.gpkg (layers: cores_qc, segments_raw)")
  written[[length(written) + 1L]] <- gpkg

  # Shapefile, zipped, for Earth Engine upload. One shapefile per layer:
  # the format has no concept of multiple layers in one file.
  shp_dir <- file.path(out_dir, "shp")
  dir.create(shp_dir, showWarnings = FALSE)
  sf::st_write(sf_cores, file.path(shp_dir, "fort_severn_cores_qc.shp"),
               quiet = TRUE, delete_dsn = TRUE)
  sf::st_write(sf_seg, file.path(shp_dir, "fort_severn_segments_raw.shp"),
               quiet = TRUE, delete_dsn = TRUE)

  zipf <- file.path(out_dir, "fort_severn_cores_shp.zip")
  if (file.exists(zipf)) unlink(zipf)
  old <- setwd(shp_dir)
  ok_zip <- tryCatch({
    utils::zip(zipf, list.files(".") , flags = "-q"); TRUE
  }, error = function(e) { log_warn("zip failed: ", conditionMessage(e)); FALSE })
  setwd(old)
  if (ok_zip && file.exists(zipf)) {
    log_ok("wrote fort_severn_cores_shp.zip (Earth Engine upload)")
    written[[length(written) + 1L]] <- zipf
  } else {
    log_warn("zip unavailable; the unzipped shapefiles are in ", shp_dir)
  }

  # Read the GeoPackage back and confirm what a GIS will actually see.
  chk <- sf::st_read(gpkg, layer = "cores_qc", quiet = TRUE)
  epsg_back <- sf::st_crs(chk)$epsg
  if (!identical(as.integer(epsg_back), 3978L)) {
    fail_loudly("The GeoPackage did not record EPSG:3978",
                paste0("It reports EPSG:", epsg_back),
                remedy = "Check the sf/GDAL installation before shipping.")
  }
  if (nrow(chk) != nrow(cores_out)) {
    fail_loudly("The GeoPackage lost rows on write",
                sprintf("wrote %d, read back %d", nrow(cores_out), nrow(chk)))
  }
  log_ok("read-back check: EPSG:", epsg_back, ", ", nrow(chk), " features")
}

# --- formats that always work ------------------------------------------------
write_csv_logged(cores_out, file.path(out_dir, "fort_severn_cores_qc.csv"),
                 "core level, with EPSG:3978 coordinates")
write_csv_logged(seg_out, file.path(out_dir, "fort_severn_segments_raw.csv"),
                 "segment level, with EPSG:3978 coordinates")

# GeoJSON is EPSG:4326 by specification, so it carries the geographic
# coordinates; the projected values travel as attributes.
write_geojson_points(cores_out, file.path(out_dir, "fort_severn_cores_qc.geojson"))
write_geojson_points(seg_out, file.path(out_dir, "fort_severn_segments_raw.geojson"))
log_ok("wrote GeoJSON for both layers (EPSG:4326 per the GeoJSON spec)")

# =============================================================================
# 5. README for the recipient
# =============================================================================

log_step("17e  README")

fl <- read_csv_verbatim(file.path(CFG$dir_tables, "01_flag_dictionary.csv"))
raised <- sort(unique(unlist(strsplit(seg$qc_flags[nzchar(seg$qc_flags)], ";"))))
fl_used <- fl[fl$flag %in% raised, ]

rd <- c(
  "# Fort Severn community soil cores — GIS layers",
  "",
  paste0("Prepared ", format(Sys.Date(), "%d %B %Y"), " for the Fort Severn ",
         "First Nation governing council."),
  "",
  "## What this is",
  "",
  paste0("The complete record of ", nrow(cores_out), " soil cores collected by ",
         "community members across two field seasons, as GIS point layers. ",
         "Nothing has been removed. Values that looked doubtful are kept and ",
         "marked with a flag, so you can see exactly what was measured and ",
         "what we thought about it."),
  "",
  "## Files",
  "",
  "| File | What it is |",
  "|------|------------|",
  "| `fort_severn_cores.gpkg` | GeoPackage, opens in ArcGIS and QGIS. Two layers: `cores_qc` (one point per core) and `segments_raw` (one point per sampled slice). |",
  "| `fort_severn_cores_shp.zip` | The same data as Shapefiles, zipped for uploading to Google Earth Engine. |",
  "| `*.csv` | Plain tables. Open in Excel or anything else. |",
  "| `*.geojson` | Web-mapping format. |",
  "| `FIELD_DICTIONARY.csv` | Every column explained, with its units. |",
  "",
  paste0("**Coordinate system: EPSG:3978 (NAD83 / Canada Atlas Lambert).** ",
         "Latitude and longitude are also included in every file."),
  "",
  "## The two layers",
  "",
  paste0("**`cores_qc`** — one point per core (", nrow(cores_out), " points). ",
         "Use this for maps and for totals."),
  "",
  paste0("**`segments_raw`** — one point per slice cut from a core (",
         nrow(seg_out), " points). Several slices share a location because ",
         "they came from the same hole, one below the other. Use this to look ",
         "at how carbon changes with depth."),
  "",
  "## Two numbers for carbon, and why",
  "",
  paste0("Cores stopped at different depths, between ",
         sprintf("%.1f", min(cores_out$core_cm)), " and ",
         sprintf("%.1f", max(cores_out$core_cm)), " cm. To compare them ",
         "fairly you must compare the same depth in each."),
  "",
  paste0("- **`soc_shal`** — carbon over the top ",
         sprintf("%.1f", cores_out$shal_cm[1]), " cm, the deepest layer ",
         "*every* core reached. Every core has a real number here, so this is ",
         "the fair comparison."),
  paste0("- **`soc_30cm`** — carbon over the top 30 cm, the depth published ",
         "maps use. Two cores never got that deep. For those, the number is a ",
         "**minimum**, not a measurement, and `soc30_lb = 1` marks them."),
  "",
  "Never mix the two columns in one total.",
  "",
  "## What the flags mean",
  "",
  "| Flag | Meaning |",
  "|------|---------|",
  paste0("| `", fl_used$flag, "` | ", fl_used$meaning, " |"),
  "",
  "## Two things you should know about the data",
  "",
  paste0("**1. The longitudes were recorded without a minus sign.** As ",
         "supplied, the coordinates placed these cores in Siberia, about ",
         "7,540 km away. This was corrected on loading, and every corrected ",
         "row is marked `lon_fixed = 1`. After correction the cores sit ",
         sprintf("%.1f to %.1f km", min(cores_out$fs_dist_km),
                 max(cores_out$fs_dist_km)),
         " from Fort Severn, which confirms the correction was right."),
  "",
  paste0("**2. The depth ordering is an assumption that needs a field check.** ",
         "The `Depth` column records the THICKNESS of each slice, not its ",
         "depth below the surface. We established this from the data itself: ",
         "read as depth-below-surface the values would decrease going down two ",
         "cores, which is impossible. Slice order comes from the number at the ",
         "end of the sample name."),
  "",
  paste0("> **Please confirm against the field notes** that slice ",
         "`PM-03-A-1` is the surface slice and `PM-03-A-4` the deepest, and ",
         "likewise for the other cores. Every depth in these files, and every ",
         "carbon total, depends on that ordering being right."),
  "",
  "## How carbon was calculated",
  "",
  "```",
  "carbon in a slice  =  bulk density  x  carbon %  /  100  x  slice thickness",
  "     (kg C/m2)          (g/cm3)                              (cm)",
  "",
  "core total  =  sum of all its slices",
  "```",
  "",
  paste0("One caution on the carbon percentages: the two field seasons used ",
         "different laboratory conversion factors to turn organic matter into ",
         "organic carbon (0.58 for the 2024 peat cores, 0.554 for the 2025 ",
         "mineral cores). This means a peat-versus-mineral comparison is ",
         "partly a laboratory difference, not only a difference in the ground. ",
         "Using one method for both, or measuring carbon directly on a few ",
         "samples, would settle it."),
  "",
  "## Questions",
  "",
  paste0("Every number in these files is produced by open scripts from the ",
         "original field data, and can be regenerated from scratch. If a value ",
         "looks wrong, it can be traced back to the exact measurement it came ",
         "from."))

writeLines(rd, file.path(out_dir, "README.md"))
log_ok("wrote README.md (", length(rd), " lines)")

# =============================================================================
# 6. Manifest
# =============================================================================

log_step("17f  DELIVERABLE MANIFEST")

files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
man <- data.frame(
  file = sub(paste0(out_dir, "/?"), "", files),
  size_kb = round(file.size(files) / 1024, 1),
  stringsAsFactors = FALSE)
man <- man[order(man$file), ]
print_table(man)
write_csv_logged(man, file.path(out_dir, "MANIFEST.csv"),
                 "everything in this delivery")

log_step("17  STEP 1 STATUS")
log_ok("point layers exported: GeoPackage",
       if (have_sf) " + zipped Shapefile" else " (SKIPPED - no sf)",
       " + CSV + GeoJSON")
log_ok("CRS confirmed EPSG:3978",
       if (pc$checked) sprintf(" (sf and base-R agree to %.4f m)", pc$worst_m)
       else " (base-R only; install sf for the cross-check)")
log_warn("STILL OPEN: the depth-ordering assumption needs confirming against ",
         "the field notes. It is flagged in the README, and every depth and ",
         "carbon total depends on it.")
log_ok("17 complete -- deliverables in ", out_dir)
