# =============================================================================
# tests/test_functions.R  --  Unit tests for the pure-function layer.
#
# Base R only. Run with:  Rscript tests/test_functions.R
# Exits non-zero if anything fails.
# =============================================================================

.root <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  d <- if (length(f)) dirname(f[[1]]) else getwd()
  while (!file.exists(file.path(d, "config.R"))) {
    p <- dirname(d)
    if (identical(p, d)) stop("config.R not found", call. = FALSE)
    d <- p
  }
  normalizePath(d)
})
source(file.path(.root, "config.R"))

# ---- micro test framework ----------------------------------------------------
.pass <- 0L; .fail <- 0L; .failures <- character(0)

ok <- function(cond, what) {
  if (isTRUE(cond)) {
    .pass <<- .pass + 1L
    cat(sprintf("  ok   %s\n", what))
  } else {
    .fail <<- .fail + 1L
    .failures <<- c(.failures, what)
    cat(sprintf("  FAIL %s\n", what))
  }
}
eq <- function(a, b, what, tol = 1e-9) {
  ok(length(a) == length(b) && all(abs(a - b) < tol), what)
}
throws <- function(expr, what) {
  ok(inherits(try(force(expr), silent = TRUE), "try-error"), what)
}
section <- function(x) cat("\n", x, "\n", sep = "")

# =============================================================================
section("trailing_int")
# =============================================================================
eq(trailing_int("PM-03-A-2"), 2L, "extracts trailing integer")
eq(trailing_int("FS-05-1"), 1L, "extracts single digit")
eq(trailing_int("PM-02-A-10"), 10L, "extracts multi-digit")
ok(is.na(trailing_int("PM-01-A")), "returns NA when no trailing integer")
eq(trailing_int(c("A-1", "B-22")), c(1L, 22L), "vectorised")

# =============================================================================
section("carbon_density_gcm3 -- the first-principles definition")
# =============================================================================
eq(carbon_density_gcm3(0.0929, 22.7), 0.0210883, "matches the supplied column exactly")
eq(carbon_density_gcm3(2.0, 0), 0, "zero carbon gives zero density")
eq(carbon_density_gcm3(1.0, 100), 1.0, "100 percent carbon gives BD")

# =============================================================================
section("integrate_window")
# =============================================================================
# One 10 cm segment at 0.05 g/cm3 over 0-10 => 0.5 g/cm2 => 5 kg/m2
r <- integrate_window(0, 10, 0.05, 0, 10)
eq(r$stock_gcm2, 0.5, "single segment integrates to density x thickness")
eq(r$stock_kgm2, 5.0, "g/cm2 -> kg/m2 conversion is x10")
ok(!r$is_lower_bound, "full coverage is not a lower bound")
eq(r$coverage_frac, 1, "coverage fraction is 1")

# Truncation: a segment extending past the window contributes only the overlap.
r <- integrate_window(c(0, 10), c(10, 25), c(0.1, 0.2), 0, 20)
eq(r$stock_gcm2, 0.1 * 10 + 0.2 * 10, "straddling segment is cut at the boundary")
eq(r$covered_cm, 20, "covered depth equals the window")
ok(!r$is_lower_bound, "truncated but complete is not a lower bound")

# Shortfall: core stops before the window closes.
r <- integrate_window(0, 14.5, 0.1, 0, 30)
ok(r$is_lower_bound, "short core is flagged as a lower bound")
eq(r$covered_cm, 14.5, "covered depth reports the shortfall")
eq(r$stock_gcm2, 1.45, "shortfall stock is NOT rescaled to the full window")

# Window not starting at zero.
r <- integrate_window(c(0, 10), c(10, 20), c(1, 2), 5, 15)
eq(r$stock_gcm2, 1 * 5 + 2 * 5, "window offset from the surface integrates correctly")

# A segment entirely outside the window contributes nothing.
r <- integrate_window(c(0, 40), c(10, 50), c(1, 99), 0, 10)
eq(r$stock_gcm2, 10, "segment outside the window is excluded")
eq(r$n_segments_used, 1, "only overlapping segments are counted")

# =============================================================================
section("add_depth_bounds")
# =============================================================================
d <- data.frame(core_id = c("A", "A", "A"), seg_index = c(1, 2, 3),
                thickness_cm = c(8.625, 8.75, 7.85), stringsAsFactors = FALSE)
b <- add_depth_bounds(d)
eq(b$depth_top_cm, c(0, 8.625, 17.375), "tops accumulate from zero")
eq(b$depth_bottom_cm, c(8.625, 17.375, 25.225), "bottoms accumulate thickness")

# Order must follow seg_index, not row order.
d2 <- d[c(3, 1, 2), ]
b2 <- add_depth_bounds(d2)
eq(b2$depth_bottom_cm, c(8.625, 17.375, 25.225), "reorders by segment index")

# =============================================================================
section("common_support_depth")
# =============================================================================
cs <- common_support_depth(c("A", "A", "B", "C"), c(10, 20, 14.5, 30))
eq(cs$depth_cm, 14.5, "common support is the shallowest core")
ok(cs$limiting_core == "B", "identifies the limiting core")

# =============================================================================
section("infer_depth_semantics")
# =============================================================================
s <- infer_depth_semantics(c("A", "A"), c(1, 2), c(15.2, 5.5))
ok(s$verdict == "thickness_cm", "a strict decrease falsifies a cumulative reading")
ok(s$n_decreases == 1, "counts the decrease")

s <- infer_depth_semantics(c("A", "A", "A"), 1:3, c(10, 10, 10))
ok(s$verdict == "thickness_cm", "repeated values also falsify a cumulative reading")
ok(s$n_ties == 2 && s$n_decreases == 0, "classifies ties separately from decreases")

s <- infer_depth_semantics(c("A", "A"), c(1, 2), c(10, 20))
ok(s$verdict == "ambiguous", "monotonic increase alone is genuinely ambiguous")

# =============================================================================
section("extrapolate_window")
# =============================================================================
e <- extrapolate_window(0, 10, 0.1, 0, 30)
eq(e$gap_cm, 20, "gap is the uncovered depth")
eq(e$stock_extrap_kgm2, (0.1 * 10 + 0.1 * 20) * 10, "carries the deepest density down")
eq(e$extrap_share_of_stock, 2 / 3, "reports the invented share of the stock")

e <- extrapolate_window(0, 30, 0.1, 0, 30)
eq(e$gap_cm, 0, "no gap when the core is complete")
ok(is.na(e$density_carried_gcm3), "nothing is carried when there is no gap")

# =============================================================================
section("gate_geolocation")
# =============================================================================
bbox <- c(xmin = -88.6, ymin = 55.5, xmax = -86.9, ymax = 56.6)
anchor <- c(lon = -87.6333, lat = 56.0167)

g <- gate_geolocation(c(87.68, 87.65), c(56.09, 56.10), anchor, bbox, TRUE)
eq(g$longitude, c(-87.68, -87.65), "repairs a positive western longitude")
eq(g$n_repaired, 2, "counts repairs")
ok(all(vapply(g$flags, function(f) "LON_SIGN_REPAIRED" %in% f, logical(1))),
   "flags every repaired row")

throws(gate_geolocation(c(87.68), c(56.09), anchor, bbox, FALSE),
       "STOPS when the sign is wrong and repair is not authorised")

g <- gate_geolocation(c(-87.68), c(56.09), anchor, bbox, TRUE)
eq(g$n_repaired, 0, "leaves an already-correct longitude alone")
ok(length(g$flags[[1]]) == 0, "raises no flag for a valid coordinate")

g <- gate_geolocation(c(-120.0), c(49.0), anchor, bbox, TRUE)
ok("COORD_OUT_OF_REGION" %in% g$flags[[1]],
   "flags a coordinate that is wrong in a way sign-flipping cannot fix")
eq(g$n_repaired, 0, "does not repair what a flip would not fix")

g <- gate_geolocation(c(NA_real_), c(56.0), anchor, bbox, TRUE)
ok("COORD_MISSING" %in% g$flags[[1]], "flags a missing coordinate")

# =============================================================================
section("gate_ocd_reconciliation")
# =============================================================================
r <- gate_ocd_reconciliation(c(1, 2, 3), c(1, 2, 3), 1e-6)
eq(r$n_bad, 0, "identical values reconcile")

r <- gate_ocd_reconciliation(c(1, 2, 3), c(1, 2, 3.0000001), 1e-6)
eq(r$n_bad, 0, "tolerance absorbs floating-point noise")

# A single bad row is flagged, not fatal: it may be one transcription error.
r <- gate_ocd_reconciliation(c(1, 2, 99), c(1, 2, 3), 1e-6)
eq(r$n_bad, 1, "one disagreement is reported")
ok(!r$ok[3], "the offending row is marked")

# A systematic failure means the column means something else, and must stop.
throws(gate_ocd_reconciliation(c(10, 20, 30), c(1, 2, 3), 1e-6),
       "STOPS when most rows disagree (the column means something else)")

# =============================================================================
section("identify_om_carbon_factor")
# =============================================================================
om  <- c(10, 20, 30, 40)
f <- identify_om_carbon_factor(om, om * 0.58, rep("peat", 4))
eq(f$factor_est, 0.58, "recovers a known factor exactly", tol = 1e-12)
ok(name_om_carbon_factor(f$factor_est) == "van Bemmelen (0.58 = 1/1.724)",
   "names the van Bemmelen factor")

f <- identify_om_carbon_factor(c(om, om), c(om * 0.58, om * 0.5),
                               rep(c("a", "b"), each = 4))
eq(sort(f$factor_est), c(0.5, 0.58), "separates factors by group")

# =============================================================================
section("flag_cross_core_duplicates")
# =============================================================================
f <- flag_cross_core_duplicates(c("A", "B", "C"), c(39.1, 39.1, 1), c(22.7, 22.7, 2))
ok("DUPLICATE_LAB_VALUES" %in% f[[1]] && "DUPLICATE_LAB_VALUES" %in% f[[2]],
   "flags an identical OM/SOC pair shared across two cores")
ok(length(f[[3]]) == 0, "leaves a unique pair unflagged")

f <- flag_cross_core_duplicates(c("A", "A"), c(39.1, 39.1), c(22.7, 22.7))
ok(length(f[[1]]) == 0,
   "does NOT flag repeats within a single core (legitimately possible)")

# =============================================================================
section("metrics")
# =============================================================================
eq(rmse(c(1, 2, 3), c(1, 2, 3)), 0, "rmse of a perfect fit is zero")
eq(r2_cv(c(1, 2, 3), c(1, 2, 3)), 1, "r2 of a perfect fit is one")
ok(r2_cv(c(1, 2, 3), c(3, 2, 1)) < 0,
   "r2 goes NEGATIVE when the model is worse than the mean")
eq(r2_cv(c(1, 2, 3), rep(2, 3)), 0, "predicting the mean gives r2 of zero")
eq(skill_vs_null(c(1, 2, 3), rep(2, 3)), 0, "null model has zero skill")

# =============================================================================
section("cross-validation grouping -- the leakage test")
# =============================================================================
# Two cores, three segments each. The covariate is constant within a core, so
# it cannot distinguish segments -- exactly the real situation at 30 m.
set.seed(1)
d <- data.frame(x = rep(c(0, 1), each = 3))
yy <- rep(c(10, 20), each = 3) + rnorm(6, 0, 0.01)
grp_core <- rep(c("A", "B"), each = 3)
grp_seg  <- as.character(1:6)

cv_seg  <- cv_leave_one_group_out(d, yy, grp_seg,  fit_ols, predict_ols)
cv_core <- cv_leave_one_group_out(d, yy, grp_core, fit_ols, predict_ols)

ok(cv_seg$metrics$r2_cv > 0.99,
   "segment-level CV reports near-perfect skill (this is the leak)")
ok(cv_core$metrics$r2_cv < cv_seg$metrics$r2_cv,
   "core-level CV is strictly less optimistic than segment-level")
ok(cv_core$n_folds == 2 && cv_seg$n_folds == 6,
   "fold counts follow the grouping unit")

# =============================================================================
section("boot_mean_ci determinism")
# =============================================================================
a <- boot_mean_ci(c(1, 2, 3, 4, 5), B = 500, seed = 42)
b <- boot_mean_ci(c(1, 2, 3, 4, 5), B = 500, seed = 42)
eq(c(a$lo, a$hi), c(b$lo, b$hi), "same seed gives identical intervals")
ok(!boot_mean_ci(c(1, 2, 3), B = 500, seed = 1)$trustworthy,
   "n=3 interval is marked untrustworthy")
ok(a$lo <= a$mean && a$mean <= a$hi, "interval brackets the mean")

# =============================================================================
section("haversine")
# =============================================================================
eq(haversine_km(0, 0, 0, 0), 0, "zero distance to self")
ok(abs(haversine_km(-87.6, 56.0, -87.6, 57.0) - 111.2) < 1,
   "one degree of latitude is about 111 km")
ok(haversine_km(87.68, 56.1, -87.68, 56.1) > 5000,
   "the sign error really is thousands of km (why it is a hard gate)")

# =============================================================================
section("near_rel")
# =============================================================================
ok(near_rel(1.0000001, 1, 1e-6), "accepts within relative tolerance")
ok(!near_rel(1.1, 1, 1e-6), "rejects outside relative tolerance")
ok(near_rel(0, 0, 1e-6), "handles zero reference")

# =============================================================================
section("make_grid")
# =============================================================================
bb <- c(xmin = -88.05, ymin = 55.88, xmax = -87.45, ymax = 56.22)
g <- make_grid(bb, 100)
ok(g$nx > 0 && g$ny > 0, "grid has positive dimensions")
ok(abs(g$x[1] - (bb[["xmin"]] + g$xres / 2)) < 1e-12,
   "first column is a cell CENTRE, half a cell east of the western edge")
ok(abs(g$y[1] - (bb[["ymax"]] - g$yres / 2)) < 1e-12,
   "first row is a cell centre, half a cell south of the northern edge")
ok(all(diff(g$y) < 0), "rows run north to south (raster orientation)")
ok(all(diff(g$x) > 0), "columns run west to east")
# Cells should be roughly square on the ground at this latitude.
ew <- haversine_km(g$x[1], g$lat_mid, g$x[2], g$lat_mid) * 1000
ns <- haversine_km(0, g$y[1], 0, g$y[2]) * 1000
ok(abs(ew - 100) < 5 && abs(ns - 100) < 5,
   "cells are ~100 m on the ground in both directions")

# =============================================================================
section("write_geotiff -- structural checks (no GDAL needed)")
# =============================================================================
tf <- tempfile(fileext = ".tif")
zz <- matrix(as.numeric(1:35), nrow = 5, ncol = 7, byrow = TRUE)
write_geotiff(zz, tf, xmin = -88.05, ymax = 56.22, xres = 0.01, yres = 0.005)
ok(file.exists(tf), "writes a file")

raw_hdr <- readBin(tf, "raw", n = 8)
ok(rawToChar(raw_hdr[1:2]) == "II", "little-endian byte order marker")
ok(readBin(raw_hdr[3:4], "integer", size = 2, endian = "little") == 42L,
   "TIFF magic number 42")
ifd_off <- readBin(raw_hdr[5:8], "integer", size = 4, endian = "little")
ok(ifd_off == 8L, "IFD begins immediately after the header")

all_raw <- readBin(tf, "raw", n = file.size(tf))
# TIFF SHORTs are UNSIGNED. R's readBin defaults to signed, which would turn
# every tag above 32767 (33550, 34735, 42113 ...) into a negative number.
u16 <- function(r) readBin(r, "integer", size = 2, endian = "little", signed = FALSE)
n_ent <- u16(all_raw[9:10])
ok(n_ent == 16L, "expected number of IFD entries")

tags <- vapply(seq_len(n_ent), function(i) {
  o <- 10 + (i - 1) * 12
  u16(all_raw[(o + 1):(o + 2)])
}, integer(1))
ok(!is.unsorted(tags), "IFD entries are in ascending tag order (TIFF requires it)")
ok(all(c(256L, 257L, 258L, 273L, 339L, 33550L, 33922L, 34735L) %in% tags),
   "all mandatory TIFF and GeoTIFF tags are present")

# File must be exactly header + IFD + payloads + one strip of 4-byte samples.
ok(file.size(tf) > 35 * 4, "file is at least as large as its pixel payload")
unlink(tf)

# Round-trip the pixel payload: read the last 35*4 bytes back as float32.
tf2 <- tempfile(fileext = ".tif")
write_geotiff(zz, tf2, xmin = 0, ymax = 10, xres = 1, yres = 1)
sz <- file.size(tf2)
con <- file(tf2, "rb"); invisible(seek(con, sz - 35 * 4))
back <- readBin(con, "double", n = 35, size = 4, endian = "little"); close(con)
eq(back, as.vector(t(zz)), "pixels round-trip in row-major order from the NW corner")
unlink(tf2)

# =============================================================================
section(".geo_key_directory")
# =============================================================================
gk <- .geo_key_directory(4326, geographic = TRUE)
ok(gk[1] == 1L && gk[2] == 1L, "GeoKeyDirectory version header")
ok(gk[4] == (length(gk) - 4) / 4, "declared key count matches the payload")
keyids <- gk[seq(5, length(gk), by = 4)]
ok(!is.unsorted(keyids), "geo keys are in ascending KeyID order")
ok(2048L %in% keyids, "GeographicTypeGeoKey present")
ok(gk[which(keyids == 2048L) * 4 + 4] == 4326L, "EPSG code is written as 4326")

# =============================================================================
section("GeoJSON writers")
# =============================================================================
gj <- tempfile(fileext = ".geojson")
dd <- data.frame(longitude = c(-87.6, -87.7), latitude = c(56.0, 56.1),
                 core_id = c("A", "B"), val = c(1.5, NA),
                 flag = c(TRUE, FALSE), stringsAsFactors = FALSE)
write_geojson_points(dd, gj)
txt <- paste(readLines(gj), collapse = "")
ok(grepl('"type": *"FeatureCollection"', txt), "writes a FeatureCollection")
ok(lengths(regmatches(txt, gregexpr('"type": *"Point"', txt))) == 2,
   "one Point feature per row")
ok(grepl("-87.6", txt) && grepl("56.1", txt), "coordinates are present")
ok(grepl('"val": *null', txt), "NA becomes JSON null, not the string NA")
ok(grepl('"flag": *true', txt), "logicals become JSON booleans")
unlink(gj)

gb <- tempfile(fileext = ".geojson")
write_geojson_bbox(bb, gb, props = list(name = "test"))
tb <- paste(readLines(gb), collapse = "")
ok(grepl('"Polygon"', tb), "bbox writes a Polygon")
ok(lengths(regmatches(tb, gregexpr("\\[-88\\.05", tb))) >= 2,
   "polygon ring is closed (first vertex repeats)")
unlink(gb)

# ---- summary -----------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(sprintf("passed %d   failed %d\n", .pass, .fail))
if (.fail > 0L) {
  cat("failures:\n"); for (f in .failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("all tests passed\n")
