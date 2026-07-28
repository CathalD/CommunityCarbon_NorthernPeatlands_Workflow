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
section("fit_ols degeneracy handling")
# =============================================================================
# A covariate that is constant within the training fold makes the design matrix
# rank-deficient. R would return a fit and predict() would warn about "doubtful
# cases" while still handing back a number. fit_ols must drop the column.
xc <- data.frame(const = rep(2, 5), good = c(1, 2, 3, 4, 5))
yc <- c(2, 4, 6, 8, 10)
mc <- fit_ols(xc, yc)
ok(identical(attr(mc, "dropped_constant"), "const"),
   "drops a constant predictor and records which one")
eq(as.numeric(predict_ols(mc, xc)), yc,
   "still fits perfectly on the informative predictor")

# Every predictor constant: nothing usable is left.
xd <- data.frame(a = rep(1, 4), b = rep(9, 4))
yd <- c(1, 2, 3, 4)
md <- fit_ols(xd, yd)
ok(inherits(md, "ols_degenerate"), "falls back to intercept-only when nothing varies")
eq(as.numeric(predict_ols(md, xd)), rep(mean(yd), 4),
   "degenerate model predicts the training mean")

# The whole point: no warning escapes to the user.
wr <- withCallingHandlers({
  invisible(predict_ols(fit_ols(xc, yc), xc)); "clean"
}, warning = function(w) { "warned" })
ok(wr == "clean", "no rank-deficiency warning reaches the caller")

# And it survives a real leave-one-out loop where a fold goes constant.
xe <- data.frame(v = c(5, 5, 5, 5, 1, 2))   # constant in most LOO folds
ye <- c(1, 2, 3, 4, 5, 6)
ge <- as.character(1:6)
cve <- withCallingHandlers(
  cv_leave_one_group_out(xe, ye, ge, fit_ols, predict_ols),
  warning = function(w) { .fail <<- .fail + 1L; invokeRestart("muffleWarning") })
ok(all(is.finite(cve$predictions$pred)), "every fold still returns a prediction")

# =============================================================================
section("skill_vs_reference")
# =============================================================================
o <- c(1, 2, 3, 4)
ok(abs(skill_vs_reference(o, o, rep(mean(o), 4)) - 1) < 1e-12,
   "perfect prediction has skill 1 against any reference")
eq(skill_vs_reference(o, rep(mean(o), 4), rep(mean(o), 4)), 0,
   "identical predictions have zero skill against each other")
ok(skill_vs_reference(o, c(4, 3, 2, 1), rep(mean(o), 4)) < 0,
   "a worse-than-reference prediction has negative skill")

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

# =============================================================================
section("detect_peat_depth")
# =============================================================================
d <- detect_peat_depth(c(0, 15.2), c(15.2, 20.7), c(39.1, 2.71))
eq(d$peat_depth_cm, 15.2, "contact placed at the top of the first mineral segment")
ok(d$status == "contact_observed", "observed contact is reported as such")
ok(!d$is_lower_bound, "an observed contact is not a lower bound")

d <- detect_peat_depth(0, 14.5, 39.1)
eq(d$peat_depth_cm, 14.5, "core ending in peat gives the core bottom")
ok(d$status == "still_peat" && d$is_lower_bound,
   "core ending in peat is flagged as a LOWER BOUND on peat depth")

d <- detect_peat_depth(c(0, 10), c(10, 20), c(5, 3))
ok(d$status == "no_peat" && d$peat_depth_cm == 0,
   "a wholly mineral core has zero peat depth")

# =============================================================================
section("kernel functions")
# =============================================================================
eq(kernel_gaussian(0, 2.5), 1, "weight is 1 at zero distance")
eq(kernel_gaussian(2.5, 2.5), exp(-0.5), "weight is exp(-1/2) at the length scale")
ok(kernel_gaussian(10, 2.5) < 0.001, "weight is negligible at 4 length scales")
ok(all(diff(kernel_gaussian(seq(0, 20, 1), 2.5)) < 0), "weight decreases monotonically")
throws(kernel_gaussian(1, 0), "rejects a non-positive length scale")

w <- kernel_truncated(c(0, 2, 6, 12), 2.5, max_km = 7.5)
ok(all(w[1:3] > 0) && w[4] == 0, "truncation zeroes weights beyond max_km")
ok(all(w >= 0 & w <= 1), "truncated weights stay within [0, 1]")

# =============================================================================
section("bayes_posterior")
# =============================================================================
# With zero weight everywhere the posterior must equal the prior exactly.
W0 <- matrix(0, 3, 2)
p <- bayes_posterior(c(10, 20, 30), rep(5, 3), c(1, 2), c(1, 1), W0)
eq(p$mean, c(10, 20, 30), "no core influence leaves the prior untouched")
eq(p$sd, rep(5, 3), "no core influence leaves the prior sd untouched")
eq(p$info_frac, rep(0, 3), "info_frac is 0 where cores have no weight")
eq(p$shift, rep(0, 3), "shift is 0 where cores have no weight")

# Full weight, equal variances: posterior mean is the midpoint.
W1 <- matrix(1, 1, 1)
p <- bayes_posterior(10, 1, 20, 1, W1)
eq(p$mean, 15, "equal precisions give the midpoint of prior and observation")
eq(p$sd, 1 / sqrt(2), "two equal precisions halve the variance")
eq(p$info_frac, 0.5, "info_frac is 0.5 when prior and data are equally precise")

# A very precise observation should dominate a vague prior.
p <- bayes_posterior(10, 100, 20, 0.01, matrix(1, 1, 1))
ok(abs(p$mean - 20) < 1e-6, "a precise observation overrides a vague prior")
ok(p$info_frac > 0.999, "info_frac approaches 1 when data dominate")

# Sample size: n identical cores must beat one, and tighten as sqrt(n).
p1 <- bayes_posterior(10, 1, 20, 1, matrix(1, 1, 1))
p4 <- bayes_posterior(10, 1, rep(20, 4), rep(1, 4), matrix(1, 1, 4))
ok(p4$mean > p1$mean, "four cores pull further from the prior than one")
eq(p4$sd, 1 / sqrt(5), "precision adds across cores (1 prior + 4 cores)")
ok(p4$info_frac > p1$info_frac, "info_frac rises with the number of cores")

# Distance: a closer core must move the posterior further.
near <- bayes_posterior(10, 1, 20, 1, matrix(kernel_gaussian(0.5, 2.5), 1, 1))
far  <- bayes_posterior(10, 1, 20, 1, matrix(kernel_gaussian(8.0, 2.5), 1, 1))
ok(near$mean > far$mean, "a nearer core moves the posterior further")
ok(abs(far$mean - 10) < 0.1, "a distant core barely moves the posterior")

# Posterior must never leave the interval spanned by prior and observations.
set.seed(11)
mu0 <- runif(50, 0, 100); s0 <- runif(50, 0.5, 10)
yv  <- c(5, 60); sy <- c(2, 3)
Wr  <- matrix(runif(100), 50, 2)
pr  <- bayes_posterior(mu0, s0, yv, sy, Wr)
lo <- pmin(mu0, min(yv)); hi <- pmax(mu0, max(yv))
ok(all(pr$mean >= lo - 1e-9 & pr$mean <= hi + 1e-9),
   "posterior mean always lies between the prior and the observations")
ok(all(pr$sd <= s0 + 1e-9), "posterior is never less certain than the prior")
ok(all(pr$info_frac >= 0 & pr$info_frac <= 1), "info_frac stays in [0, 1]")

# A missing prior must not fabricate a value.
p <- bayes_posterior(c(NA_real_, NA_real_), c(NA_real_, NA_real_),
                     20, 1, matrix(c(1, 0), 2, 1))
eq(p$mean[1], 20, "with no prior the posterior rests entirely on the cores")
ok(is.na(p$mean[2]), "no prior and no core weight yields NA, not a number")
ok(is.na(p$shift[1]), "shift is NA when there was no prior to shift from")

throws(bayes_posterior(10, 1, 20, 0, matrix(1, 1, 1)),
       "rejects a zero observation sd")
throws(bayes_posterior(10, 1, 20, 1, matrix(2, 1, 1)),
       "rejects kernel weights outside [0, 1]")

# =============================================================================
section("build_weight_matrix stratum gating")
# =============================================================================
D <- matrix(c(0.1, 0.1), 1, 2)            # one pixel, two very close cores
Wg <- build_weight_matrix(D, 2.5, 10, core_stratum = c("peat", "mineral"),
                          pixel_stratum = "peat")
ok(Wg[1, 1] > 0.9, "a same-stratum core keeps its weight")
eq(Wg[1, 2], 0, "a different-stratum core is gated to zero")

Wu <- build_weight_matrix(D, 2.5, 10)
ok(all(Wu > 0.9), "without gating both cores contribute")

# =============================================================================
section("core_sigma")
# =============================================================================
stock <- c(8.2, 5.2, 2.6, 6.9, 3.4, 3.1, 2.7, 3.5)
strat <- rep(c("mineral", "peat"), c(5, 3))
cs <- core_sigma(stock, strat)
ok(all(cs$sigma_used >= cs$floor_sd - 1e-9),
   "no core is given an sd below the pooled floor")
ok(any(cs$floored), "the small-n peat sd is floored")
ok(all(cs$sigma_used[strat == "mineral"] >= cs$sigma_raw[strat == "mineral"] - 1e-9),
   "the larger mineral sd is kept as-is")

# =============================================================================
section("full-column split and recombine")
# =============================================================================
sp <- split_full_column(100, 0.13)
eq(sp$shallow, 13, "shallow share is taken correctly")
eq(sp$deep, 87, "deep remainder is the complement")
eq(sp$shallow + sp$deep, 100, "the split conserves total carbon")
throws(split_full_column(100, 0), "rejects a zero shallow fraction")
throws(split_full_column(100, 1), "rejects a shallow fraction of one")

rc <- recombine_full_column(15, 3, 87, 4)
eq(rc$mean, 102, "recombination adds the layers")
eq(rc$sd, 5, "uncertainties add in quadrature (3-4-5)")

eq(uniform_shallow_fraction(30, 184), 30 / 184,
   "uniform-density fraction is depth ratio")
ok(uniform_shallow_fraction(30, 20) < 1,
   "fraction is capped below 1 when peat is shallower than the window")

# ---- summary -----------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(sprintf("passed %d   failed %d\n", .pass, .fail))
if (.fail > 0L) {
  cat("failures:\n"); for (f in .failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("all tests passed\n")
