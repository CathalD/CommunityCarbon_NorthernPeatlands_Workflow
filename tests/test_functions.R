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

# =============================================================================
section("oriented AOI geometry")
# =============================================================================
# Points along a known bearing must recover that bearing.
b <- 45
d <- seq(-10, 10, 2) * 1000
lat0 <- 56; lon0 <- -87.7
th <- (90 - b) * pi / 180
xx <- d * cos(th); yy <- d * sin(th)
lo <- lon0 + xx / (111320 * cos(lat0 * pi / 180)); la <- lat0 + yy / 111320
pa <- principal_axis_bearing(lo, la)
ok(abs(pa$bearing_deg - b) < 1, "recovers a known 45 degree bearing")
ok(pa$var_explained > 0.99, "collinear points put all variance on axis 1")

# A square must be built in METRES, not degrees: at 56 N one degree of
# longitude is only 0.56 of a degree of latitude on the ground, so a
# degree-square is a north-south rectangle and would fail for the wrong reason.
sq_m <- 10000
sq_lon <- lon0 + c(-1, 1, -1, 1) * sq_m / (111320 * cos(lat0 * pi / 180))
sq_lat <- lat0 + c(-1, -1, 1, 1) * sq_m / 111320
pa2 <- principal_axis_bearing(sq_lon, sq_lat)
ok(pa2$var_explained < 0.75, "a square arrangement (in metres) has no dominant axis")

rect <- oriented_rectangle(lo, la, b, along_buffer_km = 5, across_buffer_km = 2)
ok(all(point_in_oriented_rect(lo, la, rect)),
   "every input point falls inside the rectangle")
ok(abs(rect$width_km - 4) < 0.2, "width is twice the across-buffer for collinear points")
ok(abs(rect$length_km - 30) < 0.5, "length is the span plus both along-buffers")
ok(length(rect$lon) == 5 && rect$lon[1] == rect$lon[5],
   "polygon ring is closed")
ok(rect$area_km2 > 0, "area is positive")

# A point far outside must be rejected.
ok(!point_in_oriented_rect(-80, 56, rect), "a distant point is outside")
# A point displaced perpendicular to the axis beyond the 2 km across-buffer
# must be outside. Build it by inverting the rotation used in the module:
#   x = u cos(th) - v sin(th),  y = u sin(th) + v cos(th)
across_m <- 5000                     # well beyond the 2 km buffer
fx <- -across_m * sin(rect$theta)
fy <-  across_m * cos(rect$theta)
far <- .from_local(fx, fy, rect$lon0, rect$lat0)
ok(!point_in_oriented_rect(far$lon, far$lat, rect),
   "a point beyond the across-buffer is outside")
# ...and one just inside it must be inside, so the test is not passing trivially.
near_m <- 1000
nx <- -near_m * sin(rect$theta); ny <- near_m * cos(rect$theta)
nr <- .from_local(nx, ny, rect$lon0, rect$lat0)
ok(point_in_oriented_rect(nr$lon, nr$lat, rect),
   "a point within the across-buffer is inside")

oc <- oriented_coordinates(lo, la, rect)
ok(nrow(oc) == length(lo), "coordinates returned for every point")
ok(all(abs(oc$across_km) < 0.5), "collinear points have ~zero across-axis spread")
ok(abs(diff(range(oc$along_km)) - 20) < 0.5, "along-axis span matches the input")

# Envelope must contain the oriented rectangle.
ok(rect$bbox[["xmin"]] <= min(rect$lon) && rect$bbox[["xmax"]] >= max(rect$lon),
   "bbox envelopes the polygon in longitude")
ok(rect$bbox[["ymin"]] <= min(rect$lat) && rect$bbox[["ymax"]] >= max(rect$lat),
   "bbox envelopes the polygon in latitude")

sn <- ee_polygon_snippet(rect, "aoi")
ok(grepl("ee\\.Geometry\\.Polygon", sn), "emits an Earth Engine polygon literal")
ok(lengths(regmatches(sn, gregexpr("\\[-?[0-9]", sn))) >= 5,
   "snippet carries all vertices")

gp <- tempfile(fileext = ".geojson")
write_geojson_polygon(rect$lon, rect$lat, gp, props = list(name = "test"))
gt <- paste(readLines(gp), collapse = "")
ok(grepl('"Polygon"', gt), "writes a GeoJSON polygon")
ok(grepl('"name": "test"', gt), "carries properties")
unlink(gp)

# =============================================================================
section("read_geojson_polygon -- the parser that bit us")
# =============================================================================
# The failure this guards against: an earlier reader scraped every decimal
# number out of the GeoJSON, including PROPERTY values (bearing 142.91, area
# 2758.4, buffers ...), and paired them as coordinates. It produced a ring that
# looked plausible and was wrong. So the fixture below deliberately carries
# properties whose values would corrupt a naive parse.
gp2 <- tempfile(fileext = ".geojson")
write_geojson_polygon(rect$lon, rect$lat, gp2,
                      props = list(name = "test", bearing_deg = 142.91,
                                   area_km2 = 2758.4, width_km = 33.0))
ring <- read_geojson_polygon(gp2)
ok(nrow(ring) == 5, "reads exactly the ring vertices, not the properties")
# The writer emits %.8f, i.e. ~1 mm on the ground at this latitude, so the
# round-trip tolerance is set to the write precision rather than to zero.
eq(ring$lon, rect$lon, "longitudes round-trip to write precision", tol = 1e-7)
eq(ring$lat, rect$lat, "latitudes round-trip to write precision", tol = 1e-7)
ok(max(abs(ring$lon - rect$lon)) * 111320 * cos(56 * pi / 180) < 0.01,
   "longitude round-trip error is under a centimetre on the ground")
ok(ring$lon[1] == ring$lon[nrow(ring)] && ring$lat[1] == ring$lat[nrow(ring)],
   "returned ring is closed")

# The naive approach, reproduced, must disagree -- proving the fixture is a
# real trap and not a vacuous test.
naive_txt <- paste(readLines(gp2, warn = FALSE), collapse = " ")
naive_nums <- as.numeric(regmatches(naive_txt,
                gregexpr("-?[0-9]+\\.[0-9]+", naive_txt))[[1]])
ok(length(naive_nums) > 2 * nrow(ring),
   "a whole-file number scrape picks up more numbers than there are coordinates")

ok(all(abs(ring$lon) <= 180) && all(abs(ring$lat) <= 90),
   "parsed coordinates are geographically possible")
unlink(gp2)

# Malformed input must fail loudly, not return something usable.
bad <- tempfile(fileext = ".geojson")
writeLines('{"type":"FeatureCollection","features":[]}', bad)
throws(read_geojson_polygon(bad), "STOPS when there is no coordinates key")
writeLines('{"coordinates": [[[0,0],[1,1]]]}', bad)
throws(read_geojson_polygon(bad), "STOPS when the ring has too few vertices")
writeLines('{"coordinates": [[[0,0],[999,0],[999,999],[0,999],[0,0]]]}', bad)
throws(read_geojson_polygon(bad), "STOPS on impossible coordinates")
unlink(bad)

# ring_to_ee_coords shape.
ec <- ring_to_ee_coords(data.frame(lon = c(1, 2, 3), lat = c(4, 5, 6)))
ok(is.list(ec) && length(ec) == 3 && identical(ec[[1]], c(1, 4)),
   "ring converts to the nested list Earth Engine expects")

# =============================================================================
section("Earth Engine download sizing and tiling")
# =============================================================================
# The real AOI envelope, which is what getDownloadURL actually rasterises.
aoi_bb <- c(xmin = -88.317, ymin = 55.6799, xmax = -87.0821, ymax = 56.4576)

# Grid dimensions must reflect BOTH traps that made the first estimate wrong.
d <- gee_grid_dims(aoi_bb, 30, geographic = TRUE)
ok(abs(d$nx - 4583) < 5 && abs(d$ny - 2887) < 5,
   "geographic grid dims match Earth Engine's equatorial degree step")
ok(d$nx / d$ny > 1.5,
   "longitude is oversampled at 56 N, as EPSG:4326 + scale forces")
dp <- gee_grid_dims(c(xmin = 0, ymin = 0, xmax = 3000, ymax = 3000), 30,
                    geographic = FALSE)
ok(dp$nx == 100 && dp$ny == 100, "projected grid dims are a plain division")

# The estimate must now REACH the size the server reported, not undercut it.
# Server refused one Float32 band over this envelope at 66,155,605 bytes.
e1 <- gee_estimate_bytes(aoi_bb, 30, 1)
ok(e1 > 50e6, "one Float32 band over the envelope is correctly seen as too large")
ok(abs(e1 - 66155605) / 66155605 < 0.35,
   "estimate is within 35% of the size Earth Engine actually reported")
ok(gee_estimate_bytes(aoi_bb, 30, 19) > gee_estimate_bytes(aoi_bb, 30, 1),
   "more bands means more bytes")
ok(gee_estimate_bytes(aoi_bb, 250, 1) < gee_estimate_bytes(aoi_bb, 30, 1),
   "a coarser scale means fewer bytes")

# TILING is the fix band-chunking could not provide.
tl <- gee_tile_bbox(aoi_bb, 30, 1)
ok(!is.null(tl) && length(tl) > 1L,
   "a single band too large for one request IS tiled rather than refused")
ok(all(vapply(tl, function(b)
       gee_estimate_bytes(b, 30, 1) <= 24 * 1024^2, logical(1))),
   "every tile fits under the request ceiling")
# Tiles must cover the box exactly, with no gaps or overlap.
ok(abs(min(vapply(tl, `[[`, numeric(1), "xmin")) - aoi_bb[["xmin"]]) < 1e-9 &&
   abs(max(vapply(tl, `[[`, numeric(1), "xmax")) - aoi_bb[["xmax"]]) < 1e-9,
   "tiles span the full longitude range")
ok(abs(min(vapply(tl, `[[`, numeric(1), "ymin")) - aoi_bb[["ymin"]]) < 1e-9 &&
   abs(max(vapply(tl, `[[`, numeric(1), "ymax")) - aoi_bb[["ymax"]]) < 1e-9,
   "tiles span the full latitude range")
areas <- sum(vapply(tl, function(b)
  (b[["xmax"]] - b[["xmin"]]) * (b[["ymax"]] - b[["ymin"]]), numeric(1)))
full <- (aoi_bb[["xmax"]] - aoi_bb[["xmin"]]) * (aoi_bb[["ymax"]] - aoi_bb[["ymin"]])
ok(abs(areas - full) / full < 1e-9, "tiles partition the box exactly")

# Something already small enough must NOT be tiled.
ok(length(gee_tile_bbox(aoi_bb, 250, 1)) == 1L,
   "a request that already fits is left as one tile")
# More bands needs more tiles.
ok(length(gee_tile_bbox(aoi_bb, 30, 19)) >= length(gee_tile_bbox(aoi_bb, 30, 1)),
   "a 19-band stack needs at least as many tiles as one band")
# An impossible request returns NULL rather than an absurd tile count.
ok(is.null(gee_tile_bbox(aoi_bb, 1, 64, max_tiles = 64L)),
   "refuses rather than emitting an unreasonable number of tiles")

# Band chunking still applies WITHIN a tile.
tb <- tl[[1]]
ch <- gee_band_chunks(paste0("b", 1:19), tb, 30)
ok(!is.null(ch), "bands can be chunked within a tile")
ok(identical(sort(unlist(ch, use.names = FALSE)), sort(paste0("b", 1:19))),
   "chunking preserves every band exactly once")
ok(all(vapply(ch, function(g)
       gee_estimate_bytes(tb, 30, length(g)) <= 24 * 1024^2, logical(1))),
   "every band chunk within a tile fits")

# =============================================================================
section("Lambert Conformal Conic (EPSG:3978)")
# =============================================================================
p78 <- lcc_params(3978)
ok(p78$epsg == 3978L && p78$datum == "NAD83", "3978 is NAD83")
ok(lcc_params(3979)$datum == "NAD83(CSRS)", "3979 is the CSRS realisation")
ok(identical(p78$lat_1, lcc_params(3979)$lat_1) &&
   identical(p78$lon_0, lcc_params(3979)$lon_0),
   "3978 and 3979 share identical projection parameters")
throws(lcc_params(4326), "rejects an EPSG code it does not know")

# At the projection origin, easting must be exactly the false easting.
o <- lonlat_to_lcc(p78$lon_0, p78$lat_0, p78)
eq(o$x, p78$x_0, "origin longitude maps to the false easting", tol = 1e-6)
eq(o$y, p78$y_0, "origin latitude maps to the false northing", tol = 1e-6)

# The central meridian must map to x = 0 at every latitude.
cm <- lonlat_to_lcc(rep(p78$lon_0, 4), c(50, 60, 70, 77), p78)
ok(all(abs(cm$x) < 1e-6), "the central meridian is x = 0 at all latitudes")

# Symmetry about the central meridian.
l <- lonlat_to_lcc(-100, 56, p78); r <- lonlat_to_lcc(-90, 56, p78)
ok(abs(l$x + r$x) < 1e-6, "points equidistant from the central meridian mirror in x")
eq(l$y, r$y, "and share the same northing", tol = 1e-6)

# Standard parallels: the LCC scale factor is exactly 1 there, so a degree of
# longitude must project to its true length along that parallel.
#
# The reference has to be the ELLIPSOIDAL parallel arc. Using haversine_km()
# here would be wrong and would fail for the wrong reason: it assumes a sphere,
# and at 49 N the spherical and ellipsoidal parallel radii differ by 0.52% --
# far more than any error worth detecting. The formula below is the standard
# geodetic one and uses only the ellipsoid constants, not the projection
# machinery, so it is not circular.
parallel_arc_m <- function(lat_deg, dlon_deg, a, f) {
  e2  <- 2 * f - f^2
  phi <- lat_deg * pi / 180
  N   <- a / sqrt(1 - e2 * sin(phi)^2)       # radius of curvature in the prime vertical
  N * cos(phi) * (dlon_deg * pi / 180)
}
for (sp in c(p78$lat_1, p78$lat_2)) {
  a1 <- lonlat_to_lcc(p78$lon_0, sp, p78)
  b1 <- lonlat_to_lcc(p78$lon_0 + 1, sp, p78)
  proj_m <- sqrt((b1$x - a1$x)^2 + (b1$y - a1$y)^2)
  true_m <- parallel_arc_m(sp, 1, p78$a, p78$f)
  ok(abs(proj_m - true_m) / true_m < 2e-4,
     sprintf("scale is true at the standard parallel %g N", sp))
}
# ...and NOT true away from them: LCC compresses between the parallels.
mid <- (p78$lat_1 + p78$lat_2) / 2
am <- lonlat_to_lcc(p78$lon_0, mid, p78); bm <- lonlat_to_lcc(p78$lon_0 + 1, mid, p78)
ok(sqrt((bm$x - am$x)^2 + (bm$y - am$y)^2) <
     parallel_arc_m(mid, 1, p78$a, p78$f),
   "between the standard parallels the projection compresses, as LCC must")

# Round trip must return the original coordinates, not merely be self-consistent.
set.seed(3)
tl <- runif(50, -110, -80); ta <- runif(50, 48, 62)
fw <- lonlat_to_lcc(tl, ta, p78)
bk <- lcc_to_lonlat(fw$x, fw$y, p78)
ok(max(abs(bk$lon - tl)) < 1e-9, "round trip recovers longitude")
ok(max(abs(bk$lat - ta)) < 1e-9, "round trip recovers latitude")

# The Fort Severn cores, against values verified independently with PROJ.
fs <- lonlat_to_lcc(-87.62820, 55.9659990, p78)
ok(abs(fs$x - 449599.9948) < 0.01 && abs(fs$y - 792274.1280) < 0.01,
   "FS-01 matches the coordinate PROJ produces, to the centimetre")

throws(lonlat_to_lcc(-95, 90, p78), "rejects a pole, where LCC is undefined")

cmp <- compare_projections(c(0, 100), c(0, 100), c(0, 103), c(0, 104),
                           id = c("a", "b"))
eq(cmp$offset_m, c(0, 5), "compare_projections measures the offset (3-4-5)")

# =============================================================================
section("isTRUE_vec")
# =============================================================================
ok(identical(isTRUE_vec(c(TRUE, FALSE, NA)), c(TRUE, FALSE, FALSE)),
   "logical vector with NA treated as FALSE")
ok(identical(isTRUE_vec(c("TRUE", "true", "False", "")),
             c(TRUE, TRUE, FALSE, FALSE)),
   "character TRUE/true recognised, others FALSE")
ok(length(isTRUE_vec(c(TRUE, TRUE))) == 2,
   "returns one value per element, unlike isTRUE()")

# =============================================================================
section("r2_cv and skill_vs_null -- missing predictions must not create skill")
# =============================================================================
local({
  obs <- c(1, 2, 3, 4, 5, 6)
  # A model that predicted only ONE core, and got it exactly right, must not be
  # rewarded for the five folds it never ran. Before the fix, the error sum was
  # taken over the predicted row while the total sum was taken over all six,
  # which handed a mean-only null r2 = +0.95 on the 0-30 cm window.
  pred_one <- c(NA, NA, NA, NA, NA, 6)
  ok(is.na(r2_cv(obs, pred_one)),
     "r2 from a single surviving fold is NA, not near 1")

  pred_two <- c(NA, NA, NA, NA, 5, 6)
  # On the two rows it did predict, the model is perfect, so r2 on those rows
  # is 1. The point is that the denominator is now built from the same rows.
  eq(r2_cv(obs, pred_two), 1, "r2 is computed over the predicted rows only")

  # The mean-only null, scored on a subset, must come out at 0 -- not positive.
  pred_null <- c(NA, NA, NA, mean(obs), mean(obs), mean(obs))
  sub <- obs[4:6]
  eq(r2_cv(obs, pred_null), 1 - sum((sub - mean(obs))^2) / sum((sub - mean(sub))^2),
     "a null scored on a subset uses that subset's own variance")

  ok(is.na(skill_vs_null(obs, rep(NA_real_, 6))),
     "skill against the null is NA when nothing was predicted")

  m <- metric_set(obs, pred_two)
  eq(m$n_predicted, 2, "metric_set counts predicted rows")
  eq(m$n_obs, 6, "metric_set counts observed rows")
  eq(m$coverage, 2 / 6, "metric_set reports coverage so partial CV is visible")

  full <- metric_set(obs, obs)
  eq(full$coverage, 1, "complete cross-validation reports coverage 1")
  eq(full$r2_cv, 1, "a perfect complete prediction still scores r2 = 1")
})

# =============================================================================
section("random forest (backend-agnostic)")
# =============================================================================
local({
  set.seed(11)
  n <- 120
  X <- data.frame(sig = runif(n), noise = runif(n), flat = rep(2, n))
  y <- 4 * X$sig + rnorm(n, sd = 0.2)

  m <- fit_rf(X, y, num_trees = 60L, seed = 3L)
  ok(inherits(m, "ccnp_rf"), "fit_rf returns a tagged object")
  ok(!("flat" %in% m$vars), "a zero-variance covariate is dropped, not fitted")

  tr <- 1:90; te <- 91:120
  mt <- fit_rf(X[tr, ], y[tr], num_trees = 100L, seed = 3L)
  ok(r2_cv(y[te], predict_rf(mt, X[te, ])) > 0.7,
     "recovers a strong signal on held-out rows")

  # Reproducibility is not optional: an unseeded forest hands the council a
  # different map on every run.
  a <- predict_rf(fit_rf(X, y, num_trees = 40L, seed = 9L), X)
  b <- predict_rf(fit_rf(X, y, num_trees = 40L, seed = 9L), X)
  ok(identical(a, b), "the same seed gives the identical forest")

  pt <- predict_rf(m, X[1:5, ], per_tree = TRUE)
  eq(nrow(pt$trees), 5, "per-tree predictions have one row per point")
  eq(length(pt$pred), 5, "and an ensemble mean per point")
  eq(pt$pred, rowMeans(pt$trees), "the ensemble mean IS the mean of the trees")
  ok(all(rf_predict_sd(m, X[1:5, ]) >= 0), "ensemble spread is non-negative")

  # A forest handed nothing usable must return the mean, not fail or invent.
  dg <- fit_rf(data.frame(z = rep(1, 5)), c(1, 2, 3, 4, 5), seed = 1L)
  eq(predict_rf(dg, data.frame(z = rep(1, 2))), c(3, 3),
     "with no usable covariate the forest degrades to the mean")

  throws(predict_rf(m, X[, "noise", drop = FALSE]),
         "refuses to predict when a fitted covariate is absent")
})

# =============================================================================
section("rf_importance -- ranking and permutation null")
# =============================================================================
local({
  set.seed(21)
  n <- 8
  X <- data.frame(sig = c(1, 2, 3, 4, 5, 6, 7, 8) / 8,
                  noise = runif(n), noise2 = runif(n))
  y <- 6 * X$sig + rnorm(n, sd = 0.2)
  imp <- rf_importance(X, y, paste0("c", 1:n), n_perm = 2L, n_null = 5L,
                       seed = 4L, fit_args = list(num_trees = 40L),
                       verbose = FALSE)
  ok(imp$covariate[1] == "sig", "the real signal ranks first")
  ok(imp$importance[1] > imp$importance[nrow(imp)],
     "importance is ordered descending")
  ok(all(imp$p_value >= 1 / (1 + 5)),
     "no p-value below the floor the null draw count allows")
  ok(all(imp$rank == seq_len(nrow(imp))), "rank is assigned after ordering")

  # With no null draws affordable there is no evidence, so no p-value. Reporting
  # 1 or 0 would be a claim about a test that was never run.
  imp0 <- rf_importance(X, y, paste0("c", 1:n), n_perm = 2L, n_null = 0L,
                        seed = 4L, fit_args = list(num_trees = 20L),
                        verbose = FALSE)
  ok(all(is.na(imp0$p_value)), "no null draws means NA p-values, not 1")
  ok(!any(imp0$beats_chance), "and nothing is marked as beating chance")
})

# =============================================================================
section("cv_spatial_buffer")
# =============================================================================
local({
  # Two tight clusters ~40 km apart, mimicking the peat/mineral design.
  lon <- c(-87.7, -87.701, -87.702, -87.0, -87.001, -87.002)
  lat <- c(56.00, 56.001, 56.002, 56.30, 56.301, 56.302)
  X <- data.frame(v = c(1, 1.1, 0.9, 5, 5.1, 4.9))
  y <- c(1, 1.1, 0.9, 5, 5.1, 4.9)
  g <- paste0("p", 1:6)

  # A buffer of zero is ordinary leave-one-out: nothing extra is removed.
  z <- cv_spatial_buffer(X, y, g, lon, lat, 0, fit_null, predict_null)
  eq(sum(z$folds$n_excluded_by_buffer), 0, "a zero buffer excludes nothing")
  eq(z$n_folds_used, 6, "and every fold runs")

  # A 1 km buffer removes each point's own cluster-mates, so each fold is
  # trained only on the OTHER cluster -- the honest, and much harder, test.
  b <- cv_spatial_buffer(X, y, g, lon, lat, 1, fit_null, predict_null)
  eq(unique(b$folds$n_excluded_by_buffer), 2,
     "a 1 km buffer removes the two cluster-mates of each held-out point")
  eq(unique(b$folds$n_train), 3, "leaving only the far cluster to train on")
  ok(b$metrics$r2_cv < 0,
     "predicting one cluster from the other has no skill, as it must not")

  # A buffer wide enough to empty the training set must SKIP folds and say so,
  # never quietly score whichever folds happened to survive.
  w <- cv_spatial_buffer(X, y, g, lon, lat, 500, fit_null, predict_null)
  eq(w$n_folds_used, 0, "an all-consuming buffer runs no folds")
  eq(w$n_folds_skipped, 6, "and reports every skip")
  ok(all(is.na(w$predictions$pred)), "producing no predictions at all")
  ok(is.na(w$metrics$r2_cv), "and refusing to report an r2")
})

# =============================================================================
section("warp to EPSG:3978")
# =============================================================================
local({
  bbox <- c(xmin = -88.4, xmax = -87.0, ymin = 55.8, ymax = 56.4)
  src  <- make_grid(bbox, 300)
  # A surface that is exactly linear in lon/lat is reproduced exactly by
  # bilinear interpolation, so any error here is the projection's, not the
  # interpolator's.
  z <- outer(src$y, src$x, function(la, lo) 10 * lo + 3 * la)

  ring_lon <- c(-88.2, -87.2, -87.2, -88.2)
  ring_lat <- c(55.90, 55.90, 56.30, 56.30)
  g <- lcc_target_grid(ring_lon, ring_lat, 300)

  eq(g$epsg, 3978, "target grid carries the deliverable EPSG")
  ok(g$xmin %% 300 == 0 && g$ymax %% 300 == 0,
     "grid origin snaps to whole cells so stacked rasters align")
  ok(g$nx > 0 && g$ny > 0, "grid has positive extent")

  w <- warp_to_lcc(z, src, g)
  gx <- rep(g$x, times = g$ny); gy <- rep(g$y, each = g$nx)
  ll <- lcc_to_lonlat(gx, gy)
  expect <- 10 * ll$lon + 3 * ll$lat
  got <- as.vector(t(w))
  fin <- is.finite(got) & is.finite(expect)
  ok(sum(fin) > 1000, "the warp fills the target grid")
  ok(max(abs(got[fin] - expect[fin])) < 1e-6,
     "warped values match the analytic surface at the inverse-projected point")

  # The grid must be built for the CRS it is warped into.
  throws(warp_to_lcc(z, src, g, lcc_params(3979)),
         "refuses to warp a 3978 grid using 3979 parameters")

  # Nearest neighbour must return values that EXIST in the source. Averaging
  # class codes would invent classes nothing downstream would question.
  zc <- matrix(sample(c(180, 181, 188), length(z), TRUE), nrow = nrow(z))
  wn <- warp_to_lcc(zc, src, g, method = "nearest")
  ok(all(wn[is.finite(wn)] %in% c(180, 181, 188)),
     "nearest-neighbour warping never invents a class code")
  wb <- warp_to_lcc(zc, src, g, method = "bilinear")
  ok(!all(wb[is.finite(wb)] %in% c(180, 181, 188)),
     "whereas bilinear would -- which is why class layers must use nearest")
})

# =============================================================================
section("point_in_ring and masking")
# =============================================================================
local({
  sq_x <- c(0, 10, 10, 0); sq_y <- c(0, 0, 10, 10)
  ok(point_in_ring(5, 5, sq_x, sq_y), "centre is inside")
  ok(!point_in_ring(-1, 5, sq_x, sq_y), "west of the square is outside")
  ok(!point_in_ring(11, 5, sq_x, sq_y), "east of the square is outside")
  ok(!point_in_ring(5, 11, sq_x, sq_y), "north of the square is outside")
  eq(sum(point_in_ring(c(5, 5, 20), c(5, 2, 20), sq_x, sq_y)), 2,
     "vectorised over points")

  # A closed ring (first vertex repeated) must give the same answer as an open
  # one; a duplicated vertex would otherwise be counted twice by the ray test.
  ok(identical(point_in_ring(5, 5, c(sq_x, 0), c(sq_y, 0)),
               point_in_ring(5, 5, sq_x, sq_y)),
     "a repeated closing vertex does not change the result")

  ring_lon <- c(-88.2, -87.2, -87.2, -88.2)
  ring_lat <- c(55.90, 55.90, 56.30, 56.30)
  g <- lcc_target_grid(ring_lon, ring_lat, 1000)
  m <- matrix(1, nrow = g$ny, ncol = g$nx)
  mm <- mask_grid_to_ring(m, g, ring_lon, ring_lat)
  ok(sum(is.finite(mm)) < length(mm),
     "masking removes the corners the rotated AOI does not cover")
  ok(sum(is.finite(mm)) > 0.5 * length(mm),
     "but keeps the bulk of the AOI")
})

# =============================================================================
section("densify_ring")
# =============================================================================
local({
  r <- densify_ring(c(0, 1, 1, 0), c(0, 0, 1, 1), n_per_edge = 10L)
  eq(nrow(r), 4 * 9 + 1, "each edge contributes n-1 points, plus the closure")
  eq(c(r$lon[1], r$lat[1]), c(r$lon[nrow(r)], r$lat[nrow(r)]),
     "the densified ring is closed")
  ok(all(r$lon >= -1e-12 & r$lon <= 1 + 1e-12),
     "densified points stay on the ring, never outside it")
  # An already-closed ring must not gain a duplicate edge.
  r2 <- densify_ring(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0), n_per_edge = 10L)
  eq(nrow(r2), nrow(r), "closing the input explicitly changes nothing")
})

# =============================================================================
section("projected GeoTIFF (EPSG:3978)")
# =============================================================================
local({
  ok(!.epsg_is_geographic(3978), "3978 is recognised as projected")
  ok(.epsg_is_geographic(4326), "4326 is recognised as geographic")
  throws(.epsg_is_geographic(31370),
         "an unknown EPSG is refused rather than guessed at")

  kp <- .geo_key_directory(3978, geographic = FALSE, citation_len = 30L)
  ok(kp[1 + 4 * 1 + 3] == 1L || 1024L %in% kp,
     "projected keys include GTModelType")
  ok(3072L %in% kp, "projected keys carry ProjectedCSTypeGeoKey")
  ok(!(2048L %in% kp), "and NOT GeographicTypeGeoKey")
  kg <- .geo_key_directory(4326, geographic = TRUE, citation_len = 7L)
  ok(2048L %in% kg && !(3072L %in% kg), "geographic keys are the mirror image")

  # The citation length must follow the string. A fixed count would make a
  # reader run off the end of a longer CRS name.
  long <- .geo_key_directory(3978, geographic = FALSE, citation_len = 30L)
  ok(30L %in% long, "citation length is carried into the key, not hardcoded")

  p <- tempfile(fileext = ".tif")
  z <- matrix(as.numeric(1:12), nrow = 3, ncol = 4, byrow = TRUE)
  write_geotiff(z, p, xmin = 400000, ymax = 800000, xres = 500, yres = 500,
                epsg = 3978L)
  ok(file.exists(p) && file.size(p) > 0, "a projected GeoTIFF is written")
  # Read the tag set back to confirm the projected geokeys really landed.
  con <- file(p, "rb"); on.exit(close(con), add = TRUE)
  hdr <- readBin(con, "raw", 8)
  off <- readBin(hdr[5:8], "integer", size = 4L, endian = "little")
  seek(con, off)
  n_ent <- readBin(con, "integer", size = 2L, endian = "little", signed = FALSE)
  tags <- integer(n_ent)
  for (i in seq_len(n_ent)) {
    tags[i] <- readBin(con, "integer", size = 2L, endian = "little", signed = FALSE)
    readBin(con, "raw", 10)
  }
  ok(34735L %in% tags, "the GeoKeyDirectory tag is present")
  ok(33550L %in% tags && 33922L %in% tags, "pixel scale and tiepoint are present")
  unlink(p)
})

# =============================================================================
section("ringfence_categorical -- surviving a layer rename")
# =============================================================================
local({
  # The real failure this guards against: the Earth Engine export was renamed
  # gwl_fcs30 -> gwl_class between runs, a literal exclusion list stopped
  # matching, and a wetland CLASS CODE came third in 06's covariate ranking.
  nms <- c("elevation", "s1_vv", "gwl_class", "gwl_wetland", "worldcover",
           "tpi_2km")
  rf <- ringfence_categorical(nms)
  ok(all(c("gwl_class", "gwl_wetland", "worldcover") %in% rf$categorical),
     "class-code layers are ring-fenced")
  ok(all(c("elevation", "s1_vv", "tpi_2km") %in% rf$numeric),
     "genuine covariates are left alone")
  ok(!any(c("gwl_class", "worldcover") %in% rf$numeric),
     "and never appear on both sides")

  # The old name must still be caught, so re-running an older export is safe.
  ok("gwl_fcs30" %in% ringfence_categorical(c("gwl_fcs30", "ndvi"))$categorical,
     "the pre-rename name is caught too")

  # A layer that LOOKS categorical but matches no pattern must be surfaced,
  # not silently passed through as a predictor.
  sus <- ringfence_categorical(c("elevation", "soil_type_code"))
  ok("soil_type_code" %in% sus$suspected,
     "an unmatched but suspicious name is reported for the caller to act on")
  ok(!length(ringfence_categorical(c("elevation", "ndvi"))$suspected),
     "ordinary covariate names raise no suspicion")
})

# =============================================================================
section("pick_prior_column -- preference order and honest failure")
# =============================================================================
local({
  d <- data.frame(sothe2022_soc_0_30_kgm2 = c(1, 2, 3, 4),
                  soilgrids_ocs_0_30_kgm2 = c(5, 6, 7, 8),
                  li2025_peat_carbon_kgm2 = c(9, 9, 9, 9))
  p <- pick_prior_column(d, c("sothe.*0_30", "soilgrids.*0_30"))
  eq(p$column == "sothe2022_soc_0_30_kgm2", TRUE,
     "the first preference wins when it is usable")

  # Too few finite values: fall through to the next preference rather than
  # returning a column that cannot support a comparison.
  d2 <- d; d2$sothe2022_soc_0_30_kgm2 <- c(1, NA, NA, NA)
  p2 <- pick_prior_column(d2, c("sothe.*0_30", "soilgrids.*0_30"))
  ok(p2$column == "soilgrids_ocs_0_30_kgm2",
     "a mostly-missing first preference falls through to the second")

  # Li is a full peat column and must never be selected as a 0-30 cm anchor.
  p3 <- pick_prior_column(d, c("sothe.*0_30"))
  ok(!identical(p3$column, "li2025_peat_carbon_kgm2"),
     "a layer outside the preference list is never selected")

  p4 <- pick_prior_column(d, c("nothing_matches_this"))
  ok(is.na(p4$column), "no match returns NA rather than an arbitrary column")
  ok(length(p4$considered) == 3,
     "and reports what WAS present, so the preference list can be fixed")
})

# ---- summary -----------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(sprintf("passed %d   failed %d\n", .pass, .fail))
if (.fail > 0L) {
  cat("failures:\n"); for (f in .failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("all tests passed\n")
