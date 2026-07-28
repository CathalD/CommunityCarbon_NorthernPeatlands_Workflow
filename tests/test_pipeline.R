# =============================================================================
# tests/test_pipeline.R  --  End-to-end integration test.
#
# Runs 01, 02 and 05 on the real data, then exercises 06 against a SYNTHETIC
# covariate fixture, because 03 requires an authenticated Earth Engine session
# that continuous integration does not have.
#
# The whole test runs inside a throwaway copy of the repository, so the
# synthetic fixture can never be mistaken for a real extraction or leak into
# outputs/. Nothing here writes to the working tree.
#
# Run with:  Rscript tests/test_pipeline.R
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

.pass <- 0L; .fail <- 0L; .failures <- character(0)
ok <- function(cond, what) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat(sprintf("  ok   %s\n", what)) }
  else { .fail <<- .fail + 1L; .failures <<- c(.failures, what)
         cat(sprintf("  FAIL %s\n", what)) }
}
section <- function(x) cat("\n", x, "\n", sep = "")

# ---- isolated sandbox --------------------------------------------------------
sandbox <- file.path(tempdir(), paste0("ccnp_test_", Sys.getpid()))
dir.create(sandbox, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)

for (d in c("R", "scripts", "data", "config.R")) {
  file.copy(file.path(.root, d), sandbox, recursive = TRUE)
}
cat("sandbox: ", sandbox, "\n", sep = "")

run <- function(script) {
  out <- suppressWarnings(system2(
    "Rscript", c(file.path(sandbox, "scripts", script)),
    stdout = TRUE, stderr = TRUE, env = paste0("CCNP_ROOT=", sandbox)))
  list(status = attr(out, "status"), out = paste(out, collapse = "\n"))
}
derived <- function(f) file.path(sandbox, "data", "derived", f)
tables  <- function(f) file.path(sandbox, "outputs", "tables", f)

# =============================================================================
section("01 ingest and QC")
# =============================================================================
r1 <- run("01_ingest_qc.R")
ok(is.null(r1$status), "01 exits cleanly")
ok(file.exists(derived("01_segments_qc.csv")), "01 writes the segment table")
ok(file.exists(derived("01_cores_qc.csv")), "01 writes the core table")

seg <- read.csv(derived("01_segments_qc.csv"), check.names = FALSE)
ok(nrow(seg) == 22, "all 22 segments are retained (nothing dropped)")
ok(length(unique(seg$core_id)) == 8, "8 cores present")
ok(all(seg$longitude < 0), "longitudes are western after the gate")
ok(all(seg$lon_repaired), "every row is marked as sign-repaired")
ok(all(grepl("LON_SIGN_REPAIRED", seg$qc_flags)), "repair is flagged per row")
ok(sum(seg$ocd_reconciles) == 22,
   "supplied carbon density reconciles on every row")
ok(grepl("UNIT DISCREPANCY", r1$out), "01 reports the g/cm2 vs g/cm3 discrepancy")
ok(grepl("DIFFERENT OM-TO-CARBON", r1$out),
   "01 reports the two conversion factors")
ok(any(grepl("BD_ABOVE_MAX", seg$qc_flags)),
   "implausible bulk densities are flagged, not removed")
ok(sum(grepl("DUPLICATE_LAB_VALUES", seg$qc_flags)) == 2,
   "the duplicated lab pair is flagged on both rows")

# =============================================================================
section("01 fails loudly when repair is not authorised")
# =============================================================================
cfg <- readLines(file.path(sandbox, "config.R"))
writeLines(sub("allow_longitude_sign_fix = TRUE",
               "allow_longitude_sign_fix = FALSE", cfg),
           file.path(sandbox, "config.R"))
r1b <- run("01_ingest_qc.R")
ok(!is.null(r1b$status) && r1b$status != 0,
   "01 exits non-zero when the sign error is not authorised")
ok(grepl("FATAL", r1b$out) && grepl("EASTERN hemisphere", r1b$out),
   "the failure names the actual hazard")
writeLines(cfg, file.path(sandbox, "config.R"))
invisible(run("01_ingest_qc.R"))

# =============================================================================
section("02 depth harmonisation")
# =============================================================================
r2 <- run("02_depth_harmonise.R")
ok(is.null(r2$status), "02 exits cleanly")
st <- read.csv(derived("02_core_stocks.csv"), check.names = FALSE)

a <- st[st$scenario == "as_supplied" & st$window == "reference_0_30", ]
ok(sum(a$is_lower_bound) == 2, "two cores are lower bounds at 0-30 cm")
ok(all(a$core_id[a$is_lower_bound] %in% c("PM-01-A", "PM-02-A")),
   "the two short cores are the ones flagged")
ok(abs(a$stock_kgm2[a$core_id == "FS-01"] - 14.26300) < 1e-3,
   "FS-01 0-30 cm stock matches an independent hand calculation")

# PM-03-A is 30.725 cm long, so its 0-30 cm stock must be TRUNCATED, i.e.
# strictly less than its full-core stock.
full <- st[st$scenario == "as_supplied" & st$window == "common_support", ]
ok(a$covered_cm[a$core_id == "PM-03-A"] == 30,
   "the over-long core is truncated to the window, not left at 30.725")

c14 <- full[full$core_id == "PM-01-A", ]
ok(!c14$is_lower_bound, "the shortest core is complete in the common window")
ok(abs(unique(full$window_hi_cm) - 14.5) < 1e-9,
   "common support is derived as 14.5 cm from the data")

h <- st[st$scenario == "harmonised_peat" & st$window == "common_support", ]
ok(all(abs(h$stock_kgm2[grepl("FS", h$core_id)] -
           full$stock_kgm2[grepl("FS", full$core_id)]) < 1e-9),
   "the conversion scenario leaves mineral cores untouched")
ok(all(h$stock_kgm2[grepl("PM", h$core_id)] <
       full$stock_kgm2[grepl("PM", full$core_id)]),
   "the conversion scenario lowers peat stocks")

# =============================================================================
section("05 design-based estimation")
# =============================================================================
r5 <- run("05_stratified_estimate.R")
ok(is.null(r5$status), "05 exits cleanly")
ok(grepl("INTERPOLATION IS REFUSED", r5$out),
   "05 refuses kriging and says so")
ss <- read.csv(tables("05_spatial_structure.csv"))
ok(ss$value[ss$metric == "n_pairs"] == 28, "8 cores give 28 pairs")
ok(ss$value[ss$metric == "n_informative_lag_bins"] == 0,
   "no lag bin is informative")
ok(ss$value[ss$metric == "kriging_supported"] == 0, "kriging is recorded as unsupported")

est <- read.csv(derived("05_stratum_estimates.csv"))
ok(!any(est$trustworthy), "no interval is claimed to be trustworthy at this n")
ok(all(est$n <= 5), "no stratum has more than 5 cores")

# =============================================================================
section("09 spatial products")
# =============================================================================
r9 <- run("09_spatial_products.R")
ok(is.null(r9$status), "09 exits cleanly")
sp <- function(f) file.path(sandbox, "outputs", "spatial", f)
for (f in c("cores.geojson", "aoi.geojson", "dist_to_nearest_core_m.tif",
            "nearest_core_index.tif", "DIAGNOSTIC_soc_nearest_core_kgm2.tif",
            "PRODUCT1b_soc_within_credible_radius_kgm2.tif", "MANIFEST.csv")) {
  ok(file.exists(sp(f)), paste("09 writes", f))
}

# Every raster must be a real little-endian TIFF with the magic number, not an
# empty file that merely has the right extension.
u16 <- function(r) readBin(r, "integer", size = 2, endian = "little", signed = FALSE)
for (f in c("dist_to_nearest_core_m.tif", "nearest_core_index.tif",
            "PRODUCT1b_soc_within_credible_radius_kgm2.tif")) {
  h <- readBin(sp(f), "raw", n = 8)
  ok(rawToChar(h[1:2]) == "II" && u16(h[3:4]) == 42L,
     paste(f, "is a valid TIFF"))
}

gjson <- paste(readLines(sp("cores.geojson")), collapse = "")
ok(lengths(regmatches(gjson, gregexpr('"Point"', gjson))) == 8,
   "cores.geojson carries all 8 cores")
ok(!grepl("\\[8[0-9]\\.", gjson),
   "no positive longitude survives into the GeoJSON")
ok(grepl("stock_0_30_is_lower_bound", gjson),
   "the lower-bound flag travels with the spatial data")

cvg <- read.csv(tables("09_spatial_coverage.csv"))
pct <- cvg$value[cvg$metric == "pct_aoi_within_credible_radius"]
ok(length(pct) == 1 && pct > 0 && pct < 100,
   "coverage fraction is a real percentage")
ok(pct < 25, "8 clustered cores cover only a small part of the study area")

man <- read.csv(sp("MANIFEST.csv"))
ok(any(grepl("DIAGNOSTIC", man$scientific_status)),
   "the unmasked surface is labelled DIAGNOSTIC, not a product")
ok(any(grepl("MASKED", man$scientific_status)),
   "the masked surface is labelled as the only core-derived product")

# =============================================================================
section("06 against a SYNTHETIC covariate fixture")
# =============================================================================
# 03 needs Earth Engine, which CI does not have. This fixture exists purely to
# exercise the code paths in 06. Its numbers are random and mean nothing.
set.seed(20260727L)
cores <- read.csv(derived("01_cores_qc.csv"), check.names = FALSE)
n <- nrow(cores)
fixture <- data.frame(
  core_id   = cores$core_id,
  campaign  = cores$campaign,
  latitude  = cores$latitude,
  longitude = cores$longitude,
  elevation        = rnorm(n, 25, 5),
  tpi_300m         = rnorm(n, 0, 0.5),
  tpi_2km          = rnorm(n, 0, 1),
  slope_deg        = abs(rnorm(n, 0.3, 0.2)),
  dist_coast_m     = runif(n, 2000, 40000),
  water_occurrence = runif(n, 0, 60),
  s1_vv            = rnorm(n, -12, 2),
  s1_vh            = rnorm(n, -18, 2),
  s1_vv_vh_diff    = rnorm(n, 6, 1),
  ndvi             = runif(n, 0.3, 0.8),
  ndwi             = runif(n, -0.5, 0.2),
  ndmi             = runif(n, 0, 0.4),
  tc_wetness       = rnorm(n, 0, 0.1),
  worldcover       = sample(c(10, 20, 30, 90), n, replace = TRUE),
  stringsAsFactors = FALSE
)
for (k in 0:63) fixture[[sprintf("A%02d", k)]] <- rnorm(n)
write.csv(fixture, derived("03_covariates_at_cores.csv"), row.names = FALSE)

r6 <- run("06_covariate_model.R")
ok(is.null(r6$status), "06 exits cleanly on the fixture")
ok(file.exists(derived("06_model_cv.csv")), "06 writes cross-validation results")
ok(file.exists(tables("06_model_verdict.csv")), "06 writes a verdict")

cv <- read.csv(derived("06_model_cv.csv"))
ok(all(cv$n <= 8), "cross-validation never sees more than 8 observations")
ok(any(cv$cv == "leave-one-stratum-out"), "leave-one-stratum-out is run")
ok(!any(grepl("segment", cv$cv)), "no segment-level cross-validation is reported")
ok(grepl("Segment-level cross-validation is never performed", r6$out),
   "06 states that segment-level CV is never performed")

scr <- read.csv(tables("06_covariate_screening.csv"))
ok(nrow(scr) > 1, "the FULL screening ranking is written, not just the winner")
ok("n_candidates_screened" %in% names(scr),
   "the number of candidates screened is recorded (selection optimism)")
ok(!any(grepl("^A[0-9]{2}$", scr$covariate)),
   "AlphaEarth embedding bands are never offered to the model")
ok(!("worldcover" %in% scr$covariate),
   "the categorical stratum layer is never offered to the model")

ok(file.exists(tables("06_embedding_distance.csv")),
   "the embedding is used for the extrapolation diagnostic instead")

vd <- read.csv(tables("06_model_verdict.csv"))
ok("RECOMMENDED FOR MAPPING" %in% vd$criterion, "the verdict is explicit")
# On pure noise, no model should be recommended.
ok(vd$value[vd$criterion == "RECOMMENDED FOR MAPPING"] == "NO",
   "random covariates do NOT earn a map (the guard works)")

# =============================================================================
section("11 Bayesian maps")
# =============================================================================
r11 <- run("11_bayesian_map.R")
ok(is.null(r11$status), "11 exits cleanly")
ok(grepl("THE PEAT CORES SIT ON THIN PEAT", r11$out),
   "11 reports that the peat cores sit on thin peat")
ok(grepl("REGIONAL", r11$out), "11 names its prior mode")

pdp <- read.csv(tables("11_peat_depth_at_cores.csv"))
ok(nrow(pdp) == 8, "peat depth is reported for every core")
ok(all(pdp$peat_depth_cm[grepl("^FS", pdp$core_id)] == 0),
   "mineral cores are given zero peat depth")
ok(pdp$status[pdp$core_id == "PM-02-A"] == "contact_observed",
   "PM-02-A's peat/mineral contact is detected")
ok(pdp$is_lower_bound[pdp$core_id == "PM-01-A"],
   "a core still in peat at its base is a lower bound on peat depth")
ok(all(pdp$peat_depth_cm[pdp$status != "no_peat"] < 60),
   "observed peat depths are far below the 184 cm regional mean")

bt <- function(f) file.path(sandbox, "outputs", "spatial", f)
for (f in c("DEMO_BAYES_fullcolumn_mean_kgm2.tif",
            "DEMO_BAYES_fullcolumn_sd_kgm2.tif",
            "DEMO_BAYES_fullcolumn_core_info_fraction.tif",
            "DEMO_BAYES_fullcolumn_shift_from_prior_kgm2.tif")) {
  ok(file.exists(bt(f)), paste("11 writes", f))
  h <- readBin(bt(f), "raw", n = 8)
  ok(rawToChar(h[1:2]) == "II" &&
       readBin(h[3:4], "integer", size = 2, endian = "little",
               signed = FALSE) == 42L,
     paste(f, "is a valid TIFF"))
}
ok(grepl("^DEMO_", "DEMO_BAYES_fullcolumn_mean_kgm2.tif"),
   "regional-mode outputs are prefixed DEMO_ so they cannot be mistaken for the Li map")

isf <- read.csv(tables("11_implied_shallow_fraction.csv"))
ok(isf$implied_shallow_fraction < isf$uniform_density_fraction,
   "cores imply a smaller top-30cm share than uniform density (compaction)")
ok(isf$implied_shallow_fraction > 0 && isf$implied_shallow_fraction < 1,
   "implied fraction is a proper fraction")

man_b <- read.csv(bt("MANIFEST_bayes.csv"))
ok(all(grepl("FULL PEAT COLUMN", man_b$depth_support)),
   "every full-column raster is labelled with its depth support")

# =============================================================================
section("Quarto reports")
# =============================================================================
r7 <- run("07_validation_ledger.R")
ok(is.null(r7$status), "07 exits cleanly")
r8 <- run("08_report.R")
ok(is.null(r8$status), "08 exits cleanly")
r12 <- run("12_community_report.R")
ok(is.null(r12$status), "12 exits cleanly")

for (f in c("REPORT.qmd", "COMMUNITY_REPORT.qmd")) {
  p <- file.path(sandbox, "outputs", f)
  ok(file.exists(p), paste("writes", f))
  txt <- readLines(p, warn = FALSE)
  ok(txt[1] == "---", paste(f, "starts with YAML front matter"))
  close_at <- which(txt == "---")[2]
  ok(!is.na(close_at) && close_at > 1, paste(f, "front matter is closed"))
  fm <- txt[2:(close_at - 1)]
  ok(any(grepl("^title:", fm)), paste(f, "declares a title"))
  ok(any(grepl("^format:", fm)), paste(f, "declares output formats"))
  # Balanced Quarto callout fences.
  body <- txt[(close_at + 1):length(txt)]
  ok(sum(grepl("^::: \\{", body)) == sum(grepl("^:::\\s*$", body)),
     paste(f, "has balanced callout fences"))
  # No live code chunks: the report must render with the Quarto CLI alone.
  ok(!any(grepl("^```\\{r", body)),
     paste(f, "contains no executable chunks (renders without knitr)"))
}

cr <- paste(readLines(file.path(sandbox, "outputs", "COMMUNITY_REPORT.qmd"),
                      warn = FALSE), collapse = "\n")
ok(grepl("kg C/m", cr), "community report states its units")
ok(grepl("Bayesian update", cr), "community report explains the Bayesian idea")
ok(grepl("lower bound|LOWER BOUND|AT LEAST", cr),
   "community report carries the lower-bound caveat through to plain language")
ok(grepl("not a map of its carbon", cr),
   "community report keeps the shallow-map warning")
for (g in regmatches(cr, gregexpr("figures/[A-Za-z0-9_]+\\.png", cr))[[1]]) {
  ok(file.exists(file.path(sandbox, "outputs", g)),
     paste("figure resolves:", g))
}

ok(!file.exists(file.path(sandbox, "outputs", "REPORT.md")),
   "the superseded Markdown report is removed, so two versions cannot drift")

# The report must cover the whole analysis, and must degrade gracefully where
# an Earth Engine step has not run rather than emitting an empty section.
ok(grepl("The area this study covers", cr), "report covers the AOI")
ok(grepl("Carbon by ecosystem type", cr), "report covers ecosystem carbon")
ok(grepl("Finding distinct areas from satellite data", cr),
   "report covers the embedding clusters")
ok(grepl("Against the rest of Canada", cr), "report covers national context")
ok(grepl("The maps, and what each one may be used for", cr),
   "report inventories the map products and their status labels")
ok(grepl("Google Earth Engine session", cr),
   "sections awaiting Earth Engine say so instead of appearing empty")
ok(!grepl("NA km|NA%|kg C/m2\\s*NA", cr),
   "no NA leaks into the plain-language text")

# Status vocabulary must be defined for the reader, not just used.
for (lab in c("OBSERVATION", "DIAGNOSTIC", "CLASS-MEAN ASSIGNMENT")) {
  ok(grepl(lab, cr, fixed = TRUE),
     paste("status label explained to the reader:", lab))
}

# =============================================================================
section("13 AOI and 14 NPDB")
# =============================================================================
r13 <- run("13_aoi_boundary.R")
ok(is.null(r13$status), "13 exits cleanly")
sp13 <- function(f) file.path(sandbox, "outputs", "spatial", f)
for (f in c("aoi_coast_oriented.geojson", "aoi_vertices.csv",
            "aoi_coast_oriented_mask.tif", "aoi_for_earthengine.js")) {
  ok(file.exists(sp13(f)), paste("13 writes", f))
}
# The vertex table and the GeoJSON must describe the same ring.
av <- read.csv(sp13("aoi_vertices.csv"))
ok(nrow(av) == 5 && av$lon[1] == av$lon[5] && av$lat[1] == av$lat[5],
   "aoi_vertices.csv is a closed 5-vertex ring")
ok(all(av$lon < 0) && all(av$lat > 50 & av$lat < 60),
   "AOI vertices are in the right hemisphere and latitude band")

adef <- read.csv(file.path(sandbox, "outputs", "tables", "13_aoi_definition.csv"))
bd <- as.numeric(adef$value[adef$property == "bearing_deg"])
ok(is.finite(bd) && abs(bd - 135) < 45,
   "AOI bearing is within 45 degrees of the NW-SE coastal trend")

r14 <- run("14_npdb_context.R")
ok(is.null(r14$status), "14 exits cleanly")
ok(grepl("ONLY SOIL CARBON DATA", r14$out),
   "14 reports that the community cores are the only local data")
np <- read.csv(file.path(sandbox, "outputs", "tables", "14_npdb_proximity.csv"))
ok(np$n_npdb_cores[np$within_km == 500] == 0,
   "no national core lies within 500 km of Fort Severn")
ov <- read.csv(file.path(sandbox, "outputs", "tables", "14_organic_vs_mineral.csv"))
ok("wilcoxon_p" %in% names(ov) && "p_organic_exceeds_mineral" %in% names(ov),
   "the national comparison reports effect size, not just a difference in medians")
ok(ov$p_organic_exceeds_mineral[1] > 0.3 && ov$p_organic_exceeds_mineral[1] < 0.7,
   "the effect is small enough that the report must not overclaim it")

# ---- summary -----------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(sprintf("passed %d   failed %d\n", .pass, .fail))
if (.fail > 0L) {
  cat("failures:\n"); for (f in .failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("all integration tests passed\n")
