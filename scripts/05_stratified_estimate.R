# =============================================================================
# 05_stratified_estimate.R
#
# The sample-based products. These are the only carbon numbers in this
# workflow that rest entirely on the community cores.
#
# PRODUCT 1  Stratified design-based estimate of mean SOC stock per stratum.
#            Status: SAMPLE-BASED ESTIMATE. Not a map.
#
# PRODUCT 2  Stratum-mean assignment map: each stratum's mean stock painted
#            onto a stratum map.
#            Status: CLASS-MEAN ASSIGNMENT. This is NOT interpolation -- no
#            distance-weighting of samples occurs -- and NOT prediction from
#            covariates. It is a design-based mean projected onto a
#            categorical layer, and it is only as good as the assumption that
#            the stratum is internally homogeneous.
#
# Also produced here: the evidence for why geostatistical interpolation
# (kriging, IDW with a fitted structure) is REFUSED for this dataset. That
# refusal is a finding, so it is demonstrated rather than asserted.
#
# INPUT   data/derived/01_cores_qc.csv, 02_core_stocks.csv
# OUTPUT  data/derived/05_stratum_estimates.csv
#         outputs/tables/05_spatial_structure.csv
#         outputs/tables/05_variogram_binned.csv
#         outputs/tables/05_product2_map_lookup.csv
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

log_step("05  DESIGN-BASED ESTIMATION AND THE CASE AGAINST INTERPOLATION")

stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")
cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")

set.seed(CFG$seed)
log_info("seed: ", CFG$seed)

# =============================================================================
# 1. Spatial structure: can these points support interpolation at all?
# =============================================================================

log_step("05a  SPATIAL STRUCTURE DIAGNOSTIC")

nn <- nn_distance_summary(cores$longitude, cores$latitude, cores$core_id)
nn$campaign <- cores$campaign[match(nn$id, cores$core_id)]
log_info("nearest-neighbour distance for each core (km):")
print_table(nn[order(nn$nn_dist_km), ], digits = 2)

d <- pairwise_km(cores$longitude, cores$latitude)
extent_km <- max(d)
log_info(sprintf("study-area extent spanned by cores: %.1f km", extent_km))
log_info(sprintf("nearest-neighbour distance: median %.1f km, range %.1f-%.1f km",
                 median(nn$nn_dist_km), min(nn$nn_dist_km), max(nn$nn_dist_km)))

# Variogram on the window where every core is valid.
cs <- stocks[stocks$scenario == "as_supplied" &
             stocks$window == "common_support", ]
cs <- cs[match(cores$core_id, cs$core_id), ]

vc <- variogram_cloud(cores$longitude, cores$latitude, cs$stock_kgm2)
vb <- variogram_binned(vc, n_bins = 5L)
log_info(sprintf("variogram: %d cores give %d unique pairs", nrow(cores), nrow(vc)))
print_table(vb, digits = 3)

n_informative <- sum(vb$informative)
log_warn(strrep("-", 68))
log_warn("GEOSTATISTICAL INTERPOLATION IS REFUSED FOR THIS DATASET.")
log_warn(sprintf("  %d cores yield %d pairs in total.", nrow(cores), nrow(vc)))
log_warn(sprintf("  %d of %d lag bins reach the ~30 pairs conventionally",
                 n_informative, nrow(vb)))
log_warn("  required before a bin is treated as informative.")
log_warn(sprintf("  Nearest-neighbour spacing (%.1f-%.1f km) is of the same",
                 min(nn$nn_dist_km), max(nn$nn_dist_km)))
log_warn(sprintf("  order as the whole sampled extent (%.1f km), so there is",
                 extent_km))
log_warn("  no short-lag information from which to estimate a nugget or a")
log_warn("  range. A fitted variogram here would be an artefact of the")
log_warn("  fitting routine, not a measurement of spatial structure.")
log_warn("")
log_warn("  Worse, the two campaigns form two spatial clusters that differ in")
log_warn("  landscape type. Any apparent spatial correlation would largely be")
log_warn("  the peat/mineral contrast reappearing as a distance effect.")
log_warn("")
log_warn("  Consequence: NO kriged or IDW surface is produced. The sample-")
log_warn("  based product is a stratified design-based estimate instead,")
log_warn("  which makes its homogeneity assumption explicit rather than")
log_warn("  hiding it inside a covariance model.")
log_warn(strrep("-", 68))

spatial_structure <- data.frame(
  metric = c("n_cores", "n_pairs", "extent_km", "nn_median_km", "nn_min_km",
             "nn_max_km", "n_informative_lag_bins", "kriging_supported"),
  value  = c(nrow(cores), nrow(vc), round(extent_km, 2),
             round(median(nn$nn_dist_km), 2), round(min(nn$nn_dist_km), 2),
             round(max(nn$nn_dist_km), 2), n_informative, 0),
  stringsAsFactors = FALSE
)

# =============================================================================
# 2. PRODUCT 1 -- stratified design-based estimates
# =============================================================================

log_step("05b  PRODUCT 1: STRATIFIED DESIGN-BASED ESTIMATE")

# Area weights are needed to form a single study-area mean. They come from a
# stratum raster exported by 03. GWL_FCS30 is preferred over WorldCover because
# WorldCover cannot separate treed bog on peat from upland forest on mineral
# soil, which is exactly the distinction the two campaigns represent.
#
# The class-to-stratum mapping below is a modelling decision, not a fact, and
# it is stated here rather than buried: see `lookup` in section 3 for the
# caveats that travel with it.
area_weights <- NULL
strat_sources <- list(
  list(path = file.path(CFG$dir_gee, "ccnp_gwl_fcs30_30m.tif"),
       name = "GWL_FCS30",
       peat = CFG$gee$gwl_wetland_codes,     # wetland classes -> peat stratum
       mineral = c(180L)),                   # non-wetland      -> mineral stratum
  list(path = file.path(CFG$dir_gee, "ccnp_worldcover_30m.tif"),
       name = "ESA WorldCover",
       peat = c(30L, 90L),                   # grassland (wet), herbaceous wetland
       mineral = c(10L, 20L))                # tree cover, shrubland
)
strat <- Filter(function(s) file.exists(s$path), strat_sources)

if (length(strat) && requireNamespace("terra", quietly = TRUE)) {
  src <- strat[[1]]
  log_info("stratum areas read from ", basename(src$path), " (", src$name, ")")
  r  <- terra::rast(src$path)[[1]]
  tb <- terra::freq(r)
  n_peat    <- sum(tb$count[tb$value %in% src$peat],    na.rm = TRUE)
  n_mineral <- sum(tb$count[tb$value %in% src$mineral], na.rm = TRUE)
  n_other   <- sum(tb$count, na.rm = TRUE) - n_peat - n_mineral

  if (n_peat > 0 && n_mineral > 0) {
    area_weights <- c(PM_2024_peat = n_peat, FS_2025_mineral = n_mineral)
    log_ok(sprintf("stratum composition: peat %.1f%%, mineral %.1f%%, unclassified %.1f%%",
                   100 * n_peat / sum(tb$count), 100 * n_mineral / sum(tb$count),
                   100 * n_other / sum(tb$count)))
    log_info("a study-area mean IS formed, weighted by mapped AREA rather than ",
             "by how many cores happened to fall in each stratum")
    if (n_other > 0.2 * sum(tb$count)) {
      log_warn(sprintf(paste0("%.0f%% of the area falls in neither stratum and ",
                              "is EXCLUDED from the weighted mean; the study-",
                              "area figure therefore describes only the mapped ",
                              "portion."), 100 * n_other / sum(tb$count)))
    }
  } else {
    log_warn("stratum raster present but one stratum has zero mapped area; ",
             "no study-area mean is formed")
  }
} else {
  lc_path <- strat_sources[[1]]$path
  log_warn("No stratum raster found at ", lc_path)
  log_warn("Per-stratum estimates are reported. NO study-area mean is formed:")
  log_warn("weighting strata by how many cores fell in them would substitute")
  log_warn("field logistics for landscape composition, which is exactly the")
  log_warn("error that makes a convenience sample look like a survey.")
}

est_list <- list()
for (sc in unique(stocks$scenario)) {
  for (w in unique(stocks$window)) {
    d0 <- stocks[stocks$scenario == sc & stocks$window == w, ]

    # Cores whose stock is a lower bound cannot enter a mean as if it were a
    # point estimate. They are excluded from the mean and counted separately.
    usable <- d0[!d0$is_lower_bound, ]
    n_excl <- sum(d0$is_lower_bound)

    e <- stratified_estimate(usable$stock_kgm2, usable$campaign,
                             area_weight = area_weights,
                             B = CFG$boot_n, conf = CFG$boot_conf,
                             seed = CFG$seed)
    e$scenario <- sc
    e$window   <- w
    e$n_excluded_lower_bound <- n_excl
    e$excluded_cores <- paste(d0$core_id[d0$is_lower_bound], collapse = ";")
    est_list[[length(est_list) + 1L]] <- e
  }
}
est <- do.call(rbind, est_list)
rownames(est) <- NULL

log_info("stratum means, kg C/m2 (lower-bound cores excluded from the mean):")
for (sc in unique(est$scenario)) {
  for (w in unique(est$window)) {
    e <- est[est$scenario == sc & est$window == w, ]
    log_info(sprintf("  scenario=%s window=%s  (excluded as lower bounds: %s)",
                     sc, w,
                     ifelse(nzchar(e$excluded_cores[1]), e$excluded_cores[1], "none")))
    print_table(e[, c("stratum", "n", "mean", "sd", "se", "cv_pct",
                      "boot_lo", "boot_hi", "t_lo", "t_hi", "trustworthy")],
                digits = 3)
  }
}

log_warn(strrep("-", 68))
log_warn("READING THESE INTERVALS")
log_warn("  Every stratum here has n between 1 and 5. A bootstrap resampling")
log_warn("  three numbers can only produce a handful of distinct means; the")
log_warn("  interval it returns describes those three numbers, not the")
log_warn("  population they came from. `trustworthy` is FALSE throughout for")
log_warn("  that reason, and the t interval is printed alongside so the")
log_warn("  dependence on method stays visible.")
log_warn("  These are honest summaries of eight cores. They are not")
log_warn("  precise estimates of landscape means, and no amount of")
log_warn("  resampling can make them so.")
log_warn(strrep("-", 68))

# =============================================================================
# 3. PRODUCT 2 -- stratum-mean assignment map
# =============================================================================

log_step("05c  PRODUCT 2: STRATUM-MEAN ASSIGNMENT MAP")

# The two sampled strata correspond to a peat/wetland class and a
# forest/mineral class. Mapping them onto ESA WorldCover requires an explicit,
# arguable correspondence, which is stated here rather than buried.
lookup <- data.frame(
  stratum = c("PM_2024_peat", "FS_2025_mineral"),
  worldcover_classes = c("90 (herbaceous wetland), 30 (grassland) where wet",
                         "10 (tree cover), 20 (shrubland)"),
  correspondence_confidence = c("moderate", "moderate"),
  caveat = c(
    "WorldCover does not distinguish peat depth or peatland type; a herbaceous wetland class covers both shallow fen and deep bog, which differ greatly in carbon.",
    "Tree cover in the Hudson Bay Lowlands includes treed bog on peat as well as upland forest on mineral soil. This class is NOT a clean mineral-soil mask."
  ),
  stringsAsFactors = FALSE
)

ref <- est[est$scenario == "as_supplied" &
           est$window == "common_support" &
           est$stratum %in% lookup$stratum, ]
lookup$mean_stock_kgm2 <- ref$mean[match(lookup$stratum, ref$stratum)]
lookup$se_kgm2         <- ref$se[match(lookup$stratum, ref$stratum)]
lookup$n_cores         <- ref$n[match(lookup$stratum, ref$stratum)]
lookup$window          <- "common_support (0-14.5 cm)"
lookup$scientific_status <- "CLASS-MEAN ASSIGNMENT (design-based). Not interpolation; not covariate prediction."

log_info("map legend -- the product IS this lookup applied to a stratum raster:")
print_table(lookup[, c("stratum", "n_cores", "mean_stock_kgm2", "se_kgm2",
                       "worldcover_classes")], digits = 3)

log_warn("Two limits travel with this product and must travel with the map:")
log_warn("  1. Within-stratum variation is discarded by construction.")
for (st in lookup$stratum) {
  v <- cs$stock_kgm2[cs$campaign == st]
  if (length(v) > 1L) {
    log_warn(sprintf(paste0("     %-16s ranges %.2f-%.2f kg/m2 across %d ",
                            "cores (%.1f-fold);"),
                     st, min(v), max(v), length(v), max(v) / min(v)))
    log_warn(sprintf("                      one mean of %.2f replaces all of it.",
                     mean(v)))
  }
}
log_warn("  2. Campaign, landscape type, year, operator and laboratory")
log_warn("     convention are perfectly confounded in this design. The")
log_warn("     peat/mineral difference this map draws CANNOT be separated")
log_warn("     from a between-campaign difference. It is a difference")
log_warn("     between two groups of cores, labelled by landscape.")

# Rendering the map needs terra only to READ the exported stratum raster.
# Reading arbitrary GeoTIFFs (tiled, compressed, however GEE wrote them) is
# well outside what R/geotiff.R attempts, so this step is the one place the
# workflow genuinely benefits from GDAL.
if (length(strat) && requireNamespace("terra", quietly = TRUE)) {
  src <- strat[[1]]
  log_info("stratum raster present (", src$name, "); writing the assignment map")

  peat_mean    <- lookup$mean_stock_kgm2[lookup$stratum == "PM_2024_peat"]
  mineral_mean <- lookup$mean_stock_kgm2[lookup$stratum == "FS_2025_mineral"]
  rcl <- rbind(
    cbind(src$peat,    rep(peat_mean,    length(src$peat))),
    cbind(src$mineral, rep(mineral_mean, length(src$mineral)))
  )

  r <- terra::rast(src$path)[[1]]
  m <- terra::classify(r, rcl, others = NA)
  names(m) <- "soc_kgm2_class_mean"
  out <- file.path(CFG$dir_gee, "PRODUCT2_stratum_mean_soc_kgm2.tif")
  terra::writeRaster(m, out, overwrite = TRUE)
  log_ok("wrote ", basename(out), " using the ", src$name, " class mapping")
  log_warn("Classes outside the mapping are NoData, not zero. A pixel this ",
           "map leaves blank is one the cores say nothing about.")
} else if (length(strat)) {
  log_warn("stratum raster found but `terra` is not installed, so the raster ",
           "cannot be read; only the legend is written.")
} else {
  log_info("stratum raster not present, so only the legend is written.")
  log_info("Run 03_gee_covariates.R, retrieve ccnp_gwl_fcs30_30m.tif (preferred) ",
           "or ccnp_worldcover_30m.tif into ", CFG$dir_gee,
           ", then re-run this script to render the map.")
}

# =============================================================================
# 4. Write
# =============================================================================

log_step("05d  OUTPUTS")

write_csv_logged(est, file.path(CFG$dir_derived, "05_stratum_estimates.csv"),
                 "PRODUCT 1: design-based stratum means")
write_csv_logged(spatial_structure,
                 file.path(CFG$dir_tables, "05_spatial_structure.csv"),
                 "evidence that kriging is not supportable")
write_csv_logged(vb, file.path(CFG$dir_tables, "05_variogram_binned.csv"))
write_csv_logged(lookup,
                 file.path(CFG$dir_tables, "05_product2_map_lookup.csv"),
                 "PRODUCT 2: stratum-mean map legend")

log_ok("05 complete")
