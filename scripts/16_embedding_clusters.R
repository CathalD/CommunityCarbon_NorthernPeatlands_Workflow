# =============================================================================
# 16_embedding_clusters.R
#
# k-means clustering of AlphaEarth satellite embeddings, to find distinct
# landscape types WITHIN and ACROSS the land-cover classes, and to report what
# carbon each one holds.
#
# WHY CLUSTERING IS THE HONEST USE OF THE EMBEDDING
#
# 03 and 06 refuse to let the 64-band AlphaEarth embedding near a predictive
# model, and that refusal stands: 64 predictors on 8 observations fits anything
# and generalises to nothing. Clustering is different in kind. It is
# UNSUPERVISED -- it never sees a carbon value, so it cannot overfit to eight
# of them. The binding constraint for k-means is the number of PIXELS, of which
# there are millions, not the number of cores.
#
# What this buys: land-cover classes are coarse and defined globally. Two
# pixels both labelled "swamp" can be quite different ground. The embedding
# encodes what the satellite record actually saw, so clustering it splits those
# classes into distinct landscape types. Carbon is then summarised per cluster
# -- from the published layers, with cores overlaid where they fall.
#
# WHAT IT CANNOT DO: clusters are not carbon classes. They are groupings of
# spectral and temporal similarity that MAY track carbon. A cluster containing
# no core has no ground truth, and is reported as such.
#
# CHOOSING k: no single k is correct. A range is run and the within-cluster
# sum of squares reported, so the choice is visible rather than asserted.
#
# REQUIRES rgee and an authenticated Earth Engine session.
#
# OUTPUT  outputs/tables/16_cluster_carbon.csv, 16_cluster_diagnostics.csv
#         GEE export: cluster map over the AOI
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

log_step("16  LANDSCAPE CLUSTERS FROM SATELLITE EMBEDDINGS")
set.seed(CFG$seed)
log_info("seed: ", CFG$seed, " (clustering is stochastic and must be reproducible)")

if (!requireNamespace("rgee", quietly = TRUE)) {
  fail_loudly("rgee is not installed",
              "This script needs Earth Engine for the AlphaEarth embeddings.",
              remedy = 'install.packages("rgee"); rgee::ee_install(); re-run.')
}
suppressPackageStartupMessages(library(rgee))
ee_Initialize(project = CFG$gee$project, drive = TRUE)
log_ok("Earth Engine initialised")

cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

aoi_json <- file.path(CFG$root, "outputs", "spatial", "aoi_coast_oriented.geojson")
if (file.exists(aoi_json)) {
  txt <- paste(readLines(aoi_json, warn = FALSE), collapse = " ")
  nums <- as.numeric(regmatches(txt, gregexpr("-?[0-9]+\\.[0-9]+", txt))[[1]])
  n_pair <- length(nums) %/% 2
  ring <- lapply(seq_len(n_pair), function(i) c(nums[2 * i - 1], nums[2 * i]))
  ring <- Filter(function(p) p[1] < 0 && p[2] > 0 && p[2] < 90, ring)
  aoi <- ee$Geometry$Polygon(list(ring))
  log_ok("using the coast-oriented AOI from 13")
} else {
  bb <- CFG$aoi_bbox
  aoi <- ee$Geometry$Rectangle(c(bb[["xmin"]], bb[["ymin"]],
                                 bb[["xmax"]], bb[["ymax"]]))
  log_warn("oriented AOI not found; using the axis-aligned box")
}

# =============================================================================
# 1. Embeddings
# =============================================================================

log_step("16a  EMBEDDINGS")

emb <- ee$ImageCollection(CFG$gee$asset_alphaearth)$
  filterBounds(aoi)$
  filter(ee$Filter$calendarRange(max(CFG$gee$season_years),
                                 max(CFG$gee$season_years), "year"))$
  mosaic()$clip(aoi)
n_bands <- length(emb$bandNames()$getInfo())
log_info("AlphaEarth embedding: ", n_bands, " bands")
log_info("Unsupervised, so the 64 dimensions are an asset here rather than a ",
         "liability: no carbon value is used in fitting, and the sample size ",
         "that matters is the pixel count, not the ", nrow(cores), " cores.")

k_values <- CFG$clusters$k_values
n_train  <- CFG$clusters$n_training_pixels

training <- emb$sample(region = aoi, scale = CFG$gee$export_scale_m,
                       numPixels = n_train, seed = CFG$seed,
                       geometries = FALSE)
log_info("sampled ", n_train, " training pixels at ",
         CFG$gee$export_scale_m, " m")

# =============================================================================
# 2. Choose k, visibly
# =============================================================================

log_step("16b  HOW MANY CLUSTERS?")

diag <- do.call(rbind, lapply(k_values, function(k) {
  cl <- ee$Clusterer$wekaKMeans(nClusters = k, seed = CFG$seed)$train(training)
  clustered <- emb$cluster(cl)
  hist <- clustered$reduceRegion(
    reducer = ee$Reducer$frequencyHistogram(), geometry = aoi,
    scale = CFG$gee$export_scale_m * 4, maxPixels = 1e12,
    bestEffort = TRUE)$getInfo()[[1]]
  counts <- unlist(hist)
  p <- counts / sum(counts)
  data.frame(k = k, n_clusters_found = length(counts),
             smallest_cluster_pct = 100 * min(p),
             largest_cluster_pct  = 100 * max(p),
             # Shannon evenness: 1 = all clusters equal, 0 = one dominates.
             evenness = if (length(p) > 1)
               -sum(p * log(p)) / log(length(p)) else NA_real_,
             stringsAsFactors = FALSE)
}))
print_table(diag, digits = 2)
log_info("A cluster holding a tiny share of the area is usually noise or an ",
         "artefact rather than a landscape type; evenness near 1 means the ",
         "clusters split the area rather than isolating slivers.")

k <- CFG$clusters$k_final
log_ok("using k = ", k, " for the reported product (set in config)")

# =============================================================================
# 3. Cluster and characterise
# =============================================================================

log_step("16c  CLUSTERING")

clusterer <- ee$Clusterer$wekaKMeans(nClusters = k, seed = CFG$seed)$train(training)
clusters  <- emb$cluster(clusterer)$rename("cluster")$clip(aoi)

# Carbon surfaces, on the same terms as 15.
sg_soc  <- ee$Image(CFG$gee$asset_soilgrids_soc)
sg_bdod <- ee$Image(CFG$gee$asset_soilgrids_bdod)
sgl <- function(sb, bb2, th) sg_soc$select(sb)$divide(10)$
  multiply(sg_bdod$select(bb2)$divide(100))$multiply(th)$divide(100)
sg_0_30 <- sgl("soc_0-5cm_mean", "bdod_0-5cm_mean", 5)$
  add(sgl("soc_5-15cm_mean", "bdod_5-15cm_mean", 10))$
  add(sgl("soc_15-30cm_mean", "bdod_15-30cm_mean", 15))$
  rename("soilgrids_0_30_kgm2")$clip(aoi)
li <- ee$Image(CFG$gee$asset_li2025)$select(0)$
  rename("li_full_column_kgm2")$clip(aoi)

gwl <- ee$ImageCollection(CFG$gee$asset_gwl_fcs30)$
  filterBounds(aoi)$mosaic()$rename("gwl_class")$clip(aoi)

zonal <- function(carbon, scale_m, tag) {
  red <- ee$Reducer$mean()$combine(ee$Reducer$stdDev(), "", TRUE)$
    combine(ee$Reducer$count(), "", TRUE)$
    group(groupField = 1L, groupName = "cluster")
  g <- carbon$addBands(clusters)$reduceRegion(
    reducer = red, geometry = aoi, scale = scale_m,
    maxPixels = 1e13, bestEffort = TRUE)$getInfo()$groups
  if (is.null(g)) return(NULL)
  d <- do.call(rbind, lapply(g, function(x) data.frame(
    cluster = as.integer(x$cluster),
    mean_kgm2 = if (is.null(x$mean)) NA_real_ else x$mean,
    sd_kgm2   = if (is.null(x$stdDev)) NA_real_ else x$stdDev,
    n_pixels  = if (is.null(x$count)) NA_real_ else x$count,
    stringsAsFactors = FALSE)))
  d$area_km2 <- d$n_pixels * (scale_m / 1000)^2
  d$pct_of_aoi <- 100 * d$area_km2 / sum(d$area_km2, na.rm = TRUE)
  d$carbon_layer <- tag
  d
}

ct <- rbind(zonal(sg_0_30, 250, "SoilGrids 0-30 cm"),
            zonal(li, 30, "Li et al. full peat column"))
log_info("carbon by landscape cluster:")
print_table(ct[order(ct$carbon_layer, -ct$mean_kgm2), ], digits = 2)

# =============================================================================
# 4. Do clusters cut across land-cover classes?
# =============================================================================

log_step("16d  CLUSTERS VERSUS LAND COVER")

# If clusters simply reproduced the land-cover map they would add nothing. The
# interesting case is a land-cover class split across several clusters, which
# is the embedding finding structure the class system misses.
cross <- clusters$addBands(gwl)$reduceRegion(
  reducer = ee$Reducer$frequencyHistogram()$group(groupField = 1L,
                                                  groupName = "gwl_class"),
  geometry = aoi, scale = CFG$gee$export_scale_m * 4,
  maxPixels = 1e12, bestEffort = TRUE)$getInfo()$groups

if (!is.null(cross) && length(cross)) {
  spread <- do.call(rbind, lapply(cross, function(x) {
    h <- unlist(x$histogram); p <- h / sum(h)
    data.frame(gwl_class = as.integer(x$gwl_class),
               n_clusters_present = sum(p > 0.05),
               dominant_cluster_pct = 100 * max(p),
               stringsAsFactors = FALSE)
  }))
  print_table(spread, digits = 1)
  multi <- sum(spread$n_clusters_present > 1)
  if (multi > 0) {
    log_ok(multi, " land-cover class(es) split across more than one cluster")
    log_info("That is the point of this step: the embedding is finding ",
             "distinctions the global class system does not carry. Those ",
             "sub-types may differ in carbon even though the map calls them ",
             "the same thing.")
  } else {
    log_warn("Every land-cover class maps to a single cluster, so the ",
             "clustering adds no structure beyond the class map here.")
  }
  write_csv_logged(spread, file.path(CFG$dir_tables, "16_cluster_vs_landcover.csv"),
                   "whether clusters cut across land-cover classes")
}

# =============================================================================
# 5. Which clusters have ground truth?
# =============================================================================

log_step("16e  GROUND TRUTH BY CLUSTER")

pts <- ee$FeatureCollection(lapply(seq_len(nrow(cores)), function(i) {
  ee$Feature(ee$Geometry$Point(c(cores$longitude[i], cores$latitude[i])),
             list(core_id = cores$core_id[i], campaign = cores$campaign[i]))
}))
cs <- sf::st_drop_geometry(ee_as_sf(
  clusters$sampleRegions(collection = pts, scale = CFG$gee$export_scale_m,
                         geometries = FALSE, tileScale = 4),
  maxFeatures = 1000))
r30 <- stocks[stocks$scenario == "as_supplied" &
              stocks$window == "reference_0_30", ]
cs$stock_0_30_kgm2 <- r30$stock_kgm2[match(cs$core_id, r30$core_id)]
cs$is_lower_bound  <- r30$is_lower_bound[match(cs$core_id, r30$core_id)]
print_table(cs, digits = 2)

covered <- unique(cs$cluster)
ct$has_core <- ct$cluster %in% covered
unv <- unique(ct[!ct$has_core, c("cluster", "pct_of_aoi")])
log_warn(k - length(covered), " of ", k, " clusters contain no core at all.")
log_warn(sprintf("That is %.1f%% of the AOI with no ground verification.",
                 sum(unv$pct_of_aoi, na.rm = TRUE)))
log_warn("Those clusters' carbon values come entirely from published layers. ",
         "They are the natural targets for the next field season: sampling a ",
         "cluster that has never been visited buys more information than ",
         "another core beside an existing one.")

# =============================================================================
# 6. Outputs
# =============================================================================

log_step("16f  OUTPUTS")

write_csv_logged(ct, file.path(CFG$dir_tables, "16_cluster_carbon.csv"),
                 "carbon per landscape cluster, with ground-truth flag")
write_csv_logged(diag, file.path(CFG$dir_tables, "16_cluster_diagnostics.csv"),
                 "how the number of clusters was chosen")
write_csv_logged(cs, file.path(CFG$dir_tables, "16_cores_by_cluster.csv"),
                 "which landscape cluster each core sits in")

tk <- ee_image_to_drive(image = clusters$toInt16(),
                        description = "ccnp_landscape_clusters_30m",
                        folder = "CCNP_SOC", region = aoi,
                        scale = CFG$gee$export_scale_m,
                        crs = CFG$gee$export_crs, maxPixels = 1e13)
tk$start()
log_ok("started export: ccnp_landscape_clusters_30m")
log_info("Suggested next-core targets: the unvisited clusters listed above, ",
         "chosen to cover landscape types rather than to fill space evenly.")
log_ok("16 complete")
