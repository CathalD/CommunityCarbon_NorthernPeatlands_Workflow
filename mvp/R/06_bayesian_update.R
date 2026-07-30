# =============================================================================
# 06_bayesian_update.R
#
# Fuse the prior map with the core residuals from step 5: each core pulls
# nearby pixels toward what was measured, fading with distance, and multiple
# cores combine by adding precision. Far from every core, the output equals
# the prior -- this step never invents structure the cores don't support.
#
# INPUT   mvp/outputs/current/prior_shallow_equiv.tif   (from 05)
#         mvp/outputs/current/prediction_residuals.csv  (from 05)
# OUTPUT  mvp/outputs/current/carbon_posterior_mean.tif
#         mvp/outputs/current/carbon_posterior_sd.tif
#         mvp/outputs/current/carbon_difference_from_prior.tif
#
# REQUIRES: terra
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  # Works three ways: Rscript (finds --file=), source() from ANY working
  # directory (finds the sourced file's own path), and run_all.R (sets
  # MVP_R_DIR). Only falls back to getwd() if all three fail.
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))
source(file.path(.this_dir, "bayes_core.R"))

library(terra)

msg("06  Bayesian fusion")

prior_path <- file.path(CFG$dir_current, "prior_shallow_equiv.tif")
resid_path <- file.path(CFG$dir_current, "prediction_residuals.csv")
if (!file.exists(prior_path)) stop("Missing ", prior_path, " -- run 05 first.")

prior_mean <- rast(prior_path)
residuals  <- read.csv(resid_path, stringsAsFactors = FALSE)

cores <- data.frame(
  core_id  = residuals$core_id,
  lon      = residuals$lon,
  lat      = residuals$lat,
  residual = residuals$residual,
  obs_sd   = CFG$bayes$core_sd_kgm2
)

# A residual needs a prior to be a residual OF. Where the published map is
# NoData at a core's location -- open water, or ground outside a peat-specific
# product's mask -- there is nothing to correct, so that core cannot enter the
# fusion. Excluded explicitly and by name, because leaving it in turns every
# pixel of the output NaN.
unusable <- !is.finite(cores$residual)
if (any(unusable)) {
  msg("EXCLUDED ", sum(unusable), " of ", nrow(cores),
     " core(s) with no prior value at their location: ",
     paste(cores$core_id[unusable], collapse = ", "))
  cores <- cores[!unusable, , drop = FALSE]
}
if (!nrow(cores)) {
  stop("No core has a finite residual -- the prior raster is NoData at every ",
      "core location. Check prior_shallow_equiv.tif against cores_clean.geojson ",
      "in QGIS; you may need a different external_prior in config.R.")
}
msg(nrow(cores), " core(s) contributing to the update")

length_scale_km <- CFG$bayes$length_scale_km
if (is.null(length_scale_km)) {
  length_scale_km <- median_nn_distance_km(cores$lon, cores$lat)
  msg("length scale not set in config; derived from cores: ",
     round(length_scale_km, 2), " km")
}

msg(sprintf("kernel length scale %.2f km, truncated beyond %.1f km",
           length_scale_km, CFG$bayes$max_influence_km))
msg(sprintf("prior sd %.1f kg/m2 (constant), core observation sd %.2f kg/m2",
           CFG$bayes$prior_sd_kgm2, CFG$bayes$core_sd_kgm2))

result <- bayes_update_raster(
  prior_mean_r = prior_mean,
  prior_sd_r   = CFG$bayes$prior_sd_kgm2,
  cores        = cores,
  length_scale_km  = length_scale_km,
  max_influence_km = CFG$bayes$max_influence_km
)

# A carbon stock cannot be negative. Where a large negative residual sits over
# a small prior value, the weighted correction can push a pixel below zero --
# arithmetically fine, physically impossible, and indefensible in something
# handed to a council. Clamp at zero, recompute the difference layer so it stays
# consistent with the clamped mean, and say how many cells it touched: a big
# number here means the residuals are fighting the prior hard enough that the
# depth assumption (shallow_fraction_of_column) deserves another look.
n_neg <- global(result$mean < 0, "sum", na.rm = TRUE)[1, 1]
if (!is.na(n_neg) && n_neg > 0) {
  n_valid_pre <- global(!is.na(result$mean), "sum", na.rm = TRUE)[1, 1]
  msg("clamped ", as.integer(n_neg), " cell(s) (",
     sprintf("%.2f%%", 100 * n_neg / n_valid_pre),
     ") from negative to 0 kg C/m2 -- a negative carbon stock is impossible")
  result$mean <- clamp(result$mean, lower = 0, values = TRUE)
  result$difference <- result$mean - prior_mean
}

names(result$mean) <- "carbon_posterior_mean_kgm2"
names(result$sd) <- "carbon_posterior_sd_kgm2"
names(result$difference) <- "difference_from_prior_kgm2"
names(result$info_frac) <- "core_info_fraction"
names(result$sd_reduction) <- "sd_reduction_kgm2"

writeRaster(result$mean, file.path(CFG$dir_current, "carbon_posterior_mean.tif"), overwrite = TRUE)
writeRaster(result$sd, file.path(CFG$dir_current, "carbon_posterior_sd.tif"), overwrite = TRUE)
writeRaster(result$difference, file.path(CFG$dir_current, "carbon_difference_from_prior.tif"), overwrite = TRUE)
writeRaster(result$info_frac, file.path(CFG$dir_current, "core_info_fraction.tif"), overwrite = TRUE)
writeRaster(result$sd_reduction, file.path(CFG$dir_current, "sd_reduction.tif"), overwrite = TRUE)

msg("wrote carbon_posterior_mean.tif, carbon_posterior_sd.tif, ",
   "carbon_difference_from_prior.tif, core_info_fraction.tif, sd_reduction.tif")

# ---- what the cores actually contributed, as numbers ------------------------
# This is the table to quote when asked "how much did these eight cores add to
# what was already known". Every figure in it is measured, not asserted.

cell_km2 <- cellSize(result$mean, unit = "km")
valid    <- !is.na(result$mean)
gsum <- function(r) global(r, "sum", na.rm = TRUE)[1, 1]
gmean <- function(r) global(r, "mean", na.rm = TRUE)[1, 1]

aoi_km2   <- gsum(valid * cell_km2)
inf       <- result$info_frac
prior_mean_kgm2 <- gmean(prior_mean)
post_mean_kgm2  <- gmean(result$mean)

contribution <- data.frame(
  metric = c(
    "cores contributing to the update",
    "mapped area with prior coverage (km2)",
    "AOI mean carbon, PRIOR (kg C/m2)",
    "AOI mean carbon, POSTERIOR (kg C/m2)",
    "change in AOI mean (kg C/m2)",
    "change in AOI mean (%)",
    "AOI mean uncertainty, PRIOR (kg C/m2)",
    "AOI mean uncertainty, POSTERIOR (kg C/m2)",
    "uncertainty reduction (kg C/m2)",
    "mean core information fraction (%)",
    "area where cores supply >1% of the information (km2)",
    "area where cores supply >10% of the information (km2)",
    "area where the map moved more than 0.1 kg C/m2 (km2)"),
  value = c(
    nrow(cores),
    round(aoi_km2, 1),
    round(prior_mean_kgm2, 3),
    round(post_mean_kgm2, 3),
    round(post_mean_kgm2 - prior_mean_kgm2, 4),
    round(100 * (post_mean_kgm2 / prior_mean_kgm2 - 1), 3),
    round(CFG$bayes$prior_sd_kgm2, 2),
    round(gmean(result$sd), 3),
    round(CFG$bayes$prior_sd_kgm2 - gmean(result$sd), 3),
    round(100 * gmean(inf), 3),
    round(gsum((inf > 0.01) * cell_km2), 1),
    round(gsum((inf > 0.10) * cell_km2), 1),
    round(gsum((abs(result$difference) > 0.1) * cell_km2), 1))
)
write.csv(contribution, file.path(CFG$dir_current, "core_contribution.csv"), row.names = FALSE)
msg("--- what the cores contributed ---")
print(contribution, right = FALSE)
msg("wrote core_contribution.csv")

# Report with global(na.rm=TRUE) rather than minmax(), which returns NaN on a
# raster carrying NoData.
lo <- global(result$mean, "min", na.rm = TRUE)[1, 1]
hi <- global(result$mean, "max", na.rm = TRUE)[1, 1]
msg(sprintf("posterior mean range: %.2f - %.2f kg C/m2", lo, hi))
msg(sprintf("mean |difference from prior|: %.3f kg C/m2",
           global(abs(result$difference), "mean", na.rm = TRUE)[1, 1]))

# How much of the AOI the prior actually covers. A peat-specific product is
# NoData over mineral ground and open water, and the posterior inherits that
# -- so this number tells you what fraction of the map exists at all.
n_cells <- ncell(result$mean)
n_valid <- global(!is.na(result$mean), "sum", na.rm = TRUE)[1, 1]
msg(sprintf("prior/posterior coverage: %.1f%% of the AOI (%d of %d cells)",
           100 * n_valid / n_cells, as.integer(n_valid), n_cells))
if (n_valid / n_cells < 0.5) {
  msg("NOTE: the prior is NoData over most of this AOI. That is expected for ",
     "a peat-specific product like Li et al., which maps peatland only -- but ",
     "it means the 'map' is blank wherever peat isn't mapped. Switch ",
     "CFG$bayes$external_prior to 'sothe_0_30' for a wall-to-wall surface.")
}

msg("06 complete")
