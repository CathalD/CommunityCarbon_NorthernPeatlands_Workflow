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

names(result$mean) <- "carbon_posterior_mean_kgm2"
names(result$sd) <- "carbon_posterior_sd_kgm2"
names(result$difference) <- "difference_from_prior_kgm2"

writeRaster(result$mean, file.path(CFG$dir_current, "carbon_posterior_mean.tif"), overwrite = TRUE)
writeRaster(result$sd, file.path(CFG$dir_current, "carbon_posterior_sd.tif"), overwrite = TRUE)
writeRaster(result$difference, file.path(CFG$dir_current, "carbon_difference_from_prior.tif"), overwrite = TRUE)

msg("wrote carbon_posterior_mean.tif, carbon_posterior_sd.tif, carbon_difference_from_prior.tif")

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
