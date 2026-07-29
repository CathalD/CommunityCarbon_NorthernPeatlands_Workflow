# =============================================================================
# bayes_core.R  --  conjugate Bayesian fusion of a prior carbon map with new
# core observations. The one piece of genuinely custom math in this pipeline
# -- no package implements "prior raster + point observations + a distance
# kernel" as a single call. Trimmed from the original repo's R/bayes.R
# (which does the same thing at ~5x the length, for a research audience).
#
# THE MODEL
#
#   The prior gives, at every pixel, a value and an uncertainty. Each new
#   core observation is compared against the prior's value AT THAT POINT --
#   the residual -- and that residual is spread across nearby pixels with a
#   Gaussian kernel that fades with distance and is truncated beyond
#   max_influence_km. Multiple cores combine by adding precisions (1/sd^2),
#   so more/closer cores pull harder and the posterior uncertainty shrinks
#   correctly with sample size. Far from every core, the correction is ~0 and
#   the output IS the prior, unchanged -- this is what makes the fusion safe
#   to run again next season with only a couple of new cores: it never
#   invents structure the data doesn't support.
#
#   updated_map = prior_map + precision_weighted_mean(residuals, kernel)
# =============================================================================

#' Gaussian distance kernel. weight=1 at distance 0, exp(-0.5)=0.61 at the
#' length scale, ~0 beyond max_km.
kernel_weight <- function(dist_km, length_scale_km, max_km) {
  w <- exp(-0.5 * (dist_km / length_scale_km)^2)
  w[dist_km > max_km] <- 0
  w
}

#' Fuse a prior raster with core residuals. Returns posterior mean, sd, and
#' the correction (= posterior - prior, i.e. the "difference from prior" map).
#'
#' @param prior_mean_r  SpatRaster, the prior's expected value
#' @param prior_sd_r    SpatRaster or single number, the prior's uncertainty
#' @param cores         data.frame with columns lon, lat, residual, obs_sd
#' @param length_scale_km,max_influence_km  kernel parameters
bayes_update_raster <- function(prior_mean_r, prior_sd_r, cores,
                                length_scale_km, max_influence_km) {
  if (is.numeric(prior_sd_r) && length(prior_sd_r) == 1) {
    prior_sd_r <- terra::setValues(prior_mean_r, prior_sd_r)
  }
  tau_prior <- 1 / (prior_sd_r^2)

  tau_data <- terra::setValues(prior_mean_r, 0)
  weighted_resid <- terra::setValues(prior_mean_r, 0)

  for (i in seq_len(nrow(cores))) {
    pt <- terra::vect(cbind(cores$lon[i], cores$lat[i]),
                      type = "points", crs = terra::crs(prior_mean_r))
    dist_km <- terra::distance(prior_mean_r, pt) / 1000
    w <- kernel_weight(dist_km, length_scale_km, max_influence_km)
    tau_i <- w / (cores$obs_sd[i]^2)
    tau_data <- tau_data + tau_i
    weighted_resid <- weighted_resid + tau_i * cores$residual[i]
  }

  tau_post <- tau_prior + tau_data
  correction <- weighted_resid / tau_post   # 0 where no core has influence
  sd_post <- 1 / sqrt(tau_post)

  list(
    mean       = prior_mean_r + correction,
    sd         = sd_post,
    difference = correction
  )
}

#' Median nearest-neighbour distance between points, km. Used as the default
#' kernel length scale when none is configured.
median_nn_distance_km <- function(lon, lat) {
  n <- length(lon)
  if (n < 2) return(5)  # arbitrary small default for a single core
  d <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    d[i, ] <- terra::distance(
      cbind(lon[i], lat[i]), cbind(lon, lat), lonlat = TRUE) / 1000
  }
  diag(d) <- NA
  median(apply(d, 1, min, na.rm = TRUE))
}
