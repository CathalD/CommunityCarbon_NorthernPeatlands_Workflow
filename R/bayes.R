# =============================================================================
# R/bayes.R  --  Conjugate Bayesian updating of a published carbon map with
#                community core observations.
#
# THE MODEL, IN WORDS
#
# A published map gives, at every pixel, a best guess at the carbon stock and an
# uncertainty around it. That is the PRIOR. The community cores are direct
# measurements at eight points. Where a core is close to a pixel it should pull
# that pixel's value strongly towards what was measured; far away it should not
# move it at all, and the pixel should keep the published value.
#
# Written as a conjugate Gaussian update with known variances:
#
#   prior            S(x) ~ Normal( mu0(x), sigma0(x)^2 )
#   observation i    y_i  ~ Normal( S(x),   sigma_y_i^2 / w_i(x) )
#
# where w_i(x) is a distance kernel in [0, 1]. Dividing the observation
# variance by w is what makes a distant core uninformative: as w -> 0 its
# variance -> infinity and it contributes nothing. Precisions then simply add:
#
#   tau_post(x) = 1/sigma0(x)^2 + SUM_i w_i(x)/sigma_y_i^2
#   mu_post(x)  = [ mu0(x)/sigma0(x)^2 + SUM_i w_i(x) y_i / sigma_y_i^2 ] / tau_post(x)
#   sigma_post(x) = 1 / sqrt( tau_post(x) )
#
# SAMPLE SIZE enters automatically and correctly: each core contributes its own
# precision term, so n nearby cores contribute n times the precision of one,
# and the posterior tightens as sqrt(n). Nothing extra is needed for that.
#
# The most useful output is not the posterior mean but `info_frac`: the share
# of the posterior precision that came from the cores rather than the prior. It
# maps exactly where the community data actually changed the answer, and it is
# near zero across most of a study area sampled by eight points.
#
# WHAT THIS IS NOT
#
# This is not kriging and does not pretend to estimate a spatial covariance
# structure from the data -- 05 shows why that is not possible here. The kernel
# length scale is supplied, not fitted, and the whole construction is an
# explicit, auditable weighting rather than an inferred one.
#
# All functions are pure.
# =============================================================================

#' Gaussian distance kernel, returning weights in (0, 1].
#'
#' At distance 0 the weight is 1; at the length scale it is exp(-0.5) = 0.61;
#' at twice the length scale, 0.14. Chosen over an exponential or a hard radius
#' because it is smooth (no ring artefacts in the output raster) and has one
#' interpretable parameter.
kernel_gaussian <- function(dist_km, length_scale_km) {
  if (length_scale_km <= 0) stop("length_scale_km must be positive", call. = FALSE)
  exp(-0.5 * (dist_km / length_scale_km)^2)
}

#' Kernel truncated to zero beyond `max_km`.
#'
#' Without truncation every core has some influence everywhere, which produces
#' a map that looks informed across its whole extent when it is not. Truncation
#' makes the honest statement: beyond this distance, the cores say nothing and
#' the map is the published product unaltered.
kernel_truncated <- function(dist_km, length_scale_km, max_km) {
  w <- kernel_gaussian(dist_km, length_scale_km)
  w[dist_km > max_km] <- 0
  w
}

#' Conjugate Gaussian posterior for a grid of pixels and a set of point
#' observations.
#'
#' @param mu0,sigma0  numeric vectors, length n_pixel: prior mean and sd.
#' @param y,sigma_y   numeric vectors, length n_core: observed stock and its sd.
#' @param W           n_pixel x n_core matrix of kernel weights in [0, 1].
#'
#' @return data frame with one row per pixel:
#'   mean, sd, info_frac (share of precision from cores), shift (mean - mu0).
#'
#' Pixels where the prior is missing are handled: with no prior, the posterior
#' rests entirely on the cores, and where there is also no core weight the
#' result is NA rather than a fabricated number.
bayes_posterior <- function(mu0, sigma0, y, sigma_y, W) {
  n_pix <- length(mu0)
  stopifnot(length(sigma0) == n_pix,
            is.matrix(W), nrow(W) == n_pix, ncol(W) == length(y),
            length(sigma_y) == length(y))
  if (any(sigma_y <= 0, na.rm = TRUE)) {
    stop("sigma_y must be strictly positive", call. = FALSE)
  }
  if (any(W < 0 | W > 1, na.rm = TRUE)) {
    stop("kernel weights must lie in [0, 1]", call. = FALSE)
  }

  # Prior precision. A missing or non-positive prior sd contributes nothing.
  tau0 <- ifelse(is.finite(mu0) & is.finite(sigma0) & sigma0 > 0,
                 1 / sigma0^2, 0)
  mu0_safe <- ifelse(is.finite(mu0), mu0, 0)

  # Per-core precision at each pixel: w / sigma_y^2.
  Wc <- W
  Wc[!is.finite(Wc)] <- 0
  prec <- sweep(Wc, 2L, sigma_y^2, "/")

  tau_data <- rowSums(prec)
  num_data <- as.vector(prec %*% y)

  tau_post <- tau0 + tau_data
  mu_post  <- (mu0_safe * tau0 + num_data) / tau_post

  # No prior and no data: nothing is known, and the result must say so.
  none <- tau_post <= 0
  mu_post[none]  <- NA_real_
  sd_post <- ifelse(none, NA_real_, 1 / sqrt(tau_post))

  data.frame(
    mean      = mu_post,
    sd        = sd_post,
    info_frac = ifelse(none, NA_real_, tau_data / tau_post),
    shift     = ifelse(none | !is.finite(mu0), NA_real_, mu_post - mu0),
    stringsAsFactors = FALSE
  )
}

#' Build the pixel-by-core weight matrix, optionally gated by stratum.
#'
#' STRATUM GATING is on by default and matters. A mineral-soil core carries no
#' information about a peat pixel: the two differ in the property being mapped,
#' not merely in location. Without gating, the five mineral cores -- which hold
#' MORE carbon in the top 30 cm than the peat cores -- would pull peat pixels
#' upward purely because they are nearby. Gating sets the weight to zero
#' wherever the pixel's stratum and the core's stratum disagree.
#'
#' @param dist_km        n_pixel x n_core distance matrix.
#' @param core_stratum   length n_core character vector.
#' @param pixel_stratum  length n_pixel character vector, or NULL to disable
#'                       gating (in which case a warning-worthy assumption is
#'                       being made and the caller should say so).
build_weight_matrix <- function(dist_km, length_scale_km, max_km,
                                core_stratum = NULL, pixel_stratum = NULL) {
  W <- kernel_truncated(dist_km, length_scale_km, max_km)
  dim(W) <- dim(dist_km)

  if (!is.null(core_stratum) && !is.null(pixel_stratum)) {
    for (j in seq_along(core_stratum)) {
      mismatch <- !is.na(pixel_stratum) & pixel_stratum != core_stratum[j]
      W[mismatch, j] <- 0
    }
  }
  W
}

#' Observation uncertainty for each core, with a floor.
#'
#' The natural choice is the within-stratum standard deviation: it measures how
#' much stock varies between nearby cores of the same kind, which is exactly
#' what limits how far one core's value can be trusted.
#'
#' It needs a floor. With three peat cores the sample sd is 0.37 kg/m2, and
#' taken at face value that would give each peat core roughly forty times the
#' precision of a mineral core and let three points overwrite the published map
#' across their whole neighbourhood. An sd estimated from three observations is
#' not a reliable small number; it is an unreliable one that happens to be
#' small. The floor is the pooled sd across all cores, which cannot be smaller
#' than the true local variability by construction.
core_sigma <- function(stock, stratum, floor_sd = NULL) {
  pooled <- stats::sd(stock[is.finite(stock)])
  if (is.null(floor_sd)) floor_sd <- pooled
  sds <- tapply(stock, stratum, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2L) NA_real_ else stats::sd(v)
  })
  raw <- as.numeric(sds[match(stratum, names(sds))])
  data.frame(
    stratum        = stratum,
    sigma_raw      = raw,
    sigma_used     = pmax(raw, floor_sd, na.rm = TRUE),
    floored        = !is.na(raw) & raw < floor_sd,
    floor_sd       = floor_sd,
    stringsAsFactors = FALSE
  )
}

# ---- depth partitioning for full-column priors -------------------------------

#' Split a full-peat-column stock into a shallow layer and everything beneath.
#'
#' Required whenever a full-column product is to be updated with cores that
#' only reach 30 cm. There is no way to avoid an assumption here, so it is made
#' explicit and isolated in one function.
#'
#' `shallow_fraction` is the share of full-column carbon lying in the top
#' `shallow_cm`. Under uniform carbon density it would be shallow_cm/peat_depth
#' (30/184 = 0.163 for the Hudson Bay Lowlands mean). Real peat compacts with
#' depth -- these cores show bulk density rising more than fourfold between the
#' surface and 30 cm -- so the true fraction is LOWER than uniform. The default
#' should be set from that reasoning and varied in sensitivity analysis.
#'
#' Pure. Returns shallow and deep components in the prior's units.
split_full_column <- function(stock_full, shallow_fraction) {
  if (any(shallow_fraction <= 0 | shallow_fraction >= 1, na.rm = TRUE)) {
    stop("shallow_fraction must lie strictly between 0 and 1", call. = FALSE)
  }
  data.frame(
    shallow = stock_full * shallow_fraction,
    deep    = stock_full * (1 - shallow_fraction),
    stringsAsFactors = FALSE
  )
}

#' Recombine an updated shallow layer with the untouched deep layer.
#'
#' The deep component is carried through unchanged because nothing in these
#' cores observes it. Recombining is therefore addition, and the resulting
#' uncertainty is the shallow posterior's uncertainty plus the prior's
#' uncertainty on the deep part, added in quadrature as independent terms.
recombine_full_column <- function(shallow_mean, shallow_sd, deep_mean, deep_sd) {
  data.frame(
    mean = shallow_mean + deep_mean,
    sd   = sqrt(shallow_sd^2 + deep_sd^2),
    stringsAsFactors = FALSE
  )
}

#' Uniform-density shallow fraction, for reference and for sensitivity bounds.
uniform_shallow_fraction <- function(shallow_cm, peat_depth_cm) {
  pmin(1 - 1e-9, shallow_cm / peat_depth_cm)
}
