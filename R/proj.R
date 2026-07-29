# =============================================================================
# R/proj.R  --  Lambert Conformal Conic projection, in base R.
#
# WHY THIS EXISTS
#
# The community deliverables are specified in EPSG:3978 (NAD83 / Canada Atlas
# Lambert). Reprojecting is normally `sf::st_transform()`, and that is what the
# packaging script uses. This module is the INDEPENDENT CHECK on it.
#
# Eight core locations are the foundation of everything downstream. A silent
# projection error would move them, and unlike the longitude sign error in the
# raw data it would move them plausibly -- a few hundred metres, still in the
# right province, with nothing to notice. So the coordinates are computed twice,
# by two independent routes, and 17_package_for_community.R stops if they
# disagree.
#
# The closed-form ellipsoidal formulas are Snyder (1987), "Map Projections: A
# Working Manual", USGS Professional Paper 1395, pp. 107-109.
#
# Pure functions throughout. Forward projection only; that is all this needs.
# =============================================================================

#' Parameters for the Canadian Lambert Conformal Conic projections.
#'
#' EPSG:3978 and EPSG:3979 share identical projection parameters and differ
#' ONLY in datum realisation (NAD83 versus NAD83(CSRS)), which is under a metre
#' in Canada. Both are listed so the distinction is on the record rather than
#' assumed away -- the two are frequently confused, and the deliverable
#' specifies 3978.
lcc_params <- function(epsg = 3978) {
  base <- list(
    lat_1 = 49, lat_2 = 77,     # standard parallels
    lat_0 = 49, lon_0 = -95,    # origin
    x_0 = 0, y_0 = 0,           # false easting / northing
    a = 6378137,                # GRS80 semi-major axis
    f = 1 / 298.257222101       # GRS80 flattening
  )
  if (epsg == 3978) {
    c(base, list(epsg = 3978L, datum = "NAD83",
                 name = "NAD83 / Canada Atlas Lambert"))
  } else if (epsg == 3979) {
    c(base, list(epsg = 3979L, datum = "NAD83(CSRS)",
                 name = "NAD83(CSRS) / Canada Atlas Lambert"))
  } else {
    stop("lcc_params() knows only EPSG:3978 and EPSG:3979", call. = FALSE)
  }
}

#' Forward Lambert Conformal Conic (2SP), ellipsoidal. Pure and vectorised.
#'
#' @param lon,lat degrees.
#' @param p parameter list from lcc_params().
#' @return data frame of x, y in metres.
lonlat_to_lcc <- function(lon, lat, p = lcc_params(3978)) {
  if (any(abs(lat) >= 90, na.rm = TRUE)) {
    stop("latitude must be strictly between -90 and 90", call. = FALSE)
  }
  d <- pi / 180
  e <- sqrt(2 * p$f - p$f^2)

  # Snyder eq. 14-15: m is the radius of the parallel; t is the isometric
  # latitude term. Both reduce to their spherical forms as e -> 0.
  m <- function(phi) cos(phi) / sqrt(1 - e^2 * sin(phi)^2)
  t <- function(phi) {
    tan(pi / 4 - phi / 2) / ((1 - e * sin(phi)) / (1 + e * sin(phi)))^(e / 2)
  }

  phi1 <- p$lat_1 * d; phi2 <- p$lat_2 * d; phi0 <- p$lat_0 * d
  m1 <- m(phi1); m2 <- m(phi2)
  t1 <- t(phi1); t2 <- t(phi2); t0 <- t(phi0)

  # Cone constant. With two distinct standard parallels this is well defined.
  n <- (log(m1) - log(m2)) / (log(t1) - log(t2))
  bigF <- m1 / (n * t1^n)
  rho0 <- p$a * bigF * t0^n

  phi <- lat * d
  rho <- p$a * bigF * t(phi)^n
  theta <- n * ((lon - p$lon_0) * d)

  data.frame(
    x = p$x_0 + rho * sin(theta),
    y = p$y_0 + rho0 - rho * cos(theta),
    stringsAsFactors = FALSE
  )
}

#' Inverse Lambert Conformal Conic, by iteration on Snyder eq. 7-9. Pure.
#'
#' Present so the forward transform can be round-tripped in tests: a forward
#' formula that is self-consistent but wrong would pass a comparison against
#' itself, so the round trip is checked against the ORIGINAL coordinates.
lcc_to_lonlat <- function(x, y, p = lcc_params(3978), tol = 1e-12,
                          max_iter = 30L) {
  d <- pi / 180
  e <- sqrt(2 * p$f - p$f^2)
  m <- function(phi) cos(phi) / sqrt(1 - e^2 * sin(phi)^2)
  tt <- function(phi) {
    tan(pi / 4 - phi / 2) / ((1 - e * sin(phi)) / (1 + e * sin(phi)))^(e / 2)
  }
  phi1 <- p$lat_1 * d; phi2 <- p$lat_2 * d; phi0 <- p$lat_0 * d
  n <- (log(m(phi1)) - log(m(phi2))) / (log(tt(phi1)) - log(tt(phi2)))
  bigF <- m(phi1) / (n * tt(phi1)^n)
  rho0 <- p$a * bigF * tt(phi0)^n

  xs <- x - p$x_0; ys <- rho0 - (y - p$y_0)
  rho <- sign(n) * sqrt(xs^2 + ys^2)
  theta <- atan2(sign(n) * xs, sign(n) * ys)
  t_val <- (rho / (p$a * bigF))^(1 / n)

  phi <- pi / 2 - 2 * atan(t_val)
  for (i in seq_len(max_iter)) {
    prev <- phi
    phi <- pi / 2 - 2 * atan(t_val * ((1 - e * sin(phi)) /
                                      (1 + e * sin(phi)))^(e / 2))
    if (all(abs(phi - prev) < tol, na.rm = TRUE)) break
  }
  data.frame(lon = p$lon_0 + (theta / n) / d, lat = phi / d,
             stringsAsFactors = FALSE)
}

#' Compare two sets of projected coordinates and summarise the disagreement.
#'
#' Used to check an independent implementation (sf/PROJ) against this one.
#' Returns the per-point distance in metres and the worst case. Pure.
compare_projections <- function(x1, y1, x2, y2, id = NULL) {
  dxy <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  data.frame(
    id = if (is.null(id)) seq_along(dxy) else id,
    x_base_r = x1, y_base_r = y1, x_other = x2, y_other = y2,
    offset_m = dxy, stringsAsFactors = FALSE
  )
}
