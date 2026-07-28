# =============================================================================
# R/aoi.R  --  Oriented (rotated) area-of-interest geometry.
#
# An axis-aligned bounding box is a poor fit for a coastal study area. The
# Hudson Bay shoreline at Fort Severn runs north-west to south-east, and the
# community sampling follows it, so a north-up rectangle either clips the ends
# of the sampled transect or wastes most of its area over open water and
# far inland ground that was never visited.
#
# This module builds a rectangle ALIGNED TO THE SAMPLING AXIS instead, with a
# generous buffer, and provides the point-in-rectangle test needed to burn it
# into a raster mask.
#
# All functions are pure. Geometry is done in a local equirectangular
# projection about the AOI centroid: over a study area tens of kilometres
# across at 56 N the distortion is far below the precision this AOI needs, and
# it avoids a projection dependency.
# =============================================================================

#' Local projection to metres about an origin, and its inverse. Pure.
.to_local <- function(lon, lat, lon0, lat0) {
  list(x = (lon - lon0) * 111320 * cos(lat0 * pi / 180),
       y = (lat - lat0) * 111320)
}
.from_local <- function(x, y, lon0, lat0) {
  list(lon = lon0 + x / (111320 * cos(lat0 * pi / 180)),
       lat = lat0 + y / 111320)
}

#' Principal axis of a set of points, as a compass bearing in [0, 180).
#'
#' Found by principal component analysis on the locally projected coordinates:
#' the first component is the direction along which the points are most spread
#' out. For a transect that follows a shoreline, that direction IS the
#' shoreline's local trend, so the orientation is derived from the sampling
#' rather than eyeballed off a map.
#'
#' Returns the bearing plus the share of variance it explains, which is the
#' honest measure of whether the points are elongated enough for an oriented
#' rectangle to be meaningful. Below roughly 60 % the arrangement is closer to
#' a blob than a line and an axis-aligned box is the more honest choice.
principal_axis_bearing <- function(lon, lat) {
  lon0 <- mean(lon); lat0 <- mean(lat)
  p <- .to_local(lon, lat, lon0, lat0)
  pc <- stats::prcomp(cbind(p$x, p$y), center = TRUE)
  v <- pc$rotation[, 1]
  bearing <- (90 - atan2(v[2], v[1]) * 180 / pi) %% 180
  scores <- pc$x
  list(
    bearing_deg      = bearing,
    var_explained    = pc$sdev[1]^2 / sum(pc$sdev^2),
    extent_along_km  = diff(range(scores[, 1])) / 1000,
    extent_across_km = diff(range(scores[, 2])) / 1000,
    lon0 = lon0, lat0 = lat0
  )
}

#' Build a rectangle oriented along a bearing, sized to contain a set of points
#' plus a buffer on each side.
#'
#' @param bearing_deg     compass bearing of the long axis, degrees from north.
#' @param along_buffer_km buffer added beyond the points at EACH end.
#' @param across_buffer_km buffer added beyond the points on EACH side.
#'
#' Returns corner coordinates in lon/lat, closed (first vertex repeated), plus
#' the parameters needed to test membership later. Pure.
oriented_rectangle <- function(lon, lat, bearing_deg,
                               along_buffer_km, across_buffer_km) {
  lon0 <- mean(lon); lat0 <- mean(lat)
  p <- .to_local(lon, lat, lon0, lat0)

  # Rotate into the frame where the long axis is +u.
  th <- (90 - bearing_deg) * pi / 180
  u <-  p$x * cos(th) + p$y * sin(th)
  v <- -p$x * sin(th) + p$y * cos(th)

  u_lo <- min(u) - along_buffer_km * 1000
  u_hi <- max(u) + along_buffer_km * 1000
  v_lo <- min(v) - across_buffer_km * 1000
  v_hi <- max(v) + across_buffer_km * 1000

  cu <- c(u_lo, u_hi, u_hi, u_lo, u_lo)
  cv <- c(v_lo, v_lo, v_hi, v_hi, v_lo)
  cx <-  cu * cos(th) - cv * sin(th)
  cy <-  cu * sin(th) + cv * cos(th)
  g  <- .from_local(cx, cy, lon0, lat0)

  list(
    lon = g$lon, lat = g$lat,
    bearing_deg = bearing_deg,
    length_km = (u_hi - u_lo) / 1000,
    width_km  = (v_hi - v_lo) / 1000,
    area_km2  = ((u_hi - u_lo) / 1000) * ((v_hi - v_lo) / 1000),
    u_lo = u_lo, u_hi = u_hi, v_lo = v_lo, v_hi = v_hi,
    lon0 = lon0, lat0 = lat0, theta = th,
    # Axis-aligned envelope, for callers (such as Earth Engine exports) that
    # need a plain rectangle. Always at least as large as the oriented one.
    bbox = c(xmin = min(g$lon), ymin = min(g$lat),
             xmax = max(g$lon), ymax = max(g$lat))
  )
}

#' Test whether points fall inside an oriented rectangle. Pure and vectorised.
point_in_oriented_rect <- function(lon, lat, rect) {
  p <- .to_local(lon, lat, rect$lon0, rect$lat0)
  u <-  p$x * cos(rect$theta) + p$y * sin(rect$theta)
  v <- -p$x * sin(rect$theta) + p$y * cos(rect$theta)
  u >= rect$u_lo & u <= rect$u_hi & v >= rect$v_lo & v <= rect$v_hi
}

#' Distance along and across the rectangle's axes, in km. Useful as a
#' coast-parallel / coast-perpendicular coordinate system for interpreting
#' gradients: in an isostatically rebounding landscape, distance ACROSS the
#' shoreline trend is a proxy for time since marine emergence. Pure.
oriented_coordinates <- function(lon, lat, rect) {
  p <- .to_local(lon, lat, rect$lon0, rect$lat0)
  data.frame(
    along_km  = ( p$x * cos(rect$theta) + p$y * sin(rect$theta)) / 1000,
    across_km = (-p$x * sin(rect$theta) + p$y * cos(rect$theta)) / 1000,
    stringsAsFactors = FALSE
  )
}

#' Write a closed ring of lon/lat as a GeoJSON polygon. Base R only.
write_geojson_polygon <- function(lon, lat, path, props = list()) {
  if (lon[1] != lon[length(lon)] || lat[1] != lat[length(lat)]) {
    lon <- c(lon, lon[1]); lat <- c(lat, lat[1])
  }
  ring <- paste0("[", paste(sprintf("[%.8f, %.8f]", lon, lat), collapse = ", "), "]")
  p <- if (length(props)) {
    paste0(vapply(names(props), function(k)
      paste0("\"", k, "\": ", .json_val(props[[k]])), character(1)),
      collapse = ", ")
  } else ""
  writeLines(c(
    '{', '  "type": "FeatureCollection",',
    '  "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},',
    '  "features": [',
    sprintf('    {"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [%s]}, "properties": {%s}}',
            ring, p),
    '  ]', '}'
  ), path)
  invisible(path)
}

#' Render an oriented rectangle as an Earth Engine geometry literal, so the AOI
#' used locally and the AOI used in Earth Engine are provably the same shape.
ee_polygon_snippet <- function(rect, var_name = "aoi") {
  coords <- paste0("[", paste(sprintf("[%.6f, %.6f]", rect$lon, rect$lat),
                              collapse = ",\n     "), "]")
  paste0("var ", var_name, " = ee.Geometry.Polygon([\n    ", coords, "\n]);")
}
