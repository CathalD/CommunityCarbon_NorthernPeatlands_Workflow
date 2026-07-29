# =============================================================================
# R/warp.R  --  Resample a lon/lat surface onto a projected (EPSG:3978) grid.
#
# WHY THIS EXISTS
#
# The council deliverables are specified in EPSG:3978. Every raster this
# workflow produces is built on a lon/lat grid, because that is the grid Earth
# Engine hands back. Something has to move one to the other.
#
# `terra::project()` is the right tool and is used when terra is installed.
# This module is the fallback and the cross-check. It matters for two reasons:
#
#   1. The GeoTIFF writer in R/geotiff.R is base R, so a machine without GDAL
#      can still produce the deliverable in the CRS the council asked for.
#   2. Reprojection is the step where a silent error is least visible. A map
#      warped with the wrong sign, the wrong parallel, or an off-by-one-cell
#      origin still looks like a map of Hudson Bay. Computing it twice and
#      comparing is the only cheap defence.
#
# METHOD: inverse mapping. For each cell of the TARGET grid, take its centre in
# LCC metres, invert to lon/lat with lcc_to_lonlat(), and sample the source
# there. Inverse mapping is used rather than forward mapping because forward
# mapping leaves holes -- projected source cells do not tile the target grid --
# whereas inverse mapping visits every target cell exactly once by construction.
#
# Pure functions throughout; nothing here touches the filesystem.
# =============================================================================

#' Densify a closed lon/lat ring by linear interpolation between vertices.
#'
#' Needed before projecting an AOI to find its extent. A straight line in
#' lon/lat is a CURVE in Lambert Conformal Conic, so the bounding box of the
#' projected VERTICES is not the bounding box of the projected RING -- the
#' curve bulges outside it. Projecting intermediate points along each edge
#' catches that bulge. With `n_per_edge = 50` the residual error over an edge a
#' few tens of kilometres long is far below one cell.
#'
#' @param lon,lat closed ring (first vertex repeated at the end, or not).
#' @param n_per_edge points to place along each edge, including its start.
#' @return data frame of lon, lat. Pure.
densify_ring <- function(lon, lat, n_per_edge = 50L) {
  stopifnot(length(lon) == length(lat), length(lon) >= 2L)
  n_per_edge <- max(2L, as.integer(n_per_edge))
  # Close the ring if the caller did not.
  if (lon[1] != lon[length(lon)] || lat[1] != lat[length(lat)]) {
    lon <- c(lon, lon[1]); lat <- c(lat, lat[1])
  }
  out_lon <- numeric(0); out_lat <- numeric(0)
  for (i in seq_len(length(lon) - 1L)) {
    f <- seq(0, 1, length.out = n_per_edge)[-n_per_edge]   # drop the shared end
    out_lon <- c(out_lon, lon[i] + f * (lon[i + 1L] - lon[i]))
    out_lat <- c(out_lat, lat[i] + f * (lat[i + 1L] - lat[i]))
  }
  data.frame(lon = c(out_lon, lon[1]), lat = c(out_lat, lat[1]),
             stringsAsFactors = FALSE)
}

#' Define a projected target grid covering a lon/lat ring.
#'
#' The extent is snapped OUTWARD to whole multiples of `cell_m`, so the grid
#' origin is a round number in projected metres and two grids built at the same
#' resolution from different AOIs still align cell-for-cell. That alignment is
#' what lets the prediction, uncertainty and core-contribution rasters be
#' stacked in ArcGIS without resampling.
#'
#' @param lon,lat AOI ring in degrees.
#' @param cell_m  cell size in projected metres.
#' @param p       projection parameters from lcc_params().
#' @param pad_cells cells of margin added on every side.
#' @return a grid definition list. Pure.
lcc_target_grid <- function(lon, lat, cell_m, p = lcc_params(3978),
                            pad_cells = 1L) {
  stopifnot(cell_m > 0)
  ring <- densify_ring(lon, lat)
  xy   <- lonlat_to_lcc(ring$lon, ring$lat, p)

  pad <- pad_cells * cell_m
  xmin <- floor((min(xy$x) - pad) / cell_m) * cell_m
  xmax <- ceiling((max(xy$x) + pad) / cell_m) * cell_m
  ymin <- floor((min(xy$y) - pad) / cell_m) * cell_m
  ymax <- ceiling((max(xy$y) + pad) / cell_m) * cell_m

  nx <- as.integer(round((xmax - xmin) / cell_m))
  ny <- as.integer(round((ymax - ymin) / cell_m))

  list(
    epsg = p$epsg, nx = nx, ny = ny, cell_m = cell_m,
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    # Cell CENTRES. Rows run north to south, matching raster orientation and
    # the row order write_geotiff() expects.
    x = xmin + (seq_len(nx) - 0.5) * cell_m,
    y = ymax - (seq_len(ny) - 0.5) * cell_m
  )
}

#' Sample a lon/lat raster at arbitrary lon/lat points. Pure.
#'
#' @param z    matrix; z[1, 1] is the north-west cell, rows run north to south.
#' @param src  source grid definition: xmin, ymax, xres, yres (as make_grid()).
#' @param lon,lat sample locations.
#' @param method "bilinear" or "nearest".
#' @return numeric vector, NA outside the source extent or where the source is NA.
#'
#' Bilinear is the default for continuous carbon surfaces. Nearest must be used
#' for class rasters -- averaging class code 180 and class code 188 produces
#' 184, which is a different class, and nothing downstream would notice.
sample_lonlat <- function(z, src, lon, lat, method = c("bilinear", "nearest")) {
  method <- match.arg(method)
  nr <- nrow(z); nc <- ncol(z)

  # Continuous cell-index coordinates, measured from the raster's outer edge.
  # Subtracting 0.5 puts integer values at cell CENTRES, which is what the
  # interpolation below assumes. Getting this half-cell wrong shifts the whole
  # map by half a pixel, which is invisible at a glance and wrong everywhere.
  cx <- (lon - src$xmin) / src$xres - 0.5
  cy <- (src$ymax - lat) / src$yres - 0.5

  if (method == "nearest") {
    j <- round(cx) + 1L
    i <- round(cy) + 1L
    ok <- is.finite(i) & is.finite(j) & i >= 1L & i <= nr & j >= 1L & j <= nc
    out <- rep(NA_real_, length(lon))
    out[ok] <- z[cbind(i[ok], j[ok])]
    return(out)
  }

  j0 <- floor(cx); i0 <- floor(cy)
  fx <- cx - j0;   fy <- cy - i0

  # Clamp to the outermost cell centres. Half a cell of the raster lies outside
  # the centre lattice on each side; clamping treats that rim as constant rather
  # than discarding it, so the warped map does not acquire an NA fringe that the
  # source did not have. Anything genuinely off the extent is still NA.
  gv <- function(i, j) {
    ok <- is.finite(i) & is.finite(j)
    ii <- pmin(pmax(i + 1L, 1L), nr)
    jj <- pmin(pmax(j + 1L, 1L), nc)
    v <- rep(NA_real_, length(i))
    v[ok] <- z[cbind(ii[ok], jj[ok])]
    v
  }
  v00 <- gv(i0, j0);      v01 <- gv(i0, j0 + 1L)
  v10 <- gv(i0 + 1L, j0); v11 <- gv(i0 + 1L, j0 + 1L)

  out <- (1 - fy) * ((1 - fx) * v00 + fx * v01) +
         fy       * ((1 - fx) * v10 + fx * v11)

  # A bilinear cell touching NA would otherwise poison four target cells. Fall
  # back to nearest there instead, so coastlines and masked ground keep their
  # edge rather than eroding by one cell per warp.
  bad <- !is.finite(out)
  if (any(bad)) {
    out[bad] <- gv(round(cy[bad]), round(cx[bad]))
  }

  outside <- !is.finite(cx) | !is.finite(cy) |
             cx < -0.5 | cx > (nc - 0.5) | cy < -0.5 | cy > (nr - 0.5)
  out[outside] <- NA_real_
  out
}

#' Warp a lon/lat raster onto a projected LCC grid by inverse mapping. Pure.
#'
#' @param z    source matrix on the lon/lat grid `src`.
#' @param src  source grid definition (xmin, ymax, xres, yres).
#' @param grid target grid from lcc_target_grid().
#' @param p    projection parameters; must match `grid`.
#' @param method passed to sample_lonlat().
#' @return matrix with grid$ny rows and grid$nx columns, north-west first.
warp_to_lcc <- function(z, src, grid, p = lcc_params(3978),
                        method = c("bilinear", "nearest")) {
  method <- match.arg(method)
  if (!identical(as.integer(grid$epsg), as.integer(p$epsg))) {
    stop("grid EPSG (", grid$epsg, ") does not match projection EPSG (",
         p$epsg, "). Refusing to warp into a CRS the grid was not built for.",
         call. = FALSE)
  }
  # Every target cell centre, row-major from the north-west.
  gx <- rep(grid$x, times = grid$ny)
  gy <- rep(grid$y, each  = grid$nx)

  ll <- lcc_to_lonlat(gx, gy, p)
  v  <- sample_lonlat(z, src, ll$lon, ll$lat, method = method)
  matrix(v, nrow = grid$ny, ncol = grid$nx, byrow = TRUE)
}

#' Is a point inside a closed ring? Ray casting, vectorised over points. Pure.
#'
#' Used to clip a warped raster to the AOI, since the target grid is a
#' rectangle in projected space and the AOI is not.
point_in_ring <- function(px, py, ring_x, ring_y) {
  n <- length(ring_x)
  if (ring_x[1] == ring_x[n] && ring_y[1] == ring_y[n]) {
    ring_x <- ring_x[-n]; ring_y <- ring_y[-n]; n <- n - 1L
  }
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    xi <- ring_x[i]; yi <- ring_y[i]
    xj <- ring_x[j]; yj <- ring_y[j]
    # Half-open test on y: a vertex exactly on the ray is counted once, not
    # twice, which is what keeps points level with a vertex from flipping.
    crosses <- (yi > py) != (yj > py)
    xcross  <- xi + (py - yi) * (xj - xi) / (yj - yi)
    hit     <- crosses & (px < xcross)
    inside[hit] <- !inside[hit]
    j <- i
  }
  inside
}

#' Set every cell of a projected grid outside the AOI ring to NA. Pure.
#'
#' The ring is given in lon/lat and projected here, rather than the grid being
#' unprojected, because the ring has a few dozen vertices and the grid has
#' hundreds of thousands of cells.
mask_grid_to_ring <- function(m, grid, ring_lon, ring_lat,
                              p = lcc_params(3978)) {
  ring <- densify_ring(ring_lon, ring_lat)
  rxy  <- lonlat_to_lcc(ring$lon, ring$lat, p)
  gx <- rep(grid$x, times = grid$ny)
  gy <- rep(grid$y, each  = grid$nx)
  keep <- point_in_ring(gx, gy, rxy$x, rxy$y)
  m[!matrix(keep, nrow = grid$ny, ncol = grid$nx, byrow = TRUE)] <- NA_real_
  m
}
