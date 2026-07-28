# =============================================================================
# R/gee_io.R  --  Getting Earth Engine rasters onto local disk.
#
# LOCAL ONLY. There is no Drive export here and no Drive authentication.
# Rasters are fetched straight from Earth Engine over HTTPS with
# ee$Image$getDownloadURL() and written into outputs/gee/.
#
# Why this rather than an export task: a Drive export means starting a task,
# polling it for minutes, then a second authenticated hop to pull the file
# back. Every one of those steps can fail on its own, and the failure that bit
# this project was in the last of them. A signed download URL is one request,
# needs no Drive scope, and either produces the file or says why not.
#
# THE ONE CONSTRAINT: getDownloadURL has a hard request-size limit (~48 MB).
# A single Float32 band at 30 m over this study area is about 12 MB, so single
# layers are comfortable, but a nineteen-band predictor stack is not. The
# downloader therefore estimates the size up front and splits the request into
# band groups that fit, rather than discovering the limit as an opaque server
# error partway through.
#
# The local copy is what the reports embed, what terra reads, and what
# 11_bayesian_map.R needs in order to drop its DEMO_ fallback and use the real
# published priors.
# =============================================================================

#' Is a package available? Small wrapper so callers read cleanly.
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# Earth Engine rejects a getDownloadURL request above 50331648 bytes. Aim well
# under it: the estimate is approximate and a rejected request costs a round
# trip for nothing.
.GEE_HARD_LIMIT <- 50331648
.GEE_MAX_BYTES  <- 24 * 1024^2

# Observed overhead between a naive pixels x bands x bytes count and the size
# Earth Engine actually reports. Measured at ~1.25 across several requests.
.GEE_OVERHEAD <- 1.30

#' Pixel dimensions Earth Engine will actually produce for a region and scale.
#'
#' TWO THINGS THAT MAKE THE NAIVE AREA CALCULATION WRONG, both learned the hard
#' way when requests five times larger than predicted were refused:
#'
#'  1. getDownloadURL rasterises the region's BOUNDING BOX, not the polygon. A
#'     rotated rectangle costs its envelope. Here that is 6629 km2 for a
#'     2758 km2 AOI -- a factor of 2.4 before anything else.
#'
#'  2. In a geographic CRS, `scale` is metres AT THE EQUATOR converted to a
#'     single degree step used on BOTH axes. At 56 N a 30 m step in longitude
#'     covers only 16.8 m of ground, so the grid is oversampled by 1/cos(lat)
#'     -- another factor of 1.8 here.
#'
#' Pure.
gee_grid_dims <- function(bbox, scale_m, geographic = TRUE) {
  if (geographic) {
    deg <- scale_m / 111319.49      # Earth Engine's equatorial metres-per-degree
    nx <- ceiling((bbox[["xmax"]] - bbox[["xmin"]]) / deg)
    ny <- ceiling((bbox[["ymax"]] - bbox[["ymin"]]) / deg)
  } else {
    nx <- ceiling((bbox[["xmax"]] - bbox[["xmin"]]) / scale_m)
    ny <- ceiling((bbox[["ymax"]] - bbox[["ymin"]]) / scale_m)
  }
  list(nx = as.numeric(nx), ny = as.numeric(ny), n_px = as.numeric(nx) * as.numeric(ny))
}

#' Estimated request size in bytes for a bounding box. Pure.
gee_estimate_bytes <- function(bbox, scale_m, n_bands, bytes_per_px = 4,
                               geographic = TRUE) {
  d <- gee_grid_dims(bbox, scale_m, geographic)
  d$n_px * n_bands * bytes_per_px * .GEE_OVERHEAD
}

#' Split a bounding box into a grid of tiles, each small enough to request.
#'
#' This is the fix that band-chunking could not provide. Chunking by band only
#' helps when a single band already fits; when one band is itself too large the
#' only way through is to cut the AREA. Tiles are returned as bounding boxes to
#' be requested separately and mosaicked afterwards.
#'
#' Pure. Returns a list of bboxes, in row-major order.
gee_tile_bbox <- function(bbox, scale_m, n_bands, bytes_per_px = 4,
                          max_bytes = .GEE_MAX_BYTES, geographic = TRUE,
                          max_tiles = 256L) {
  total <- gee_estimate_bytes(bbox, scale_m, n_bands, bytes_per_px, geographic)
  if (total <= max_bytes) {
    return(list(bbox))
  }
  # Aim for roughly square tiles: split each axis by the same factor.
  need <- ceiling(total / max_bytes)
  side <- ceiling(sqrt(need))

  repeat {
    if (side^2 > max_tiles) return(NULL)      # would need an unreasonable count
    dx <- (bbox[["xmax"]] - bbox[["xmin"]]) / side
    dy <- (bbox[["ymax"]] - bbox[["ymin"]]) / side
    probe <- c(xmin = bbox[["xmin"]], ymin = bbox[["ymin"]],
               xmax = bbox[["xmin"]] + dx, ymax = bbox[["ymin"]] + dy)
    if (gee_estimate_bytes(probe, scale_m, n_bands, bytes_per_px,
                           geographic) <= max_bytes) break
    side <- side + 1L
  }

  out <- list()
  for (j in seq_len(side)) {
    for (i in seq_len(side)) {
      out[[length(out) + 1L]] <- c(
        xmin = bbox[["xmin"]] + (i - 1) * dx,
        ymin = bbox[["ymin"]] + (j - 1) * dy,
        xmax = bbox[["xmin"]] + i * dx,
        ymax = bbox[["ymin"]] + j * dy)
    }
  }
  out
}

#' Split band names into groups small enough to download in one request. Pure.
#'
#' Used together with tiling: tiles cut the area, band groups cut the depth.
gee_band_chunks <- function(bands, bbox, scale_m, max_bytes = .GEE_MAX_BYTES,
                            bytes_per_px = 4, geographic = TRUE) {
  per_band <- gee_estimate_bytes(bbox, scale_m, 1, bytes_per_px, geographic)
  if (per_band >= max_bytes) return(NULL)     # tile first, then call again
  per_chunk <- max(1L, floor(max_bytes / per_band))
  split(bands, ceiling(seq_along(bands) / per_chunk))
}

#' Download an Earth Engine image to a local GeoTIFF. No Drive involved.
#'
#' @param image     ee$Image.
#' @param name      file stem. Downstream scripts match on these names, so keep
#'                  them stable.
#' @param region    ee$Geometry to clip to.
#' @param scale     metres per pixel.
#' @param crs       output CRS. EPSG:4326 keeps these stackable with every
#'                  other product this workflow writes and avoids the
#'                  reprojection that a projected CRS forces on the server.
#' @param dir_local destination directory.
#' @param overwrite re-download a file that already exists.
#'
#' @return a one-row data frame for the manifest.
gee_download_image <- function(image, name, region, scale, dir_local,
                               crs = "EPSG:4326", overwrite = FALSE,
                               bytes_per_px = 4, mosaic = TRUE) {
  if (!has_pkg("rgee")) stop("rgee is required for gee_download_image()", call. = FALSE)
  if (!dir.exists(dir_local)) dir.create(dir_local, recursive = TRUE,
                                         showWarnings = FALSE)
  dest <- file.path(dir_local, paste0(name, ".tif"))

  if (!overwrite && file.exists(dest)) {
    log_ok("already on disk, skipping: ", basename(dest))
    return(data.frame(name = name, local_path = dest, status = "cached",
                      n_files = 1L, n_tiles = NA_integer_,
                      note = gee_verify_local(dest), stringsAsFactors = FALSE))
  }

  bands <- tryCatch(image$bandNames()$getInfo(), error = function(e) NULL)
  if (is.null(bands) || !length(bands)) {
    log_warn("could not read band names for ", name, "; skipping")
    return(data.frame(name = name, local_path = NA_character_,
                      status = "failed", n_files = 0L, n_tiles = 0L,
                      note = "band names unavailable", stringsAsFactors = FALSE))
  }

  # The BOUNDING BOX is what gets rasterised, so that is what must be sized.
  bb <- tryCatch({
    co <- region$bounds(maxError = 10)$coordinates()$getInfo()[[1]]
    xs <- vapply(co, `[`, numeric(1), 1); ys <- vapply(co, `[`, numeric(1), 2)
    c(xmin = min(xs), ymin = min(ys), xmax = max(xs), ymax = max(ys))
  }, error = function(e) NULL)
  if (is.null(bb)) {
    log_warn("could not read region bounds for ", name, "; skipping")
    return(data.frame(name = name, local_path = NA_character_,
                      status = "failed", n_files = 0L, n_tiles = 0L,
                      note = "region bounds unavailable", stringsAsFactors = FALSE))
  }

  geographic <- grepl("4326", crs, fixed = TRUE)
  dims <- gee_grid_dims(bb, scale, geographic)
  est  <- gee_estimate_bytes(bb, scale, length(bands), bytes_per_px, geographic)
  log_info(sprintf("%s: %d band(s), %.0f x %.0f px over the bounding box, est %.0f MB",
                   name, length(bands), dims$nx, dims$ny, est / 1024^2))

  tiles <- gee_tile_bbox(bb, scale, length(bands), bytes_per_px,
                         geographic = geographic)
  if (is.null(tiles)) {
    log_warn(strrep("-", 68))
    log_warn(name, " cannot be tiled into a reasonable number of requests at ",
             scale, " m.")
    log_warn(sprintf("  Estimated total: %.0f MB.", est / 1024^2))
    log_warn("  Coarsen the scale or shrink the area. Nothing is written,")
    log_warn("  rather than a partial file that looks complete.")
    log_warn(strrep("-", 68))
    return(data.frame(name = name, local_path = NA_character_,
                      status = "too_large", n_files = 0L, n_tiles = 0L,
                      note = sprintf("est %.0f MB at %g m", est / 1024^2, scale),
                      stringsAsFactors = FALSE))
  }
  if (length(tiles) > 1L) {
    log_info("tiling into ", length(tiles), " spatial tile(s) to stay under ",
             "the ", round(.GEE_MAX_BYTES / 1024^2), " MB request ceiling")
  }

  paths <- character(0)
  for (ti in seq_along(tiles)) {
    tb <- tiles[[ti]]
    chunks <- gee_band_chunks(bands, tb, scale, bytes_per_px = bytes_per_px,
                              geographic = geographic)
    if (is.null(chunks)) { log_warn("tile ", ti, " still too large; skipping"); next }
    treg <- ee$Geometry$Rectangle(c(tb[["xmin"]], tb[["ymin"]],
                                    tb[["xmax"]], tb[["ymax"]]))
    for (ci in seq_along(chunks)) {
      stem <- name
      if (length(tiles) > 1L)  stem <- sprintf("%s_t%03d", stem, ti)
      if (length(chunks) > 1L) stem <- sprintf("%s_b%02d", stem, ci)
      out <- file.path(dir_local, paste0(stem, ".tif"))

      url <- tryCatch(
        image$select(chunks[[ci]])$getDownloadURL(list(
          name = stem, region = treg, scale = scale, crs = crs,
          filePerBand = FALSE, format = "GEO_TIFF")),
        error = function(e) { log_warn("URL request failed for ", stem, ": ",
                                       conditionMessage(e)); NULL })
      if (is.null(url)) next

      ok_dl <- tryCatch({
        suppressWarnings(utils::download.file(url, destfile = out, mode = "wb",
                                              quiet = TRUE))
        TRUE
      }, error = function(e) { log_warn("download failed for ", stem, ": ",
                                        conditionMessage(e)); FALSE })
      if (!ok_dl || !file.exists(out) || file.size(out) < 1024) {
        log_warn("no usable file produced for ", stem)
        if (file.exists(out)) unlink(out)
        next
      }
      paths <- c(paths, out)
    }
  }

  if (!length(paths)) {
    return(data.frame(name = name, local_path = NA_character_,
                      status = "failed", n_files = 0L, n_tiles = length(tiles),
                      note = "no file downloaded", stringsAsFactors = FALSE))
  }
  log_ok(sprintf("downloaded %d file(s), %.1f MB total", length(paths),
                 sum(file.size(paths)) / 1024^2))

  # Stitch the tiles back into one raster so downstream code sees a single
  # file, which is what every other part of this workflow expects.
  if (length(paths) > 1L && mosaic && has_pkg("terra")) {
    merged <- tryCatch({
      r <- terra::vrt(paths, overwrite = TRUE)
      terra::writeRaster(r, dest, overwrite = TRUE)
      dest
    }, error = function(e) { log_warn("mosaic failed for ", name, ": ",
                                      conditionMessage(e)); NULL })
    if (!is.null(merged)) {
      log_ok("mosaicked ", length(paths), " tiles into ", basename(dest))
      unlink(paths)                       # tiles superseded by the mosaic
      return(data.frame(name = name, local_path = dest, status = "local",
                        n_files = 1L, n_tiles = length(tiles),
                        note = gee_verify_local(dest), stringsAsFactors = FALSE))
    }
  }
  if (length(paths) > 1L) {
    log_warn(name, " is ", length(paths), " tile files; install `terra` to ",
             "have them mosaicked automatically, or use terra::vrt() yourself.")
  }

  data.frame(name = name, local_path = paths[[1]],
             status = if (length(paths) > 1L) "local_multipart" else "local",
             n_files = length(paths), n_tiles = length(tiles),
             note = gee_verify_local(paths[[1]]), stringsAsFactors = FALSE)
}

#' Read a downloaded GeoTIFF's basic properties, to confirm it is usable.
#'
#' Verification rather than assumption: a download can succeed and still leave
#' a file that is empty, in the wrong projection, or entirely NoData.
gee_verify_local <- function(path) {
  if (!has_pkg("terra")) return("not verified (terra not installed)")
  r <- tryCatch(terra::rast(path), error = function(e) NULL)
  if (is.null(r)) return("UNREADABLE by terra")
  v <- tryCatch(terra::global(r[[1]], "mean", na.rm = TRUE)[[1]],
                error = function(e) NA_real_)
  msg <- sprintf("%d x %d x %d layer(s), %s, mean %.4g",
                 terra::nrow(r), terra::ncol(r), terra::nlyr(r),
                 terra::crs(r, describe = TRUE)$code, v)
  if (!is.finite(v)) {
    log_warn(basename(path), " reads as entirely NoData. The region may not ",
             "overlap the asset, or the asset may be masked there.")
    msg <- paste0(msg, " -- ALL NoData")
  }
  msg
}

#' Sample an Earth Engine image at points and return a plain data frame.
#'
#' sampleRegions DROPS a feature entirely when ANY band is masked there, which
#' silently loses cores rather than returning NA for the missing band. The
#' sentinel below keeps every point and is converted back to NA afterwards, so
#' a masked layer costs one value instead of one core.
gee_sample_points <- function(image, points, scale, sentinel = -9999,
                              max_features = 10000) {
  if (!has_pkg("rgee")) stop("rgee is required", call. = FALSE)
  samp <- image$unmask(sentinel, FALSE)$sampleRegions(
    collection = points, scale = scale, geometries = FALSE, tileScale = 4)
  df <- rgee::ee_as_sf(samp, maxFeatures = max_features)
  if (has_pkg("sf")) df <- sf::st_drop_geometry(df)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) {
    x[is.finite(x) & abs(x - sentinel) < 1e-6] <- NA_real_
    x
  })
  df
}

#' Write the download manifest and summarise what landed.
gee_write_manifest <- function(rows, path) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(invisible(NULL))
  m <- do.call(rbind, rows)
  write_csv_logged(m, path, "what was downloaded, and where it landed")

  ok_local <- sum(m$status %in% c("local", "local_multipart", "cached"))
  log_info(ok_local, " of ", nrow(m), " products are on local disk")
  bad <- m[!m$status %in% c("local", "local_multipart", "cached"), ]
  if (nrow(bad)) {
    log_warn(nrow(bad), " product(s) did not download:")
    for (i in seq_len(nrow(bad))) {
      log_warn("   ", bad$name[i], " -- ", bad$status[i], ": ", bad$note[i])
    }
  }
  invisible(m)
}
