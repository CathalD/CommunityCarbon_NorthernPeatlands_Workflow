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

# Earth Engine rejects a getDownloadURL request above roughly 48 MB. Aim well
# under it: the estimate below is approximate, and a failed request costs a
# round trip.
.GEE_MAX_BYTES <- 32 * 1024^2

#' Estimate the download size of an image over a region. Pure arithmetic.
#'
#' @param area_km2 area of the export region.
#' @param scale_m  metres per pixel.
#' @param n_bands  number of bands.
#' @param bytes_per_px bytes per band per pixel (4 for Float32/Int32).
gee_estimate_bytes <- function(area_km2, scale_m, n_bands, bytes_per_px = 4) {
  n_px <- (area_km2 * 1e6) / (scale_m^2)
  n_px * n_bands * bytes_per_px
}

#' Split band names into groups small enough to download in one request. Pure.
gee_band_chunks <- function(bands, area_km2, scale_m, max_bytes = .GEE_MAX_BYTES,
                            bytes_per_px = 4) {
  per_band <- gee_estimate_bytes(area_km2, scale_m, 1, bytes_per_px)
  if (per_band >= max_bytes) return(NULL)     # even one band is too large
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
                               bytes_per_px = 4) {
  if (!has_pkg("rgee")) stop("rgee is required for gee_download_image()", call. = FALSE)
  if (!dir.exists(dir_local)) dir.create(dir_local, recursive = TRUE,
                                         showWarnings = FALSE)
  dest <- file.path(dir_local, paste0(name, ".tif"))

  if (!overwrite && file.exists(dest)) {
    log_ok("already on disk, skipping: ", basename(dest))
    return(data.frame(name = name, local_path = dest, status = "cached",
                      n_files = 1L, note = gee_verify_local(dest),
                      stringsAsFactors = FALSE))
  }

  bands <- tryCatch(image$bandNames()$getInfo(),
                    error = function(e) NULL)
  if (is.null(bands)) {
    log_warn("could not read band names for ", name, "; skipping")
    return(data.frame(name = name, local_path = NA_character_,
                      status = "failed", n_files = 0L,
                      note = "band names unavailable", stringsAsFactors = FALSE))
  }
  area_km2 <- tryCatch(region$area(maxError = 100)$getInfo() / 1e6,
                       error = function(e) NA_real_)
  if (!is.finite(area_km2)) area_km2 <- 3000   # conservative fallback

  est <- gee_estimate_bytes(area_km2, scale, length(bands), bytes_per_px)
  log_info(sprintf("%s: %d band(s), %.0f km2 at %g m, estimated %.1f MB",
                   name, length(bands), area_km2, scale, est / 1024^2))

  chunks <- gee_band_chunks(bands, area_km2, scale, bytes_per_px = bytes_per_px)
  if (is.null(chunks)) {
    log_warn(strrep("-", 68))
    log_warn("A SINGLE BAND of ", name, " exceeds the download limit at ",
             scale, " m.")
    log_warn(sprintf("  One band is about %.0f MB; the ceiling is %.0f MB.",
                     gee_estimate_bytes(area_km2, scale, 1, bytes_per_px) / 1024^2,
                     .GEE_MAX_BYTES / 1024^2))
    log_warn("  Coarsen the scale, or shrink the area, and re-run. Nothing is")
    log_warn("  written rather than a partial file that looks complete.")
    log_warn(strrep("-", 68))
    return(data.frame(name = name, local_path = NA_character_,
                      status = "too_large", n_files = 0L,
                      note = sprintf("one band ~%.0f MB at %g m",
                                     gee_estimate_bytes(area_km2, scale, 1,
                                                        bytes_per_px) / 1024^2,
                                     scale),
                      stringsAsFactors = FALSE))
  }
  if (length(chunks) > 1L) {
    log_info("splitting into ", length(chunks),
             " band group(s) to stay under the request-size limit")
  }

  paths <- character(0)
  for (i in seq_along(chunks)) {
    bn <- chunks[[i]]
    stem <- if (length(chunks) == 1L) name else sprintf("%s_part%02d", name, i)
    out  <- file.path(dir_local, paste0(stem, ".tif"))
    url <- tryCatch(
      image$select(bn)$getDownloadURL(list(
        name = stem, region = region, scale = scale, crs = crs,
        filePerBand = FALSE, format = "GEO_TIFF")),
      error = function(e) { log_warn("URL request failed for ", stem, ": ",
                                     conditionMessage(e)); NULL })
    if (is.null(url)) next

    ok_dl <- tryCatch({
      utils::download.file(url, destfile = out, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) { log_warn("download failed for ", stem, ": ",
                                      conditionMessage(e)); FALSE })
    if (!ok_dl || !file.exists(out) || file.size(out) < 1024) {
      log_warn("no usable file produced for ", stem)
      if (file.exists(out)) unlink(out)
      next
    }
    log_ok(sprintf("wrote %s (%.1f MB)", basename(out), file.size(out) / 1024^2))
    paths <- c(paths, out)
  }

  if (!length(paths)) {
    return(data.frame(name = name, local_path = NA_character_,
                      status = "failed", n_files = 0L,
                      note = "no file downloaded", stringsAsFactors = FALSE))
  }
  if (length(paths) > 1L) {
    log_warn(name, " arrived as ", length(paths), " band-group files. Stack ",
             "them with terra::rast(c(...)) before use.")
  }

  data.frame(name = name, local_path = paths[[1]],
             status = if (length(paths) > 1L) "local_multipart" else "local",
             n_files = length(paths), note = gee_verify_local(paths[[1]]),
             stringsAsFactors = FALSE)
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
