# =============================================================================
# R/gee_io.R  --  Getting Earth Engine results out of Earth Engine.
#
# Every raster this workflow builds in Earth Engine is wanted in two places:
#
#   GOOGLE DRIVE  so it can be shared, opened in Earth Engine again, or handed
#                 to someone without an R session.
#   LOCAL DISK    so the reports can embed it, `terra` can read it, and the
#                 downstream scripts that expect a real prior raster -- notably
#                 11_bayesian_map.R -- can find one.
#
# The second is the one that changes the analysis. 11 falls back to a spatially
# constant regional prior and marks its outputs DEMO_ when no reference raster
# is on disk. Land the GeoTIFFs in outputs/gee/ under the names 11 looks for
# and it silently upgrades to the real published surfaces.
#
# DESIGN NOTE: a failed download must never lose work. The Drive export is
# started first and is allowed to succeed on its own; downloading is a separate,
# recoverable step. If the download fails the task id and filename are reported
# so the file can be fetched by hand, and the script carries on.
#
# These functions require rgee, and terra only for the optional verification
# step. Nothing else in the workflow depends on either.
# =============================================================================

#' Is a package available? Small wrapper so callers read cleanly.
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

#' Export an Earth Engine image to Drive and, optionally, download it locally.
#'
#' @param image      ee$Image to export.
#' @param name       export description AND the local file stem. Keep it stable:
#'                   downstream scripts match on these names.
#' @param region     ee$Geometry to clip to.
#' @param scale      metres per pixel.
#' @param crs        output CRS, e.g. "EPSG:4326".
#' @param folder     Drive folder.
#' @param dir_local  local directory; NULL to skip the download.
#' @param wait       block until the export finishes (needed to download).
#' @param max_minutes give up waiting after this long. The Drive export keeps
#'                   running regardless; only the local copy is abandoned.
#'
#' @return a one-row data frame recording what happened, for the manifest.
gee_export_image <- function(image, name, region, scale, crs,
                             folder = "CCNP_SOC", dir_local = NULL,
                             wait = TRUE, max_minutes = 30,
                             overwrite = TRUE) {
  if (!has_pkg("rgee")) stop("rgee is required for gee_export_image()", call. = FALSE)

  local_path <- if (is.null(dir_local)) NA_character_ else
    file.path(dir_local, paste0(name, ".tif"))

  if (!is.null(dir_local) && !dir.exists(dir_local)) {
    dir.create(dir_local, recursive = TRUE, showWarnings = FALSE)
  }

  # Skip work that has already been done, unless told otherwise.
  if (!overwrite && !is.na(local_path) && file.exists(local_path)) {
    log_ok("already present locally, skipping export: ", basename(local_path))
    return(data.frame(name = name, drive_folder = folder,
                      local_path = local_path, status = "cached",
                      note = "local file already existed",
                      stringsAsFactors = FALSE))
  }

  log_info("exporting ", name, " at ", scale, " m to Drive/", folder)
  task <- rgee::ee_image_to_drive(
    image = image, description = name, folder = folder, region = region,
    scale = scale, crs = crs, maxPixels = 1e13, fileFormat = "GeoTIFF")
  task$start()

  if (!wait) {
    log_ok("started (not waiting): ", name)
    return(data.frame(name = name, drive_folder = folder,
                      local_path = NA_character_, status = "started",
                      note = "wait = FALSE; download separately",
                      stringsAsFactors = FALSE))
  }

  log_info("waiting for Earth Engine to finish ", name,
           " (up to ", max_minutes, " min)")
  waited <- tryCatch({
    rgee::ee_monitoring(task, task_time = 10,
                        max_attempts = ceiling(max_minutes * 6))
    TRUE
  }, error = function(e) {
    log_warn("monitoring stopped for ", name, ": ", conditionMessage(e))
    FALSE
  })

  if (!waited) {
    log_warn("The Drive export for ", name, " may still be running. Check ",
             "rgee::ee_monitoring() and fetch it by hand if needed.")
    return(data.frame(name = name, drive_folder = folder,
                      local_path = NA_character_, status = "drive_pending",
                      note = "export started; wait timed out or errored",
                      stringsAsFactors = FALSE))
  }

  if (is.null(dir_local)) {
    log_ok("exported to Drive: ", name)
    return(data.frame(name = name, drive_folder = folder,
                      local_path = NA_character_, status = "drive_only",
                      note = "no local directory requested",
                      stringsAsFactors = FALSE))
  }

  # The Drive copy exists at this point. Downloading is separately recoverable.
  got <- tryCatch(
    rgee::ee_drive_to_local(task = task, dsn = local_path,
                            overwrite = TRUE, consider = FALSE),
    error = function(e) { log_warn("download failed for ", name, ": ",
                                   conditionMessage(e)); NULL })

  if (is.null(got) || !length(got)) {
    log_warn("Drive copy of ", name, " is fine; the LOCAL copy failed. ",
             "Download it from Drive/", folder, " into ", dir_local,
             " and re-run the downstream script.")
    return(data.frame(name = name, drive_folder = folder,
                      local_path = NA_character_, status = "drive_only",
                      note = "export succeeded; local download failed",
                      stringsAsFactors = FALSE))
  }

  # Earth Engine splits large exports across several files. Record that rather
  # than silently using only the first.
  n_files <- length(got)
  if (n_files > 1L) {
    log_warn(name, " came back as ", n_files, " tiles. `terra` can mosaic ",
             "them with terra::merge(); downstream scripts expecting a single ",
             "file will need that step first.")
  }
  log_ok("wrote ", basename(got[[1]]),
         if (n_files > 1L) paste0(" (+", n_files - 1L, " more tiles)") else "")

  info <- gee_verify_local(got[[1]])
  data.frame(name = name, drive_folder = folder,
             local_path = got[[1]], status = "drive_and_local",
             note = info, stringsAsFactors = FALSE)
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
  nlyr <- terra::nlyr(r)
  msg <- sprintf("%d x %d x %d layer(s), %s, mean %.4g",
                 terra::nrow(r), terra::ncol(r), nlyr,
                 terra::crs(r, describe = TRUE)$code, v)
  if (!is.finite(v)) {
    log_warn(basename(path), " reads as entirely NoData. The export region ",
             "may not overlap the asset, or the asset may be masked there.")
    msg <- paste0(msg, " -- ALL NoData")
  }
  msg
}

#' Export an Earth Engine FeatureCollection to Drive and read it locally.
#'
#' Used for point extractions, which are small enough that the direct route is
#' reliable; Drive is used only so the same artefact is shareable.
gee_export_table <- function(fc, name, dir_local, folder = "CCNP_SOC",
                             max_features = 10000) {
  if (!has_pkg("rgee")) stop("rgee is required", call. = FALSE)
  if (!dir.exists(dir_local)) dir.create(dir_local, recursive = TRUE,
                                         showWarnings = FALSE)
  task <- rgee::ee_table_to_drive(collection = fc, description = name,
                                  folder = folder, fileFormat = "CSV")
  task$start()
  log_ok("started table export to Drive: ", name)

  local <- tryCatch({
    df <- rgee::ee_as_sf(fc, maxFeatures = max_features)
    if (has_pkg("sf")) df <- sf::st_drop_geometry(df)
    p <- file.path(dir_local, paste0(name, ".csv"))
    utils::write.csv(df, p, row.names = FALSE)
    log_ok("wrote ", basename(p), " (", nrow(df), " rows)")
    p
  }, error = function(e) {
    log_warn("local table read failed for ", name, ": ", conditionMessage(e))
    NA_character_
  })
  data.frame(name = name, drive_folder = folder, local_path = local,
             status = if (is.na(local)) "drive_only" else "drive_and_local",
             note = "feature collection", stringsAsFactors = FALSE)
}

#' Write the export manifest, and say plainly what it unlocks.
gee_write_manifest <- function(rows, path) {
  if (!length(rows)) return(invisible(NULL))
  m <- do.call(rbind, Filter(Negate(is.null), rows))
  write_csv_logged(m, path, "what was exported, and where it landed")

  ok_local <- sum(m$status %in% c("drive_and_local", "cached"))
  log_info(ok_local, " of ", nrow(m), " products are now on local disk")
  if (ok_local < nrow(m)) {
    log_warn(nrow(m) - ok_local, " product(s) reached Drive but not local disk. ",
             "They are listed with status 'drive_only' or 'drive_pending'.")
  }
  invisible(m)
}
