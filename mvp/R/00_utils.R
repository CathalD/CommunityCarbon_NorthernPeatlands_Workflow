# =============================================================================
# 00_utils.R  --  the only shared helpers in the MVP. Everything else uses
# packages directly rather than a custom wrapper.
# =============================================================================

#' Locate the mvp/ directory regardless of the caller's working directory.
find_mvp_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  d <- if (length(f)) dirname(f[[1]]) else getwd()
  d <- normalizePath(d, mustWork = FALSE)
  for (i in 1:5) {
    if (file.exists(file.path(d, "config.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not find mvp/config.R walking up from '", getwd(), "'. ",
       "Run scripts from inside the mvp/ directory or its R/ subfolder.",
       call. = FALSE)
}

#' One-line timestamped progress message. Not a logging framework -- just
#' enough to see where a run is without a package dependency for it.
msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' Path to the most recently archived version (e.g. .../versions/carbon_map_v3),
#' or NULL if none exist yet. Version folders are named carbon_map_vN.
latest_version_dir <- function(dir_versions) {
  if (!dir.exists(dir_versions)) return(NULL)
  vdirs <- list.dirs(dir_versions, recursive = FALSE)
  vnums <- suppressWarnings(as.integer(sub(".*carbon_map_v", "", basename(vdirs))))
  ok <- !is.na(vnums)
  if (!any(ok)) return(NULL)
  vdirs[ok][which.max(vnums[ok])]
}

#' Next version number to write, given what's already archived.
next_version_number <- function(dir_versions) {
  latest <- latest_version_dir(dir_versions)
  if (is.null(latest)) return(1L)
  as.integer(sub(".*carbon_map_v", "", basename(latest))) + 1L
}
