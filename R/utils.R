# =============================================================================
# R/utils.R  --  Logging, small pure helpers, and I/O with loud failure.
# Base R only, by design: the scientifically load-bearing steps of this
# workflow must run on any R install without a package-install ordeal.
# =============================================================================

# ---- progress reporting ------------------------------------------------------

.ts <- function() format(Sys.time(), "%H:%M:%S")

log_step <- function(...) {
  message("\n", strrep("=", 74))
  message("[", .ts(), "] ", ...)
  message(strrep("=", 74))
  invisible(NULL)
}

log_info <- function(...) {
  message("  [", .ts(), "] ", ...)
  invisible(NULL)
}

log_ok <- function(...) {
  message("  [", .ts(), "] OK   ", ...)
  invisible(NULL)
}

log_flag <- function(...) {
  message("  [", .ts(), "] FLAG ", ...)
  invisible(NULL)
}

log_warn <- function(...) {
  message("  [", .ts(), "] WARN ", ...)
  invisible(NULL)
}

#' Abort with a clearly framed, actionable message.
#'
#' Used for conditions that would silently corrupt the science if allowed
#' through: geolocation errors, unit mismatches, missing columns.
fail_loudly <- function(title, ..., remedy = NULL) {
  bar <- strrep("!", 74)
  msg <- paste0(
    "\n", bar, "\n",
    "FATAL: ", title, "\n", bar, "\n",
    paste0("  ", c(...), collapse = "\n"),
    if (!is.null(remedy)) paste0("\n\n  REMEDY: ", remedy) else "",
    "\n", bar, "\n"
  )
  stop(msg, call. = FALSE)
}

# ---- pure helpers ------------------------------------------------------------

#' Are two numeric vectors equal to within a relative tolerance?
#' Falls back to absolute comparison where the reference is ~zero.
near_rel <- function(x, y, tol = 1e-6) {
  d <- abs(x - y)
  s <- abs(y)
  ifelse(s < .Machine$double.eps^0.5, d < tol, d / s < tol)
}

#' Trailing integer of an identifier, e.g. "PM-03-A-2" -> 2L.
#' Returns NA_integer_ when the identifier has no trailing integer.
#' Pure.
trailing_int <- function(x) {
  m <- regmatches(x, regexpr("[0-9]+$", x))
  out <- rep(NA_integer_, length(x))
  hit <- regexpr("[0-9]+$", x) > 0
  out[hit] <- as.integer(m)
  out
}

#' Vectorised, NA-safe truth test.
#'
#' Logical columns read from CSV can arrive as character ("TRUE"/"true") or
#' with NAs. `isTRUE()` is scalar-only and returns FALSE for a vector, which
#' silently drops every row when used to subset. This does what the name
#' suggests, elementwise, treating NA as FALSE.
isTRUE_vec <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  !is.na(x) & tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

#' Coefficient of variation (%), NA-safe.
cv_pct <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || mean(x) == 0) return(NA_real_)
  100 * stats::sd(x) / mean(x)
}

#' Standard error of the mean.
se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

#' Great-circle distance in km between lon/lat pairs (haversine). Pure.
haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371.0088
  p <- pi / 180
  dlat <- (lat2 - lat1) * p
  dlon <- (lon2 - lon1) * p
  a <- sin(dlat / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

#' Full pairwise great-circle distance matrix (km) for a set of points. Pure.
pairwise_km <- function(lon, lat) {
  n <- length(lon)
  m <- matrix(0, n, n)
  for (i in seq_len(n)) {
    m[i, ] <- haversine_km(lon[i], lat[i], lon, lat)
  }
  m
}

# ---- I/O ---------------------------------------------------------------------

#' Read a CSV with verbatim headers preserved and no type guessing surprises.
#' check.names = FALSE keeps "Core Id" and "Organic Carbon Density (g/cm^2)"
#' exactly as supplied, which matters for the reconciliation audit.
read_csv_verbatim <- function(path) {
  if (!file.exists(path)) {
    fail_loudly(
      "Input file not found",
      paste0("Expected: ", path),
      remedy = "Check the paths block in config.R."
    )
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NaN"))
}

#' Write a CSV and report it.
write_csv_logged <- function(df, path, what = NULL) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
  log_ok(sprintf("wrote %-52s (%d rows)%s",
                 basename(path), nrow(df),
                 if (is.null(what)) "" else paste0("  <- ", what)))
  invisible(path)
}

#' Load a derived artefact produced by an earlier script, or stop with an
#' instruction naming the script that produces it.
require_artifact <- function(path, produced_by) {
  if (!file.exists(path)) {
    fail_loudly(
      "Required upstream artefact is missing",
      paste0("Missing: ", path),
      paste0("Produced by: ", produced_by),
      remedy = paste0("Run ", produced_by, " first (or run_all.R).")
    )
  }
  read_csv_verbatim(path)
}

#' Fixed-width console table for the run log. Returns the data invisibly.
print_table <- function(df, digits = 4, max_rows = 60L) {
  d <- utils::head(df, max_rows)
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) round(v, digits))
  print(as.data.frame(d), row.names = FALSE)
  if (nrow(df) > max_rows) {
    message("  ... ", nrow(df) - max_rows, " further rows not shown")
  }
  invisible(df)
}
