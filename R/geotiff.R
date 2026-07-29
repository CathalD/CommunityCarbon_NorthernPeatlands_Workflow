# =============================================================================
# R/geotiff.R  --  Minimal GeoTIFF and GeoJSON writers in pure base R.
#
# WHY THIS EXISTS
#
# The scientifically load-bearing parts of this workflow deliberately depend on
# nothing but base R, so that a partner organisation can re-run them on any
# machine. Spatial OUTPUT would normally break that promise, because writing a
# GeoTIFF usually means `terra`, which means GDAL, which means a build
# toolchain. For a community-science pipeline that is a real barrier.
#
# A single-band, uncompressed, strip-organised GeoTIFF is a simple enough
# format to write directly, so that is what this does. The output is a valid
# GeoTIFF that QGIS, ArcGIS, GDAL and terra all read normally.
#
# SCOPE AND LIMITS -- read before extending
#   - Little-endian classic TIFF (not BigTIFF): files must stay under ~4 GB.
#   - Single band, uncompressed, one strip per raster.
#   - Float32 or Int32 samples.
#   - Geographic (EPSG:4326) or any EPSG code written as a projected key. Only
#     4326 is exercised by this workflow and by the tests.
#   - Axis-aligned, north-up grids only (no rotation or shear).
# Anything beyond that should use GDAL rather than be bolted on here.
#
# Reference: TIFF 6.0 specification; OGC GeoTIFF 1.1 (OGC 19-008r4).
# =============================================================================

# ---- low-level TIFF helpers (pure) ------------------------------------------

.tiff_types <- c(BYTE = 1L, ASCII = 2L, SHORT = 3L, LONG = 4L, RATIONAL = 5L,
                 SBYTE = 6L, UNDEFINED = 7L, SSHORT = 8L, SLONG = 9L,
                 SRATIONAL = 10L, FLOAT = 11L, DOUBLE = 12L)

.tiff_type_size <- c("1" = 1L, "2" = 1L, "3" = 2L, "4" = 4L, "5" = 8L,
                     "6" = 1L, "7" = 1L, "8" = 2L, "9" = 4L, "10" = 8L,
                     "11" = 4L, "12" = 8L)

#' Serialise a vector to raw bytes for a given TIFF type. Pure.
.tiff_raw <- function(value, type) {
  switch(as.character(type),
    "2"  = c(charToRaw(paste(value, collapse = "")), as.raw(0L)),
    "3"  = writeBin(as.integer(value), raw(), size = 2L, endian = "little"),
    "4"  = writeBin(as.integer(value), raw(), size = 4L, endian = "little"),
    "9"  = writeBin(as.integer(value), raw(), size = 4L, endian = "little"),
    "11" = writeBin(as.double(value),  raw(), size = 4L, endian = "little"),
    "12" = writeBin(as.double(value),  raw(), size = 8L, endian = "little"),
    stop("unsupported TIFF type: ", type)
  )
}

#' Build the GeoKeyDirectory for a given EPSG code. Pure.
#'
#' Keys must appear in ascending KeyID order; each key is four SHORTs:
#'   (KeyID, TIFFTagLocation, Count, Value_or_Offset).
#' TIFFTagLocation 0 means the value is the Value_or_Offset field itself.
.geo_key_directory <- function(epsg, geographic = TRUE, citation_len = 7L) {
  # 1024 GTModelTypeGeoKey    : 2 = geographic, 1 = projected
  # 1025 GTRasterTypeGeoKey   : 1 = RasterPixelIsArea
  # 2048 GeographicTypeGeoKey / 3072 ProjectedCSTypeGeoKey
  # 2049 GeogCitationGeoKey / 3073 PCSCitationGeoKey : ASCII, held in tag 34737
  #
  # The citation key must declare the true length of the ASCII payload. A fixed
  # count here would make readers run off the end of the string, or truncate a
  # CRS name, so it is passed in rather than assumed.
  cl <- as.integer(citation_len)
  keys <- if (geographic) {
    rbind(c(1024L, 0L, 1L, 2L),
          c(1025L, 0L, 1L, 1L),
          c(2048L, 0L, 1L, as.integer(epsg)),
          c(2049L, 34737L, cl, 0L))
  } else {
    rbind(c(1024L, 0L, 1L, 1L),
          c(1025L, 0L, 1L, 1L),
          c(3072L, 0L, 1L, as.integer(epsg)),
          c(3073L, 34737L, cl, 0L))
  }
  keys <- keys[order(keys[, 1]), , drop = FALSE]
  header <- c(1L, 1L, 0L, nrow(keys))   # version, revision, minor, n keys
  c(header, as.integer(t(keys)))
}

#' Is an EPSG code a geographic (lon/lat) CRS, as far as this workflow needs?
#'
#' Only the codes this workflow actually writes are listed. Anything else is
#' refused rather than guessed at: silently mislabelling a projected raster as
#' geographic would place it in the Gulf of Guinea, and the failure would be
#' obvious; mislabelling one Canadian projection as another would not be.
.epsg_is_geographic <- function(epsg) {
  epsg <- as.integer(epsg)
  if (epsg %in% c(4326L, 4269L, 4617L)) return(TRUE)
  if (epsg %in% c(3978L, 3979L, 3573L, 3347L)) return(FALSE)
  stop("write_geotiff() does not know whether EPSG:", epsg,
       " is geographic or projected. Add it to .epsg_is_geographic().",
       call. = FALSE)
}

#' Write a single-band GeoTIFF from a matrix. Pure with respect to `z`.
#'
#' @param z       matrix; z[1, 1] is the NORTH-WEST corner, rows run north to
#'                south and columns west to east (standard raster orientation).
#' @param path    output file.
#' @param xmin    western edge of the raster (not the centre of the first cell).
#' @param ymax    northern edge of the raster.
#' @param xres,yres  cell size in CRS units; both POSITIVE. For a projected CRS
#'                these are METRES, not degrees.
#' @param epsg    EPSG code. 4326 (lon/lat) and 3978 (NAD83 / Canada Atlas
#'                Lambert, the community deliverable CRS) are both exercised.
#' @param datatype "Float32" or "Int32".
#' @param nodata  value written to the GDAL_NODATA tag. Defaults to NaN for
#'                Float32 and -9999 for Int32.
#' @param citation ASCII CRS name. Defaults to a name matching `epsg` rather
#'                than to "WGS 84", so a projected raster cannot inherit a
#'                geographic citation by forgetting an argument.
#'
#' @return the path, invisibly.
write_geotiff <- function(z, path, xmin, ymax, xres, yres,
                          epsg = 4326L, datatype = c("Float32", "Int32"),
                          nodata = NULL, citation = NULL) {
  datatype <- match.arg(datatype)
  geographic <- .epsg_is_geographic(epsg)
  if (is.null(citation)) {
    citation <- switch(as.character(as.integer(epsg)),
                       "4326" = "WGS 84|",
                       "4269" = "NAD83|",
                       "4617" = "NAD83(CSRS)|",
                       "3978" = "NAD83 / Canada Atlas Lambert|",
                       "3979" = "NAD83(CSRS) / Canada Atlas Lambert|",
                       paste0("EPSG:", epsg, "|"))
  }
  if (!is.matrix(z)) stop("`z` must be a matrix", call. = FALSE)
  if (xres <= 0 || yres <= 0) stop("`xres` and `yres` must be positive", call. = FALSE)

  nrow_z <- nrow(z); ncol_z <- ncol(z)
  is_float <- identical(datatype, "Float32")
  if (is.null(nodata)) nodata <- if (is_float) NaN else -9999L

  sample_format <- if (is_float) 3L else 2L      # 3 = IEEE float, 2 = signed int
  bits <- 32L
  bytes_per_px <- 4L
  strip_bytes <- as.integer(nrow_z) * as.integer(ncol_z) * bytes_per_px

  nodata_str <- if (is.nan(nodata) || (is.numeric(nodata) && is.na(nodata))) {
    "nan"
  } else {
    format(nodata, scientific = FALSE, trim = TRUE)
  }

  geokeys  <- .geo_key_directory(epsg, geographic = geographic,
                                 citation_len = nchar(citation))
  pixscale <- c(xres, yres, 0)
  # Tiepoint maps raster (0,0) -- the OUTER corner of the first cell -- to the
  # world coordinate (xmin, ymax). RasterPixelIsArea, so this is an edge.
  tiepoint <- c(0, 0, 0, xmin, ymax, 0)

  # --- entries, in ascending tag order (a TIFF requirement) -------------------
  ent <- list(
    list(tag = 256L,   type = 4L,  count = 1L, val = ncol_z),        # ImageWidth
    list(tag = 257L,   type = 4L,  count = 1L, val = nrow_z),        # ImageLength
    list(tag = 258L,   type = 3L,  count = 1L, val = bits),          # BitsPerSample
    list(tag = 259L,   type = 3L,  count = 1L, val = 1L),            # Compression: none
    list(tag = 262L,   type = 3L,  count = 1L, val = 1L),            # BlackIsZero
    list(tag = 273L,   type = 4L,  count = 1L, val = NA_integer_),   # StripOffsets (patched)
    list(tag = 277L,   type = 3L,  count = 1L, val = 1L),            # SamplesPerPixel
    list(tag = 278L,   type = 4L,  count = 1L, val = nrow_z),        # RowsPerStrip
    list(tag = 279L,   type = 4L,  count = 1L, val = strip_bytes),   # StripByteCounts
    list(tag = 284L,   type = 3L,  count = 1L, val = 1L),            # PlanarConfig
    list(tag = 339L,   type = 3L,  count = 1L, val = sample_format), # SampleFormat
    list(tag = 33550L, type = 12L, count = 3L, val = pixscale),      # ModelPixelScale
    list(tag = 33922L, type = 12L, count = 6L, val = tiepoint),      # ModelTiepoint
    list(tag = 34735L, type = 3L,  count = length(geokeys), val = geokeys),
    list(tag = 34737L, type = 2L,  count = nchar(citation), val = citation),
    list(tag = 42113L, type = 2L,  count = nchar(nodata_str), val = nodata_str)
  )
  n_ent <- length(ent)

  # --- offset arithmetic ------------------------------------------------------
  ifd_offset  <- 8L
  ifd_bytes   <- 2L + 12L * n_ent + 4L
  data_offset <- ifd_offset + ifd_bytes

  payloads <- vector("list", n_ent)
  offsets  <- integer(n_ent)
  cursor   <- data_offset

  for (i in seq_len(n_ent)) {
    e <- ent[[i]]
    if (e$tag == 273L) next                        # patched once the data area is sized
    sz <- .tiff_type_size[[as.character(e$type)]] * e$count
    if (e$type == 2L) sz <- e$count + 1L           # ASCII carries a trailing NUL
    if (sz > 4L) {
      raw_i <- .tiff_raw(e$val, e$type)
      payloads[[i]] <- raw_i
      offsets[i] <- cursor
      cursor <- cursor + length(raw_i)
      if (cursor %% 2L != 0L) {                    # TIFF wants word alignment
        payloads[[i]] <- c(payloads[[i]], as.raw(0L))
        cursor <- cursor + 1L
      }
    }
  }
  strip_offset <- cursor

  # --- write ------------------------------------------------------------------
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)

  writeBin(charToRaw("II"), con)                                    # little-endian
  writeBin(42L, con, size = 2L, endian = "little")                  # magic
  writeBin(ifd_offset, con, size = 4L, endian = "little")

  writeBin(as.integer(n_ent), con, size = 2L, endian = "little")
  for (i in seq_len(n_ent)) {
    e <- ent[[i]]
    writeBin(as.integer(e$tag),   con, size = 2L, endian = "little")
    writeBin(as.integer(e$type),  con, size = 2L, endian = "little")
    writeBin(as.integer(e$count), con, size = 4L, endian = "little")

    if (e$tag == 273L) {
      writeBin(strip_offset, con, size = 4L, endian = "little")
    } else if (!is.null(payloads[[i]])) {
      writeBin(as.integer(offsets[i]), con, size = 4L, endian = "little")
    } else {
      # Fits in the 4-byte value field. SHORTs occupy the low 2 bytes.
      pad <- .tiff_raw(e$val, e$type)
      pad <- c(pad, raw(4L - length(pad)))
      writeBin(pad[1:4], con)
    }
  }
  writeBin(0L, con, size = 4L, endian = "little")                   # no next IFD

  for (i in seq_len(n_ent)) {
    if (!is.null(payloads[[i]])) writeBin(payloads[[i]], con)
  }

  # Image data: row-major from the north-west corner.
  v <- as.vector(t(z))
  if (is_float) {
    v[is.na(v)] <- NaN
    writeBin(as.double(v), con, size = 4L, endian = "little")
  } else {
    v[is.na(v)] <- as.integer(nodata)
    writeBin(as.integer(v), con, size = 4L, endian = "little")
  }

  invisible(path)
}

# ---- grid construction (pure) ------------------------------------------------

#' Build a regular lon/lat grid definition for an AOI bounding box.
#'
#' Cell size is given in metres and converted to degrees at the AOI's centre
#' latitude, so cells are approximately square on the ground. The grid is
#' returned as edges plus cell-centre coordinate vectors.
make_grid <- function(bbox, cell_m) {
  lat_mid <- (bbox[["ymin"]] + bbox[["ymax"]]) / 2
  deg_lat <- cell_m / 111320                       # metres per degree latitude
  deg_lon <- cell_m / (111320 * cos(lat_mid * pi / 180))

  nx <- as.integer(ceiling((bbox[["xmax"]] - bbox[["xmin"]]) / deg_lon))
  ny <- as.integer(ceiling((bbox[["ymax"]] - bbox[["ymin"]]) / deg_lat))

  list(
    nx = nx, ny = ny,
    xres = deg_lon, yres = deg_lat,
    xmin = bbox[["xmin"]], ymax = bbox[["ymax"]],
    # Cell CENTRES. Rows run north to south to match raster orientation.
    x = bbox[["xmin"]] + (seq_len(nx) - 0.5) * deg_lon,
    y = bbox[["ymax"]] - (seq_len(ny) - 0.5) * deg_lat,
    cell_m = cell_m, lat_mid = lat_mid
  )
}

# ---- GeoJSON (pure) ----------------------------------------------------------

#' Escape a string for JSON. Pure.
.json_str <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

#' Render one R value as a JSON scalar. Pure.
.json_val <- function(v) {
  if (length(v) != 1L || is.na(v)) return("null")
  if (is.logical(v)) return(if (v) "true" else "false")
  if (is.numeric(v)) {
    if (!is.finite(v)) return("null")
    return(format(v, scientific = FALSE, trim = TRUE, digits = 12))
  }
  .json_str(as.character(v))
}

#' Write a point FeatureCollection as GeoJSON. Base R only.
#'
#' @param df    data frame of attributes.
#' @param lon,lat  column names holding coordinates.
#' @param props    attribute columns to carry (defaults to all others).
write_geojson_points <- function(df, path, lon = "longitude", lat = "latitude",
                                 props = NULL) {
  if (is.null(props)) props <- setdiff(names(df), c(lon, lat))
  feats <- vapply(seq_len(nrow(df)), function(i) {
    p <- paste0(vapply(props, function(k)
      paste0(.json_str(k), ": ", .json_val(df[[k]][i])), character(1)),
      collapse = ", ")
    sprintf('    {"type": "Feature", "geometry": {"type": "Point", "coordinates": [%s, %s]}, "properties": {%s}}',
            format(df[[lon]][i], digits = 12, trim = TRUE),
            format(df[[lat]][i], digits = 12, trim = TRUE), p)
  }, character(1))

  writeLines(c(
    '{',
    '  "type": "FeatureCollection",',
    '  "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},',
    '  "features": [',
    paste(feats, collapse = ",\n"),
    '  ]',
    '}'
  ), path)
  invisible(path)
}

#' Write a bounding box as a GeoJSON polygon. Pure.
write_geojson_bbox <- function(bbox, path, props = list()) {
  ring <- sprintf("[[%s, %s], [%s, %s], [%s, %s], [%s, %s], [%s, %s]]",
                  bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymin"]],
                  bbox[["xmax"]], bbox[["ymax"]], bbox[["xmin"]], bbox[["ymax"]],
                  bbox[["xmin"]], bbox[["ymin"]])
  p <- if (length(props)) {
    paste0(vapply(names(props), function(k)
      paste0(.json_str(k), ": ", .json_val(props[[k]])), character(1)),
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
