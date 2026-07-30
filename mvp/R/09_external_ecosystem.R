# =============================================================================
# 09_external_ecosystem.R      ADD-ON, needs Earth Engine.
#
# Attach an ecosystem class to every external profile AND to our own eight
# cores, so the comparison can be grouped by what kind of ground each sample
# sits on rather than only by which database it came from.
#
# Two class layers, because they disagree in a way that is itself informative:
#
#   GWL_FCS30 (Zhang et al. 2023)  a 30 m WETLAND classification -- swamp,
#       marsh, flooded flat, saline, tidal flat. Preferred here, because ESA
#       WorldCover cannot separate treed bog growing on metres of peat from
#       upland spruce on mineral soil: both read as tree cover. In the Hudson
#       Bay Lowlands that is exactly the distinction that decides carbon.
#
#   ESA WorldCover v200  a 10 m general LAND COVER classification. Reported
#       alongside so the two can be compared; where they disagree, the
#       disagreement is worth looking at.
#
# ~11,500 points, so requests are chunked (CFG$ext$gee_chunk_size). This is
# well within Earth Engine's limits once chunked, but a single call with every
# point would be refused.
#
# INPUT   mvp/outputs/current/external_profiles.csv   (from 08)
#         mvp/outputs/current/cores_clean.csv         (from 01)
# OUTPUT  mvp/outputs/current/ecosystem_classes.csv
#
# REQUIRES: readr, dplyr, sf, rgee
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

library(readr)
library(dplyr)
library(sf)
library(rgee)

msg("09  ecosystem classes at every profile (Earth Engine)")

E <- CFG$ext

ext_path <- file.path(CFG$dir_current, "external_profiles.csv")
if (!file.exists(ext_path)) stop("Missing ", ext_path, " -- run 08 first.")

ext <- read_csv(ext_path, show_col_types = FALSE) %>%
  transmute(dataset, profile_id = as.character(profile_id), latitude, longitude)

ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE) %>%
  transmute(dataset = "Fort Severn", profile_id = core_id, latitude, longitude)

pts <- bind_rows(ours, ext) %>% filter(!is.na(latitude), !is.na(longitude))
msg(nrow(pts), " points to classify (", nrow(ours), " ours + ", nrow(ext), " external)")

ee_Initialize(project = CFG$gee$project)
msg("Earth Engine initialised, project: ", CFG$gee$project)

# ---- the class image -------------------------------------------------------
# Both layers are categorical, so sampling MUST use a mode/first reducer, never
# a mean: the average of class codes 182 and 186 is 184, a different class
# entirely and a completely fabricated one.

gwl <- ee$ImageCollection(E$gee_asset_gwl_fcs30)$mosaic()$rename("gwl_class")
wc  <- ee$ImageCollection(E$gee_asset_worldcover)$first()$
  select("Map")$rename("worldcover")

# UNMASK WITH A SENTINEL BEFORE SAMPLING. sampleRegions() drops a point
# entirely if ANY band is masked there, so a single gap in one class layer
# would silently delete that profile from the comparison rather than return a
# missing class for it. 0 is not a valid code in either layer, so it is
# unambiguous as "no data" and is converted back to NA below.
NODATA <- 0L
class_img <- gwl$addBands(wc)$toInt()$unmask(NODATA)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- chunked sampling ------------------------------------------------------

chunk_size <- E$gee_chunk_size
n_chunks <- ceiling(nrow(pts) / chunk_size)
msg("sampling in ", n_chunks, " chunk(s) of up to ", chunk_size, " points")

sample_chunk <- function(idx) {
  sub <- pts[idx, ]
  fc <- ee$FeatureCollection(lapply(seq_len(nrow(sub)), function(i) {
    ee$Feature(ee$Geometry$Point(c(sub$longitude[i], sub$latitude[i])),
               list(row_key = paste(sub$dataset[i], sub$profile_id[i], sep = "||")))
  }))
  # scale = 30 m matches GWL_FCS30's native grid. WorldCover is 10 m and gets
  # aggregated up; for a categorical layer that means the dominant class, which
  # is why sampling uses first()/mode semantics and never a mean -- the average
  # of class codes 182 and 186 is 184, a different and entirely invented class.
  samp <- class_img$sampleRegions(collection = fc,
                                 properties = list("row_key"),
                                 scale = 30,
                                 geometries = FALSE)
  feats <- samp$getInfo()$features
  if (!length(feats)) return(NULL)
  do.call(rbind, lapply(feats, function(f) {
    p <- f$properties
    data.frame(row_key    = as.character(p$row_key %||% NA),
              gwl_class  = as.integer(p$gwl_class %||% NA),
              worldcover = as.integer(p$worldcover %||% NA),
              stringsAsFactors = FALSE)
  }))
}

results <- list()
for (k in seq_len(n_chunks)) {
  idx <- ((k - 1) * chunk_size + 1):min(k * chunk_size, nrow(pts))
  res <- tryCatch(sample_chunk(idx), error = function(e) {
    msg("  chunk ", k, "/", n_chunks, " FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    results[[length(results) + 1L]] <- res
    msg("  chunk ", k, "/", n_chunks, " ok (", nrow(res), " of ", length(idx), " points)")
  }
}

if (!length(results)) {
  stop("Every chunk failed. Check the GWL_FCS30 asset is readable by project ",
      CFG$gee$project, " -- it is a community asset and may need to be added.")
}

sampled <- bind_rows(results)

# ---- decode to readable labels ---------------------------------------------

gwl_labels <- c("180" = "non-wetland", "181" = "permanent water", "182" = "swamp",
               "183" = "marsh", "184" = "flooded flat", "185" = "saline",
               "186" = "mangrove", "187" = "salt marsh", "188" = "tidal flat")
wc_labels <- c("10" = "tree cover", "20" = "shrubland", "30" = "grassland",
              "40" = "cropland", "50" = "built-up", "60" = "bare/sparse",
              "70" = "snow and ice", "80" = "permanent water",
              "90" = "herbaceous wetland", "95" = "mangrove",
              "100" = "moss and lichen")

out <- sampled %>%
  mutate(
    dataset    = sub("\\|\\|.*$", "", row_key),
    profile_id = sub("^.*\\|\\|", "", row_key),
    # Undo the sampling sentinel: 0 means "no class here", not class zero.
    gwl_class  = ifelse(gwl_class  == NODATA, NA_integer_, gwl_class),
    worldcover = ifelse(worldcover == NODATA, NA_integer_, worldcover),
    gwl_label  = unname(gwl_labels[as.character(gwl_class)]),
    wc_label   = unname(wc_labels[as.character(worldcover)]),
    gwl_is_wetland = gwl_class %in% E$gwl_wetland_codes
  ) %>%
  select(dataset, profile_id, gwl_class, gwl_label, gwl_is_wetland,
        worldcover, wc_label)

msg("classified ", nrow(out), " of ", nrow(pts), " points")
if (nrow(out) < nrow(pts)) {
  msg("  ", nrow(pts) - nrow(out), " point(s) returned nothing -- either a ",
     "failed chunk above, or the point falls outside both class layers.")
}

msg("--- GWL_FCS30 wetland class ---")
print(as.data.frame(count(out, gwl_label, sort = TRUE)))
msg("--- ESA WorldCover ---")
print(as.data.frame(count(out, wc_label, sort = TRUE)))

msg("--- our own cores ---")
print(as.data.frame(out[out$dataset == "Fort Severn", ]))

write_csv(out, file.path(CFG$dir_current, "ecosystem_classes.csv"))
msg("wrote ecosystem_classes.csv")
msg("09 complete")
