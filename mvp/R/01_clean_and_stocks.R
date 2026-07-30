# =============================================================================
# 01_clean_and_stocks.R
#
# Read the raw core CSV, fix the one known data error (longitude sign),
# flag (never drop) implausible bulk density values, and compute per-core
# carbon stock.
#
# No Earth Engine needed. Run this one first.
#
# INPUT   data/raw/community_soil_cores.csv          (repo root, shared with
#                                                      the original workflow)
# OUTPUT  mvp/outputs/current/cores_clean.geojson     one point per core
#         mvp/outputs/current/cores_clean.csv         same data, flat table
#
# REQUIRES: readr, dplyr, sf
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  # Works three ways: Rscript (finds --file=), source() from ANY working
  # directory (finds the sourced file's own path), and run_all.R (sets
  # MVP_R_DIR). Only falls back to getwd() if all three fail.
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

library(readr)
library(dplyr)
library(sf)

msg("01  clean & compute stocks")

# ---- 1. read --------------------------------------------------------------

raw <- read_csv(CFG$file_cores_raw, show_col_types = FALSE)
msg("read ", nrow(raw), " segment rows from ", basename(CFG$file_cores_raw))

seg <- raw %>%
  transmute(
    core_id           = trimws(`Core Id`),
    sample_id         = trimws(`Sample Id`),
    latitude          = Latitude,
    longitude_raw     = Longitude,
    # NOTE: `Depth` in the raw file is segment THICKNESS, not depth below
    # surface. Read as depth-below-surface the values would decrease going
    # down two of the cores, which is impossible -- that is how the original
    # workflow established the semantics, and it holds here.
    thickness_cm      = Depth,
    bulk_density_gcm3 = `Bulk Density`,
    soc_pct           = SOC
  ) %>%
  mutate(
    # The raw file stores western longitude as POSITIVE, which places every
    # core in Siberia. This is a known, one-directional error -- fix it and
    # flag every row it touched.
    lon_repaired = longitude_raw > 0,
    longitude    = if_else(lon_repaired, -longitude_raw, longitude_raw),
    campaign     = case_when(
      grepl("^PM", core_id) ~ "peat",
      grepl("^FS", core_id) ~ "mineral",
      TRUE ~ NA_character_
    ),
    bd_flagged   = bulk_density_gcm3 < CFG$qc$bd_min |
                   bulk_density_gcm3 > CFG$qc$bd_max
  )

if (any(is.na(seg$campaign))) {
  stop("Core id(s) matched neither the 'PM' nor 'FS' prefix: ",
       paste(unique(seg$core_id[is.na(seg$campaign)]), collapse = ", "))
}

n_flagged <- sum(seg$bd_flagged)
if (n_flagged > 0) {
  msg("flagged (not dropped) ", n_flagged, " segment(s) with implausible bulk density")
}
msg("repaired longitude sign on ", sum(seg$lon_repaired), " of ", nrow(seg), " rows")

# ---- 2. carbon stock per segment, then per core ----------------------------
#
#   stock (kg C / m2) = bulk_density (g/cm3) x SOC (fraction) x thickness (cm) x 10
#
# The x10 converts g/cm2 (bulk_density x thickness) x fraction to kg/m2.

seg <- seg %>%
  mutate(stock_kgm2_segment = bulk_density_gcm3 * (soc_pct / 100) * thickness_cm * 10,
        # Carbon per centimetre of profile. This, not the segment total, is
        # what makes profiles with unequal segment thickness comparable, so it
        # is what the raw profile plot shows.
        carbon_density_kgm2_per_cm = bulk_density_gcm3 * (soc_pct / 100) * 10)

# ---- 2b. depth intervals, from the segment order --------------------------
# Segment order comes from the trailing integer in Sample Id (PM-03-A-1 ..
# PM-03-A-4). Single-segment cores carry no trailing integer and get index 1.
# Depth intervals are then the running sum of thickness within each core.
seg <- seg %>%
  mutate(seg_index = suppressWarnings(as.integer(sub(".*-(\\d+)$", "\\1", sample_id))),
        seg_index = ifelse(is.na(seg_index), 1L, seg_index)) %>%
  arrange(core_id, seg_index) %>%
  group_by(core_id) %>%
  mutate(depth_bottom_cm = cumsum(thickness_cm),
        depth_top_cm    = depth_bottom_cm - thickness_cm) %>%
  ungroup()

cores <- seg %>%
  group_by(core_id, campaign) %>%
  summarise(
    latitude    = first(latitude),
    longitude   = first(longitude),
    lon_repaired = first(lon_repaired),
    core_depth_cm = sum(thickness_cm),
    stock_kgm2    = sum(stock_kgm2_segment),
    any_bd_flagged = any(bd_flagged),
    n_segments     = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # A core that never reached the reference depth is a LOWER BOUND, not a
    # complete measurement. Kept in the dataset, never rescaled to look
    # complete -- just flagged so downstream steps (and you) know.
    stock_is_lower_bound = core_depth_cm < (CFG$reference_depth_cm - 0.1)
  )

msg(nrow(cores), " cores: ", sum(cores$campaign == "peat"), " peat, ",
   sum(cores$campaign == "mineral"), " mineral")
if (any(cores$stock_is_lower_bound)) {
  msg("lower-bound cores (core depth < ", CFG$reference_depth_cm, " cm): ",
     paste(cores$core_id[cores$stock_is_lower_bound], collapse = ", "))
}
print(cores %>% select(core_id, campaign, core_depth_cm, stock_kgm2,
                       stock_is_lower_bound), n = Inf)

# ---- 3. write ---------------------------------------------------------------

ensure_dir(CFG$dir_current)

cores_sf <- st_as_sf(cores, coords = c("longitude", "latitude"),
                     crs = CFG$crs_geographic, remove = FALSE)

st_write(cores_sf, file.path(CFG$dir_current, "cores_clean.geojson"),
         delete_dsn = TRUE, quiet = TRUE)
write_csv(cores, file.path(CFG$dir_current, "cores_clean.csv"))
# Segment level, with depth intervals -- the input for the profile plots.
write_csv(seg, file.path(CFG$dir_current, "segments_clean.csv"))

msg("wrote cores_clean.geojson, cores_clean.csv, segments_clean.csv")
msg("01 complete  --  next: 01b_plot_profiles.R to eyeball the profiles")
