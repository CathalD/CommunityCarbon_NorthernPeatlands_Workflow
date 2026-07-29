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
  .args <- commandArgs(trailingOnly = FALSE)
  .this_file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  .this_dir <- if (length(.this_file)) dirname(.this_file[[1]]) else getwd()
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
    latitude          = Latitude,
    longitude_raw     = Longitude,
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
  mutate(stock_kgm2_segment = bulk_density_gcm3 * (soc_pct / 100) * thickness_cm * 10)

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

msg("wrote ", file.path(CFG$dir_current, "cores_clean.geojson"))
msg("wrote ", file.path(CFG$dir_current, "cores_clean.csv"))
msg("01 complete")
