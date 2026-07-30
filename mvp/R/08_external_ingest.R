# =============================================================================
# 08_external_ingest.R      ADD-ON. Steps 01-07 do not depend on this.
#
# Harmonise four external soil/sediment databases onto one schema so the Fort
# Severn cores can be compared against them, and derive the numbers behind the
# data-gap narrative.
#
#   CanPeat        1,217 profiles / 37,072 layers   Canadian peatlands
#   NPDB           9,017 profiles / 48,372 layers   AAFC National Pedon Database
#   Janousek       1,284 profiles / 23,018 layers   US Pacific + Gulf tidal wetlands
#   WOSIS Canada      29 profiles /    124 layers   Canadian subset only
#
# This extends the supplied data/reference_scripts/04_combine_all.R, which
# harmonises CanPeat, Janousek and WOSIS but not NPDB. NPDB is the only source
# with a mineral/organic flag and a standardised 0-30 cm stock, so it carries
# most of the mineral-soil comparison and cannot be left out.
#
# TWO STOCK COLUMNS, ALWAYS, AND NEVER INTERCHANGEABLE
#
#   stock_kgm2_total  carbon over whatever depth the profile actually covers.
#                     CanPeat's median profile is 243 cm; ours are 30 cm.
#   stock_kgm2_0_30   carbon over 0-30 cm only, integrated from the layers with
#                     partial layers apportioned by overlap.
#
#   Comparing our 30 cm total against a 243 cm total is a category error, so
#   both columns exist and every figure states which it uses. Profiles that do
#   not reach 30 cm are flagged rather than rescaled.
#
# Profile totals are DERIVED HERE from each layer table rather than read from
# any pre-summarised file, so all four datasets are integrated identically.
#
# INPUT   see CFG$ext in config.R
# OUTPUT  mvp/outputs/current/external_layers.csv
#         mvp/outputs/current/external_profiles.csv
#         mvp/outputs/current/narrative_stats.csv
#
# REQUIRES: readr, dplyr, sf
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

msg("08  external dataset ingest")

E <- CFG$ext
FS_LON <- CFG$site_anchor[["lon"]]
FS_LAT <- CFG$site_anchor[["lat"]]

# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------

#' Great-circle distance to Fort Severn, km.
dist_to_site_km <- function(lon, lat) {
  R <- 6371
  p1 <- FS_LAT * pi / 180; p2 <- lat * pi / 180
  a <- sin((p2 - p1) / 2)^2 +
       cos(p1) * cos(p2) * sin((lon - FS_LON) * pi / 180 / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

#' The standard layer schema every dataset is mapped onto.
LAYER_COLS <- c("dataset", "profile_id", "layer_id", "latitude", "longitude",
               "upper_depth", "lower_depth", "layer_thickness_cm",
               "BDOD", "OrgC_pct", "stock_kgm2_layer", "year")

standardise <- function(df) {
  for (cc in setdiff(LAYER_COLS, names(df))) df[[cc]] <- NA
  df[, LAYER_COLS]
}

# -----------------------------------------------------------------------------
# 1. CanPeat -- raw layer table, needs unit work
# -----------------------------------------------------------------------------
# OrgC and TOTC are already %; BDOD already g/cm3. Layer stock is computed from
# first principles (OrgC% x BDOD x thickness / 10 = kg C/m2) in preference to
# converting the supplied TOTC_Stock from Mg C/ha, because the from-scratch
# route uses the same arithmetic applied to our own cores in step 01. Spot
# checks agree: core 1 layer 1 gives 2.70 kg/m2 either way.

canpeat <- read_csv(E$file_canpeat_layers, show_col_types = FALSE) %>%
  transmute(
    dataset            = "CanPeat",
    profile_id         = as.character(CORE_ID),
    layer_id           = paste(CORE_ID, SAMPLE_NO, sep = "_"),
    latitude           = Latitude,
    longitude          = Longitude,
    upper_depth        = as.numeric(upper_depth),
    lower_depth        = as.numeric(lower_depth),
    layer_thickness_cm = as.numeric(layer_thickness_cm),
    BDOD               = as.numeric(BDOD),
    OrgC_pct           = as.numeric(OrgC),
    stock_kgm2_layer   = case_when(
      !is.na(OrgC_pct) & !is.na(BDOD) & !is.na(layer_thickness_cm) ~
        OrgC_pct * BDOD * layer_thickness_cm / 10,
      !is.na(as.numeric(TOTC_Stock)) ~ as.numeric(TOTC_Stock) * 0.1,  # Mg/ha -> kg/m2
      TRUE ~ NA_real_),
    year               = suppressWarnings(as.integer(Year))
  ) %>% standardise()
msg("CanPeat:       ", nrow(canpeat), " layers, ",
   n_distinct(canpeat$profile_id), " profiles")

# -----------------------------------------------------------------------------
# 2. Janousek -- already harmonised upstream
# -----------------------------------------------------------------------------

janousek <- read_csv(E$file_janousek_layers, show_col_types = FALSE) %>%
  transmute(
    dataset            = "Janousek",
    profile_id         = as.character(profile_id),
    layer_id           = as.character(layer_id),
    latitude, longitude,
    upper_depth        = as.numeric(upper_depth),
    lower_depth        = as.numeric(lower_depth),
    layer_thickness_cm = as.numeric(layer_thickness_cm),
    BDOD               = as.numeric(BDOD),
    OrgC_pct           = as.numeric(OrgC_pct),
    stock_kgm2_layer   = as.numeric(OrgC_Stock_kgm2),
    year               = NA_integer_
  ) %>% standardise()
msg("Janousek:      ", nrow(janousek), " layers, ",
   n_distinct(janousek$profile_id), " profiles")

# -----------------------------------------------------------------------------
# 3. WOSIS -- CANADA SUBSET ONLY
# -----------------------------------------------------------------------------
# OrgC arrives in g/kg, so divide by 10 for %. OrgC_Stock is already kg/m2.
# The global WOSIS table is deliberately excluded: 92% of it is United States
# and it would swamp every comparison with irrelevant geography.

wosis <- read_csv(E$file_wosis_layers, show_col_types = FALSE) %>%
  transmute(
    dataset            = "WOSIS Canada",
    profile_id         = as.character(profile_id),
    layer_id           = as.character(layer_id),
    latitude, longitude,
    upper_depth        = as.numeric(upper_depth),
    lower_depth        = as.numeric(lower_depth),
    layer_thickness_cm = as.numeric(layer_thickness_cm),
    BDOD               = as.numeric(BDOD),
    OrgC_pct           = as.numeric(OrgC) / 10,            # g/kg -> %
    stock_kgm2_layer   = as.numeric(OrgC_Stock),           # already kg/m2
    year               = suppressWarnings(as.integer(substr(as.character(date), 1, 4)))
  ) %>% standardise()
msg("WOSIS Canada:  ", nrow(wosis), " layers, ",
   n_distinct(wosis$profile_id), " profiles")

# -----------------------------------------------------------------------------
# 4. NPDB -- the source 04_combine_all.R was missing
# -----------------------------------------------------------------------------
# stock_Mgha is Mg C/ha; x 0.1 gives kg C/m2. carbon_pct is already %.

npdb_l <- read_csv(E$file_npdb_layers, show_col_types = FALSE) %>%
  transmute(
    dataset            = "NPDB",
    profile_id         = as.character(pedon_id),
    layer_id           = paste(pedon_id, layer_id, sep = "_"),
    latitude, longitude,
    upper_depth        = as.numeric(top_depth_cm),
    lower_depth        = as.numeric(bot_depth_cm),
    layer_thickness_cm = as.numeric(thickness_cm),
    BDOD               = as.numeric(bulk_density),
    OrgC_pct           = as.numeric(carbon_pct),
    stock_kgm2_layer   = as.numeric(stock_Mgha) * 0.1,     # Mg/ha -> kg/m2
    year               = suppressWarnings(as.integer(cal_year))
  ) %>% standardise()
msg("NPDB:          ", nrow(npdb_l), " layers, ",
   n_distinct(npdb_l$profile_id), " profiles")

layers <- bind_rows(canpeat, janousek, wosis, npdb_l) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(upper_depth))

# -----------------------------------------------------------------------------
# 5. Integrate each profile: full column AND 0-30 cm
# -----------------------------------------------------------------------------
# Partial layers straddling 30 cm are apportioned by the fraction of their
# thickness that falls inside the window, which is the same treatment step 01
# applies to our own cores.

ref <- CFG$reference_depth_cm

layers <- layers %>%
  mutate(
    overlap_cm  = pmax(0, pmin(ref, lower_depth) - pmax(0, upper_depth)),
    frac_in_30  = ifelse(is.na(layer_thickness_cm) | layer_thickness_cm <= 0,
                        NA_real_, overlap_cm / layer_thickness_cm),
    stock_in_30 = stock_kgm2_layer * pmin(1, pmax(0, frac_in_30))
  )

profiles <- layers %>%
  group_by(dataset, profile_id) %>%
  summarise(
    latitude          = first(latitude),
    longitude         = first(longitude),
    year              = first(year),
    n_layers          = n(),
    total_depth_cm    = suppressWarnings(max(lower_depth, na.rm = TRUE)),
    stock_kgm2_total  = sum(stock_kgm2_layer, na.rm = TRUE),
    stock_kgm2_0_30   = sum(stock_in_30, na.rm = TRUE),
    surface_OrgC_pct  = OrgC_pct[which.min(upper_depth)][1],
    mean_BDOD         = mean(BDOD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_depth_cm = ifelse(is.finite(total_depth_cm), total_depth_cm, NA_real_),
    # A profile that stopped short of 30 cm gives a LOWER BOUND, not a
    # measurement. Same rule as our own short cores in step 01.
    reaches_30cm   = !is.na(total_depth_cm) & total_depth_cm >= ref - 0.1,
    stock_0_30_is_lower_bound = !reaches_30cm,
    dist_fort_severn_km = dist_to_site_km(longitude, latitude)
  )

# -----------------------------------------------------------------------------
# 6. Organic vs mineral
# -----------------------------------------------------------------------------
# Each dataset's own flag where it has one; the conventional 17% surface
# organic carbon threshold otherwise. Recorded per profile so the basis of the
# classification is never lost.

npdb_organic <- read_csv(E$file_npdb_profiles, show_col_types = FALSE) %>%
  transmute(profile_id = as.character(pedon_id),
           npdb_organic = as.logical(has_organic_horizon))

profiles <- profiles %>%
  left_join(npdb_organic, by = "profile_id") %>%
  mutate(
    soil_type = case_when(
      dataset == "NPDB" & !is.na(npdb_organic) & npdb_organic  ~ "organic",
      dataset == "NPDB" & !is.na(npdb_organic) & !npdb_organic ~ "mineral",
      dataset == "CanPeat"                                     ~ "organic",
      !is.na(surface_OrgC_pct) & surface_OrgC_pct >= E$organic_orgc_pct ~ "organic",
      !is.na(surface_OrgC_pct)                                 ~ "mineral",
      TRUE                                                     ~ "unknown"
    ),
    soil_type_basis = case_when(
      dataset == "NPDB" & !is.na(npdb_organic) ~ "NPDB has_organic_horizon flag",
      dataset == "CanPeat"                     ~ "peat database by construction",
      !is.na(surface_OrgC_pct) ~ sprintf("surface OrgC >= %g%%", E$organic_orgc_pct),
      TRUE                                     ~ "no basis available"
    )
  ) %>%
  select(-npdb_organic)

# -----------------------------------------------------------------------------
# 7. Inside the Hudson & James Bay Lowlands?
# -----------------------------------------------------------------------------

hbl <- st_union(st_geometry(st_read(E$file_hbl, quiet = TRUE)))
prof_sf <- st_as_sf(profiles, coords = c("longitude", "latitude"),
                    crs = CFG$crs_geographic, remove = FALSE)
profiles$in_hbl <- as.vector(apply(st_intersects(prof_sf, hbl, sparse = FALSE), 1, any))

msg(nrow(profiles), " external profiles harmonised; ",
   sum(profiles$in_hbl), " inside the Hudson & James Bay Lowlands")

# -----------------------------------------------------------------------------
# 8. QC
# -----------------------------------------------------------------------------

msg("--- per dataset ---")
qc <- profiles %>%
  group_by(dataset) %>%
  summarise(n = n(),
           med_depth_cm = median(total_depth_cm, na.rm = TRUE),
           med_total    = round(median(stock_kgm2_total, na.rm = TRUE), 2),
           med_0_30     = round(median(stock_kgm2_0_30, na.rm = TRUE), 2),
           reach_30     = sum(reaches_30cm),
           organic      = sum(soil_type == "organic"),
           mineral      = sum(soil_type == "mineral"),
           in_hbl       = sum(in_hbl), .groups = "drop")
print(as.data.frame(qc))

# A 0-30 cm stock far above ~40 kg/m2 is not physically credible for 30 cm of
# soil and usually means a unit slipped somewhere. Surface it rather than let
# it quietly widen every distribution.
implausible <- sum(profiles$stock_kgm2_0_30 > 40, na.rm = TRUE)
if (implausible > 0) {
  msg("WARNING: ", implausible, " profile(s) have a 0-30 cm stock above ",
     "40 kg/m2, which is not physically credible for 30 cm of soil. ",
     "Check units in the source before quoting those distributions.")
}

# -----------------------------------------------------------------------------
# 9. Narrative statistics
# -----------------------------------------------------------------------------
# Written to a file rather than only printed, so the numbers quoted in a report
# are reproducible and traceable to this run.

ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE)

usable <- profiles %>% filter(reaches_30cm, stock_kgm2_0_30 > 0)
usable_min <- usable %>% filter(soil_type == "mineral")
in_hbl <- profiles %>% filter(in_hbl)

ring_counts <- vapply(E$rings_km, function(r) sum(profiles$dist_fort_severn_km <= r), integer(1))

stats <- tibble(
  statistic = c(
    "our cores",
    "our mineral cores",
    "external profiles harmonised",
    "external profiles inside the Hudson & James Bay Lowlands",
    "  of which MINERAL",
    "  of which organic",
    paste0("external profiles within ", E$rings_km, " km of Fort Severn"),
    "nearest external profile of any kind (km)",
    "nearest external profile with a complete 0-30 cm stock (km)",
    "nearest MINERAL external profile with a complete 0-30 cm stock (km)",
    "external profiles with a complete 0-30 cm stock, nationally",
    "  of which mineral"
  ),
  value = c(
    nrow(ours),
    sum(ours$campaign == "mineral"),
    nrow(profiles),
    nrow(in_hbl),
    sum(in_hbl$soil_type == "mineral"),
    sum(in_hbl$soil_type == "organic"),
    as.numeric(ring_counts),
    round(min(profiles$dist_fort_severn_km, na.rm = TRUE), 1),
    if (nrow(usable)) round(min(usable$dist_fort_severn_km), 1) else NA_real_,
    if (nrow(usable_min)) round(min(usable_min$dist_fort_severn_km), 1) else NA_real_,
    nrow(usable),
    nrow(usable_min)
  )
)

msg("--- narrative statistics ---")
print(as.data.frame(stats), right = FALSE)

# -----------------------------------------------------------------------------
# 10. Write
# -----------------------------------------------------------------------------

ensure_dir(CFG$dir_current)
write_csv(layers, file.path(CFG$dir_current, "external_layers.csv"))
write_csv(profiles, file.path(CFG$dir_current, "external_profiles.csv"))
write_csv(stats, file.path(CFG$dir_current, "narrative_stats.csv"))
msg("wrote external_layers.csv (", nrow(layers), " rows), external_profiles.csv (",
   nrow(profiles), " rows), narrative_stats.csv")
msg("08 complete  --  next: 09 for ecosystem classes (needs GEE), or 10 for outputs")
