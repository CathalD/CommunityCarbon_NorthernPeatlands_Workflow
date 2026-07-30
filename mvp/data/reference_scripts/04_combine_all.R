# =============================================================================
# SCRIPT 4: Combine All Three Databases — WOSIS, CANPEAT, Janousek
# =============================================================================
# Applies unit conversions to WOSIS and CANPEAT outputs, standardizes
# column names across all three, and row-binds into combined outputs.
#
# Prerequisites: Run scripts 01, 02, and 03 first to generate all 6 inputs.
#
# Unit conversions applied:
#   WOSIS:   OrgC g/kg → % (÷ 10), TOTC g/kg → % (÷ 10)
#            OrgC_Stock already in kg/m² — no change
#            mean_BDOD and mean_OrgC_pct absent from wosis_profiles.csv —
#            computed here from wosis_layers.csv
#   CANPEAT: TOTC_Stock Mg C/ha → kg/m² (× 0.1) [layer-level]
#            Total_C_MgCha Mg C/ha → kg/m² (× 0.1) [profile-level]
#            OrgC/TOTC already in % — no change
#            BDOD already in g/cm³ — no change
#            country_name: hardcoded "Canada"
#   JANOUSEK: all units already harmonized in script 03 — no change
#
# Outputs:
#   - combined_layers.csv   : one row per layer across all three databases
#   - combined_profiles.csv : one row per profile/core across all three databases
#
# Combined layer columns:
#   dataset, profile_id, layer_id, latitude, longitude,
#   upper_depth, lower_depth, layer_thickness_cm,
#   BDOD, OrgC_pct, TOTC_pct, OrgC_Stock_kgm2,
#   country_name, year
#
# Combined profile columns:
#   dataset, profile_id, latitude, longitude,
#   total_depth_cm, n_layers, sum_OrgC_Stock_kgm2,
#   mean_BDOD, mean_OrgC_pct,
#   country_name, year
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 0. USER SETTINGS — update paths as needed
# -----------------------------------------------------------------------------

PATH_WOSIS_LAYERS    <- "output/wosis_layers.csv"
PATH_WOSIS_PROFILES  <- "output/wosis_profiles.csv"
PATH_PEAT_LAYERS     <- "output/peat_layers.csv"
PATH_PEAT_PROFILES   <- "output/peat_profiles_gee.csv"
PATH_JAN_LAYERS      <- "output/janousek_layers.csv"
PATH_JAN_PROFILES    <- "output/janousek_profiles.csv"

PATH_OUT_LAYERS      <- "output/combined_layers.csv"
PATH_OUT_PROFILES    <- "output/combined_profiles.csv"

# -----------------------------------------------------------------------------
# 1. LOAD ALL SIX INPUT FILES
# -----------------------------------------------------------------------------
message("--- Loading all input files ---")

wosis_l  <- read_csv(PATH_WOSIS_LAYERS,   show_col_types = FALSE)
wosis_p  <- read_csv(PATH_WOSIS_PROFILES,  show_col_types = FALSE)
peat_l   <- read_csv(PATH_PEAT_LAYERS,     show_col_types = FALSE)
peat_p   <- read_csv(PATH_PEAT_PROFILES,   show_col_types = FALSE)
jan_l    <- read_csv(PATH_JAN_LAYERS,      show_col_types = FALSE)
jan_p    <- read_csv(PATH_JAN_PROFILES,    show_col_types = FALSE)

message(sprintf("WOSIS layers:     %d rows", nrow(wosis_l)))
message(sprintf("WOSIS profiles:   %d rows", nrow(wosis_p)))
message(sprintf("CANPEAT layers:   %d rows", nrow(peat_l)))
message(sprintf("CANPEAT profiles: %d rows", nrow(peat_p)))
message(sprintf("Janousek layers:  %d rows", nrow(jan_l)))
message(sprintf("Janousek profiles:%d rows", nrow(jan_p)))

# -----------------------------------------------------------------------------
# 2. HARMONIZE WOSIS LAYERS
# -----------------------------------------------------------------------------
message("\n--- Harmonizing WOSIS layers ---")

# OrgC/TOTC in g/kg → divide by 10 to get %
# OrgC_Stock already in kg/m²
# date is "YYYY-MM-D" string → extract 4-char year prefix

wosis_l_harm <- wosis_l %>%
  mutate(
    profile_id      = as.character(profile_id),
    layer_id        = as.character(layer_id),
    OrgC_pct        = OrgC / 10,           # g/kg → %
    TOTC_pct        = TOTC / 10,           # g/kg → %
    OrgC_Stock_kgm2 = OrgC_Stock,          # already kg/m²
    year            = as.integer(substr(as.character(date), 1, 4))
  ) %>%
  select(
    dataset,
    profile_id,
    layer_id,
    latitude,
    longitude,
    upper_depth,
    lower_depth,
    layer_thickness_cm,
    BDOD,
    OrgC_pct,
    TOTC_pct,
    OrgC_Stock_kgm2,
    country_name,
    year
  )

cat("\nWOSIS OrgC_pct summary (should be 0-13%, not 0-130 g/kg):\n")
print(summary(wosis_l_harm$OrgC_pct))
cat("WOSIS OrgC_Stock_kgm2 summary:\n")
print(summary(wosis_l_harm$OrgC_Stock_kgm2))

# -----------------------------------------------------------------------------
# 3. HARMONIZE CANPEAT LAYERS
# -----------------------------------------------------------------------------
message("\n--- Harmonizing CANPEAT layers ---")

# OrgC/TOTC already in %; BDOD already in g/cm³
# TOTC_Stock is Mg C/ha; recalculate OrgC_Stock from scratch where possible,
# fall back to TOTC_Stock × 0.1 (Mg C/ha → kg/m²)
# country_name: hardcoded "Canada" (all CANPEAT cores are Canadian)

peat_l_harm <- peat_l %>%
  mutate(
    profile_id      = as.character(CORE_ID),
    layer_id        = paste(CORE_ID, SAMPLE_NO, sep = "_"),
    latitude        = Latitude,
    longitude       = Longitude,
    OrgC_pct        = OrgC,                # already %
    TOTC_pct        = TOTC,                # already %
    # Prefer from-scratch calculation; fall back to pre-calculated conversion
    OrgC_Stock_kgm2 = case_when(
      !is.na(OrgC) & !is.na(BDOD) & !is.na(layer_thickness_cm) ~
        OrgC * BDOD * layer_thickness_cm / 10,
      !is.na(TOTC_Stock) ~ TOTC_Stock * 0.1,   # Mg C/ha → kg/m²
      TRUE ~ NA_real_
    ),
    country_name    = "Canada",
    year            = as.integer(Year)
  ) %>%
  select(
    dataset,
    profile_id,
    layer_id,
    latitude,
    longitude,
    upper_depth,
    lower_depth,
    layer_thickness_cm,
    BDOD,
    OrgC_pct,
    TOTC_pct,
    OrgC_Stock_kgm2,
    country_name,
    year
  )

cat("\nCANPEAT OrgC_pct summary (should be in %):\n")
print(summary(peat_l_harm$OrgC_pct))
cat("CANPEAT OrgC_Stock_kgm2 summary:\n")
print(summary(peat_l_harm$OrgC_Stock_kgm2))

# -----------------------------------------------------------------------------
# 4. HARMONIZE JANOUSEK LAYERS
# -----------------------------------------------------------------------------
message("\n--- Harmonizing Janousek layers ---")

# Script 03 already produced harmonized column names and units
# Just enforce character types and ensure NA columns are present

jan_l_harm <- jan_l %>%
  mutate(
    profile_id   = as.character(profile_id),
    layer_id     = as.character(layer_id),
    country_name = NA_character_,
    year         = NA_integer_
  ) %>%
  select(
    dataset,
    profile_id,
    layer_id,
    latitude,
    longitude,
    upper_depth,
    lower_depth,
    layer_thickness_cm,
    BDOD,
    OrgC_pct,
    TOTC_pct,
    OrgC_Stock_kgm2,
    country_name,
    year
  )

cat("\nJanousek OrgC_pct summary:\n")
print(summary(jan_l_harm$OrgC_pct))
cat("Janousek OrgC_Stock_kgm2 summary:\n")
print(summary(jan_l_harm$OrgC_Stock_kgm2))

# -----------------------------------------------------------------------------
# 5. ROW-BIND LAYERS + QC + EXPORT
# -----------------------------------------------------------------------------
message("\n--- Combining all layer tables ---")

combined_layers <- bind_rows(wosis_l_harm, peat_l_harm, jan_l_harm)

message(sprintf("Combined layer table: %d rows", nrow(combined_layers)))
message(sprintf("  WOSIS:    %d", nrow(wosis_l_harm)))
message(sprintf("  CANPEAT:  %d", nrow(peat_l_harm)))
message(sprintf("  Janousek: %d", nrow(jan_l_harm)))

# QC: verify row count matches sum of inputs
expected_rows <- nrow(wosis_l_harm) + nrow(peat_l_harm) + nrow(jan_l_harm)
if (nrow(combined_layers) != expected_rows) {
  warning(sprintf("Row count mismatch: expected %d, got %d",
                  expected_rows, nrow(combined_layers)))
} else {
  message("OK: Combined row count matches sum of inputs.")
}

# QC: row counts by dataset
cat("\nLayer rows per dataset:\n")
print(count(combined_layers, dataset))

# QC: NA rates per key column by dataset
cat("\nNA counts per key column by dataset:\n")
combined_layers %>%
  group_by(dataset) %>%
  summarise(
    n_rows        = n(),
    na_BDOD       = sum(is.na(BDOD)),
    na_OrgC_pct   = sum(is.na(OrgC_pct)),
    na_TOTC_pct   = sum(is.na(TOTC_pct)),
    na_Stock      = sum(is.na(OrgC_Stock_kgm2)),
    na_country    = sum(is.na(country_name)),
    na_year       = sum(is.na(year)),
    .groups = "drop"
  ) %>%
  print()

# QC: unit sanity — median OrgC_pct should be comparable across datasets (not 10× off)
cat("\nUnit sanity check — median OrgC_pct and OrgC_Stock by dataset:\n")
combined_layers %>%
  group_by(dataset) %>%
  summarise(
    median_OrgC_pct       = median(OrgC_pct, na.rm = TRUE),
    median_OrgC_Stock_kgm2 = median(OrgC_Stock_kgm2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# QC: check for duplicate layer_id within a dataset
dup_layers <- combined_layers %>%
  count(dataset, layer_id) %>%
  filter(n > 1)
if (nrow(dup_layers) > 0) {
  warning(sprintf("%d duplicate dataset+layer_id combinations found", nrow(dup_layers)))
} else {
  message("OK: No duplicate layer_id within any dataset.")
}

write_csv(combined_layers, PATH_OUT_LAYERS)
message(sprintf("Saved: %s  (%d rows)", PATH_OUT_LAYERS, nrow(combined_layers)))

# -----------------------------------------------------------------------------
# 6. HARMONIZE WOSIS PROFILES
# -----------------------------------------------------------------------------
message("\n--- Harmonizing WOSIS profiles ---")

# wosis_profiles.csv does not contain mean_BDOD or mean_OrgC_pct
# Compute them from the harmonized layer table built in Section 2

wosis_profile_means <- wosis_l_harm %>%
  group_by(profile_id) %>%
  summarise(
    mean_BDOD     = mean(BDOD, na.rm = TRUE),
    mean_OrgC_pct = mean(OrgC_pct, na.rm = TRUE),
    .groups = "drop"
  )

wosis_p_harm <- wosis_p %>%
  mutate(
    profile_id          = as.character(profile_id),
    sum_OrgC_Stock_kgm2 = sum_OrgC_Stock,          # already kg/m²
    year                = as.integer(substr(as.character(date), 1, 4))
  ) %>%
  left_join(wosis_profile_means, by = "profile_id") %>%
  select(
    dataset,
    profile_id,
    latitude,
    longitude,
    total_depth_cm,
    n_layers,
    sum_OrgC_Stock_kgm2,
    mean_BDOD,
    mean_OrgC_pct,
    country_name,
    year
  )

message(sprintf("WOSIS profiles harmonized: %d rows", nrow(wosis_p_harm)))

# -----------------------------------------------------------------------------
# 7. HARMONIZE CANPEAT PROFILES
# -----------------------------------------------------------------------------
message("\n--- Harmonizing CANPEAT profiles ---")

# peat_profiles_gee.csv has Total_C_MgCha (Mg C/ha) → convert to kg/m² (× 0.1)
# n_samples → n_layers; Latitude/Longitude → lowercase

peat_p_harm <- peat_p %>%
  mutate(
    profile_id          = as.character(CORE_ID),
    latitude            = Latitude,
    longitude           = Longitude,
    total_depth_cm      = Total_depth_cm,
    n_layers            = n_samples,
    sum_OrgC_Stock_kgm2 = Total_C_MgCha * 0.1,    # Mg C/ha → kg/m²
    country_name        = "Canada",
    year                = as.integer(Year)
  ) %>%
  select(
    dataset,
    profile_id,
    latitude,
    longitude,
    total_depth_cm,
    n_layers,
    sum_OrgC_Stock_kgm2,
    mean_BDOD,
    mean_OrgC_pct,
    country_name,
    year
  )

message(sprintf("CANPEAT profiles harmonized: %d rows", nrow(peat_p_harm)))

# -----------------------------------------------------------------------------
# 8. HARMONIZE JANOUSEK PROFILES
# -----------------------------------------------------------------------------
message("\n--- Harmonizing Janousek profiles ---")

# Script 03 already produced harmonized column names
# Enforce types; add NA columns for country_name and year

jan_p_harm <- jan_p %>%
  mutate(
    profile_id   = as.character(profile_id),
    country_name = NA_character_,
    year         = NA_integer_
  ) %>%
  select(
    dataset,
    profile_id,
    latitude,
    longitude,
    total_depth_cm,
    n_layers,
    sum_OrgC_Stock_kgm2,
    mean_BDOD,
    mean_OrgC_pct,
    country_name,
    year
  )

message(sprintf("Janousek profiles harmonized: %d rows", nrow(jan_p_harm)))

# -----------------------------------------------------------------------------
# 9. ROW-BIND PROFILES + QC + EXPORT
# -----------------------------------------------------------------------------
message("\n--- Combining all profile tables ---")

combined_profiles <- bind_rows(wosis_p_harm, peat_p_harm, jan_p_harm)

message(sprintf("Combined profile table: %d rows", nrow(combined_profiles)))
message(sprintf("  WOSIS:    %d", nrow(wosis_p_harm)))
message(sprintf("  CANPEAT:  %d", nrow(peat_p_harm)))
message(sprintf("  Janousek: %d", nrow(jan_p_harm)))

# QC: verify row count
expected_profiles <- nrow(wosis_p_harm) + nrow(peat_p_harm) + nrow(jan_p_harm)
if (nrow(combined_profiles) != expected_profiles) {
  warning(sprintf("Profile row count mismatch: expected %d, got %d",
                  expected_profiles, nrow(combined_profiles)))
} else {
  message("OK: Combined profile count matches sum of inputs.")
}

# QC: row counts by dataset
cat("\nProfile rows per dataset:\n")
print(count(combined_profiles, dataset))

# QC: stock and mean values by dataset
cat("\nstock/means sanity check by dataset:\n")
combined_profiles %>%
  group_by(dataset) %>%
  summarise(
    n               = n(),
    median_stock    = median(sum_OrgC_Stock_kgm2, na.rm = TRUE),
    median_BDOD     = median(mean_BDOD, na.rm = TRUE),
    median_OrgC_pct = median(mean_OrgC_pct, na.rm = TRUE),
    na_stock        = sum(is.na(sum_OrgC_Stock_kgm2)),
    na_BDOD         = sum(is.na(mean_BDOD)),
    .groups = "drop"
  ) %>%
  print()

write_csv(combined_profiles, PATH_OUT_PROFILES)
message(sprintf("Saved: %s  (%d profiles)", PATH_OUT_PROFILES, nrow(combined_profiles)))

message("\n=== Combination complete ===")
message(sprintf("combined_layers.csv:   %d rows across %d databases",
                nrow(combined_layers),
                n_distinct(combined_layers$dataset)))
message(sprintf("combined_profiles.csv: %d rows across %d databases",
                nrow(combined_profiles),
                n_distinct(combined_profiles$dataset)))
