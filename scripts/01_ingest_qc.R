# =============================================================================
# 01_ingest_qc.R
#
# Read the community core data, recompute every derived quantity from first
# principles, reconcile it against what was supplied, and flag what is
# implausible. Nothing is dropped.
#
# Two classes of problem stop the run outright: a coordinate that would place
# the cores on the wrong continent, and a supplied column that does not mean
# what its header says.
#
# INPUT   data/raw/community_soil_cores.csv
# OUTPUT  data/derived/01_segments_qc.csv        one row per segment, flagged
#         data/derived/01_cores_qc.csv           one row per core
#         outputs/tables/01_flag_dictionary.csv
#         outputs/tables/01_flag_summary.csv
#         outputs/tables/01_reconciliation.csv
#         outputs/tables/01_om_carbon_factor.csv
#         outputs/tables/01_prior_pipeline_audit.csv
# =============================================================================

## --- bootstrap ---------------------------------------------------------------
.root <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  d <- if (length(f)) dirname(f[[1]]) else getwd()
  while (!file.exists(file.path(d, "config.R"))) {
    p <- dirname(d)
    if (identical(p, d)) stop("config.R not found above '", getwd(), "'", call. = FALSE)
    d <- p
  }
  normalizePath(d)
})
source(file.path(.root, "config.R"))

log_step("01  INGEST AND QUALITY CONTROL")

# =============================================================================
# 1. Read, with headers exactly as supplied
# =============================================================================

raw <- read_csv_verbatim(CFG$file_cores_raw)
log_info("read ", nrow(raw), " rows x ", ncol(raw), " columns from ",
         basename(CFG$file_cores_raw))

gate_columns(raw, CFG$raw_cols)
log_ok("all ", length(CFG$raw_cols), " expected columns present, names verbatim")

# Rename to safe internal names. The originals are preserved alongside so the
# reconciliation audit can quote what was actually supplied.
seg <- data.frame(
  year               = as.integer(raw[["year"]]),
  core_id            = trimws(as.character(raw[["Core Id"]])),
  sample_id          = trimws(as.character(raw[["Sample Id"]])),
  latitude           = as.numeric(raw[["Latitude"]]),
  longitude_supplied = as.numeric(raw[["Longitude"]]),
  thickness_cm       = as.numeric(raw[["Depth"]]),
  bulk_density_gcm3  = as.numeric(raw[["Bulk Density"]]),
  om_pct             = as.numeric(raw[["OM"]]),
  soc_pct            = as.numeric(raw[["SOC"]]),
  ocd_supplied       = as.numeric(raw[["Organic Carbon Density (g/cm^2)"]]),
  stringsAsFactors   = FALSE
)

# Campaign stratum. This is carried through every later step and never pooled
# away: the two campaigns differ in year, landscape type, operator and lab
# batch simultaneously.
seg$campaign <- ifelse(grepl("^PM", seg$core_id), "PM_2024_peat",
                ifelse(grepl("^FS", seg$core_id), "FS_2025_mineral", NA))
if (any(is.na(seg$campaign))) {
  fail_loudly("Core identifier matches no known campaign prefix",
              paste0("Unrecognised: ",
                     paste(unique(seg$core_id[is.na(seg$campaign)]), collapse = ", ")),
              remedy = "Extend the campaign rule in 01_ingest_qc.R.")
}
log_ok("assigned campaign stratum: ",
       paste(sprintf("%s n=%d", names(table(seg$campaign)),
                     as.integer(table(seg$campaign))), collapse = "  |  "))

# =============================================================================
# 2. HARD GATE: geolocation
# =============================================================================

log_step("01a  GEOLOCATION GATE")
log_info("stated study area: ", CFG$study_area)
log_info(sprintf("supplied longitude range: %+.5f to %+.5f",
                 min(seg$longitude_supplied), max(seg$longitude_supplied)))
log_info(sprintf("supplied latitude  range: %+.5f to %+.5f",
                 min(seg$latitude), max(seg$latitude)))

geo <- gate_geolocation(seg$longitude_supplied, seg$latitude,
                        anchor   = CFG$site_anchor,
                        bbox     = CFG$plausible_bbox,
                        allow_fix = CFG$allow_longitude_sign_fix)

seg$longitude   <- geo$longitude
seg$lon_repaired <- geo$repaired
flags_geo <- geo$flags

# Distance from the stated anchor, as an independent confirmation that the
# repaired coordinates land where the study area says they should.
seg$dist_from_fort_severn_km <- haversine_km(
  seg$longitude, seg$latitude,
  CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]]
)
log_ok(sprintf("after repair, cores lie %.1f-%.1f km from Fort Severn",
               min(seg$dist_from_fort_severn_km),
               max(seg$dist_from_fort_severn_km)))
if (max(seg$dist_from_fort_severn_km) > 80) {
  log_flag("a core sits further than 80 km from the stated anchor")
}

# =============================================================================
# 3. Segment order and depth geometry
# =============================================================================

log_step("01b  SEGMENT ORDER AND DEPTH SEMANTICS")

seg$seg_index <- trailing_int(seg$sample_id)
absent <- is.na(seg$seg_index)
if (any(absent)) {
  log_flag(sum(absent), " sample id(s) carry no trailing integer: ",
           paste(seg$sample_id[absent], collapse = ", "))
  log_info("assigning index 1 to these; each belongs to a single-segment core")
  # Only safe when the core truly has one segment. Verify before assigning.
  for (i in which(absent)) {
    n_in_core <- sum(seg$core_id == seg$core_id[i])
    if (n_in_core > 1L) {
      fail_loudly(
        "Cannot determine segment order",
        paste0("Sample '", seg$sample_id[i], "' has no trailing integer, but ",
               "core '", seg$core_id[i], "' has ", n_in_core, " segments."),
        remedy = "Segment order drives depth assignment. Obtain the ordering from the field records."
      )
    }
  }
  seg$seg_index[absent] <- 1L
}
flags_idx <- replicate(nrow(seg), character(0), simplify = FALSE)
for (i in which(absent)) flags_idx[[i]] <- "SEGMENT_INDEX_ABSENT"

# Is `Depth` a thickness or a cumulative bottom depth? Decide from evidence.
sem <- infer_depth_semantics(seg$core_id, seg$seg_index, seg$thickness_cm)
log_info("depth column semantics test: ", sem$reasoning)
if (sem$n_decreases > 0L) {
  dec <- sem$monotonic_violations[sem$monotonic_violations$kind == "decrease", ]
  for (k in seq_len(nrow(dec))) {
    log_info(sprintf("   decisive: %-8s %.3f -> %.3f", dec$core_id[k],
                     dec$from[k], dec$to[k]))
  }
}
if (sem$n_ties > 0L) {
  log_info(sprintf("   plus %d step(s) repeating the same value, in %d core(s)",
                   sem$n_ties,
                   length(unique(sem$monotonic_violations$core_id[
                     sem$monotonic_violations$kind == "tie"]))))
}
if (sem$verdict != CFG$depth_col_semantics) {
  if (sem$verdict == "ambiguous") {
    log_warn("monotonicity alone cannot settle the reading; ",
             "proceeding on the configured assumption '",
             CFG$depth_col_semantics, "' and corroborating below")
  } else {
    fail_loudly("Depth column semantics contradict the configuration",
                paste0("Evidence says: ", sem$verdict),
                paste0("config.R says: ", CFG$depth_col_semantics),
                remedy = "Resolve before computing any stock.")
  }
}

seg <- add_depth_bounds(seg)

# Corroboration: do the accumulated thicknesses reproduce sensible core lengths?
core_len <- tapply(seg$depth_bottom_cm, seg$core_id, max)
log_ok(sprintf("core lengths from accumulated thickness: %.1f-%.1f cm",
               min(core_len), max(core_len)))
log_info("  ", paste(sprintf("%s=%.3f", names(core_len), core_len), collapse = "  "))

# =============================================================================
# 4. Recompute every derived quantity from first principles
# =============================================================================

log_step("01c  RECOMPUTATION AND RECONCILIATION")

seg$carbon_density_gcm3 <- carbon_density_gcm3(seg$bulk_density_gcm3, seg$soc_pct)

rec <- gate_ocd_reconciliation(seg$ocd_supplied, seg$carbon_density_gcm3,
                               tol = CFG$qc$ocd_rel_tol)
seg$ocd_reconciles  <- rec$ok
seg$ocd_ratio       <- seg$ocd_supplied / seg$carbon_density_gcm3

if (rec$n_bad == 0) {
  log_ok("supplied 'Organic Carbon Density' reconciles EXACTLY with ",
         "bulk density x SOC/100 on all ", nrow(seg), " rows (rel. tol ",
         CFG$qc$ocd_rel_tol, ")")
} else {
  log_flag(rec$n_bad, " row(s) fail reconciliation")
}

# --- the unit finding --------------------------------------------------------
# The column reconciles as BD [g/cm3] * SOC [%] / 100, which is dimensionally
# g/cm3: a carbon CONCENTRATION per unit volume of soil. The supplied header
# says g/cm2, which is an AREAL density. These differ by one length dimension.
# Summing this column down a core -- as the prior pipeline did -- adds
# concentrations and yields a number with no physical meaning unless every
# segment happens to be 1 cm thick.
log_warn(strrep("-", 68))
log_warn("UNIT DISCREPANCY in supplied column 'Organic Carbon Density (g/cm^2)'")
log_warn("  reconciles exactly as: bulk density [g/cm3] * SOC [%] / 100")
log_warn("  that quantity has units g/cm3 (carbon per unit VOLUME),")
log_warn("  not g/cm2 (carbon per unit AREA) as the header states.")
log_warn("  An areal stock requires multiplication by segment thickness.")
log_warn("  This workflow treats the column as g/cm3 and integrates properly;")
log_warn("  the supplied values are retained unaltered for audit.")
log_warn(strrep("-", 68))

seg$ocd_supplied_units_asserted <- "g/cm^2"
seg$ocd_units_evidenced         <- "g/cm^3"

reconciliation <- data.frame(
  sample_id            = seg$sample_id,
  bulk_density_gcm3    = seg$bulk_density_gcm3,
  soc_pct              = seg$soc_pct,
  ocd_supplied         = seg$ocd_supplied,
  ocd_recomputed_gcm3  = seg$carbon_density_gcm3,
  abs_diff             = abs(seg$ocd_supplied - seg$carbon_density_gcm3),
  ratio                = seg$ocd_ratio,
  reconciles           = seg$ocd_reconciles,
  stringsAsFactors     = FALSE
)

# --- organic matter to carbon factor forensics -------------------------------
fac <- identify_om_carbon_factor(seg$om_pct, seg$soc_pct, seg$campaign)
fac$factor_name <- vapply(fac$factor_est, name_om_carbon_factor,
                          character(1), tol = CFG$factor_match_tol)
log_info("organic-matter-to-carbon factor implied by the supplied SOC column:")
for (k in seq_len(nrow(fac))) {
  log_info(sprintf("   %-16s n=%2d  factor=%.4f  (1/%.3f)  -> %s",
                   fac$group[k], fac$n[k], fac$factor_est[k],
                   fac$implied_divisor[k], fac$factor_name[k]))
}
if (length(unique(round(fac$factor_est, 2))) > 1L) {
  log_warn(strrep("-", 68))
  log_warn("THE TWO CAMPAIGNS USED DIFFERENT OM-TO-CARBON CONVERSION FACTORS.")
  log_warn("  SOC is therefore not a single measured quantity across the")
  log_warn("  dataset; it is two differently-derived quantities sharing a")
  log_warn("  column. Campaign must be preserved as a stratum, and any")
  log_warn("  peat-versus-mineral contrast is partly a lab-convention")
  log_warn("  artefact rather than a property of the landscape.")
  log_warn("  The peat cores carry the larger factor, which is the direction")
  log_warn("  that OVERSTATES carbon in peat: measured C:OM for northern peat")
  log_warn("  is nearer ", CFG$factor_peat_lit, " than ", CFG$factor_van_bemmelen, ".")
  log_warn("  Sensitivity to this is quantified in 05.")
  log_warn(strrep("-", 68))
}

# =============================================================================
# 5. Value flags -- surfaced, never dropped
# =============================================================================

log_step("01d  VALUE PLAUSIBILITY FLAGS")

flags_val <- flag_values(seg, CFG$qc)
flags_dup <- flag_cross_core_duplicates(seg$core_id, seg$om_pct, seg$soc_pct)

flags_ocd <- replicate(nrow(seg), character(0), simplify = FALSE)
for (i in which(!seg$ocd_reconciles)) flags_ocd[[i]] <- "OCD_RECONCILE_FAIL"

n_seg_per_core <- table(seg$core_id)
flags_single <- replicate(nrow(seg), character(0), simplify = FALSE)
for (i in which(n_seg_per_core[seg$core_id] == 1L)) {
  flags_single[[i]] <- "SINGLE_SEGMENT_CORE"
}

flags <- merge_flags(flags_geo, flags_idx, flags_val, flags_dup,
                     flags_ocd, flags_single)
seg$qc_flags   <- flags_to_string(flags)
seg$qc_n_flags <- lengths(flags)

flag_summary <- summarise_flags(flags)
dict <- qc_flag_dictionary()
flag_summary$meaning <- dict$meaning[match(flag_summary$flag, dict$flag)]

if (nrow(flag_summary)) {
  log_info("flags raised (rows are kept, not removed):")
  for (k in seq_len(nrow(flag_summary))) {
    log_flag(sprintf("%-22s %2d row(s)", flag_summary$flag[k],
                     flag_summary$n_rows[k]))
  }
} else {
  log_ok("no value flags raised")
}
log_info(sum(seg$qc_n_flags > 0), " of ", nrow(seg),
         " segments carry at least one flag; all retained")

# Narrate the specific flags that change how a reader should treat the data.
for (i in which(grepl("BD_ABOVE_MAX", seg$qc_flags))) {
  log_flag(sprintf("%-12s bulk density %.3f g/cm3 exceeds the %.2f ceiling",
                   seg$sample_id[i], seg$bulk_density_gcm3[i], CFG$qc$bd_max))
}
for (i in which(grepl("DUPLICATE_LAB_VALUES", seg$qc_flags))) {
  log_flag(sprintf("%-12s OM=%.2f SOC=%.2f identical to another core's sample",
                   seg$sample_id[i], seg$om_pct[i], seg$soc_pct[i]))
}

# =============================================================================
# 6. Core-level table
# =============================================================================

agg <- function(f, col) tapply(seg[[col]], seg$core_id, f)
cores <- data.frame(
  core_id          = names(agg(length, "thickness_cm")),
  campaign         = agg(function(x) x[1], "campaign"),
  year             = agg(function(x) x[1], "year"),
  latitude         = agg(function(x) x[1], "latitude"),
  longitude        = agg(function(x) x[1], "longitude"),
  n_segments       = as.integer(agg(length, "thickness_cm")),
  core_length_cm   = agg(max, "depth_bottom_cm"),
  bd_min_gcm3      = agg(min, "bulk_density_gcm3"),
  bd_max_gcm3      = agg(max, "bulk_density_gcm3"),
  om_min_pct       = agg(min, "om_pct"),
  om_max_pct       = agg(max, "om_pct"),
  n_flagged_segs   = as.integer(tapply(seg$qc_n_flags > 0, seg$core_id, sum)),
  stringsAsFactors = FALSE
)
cores$dist_from_fort_severn_km <- haversine_km(
  cores$longitude, cores$latitude,
  CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]])
cores <- cores[order(cores$campaign, cores$core_id), ]
rownames(cores) <- NULL

# Within-core coordinate consistency: every segment of a core must share one
# location. A disagreement means the core was mis-keyed.
for (cc in unique(seg$core_id)) {
  s <- seg[seg$core_id == cc, ]
  if (length(unique(s$latitude)) > 1L || length(unique(s$longitude)) > 1L) {
    fail_loudly("Segments of one core carry different coordinates",
                paste0("Core: ", cc),
                remedy = "Resolve the keying error before any spatial step.")
  }
}
log_ok("every core has internally consistent coordinates")

log_info("core inventory:")
print_table(cores[, c("core_id", "campaign", "n_segments", "core_length_cm",
                      "bd_min_gcm3", "bd_max_gcm3", "om_max_pct",
                      "n_flagged_segs")], digits = 3)

# =============================================================================
# 7. Audit the prior pipeline output, if present
# =============================================================================

prior_audit <- NULL
if (file.exists(CFG$file_prior_pipeline)) {
  log_step("01e  AUDIT OF PRIOR PIPELINE OUTPUT")
  log_info("this file is NEVER used as an analysis input; it is checked only")
  log_info("so that errors it introduced are not repeated silently")
  p <- read_csv_verbatim(CFG$file_prior_pipeline)

  # Its 'total_ocd_naive_sum' adds the g/cm3 column down each core. Compare
  # that against a properly integrated areal stock over the same interval.
  proper <- integrate_all_cores(seg, 0, max(seg$depth_bottom_cm))
  m <- match(p$core_id, proper$core_id)
  prior_audit <- data.frame(
    core_id                = p$core_id,
    prior_longitude        = p$longitude,
    our_longitude          = cores$longitude[match(p$core_id, cores$core_id)],
    prior_total_ocd_naive  = p$total_ocd_naive_sum,
    proper_stock_gcm2      = proper$stock_gcm2[m],
    ratio_proper_to_naive  = proper$stock_gcm2[m] / p$total_ocd_naive_sum,
    core_length_cm         = proper$covered_cm[m],
    stringsAsFactors       = FALSE
  )
  log_info("prior 'total_ocd_naive_sum' versus a correctly integrated stock:")
  print_table(prior_audit, digits = 4)
  log_warn("The prior sum adds carbon CONCENTRATIONS (g/cm3) down the core ",
           "instead of integrating concentration x thickness.")
  log_warn("Where every segment is 10 cm the error is a clean factor of 10.")
  log_warn("Where segment thickness varies -- all three PM peat cores -- the ",
           "error varies between cores and cannot be undone by rescaling.")
  log_warn("The prior file also silently negated the longitude without ",
           "recording that it had done so, and then labelled the result ",
           "'lat_lon_consistent = TRUE'.")
}

# =============================================================================
# 8. Write
# =============================================================================

log_step("01f  OUTPUTS")

seg_out <- seg[order(seg$campaign, seg$core_id, seg$seg_index), ]
rownames(seg_out) <- NULL

write_csv_logged(seg_out, file.path(CFG$dir_derived, "01_segments_qc.csv"),
                 "segment level, QC flagged, depths resolved")
write_csv_logged(cores, file.path(CFG$dir_derived, "01_cores_qc.csv"),
                 "core level inventory")
write_csv_logged(dict, file.path(CFG$dir_tables, "01_flag_dictionary.csv"),
                 "what every flag means")
write_csv_logged(flag_summary, file.path(CFG$dir_tables, "01_flag_summary.csv"))
write_csv_logged(reconciliation,
                 file.path(CFG$dir_tables, "01_reconciliation.csv"),
                 "supplied vs first-principles carbon density")
write_csv_logged(fac, file.path(CFG$dir_tables, "01_om_carbon_factor.csv"),
                 "which OM-to-C factor each campaign used")
if (!is.null(prior_audit)) {
  write_csv_logged(prior_audit,
                   file.path(CFG$dir_tables, "01_prior_pipeline_audit.csv"))
}

log_ok("01 complete")
