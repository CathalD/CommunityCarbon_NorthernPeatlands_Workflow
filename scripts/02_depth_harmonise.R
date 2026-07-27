# =============================================================================
# 02_depth_harmonise.R
#
# Turn flagged segments into core-level carbon stocks that are only ever
# compared on a like-for-like basis.
#
# Two depth windows are produced and they are NOT interchangeable:
#
#   COMMON SUPPORT  0 to the deepest depth EVERY core reaches. Derived from the
#                   data. No extrapolation anywhere; all eight cores are
#                   directly comparable. This is the window in which the
#                   sample-based products are honest.
#
#   REFERENCE       0-30 cm, chosen to match the published layers. Two cores do
#                   not reach 30 cm, so their stocks are LOWER BOUNDS and are
#                   flagged. They are never rescaled to look complete.
#
# Two carbon-conversion scenarios are also produced, because 01 established
# that the two campaigns used different organic-matter-to-carbon factors.
#
# INPUT   data/derived/01_segments_qc.csv
# OUTPUT  data/derived/02_core_stocks.csv
#         outputs/tables/02_depth_support.csv
#         outputs/tables/02_stock_scenarios.csv
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

log_step("02  DEPTH HARMONISATION AND CARBON STOCKS")

seg <- require_artifact(file.path(CFG$dir_derived, "01_segments_qc.csv"),
                        "scripts/01_ingest_qc.R")
log_info("loaded ", nrow(seg), " QC'd segments across ",
         length(unique(seg$core_id)), " cores")

# =============================================================================
# 1. Depth support: what can actually be compared?
# =============================================================================

log_step("02a  DEPTH SUPPORT")

cs <- common_support_depth(seg$core_id, seg$depth_bottom_cm)
z_common <- if (is.null(CFG$window_common_cm)) c(0, cs$depth_cm) else CFG$window_common_cm
z_ref    <- CFG$window_reference_cm

log_info("core lengths, shallowest first:")
for (nm in names(cs$core_totals_cm)) {
  log_info(sprintf("   %-8s %6.3f cm", nm, cs$core_totals_cm[[nm]]))
}
log_ok(sprintf("common support = 0-%.3f cm, limited by core %s",
               cs$depth_cm, cs$limiting_core))
log_info(sprintf("reference window = %.0f-%.0f cm (to match the published layers)",
                 z_ref[1], z_ref[2]))

reach_ref <- cs$core_totals_cm >= z_ref[2] - 1e-8
log_info(sprintf("%d of %d cores reach %.0f cm; %d do not and will be flagged ",
                 sum(reach_ref), length(reach_ref), z_ref[2], sum(!reach_ref)),
         "as lower bounds: ",
         paste(names(cs$core_totals_cm)[!reach_ref], collapse = ", "))

depth_support <- data.frame(
  core_id           = names(cs$core_totals_cm),
  core_length_cm    = as.numeric(cs$core_totals_cm),
  reaches_reference = as.logical(reach_ref),
  common_support_cm = cs$depth_cm,
  stringsAsFactors  = FALSE
)
depth_support$campaign <- seg$campaign[match(depth_support$core_id, seg$core_id)]

# =============================================================================
# 2. Carbon-conversion scenarios
# =============================================================================

log_step("02b  CARBON CONVERSION SCENARIOS")

# Recover, per campaign, the factor actually used, then build a scenario in
# which the peat cores are converted with a peat-appropriate factor instead.
fac <- identify_om_carbon_factor(seg$om_pct, seg$soc_pct, seg$campaign)
log_info("factors implied by the supplied SOC column:")
for (k in seq_len(nrow(fac))) {
  log_info(sprintf("   %-16s %.4f", fac$group[k], fac$factor_est[k]))
}

#' Build a segment table under a named carbon-conversion scenario. Pure.
apply_carbon_scenario <- function(seg, scenario, peat_factor, fac_tbl) {
  s <- seg
  if (identical(scenario, "as_supplied")) {
    s$soc_pct_used <- s$soc_pct
  } else if (identical(scenario, "harmonised_peat")) {
    # Rescale ONLY the peat campaign, from the factor its lab used to a
    # factor measured for northern peat. Mineral cores are left alone: their
    # implied factor is a defensible convention for mineral soil.
    b_used <- fac_tbl$factor_est[match(s$campaign, fac_tbl$group)]
    is_peat <- grepl("peat", s$campaign)
    s$soc_pct_used <- ifelse(is_peat, s$soc_pct * (peat_factor / b_used),
                             s$soc_pct)
  } else {
    stop("unknown scenario: ", scenario)
  }
  s$carbon_density_gcm3 <- carbon_density_gcm3(s$bulk_density_gcm3, s$soc_pct_used)
  s$scenario <- scenario
  s
}

scenarios <- c("as_supplied", "harmonised_peat")
log_info("scenario 'as_supplied'     : SOC column used exactly as provided")
log_info(sprintf(paste0("scenario 'harmonised_peat': peat cores re-converted at ",
                        "%.2f instead of %.4f;"), CFG$factor_peat_lit,
                 fac$factor_est[grepl("peat", fac$group)]))
log_info("                            mineral cores unchanged")

# =============================================================================
# 3. Integrate
# =============================================================================

log_step("02c  INTEGRATION")

windows <- list(
  list(name = "common_support", z = z_common,
       status = "interpolation-ready: every core fully covers this window"),
  list(name = "reference_0_30", z = z_ref,
       status = "reference-comparable: short cores are lower bounds")
)

res <- list()
for (sc in scenarios) {
  s <- apply_carbon_scenario(seg, sc, CFG$factor_peat_lit, fac)
  for (w in windows) {
    r <- integrate_all_cores(s, w$z[1], w$z[2],
                             extrapolate = CFG$report_extrapolated_variant)
    r$scenario    <- sc
    r$window      <- w$name
    r$window_note <- w$status
    res[[length(res) + 1L]] <- r
  }
}
stocks <- do.call(rbind, res)

# Attach core metadata.
meta <- unique(seg[, c("core_id", "campaign", "year", "latitude", "longitude")])
stocks <- merge(stocks, meta, by = "core_id", all.x = TRUE)
stocks$stock_tha <- kgm2_to_tha(stocks$stock_kgm2)

# A stock is usable as a plain number only when it is not a lower bound.
stocks$usable_as_point_estimate <- !stocks$is_lower_bound

stocks <- stocks[order(stocks$scenario, stocks$window, stocks$campaign,
                       stocks$core_id), ]
rownames(stocks) <- NULL

for (sc in scenarios) {
  for (w in vapply(windows, `[[`, character(1), "name")) {
    d <- stocks[stocks$scenario == sc & stocks$window == w, ]
    log_info(sprintf("scenario=%s  window=%s", sc, w))
    print_table(d[, c("core_id", "campaign", "covered_cm", "stock_kgm2",
                      "is_lower_bound")], digits = 3)
    nlb <- sum(d$is_lower_bound)
    if (nlb) {
      log_flag(nlb, " core(s) are LOWER BOUNDS in this window: ",
               paste(d$core_id[d$is_lower_bound], collapse = ", "))
    }
  }
}

# =============================================================================
# 4. What the numbers say
# =============================================================================

log_step("02d  READING THE RESULT")

# The comparison must be made in the window where every core is valid. In the
# 0-30 cm window only one peat core is complete, so a peat-versus-mineral
# contrast there would rest on a single peat observation.
d30 <- stocks[stocks$scenario == "as_supplied" &
              stocks$window == "reference_0_30", ]
n_complete_30 <- tapply(!d30$is_lower_bound, d30$campaign, sum)
log_info("0-30 cm window, complete cores per campaign: ",
         paste(sprintf("%s=%d", names(n_complete_30),
                       as.integer(n_complete_30)), collapse = "  "))
if (any(n_complete_30 < 2L)) {
  log_warn("At 0-30 cm, ",
           paste(names(n_complete_30)[n_complete_30 < 2L], collapse = ", "),
           " has fewer than two complete cores.")
  log_warn("No campaign contrast is drawn in this window; it would rest on a ",
           "single observation. The contrast below uses the common-support ",
           "window instead, where all eight cores are valid.")
}

d <- stocks[stocks$scenario == "as_supplied" &
            stocks$window == "common_support", ]
by_camp <- tapply(d$stock_kgm2, d$campaign, mean)
log_info(sprintf("mean 0-%.1f cm stock (common support, all cores valid, kg C/m2):",
                 z_common[2]))
for (nm in names(by_camp)) {
  log_info(sprintf("   %-16s %6.2f  (n=%d)", nm, by_camp[[nm]],
                   sum(d$campaign == nm)))
}

if (length(by_camp) == 2L && by_camp[["PM_2024_peat"]] < by_camp[["FS_2025_mineral"]]) {
  log_warn(strrep("-", 68))
  log_warn(sprintf("Over 0-%.1f cm the PEAT cores hold LESS carbon (%.2f kg/m2)",
                   z_common[2], by_camp[["PM_2024_peat"]]))
  log_warn(sprintf("than the MINERAL cores (%.2f kg/m2).",
                   by_camp[["FS_2025_mineral"]]))
  log_warn("This is not an error. Peat here has bulk density near")
  log_warn("0.05-0.27 g/cm3 against 0.75-2.18 for the mineral soils, so a")
  log_warn("high carbon PERCENTAGE still yields little carbon MASS in a fixed")
  log_warn("volume.")
  log_warn("")
  log_warn("The consequence matters for this whole exercise: a 0-30 cm")
  log_warn("accounting depth systematically under-represents peatland carbon,")
  log_warn("because Hudson Bay Lowlands peat is metres thick and almost all")
  log_warn("of its stock lies BELOW the window every published 0-30 cm")
  log_warn("product reports. A 0-30 cm map of this landscape is not a map of")
  log_warn("its carbon; it is a map of the top of its carbon.")
  log_warn(strrep("-", 68))
}

# Effect of the conversion scenario on the peat cores.
a <- stocks[stocks$scenario == "as_supplied" & stocks$window == "reference_0_30", ]
h <- stocks[stocks$scenario == "harmonised_peat" & stocks$window == "reference_0_30", ]
cmp <- merge(a[, c("core_id", "campaign", "stock_kgm2")],
             h[, c("core_id", "stock_kgm2")], by = "core_id",
             suffixes = c("_supplied", "_harmonised"))
cmp$pct_change <- 100 * (cmp$stock_kgm2_harmonised / cmp$stock_kgm2_supplied - 1)
log_info("effect of re-converting peat carbon at ", CFG$factor_peat_lit, ":")
print_table(cmp, digits = 3)
log_info(sprintf("peat cores shift by %.1f%%; mineral cores are untouched",
                 unique(round(cmp$pct_change[grepl("peat", cmp$campaign)], 1))[1]))

if (CFG$report_extrapolated_variant) {
  e <- stocks[stocks$scenario == "as_supplied" &
              stocks$window == "reference_0_30" & stocks$is_lower_bound, ]
  if (nrow(e)) {
    log_info("EXTRAPOLATED variant for the two short cores (scenario only, ",
             "never a default product):")
    print_table(e[, c("core_id", "covered_cm", "gap_cm", "stock_kgm2",
                      "stock_extrap_kgm2", "extrap_share_of_stock")], digits = 3)
    log_warn("Between ",
             sprintf("%.0f%% and %.0f%%",
                     100 * min(e$extrap_share_of_stock),
                     100 * max(e$extrap_share_of_stock)),
             " of the extrapolated stock is invented by carrying the")
    log_warn("deepest observed density downward.")
    # Whether that is a conservative or an inflationary assumption depends on
    # what each core was actually doing at its base. Read it, do not assume it.
    for (cc in e$core_id) {
      s <- seg[seg$core_id == cc, ]
      s <- s[order(s$seg_index), ]
      if (nrow(s) < 2L) {
        log_warn(sprintf("   %-8s single segment: no depth trend exists, so ",
                         cc),
                 "the carried density is an assumption with no support at all.")
      } else {
        last2 <- utils::tail(s$carbon_density_gcm3, 2)
        dir <- if (last2[2] < last2[1]) "FALLING" else "rising"
        log_warn(sprintf(paste0("   %-8s carbon density is %s at the base ",
                                "(%.5f -> %.5f g/cm3)"),
                         cc, dir, last2[1], last2[2]))
      }
    }
    log_warn("Treat the extrapolated column as a scenario, not a value.")
  }
}

# =============================================================================
# 5. Write
# =============================================================================

log_step("02e  OUTPUTS")

write_csv_logged(stocks, file.path(CFG$dir_derived, "02_core_stocks.csv"),
                 "core stocks by window and conversion scenario")
write_csv_logged(depth_support,
                 file.path(CFG$dir_tables, "02_depth_support.csv"))
write_csv_logged(cmp, file.path(CFG$dir_tables, "02_stock_scenarios.csv"),
                 "sensitivity to the OM-to-carbon factor")

log_ok("02 complete")
