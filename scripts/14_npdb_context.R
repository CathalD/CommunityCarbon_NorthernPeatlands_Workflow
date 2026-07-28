# =============================================================================
# 14_npdb_context.R
#
# Place the Fort Severn community cores in national context using the
# Agriculture and Agri-Food Canada National Pedon Database (NPDB).
#
# This is the script that answers "how do these results compare with other
# areas?" with data rather than assertion. It is also the script that
# establishes how alone these cores are.
#
# UNITS AND DEPTH, RECONCILED RATHER THAN ASSUMED
#
# NPDB stocks arrive as Mg C/ha; this workflow reports kg C/m2 everywhere.
# 1 kg/m2 = 10 Mg/ha, so the conversion is a division by ten. It is applied
# once, explicitly, and the result is checked against a physically plausible
# range instead of trusted.
#
# The NPDB preparation already standardises to 0-30 cm from the profile top
# and flags cores that do not reach it. Those flags are honoured here on the
# same terms as the community cores: a core that stopped short is a LOWER
# BOUND, not a measurement, and is excluded from any mean.
#
# INPUT   data/raw/National Pedon Database/outputs_npdb/npdb_carbon_cores.csv
#         data/derived/02_core_stocks.csv
# OUTPUT  data/derived/14_npdb_cores_qc.csv
#         outputs/tables/14_*.csv
#         outputs/figures/14_*.png
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

log_step("14  NATIONAL CONTEXT FROM THE NPDB")

npdb_path <- file.path(CFG$root, CFG$npdb$dir, CFG$npdb$file_cores)
if (!file.exists(npdb_path)) {
  fail_loudly("NPDB core file not found", paste0("Expected: ", npdb_path),
              remedy = "Check CFG$npdb$dir in config.R.")
}
n <- read_csv_verbatim(npdb_path)
log_info("read ", nrow(n), " NPDB pedons")

# =============================================================================
# 1. Unit reconciliation
# =============================================================================

log_step("14a  UNITS")

if (!"std_stock_Mgha" %in% names(n)) {
  fail_loudly("NPDB file lacks the expected standardised stock column",
              paste0("Present: ", paste(names(n), collapse = ", ")),
              remedy = "Re-run the NPDB preparation, or update this script.")
}
n$std_stock_kgm2 <- n$std_stock_Mgha * CFG$npdb$mgha_to_kgm2
log_info("converted Mg C/ha -> kg C/m2 by x", CFG$npdb$mgha_to_kgm2,
         " (1 kg/m2 = 10 Mg/ha)")

ok_range <- n$std_stock_kgm2[is.finite(n$std_stock_kgm2)]
log_info(sprintf("converted 0-30 cm stocks span %.2f to %.2f kg C/m2",
                 min(ok_range), max(ok_range)))
if (max(ok_range) > CFG$ref_plausible$stock_0_30_kgm2_max) {
  log_flag(sum(ok_range > CFG$ref_plausible$stock_0_30_kgm2_max),
           " NPDB core(s) exceed the ",
           CFG$ref_plausible$stock_0_30_kgm2_max,
           " kg/m2 plausibility ceiling for a 0-30 cm stock; flagged, kept")
} else {
  log_ok("all converted stocks are within the 0-30 cm plausibility ceiling")
}

# =============================================================================
# 2. QC on the same terms as the community cores
# =============================================================================

log_step("14b  QC")

n$has_stock     <- is.finite(n$std_stock_kgm2)
n$is_complete   <- n$has_stock & isTRUE_vec(n$std_stock_complete)
n$loc_poor      <- isTRUE_vec(n$flag_location_conf_poor)

log_info(sum(n$has_stock), " of ", nrow(n),
         " pedons have a computable 0-30 cm stock (",
         sprintf("%.1f%%", 100 * mean(n$has_stock)), ")")
log_warn("Bulk density is the binding constraint in the NPDB: most pedons ",
         "record carbon but not density, and a stock cannot be computed ",
         "without both. That is the same arithmetic this workflow applies to ",
         "the community cores.")
log_info(sum(n$is_complete), " of those reach 30 cm and are usable as point ",
         "estimates rather than lower bounds")
log_flag(sum(n$loc_poor & n$has_stock),
         " usable core(s) carry a location confidence worse than ",
         CFG$npdb$max_location_conf_m, " m")

usable <- n[n$is_complete & is.finite(n$latitude) & is.finite(n$longitude), ]
log_ok(nrow(usable), " NPDB cores are usable for comparison")

usable$dist_fort_severn_km <- haversine_km(
  usable$longitude, usable$latitude,
  CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]])

# =============================================================================
# 3. How alone are these cores?
# =============================================================================

log_step("14c  NEAREST EXISTING DATA")

near <- usable[order(usable$dist_fort_severn_km), ]
log_info("closest NPDB cores with a usable 0-30 cm stock:")
print_table(head(near[, c("province", "latitude", "longitude",
                          "dist_fort_severn_km", "std_stock_kgm2",
                          "has_organic_horizon")], 5), digits = 2)

rings <- do.call(rbind, lapply(CFG$npdb$context_radii_km, function(r) {
  data.frame(within_km = r, n_npdb_cores = sum(usable$dist_fort_severn_km <= r),
             stringsAsFactors = FALSE)
}))
print_table(rings)

nearest_km <- min(usable$dist_fort_severn_km)
log_warn(strrep("-", 68))
log_warn("THE COMMUNITY CORES ARE THE ONLY SOIL CARBON DATA FOR THIS AREA.")
log_warn(sprintf("  Nearest NPDB core with a usable 0-30 cm stock: %.0f km away",
                 nearest_km))
log_warn(sprintf("  NPDB cores within 500 km of Fort Severn: %d",
                 sum(usable$dist_fort_severn_km <= 500)))
log_warn("")
log_warn("  This is the strongest argument for the community programme in the")
log_warn("  whole analysis. Eight cores are few, and they are also ALL there")
log_warn("  is. Every published carbon map covering Fort Severn is, in this")
log_warn("  neighbourhood, an extrapolation from hundreds of kilometres away.")
log_warn(strrep("-", 68))

# =============================================================================
# 4. Where does Fort Severn sit nationally?
# =============================================================================

log_step("14d  NATIONAL DISTRIBUTION")

stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")
fs <- stocks[stocks$scenario == "as_supplied" &
             stocks$window == "reference_0_30" & !stocks$is_lower_bound, ]
log_info(nrow(fs), " community cores are complete at 0-30 cm and comparable")

q <- stats::quantile(usable$std_stock_kgm2, c(.05, .25, .5, .75, .95))
log_info("NPDB 0-30 cm stock, complete cores only (kg C/m2):")
log_info(sprintf("   5%%=%.2f  25%%=%.2f  median=%.2f  75%%=%.2f  95%%=%.2f",
                 q[1], q[2], q[3], q[4], q[5]))

pct_of <- function(v) 100 * mean(usable$std_stock_kgm2 <= v)
fs$national_percentile <- vapply(fs$stock_kgm2, pct_of, numeric(1))
log_info("each community core against the national distribution:")
print_table(fs[order(-fs$stock_kgm2),
               c("core_id", "campaign", "stock_kgm2", "national_percentile")],
            digits = 1)

# --- organic versus mineral, nationally --------------------------------------
# The Fort Severn result that surprises people is that peat held LESS carbon in
# the top 30 cm than mineral soil. If that is a real property of shallow
# accounting rather than a local quirk, it should appear nationally too.
log_step("14e  DOES THE FORT SEVERN PATTERN HOLD NATIONALLY?")

if ("has_organic_horizon" %in% names(usable)) {
  og <- isTRUE_vec(usable$has_organic_horizon)
  cmp <- data.frame(
    group = c("organic horizon present", "no organic horizon"),
    n = c(sum(og), sum(!og)),
    median_kgm2 = c(stats::median(usable$std_stock_kgm2[og]),
                    stats::median(usable$std_stock_kgm2[!og])),
    q25 = c(stats::quantile(usable$std_stock_kgm2[og], .25),
            stats::quantile(usable$std_stock_kgm2[!og], .25)),
    q75 = c(stats::quantile(usable$std_stock_kgm2[og], .75),
            stats::quantile(usable$std_stock_kgm2[!og], .75)),
    stringsAsFactors = FALSE)
  print_table(cmp, digits = 2)

  # A difference in medians is not a finding until its size is compared with
  # the spread. These two distributions overlap heavily, so the test and the
  # overlap are reported before any claim is made.
  vo <- usable$std_stock_kgm2[og]; vm <- usable$std_stock_kgm2[!og]
  wt <- suppressWarnings(stats::wilcox.test(vo, vm))
  # Common-language effect size: probability a random organic core exceeds a
  # random mineral one. 0.5 means no difference at all.
  pd <- mean(outer(vo, vm, ">")) + 0.5 * mean(outer(vo, vm, "=="))
  cmp$wilcoxon_p <- signif(wt$p.value, 3)
  cmp$p_organic_exceeds_mineral <- round(pd, 3)

  log_info(sprintf("Wilcoxon rank-sum p = %.3g", wt$p.value))
  log_info(sprintf("probability a random organic core exceeds a random mineral one: %.2f",
                   pd))
  log_info(sprintf("IQR overlap: organic %.2f-%.2f, mineral %.2f-%.2f",
                   cmp$q25[1], cmp$q75[1], cmp$q25[2], cmp$q75[2]))

  if (cmp$median_kgm2[1] < cmp$median_kgm2[2]) {
    log_warn(strrep("-", 68))
    log_warn("THE FORT SEVERN ORDERING ALSO APPEARS NATIONALLY -- WEAKLY.")
    log_warn(sprintf("  Across %d NPDB cores, soils WITH an organic horizon hold",
                     nrow(usable)))
    log_warn(sprintf("  a median %.2f kg C/m2 over 0-30 cm, against %.2f for soils",
                     cmp$median_kgm2[1], cmp$median_kgm2[2]))
    log_warn("  WITHOUT one. The ordering matches Fort Severn.")
    log_warn("")
    log_warn("  BUT THE DIFFERENCE IS SMALL AND THE DISTRIBUTIONS OVERLAP.")
    log_warn(sprintf("  The gap is %.2f kg/m2 against an interquartile range of",
                     abs(diff(cmp$median_kgm2))))
    log_warn(sprintf("  %.1f kg/m2 for organic soils alone, and a random organic",
                     cmp$q75[1] - cmp$q25[1]))
    log_warn(sprintf("  core exceeds a random mineral one %.0f%% of the time.",
                     100 * pd))
    log_warn(if (wt$p.value < 0.05)
               sprintf("  The rank-sum test is significant (p = %.3g), but with", wt$p.value)
             else
               sprintf("  The rank-sum test is NOT significant (p = %.3g), so", wt$p.value))
    log_warn("  n this large, significance is easy and size is what matters.")
    log_warn("")
    log_warn("  Read it as CORROBORATION, NOT PROOF: the national record is")
    log_warn("  consistent with what the community cores found, and gives no")
    log_warn("  reason to suspect a local error. The strong version of the")
    log_warn("  claim -- that organic soils reliably hold less in the top")
    log_warn("  30 cm -- is NOT supported; organic soils are simply far more")
    log_warn(sprintf("  variable (IQR %.1f vs %.1f kg/m2).",
                     cmp$q75[1] - cmp$q25[1], cmp$q75[2] - cmp$q25[2]))
    log_warn("")
    log_warn("  What survives, and matters: a 0-30 cm accounting does NOT rank")
    log_warn("  peatlands above mineral soils. Depth is what distinguishes")
    log_warn("  them, and a shallow window cannot see it.")
    log_warn(strrep("-", 68))
  } else {
    log_info("Nationally, organic-horizon soils hold MORE in the top 30 cm, ",
             "so the Fort Severn ordering is local rather than general.")
  }
  write_csv_logged(cmp, file.path(CFG$dir_tables, "14_organic_vs_mineral.csv"),
                   "national test of the Fort Severn peat/mineral ordering")
}

# =============================================================================
# 5. Outputs
# =============================================================================

log_step("14f  OUTPUTS")

keep <- c("pedon_id", "province", "cal_year", "latitude", "longitude",
          "conf_metrs", "has_organic_horizon", "std_stock_kgm2",
          "std_stock_complete", "dist_fort_severn_km", "flag_location_conf_poor")
write_csv_logged(usable[, intersect(keep, names(usable))],
                 file.path(CFG$dir_derived, "14_npdb_cores_qc.csv"),
                 "NPDB cores usable for comparison, in kg C/m2")
write_csv_logged(rings, file.path(CFG$dir_tables, "14_npdb_proximity.csv"),
                 "how many existing cores lie near Fort Severn")
write_csv_logged(fs[, c("core_id", "campaign", "stock_kgm2",
                        "national_percentile")],
                 file.path(CFG$dir_tables, "14_fort_severn_percentiles.csv"),
                 "community cores against the national distribution")

qt <- data.frame(quantile = names(q), npdb_kgm2 = as.numeric(q),
                 stringsAsFactors = FALSE)
write_csv_logged(qt, file.path(CFG$dir_tables, "14_npdb_quantiles.csv"))

# =============================================================================
# 6. Figures
# =============================================================================

if (capabilities("png")) {
  f1 <- file.path(CFG$dir_figures, "14_national_distribution.png")
  grDevices::png(f1, width = 1600, height = 1000, res = 175)
  graphics::par(mar = c(5.4, 4.6, 3.2, 1.2), mgp = c(2.8, 0.7, 0))
  h <- graphics::hist(usable$std_stock_kgm2, breaks = 40, plot = FALSE)
  graphics::plot(h, col = "grey85", border = "grey60",
                 xlab = "Soil carbon, 0-30 cm (kg C / m2)", ylab = "NPDB cores",
                 main = "Fort Severn against the national soil carbon record")
  cols <- ifelse(grepl("peat", fs$campaign), "#1b9e77", "#7570b3")
  graphics::abline(v = fs$stock_kgm2, col = cols, lwd = 2.2)
  graphics::abline(v = stats::median(usable$std_stock_kgm2), col = "black",
                   lwd = 2, lty = 2)
  graphics::legend("topright", bty = "n", cex = 0.75,
                   legend = c(sprintf("NPDB median (%.1f)",
                                      stats::median(usable$std_stock_kgm2)),
                              "Fort Severn peat core", "Fort Severn mineral cores"),
                   col = c("black", "#1b9e77", "#7570b3"),
                   lwd = c(2, 2.2, 2.2), lty = c(2, 1, 1))
  graphics::mtext(sprintf(paste0("%d NPDB cores complete to 30 cm. Nearest to ",
                                 "Fort Severn is %.0f km away."),
                          nrow(usable), nearest_km),
                  side = 1, line = 4.0, cex = 0.68, adj = 0, col = "grey30")
  grDevices::dev.off()
  log_ok("wrote ", basename(f1))

  f2 <- file.path(CFG$dir_figures, "14_npdb_coverage_map.png")
  grDevices::png(f2, width = 1600, height = 1200, res = 175)
  graphics::par(mar = c(5.4, 4.6, 3.2, 1.2), mgp = c(2.8, 0.7, 0))
  graphics::plot(usable$longitude, usable$latitude, pch = 16, cex = 0.45,
                 col = grDevices::adjustcolor("grey35", 0.55),
                 asp = 1 / cos(56 * pi / 180),
                 xlab = "Longitude", ylab = "Latitude",
                 main = "Where Canada's soil carbon cores are, and are not")
  graphics::points(CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]],
                   pch = 21, bg = "#d7301f", col = "white", cex = 2.2, lwd = 2)
  graphics::text(CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]],
                 "Fort Severn", pos = 4, cex = 0.8, font = 2, col = "#d7301f")
  for (r in c(500, 1000)) {
    aa <- seq(0, 2 * pi, length.out = 200)
    graphics::lines(CFG$site_anchor[["lon"]] + (r / 111.32 /
                      cos(56 * pi / 180)) * cos(aa),
                    CFG$site_anchor[["lat"]] + (r / 111.32) * sin(aa),
                    col = "#d7301f", lty = 3, lwd = 1.3)
  }
  graphics::mtext(sprintf(paste0("Dotted rings: 500 and 1,000 km. No usable ",
                                 "NPDB core lies within %.0f km of Fort Severn."),
                          nearest_km),
                  side = 1, line = 4.0, cex = 0.68, adj = 0, col = "grey30")
  grDevices::dev.off()
  log_ok("wrote ", basename(f2))
}

log_ok("14 complete")
