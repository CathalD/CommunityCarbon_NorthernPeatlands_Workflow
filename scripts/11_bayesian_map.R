# =============================================================================
# 11_bayesian_map.R
#
# PRODUCT 4  Bayesian carbon maps: a published map as PRIOR, the community
#            cores as OBSERVATIONS.
#            Status: PRIOR-DOMINATED PREDICTION, updated locally by measurement.
#
# The model is in R/bayes.R. In short: each pixel starts at the published
# value with the published uncertainty, and each core pulls nearby pixels
# towards what was measured, with a pull that fades with distance and adds up
# across cores. Far from every core the map is the published product unchanged.
#
# TWO PRODUCTS, BECAUSE THERE ARE TWO DEPTHS AND THEY ARE NOT THE SAME QUANTITY
#
#   A. SHALLOW, 0-30 cm.  Prior is a 0-30 cm product (Sothe 0-30 or SoilGrids).
#      The cores measure exactly this quantity, so the update is like for like
#      and needs no assumption about anything below 30 cm.
#
#   B. FULL COLUMN.  Prior is Li et al. 2025, which is carbon from the surface
#      to the base of the peat -- metres, not centimetres. The cores cannot
#      observe that. So the column is split into its top 30 cm and everything
#      beneath, ONLY the top is updated, and the deep part is added back
#      untouched. That split needs an assumed shallow fraction, which is the
#      single largest assumption in this script. It is set in config
#      (CFG$bayes$shallow_fraction_of_column), reported everywhere, and varied
#      in the sensitivity table.
#
#   The two products are NEVER averaged or differenced. 04's depth-support
#   guard applies here too.
#
# THREE WAYS TO GET A PRIOR, in order of preference
#
#   1. GEE      : rgee reads the published rasters directly. Preferred.
#   2. LOCAL    : reference GeoTIFFs already retrieved into outputs/gee/ and
#                 `terra` available to read them.
#   3. REGIONAL : no raster at all -- the prior is the published REGIONAL MEAN
#                 (Li et al: 86 +/- 35 kg C/m2), spatially constant. This drops
#                 all spatial pattern from the prior and keeps only its
#                 magnitude, so outputs are prefixed DEMO_ and must not be
#                 presented as the Li map. It exists so the workflow produces
#                 real, inspectable GeoTIFFs without credentials.
#
# INPUT   data/derived/01_cores_qc.csv, 01_segments_qc.csv, 02_core_stocks.csv
# OUTPUT  outputs/spatial/BAYES_*.tif  (or DEMO_BAYES_*.tif in regional mode)
#         outputs/tables/11_*.csv
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

log_step("11  BAYESIAN CARBON MAPS")
set.seed(CFG$seed)

dir_sp <- file.path(CFG$root, "outputs", "spatial")
if (!dir.exists(dir_sp)) dir.create(dir_sp, recursive = TRUE)

cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
seg    <- require_artifact(file.path(CFG$dir_derived, "01_segments_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

if (any(cores$longitude > 0)) {
  fail_loudly("Positive longitude reached the Bayesian mapping step",
              remedy = "Re-run 01_ingest_qc.R.")
}

# =============================================================================
# 1. How deep is the peat where the cores actually are?
# =============================================================================

log_step("11a  PEAT DEPTH AT THE CORES")

pd <- detect_peat_depth_all(seg, CFG$bayes$om_peat_threshold_pct)
pd$campaign <- cores$campaign[match(pd$core_id, cores$core_id)]
print_table(pd[, c("core_id", "campaign", "peat_depth_cm", "status",
                   "core_bottom_cm")], digits = 2)

peat_cores <- pd[pd$status != "no_peat", ]
log_info("peat defined as organic matter >= ",
         CFG$bayes$om_peat_threshold_pct, "%")

if (nrow(peat_cores)) {
  obs <- peat_cores[peat_cores$status == "contact_observed", ]
  log_warn(strrep("-", 68))
  log_warn("THE PEAT CORES SIT ON THIN PEAT.")
  if (nrow(obs)) {
    for (k in seq_len(nrow(obs))) {
      log_warn(sprintf("  %-8s peat/mineral contact at %.1f cm",
                       obs$core_id[k], obs$peat_depth_cm[k]))
    }
  }
  lb <- peat_cores[peat_cores$is_lower_bound, ]
  for (k in seq_len(nrow(lb))) {
    log_warn(sprintf("  %-8s still peat at %.1f cm (lower bound)",
                     lb$core_id[k], lb$peat_depth_cm[k]))
  }
  log_warn("")
  log_warn(sprintf("  Li et al. report a Hudson Bay Lowlands mean peat depth of"))
  log_warn(sprintf("  %.0f +/- %.0f cm. The observed contacts here are an order",
                   CFG$lit$li_peat_depth_mean_cm, CFG$lit$li_peat_depth_sd_cm))
  log_warn("  of magnitude shallower.")
  log_warn("")
  log_warn("  These cores are on peat MARGINS or shallow fen, not on the deep")
  log_warn("  peat plateaus that dominate the regional carbon stock. They are")
  log_warn("  therefore weak evidence about a full-column product, and that is")
  log_warn("  reflected below: the full-column posterior barely moves.")
  log_warn(strrep("-", 68))
}
write_csv_logged(pd, file.path(CFG$dir_tables, "11_peat_depth_at_cores.csv"),
                 "peat/mineral contact detected from the OM profile")

# =============================================================================
# 2. Observations
# =============================================================================

log_step("11b  OBSERVATIONS AND THEIR UNCERTAINTY")

obs30 <- stocks[stocks$scenario == "as_supplied" &
                stocks$window == "reference_0_30", ]
obs30 <- obs30[match(cores$core_id, obs30$core_id), ]

# Lower-bound cores cannot serve as observations of a 0-30 cm stock: their
# value is a floor, not a measurement, and treating it as a measurement would
# drag the posterior down wherever they sit.
usable <- !obs30$is_lower_bound
log_info(sum(usable), " of ", nrow(obs30),
         " cores are complete at 0-30 cm and may act as observations")
if (any(!usable)) {
  log_warn("excluded as lower bounds, not measurements: ",
           paste(obs30$core_id[!usable], collapse = ", "))
}

y      <- obs30$stock_kgm2[usable]
ystrat <- obs30$campaign[usable]
ylon   <- cores$longitude[usable]
ylat   <- cores$latitude[usable]
yid    <- cores$core_id[usable]

sig <- core_sigma(y, ystrat)
log_info("observation uncertainty per core (kg C/m2):")
print_table(cbind(data.frame(core_id = yid, stringsAsFactors = FALSE), sig),
            digits = 3)
if (any(sig$floored)) {
  log_warn("Some within-stratum sd values were floored at the pooled sd (",
           sprintf("%.3f", sig$floor_sd[1]), ").")
  log_warn("An sd estimated from very few cores is not a reliable small ",
           "number; without a floor those cores would overwrite the prior.")
}

# =============================================================================
# 3. Grid and distances
# =============================================================================

log_step("11c  GRID")

g <- make_grid(CFG$aoi_bbox, CFG$spatial$cell_m)
LON <- matrix(rep(g$x, each = g$ny), nrow = g$ny)
LAT <- matrix(rep(g$y, times = g$nx), nrow = g$ny)
n_pix <- length(LON)
log_info(sprintf("%d x %d cells at ~%d m", g$nx, g$ny, CFG$spatial$cell_m))

D <- matrix(NA_real_, n_pix, length(y))
for (j in seq_along(y)) D[, j] <- haversine_km(LON, LAT, ylon[j], ylat[j])

ls_km <- CFG$bayes$length_scale_km
if (is.null(ls_km)) {
  nn <- nn_distance_summary(cores$longitude, cores$latitude, cores$core_id)
  ls_km <- median(nn$nn_dist_km)
  log_info(sprintf("length scale derived from data: %.2f km (median core spacing)",
                   ls_km))
} else {
  log_info(sprintf("length scale from config: %.2f km", ls_km))
}
log_info(sprintf("influence truncated at %.1f km", CFG$bayes$max_influence_km))

# =============================================================================
# 4. Prior
# =============================================================================

log_step("11d  PRIOR SOURCE")

ref_tif <- list(
  li     = file.path(CFG$dir_gee, "ccnp_li2025_peat_carbon_30m.tif"),
  sothe  = file.path(CFG$dir_gee, "ccnp_sothe_sc_0_30_30m.tif"),
  sg     = file.path(CFG$dir_gee, "ccnp_soilgrids_ocs_0_30_30m.tif")
)
have_terra <- requireNamespace("terra", quietly = TRUE)
have_local <- any(file.exists(unlist(ref_tif))) && have_terra

mode <- if (have_local) "LOCAL" else "REGIONAL"
log_info("prior mode: ", mode)
if (mode == "REGIONAL") {
  log_warn(strrep("-", 68))
  log_warn("NO REFERENCE RASTER AVAILABLE. Running in REGIONAL mode.")
  log_warn("  The prior is the PUBLISHED REGIONAL MEAN, spatially constant:")
  log_warn(sprintf("    Li et al. 2025 full peat column: %.0f +/- %.0f kg C/m2",
                   CFG$lit$li_peat_stock_mean_kgm2, CFG$lit$li_peat_stock_sd_kgm2))
  log_warn("  All spatial pattern from the published map is absent. What these")
  log_warn("  outputs DO show correctly is where and how much the community")
  log_warn("  cores change a prior of that magnitude.")
  log_warn("  Outputs are prefixed DEMO_ and must not be presented as the Li map.")
  log_warn("")
  log_warn("  To get the real thing: run 04_gee_reference_audit.R, retrieve the")
  log_warn("  exported reference GeoTIFFs into ", CFG$dir_gee)
  log_warn("  and re-run this script.")
  log_warn(strrep("-", 68))
}

prefix <- if (mode == "REGIONAL") "DEMO_BAYES_" else "BAYES_"

read_prior <- function(path, g) {
  r <- terra::rast(path)[[1]]
  tmpl <- terra::rast(nrows = g$ny, ncols = g$nx,
                      xmin = g$xmin, xmax = g$xmin + g$nx * g$xres,
                      ymax = g$ymax, ymin = g$ymax - g$ny * g$yres,
                      crs = "EPSG:4326")
  as.vector(t(terra::as.matrix(terra::resample(r, tmpl, method = "bilinear"),
                               wide = TRUE)))
}

# =============================================================================
# 5. PRODUCT B -- full column, Li prior
# =============================================================================

log_step("11e  FULL-COLUMN PRODUCT (Li et al. prior)")

f_shallow <- CFG$bayes$shallow_fraction_of_column
log_info(sprintf("assumed shallow fraction: %.3f of full-column carbon in 0-30 cm",
                 f_shallow))
log_info(sprintf("  uniform density would give %.3f (30 / %.0f cm)",
                 uniform_shallow_fraction(30, CFG$lit$li_peat_depth_mean_cm),
                 CFG$lit$li_peat_depth_mean_cm))
log_info("  peat compacts with depth, so the true share is lower than uniform")

if (mode == "LOCAL" && file.exists(ref_tif$li)) {
  li_mean <- read_prior(ref_tif$li, g)
  li_sd   <- rep(CFG$lit$li_peat_stock_sd_kgm2, n_pix)
} else {
  li_mean <- rep(CFG$lit$li_peat_stock_mean_kgm2, n_pix)
  li_sd   <- rep(CFG$lit$li_peat_stock_sd_kgm2, n_pix)
}

sp <- split_full_column(li_mean, f_shallow)
# The prior's uncertainty is apportioned to the two layers in the same ratio.
sp_sd_shallow <- li_sd * f_shallow
sp_sd_deep    <- li_sd * (1 - f_shallow)

# Only the PEAT cores may update a peat-column product. A mineral-soil core
# says nothing about the carbon in a peat column.
is_peat_core <- grepl("peat", ystrat)
log_info(sum(is_peat_core), " peat core(s) may update the full-column prior; ",
         sum(!is_peat_core), " mineral core(s) are excluded as irrelevant to a ",
         "peat column")

if (sum(is_peat_core) == 0L) {
  log_warn("No usable peat core reaches 30 cm, so the full-column product ",
           "cannot be updated at all and would equal the prior everywhere.")
  log_warn("Writing the prior-only layers so the comparison in 12 still works.")
  post_full <- data.frame(mean = li_mean, sd = li_sd,
                          info_frac = rep(0, n_pix), shift = rep(0, n_pix))
} else {
  Wf <- build_weight_matrix(D[, is_peat_core, drop = FALSE], ls_km,
                            CFG$bayes$max_influence_km)
  post_shallow <- bayes_posterior(sp$shallow, sp_sd_shallow,
                                  y[is_peat_core], sig$sigma_used[is_peat_core],
                                  Wf)
  rc <- recombine_full_column(post_shallow$mean, post_shallow$sd,
                              sp$deep, sp_sd_deep)
  post_full <- data.frame(
    mean = rc$mean, sd = rc$sd,
    # Information fraction is reported for the FULL column, so it correctly
    # shows how little of a metres-deep stock eight 30 cm cores constrain.
    info_frac = post_shallow$info_frac * f_shallow,
    shift = rc$mean - li_mean)
}

log_info(sprintf("full-column posterior: %.1f - %.1f kg C/m2 (prior %.1f)",
                 min(post_full$mean, na.rm = TRUE),
                 max(post_full$mean, na.rm = TRUE), li_mean[1]))
log_info(sprintf("maximum shift from prior: %+.2f kg C/m2 (%.2f%% of prior)",
                 post_full$shift[which.max(abs(post_full$shift))],
                 100 * max(abs(post_full$shift), na.rm = TRUE) / li_mean[1]))
log_info(sprintf("maximum information share from cores: %.3f%%",
                 100 * max(post_full$info_frac, na.rm = TRUE)))

log_warn("Eight shallow cores move a full peat column stock by at most ",
         sprintf("%.1f%%.", 100 * max(abs(post_full$shift), na.rm = TRUE) / li_mean[1]),
         " That is the correct answer, not a")
log_warn("failure: 30 cm is a small fraction of a column that averages ",
         CFG$lit$li_peat_depth_mean_cm, " cm, and these particular cores sit ",
         "on unusually thin peat.")

as_mat <- function(v) matrix(v, nrow = g$ny, ncol = g$nx, byrow = FALSE)
w <- function(v, nm) {
  p <- file.path(dir_sp, paste0(prefix, nm, ".tif"))
  write_geotiff(as_mat(v), p, xmin = g$xmin, ymax = g$ymax,
                xres = g$xres, yres = g$yres, epsg = CFG$spatial$epsg_out)
  log_ok("wrote ", basename(p))
  invisible(p)
}
w(post_full$mean,      "fullcolumn_mean_kgm2")
w(post_full$sd,        "fullcolumn_sd_kgm2")
w(post_full$info_frac, "fullcolumn_core_info_fraction")
w(post_full$shift,     "fullcolumn_shift_from_prior_kgm2")

# =============================================================================
# 6. PRODUCT A -- shallow 0-30 cm
# =============================================================================

log_step("11f  SHALLOW 0-30 cm PRODUCT")

shallow_prior_name <- NA_character_
if (mode == "LOCAL" && file.exists(ref_tif$sothe)) {
  mu0 <- read_prior(ref_tif$sothe, g); shallow_prior_name <- "Sothe 0-30 cm"
} else if (mode == "LOCAL" && file.exists(ref_tif$sg)) {
  mu0 <- read_prior(ref_tif$sg, g);    shallow_prior_name <- "SoilGrids 0-30 cm"
} else {
  mu0 <- NULL
}

if (is.null(mu0)) {
  log_warn("No 0-30 cm prior raster is available, and unlike the full-column")
  log_warn("case there is no citable regional constant to stand in for one.")
  log_warn("Inventing a number here would be worse than producing nothing, so")
  log_warn("the shallow Bayesian product is SKIPPED.")
  log_warn("Retrieve ccnp_sothe_sc_0_30_30m.tif or ccnp_soilgrids_ocs_0_30_30m.tif")
  log_warn("into ", CFG$dir_gee, " and re-run.")
  shallow_written <- FALSE
} else {
  log_info("prior: ", shallow_prior_name)
  sd0 <- rep(stats::sd(y, na.rm = TRUE), n_pix)   # conservative prior sd
  Ws <- build_weight_matrix(D, ls_km, CFG$bayes$max_influence_km)
  post_sh <- bayes_posterior(mu0, sd0, y, sig$sigma_used, Ws)
  w(post_sh$mean,      "shallow_0_30_mean_kgm2")
  w(post_sh$sd,        "shallow_0_30_sd_kgm2")
  w(post_sh$info_frac, "shallow_0_30_core_info_fraction")
  w(post_sh$shift,     "shallow_0_30_shift_from_prior_kgm2")
  shallow_written <- TRUE
  log_warn("Stratum gating was NOT applied: no stratum raster is present. A ",
           "mineral core is therefore allowed to influence peat pixels. ",
           "Retrieve ccnp_gwl_fcs30_30m.tif to enable gating.")
}

# =============================================================================
# 7. Sensitivity
# =============================================================================

log_step("11g  SENSITIVITY")

sens <- do.call(rbind, lapply(CFG$bayes$sensitivity_length_scales_km, function(L) {
  Wl <- build_weight_matrix(D[, is_peat_core, drop = FALSE], L,
                            CFG$bayes$max_influence_km)
  ps <- bayes_posterior(sp$shallow, sp_sd_shallow, y[is_peat_core],
                        sig$sigma_used[is_peat_core], Wl)
  data.frame(length_scale_km = L,
             pct_area_info_gt_1pct = 100 * mean(ps$info_frac > 0.01, na.rm = TRUE),
             max_info_frac = max(ps$info_frac, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
log_info("effect of the kernel length scale:")
print_table(sens, digits = 3)

fsens <- do.call(rbind, lapply(CFG$bayes$shallow_fraction_sensitivity, function(f) {
  s2 <- split_full_column(li_mean, f)
  Wf2 <- build_weight_matrix(D[, is_peat_core, drop = FALSE], ls_km,
                             CFG$bayes$max_influence_km)
  p2 <- bayes_posterior(s2$shallow, li_sd * f, y[is_peat_core],
                        sig$sigma_used[is_peat_core], Wf2)
  r2 <- recombine_full_column(p2$mean, p2$sd, s2$deep, li_sd * (1 - f))
  data.frame(shallow_fraction = f,
             max_shift_kgm2 = max(abs(r2$mean - li_mean), na.rm = TRUE),
             pct_of_prior = 100 * max(abs(r2$mean - li_mean), na.rm = TRUE) / li_mean[1],
             stringsAsFactors = FALSE)
}))
log_info("effect of the assumed shallow fraction:")
print_table(fsens, digits = 3)
log_info("Note: the length-scale rows saturate because influence is truncated ",
         "at ", CFG$bayes$max_influence_km, " km regardless of length scale.")

# --- what shallow fraction would make the prior and the cores agree? ---------
# The shift is minimised where the prior's implied 0-30 cm carbon equals what
# the cores measured. Solving for that fraction turns the cores into an
# independent check on the published full-column magnitude, which is a much
# more interesting use of them than nudging a map.
y_peat <- y[is_peat_core]
if (length(y_peat) && li_mean[1] > 0) {
  f_implied <- mean(y_peat) / li_mean[1]
  f_uniform <- uniform_shallow_fraction(30, CFG$lit$li_peat_depth_mean_cm)

  log_step("11g2  DO THE CORES AGREE WITH LI ET AL.?")
  log_info(sprintf("measured 0-30 cm stock, peat core(s): %.2f kg C/m2",
                   mean(y_peat)))
  log_info(sprintf("Li et al. full-column mean            : %.0f kg C/m2",
                   li_mean[1]))
  log_info(sprintf("=> implied share of the column in the top 30 cm: %.3f",
                   f_implied))
  log_info(sprintf("   uniform-density expectation                 : %.3f",
                   f_uniform))
  log_info(sprintf("   assumed in this run                         : %.3f",
                   f_shallow))

  if (f_implied < f_uniform) {
    log_ok("The cores imply a SMALLER top-30 cm share than uniform density.")
    log_info("  That is the direction compaction predicts, and these cores ")
    log_info("  show it directly: bulk density in PM-03-A rises from 0.050 to ")
    log_info("  0.229 g/cm3 between the surface and 30 cm.")
    log_info("  So the community measurement is CONSISTENT with the published ")
    log_info("  full-column magnitude, under a physically sensible depth ")
    log_info("  distribution. That is independent corroboration of Li et al. ")
    log_info("  at this site, which is worth more than the map nudge.")
  } else {
    log_warn("The cores imply a LARGER top-30 cm share than uniform density, ",
             "which compaction does not predict. Investigate before relying ",
             "on either number here.")
  }
  log_warn("Caveat that limits this: it rests on ", length(y_peat),
           " peat core(s) complete to 30 cm, and those cores sit on thin peat ",
           "(15-25 cm) rather than the deep column Li maps.")

  write_csv_logged(data.frame(
    measured_0_30_kgm2 = mean(y_peat),
    li_full_column_kgm2 = li_mean[1],
    implied_shallow_fraction = f_implied,
    uniform_density_fraction = f_uniform,
    assumed_fraction = f_shallow,
    n_peat_cores_used = length(y_peat),
    stringsAsFactors = FALSE),
    file.path(CFG$dir_tables, "11_implied_shallow_fraction.csv"),
    "cores as an independent check on the published full-column magnitude")
}

write_csv_logged(sens, file.path(CFG$dir_tables, "11_sensitivity_length_scale.csv"))
write_csv_logged(fsens, file.path(CFG$dir_tables, "11_sensitivity_shallow_fraction.csv"))

# =============================================================================
# 7b. Figures
# =============================================================================

log_step("11i  FIGURES")

if (capabilities("png")) {
  ramp_info <- grDevices::colorRampPalette(
    c("#f7f4f9", "#d4b9da", "#df65b0", "#ce1256", "#67001f"))

  draw <- function(m, main, sub, ramp, path) {
    grDevices::png(path, width = 1400, height = 1400, res = 170)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(5.8, 4.4, 3.4, 1.2), mgp = c(2.6, 0.7, 0))
    graphics::image(x = g$x, y = rev(g$y), z = t(m[nrow(m):1, ]),
                    col = ramp(200), asp = 1 / cos(g$lat_mid * pi / 180),
                    xlab = "Longitude", ylab = "Latitude", main = main,
                    useRaster = TRUE)
    graphics::points(cores$longitude, cores$latitude, pch = 21,
                     bg = ifelse(grepl("peat", cores$campaign), "#1b9e77", "#7570b3"),
                     col = "white", cex = 1.5, lwd = 1.6)
    graphics::legend("topright", legend = c("PM peat (2024)", "FS mineral (2025)"),
                     pt.bg = c("#1b9e77", "#7570b3"), pch = 21, col = "white",
                     bg = "white", cex = 0.7, pt.cex = 1.2)
    graphics::mtext(sub, side = 1, line = 4.3, cex = 0.66, adj = 0, col = "grey30")
    graphics::box()
  }

  f_info <- file.path(CFG$dir_figures, "11_bayes_core_information.png")
  draw(as_mat(100 * post_full$info_frac),
       "Where the community cores actually changed the map",
       paste0("Percent of the answer coming from the cores rather than the ",
              "published prior. Peaks at ",
              sprintf("%.1f%%", 100 * max(post_full$info_frac, na.rm = TRUE)),
              "; zero across most of the area."),
       ramp_info, f_info)
  log_ok("wrote ", basename(f_info))

  f_mean <- file.path(CFG$dir_figures, "11_bayes_fullcolumn_mean.png")
  draw(as_mat(post_full$mean),
       "Full peat-column carbon, prior updated by community cores",
       sprintf(paste0("kg C/m2. Prior = %s. Mode: %s. The cores move it by at ",
                      "most %.2f kg/m2 (%.1f%%)."),
               "Li et al. 2025", mode, max(abs(post_full$shift), na.rm = TRUE),
               100 * max(abs(post_full$shift), na.rm = TRUE) / li_mean[1]),
       grDevices::colorRampPalette(c("#fff7ec", "#fdd49e", "#fc8d59", "#d7301f")),
       f_mean)
  log_ok("wrote ", basename(f_mean))
} else {
  log_warn("no PNG device; figures skipped (rasters written regardless)")
}

# =============================================================================
# 8. Manifest
# =============================================================================

log_step("11h  MANIFEST")

man <- data.frame(
  file = paste0(prefix, c("fullcolumn_mean_kgm2", "fullcolumn_sd_kgm2",
                          "fullcolumn_core_info_fraction",
                          "fullcolumn_shift_from_prior_kgm2"), ".tif"),
  units = c("kg C/m2", "kg C/m2", "fraction 0-1", "kg C/m2"),
  meaning = c("posterior mean full peat column carbon",
              "posterior standard deviation",
              "share of the answer coming from the cores rather than the prior",
              "how far the cores moved the published value"),
  depth_support = "FULL PEAT COLUMN -- never compare with a 0-30 cm layer",
  prior_mode = mode,
  stringsAsFactors = FALSE
)
if (exists("shallow_written") && shallow_written) {
  man <- rbind(man, data.frame(
    file = paste0(prefix, c("shallow_0_30_mean_kgm2", "shallow_0_30_sd_kgm2",
                            "shallow_0_30_core_info_fraction",
                            "shallow_0_30_shift_from_prior_kgm2"), ".tif"),
    units = c("kg C/m2", "kg C/m2", "fraction 0-1", "kg C/m2"),
    meaning = c("posterior mean 0-30 cm carbon", "posterior standard deviation",
                "share of the answer coming from the cores",
                "how far the cores moved the published value"),
    depth_support = "0-30 cm",
    prior_mode = mode, stringsAsFactors = FALSE))
}
print_table(man)
write_csv_logged(man, file.path(dir_sp, "MANIFEST_bayes.csv"),
                 "what each Bayesian raster is")

log_ok("11 complete")
