# =============================================================================
# 11_context_figures.R      ADD-ON. Steps 01-07 do not depend on this.
#
# Six side-by-side comparison figures placing the Fort Severn cores against the
# four external databases.
#
#   1  ecosystem panels, one line per study   -- do the databases agree?
#   2  paired box plots by study              -- is mineral > organic real?
#   3  bulk density vs depth                  -- WHY peat holds less shallow C
#   4  cumulative carbon + share at 30 cm     -- what 30 cm misses, MEASURED
#   5  0-30 cm vs full column, two panels     -- why our numbers look small
#   6  distance vs stock                      -- is any nearby analogue valid?
#
# Figure 4 also writes shallow_fraction_measured.csv. Step 06 currently ASSUMES
# CFG$bayes$shallow_fraction_of_column = 0.15; CanPeat's full-depth profiles can
# estimate that empirically, so this figure feeds a real number back into the
# mapping pipeline rather than only describing the data.
#
# Figure 1 needs step 09; the other five do not.
#
# INPUT   mvp/outputs/current/external_profiles.csv, external_layers.csv  (08)
#         mvp/outputs/current/ecosystem_classes.csv                       (09, optional)
#         mvp/outputs/current/cores_clean.csv, segments_clean.csv         (01)
# OUTPUT  mvp/outputs/current/figures/context_*.png
#         mvp/outputs/current/shallow_fraction_measured.csv
#
# REQUIRES: readr, dplyr
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))
source(file.path(.this_dir, "plot_helpers.R"))

library(readr)
library(dplyr)

msg("11  context comparison figures")

ensure_dir(CFG$dir_figures)
REF <- CFG$reference_depth_cm

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------

need <- c("external_profiles.csv", "external_layers.csv",
         "cores_clean.csv", "segments_clean.csv")
miss <- need[!file.exists(file.path(CFG$dir_current, need))]
if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "), " -- run 01 and 08 first.")

prof <- read_csv(file.path(CFG$dir_current, "external_profiles.csv"), show_col_types = FALSE)
lay  <- read_csv(file.path(CFG$dir_current, "external_layers.csv"), show_col_types = FALSE)
ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE)
oseg <- read_csv(file.path(CFG$dir_current, "segments_clean.csv"), show_col_types = FALSE)

eco_path <- file.path(CFG$dir_current, "ecosystem_classes.csv")
have_eco <- file.exists(eco_path)
if (have_eco) {
  eco <- read_csv(eco_path, show_col_types = FALSE)
  prof <- left_join(prof, eco, by = c("dataset", "profile_id"))
}

lay <- lay %>%
  mutate(dens_kgm2_per_cm = ifelse(is.finite(layer_thickness_cm) & layer_thickness_cm > 0,
                                  stock_kgm2_layer / layer_thickness_cm, NA_real_)) %>%
  left_join(prof %>% select(dataset, profile_id, soil_type,
                           any_of(c("gwl_label", "wc_label"))),
           by = c("dataset", "profile_id"))

ours_type <- ifelse(ours$campaign == "peat", "organic", "mineral")
xmax_dens <- quantile(lay$dens_kgm2_per_cm, 0.995, na.rm = TRUE)

# =============================================================================
# FIGURE 1  Ecosystem panels, one line per study
# =============================================================================

if (have_eco && "gwl_label" %in% names(lay)) {

  # WHICH CLASS LAYER TO PANEL BY, decided from coverage rather than preference.
  #
  # GWL_FCS30 is the better layer in principle -- it is a wetland map, and ESA
  # WorldCover famously cannot tell treed bog on peat from upland spruce on
  # mineral soil. But it only returns a class where it has coverage, and on this
  # point set that turned out to be a minority of profiles. Panelling by a layer
  # that is missing for most of the data would produce panels built from a
  # handful of profiles and imply agreement that was never tested.
  cov_gwl <- mean(!is.na(prof$gwl_label))
  cov_wc  <- if ("wc_label" %in% names(prof)) mean(!is.na(prof$wc_label)) else 0
  msg(sprintf("class coverage: GWL_FCS30 %.0f%%, WorldCover %.0f%%",
             100 * cov_gwl, 100 * cov_wc))
  class_col <- if (cov_gwl >= 0.5 || cov_gwl >= cov_wc) "gwl_label" else "wc_label"
  class_name <- if (class_col == "gwl_label") "wetland class (GWL_FCS30)" else
                  "land cover (ESA WorldCover)"
  if (cov_gwl < 0.5 && class_col == "wc_label") {
    msg("GWL_FCS30 covers under half the profiles, so the panels use WorldCover ",
       "instead. Note the tradeoff: WorldCover cannot separate treed bog on ",
       "peat from upland forest on mineral soil, which is the distinction this ",
       "landscape turns on -- read the panels with that in mind.")
  }

  classes <- names(sort(table(prof[[class_col]]), decreasing = TRUE))
  classes <- head(classes[!is.na(classes)], 6)
  our_eco <- read_csv(eco_path, show_col_types = FALSE) %>%
    filter(dataset == "Fort Severn")

  fig("context_1_ecosystem_panels", {
    op <- par(mfrow = c(2, 3), mar = c(4.0, 4.0, 3.0, 0.8), oma = c(0, 0, 3.0, 0))
    on.exit(par(op), add = TRUE)
    for (cl in classes) {
      sub <- lay[!is.na(lay[[class_col]]) & lay[[class_col]] == cl, ]
      env_ds <- envelope(sub, "dataset")
      env_ds <- env_ds[!vapply(env_ds, is.null, logical(1))]
      plot(NA, xlim = c(0, xmax_dens), ylim = c(200, 0), xlab = "", ylab = "",
           main = sprintf("%s\n(%d profiles)", cl, sum(prof[[class_col]] == cl, na.rm = TRUE)),
           cex.main = 0.95)
      grid(col = "grey93", lty = 1)
      for (ds in names(env_ds)) draw_env(env_ds[[ds]], DS_COL[[ds]], lwd = 2.2)
      # Our cores appear only in the panel matching their own class.
      if (nrow(our_eco) && cl %in% our_eco[[class_col]]) {
        keep <- our_eco$profile_id[!is.na(our_eco[[class_col]]) & our_eco[[class_col]] == cl]
        draw_our_cores(oseg[oseg$core_id %in% keep, ])
      }
      if (!length(env_ds)) text(xmax_dens / 2, 100, "fewer than 3 profiles",
                               cex = 0.8, col = "grey45")
    }
    mtext(sprintf("Carbon profiles by %s, one colour per study\nbands = 10th-90th percentile   x = kg C/m2 per cm   y = depth (cm)",
                  class_name),
          outer = TRUE, cex = 0.9, font = 2)
  }, width = 2100, height = 1350)

  fig("context_1b_study_legend", {
    par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", legend = names(DS_COL), col = unname(DS_COL), lwd = 3,
           bty = "n", cex = 1.2, title = "Study")
  }, width = 700, height = 500)
} else {
  msg("skipped figure 1 (needs step 09)")
}

# =============================================================================
# FIGURE 2  Paired box plots, organic vs mineral within each study
# =============================================================================
# Only profiles that actually reach 30 cm: a lower bound is a floor, not a
# measurement, and mixing them in would bias every box downward.

box_df <- bind_rows(
  prof %>% filter(reaches_30cm) %>%
    transmute(dataset, soil_type, stock = stock_kgm2_0_30),
  tibble(dataset = "Fort Severn", soil_type = ours_type,
        stock = ours$stock_kgm2)[!ours$stock_is_lower_bound, ]
) %>% filter(soil_type %in% c("organic", "mineral"), is.finite(stock))

studies <- c("Fort Severn", "CanPeat", "NPDB", "Janousek", "WOSIS Canada")
studies <- studies[studies %in% box_df$dataset]

fig("context_2_paired_boxplots", {
  par(mar = c(6.2, 4.8, 4.0, 1))
  dat <- list(); at <- c(); cols <- c(); i <- 0
  for (s in studies) {
    for (tp in c("mineral", "organic")) {
      v <- box_df$stock[box_df$dataset == s & box_df$soil_type == tp]
      i <- i + 1
      dat[[i]] <- if (length(v)) v else NA_real_
      at <- c(at, which(studies == s) * 3 + ifelse(tp == "mineral", -0.45, 0.45))
      cols <- c(cols, COL[[tp]])
    }
  }
  bx <- boxplot(dat, at = at, col = fade("#FFFFFF", 1), border = cols,
               outline = FALSE, xaxt = "n", boxwex = 0.75, las = 1,
               ylab = expression("Carbon stock, 0-30 cm  (kg C m"^-2*")"),
               main = paste0("Is mineral soil really above peat over 0-30 cm?",
                             "\nsame question asked of every database"))
  grid(col = "grey93", lty = 1)
  boxplot(dat, at = at, col = vapply(cols, fade, character(1), 0.22),
         border = cols, outline = FALSE, xaxt = "n", boxwex = 0.75, add = TRUE)
  axis(1, at = seq_along(studies) * 3, labels = studies, las = 2, cex.axis = 0.88)
  # n under each box
  ns <- vapply(dat, function(v) sum(is.finite(v)), integer(1))
  mtext(paste0("n=", ns), side = 1, at = at, line = 0.15, cex = 0.6, col = "grey35")
  legend("topright", legend = c("mineral", "organic"),
         fill = c(fade(COL[["mineral"]], 0.22), fade(COL[["organic"]], 0.22)),
         border = c(COL[["mineral"]], COL[["organic"]]), bty = "n", cex = 0.85)
})

# The number the figure is testing, printed so it can be quoted.
med_cmp <- box_df %>% group_by(dataset, soil_type) %>%
  summarise(n = n(), median_0_30 = round(median(stock), 2), .groups = "drop")
msg("--- median 0-30 cm stock by study and soil type ---")
print(as.data.frame(med_cmp))

# =============================================================================
# FIGURE 3  Bulk density vs depth -- the mechanism
# =============================================================================
# Log x axis: our peat cores run 0.05-0.27 g/cm3 against 0.75-2.18 for mineral,
# an order of magnitude, and a linear axis would compress the peat end to a line.

bd <- lay %>%
  filter(is.finite(BDOD), BDOD > 0, is.finite(upper_depth),
        soil_type %in% c("organic", "mineral")) %>%
  mutate(depth_mid = (upper_depth + lower_depth) / 2,
        bin = cut(depth_mid, breaks = seq(0, 300, by = 10), labels = FALSE))

bd_sum <- bd %>%
  group_by(soil_type, bin) %>%
  summarise(depth = mean(depth_mid, na.rm = TRUE), n = n(),
           p25 = quantile(BDOD, .25, na.rm = TRUE),
           p50 = median(BDOD, na.rm = TRUE),
           p75 = quantile(BDOD, .75, na.rm = TRUE), .groups = "drop") %>%
  filter(n >= 10)

fig("context_3_bulk_density", {
  par(mar = c(4.8, 4.8, 4.2, 1))
  plot(NA, xlim = c(0.02, 2.5), ylim = c(200, 0), log = "x",
       xlab = expression("Dry bulk density  (g cm"^-3*", log scale)"),
       ylab = "Depth below surface (cm)",
       main = paste0("Why peat holds less carbon in the top 30 cm",
                     "\nhigh carbon % x low bulk density = little carbon mass"))
  grid(col = "grey93", lty = 1)
  for (tp in c("organic", "mineral")) {
    e <- bd_sum[bd_sum$soil_type == tp, ]
    if (!nrow(e)) next
    polygon(c(e$p25, rev(e$p75)), c(e$depth, rev(e$depth)),
            col = fade(COL[[tp]], 0.18), border = NA)
    lines(e$p50, e$depth, col = COL[[tp]], lwd = 2.8)
  }
  # Our cores
  for (cc in unique(oseg$core_id)) {
    s <- oseg[oseg$core_id == cc, ]; s <- s[order(s$depth_top_cm), ]
    for (i in seq_len(nrow(s))) {
      lines(rep(s$bulk_density_gcm3[i], 2),
            c(s$depth_top_cm[i], s$depth_bottom_cm[i]),
            col = COL[["ours"]], lwd = 2.0)
    }
  }
  abline(h = REF, col = "grey40", lty = 3)
  legend("bottomleft",
         legend = c("organic (external, IQR)", "mineral (external, IQR)",
                    "Fort Severn cores"),
         col = c(COL[["organic"]], COL[["mineral"]], COL[["ours"]]),
         lwd = c(2.8, 2.8, 2.0), bty = "n", cex = 0.8)
})

# =============================================================================
# FIGURE 4  Cumulative carbon, and the share sitting above 30 cm
# =============================================================================
# This measures CFG$bayes$shallow_fraction_of_column, which step 06 currently
# assumes to be 0.15.

group_col <- if (have_eco && "gwl_label" %in% names(lay)) "gwl_label" else "soil_type"
cum <- cumulative_by_group(lay, group_col, ref = REF)
cum <- cum[!vapply(cum, is.null, logical(1))]
# Keep the best-supported groups so the figure stays legible.
cum <- cum[order(vapply(cum, function(x) -x$n_profiles, numeric(1)))]
cum <- head(cum, 6)
pal_cum <- setNames(hcl.colors(length(cum), "Dark 3"), names(cum))

fig("context_4_cumulative_carbon", {
  par(mar = c(4.8, 4.8, 4.4, 1))
  xm <- max(vapply(cum, function(x) max(x$curve$p50, na.rm = TRUE), numeric(1)), na.rm = TRUE)
  plot(NA, xlim = c(0, xm * 1.02), ylim = c(200, 0),
       xlab = expression("Cumulative carbon from the surface  (kg C m"^-2*")"),
       ylab = "Depth below surface (cm)",
       main = paste0("What a 0-", REF, " cm survey captures, and what it misses",
                     "\nmedian cumulative carbon; % = share of the full column above ", REF, " cm"))
  grid(col = "grey93", lty = 1)
  for (nm in names(cum)) {
    cv <- cum[[nm]]$curve
    cv <- cv[is.finite(cv$p50) & cv$n >= 3, ]
    lines(cv$p50, cv$depth, col = pal_cum[[nm]], lwd = 2.8)
  }
  abline(h = REF, col = "grey25", lty = 2, lwd = 1.4)
  labs <- vapply(names(cum), function(nm)
    sprintf("%s  (%d profiles, %.0f%% above %d cm)", nm, cum[[nm]]$n_profiles,
            100 * cum[[nm]]$share_at_ref, REF), character(1))
  legend("bottomright", legend = labs, col = pal_cum[names(cum)],
         lwd = 2.8, bty = "n", cex = 0.72)
  # Our cores' own total, as a reference tick on the 30 cm line.
  our_tot <- sum(ours$stock_kgm2[!ours$stock_is_lower_bound]) /
             sum(!ours$stock_is_lower_bound)
  points(our_tot, REF, pch = 23, bg = COL[["ours"]], col = "white", cex = 1.6, lwd = 1.5)
  text(our_tot, REF - 9, sprintf("Fort Severn mean\n%.1f kg C/m2 at %d cm", our_tot, REF),
       cex = 0.68, col = "grey20")
})

shallow_tbl <- tibble(
  group = names(cum),
  n_profiles = vapply(cum, function(x) x$n_profiles, numeric(1)),
  measured_share_above_ref = round(vapply(cum, function(x) x$share_at_ref, numeric(1)), 4),
  ref_depth_cm = REF,
  config_assumption = CFG$bayes$shallow_fraction_of_column
)
write_csv(shallow_tbl, file.path(CFG$dir_current, "shallow_fraction_measured.csv"))
msg("--- share of the full column above ", REF, " cm (MEASURED) ---")
print(as.data.frame(shallow_tbl))
msg("config currently assumes CFG$bayes$shallow_fraction_of_column = ",
   CFG$bayes$shallow_fraction_of_column,
   " -- compare against the peat-bearing rows above and update if they disagree.")

# =============================================================================
# FIGURE 5  Same groups, two depth windows
# =============================================================================

f5 <- prof %>% filter(soil_type %in% c("organic", "mineral"))
fig("context_5_two_windows", {
  op <- par(mfrow = c(1, 2), mar = c(5.6, 4.6, 4.2, 0.8), oma = c(0, 0, 2.2, 0))
  for (w in c("stock_kgm2_0_30", "stock_kgm2_total")) {
    d <- f5[f5$reaches_30cm | w == "stock_kgm2_total", ]
    vals <- split(d[[w]], d$soil_type)
    boxplot(vals, outline = FALSE, las = 1,
            col = vapply(names(vals), function(n) fade(COL[[n]], 0.22), character(1)),
            border = vapply(names(vals), function(n) COL[[n]], character(1)),
            ylab = expression("kg C m"^-2),
            main = if (w == "stock_kgm2_0_30") paste0("0-", REF, " cm only") else "full column")
    grid(col = "grey93", lty = 1)
    # Our cores are 0-30 cm in BOTH panels -- that is the entire point of the
    # figure, so the same values are plotted twice against different y scales.
    xs <- match(ours_type, names(vals))
    points(jitter(xs, amount = 0.12), ours$stock_kgm2,
           pch = 19, col = COL[["ours"]], cex = 1.2)
  }
  par(op)
  mtext(paste0("The same soils, two accounting depths -- our cores are 0-", REF,
               " cm in BOTH panels\nblack points = Fort Severn cores"),
        outer = TRUE, cex = 0.92, font = 2)
}, width = 1500, height = 1000)

# =============================================================================
# FIGURE 6  Distance from Fort Severn vs stock
# =============================================================================

f6 <- prof %>% filter(reaches_30cm, is.finite(stock_kgm2_0_30), stock_kgm2_0_30 > 0,
                     soil_type %in% c("organic", "mineral"))

nearest <- min(f6$dist_fort_severn_km, na.rm = TRUE)

fig("context_6_distance_vs_stock", {
  par(mar = c(4.8, 4.8, 4.2, 1))
  plot(pmax(f6$dist_fort_severn_km, 1), f6$stock_kgm2_0_30, log = "x",
       pch = 21, cex = 0.7, lwd = 0.4,
       bg = vapply(f6$soil_type, function(t) fade(COL[[t]], 0.45), character(1)),
       col = NA,
       xlab = "Distance from Fort Severn (km, log scale)",
       ylab = expression("Carbon stock, 0-30 cm  (kg C m"^-2*")"),
       main = paste0("There is no nearby analogue to borrow from",
                     sprintf("\nnearest of any soil type %.0f km; nearest MINERAL %s; and no trend with distance beyond",
                             nearest,
                             if (any(f6$soil_type == "mineral"))
                               sprintf("%.0f km", min(f6$dist_fort_severn_km[f6$soil_type == "mineral"], na.rm = TRUE))
                             else "none")))
  grid(col = "grey93", lty = 1)
  rect(1, par("usr")[3], nearest, 10^par("usr")[4],
       col = fade("#C1660A", 0.07), border = NA)
  abline(v = nearest, col = "#C1660A", lty = 2, lwd = 1.6)
  text(sqrt(1 * nearest), max(f6$stock_kgm2_0_30, na.rm = TRUE) * 0.94,
       sprintf("no comparable\nmeasurement\nwithin %.0f km", nearest),
       cex = 0.76, col = "#8C3A00")
  pts <- ours$stock_kgm2[!ours$stock_is_lower_bound]
  points(rep(1, length(pts)), pts, pch = 23, bg = COL[["ours"]],
         col = "white", cex = 1.4, lwd = 1.4)
  legend("topright", legend = c("organic", "mineral", "Fort Severn cores"),
         pch = c(21, 21, 23), pt.bg = c(COL[["organic"]], COL[["mineral"]], COL[["ours"]]),
         col = NA, bty = "n", cex = 0.8)
})

msg("11 complete")
