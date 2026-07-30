# =============================================================================
# 01b_plot_profiles.R
#
# Two carbon-depth profile plots, so you can see the cores rather than read
# them off a table. Carbon on x, depth on y (increasing downward, as a soil
# profile is normally drawn).
#
#   RAW          carbon density per centimetre, drawn at each segment's ACTUAL
#                depth interval. Segment thicknesses differ between cores
#                (5.5 to 15.2 cm), so this is the honest picture of what was
#                measured, gaps and all.
#
#   HARMONISED   cumulative carbon stock down the profile, every core on the
#                same 0-30 cm axis. A core that stopped short is drawn DASHED
#                and simply ends -- it is a lower bound, and nothing is
#                extrapolated to make it look complete.
#
# Base R graphics on purpose: no ggplot2 dependency for two simple plots.
#
# INPUT   mvp/outputs/current/segments_clean.csv   (from 01)
# OUTPUT  mvp/outputs/current/figures/profiles_raw.png
#         mvp/outputs/current/figures/profiles_harmonised.png
#
# REQUIRES: nothing beyond base R
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

msg("01b  carbon profile plots")

seg_path <- file.path(CFG$dir_current, "segments_clean.csv")
if (!file.exists(seg_path)) stop("Missing ", seg_path, " -- run 01 first.")
seg <- read.csv(seg_path, stringsAsFactors = FALSE)
ensure_dir(CFG$dir_figures)

cores <- sort(unique(seg$core_id))
campaign_of <- seg$campaign[match(cores, seg$core_id)]

# One distinct colour per core: warm ramp for peat, cool for mineral. Colour
# carries core identity so that line TYPE stays free to mean exactly one thing
# in the harmonised plot -- "this core is a lower bound".
peat_ids    <- cores[campaign_of == "peat"]
mineral_ids <- cores[campaign_of == "mineral"]
pal <- setNames(rep(NA_character_, length(cores)), cores)
if (length(peat_ids)) {
  pal[peat_ids] <- colorRampPalette(c("#8C3A00", "#F2A24B"))(length(peat_ids))
}
if (length(mineral_ids)) {
  pal[mineral_ids] <- colorRampPalette(c("#12385F", "#69A8DC"))(length(mineral_ids))
}

legend_labels <- paste0(cores, "  (", campaign_of, ")")
max_depth <- max(seg$depth_bottom_cm)

# =============================================================================
# 1. RAW: carbon density per cm, at each segment's real depth interval
# =============================================================================

png(file.path(CFG$dir_figures, "profiles_raw.png"),
    width = 1500, height = 1100, res = 190)
par(mar = c(4.5, 4.5, 3.5, 1))

plot(NA, xlim = c(0, max(seg$carbon_density_kgm2_per_cm) * 1.05),
     ylim = c(max_depth, 0),   # reversed: depth increases downward
     xlab = expression("Carbon density  (kg C m"^-2*" per cm of depth)"),
     ylab = "Depth below surface (cm)",
     main = "Raw core profiles\ndrawn at each segment's measured depth interval")
grid(col = "grey90", lty = 1)

for (cc in cores) {
  s <- seg[seg$core_id == cc, ]
  s <- s[order(s$depth_top_cm), ]
  # One horizontal bar per segment: constant density across its own interval.
  for (i in seq_len(nrow(s))) {
    lines(rep(s$carbon_density_kgm2_per_cm[i], 2),
          c(s$depth_top_cm[i], s$depth_bottom_cm[i]),
          col = pal[[cc]], lwd = 2.6)
  }
  # Horizontal connectors at each segment boundary, so a core reads as one
  # profile rather than a set of floating bars.
  if (nrow(s) > 1) {
    for (i in seq_len(nrow(s) - 1)) {
      lines(c(s$carbon_density_kgm2_per_cm[i], s$carbon_density_kgm2_per_cm[i + 1]),
            rep(s$depth_bottom_cm[i], 2),
            col = pal[[cc]], lwd = 1.1, lty = 3)
    }
  }
}

legend("bottomright", legend = legend_labels, col = pal[cores],
       lty = 1, lwd = 2.4, bty = "n", cex = 0.75)
dev.off()
msg("wrote figures/profiles_raw.png")

# =============================================================================
# 2. HARMONISED: cumulative stock, common 0-30 cm axis
# =============================================================================

ref <- CFG$reference_depth_cm

png(file.path(CFG$dir_figures, "profiles_harmonised.png"),
    width = 1500, height = 1100, res = 190)
par(mar = c(4.5, 4.5, 3.5, 1))

cum_max <- max(tapply(seg$stock_kgm2_segment, seg$core_id, sum))
plot(NA, xlim = c(0, cum_max * 1.05), ylim = c(ref, 0),
     xlab = expression("Cumulative carbon stock  (kg C m"^-2*")"),
     ylab = "Depth below surface (cm)",
     main = paste0("Harmonised profiles: cumulative carbon to ", ref, " cm",
                   "\ndashed = core stopped short, so its total is a LOWER BOUND"))
grid(col = "grey90", lty = 1)

for (cc in cores) {
  s <- seg[seg$core_id == cc, ]
  s <- s[order(s$depth_top_cm), ]
  # Step from (0,0) down through each segment boundary.
  depth <- c(0, s$depth_bottom_cm)
  cumc  <- c(0, cumsum(s$stock_kgm2_segment))
  short <- max(s$depth_bottom_cm) < (ref - 0.1)
  lines(cumc, depth, col = pal[[cc]], lwd = 2.6, lty = if (short) 2 else 1)
  points(utils::tail(cumc, 1), utils::tail(depth, 1),
         col = pal[[cc]], pch = if (short) 4 else 19, cex = 1.1)
}

abline(h = ref, col = "grey45", lty = 3)

is_short <- vapply(cores, function(cc)
  max(seg$depth_bottom_cm[seg$core_id == cc]) < (ref - 0.1), logical(1))

legend("bottomright",
       legend = c(paste0(legend_labels, ifelse(is_short, "  LOWER BOUND", "")),
                  "solid + dot = reached reference depth",
                  "dashed + x  = stopped short"),
       col = c(pal[cores], "grey35", "grey35"),
       lty = c(ifelse(is_short, 2, 1), 1, 2),
       pch = c(ifelse(is_short, 4, 19), NA, NA),
       lwd = 2.4, bty = "n", cex = 0.72)
dev.off()
msg("wrote figures/profiles_harmonised.png")

msg("01b complete")
