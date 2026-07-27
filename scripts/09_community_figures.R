# =============================================================================
# 09_community_figures.R
#
# Create a small, shareable figure pack and plain-language community brief around
# four questions:
#   1. What is the carbon stock around Fort Severn?
#   2. What did the new community core data teach us?
#   3. How do these results compare with other peatlands / Canada / global values?
#   4. What does this mean for conservation?
#
# Base R only. These are communication products, not new statistical products.
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

log_step("09  COMMUNITY FIGURES AND MAP BRIEF")

cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")
est    <- require_artifact(file.path(CFG$dir_derived, "05_stratum_estimates.csv"),
                           "scripts/05_stratified_estimate.R")
lookup <- require_artifact(file.path(CFG$dir_tables, "05_product2_map_lookup.csv"),
                           "scripts/05_stratified_estimate.R")

# ---- helpers ----------------------------------------------------------------
png_open <- function(file, width = 1800, height = 1200) {
  grDevices::png(file.path(CFG$dir_figures, file), width = width, height = height,
                 res = 180)
  par(family = "sans", mar = c(5, 5, 4, 1) + 0.1, las = 1)
}
finish <- function() invisible(grDevices::dev.off())
stock_row <- function(window, scenario = "as_supplied") {
  stocks[stocks$window == window & stocks$scenario == scenario, ]
}
fmt <- function(x, d = 1) format(round(x, d), nsmall = d, trim = TRUE)

# 02_core_stocks.csv already carries campaign/latitude/longitude. Keep that
# table as the source of plotting data so merge suffixes cannot accidentally
# create latitude.x/latitude.y and leave common$latitude as NULL. If an older
# derived table lacks the location columns, add them explicitly from 01_cores_qc.
common <- stock_row("common_support")
loc_cols <- c("campaign", "latitude", "longitude")
missing_loc <- setdiff(loc_cols, names(common))
if (length(missing_loc)) {
  common <- merge(
    common,
    cores[, c("core_id", missing_loc), drop = FALSE],
    by = "core_id", all.x = TRUE, sort = FALSE
  )
}
ref30  <- stock_row("reference_0_30")
stratum_common <- est[est$scenario == "as_supplied" & est$window == "common_support", ]

need_finite_columns <- function(df, cols, context) {
  bad <- cols[!vapply(cols, function(nm) {
    nm %in% names(df) && any(is.finite(df[[nm]]))
  }, logical(1))]
  if (length(bad)) {
    fail_loudly(
      paste0("Cannot draw ", context),
      paste0("Missing or non-finite columns: ", paste(bad, collapse = ", ")),
      remedy = "Re-run scripts/01_ingest_qc.R and scripts/02_depth_harmonise.R, then re-run 09_community_figures.R."
    )
  }
  invisible(TRUE)
}
need_finite_columns(common, c("longitude", "latitude", "stock_kgm2"),
                    "community figure pack")

peat_mean <- stratum_common$mean[stratum_common$stratum == "PM_2024_peat"]
min_mean  <- stratum_common$mean[stratum_common$stratum == "FS_2025_mineral"]
peat_se   <- stratum_common$se[stratum_common$stratum == "PM_2024_peat"]
min_se    <- stratum_common$se[stratum_common$stratum == "FS_2025_mineral"]

# Literature / reference comparison values, intentionally broad and labelled.
# They are context anchors for communication, not calibration data.
compare <- data.frame(
  label = c("Fort Severn peat cores\n0-14.5 cm",
            "Fort Severn mineral cores\n0-14.5 cm",
            "Fort Severn observed cores\n0-30 cm where complete",
            "SoilGrids/Sothe-style\nshallow products",
            "Hudson Bay Lowlands peat\nfull peat column (Li 2025)",
            "Northern peatlands globally\nfull peat column"),
  stock_kgm2 = c(peat_mean, min_mean,
                 mean(ref30$stock_kgm2[!ref30$is_lower_bound]),
                 mean(c(3, 14)),
                 86,
                 50),
  note = c("new cores", "new cores", "6 complete cores", "typical shallow range in this workflow", "published HBL mean", "order-of-magnitude context"),
  stringsAsFactors = FALSE
)
# ---- Figure 1: core map ------------------------------------------------------
png_open("community_01_core_locations.png")
cols <- ifelse(common$campaign == "PM_2024_peat", "#2b8cbe", "#7f3b08")
plot(common$longitude, common$latitude, pch = 21, bg = cols, col = "white", cex = 2.2,
     xlab = "Longitude", ylab = "Latitude",
     main = "Community soil cores near Fort Severn")
text(common$longitude, common$latitude, labels = common$core_id, pos = 3, cex = 0.75)
points(CFG$site_anchor[["lon"]], CFG$site_anchor[["lat"]], pch = 8, cex = 2, lwd = 2)
legend("bottomleft", pch = c(21, 21, 8), pt.bg = c("#2b8cbe", "#7f3b08", NA),
       col = c("white", "white", "black"), pt.cex = c(1.6, 1.6, 1.4), bty = "n",
       legend = c("peat / wetland cores", "mineral / forest cores", "Fort Severn"))
mtext("Eight cores in two clusters: useful for careful summaries, not enough for kriging.",
      side = 1, line = 4, las = 0, cex = 0.85)
finish()

# ---- Figure 2: new core stocks ----------------------------------------------
png_open("community_02_core_carbon_stocks.png")
o <- order(common$campaign, common$stock_kgm2)
bar_cols <- ifelse(common$campaign[o] == "PM_2024_peat", "#2b8cbe", "#7f3b08")
bp <- barplot(common$stock_kgm2[o], names.arg = common$core_id[o], col = bar_cols,
              ylab = "Soil carbon stock (kg C / m²)",
              main = "What the new cores measured in the shared 0-14.5 cm layer",
              ylim = c(0, max(common$stock_kgm2) * 1.25), cex.names = 0.8)
abline(h = c(peat_mean, min_mean), col = c("#2b8cbe", "#7f3b08"), lwd = 2, lty = 2)
legend("topleft", bty = "n", fill = c("#2b8cbe", "#7f3b08"),
       legend = c(paste0("peat mean = ", fmt(peat_mean), " kg C/m²"),
                  paste0("mineral mean = ", fmt(min_mean), " kg C/m²")))
mtext("The peat cores have high carbon %, but low bulk density; most peat carbon is deeper than this shallow window.",
      side = 1, line = 4, las = 0, cex = 0.85)
finish()

# ---- Figure 3: comparison ladder --------------------------------------------
png_open("community_03_compare_context.png", height = 1400)
y <- seq_len(nrow(compare))
plot(compare$stock_kgm2, y, pch = 21, bg = "#4daf4a", col = "white", cex = 2,
     yaxt = "n", xlab = "Soil carbon stock (kg C / m²)", ylab = "",
     main = "How the Fort Severn shallow cores compare with broader carbon contexts",
     xlim = c(0, max(compare$stock_kgm2) * 1.15))
axis(2, at = y, labels = compare$label, las = 1, cex.axis = 0.78)
abline(v = c(10, 50, 86), col = "grey75", lty = 3)
text(compare$stock_kgm2, y, labels = paste0(fmt(compare$stock_kgm2), "  ", compare$note),
     pos = 4, cex = 0.8)
mtext("Do not read shallow 0-30 cm products as whole-peatland carbon maps.",
      side = 1, line = 4, las = 0, cex = 0.85)
finish()

# ---- Figure 4: conservation message -----------------------------------------
png_open("community_04_conservation_messages.png", width = 1800, height = 1000)
plot.new(); par(xpd = NA)
text(0.5, 0.92, "Conservation meaning for Fort Severn peatlands", cex = 1.8, font = 2)
boxes <- c("Keep peat wet\nDrainage or drying can turn deep stored carbon into emissions.",
           "Avoid shallow-map mistakes\n0-30 cm maps miss most peat carbon and can undervalue peatlands.",
           "Use community data as evidence\nThese cores are among the first local peat/mineral soil samples from this area.",
           "Sample deeper next\nPeat depth to mineral contact is the key missing conservation metric.")
x <- c(0.25, 0.75, 0.25, 0.75); yb <- c(0.65, 0.65, 0.32, 0.32)
for (i in seq_along(boxes)) {
  rect(x[i] - 0.21, yb[i] - 0.12, x[i] + 0.21, yb[i] + 0.12,
       col = c("#c7e9b4", "#fdd49e", "#bdd7e7", "#dadaeb")[i], border = NA)
  text(x[i], yb[i], boxes[i], cex = 1.05)
}
finish()

# ---- plain-language brief ----------------------------------------------------
brief <- file.path(CFG$root, "outputs", "COMMUNITY_BRIEF.md")
lines <- c(
  "# Community figure brief: Fort Severn soil carbon",
  "",
  paste0("Generated by `09_community_figures.R` on ", format(Sys.time(), "%Y-%m-%d %H:%M"), "."),
  "",
  "## 1. What is the carbon stock around Fort Severn?",
  "",
  paste0("In the shared shallow layer measured by every core (0-14.5 cm), peat/wetland cores averaged **", fmt(peat_mean, 2), " kg C/m²** and mineral/forest cores averaged **", fmt(min_mean, 2), " kg C/m²**."),
  paste0("For the 0-30 cm reference layer, only ", sum(!ref30$is_lower_bound), " of ", nrow(ref30), " cores reached the full depth; the two shorter cores are lower bounds, not complete 0-30 cm measurements."),
  "",
  "## 2. What did the new data teach us?",
  "",
  "These are some of the first local peat and mineral soil samples from this area. They show that high peat carbon percentage does not automatically mean high shallow carbon mass: low peat bulk density keeps shallow stocks modest, while the major peatland carbon store sits deeper in the peat column.",
  "",
  "## 3. How do these results compare to other areas?",
  "",
  "The Fort Severn core values are shallow stocks and should be compared with other shallow products only when depth support matches. Full Hudson Bay Lowlands peat-column products are much larger (Li et al. 2025 reports an HBL mean near 86 kg C/m²) because they include metres of peat, not just the top 30 cm.",
  "",
  "## 4. Conservation perspective",
  "",
  "The conservation message is that peatlands can be undervalued by shallow maps. Protecting hydrology and avoiding disturbance matters because most carbon is below the sampled window. The next most valuable data would be deeper peat cores or peat-depth measurements to mineral contact.",
  "",
  "## Shareable files",
  "",
  "- `outputs/figures/community_01_core_locations.png`",
  "- `outputs/figures/community_02_core_carbon_stocks.png`",
  "- `outputs/figures/community_03_compare_context.png`",
  "- `outputs/figures/community_04_conservation_messages.png`"
)
writeLines(lines, brief)
log_ok("wrote COMMUNITY_BRIEF.md and 4 PNG figures")
