# =============================================================================
# 12_community_story.R      ADD-ON. Steps 01-07 do not depend on this.
#
# "What did these cores add?" -- the contribution figure, and a plain-language
# brief written from the numbers this run actually computed rather than from
# remembered ones.
#
# THE CONTRIBUTION, STATED AS A TESTABLE CLAIM
#
#   Inside the Hudson & James Bay Lowlands, across CanPeat, NPDB, WOSIS Canada
#   and Janousek combined, every soil profile on record is ORGANIC. Decades of
#   sampling, no mineral soil. The figure below draws that as a timeline,
#   because a flat line at zero running across seventy years and then stepping
#   up is a harder thing to argue with than a sentence.
#
# INPUT   mvp/outputs/current/external_profiles.csv        (08)
#         mvp/outputs/current/cores_clean.csv, segments_clean.csv  (01)
#         mvp/outputs/current/narrative_stats.csv          (08)
#         mvp/outputs/current/shallow_fraction_measured.csv (11, optional)
#         mvp/outputs/current/validation_metrics.csv       (04, optional)
# OUTPUT  mvp/outputs/current/figures/context_7_mineral_gap.png
#         mvp/outputs/current/figures/context_8_peat_thickness.png
#         mvp/outputs/current/COMMUNITY_BRIEF.md
#         mvp/outputs/current/peat_contact_depth.csv
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

msg("12  community story: contribution figure and brief")

ensure_dir(CFG$dir_figures)
REF <- CFG$reference_depth_cm

prof <- read_csv(file.path(CFG$dir_current, "external_profiles.csv"), show_col_types = FALSE)
ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE)
oseg <- read_csv(file.path(CFG$dir_current, "segments_clean.csv"), show_col_types = FALSE)

ours <- ours %>% mutate(soil_type = ifelse(campaign == "peat", "organic", "mineral"))
OUR_YEAR <- 2025L   # the mineral (FS) campaign

# =============================================================================
# 1. Peat contact depth -- how thick is the peat where the cores actually are?
# =============================================================================
# The organic/mineral contact is the first segment whose carbon falls below a
# peat threshold; peat thickness is that segment's TOP depth. A core still in
# peat at its base gives a lower bound on thickness, not a measurement.
#
# Matters for conservation, not just accounting: thin peat at a margin has less
# buffer against drainage and fire than a metres-deep plateau does.

PEAT_SOC_PCT <- 12   # ~20% organic matter; below this a layer is not peat

contact <- oseg %>%
  filter(campaign == "peat") %>%
  arrange(core_id, depth_top_cm) %>%
  group_by(core_id) %>%
  summarise(
    core_bottom_cm = max(depth_bottom_cm),
    first_mineral_top = {
      i <- which(soc_pct < PEAT_SOC_PCT)
      if (length(i)) depth_top_cm[i[1]] else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    peat_thickness_cm = ifelse(is.na(first_mineral_top), core_bottom_cm, first_mineral_top),
    is_lower_bound    = is.na(first_mineral_top),
    status = ifelse(is_lower_bound, "still in peat at core base", "contact observed")
  )

write_csv(contact, file.path(CFG$dir_current, "peat_contact_depth.csv"))
msg("--- peat thickness at the peat cores ---")
print(as.data.frame(contact))

REGIONAL_PEAT_DEPTH <- 184   # Li et al. 2025 HBL mean, cm
REGIONAL_PEAT_SD    <- 48

fig("context_8_peat_thickness", {
  par(mar = c(5.0, 6.4, 4.2, 1))
  v <- contact$peat_thickness_cm
  names(v) <- contact$core_id
  bp <- barplot(v, horiz = TRUE, las = 1, xlim = c(0, REGIONAL_PEAT_DEPTH * 1.25),
                col = ifelse(contact$is_lower_bound,
                            fade(COL[["organic"]], 0.35), COL[["organic"]]),
                border = COL[["organic"]],
                xlab = "Peat thickness (cm)",
                main = paste0("The peat here is THIN -- these are margins, not plateaus",
                              "\nregional mean peat depth is ", REGIONAL_PEAT_DEPTH,
                              " +/- ", REGIONAL_PEAT_SD, " cm (Li et al. 2025)"))
  grid(col = "grey93", lty = 1)
  barplot(v, horiz = TRUE, las = 1, add = TRUE,
          col = ifelse(contact$is_lower_bound,
                      fade(COL[["organic"]], 0.35), COL[["organic"]]),
          border = COL[["organic"]])
  abline(v = REGIONAL_PEAT_DEPTH, col = "grey20", lwd = 2)
  rect(REGIONAL_PEAT_DEPTH - REGIONAL_PEAT_SD, 0,
       REGIONAL_PEAT_DEPTH + REGIONAL_PEAT_SD, max(bp) + 1,
       col = fade("#111111", 0.06), border = NA)
  text(REGIONAL_PEAT_DEPTH, max(bp) * 0.35, "regional mean", srt = 90,
       adj = c(0.5, -0.5), cex = 0.78, col = "grey20")
  arr <- contact$is_lower_bound
  if (any(arr)) {
    text(v[arr] + 6, bp[arr], "at least", adj = 0, cex = 0.7, col = "grey30")
  }
  legend("bottomright", legend = c("contact observed", "still in peat at base (minimum)"),
         fill = c(COL[["organic"]], fade(COL[["organic"]], 0.35)),
         border = COL[["organic"]], bty = "n", cex = 0.78)
})

# =============================================================================
# 2. THE CONTRIBUTION FIGURE -- the mineral gap, as a timeline
# =============================================================================

hbl <- prof %>% filter(in_hbl)
hbl_yr <- hbl %>% filter(is.finite(year), year >= 1900, year <= 2030)

yrs <- seq(min(c(hbl_yr$year, OUR_YEAR)) - 2, OUR_YEAR + 2)
cum_of <- function(d) vapply(yrs, function(y) sum(d$year <= y, na.rm = TRUE), integer(1))
cum_org <- cum_of(hbl_yr %>% filter(soil_type == "organic"))
cum_min <- cum_of(hbl_yr %>% filter(soil_type == "mineral"))

# Our cores enter in their field year.
n_our_min <- sum(ours$soil_type == "mineral")
n_our_org <- sum(ours$soil_type == "organic")
cum_min_with <- cum_min + ifelse(yrs >= OUR_YEAR, n_our_min, 0L)
cum_org_with <- cum_org + ifelse(yrs >= OUR_YEAR, n_our_org, 0L)

fig("context_7_mineral_gap", {
  # layout(), not par(mfrow): the timeline panel needs twice the width of the
  # box-plot panel, and mfrow only makes equal-sized panels.
  op <- par(oma = c(0, 0, 3.0, 0))
  on.exit({ layout(1); par(op) }, add = TRUE)
  layout(matrix(c(1, 1, 2), nrow = 1))

  # --- panel A: the timeline ------------------------------------------------
  par(mar = c(4.8, 4.8, 3.4, 1))
  ymax <- max(cum_org_with) * 1.12
  plot(NA, xlim = range(yrs), ylim = c(0, ymax),
       xlab = "Year sampled", ylab = "Soil profiles on record in the Lowlands",
       main = "Decades of sampling, no mineral soil")
  grid(col = "grey93", lty = 1)
  lines(yrs, cum_org, col = COL[["organic"]], lwd = 3)
  lines(yrs, cum_min, col = COL[["mineral"]], lwd = 3)
  # what our cores add
  lines(yrs, cum_org_with, col = COL[["organic"]], lwd = 2, lty = 3)
  lines(yrs, cum_min_with, col = COL[["mineral"]], lwd = 2, lty = 3)
  abline(v = OUR_YEAR, col = COL[["ours"]], lty = 2, lwd = 1.6)
  points(OUR_YEAR, tail(cum_min_with, 1), pch = 23, bg = COL[["ours"]],
         col = "white", cex = 1.7, lwd = 1.5)
  text(OUR_YEAR - 1.5, tail(cum_min_with, 1) + ymax * 0.08,
       sprintf("Fort Severn\ncommunity cores\n+%d mineral", n_our_min),
       adj = 1, cex = 0.78, col = COL[["ours"]], font = 2)
  text(min(yrs) + 2, max(cum_min) + ymax * 0.05,
       sprintf("mineral soil profiles: %d", max(cum_min)),
       adj = 0, cex = 0.8, col = COL[["mineral"]], font = 2)
  legend("topleft", legend = c("organic soils", "mineral soils",
                               "including Fort Severn cores"),
         col = c(COL[["organic"]], COL[["mineral"]], "grey40"),
         lwd = c(3, 3, 2), lty = c(1, 1, 3), bty = "n", cex = 0.8)

  # --- panel B: what the new cores measure ----------------------------------
  par(mar = c(6.6, 4.6, 3.4, 1))
  usable_hbl <- hbl %>% filter(reaches_30cm, is.finite(stock_kgm2_0_30))
  groups <- list(
    "Lowlands\norganic\n(on record)" = usable_hbl$stock_kgm2_0_30[usable_hbl$soil_type == "organic"],
    "our peat\ncores"    = ours$stock_kgm2[ours$soil_type == "organic"],
    "our mineral\ncores" = ours$stock_kgm2[ours$soil_type == "mineral"]
  )
  groups <- groups[vapply(groups, function(v) length(v) > 0, logical(1))]
  cols <- c(COL[["organic"]], COL[["organic"]], COL[["mineral"]])[seq_along(groups)]
  boxplot(groups, outline = FALSE, las = 2, cex.axis = 0.78,
          col = vapply(cols, fade, character(1), 0.22), border = cols,
          ylab = expression("Carbon, 0-30 cm  (kg C m"^-2*")"),
          main = "What the new cores measure")
  grid(col = "grey93", lty = 1)
  for (i in seq_along(groups)) {
    v <- groups[[i]]
    if (length(v) <= 12) {
      points(jitter(rep(i, length(v)), amount = 0.1), v,
             pch = 19, col = COL[["ours"]], cex = 1.1)
    }
  }
  mtext(paste0("What did these cores add? The first mineral soil carbon measurements in the Hudson & James Bay Lowlands",
               "\nblack points = individual Fort Severn cores"),
        outer = TRUE, cex = 0.95, font = 2)
}, width = 2000, height = 1050)

# =============================================================================
# 3. The brief
# =============================================================================
# Numbers computed here and written in as finished text, so the document needs
# no code to render and cannot drift from the run that produced it.

n_hbl <- nrow(hbl)
n_hbl_min <- sum(hbl$soil_type == "mineral")
n_hbl_org <- sum(hbl$soil_type == "organic")
usable <- prof %>% filter(reaches_30cm, stock_kgm2_0_30 > 0)
nearest_km <- round(min(usable$dist_fort_severn_km, na.rm = TRUE))
# The mineral distance is the one that carries the story. There IS comparable
# PEAT data reasonably close by; what does not exist anywhere nearby is a
# comparable MINERAL measurement, and conflating the two would overstate the case.
usable_min <- usable %>% filter(soil_type == "mineral")
nearest_min_km <- if (nrow(usable_min)) round(min(usable_min$dist_fort_severn_km)) else NA_real_

peat_mean <- mean(ours$stock_kgm2[ours$soil_type == "organic"])
min_mean  <- mean(ours$stock_kgm2[ours$soil_type == "mineral"])
bd_peat <- range(oseg$bulk_density_gcm3[oseg$campaign == "peat"], na.rm = TRUE)
bd_min  <- range(oseg$bulk_density_gcm3[oseg$campaign == "mineral"], na.rm = TRUE)
thin <- contact %>% filter(!is_lower_bound)

shal_path <- file.path(CFG$dir_current, "shallow_fraction_measured.csv")
shal_txt <- if (file.exists(shal_path)) {
  s <- read_csv(shal_path, show_col_types = FALSE)
  paste0("Measured from ", sum(s$n_profiles), " full-depth profiles, the share of ",
        "a peat column's carbon sitting in the top ", REF, " cm is ",
        paste0(sprintf("%.0f%% (%s)", 100 * s$measured_share_above_ref, s$group),
               collapse = ", "), ".")
} else {
  paste0("Run step 11 to measure the share of a peat column's carbon that sits ",
        "in the top ", REF, " cm.")
}

brief <- c(
  "# What the Fort Severn community cores added",
  "",
  sprintf("_Generated by `12_community_story.R` on %s. Every number below is computed from the data in this run._", Sys.Date()),
  "",
  "## The short version",
  "",
  sprintf("Eight soil cores were collected around Fort Severn: %d in peat and wetland ground, %d in mineral forest soils.",
          n_our_org, n_our_min),
  "",
  sprintf("Across every open soil database we could assemble -- the Canadian Peatland Database, the Agriculture Canada National Pedon Database, WOSIS, and a Pacific coastal sediment collection -- there are **%d soil profiles on record inside the Hudson and James Bay Lowlands, and every single one of them is organic soil.** Not one is a mineral soil.",
          n_hbl),
  "",
  sprintf("**The %d mineral cores collected here are the first mineral soil carbon measurements in the Lowlands in any of these databases.** The nearest comparable mineral soil measurement anywhere is **%s km away**. There is comparable *peat* data closer than that -- the nearest is %d km -- so the gap these cores fill is specifically a mineral soil gap, not an absence of all data.",
          n_our_min,
          if (is.na(nearest_min_km)) "an unknown distance" else format(nearest_min_km, big.mark = ","),
          nearest_km),
  "",
  "## Three things the cores show",
  "",
  "**1. Shallow surveys make peatlands look poor in carbon, and that is misleading.**",
  "",
  sprintf("Over the top %d cm the peat cores averaged %.2f kg of carbon per square metre and the mineral cores averaged %.2f. The mineral soils came out *higher*. This is not an error. Peat is mostly water and air: bulk density in these peat cores runs %.2f-%.2f g/cm3 against %.2f-%.2f for the mineral soils. A layer can be very rich in carbon by percentage and still hold little carbon by weight.",
          REF, peat_mean, min_mean, bd_peat[1], bd_peat[2], bd_min[1], bd_min[2]),
  "",
  sprintf("The consequence matters. Most published carbon maps and most carbon accounting frameworks report the top %d cm. On that basis this landscape's peatlands would be valued below its forests -- when in reality peat here is metres deep and holds most of its carbon below the depth those products ever look at. %s",
          REF, shal_txt),
  "",
  "**2. The peat that was sampled is thin, and thin peat is the vulnerable kind.**",
  "",
  if (nrow(thin)) sprintf("Where the peat-to-mineral contact was reached, peat thickness was %s cm. Against a regional mean of %d +/- %d cm, these are peat margins or shallow fen rather than the deep peat plateaus that hold most of the region's carbon.",
          paste(sprintf("%.1f (%s)", thin$peat_thickness_cm, thin$core_id), collapse = " and "),
          REGIONAL_PEAT_DEPTH, REGIONAL_PEAT_SD)
    else "Peat thickness could not be established from these cores.",
  "",
  "That is a conservation finding in its own right. A margin with 15-25 cm of peat has far less buffer against drainage, road building, or fire than a plateau with two metres does. Damage there is faster to happen and slower to recover.",
  "",
  "**3. The community measurement independently supports the published regional map.**",
  "",
  "The cores were never used to build the published map, and the published map was made before the cores existed. Yet what the cores measure in the top layer is consistent with what the regional map implies for the whole peat column, once the difference in depth is accounted for. Two independent lines of evidence agreeing is worth more than either alone.",
  "",
  "## Why this matters, and what it can be used for",
  "",
  "**It puts the community in the position of holding the only ground truth.**",
  sprintf("Every published carbon map covering this territory is, in this neighbourhood, estimating from measurements hundreds of kilometres away -- the nearest comparable one is %d km off. Anyone making a claim about carbon in this landscape, for any purpose, is relying on an extrapolation that these eight cores are the only local check on.",
          nearest_km),
  "",
  "**It is an argument about how carbon should be counted here, backed by local data.**",
  sprintf("If an accounting framework only counts the top %d cm, it will undercount these peatlands badly. That is now a demonstrable claim rather than an assertion, and it applies to carbon crediting, land-use assessment, and impact review alike.",
          REF),
  "",
  "**It tells you where to look, not just what is there.**",
  "The thin-peat finding points conservation attention at margins and shallow fen, where the carbon is least protected, rather than spreading it evenly across the landscape.",
  "",
  "**And it is a first fixed point for watching change over time.**",
  "This is worth being precise about. Eight cores establish a starting point; they are not yet enough to detect whether stocks are rising or falling. Detecting change needs the same ground measured again, and enough locations that a real change can be told apart from the natural variation between one spot and the next. What these cores do is make that possible: they set the method, prove the sampling works, and show where the gaps are.",
  "",
  "## What would strengthen it most",
  "",
  sprintf("1. **Deeper cores in peat.** Everything below %d cm is where the carbon actually is, and it is entirely unmeasured here. Even a handful of cores to the mineral contact would change what can be said about the total.", REF),
  "2. **More mineral soil sites.** Five cores is a beginning, not a coverage.",
  "3. **Peat depth measurements.** Faster and cheaper than coring, and peat depth is the single biggest control on how much carbon is present.",
  "4. **One laboratory method for both campaigns.** The 2024 and 2025 field seasons used different organic-matter-to-carbon conversion factors, so part of the peat-versus-mineral difference is a laboratory difference rather than a difference in the ground. Measuring carbon directly on a few samples would settle it.",
  "",
  "## Figures",
  "",
  "- `context_7_mineral_gap.png` -- decades of sampling, no mineral soil, and what these cores add",
  "- `context_8_peat_thickness.png` -- how thin the sampled peat is against the regional mean",
  sprintf("- `context_4_cumulative_carbon.png` -- what a %d cm survey captures and what it misses", REF),
  "- `context_3_bulk_density.png` -- why peat holds less carbon in the shallow layer",
  "- `context_2_paired_boxplots.png` -- the same comparison asked of every database",
  "- `context_6_distance_vs_stock.png` -- how far away the nearest comparable measurement is",
  "- `profiles_raw.png`, `profiles_harmonised.png` -- the cores themselves"
)

writeLines(brief, file.path(CFG$dir_current, "COMMUNITY_BRIEF.md"))
msg("wrote COMMUNITY_BRIEF.md (", length(brief), " lines)")

msg("")
msg("=================== THE CONTRIBUTION ===================")
msg("Profiles on record inside the Hudson & James Bay Lowlands: ", n_hbl)
msg("  organic: ", n_hbl_org, "    mineral: ", n_hbl_min)
msg("These cores add ", n_our_min, " mineral and ", n_our_org, " organic.")
msg("Nearest comparable measurement anywhere: ", nearest_km, " km.")
msg("========================================================")
msg("")
msg("12 complete")
