# =============================================================================
# 10_comparison_outputs.R      ADD-ON. Steps 01-07 do not depend on this.
#
# The comparison deliverables: profile figures, the data-gap figures, and a
# multi-layer GeoPackage for mapping.
#
# Runs with or without step 09. If ecosystem_classes.csv is absent the
# ecosystem-grouped figure is skipped and everything else is produced, so you
# never have to authenticate Earth Engine just to see the profile comparison.
#
# READABILITY CHOICE WORTH KNOWING ABOUT
#
#   Drawing 11,547 external profiles as individual lines is unreadable. So
#   external datasets are drawn as a 10th-90th percentile ENVELOPE plus a
#   median line per group, and the eight Fort Severn cores are drawn as
#   individual bold lines on top. Same information, legible. Set
#   SHOW_RAW_LINES <- TRUE below to draw every external profile faintly
#   instead.
#
# INPUT   mvp/outputs/current/external_profiles.csv, external_layers.csv  (08)
#         mvp/outputs/current/ecosystem_classes.csv                       (09, optional)
#         mvp/outputs/current/cores_clean.csv, segments_clean.csv         (01)
# OUTPUT  mvp/outputs/current/figures/compare_*.png
#         mvp/outputs/current/comparison_cores.gpkg   (5 layers)
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

msg("10  comparison figures, GeoPackage and narrative")

E <- CFG$ext
SHOW_RAW_LINES <- FALSE   # TRUE = draw every external profile faintly

ensure_dir(CFG$dir_figures)

#' Open a PNG, run the plotting code, close the device even if it errors.
fig <- function(name, expr, width = 1600, height = 1150, res = 185) {
  png(file.path(CFG$dir_figures, paste0(name, ".png")),
      width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  msg("wrote figures/", name, ".png")
}

#' Same colour at a given opacity. Defined here rather than relying on base R's
#' colour-adjust helper, so there is nothing to get wrong about its name or
#' argument signature across R versions.
fade <- function(col, alpha) {
  v <- grDevices::col2rgb(col) / 255
  grDevices::rgb(v[1], v[2], v[3], alpha = alpha)
}

# Palette. Organic warm, mineral cool, our cores black so they always read as
# the subject rather than one series among thousands.
COL <- c(organic = "#C1660A", mineral = "#1F5FA8", unknown = "#8A8A8A",
        ours = "#111111")
DS_COL <- c("CanPeat" = "#C1660A", "NPDB" = "#1F5FA8",
           "Janousek" = "#2E8B6F", "WOSIS Canada" = "#8B5FA8",
           "Fort Severn" = "#111111")

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------

need <- c("external_profiles.csv", "external_layers.csv",
         "cores_clean.csv", "segments_clean.csv")
miss <- need[!file.exists(file.path(CFG$dir_current, need))]
if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "),
                      " -- run steps 01 and 08 first.")

prof <- read_csv(file.path(CFG$dir_current, "external_profiles.csv"), show_col_types = FALSE)
lay  <- read_csv(file.path(CFG$dir_current, "external_layers.csv"), show_col_types = FALSE)
ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE)
oseg <- read_csv(file.path(CFG$dir_current, "segments_clean.csv"), show_col_types = FALSE)

eco_path <- file.path(CFG$dir_current, "ecosystem_classes.csv")
have_eco <- file.exists(eco_path)
eco <- if (have_eco) read_csv(eco_path, show_col_types = FALSE) else NULL
if (have_eco) {
  prof <- left_join(prof, eco, by = c("dataset", "profile_id"))
  msg("ecosystem classes joined")
} else {
  msg("no ecosystem_classes.csv -- run 09 to add the ecosystem figure; ",
     "everything else still builds")
}

# Our cores on the external schema, so they can share axes and the GeoPackage.
ours_std <- ours %>%
  transmute(dataset = "Fort Severn", profile_id = core_id,
           latitude, longitude, year = NA_integer_,
           n_layers = n_segments, total_depth_cm = core_depth_cm,
           stock_kgm2_total = stock_kgm2, stock_kgm2_0_30 = stock_kgm2,
           reaches_30cm = !stock_is_lower_bound,
           stock_0_30_is_lower_bound = stock_is_lower_bound,
           dist_fort_severn_km = 0,
           soil_type = ifelse(campaign == "peat", "organic", "mineral"),
           soil_type_basis = "field campaign (peat vs mineral transect)",
           in_hbl = TRUE)
if (have_eco) {
  ours_std <- left_join(ours_std, eco, by = c("dataset", "profile_id"))
}

# -----------------------------------------------------------------------------
# 2. Profile envelopes
# -----------------------------------------------------------------------------
# Carbon density per cm on a common 1 cm depth grid, so profiles with unequal
# layer thickness are directly comparable. Same quantity plotted in 01b.

grid_cm <- 0:300

#' Step-interpolate a profile's carbon density onto the 1 cm grid. Returns NA
#' below the profile's own bottom -- profiles are never extended past what was
#' measured.
densify <- function(upper, lower, dens, grid = grid_cm) {
  out <- rep(NA_real_, length(grid))
  ok <- is.finite(upper) & is.finite(lower) & is.finite(dens) & lower > upper
  for (i in which(ok)) {
    sel <- grid >= upper[i] & grid < lower[i]
    out[sel] <- dens[i]
  }
  out
}

lay <- lay %>%
  mutate(dens_kgm2_per_cm = ifelse(is.finite(layer_thickness_cm) & layer_thickness_cm > 0,
                                  stock_kgm2_layer / layer_thickness_cm, NA_real_))

#' Per-group percentile envelope across profiles, on the common grid.
#'
#' Uses split() rather than filtering inside a loop: NPDB alone is 9,017
#' profiles over 48,372 layers, and a per-profile subset scan would be
#' quadratic and take minutes.
envelope <- function(layer_df, group_col) {
  groups <- unique(layer_df[[group_col]])
  groups <- groups[!is.na(groups)]
  lapply(setNames(groups, groups), function(g) {
    sub <- layer_df[!is.na(layer_df[[group_col]]) & layer_df[[group_col]] == g, ]
    by_prof <- split(sub, sub$profile_id)
    if (length(by_prof) < 3) return(NULL)
    m <- vapply(by_prof, function(s) {
      densify(s$upper_depth, s$lower_depth, s$dens_kgm2_per_cm)
    }, numeric(length(grid_cm)))
    qs <- function(p) apply(m, 1, quantile, probs = p, na.rm = TRUE, names = FALSE)
    data.frame(
      depth = grid_cm,
      n     = apply(m, 1, function(v) sum(is.finite(v))),
      p10   = qs(0.10), p50 = qs(0.50), p90 = qs(0.90)
    )
  })
}

#' Draw one envelope. Only where at least `min_n` profiles still contribute, so
#' the band does not narrow to a spurious line at depth as profiles drop out.
draw_env <- function(e, col, min_n = 3) {
  if (is.null(e)) return(invisible())
  e <- e[is.finite(e$p50) & e$n >= min_n, ]
  if (!nrow(e)) return(invisible())
  polygon(c(e$p10, rev(e$p90)), c(e$depth, rev(e$depth)),
          col = fade(col, 0.18), border = NA)
  lines(e$p50, e$depth, col = col, lwd = 2.6)
}

# -----------------------------------------------------------------------------
# FIGURE 1: full-depth profiles, organic vs mineral
# -----------------------------------------------------------------------------

lay_typed <- lay %>%
  left_join(prof %>% select(dataset, profile_id, soil_type), by = c("dataset", "profile_id"))
env_type <- envelope(lay_typed, "soil_type")

xmax <- quantile(lay$dens_kgm2_per_cm, 0.995, na.rm = TRUE)

fig("compare_profiles_full", {
  par(mar = c(4.6, 4.6, 3.8, 1))
  plot(NA, xlim = c(0, xmax), ylim = c(300, 0),
       xlab = expression("Carbon density  (kg C m"^-2*" per cm of depth)"),
       ylab = "Depth below surface (cm)",
       main = paste0("Full-depth carbon profiles: Fort Severn against ",
                     format(nrow(prof), big.mark = ","), " external profiles",
                     "\nbands = 10th-90th percentile, line = median"))
  grid(col = "grey92", lty = 1)

  if (SHOW_RAW_LINES) {
    for (pid in unique(lay_typed$profile_id)) {
      s <- lay_typed[lay_typed$profile_id == pid, ]
      st <- s$soil_type[1]
      if (is.na(st)) next
      lines(s$dens_kgm2_per_cm, (s$upper_depth + s$lower_depth) / 2,
            col = fade(COL[[st]], 0.05), lwd = 0.5)
    }
  } else {
    for (nm in names(env_type)) draw_env(env_type[[nm]], COL[[nm]])
  }

  # Our eight cores, individually, on top.
  for (cc in unique(oseg$core_id)) {
    s <- oseg[oseg$core_id == cc, ]
    s <- s[order(s$depth_top_cm), ]
    for (i in seq_len(nrow(s))) {
      lines(rep(s$carbon_density_kgm2_per_cm[i], 2),
            c(s$depth_top_cm[i], s$depth_bottom_cm[i]),
            col = COL[["ours"]], lwd = 2.2)
    }
  }
  abline(h = CFG$reference_depth_cm, col = "grey40", lty = 3)
  text(xmax * 0.98, CFG$reference_depth_cm - 6,
       paste0(CFG$reference_depth_cm, " cm - the full depth of our cores"),
       adj = 1, cex = 0.72, col = "grey30")

  legend("bottomright",
         legend = c("organic (external)", "mineral (external)",
                    "Fort Severn cores (n=8)"),
         col = c(COL[["organic"]], COL[["mineral"]], COL[["ours"]]),
         lwd = c(2.6, 2.6, 2.2), bty = "n", cex = 0.78)
})

# -----------------------------------------------------------------------------
# FIGURE 2: by ecosystem class (needs step 09)
# -----------------------------------------------------------------------------

if (have_eco && "gwl_label" %in% names(prof)) {
  lay_eco <- lay %>%
    left_join(prof %>% select(dataset, profile_id, gwl_label), by = c("dataset", "profile_id"))
  keep <- names(sort(table(prof$gwl_label), decreasing = TRUE))
  keep <- head(keep[!is.na(keep)], 5)
  lay_eco <- lay_eco[lay_eco$gwl_label %in% keep, ]
  env_eco <- envelope(lay_eco, "gwl_label")
  pal_eco <- setNames(hcl.colors(length(keep), "Dark 3"), keep)

  our_eco <- if ("gwl_label" %in% names(ours_std)) unique(na.omit(ours_std$gwl_label)) else character(0)

  fig("compare_profiles_ecosystem", {
    par(mar = c(4.6, 4.6, 3.8, 1))
    plot(NA, xlim = c(0, xmax), ylim = c(300, 0),
         xlab = expression("Carbon density  (kg C m"^-2*" per cm of depth)"),
         ylab = "Depth below surface (cm)",
         main = paste0("Carbon profiles by wetland class (GWL_FCS30)",
                       "\nbands = 10th-90th percentile, line = median"))
    grid(col = "grey92", lty = 1)
    for (nm in names(env_eco)) draw_env(env_eco[[nm]], pal_eco[[nm]])
    for (cc in unique(oseg$core_id)) {
      s <- oseg[oseg$core_id == cc, ]; s <- s[order(s$depth_top_cm), ]
      for (i in seq_len(nrow(s))) {
        lines(rep(s$carbon_density_kgm2_per_cm[i], 2),
              c(s$depth_top_cm[i], s$depth_bottom_cm[i]),
              col = COL[["ours"]], lwd = 2.2)
      }
    }
    legend("bottomright",
           legend = c(keep, paste0("Fort Severn cores",
                                   if (length(our_eco)) paste0(" (", paste(our_eco, collapse = ", "), ")") else "")),
           col = c(pal_eco[keep], COL[["ours"]]), lwd = 2.4, bty = "n", cex = 0.74)
  })
} else {
  msg("skipped compare_profiles_ecosystem.png (needs step 09)")
}

# -----------------------------------------------------------------------------
# FIGURE 3: the data gap, as distance rings
# -----------------------------------------------------------------------------

all_prof <- bind_rows(prof, ours_std)
rings <- c(0, E$rings_km)
ring_lab <- paste0(head(rings, -1), "-", rings[-1])

usable <- all_prof %>% filter(reaches_30cm, stock_kgm2_0_30 > 0)
tab <- t(vapply(seq_along(ring_lab), function(i) {
  sel <- usable$dist_fort_severn_km > rings[i] & usable$dist_fort_severn_km <= rings[i + 1]
  c(mineral = sum(sel & usable$soil_type == "mineral"),
    organic = sum(sel & usable$soil_type == "organic"))
}, numeric(2)))
rownames(tab) <- ring_lab

fig("compare_datagap_rings", {
  par(mar = c(4.8, 4.8, 4.2, 1))
  bp <- barplot(t(tab), beside = TRUE, las = 1,
                col = c(COL[["mineral"]], COL[["organic"]]),
                xlab = "Distance from Fort Severn (km)",
                ylab = "Profiles with a complete 0-30 cm carbon stock",
                main = paste0("Where comparable measurements exist, and where they do not",
                              "\nour cores at 0 km, then no external profile of either kind until 700 km"))
  legend("topleft", legend = c("mineral soil", "organic soil"),
         fill = c(COL[["mineral"]], COL[["organic"]]), bty = "n", cex = 0.82)
  # Name the empty rings explicitly -- an absent bar is easy to miss.
  empty <- which(rowSums(tab) == 0)
  if (length(empty)) {
    text(colMeans(bp)[empty], max(tab) * 0.04, "none", srt = 90,
         adj = 0, cex = 0.7, col = "grey25")
  }
})

# -----------------------------------------------------------------------------
# FIGURE 4: 0-30 cm stock by soil type, our cores overlaid
# -----------------------------------------------------------------------------

fig("compare_stock_boxplots", {
  par(mar = c(4.6, 4.8, 4.0, 1))
  d <- usable %>% filter(dataset != "Fort Severn", soil_type %in% c("organic", "mineral"))
  boxplot(stock_kgm2_0_30 ~ soil_type, data = d, outline = FALSE,
          col = c(fade(COL[["mineral"]], 0.25), fade(COL[["organic"]], 0.25)),
          border = c(COL[["mineral"]], COL[["organic"]]),
          ylab = expression("Carbon stock, 0-30 cm  (kg C m"^-2*")"),
          xlab = "", las = 1,
          main = paste0("0-30 cm carbon stock: our cores against every comparable profile",
                        "\nlower-bound profiles excluded from both"))
  grid(col = "grey92", lty = 1)
  # Our eight, jittered onto the matching category.
  set.seed(CFG$seed)
  om <- ours_std %>% filter(soil_type %in% c("organic", "mineral"))
  xpos <- ifelse(om$soil_type == "mineral", 1, 2) + runif(nrow(om), -0.13, 0.13)
  pch <- ifelse(om$stock_0_30_is_lower_bound, 2, 19)
  points(xpos, om$stock_kgm2_0_30, pch = pch, col = COL[["ours"]], cex = 1.25, lwd = 1.6)
  legend("topright",
         legend = c("Fort Severn core", "Fort Severn core (lower bound)"),
         pch = c(19, 2), col = COL[["ours"]], bty = "n", cex = 0.76)
})

# -----------------------------------------------------------------------------
# FIGURE 5: the map
# -----------------------------------------------------------------------------

hbl <- st_union(st_geometry(st_read(E$file_hbl, quiet = TRUE)))

fig("compare_datagap_map", {
  par(mar = c(4.2, 4.2, 3.8, 1))
  bb <- st_bbox(hbl)
  plot(NA, xlim = c(bb["xmin"] - 4, bb["xmax"] + 3), ylim = c(bb["ymin"] - 1, bb["ymax"] + 1.5),
       xlab = "Longitude", ylab = "Latitude", asp = 1/cos(56 * pi / 180),
       main = paste0("Every comparable soil profile near Fort Severn",
                     "\nHudson & James Bay Lowlands outlined"))
  plot(hbl, add = TRUE, col = fade("#2E8B6F", 0.07), border = "#2E8B6F", lwd = 1.6)
  for (ds in c("NPDB", "CanPeat", "WOSIS Canada", "Janousek")) {
    d <- prof[prof$dataset == ds, ]
    if (!nrow(d)) next
    points(d$longitude, d$latitude, pch = 21, cex = 0.5,
           bg = fade(DS_COL[[ds]], 0.5), col = NA)
  }
  points(ours$longitude, ours$latitude, pch = 21, bg = COL[["ours"]],
         col = "white", cex = 1.5, lwd = 1.4)
  legend("bottomleft",
         legend = c(paste0(c("NPDB", "CanPeat", "WOSIS Canada", "Janousek"),
                           " (", vapply(c("NPDB","CanPeat","WOSIS Canada","Janousek"),
                                        function(d) sum(prof$dataset == d), integer(1)), ")"),
                    paste0("Fort Severn cores (", nrow(ours), ")")),
         pch = 21, pt.bg = c(DS_COL[["NPDB"]], DS_COL[["CanPeat"]],
                             DS_COL[["WOSIS Canada"]], DS_COL[["Janousek"]], COL[["ours"]]),
         col = NA, bty = "n", cex = 0.74)
})

# -----------------------------------------------------------------------------
# 3. GeoPackage -- one layer per dataset
# -----------------------------------------------------------------------------
# Aggregated to one point per profile, for symbolising by carbon stock.
#
# BOTH stock columns are carried on purpose. Symbolising size by
# stock_kgm2_total makes our 30 cm cores (3-14) almost invisible beside
# CanPeat's full columns (median 119) -- not because they hold less carbon but
# because they are a tenth of the depth. Use stock_kgm2_0_30 for a like-for-like
# size comparison and stock_kgm2_total to show the depth story.

gpkg <- file.path(CFG$dir_current, "comparison_cores.gpkg")
if (file.exists(gpkg)) unlink(gpkg)

layer_of <- list(
  fort_severn_cores = ours_std,
  canpeat  = prof %>% filter(dataset == "CanPeat"),
  npdb     = prof %>% filter(dataset == "NPDB"),
  wosis    = prof %>% filter(dataset == "WOSIS Canada"),
  janousek = prof %>% filter(dataset == "Janousek")
)

for (nm in names(layer_of)) {
  d <- layer_of[[nm]]
  if (!nrow(d)) { msg("  (", nm, ": no rows, layer skipped)"); next }
  st_write(st_as_sf(d, coords = c("longitude", "latitude"),
                    crs = CFG$crs_geographic, remove = FALSE),
           gpkg, layer = nm, append = file.exists(gpkg), quiet = TRUE)
  msg("  gpkg layer '", nm, "': ", nrow(d), " profiles")
}
msg("wrote comparison_cores.gpkg (", length(layer_of), " layers)")

# -----------------------------------------------------------------------------
# 4. The narrative, in numbers
# -----------------------------------------------------------------------------

in_hbl <- prof %>% filter(in_hbl)
hbl_min <- sum(in_hbl$soil_type == "mineral")

msg("")
msg("================ THE FINDING ================")
msg("Inside the Hudson & James Bay Lowlands, across CanPeat, NPDB, WOSIS")
msg("Canada and Janousek combined: ", nrow(in_hbl), " external profiles, of which")
msg("MINERAL: ", hbl_min, ".")
msg("")
msg("Our ", sum(ours$campaign == "mineral"), " mineral cores are therefore the first ",
   if (hbl_min == 0) "mineral soil carbon" else "additional mineral")
msg("measurements in the lowlands in any of these open databases.")
msg("")
msg("Nearest external profile with a complete 0-30 cm stock: ",
   round(min(usable$dist_fort_severn_km[usable$dataset != "Fort Severn"], na.rm = TRUE), 0), " km.")
msg("=============================================")
msg("")
msg("10 complete")
