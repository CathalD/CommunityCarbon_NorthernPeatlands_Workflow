# =============================================================================
# 13_regional_comparison.R      ADD-ON, needs Earth Engine.
#
# Place the Fort Severn cores against four reference groups, which is the
# comparison a presentation actually needs:
#
#   A  MOST SIMILAR CORES     CanPeat profiles within CFG$ext$coastal_buffer_km
#                             of permanent water inside the Lowlands -- the
#                             coastal Hudson/James Bay cores most like these.
#   B  THE LOWLANDS           mean stock across the Hudson & James Bay Lowlands,
#                             at 0-30 cm AND full column, from both the measured
#                             profiles and the published rasters.
#   C  THE REST OF ONTARIO    Sothe et al. 0-30 cm and SoilGrids 0-30 cm averaged
#                             over the province.
#   D  CANADA                 published figures, read from
#                             data/literature_values.csv. Blank rows are skipped
#                             rather than guessed at.
#
# WHY DISTANCE-TO-COAST IS COMPUTED HERE RATHER THAN TAKEN FROM STEP 02
#   Step 02's dist_coast_m band only covers the mapping AOI. Group A needs the
#   same quantity at CanPeat profiles spread across the whole Lowlands, so the
#   distance transform is rebuilt over the Lowlands and sampled at those points.
#
# INPUT   mvp/outputs/current/external_profiles.csv   (08)
#         mvp/outputs/current/cores_clean.csv         (01)
#         mvp/data/hbl_boundary.geojson, literature_values.csv
# OUTPUT  mvp/outputs/current/regional_comparison.csv
#         mvp/outputs/current/canpeat_coastal.csv
#         mvp/outputs/current/figures/context_9_four_groups.png
#
# REQUIRES: readr, dplyr, sf, rgee
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
library(sf)
library(rgee)

msg("13  regional comparison: four reference groups")

E <- CFG$ext
REF <- CFG$reference_depth_cm
ensure_dir(CFG$dir_figures)

prof <- read_csv(file.path(CFG$dir_current, "external_profiles.csv"), show_col_types = FALSE)
ours <- read_csv(file.path(CFG$dir_current, "cores_clean.csv"), show_col_types = FALSE) %>%
  mutate(soil_type = ifelse(campaign == "peat", "organic", "mineral"))

ee_Initialize(project = CFG$gee$project)
msg("Earth Engine initialised, project: ", CFG$gee$project)

hbl_sf <- st_read(E$file_hbl, quiet = TRUE) %>% st_transform(CFG$crs_geographic)
hbl_ee <- sf_as_ee(st_sf(geometry = st_union(st_geometry(hbl_sf))))$geometry()

# -----------------------------------------------------------------------------
# 1. GROUP A -- the most similar cores: coastal CanPeat inside the Lowlands
# -----------------------------------------------------------------------------

gsw <- ee$Image(CFG$gee$asset_jrc_water)
sea <- gsw$select("occurrence")$gt(90)$selfMask()
dist_coast <- sea$fastDistanceTransform(4096)$sqrt()$
  multiply(ee$Image$pixelArea()$sqrt())$
  rename("dist_coast_m")$unmask(-1)

cp <- prof %>% filter(dataset == "CanPeat", in_hbl,
                     !is.na(latitude), !is.na(longitude))
msg("CanPeat profiles inside the Lowlands: ", nrow(cp))

if (nrow(cp)) {
  fc <- ee$FeatureCollection(lapply(seq_len(nrow(cp)), function(i) {
    ee$Feature(ee$Geometry$Point(c(cp$longitude[i], cp$latitude[i])),
               list(pid = as.character(cp$profile_id[i])))
  }))
  samp <- dist_coast$sampleRegions(collection = fc, properties = list("pid"),
                                  scale = 300, geometries = FALSE)
  feats <- samp$getInfo()$features
  dc <- do.call(rbind, lapply(feats, function(f)
    data.frame(profile_id = as.character(f$properties$pid),
              dist_coast_m = as.numeric(f$properties$dist_coast_m),
              stringsAsFactors = FALSE)))
  cp <- left_join(cp, dc, by = "profile_id") %>%
    mutate(dist_coast_km = ifelse(is.finite(dist_coast_m) & dist_coast_m >= 0,
                                 dist_coast_m / 1000, NA_real_),
          is_coastal = !is.na(dist_coast_km) & dist_coast_km <= E$coastal_buffer_km)
  write_csv(cp, file.path(CFG$dir_current, "canpeat_coastal.csv"))
  msg("  within ", E$coastal_buffer_km, " km of the coast: ", sum(cp$is_coastal),
     " of ", nrow(cp))
} else {
  cp$dist_coast_km <- numeric(0); cp$is_coastal <- logical(0)
}

group_a <- cp %>% filter(is_coastal)

# -----------------------------------------------------------------------------
# 2. Published rasters averaged over each region
# -----------------------------------------------------------------------------

#' SoilGrids 0-30 cm stock, kg C/m2, built from concentration x bulk density so
#' the integration matches the one applied to the cores.
soilgrids_0_30 <- function() {
  soc  <- ee$Image(E$gee_asset_soilgrids_soc)
  bdod <- ee$Image(E$gee_asset_soilgrids_bdod)
  parts <- lapply(seq_along(E$soilgrids_bands), function(i) {
    b <- E$soilgrids_bands[i]; th <- E$soilgrids_thick[i]
    # SoilGrids band names carry the variable and statistic, e.g.
    # "soc_15-30cm_mean" -- selecting the bare depth string "15-30cm" matches
    # nothing and fails with "Band pattern did not match any bands".
    soc$select(paste0("soc_", b, "_mean"))$divide(10)$
      multiply(bdod$select(paste0("bdod_", b, "_mean"))$divide(100))$
      multiply(th)$divide(100)
  })
  Reduce(function(a, b) a$add(b), parts)$rename("soilgrids_0_30_kgm2")
}

# Canada's GAUL outline has ~6.7 million edges, far past Earth Engine's 2 million
# limit, so requesting it raw fails with "Geometry has too many edges". Simplify
# before use: for a national mean at 1 km, a 5 km boundary tolerance is
# immaterial, and it is a better answer than the bounding box, which would drag
# in ocean and a slice of Alaska.
simplify_geom <- function(fc, max_error_m = 5000) {
  fc$geometry()$simplify(maxError = max_error_m)
}

regions <- list(
  "Hudson & James Bay Lowlands" = hbl_ee,
  "Ontario" = simplify_geom(ee$FeatureCollection(E$gee_asset_gaul1)$
    filter(ee$Filter$eq("ADM1_NAME", "Ontario"))),
  "Canada"  = simplify_geom(ee$FeatureCollection(E$gee_asset_gaul0)$
    filter(ee$Filter$eq("ADM0_NAME", "Canada")))
)

layers <- list(
  "SoilGrids 0-30cm" = list(img = soilgrids_0_30(), basis = "0-30 cm"),
  "Sothe 0-30cm"     = list(img = ee$Image(CFG$gee$asset_sothe_sc_0_30)$
                                    rename("sothe_0_30_kgm2"), basis = "0-30 cm"),
  "Li et al. full column" = list(img = ee$Image(CFG$gee$asset_li2025)$
                                    rename("li_full_kgm2"), basis = "full peat column")
)

#' Regional mean and sd of one image. Returns NA rather than failing, because a
#' private or region-limited asset legitimately has nothing to say about Canada.
region_stat <- function(img, geom, label, region_name) {
  out <- tryCatch({
    r <- img$reduceRegion(
      reducer = ee$Reducer$mean()$combine(ee$Reducer$stdDev(), sharedInputs = TRUE)$
        combine(ee$Reducer$count(), sharedInputs = TRUE),
      geometry = geom, scale = E$regional_scale_m,
      maxPixels = 1e12, bestEffort = TRUE)$getInfo()
    v <- unlist(r)
    c(mean = unname(v[grepl("_mean$", names(v))][1]),
      sd   = unname(v[grepl("_stdDev$", names(v))][1]),
      n    = unname(v[grepl("_count$", names(v))][1]))
  }, error = function(e) {
    msg("  ", label, " over ", region_name, ": ", conditionMessage(e))
    c(mean = NA_real_, sd = NA_real_, n = NA_real_)
  })
  out
}

raster_rows <- list()
for (rn in names(regions)) {
  for (ln in names(layers)) {
    s <- region_stat(layers[[ln]]$img, regions[[rn]], ln, rn)
    if (is.finite(s[["mean"]])) {
      msg(sprintf("  %-24s over %-28s mean %8.2f  sd %7.2f  (n=%s)",
                 ln, rn, s[["mean"]], s[["sd"]], format(s[["n"]], big.mark = ",")))
    }
    raster_rows[[length(raster_rows) + 1L]] <- data.frame(
      group = rn, source = ln, depth_basis = layers[[ln]]$basis,
      evidence = "published raster",
      mean_kgm2 = s[["mean"]], sd_kgm2 = s[["sd"]],
      n = s[["n"]], stringsAsFactors = FALSE)
  }
}
raster_tbl <- bind_rows(raster_rows) %>% filter(is.finite(mean_kgm2))

# -----------------------------------------------------------------------------
# 3. Measured-profile summaries
# -----------------------------------------------------------------------------

#' Summarise one group's stocks.
#'
#' NOTE ON LOWER BOUNDS. This keeps every profile with a positive stock,
#' including cores that stopped short of the reference depth, because for a
#' REGIONAL MEAN dropping them would discard two of our three peat cores and
#' leave the peat figure resting on one observation. Figure 2 in step 11 takes
#' the opposite decision for the same reason in reverse: a BOX PLOT of the
#' distribution should not mix floors with measurements. So the two are
#' deliberately different, and the n column is what reconciles them -- step 13
#' reports Fort Severn peat as n=3, step 11 as n=1.
summarise_profiles <- function(d, group, basis, col) {
  v <- d[[col]]
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(NULL)
  data.frame(group = group, source = "measured profiles", depth_basis = basis,
            evidence = "cores", mean_kgm2 = mean(v), sd_kgm2 = sd(v),
            n = length(v), stringsAsFactors = FALSE)
}

hbl_prof <- prof %>% filter(in_hbl)
measured <- bind_rows(
  summarise_profiles(ours, "Fort Severn cores (all)", paste0("0-", REF, " cm"), "stock_kgm2"),
  summarise_profiles(ours %>% filter(soil_type == "organic"),
                    "Fort Severn peat cores", paste0("0-", REF, " cm"), "stock_kgm2"),
  summarise_profiles(ours %>% filter(soil_type == "mineral"),
                    "Fort Severn mineral cores", paste0("0-", REF, " cm"), "stock_kgm2"),
  summarise_profiles(group_a %>% filter(reaches_30cm),
                    sprintf("Coastal CanPeat, within %d km of coast", E$coastal_buffer_km),
                    paste0("0-", REF, " cm"), "stock_kgm2_0_30"),
  summarise_profiles(group_a,
                    sprintf("Coastal CanPeat, within %d km of coast", E$coastal_buffer_km),
                    "full column", "stock_kgm2_total"),
  summarise_profiles(hbl_prof %>% filter(reaches_30cm),
                    "Hudson & James Bay Lowlands", paste0("0-", REF, " cm"), "stock_kgm2_0_30"),
  summarise_profiles(hbl_prof, "Hudson & James Bay Lowlands", "full column", "stock_kgm2_total")
)

# -----------------------------------------------------------------------------
# 4. GROUP D -- published Canada figures, read not computed
# -----------------------------------------------------------------------------

lit <- if (file.exists(E$file_literature)) {
  read_csv(E$file_literature, show_col_types = FALSE) %>%
    filter(is.finite(mean_kgm2)) %>%
    transmute(group, source = ifelse(is.na(source), "literature", source),
             depth_basis, evidence = "literature",
             mean_kgm2, sd_kgm2, n = NA_real_)
} else tibble()

if (!nrow(lit)) {
  msg("NOTE: no usable rows in ", basename(E$file_literature),
     " -- fill in the Canada figures there and re-run to add group D. ",
     "They are deliberately not guessed at here.")
} else {
  msg("literature rows used: ", nrow(lit))
}

# -----------------------------------------------------------------------------
# 5. The table
# -----------------------------------------------------------------------------

comparison <- bind_rows(measured, raster_tbl, lit) %>%
  mutate(across(c(mean_kgm2, sd_kgm2), ~ round(.x, 2))) %>%
  arrange(depth_basis, desc(mean_kgm2))

# A published raster that disagrees with the measured profiles by several fold
# over the SAME region and depth is a units or depth-support question, not a
# finding. Surface it here rather than let it into a slide unremarked.
shallow <- comparison %>% filter(depth_basis == paste0("0-", REF, " cm"))
meas_hbl <- shallow$mean_kgm2[shallow$group == "Hudson & James Bay Lowlands" &
                             shallow$evidence == "cores"]
rast_hbl <- shallow$mean_kgm2[shallow$group == "Hudson & James Bay Lowlands" &
                             shallow$evidence == "published raster"]
if (length(meas_hbl) && length(rast_hbl)) {
  ratio <- max(rast_hbl, na.rm = TRUE) / meas_hbl[1]
  if (is.finite(ratio) && (ratio > 2 || ratio < 0.5)) {
    msg("")
    msg("CHECK BEFORE PRESENTING: over the same region and the same 0-", REF,
       " cm window, the published raster mean (", round(max(rast_hbl, na.rm = TRUE), 1),
       ") is ", round(ratio, 1), "x the measured-profile mean (",
       round(meas_hbl[1], 1), ").")
    msg("  A gap that size is usually units or depth support, not ground truth. ",
       "Confirm what the raster's values actually represent before putting both ",
       "numbers on one slide.")
  }
}

write_csv(comparison, file.path(CFG$dir_current, "regional_comparison.csv"))
msg("--- regional comparison ---")
print(as.data.frame(comparison), right = FALSE)
msg("wrote regional_comparison.csv")

# -----------------------------------------------------------------------------
# 6. The figure -- two panels, because the depth bases must not be mixed
# -----------------------------------------------------------------------------

fig("context_9_four_groups", {
  op <- par(mfrow = c(1, 2), mar = c(4.6, 15.5, 3.6, 1), oma = c(0, 0, 3.0, 0))
  on.exit(par(op), add = TRUE)
  # Two panels, and the split must be exact. depth_basis carries three distinct
  # values -- "0-30 cm", "full column" (measured profiles) and "full peat
  # column" (Li et al. and the literature rows) -- and the two full-depth labels
  # belong together on the right while the shallow one stands alone on the left.
  panels <- list(
    `shallow` = list(title = paste0("0-", REF, " cm"),
                    keep = paste0("0-", REF, " cm")),
    `deep`    = list(title = "full column / full peat column",
                    keep = c("full column", "full peat column"))
  )
  for (pn in names(panels)) {
    basis <- panels[[pn]]$title
    d <- comparison %>% filter(depth_basis %in% panels[[pn]]$keep)
    if (!nrow(d)) { plot.new(); title(main = paste0("no data for ", basis)); next }
    d <- d %>% arrange(mean_kgm2)
    bar_col <- ifelse(grepl("^Fort Severn", d$group), COL[["ours"]],
                     ifelse(d$evidence == "literature", "#8A8A8A",
                            ifelse(d$evidence == "cores", COL[["organic"]], COL[["mineral"]])))
    labs <- ifelse(d$source == "measured profiles",
                  sprintf("%s  (n=%s)", d$group, format(d$n, big.mark = ",")),
                  sprintf("%s - %s", d$group, d$source))
    bp <- barplot(d$mean_kgm2, horiz = TRUE, names.arg = labs, las = 1,
                  col = bar_col, border = NA, cex.names = 0.62,
                  xlim = c(0, max(d$mean_kgm2 + ifelse(is.finite(d$sd_kgm2), d$sd_kgm2, 0),
                                  na.rm = TRUE) * 1.08),
                  xlab = expression("kg C m"^-2), main = basis)
    ok <- is.finite(d$sd_kgm2)
    if (any(ok)) arrows(d$mean_kgm2[ok], bp[ok],
                        pmax(0, d$mean_kgm2[ok] + d$sd_kgm2[ok]), bp[ok],
                        angle = 90, length = 0.03, col = "grey30", lwd = 1)
    grid(col = "grey93", lty = 1, ny = NA)
  }
  mtext(paste0("Fort Severn in context: coastal CanPeat cores, the Lowlands, Ontario and Canada",
               "\nbars = mean, whisker = +1 SD.  The two panels are DIFFERENT quantities and must not be compared across."),
        outer = TRUE, cex = 0.88, font = 2)
}, width = 2200, height = 1100)

msg("13 complete")
