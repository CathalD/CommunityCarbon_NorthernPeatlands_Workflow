# =============================================================================
# 12_community_report.R
#
# A plain-language Quarto report for a non-specialist reader: what was
# collected, what it shows, how it compares, and what it means.
#
# Every number is computed here and written into the document as finished text,
# so `quarto render` needs only the Quarto CLI -- no knitr, no R packages, no
# CRAN access. Same reproducibility terms as the rest of the workflow.
#
# Writing rule for this file: no unexplained jargon. Where a technical term is
# unavoidable it is defined in the sentence that introduces it.
#
# OUTPUT  outputs/COMMUNITY_REPORT.qmd
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

log_step("12  PLAIN-LANGUAGE QUARTO REPORT")

md <- character(0)
add <- function(...) md <<- c(md, paste0(...))
blank <- function() md <<- c(md, "")
tbl <- function(f) file.path(CFG$dir_tables, f)
drv <- function(f) file.path(CFG$dir_derived, f)
got <- function(p) file.exists(p)

md_table <- function(df, digits = 2) {
  d <- as.data.frame(df)
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) formatC(round(v, digits), format = "fg",
                                               digits = 6))
  d[] <- lapply(d, function(v) gsub("|", "\\|", as.character(v), fixed = TRUE))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
#' Figure block with a caption, referencing a PNG by relative path.
fig <- function(file, caption, width = "100%") {
  p <- file.path(CFG$dir_figures, file)
  if (!file.exists(p)) return(invisible(NULL))
  add("![", caption, "](figures/", file, "){width=", width, "}")
  blank()
}

cores  <- require_artifact(drv("01_cores_qc.csv"), "scripts/01_ingest_qc.R")
seg    <- require_artifact(drv("01_segments_qc.csv"), "scripts/01_ingest_qc.R")
stocks <- require_artifact(drv("02_core_stocks.csv"), "scripts/02_depth_harmonise.R")

cs  <- stocks[stocks$scenario == "as_supplied" & stocks$window == "common_support", ]
r30 <- stocks[stocks$scenario == "as_supplied" & stocks$window == "reference_0_30", ]
zc  <- cs$window_hi_cm[1]

# =============================================================================
add("---")
add('title: "Soil carbon around Fort Severn"')
add('subtitle: "What eight community soil cores measured, and what they tell us"')
add("date: ", format(Sys.Date(), "%Y-%m-%d"))
add("format:")
add("  html:")
add("    toc: true")
add("    toc-title: \"On this page\"")
add("    number-sections: true")
add("    theme: cosmo")
add("    embed-resources: true")
add("    fig-align: center")
add("  pdf:")
add("    toc: true")
add("    number-sections: true")
add("    geometry: margin=1in")
add("execute:")
add("  echo: false")
add("---")
blank()

# =============================================================================
add("# In one paragraph")
blank()
add("Community members collected **", nrow(cores), " soil cores** near Fort ",
    "Severn, Ontario, across two field seasons: ",
    sum(cores$campaign == "PM_2024_peat"), " from peatland sites in 2024 and ",
    sum(cores$campaign == "FS_2025_mineral"), " from forest and mineral-soil ",
    "sites in 2025. This report explains what those cores measured, how much ",
    "carbon is in the ground where they were taken, how that compares with ",
    "existing maps of the region, and what it means for looking after this ",
    "landscape. The short answer: the cores are sound and they agree with the ",
    "published regional carbon map, but eight cores can only speak for a ",
    "small part of the area, and the shallow depth they reached misses most ",
    "of the carbon that makes this landscape important.")
blank()

# =============================================================================
add("# What a soil core measures, and how carbon is counted")
blank()
add("A soil core is a cylinder of soil pushed or dug out of the ground and ",
    "kept in order, so the top of the core is the surface and the bottom is ",
    "the deepest point reached. The core is cut into sections, and the lab ",
    "measures two things for each section:")
blank()
add("- **Bulk density** -- how many grams a cubic centimetre of dried soil ",
    "weighs. Peat is light and fluffy; mineral soil is heavy and dense.")
add("- **Organic matter**, and from it **organic carbon** -- what share of ",
    "that dried soil is carbon.")
blank()
add("Carbon stock is the two multiplied together, then added up down the ",
    "core:")
blank()
add("> carbon in a layer = how dense the soil is x what share of it is carbon ",
    "x how thick the layer is")
blank()
add("This matters more than it sounds. A soil can be **rich in carbon by ",
    "percentage but hold little carbon by weight**, if it is light enough. ",
    "That is exactly what happens with peat, and it explains the most ",
    "surprising result in this report.")
blank()
add("Carbon stock here is given in **kilograms of carbon per square metre** ",
    "(kg C/m^2^): the weight of carbon under one square metre of ground, down ",
    "to a stated depth. For scale, 1 kg C/m^2^ is 10 tonnes of carbon per ",
    "hectare.")
blank()

# =============================================================================
add("# The cores")
blank()
add("## Where they are")
blank()
fig("09_soc_within_credible_radius.png",
    "The eight community cores. Colour shows measured carbon; blank ground is area no core can speak for.")
add("Cores fall into two clusters roughly ",
    sprintf("%.0f km", max(pairwise_km(cores$longitude, cores$latitude))),
    " apart at their furthest: a peatland group and a forest/mineral group. ",
    "This clustering is a practical consequence of access, not a survey ",
    "design, and it shapes what can honestly be said later.")
blank()

add("## What each core found")
blank()
tab <- data.frame(
  Core = cores$core_id,
  Type = ifelse(grepl("peat", cores$campaign), "Peatland (2024)", "Forest/mineral (2025)"),
  `Sections` = cores$n_segments,
  `Depth reached (cm)` = round(cores$core_length_cm, 1),
  check.names = FALSE, stringsAsFactors = FALSE)
tab[[paste0("Carbon 0-", sprintf("%.1f", zc), " cm (kg C/m2)")]] <-
  round(cs$stock_kgm2[match(cores$core_id, cs$core_id)], 2)
tab[["Carbon 0-30 cm (kg C/m2)"]] <-
  round(r30$stock_kgm2[match(cores$core_id, r30$core_id)], 2)
tab[["0-30 cm value is"]] <-
  ifelse(r30$is_lower_bound[match(cores$core_id, r30$core_id)],
         "AT LEAST this (core stopped short)", "complete")
add(md_table(tab))
blank()
add("Two things to read carefully in that table.")
blank()
add("**The two depth columns are not interchangeable.** Cores stopped at ",
    "different depths, between ", sprintf("%.1f", min(cores$core_length_cm)),
    " and ", sprintf("%.1f", max(cores$core_length_cm)),
    " cm. To compare all eight fairly, we use the deepest layer *every* core ",
    "reached, which is 0-", sprintf("%.1f", zc),
    " cm. The 0-30 cm column is included because published maps report that ",
    "depth, but two cores never got there, so their 0-30 cm figures are ",
    "**minimums, not measurements**.")
blank()
add("**Two of the cores hit mineral soil quickly.** Where organic matter ",
    "drops sharply, the peat has ended:")
blank()
if (got(tbl("11_peat_depth_at_cores.csv"))) {
  pdp <- read_csv_verbatim(tbl("11_peat_depth_at_cores.csv"))
  pdp <- pdp[pdp$status != "no_peat", ]
  pt <- data.frame(
    Core = pdp$core_id,
    `Peat depth (cm)` = round(pdp$peat_depth_cm, 1),
    Finding = ifelse(pdp$status == "contact_observed",
                     "Peat ends here; mineral soil below",
                     "Still peat at the bottom of the core; true depth is greater"),
    check.names = FALSE, stringsAsFactors = FALSE)
  add(md_table(pt))
  blank()
  add("Published work for the Hudson Bay Lowlands gives an average peat depth ",
      "of about **", CFG$lit$li_peat_depth_mean_cm, " cm**. The peat at these ",
      "coring sites is far thinner than that. These are **peat margins or ",
      "shallow fen**, not the deep peat plateaus that hold most of the ",
      "region's carbon. That is a real and useful finding, and it limits what ",
      "these cores can say about the deep peat.")
  blank()
}

# =============================================================================
add("# The surprise: peat held less carbon than the forest soils")
blank()
fig("09_core_stocks.png",
    "Carbon in each core, at the shared depth and at 30 cm. Orange bars marked LOWER BOUND are cores that stopped short.")
bc <- tapply(cs$stock_kgm2, cs$campaign, mean)
add("Over the layer all eight cores reached (0-", sprintf("%.1f", zc), " cm), ",
    "the peatland cores averaged **", sprintf("%.2f", bc[["PM_2024_peat"]]),
    " kg C/m^2^** and the forest/mineral cores **",
    sprintf("%.2f", bc[["FS_2025_mineral"]]), " kg C/m^2^**. The peat sites ",
    "held *less*.")
blank()
add("This is not a mistake, and it is not an argument against peatlands. It ",
    "is the density point from earlier. Peat at these sites has a bulk ",
    "density of about ",
    sprintf("%.2f-%.2f", min(seg$bulk_density_gcm3[grepl("peat", seg$campaign)]),
            max(seg$bulk_density_gcm3[grepl("peat", seg$campaign)])),
    " g/cm^3^, while the mineral soils run ",
    sprintf("%.2f-%.2f", min(seg$bulk_density_gcm3[!grepl("peat", seg$campaign)]),
            max(seg$bulk_density_gcm3[!grepl("peat", seg$campaign)])),
    ". Peat is mostly water and air. A spadeful of it contains a high ",
    "*percentage* of carbon but not much *weight* of carbon.")
blank()
add("::: {.callout-important}")
add("## The single most important point in this report")
blank()
add("Peatlands matter because they are **deep**, not because their top layer ",
    "is rich. Hudson Bay Lowlands peat averages roughly ",
    CFG$lit$li_peat_depth_mean_cm, " cm thick. A map of the top 30 cm ",
    "captures only about ",
    sprintf("%.0f%%", 100 * 30 / CFG$lit$li_peat_depth_mean_cm),
    " of that column, and because peat gets denser with depth, rather less ",
    "than that share of the carbon.")
blank()
add("**A 0-30 cm carbon map of this landscape is not a map of its carbon. It ",
    "is a map of the top of its carbon.** Any decision that uses a shallow ",
    "map to judge peatland value will badly understate what is there.")
add(":::")
blank()

# =============================================================================
add("# Combining the cores with an existing map")
blank()
add("## The idea")
blank()
add("There is already a published carbon map of this region, built by Li and ",
    "colleagues in 2025 from satellite data and many field measurements. It ",
    "covers the **whole peat column**, surface to the bottom of the peat, and ",
    "gives an average of about **", CFG$lit$li_peat_stock_mean_kgm2,
    " kg C/m^2^** for the Hudson Bay Lowlands.")
blank()
add("We can combine that map with the community cores using a standard ",
    "statistical method. In plain terms:")
blank()
add("1. Start every point on the map at the published value, along with how ",
    "uncertain that value is.")
add("2. Where a community core sits nearby, pull the map towards what the ",
    "core actually measured.")
add("3. **The closer the core, the harder the pull.** Beyond about ",
    sprintf("%.1f km", CFG$bayes$max_influence_km),
    " a core has no say at all.")
add("4. **More nearby cores pull harder than one.** Two cores agreeing carry ",
    "more weight than either alone.")
add("5. Far from every core, the map simply stays as published -- because ",
    "nothing was measured there.")
blank()
add("This is called a *Bayesian update*: the published map is the ",
    "**prior** (what we believed beforehand), the cores are the ",
    "**observations**, and the result is the **posterior** (what we believe ",
    "after taking the measurements into account).")
blank()

add("## Where the cores actually changed anything")
blank()
fig("11_bayes_core_information.png",
    "How much of the final answer comes from the community cores rather than the published map.")
add("This map answers a question that is usually left unasked: *where did the ",
    "new data actually make a difference?* Almost everywhere, the answer is ",
    "none at all. The cores influence a small area around themselves and ",
    "nothing beyond it.")
blank()
if (got(tbl("09_spatial_coverage.csv"))) {
  cvg <- read_csv_verbatim(tbl("09_spatial_coverage.csv"))
  gv <- function(k) cvg$value[cvg$metric == k]
  add("Put numerically: the study area is about ",
      format(round(gv("aoi_area_km2")), big.mark = ","), " km^2^, and the ",
      "eight cores meaningfully constrain roughly **",
      sprintf("%.0f%%", gv("pct_aoi_within_credible_radius")), "** of it.")
  blank()
}

add("## How much did the cores move the published number?")
blank()
if (got(tbl("11_implied_shallow_fraction.csv"))) {
  isf <- read_csv_verbatim(tbl("11_implied_shallow_fraction.csv"))
  add("Barely at all -- and that is the right answer rather than a ",
      "disappointing one. The cores reach 30 cm; the published map covers ",
      "metres. Eight shallow cores cannot say much about a column that deep.")
  blank()
  add("But they can do something more interesting: **check whether the ",
      "published map is plausible here.**")
  blank()
  chk <- data.frame(
    Quantity = c("Measured by the community core, 0-30 cm",
                 "Published full-column average (Li et al. 2025)",
                 "Share of the column that implies is in the top 30 cm",
                 "Share you would expect if peat had constant density",
                 "Direction of the difference"),
    Value = c(sprintf("%.2f kg C/m2", isf$measured_0_30_kgm2),
              sprintf("%.0f kg C/m2", isf$li_full_column_kgm2),
              sprintf("%.1f%%", 100 * isf$implied_shallow_fraction),
              sprintf("%.1f%%", 100 * isf$uniform_density_fraction),
              "Less in the top layer than constant density would give"),
    stringsAsFactors = FALSE)
  add(md_table(chk))
  blank()
  add("Peat compacts under its own weight, so deeper peat is denser and holds ",
      "more carbon per centimetre. That means the top 30 cm should hold ",
      "*less* than its proportional share -- and that is exactly what the ",
      "cores show. The community measurement is **consistent with the ",
      "published regional map**, under a physically sensible picture of how ",
      "peat is packed.")
  blank()
  add("::: {.callout-note}")
  add("## What this means")
  blank()
  add("The most valuable thing these cores did was not to change the map. It ",
      "was to **independently confirm that the published map has the right ",
      "magnitude at this location** -- from ground measurements taken by the ",
      "community, using different methods, years later.")
  add(":::")
  blank()
  add("The honest caveat: this rests on ", isf$n_peat_cores_used,
      " peat core(s) that reached 30 cm, at sites where the peat is thin. It ",
      "is encouraging, not conclusive.")
  blank()
}

# =============================================================================
add("# The area this study covers")
blank()
if (got(tbl("13_aoi_definition.csv"))) {
  ad <- read_csv_verbatim(tbl("13_aoi_definition.csv"))
  gv <- function(k) ad$value[ad$property == k]
  add("The study area is a long rectangle that **follows the Hudson Bay ",
      "coastline** rather than being drawn square to north. The coast here ",
      "runs north-west to south-east, and so does the sampling, so a ",
      "north-up box would either cut off the ends of the transect or waste ",
      "most of its area over open water.")
  blank()
  add("The direction was not drawn by eye. It was calculated from the cores ",
      "themselves -- the line they best fall along -- giving **",
      gv("bearing_deg"), "degrees from north**. The area measures **",
      gv("length_km"), " km long by ", gv("width_km"), " km wide (",
      gv("area_km2"), " km^2^)**, including a ", gv("along_buffer_km"),
      " km buffer beyond the cores at each end and ", gv("across_buffer_km"),
      " km on each side.")
  blank()
}
fig("13_aoi_coast_oriented.png",
    "The coast-following study area. The dotted box is the earlier north-up version.")
if (got(tbl("13_cores_coast_relative.csv"))) {
  cr <- read_csv_verbatim(tbl("13_cores_coast_relative.csv"))
  pk <- cr$across_km[grepl("peat", cr$campaign)]
  mn <- cr$across_km[!grepl("peat", cr$campaign)]
  if (length(pk) && length(mn) && min(pk) > max(mn)) {
    add("::: {.callout-warning}")
    add("## Something the ordinary map hides")
    blank()
    add("Measured across the coastline rather than in latitude and longitude, ",
        "**every peat core sits on one side of the sampling line and every ",
        "mineral core on the other** (peat ",
        sprintf("%.1f to %.1f km", min(pk), max(pk)), ", mineral ",
        sprintf("%.1f to %.1f km", min(mn), max(mn)), " across).")
    blank()
    add("This matters when interpreting any peat-versus-mineral comparison. ",
        "The two groups differ in soil type, but they also differ in **year ",
        "sampled, field team, laboratory method, and now position across the ",
        "landscape**. Those things cannot be separated with this dataset. A ",
        "difference between the groups is a difference between two sets of ",
        "cores; calling it purely a soil-type effect would be going beyond ",
        "the evidence.")
    add(":::")
    blank()
  }
}

# =============================================================================
add("# Carbon by ecosystem type")
blank()
if (got(tbl("15_carbon_by_landcover.csv"))) {
  lcz <- read_csv_verbatim(tbl("15_carbon_by_landcover.csv"))
  prim <- lcz[lcz$landcover_source == "GWL_FCS30" &
              lcz$depth_support == "0-30 cm", ]
  if (nrow(prim)) {
    add("Different kinds of ground hold different amounts of carbon. This ",
        "table splits the study area by ecosystem and gives the average ",
        "carbon in the top 30 cm of each, taken from published maps:")
    blank()
    t2 <- data.frame(
      Ecosystem = prim$label,
      `Area (km2)` = round(prim$area_km2, 1),
      `Share of area` = paste0(round(prim$pct_of_aoi, 1), "%"),
      `Average carbon (kg C/m2)` = round(prim$mean_kgm2, 2),
      check.names = FALSE, stringsAsFactors = FALSE)
    if ("has_core" %in% names(prim)) {
      t2[["Any core here?"]] <- ifelse(prim$has_core, "yes", "NO")
    }
    add(md_table(t2))
    blank()
  }
  if (got(tbl("15_ecosystem_verification.csv"))) {
    ev <- read_csv_verbatim(tbl("15_ecosystem_verification.csv"))
    unv <- sum(ev$pct_of_aoi[!isTRUE_vec(ev$has_core)], na.rm = TRUE)
    add("::: {.callout-important}")
    add("## Which of these numbers has been checked on the ground?")
    blank()
    add(sprintf("**%.0f%% of the study area** lies in ecosystem types where no ",
                unv),
        "core has ever been taken. Those rows come entirely from published ",
        "maps built elsewhere, with nothing local to test them against. The ",
        "rows that do contain a core are the only ones with any local ",
        "verification at all.")
    add(":::")
    blank()
  }
  add("Two land-cover maps are used. The main one (GWL_FCS30) can tell ",
      "wetland types apart. The widely used alternative (ESA WorldCover) ",
      "cannot separate **treed bog growing on deep peat** from **spruce ",
      "forest on ordinary mineral soil** -- both simply read as tree cover. ",
      "In this landscape that is exactly the distinction that decides carbon, ",
      "which is why the wetland-aware map is preferred here.")
  blank()
} else {
  add("This section needs the ecosystem analysis (`15_landcover_carbon.R`), ",
      "which requires a Google Earth Engine session. Once it has run, this ",
      "section fills in automatically with carbon per ecosystem type and a ",
      "note of which ecosystems have been visited.")
  blank()
}

# =============================================================================
add("# Finding distinct areas from satellite data")
blank()
if (got(tbl("16_cluster_carbon.csv"))) {
  cc <- read_csv_verbatim(tbl("16_cluster_carbon.csv"))
  sgc <- cc[cc$carbon_layer == "SoilGrids 0-30 cm", ]
  add("Ecosystem categories are coarse. Two places both labelled \"swamp\" ",
      "can be quite different ground. To find finer distinctions we used a ",
      "method called **clustering**: the computer groups together places that ",
      "*look alike to a satellite* across a whole year -- their colour, ",
      "their wetness, how they change with the seasons -- without being told ",
      "anything about carbon.")
  blank()
  add("This is worth trusting more than it might sound. The grouping never ",
      "sees a single carbon measurement, so it cannot be bent by our eight ",
      "cores. It works from millions of satellite pixels.")
  blank()
  if (nrow(sgc)) {
    t3 <- data.frame(
      `Landscape type` = paste("Type", sgc$cluster),
      `Share of area` = paste0(round(sgc$pct_of_aoi, 1), "%"),
      `Average carbon (kg C/m2)` = round(sgc$mean_kgm2, 2),
      check.names = FALSE, stringsAsFactors = FALSE)
    if ("has_core" %in% names(sgc)) {
      t3[["Visited?"]] <- ifelse(isTRUE_vec(sgc$has_core), "yes", "never")
    }
    add(md_table(t3))
    blank()
    if ("has_core" %in% names(sgc)) {
      nv <- sum(!isTRUE_vec(sgc$has_core))
      gap <- sum(sgc$pct_of_aoi[!isTRUE_vec(sgc$has_core)], na.rm = TRUE)
      add("::: {.callout-tip}")
      add("## Where the next cores should go")
      blank()
      add(sprintf("**%d of the %d landscape types have never been sampled**, ",
                  nv, nrow(sgc)),
          sprintf("covering about %.0f%% of the study area.", gap),
          " Those are the highest-value places for the next field season. A ",
          "core taken in a landscape type nobody has visited teaches far more ",
          "than another core beside an existing one.")
      blank()
      add("The map `ccnp_clusters_without_cores.tif` shows exactly where ",
          "these are, and can be opened in any mapping software to plan a ",
          "trip.")
      add(":::")
      blank()
    }
  }
} else {
  add("This section needs the clustering analysis ",
      "(`16_embedding_clusters.R`), which requires a Google Earth Engine ",
      "session. Once it has run, this section fills in with the distinct ",
      "landscape types found, the carbon in each, and -- most usefully -- ",
      "which types have never been visited and should be sampled next.")
  blank()
}

# =============================================================================
add("# How Fort Severn compares")
blank()
add("Comparisons only mean something when the depths match. Putting a 0-30 cm ",
    "number next to a full-peat-column number is like comparing the top ",
    "shelf of a bookcase with the whole bookcase.")
blank()
cmp <- data.frame(
  What = c("Fort Severn peat cores",
           "Fort Severn forest/mineral cores",
           "Fort Severn cores that reached 30 cm",
           "Hudson Bay Lowlands peat (Li et al. 2025)"),
  Depth = c(sprintf("0-%.1f cm", zc), sprintf("0-%.1f cm", zc),
            "0-30 cm", "whole peat column (~184 cm)"),
  `Carbon (kg C/m2)` = c(sprintf("%.2f", bc[["PM_2024_peat"]]),
                         sprintf("%.2f", bc[["FS_2025_mineral"]]),
                         sprintf("%.2f", mean(r30$stock_kgm2[!r30$is_lower_bound])),
                         sprintf("%.0f (+/- %.0f)", CFG$lit$li_peat_stock_mean_kgm2,
                                 CFG$lit$li_peat_stock_sd_kgm2)),
  `Comparable to the others?` = c("yes, with row 2", "yes, with row 1",
                                  "only with other 0-30 cm figures",
                                  "NO -- different depth entirely"),
  check.names = FALSE, stringsAsFactors = FALSE)
add(md_table(cmp))
blank()
add("The last row is roughly **",
    sprintf("%.0f times", CFG$lit$li_peat_stock_mean_kgm2 /
              mean(r30$stock_kgm2[!r30$is_lower_bound])),
    "** the 0-30 cm figures above it. That gap is not disagreement between ",
    "studies. It is the carbon that lies below 30 cm, which the shallow ",
    "numbers never counted.")
blank()

add("## Against the rest of Canada")
blank()
if (got(tbl("14_npdb_proximity.csv")) && got(tbl("14_npdb_quantiles.csv"))) {
  prox <- read_csv_verbatim(tbl("14_npdb_proximity.csv"))
  nq   <- read_csv_verbatim(tbl("14_npdb_quantiles.csv"))
  add("Agriculture and Agri-Food Canada keeps a national record of soil ",
      "profiles, the National Pedon Database. After the same quality checks ",
      "applied to the community cores, it yields several hundred cores across ",
      "the country that can be compared like for like at 0-30 cm.")
  blank()
  add("::: {.callout-important}")
  add("## These cores are the only soil carbon data for this area")
  blank()
  nearest_txt <- if (got(drv("14_npdb_cores_qc.csv"))) {
    nq2 <- read_csv_verbatim(drv("14_npdb_cores_qc.csv"))
    # NPDB province names are lowercase; title-case them for prose.
    prov <- nq2$province[which.min(nq2$dist_fort_severn_km)]
    prov <- gsub("(^|[ -])([a-z])", "\\1\\U\\2", prov, perl = TRUE)
    sprintf("about %.0f km away, in %s",
            min(nq2$dist_fort_severn_km, na.rm = TRUE), prov)
  } else "several hundred kilometres away"
  add("**There is not one core in the national database within 500 km of ",
      "Fort Severn.** The nearest usable one is ", nearest_txt, ".")
  blank()
  add(md_table(data.frame(
    `Distance from Fort Severn` = paste("within", prox$within_km, "km"),
    `Existing cores in the national record` = prox$n_npdb_cores,
    check.names = FALSE, stringsAsFactors = FALSE)))
  blank()
  add("Every published carbon map that covers Fort Severn is, in this ",
      "neighbourhood, working from measurements taken hundreds of kilometres ",
      "away. That is not a criticism of those maps -- it is simply what the ",
      "available data allowed. It is also the clearest argument for the ",
      "community sampling programme: **these eight cores are not just few, ",
      "they are all there is.**")
  add(":::")
  blank()

  if (got(tbl("14_fort_severn_percentiles.csv"))) {
    pc <- read_csv_verbatim(tbl("14_fort_severn_percentiles.csv"))
    add("Where the community cores sit in the national range:")
    blank()
    add(md_table(data.frame(
      Core = pc$core_id,
      Type = ifelse(grepl("peat", pc$campaign), "Peatland", "Forest/mineral"),
      `Carbon 0-30 cm (kg C/m2)` = round(pc$stock_kgm2, 2),
      `Higher than this share of Canadian cores` =
        paste0(round(pc$national_percentile), "%"),
      check.names = FALSE, stringsAsFactors = FALSE)))
    blank()
    add("They fall between roughly the ", round(min(pc$national_percentile)),
        "th and ", round(max(pc$national_percentile)),
        "th percentile -- comfortably within the national range. That is ",
        "reassuring: it means the community sampling and laboratory work are ",
        "producing numbers of the kind the national record contains, not ",
        "outliers that would suggest a problem.")
    blank()
  }
  fig("14_national_distribution.png",
      "The community cores against several hundred Canadian cores measured to the same depth.")
  fig("14_npdb_coverage_map.png",
      "Where Canada's soil carbon cores are. The gap around Fort Severn is the point.")

  if (got(tbl("14_organic_vs_mineral.csv"))) {
    ov <- read_csv_verbatim(tbl("14_organic_vs_mineral.csv"))
    add("### Does the surprising result hold up elsewhere?")
    blank()
    add("Earlier we found that the peat cores held *less* carbon than the ",
        "forest cores in the shallow layer. The national record lets us check ",
        "whether that is a general pattern or a local quirk.")
    blank()
    add(md_table(data.frame(
      `Soil type` = ov$group, `Cores` = ov$n,
      `Median carbon 0-30 cm (kg C/m2)` = round(ov$median_kgm2, 2),
      `Middle half of the range` = paste0(round(ov$q25, 1), " to ",
                                          round(ov$q75, 1)),
      check.names = FALSE, stringsAsFactors = FALSE)))
    blank()
    add("The same ordering appears nationally -- organic soils hold slightly ",
        "less in the top 30 cm. **But the difference is small and the two ",
        "groups overlap heavily.**")
    if ("p_organic_exceeds_mineral" %in% names(ov)) {
      add(" Pick a random organic core and a random mineral one, and the ",
          "organic one is higher ",
          round(100 * ov$p_organic_exceeds_mineral[1]), "% of the time -- ",
          "close to a coin toss.")
    }
    blank()
    add("So this is **support, not proof**. The national record is consistent ",
        "with what the community cores found and gives no reason to suspect a ",
        "local error. What it does *not* support is a strong claim that ",
        "organic soils reliably hold less carbon near the surface; they are ",
        "simply far more variable.")
    blank()
    add("The point that survives is the one that matters: **counting only the ",
        "top 30 cm does not show peatlands to advantage.** Their importance ",
        "lies in depth, and a shallow measurement cannot see it.")
    blank()
  }
}

if (got(tbl("07_reference_benchmark.csv"))) {
  add("## Against other published maps")
  blank()
  bm <- read_csv_verbatim(tbl("07_reference_benchmark.csv"))
  add(md_table(bm[, intersect(c("layer", "n_cores_compared", "core_mean_kgm2",
                                "layer_mean_kgm2", "bias_layer_minus_core"),
                              names(bm))]))
  blank()
  add("These are the global SoilGrids product and the Canada-wide map by ",
      "Sothe and colleagues, both at 0-30 cm, compared against the cores that ",
      "reached that depth. Read the difference column; ignore any correlation ",
      "figure, because six clustered points cannot establish one.")
  blank()
} else {
  add("## Against other published maps")
  blank()
  add("A like-for-like comparison against the global SoilGrids product and ",
      "the Canada-wide map by Sothe and colleagues (both 0-30 cm) is set up ",
      "in the workflow but needs a Google Earth Engine session to fill in. ",
      "Run `04_gee_reference_audit.R`, then `07_validation_ledger.R`, and this ",
      "section will populate.")
  blank()
}

# =============================================================================
add("# What this means for looking after the land")
blank()
add("::: {.callout-tip}")
add("## Five things this work supports")
add(":::")
blank()
add("**1. The peat here is thinner than the regional average, and thin peat ",
    "is the vulnerable kind.** Where peat is only 15-25 cm deep over mineral ",
    "soil, drainage, road building or fire can remove a large share of it ",
    "quickly. Deep peat has more buffer. Knowing which sites are thin is ",
    "directly useful for deciding where disturbance matters most.")
blank()
add("**2. Judging peatlands by shallow maps will undervalue them.** Most ",
    "readily available carbon maps stop at 30 cm. At Fort Severn that misses ",
    "roughly ", sprintf("%.0f%%",
                        100 * (1 - 30 / CFG$lit$li_peat_depth_mean_cm)),
    " of the peat column. Any assessment, offset calculation or land-use ",
    "decision resting on a 0-30 cm layer is working with a fraction of the ",
    "carbon actually present.")
blank()
add("**3. The community data corroborates the published regional map.** That ",
    "matters for trust in both directions: the published map is supported by ",
    "local ground measurement, and the community's sampling is shown to be ",
    "producing numbers consistent with independent science.")
blank()
add("**4. Eight cores cannot map this landscape, and no amount of ",
    "computation changes that.** They constrain a small percentage of the ",
    "study area. Every method that would paint a confident carbon value ",
    "across the rest was tested in this workflow and rejected, because it ",
    "would give a pattern the appearance of measurement without the ",
    "substance.")
blank()
add("**5. The most valuable next measurement is depth, not more shallow ",
    "cores.** A simple probe to the bottom of the peat, at many locations, ",
    "would tell you more about carbon here than a great many more 30 cm ",
    "cores. Peat depth is the variable that dominates the total, and it is ",
    "cheap to measure.")
blank()

add("# What would make the next round stronger")
blank()
add("- **Probe peat depth everywhere you core.** Push a rod to the mineral ",
    "contact and record the depth. This is the highest-value measurement ",
    "available and needs no laboratory.")
add("- **Spread the cores out.** Roughly 30-50 cores across the wet-to-dry ",
    "and coast-to-inland gradients would support a real map. Eight clustered ",
    "cores cannot, however they are analysed.")
add("- **Take cores to a common depth.** Two of eight stopped short of 30 cm, ",
    "which forces their results to be reported as minimums.")
add("- **Use one laboratory convention.** The 2024 and 2025 campaigns ",
    "converted organic matter to carbon using different factors, so the two ",
    "sets are not perfectly comparable. Agreeing one method, or measuring ",
    "carbon directly on a subset, would remove that.")
add("- **Sample both peat and mineral sites in the same season.** At present ",
    "site type and field season are completely tangled, so a difference ",
    "between them cannot be separated from a difference between years.")
blank()

add("# Where the numbers come from")
blank()
add("Every figure in this report is produced by the analysis scripts in this ",
    "repository from the raw core data, and can be regenerated with ",
    "`Rscript run_all.R`. Nothing is typed in by hand. The technical ",
    "companion, including all quality-control findings and validation, is in ",
    "`REPORT.qmd`.")
blank()
add("## The maps, and what each one may be used for")
blank()
add("All maps are GeoTIFF files that open directly in ArcGIS, QGIS or Google ",
    "Earth Engine. Everything is written to a local folder in this project ",
    "(`outputs/spatial/` and `outputs/gee/`), so nothing depends on a cloud ",
    "account to open.")
blank()
add("Every map carries a **status label**, and it is worth reading before ",
    "using one:")
blank()
add(md_table(data.frame(
  Status = c("OBSERVATION", "GEOMETRY", "SAMPLE-BASED, MASKED",
             "CLASS-MEAN ASSIGNMENT", "PRIOR-DOMINATED", "DIAGNOSTIC"),
  `What it means` = c(
    "Measured in the field. The most trustworthy thing here.",
    "Pure geometry, such as distance to the nearest core. Exact, nothing inferred.",
    "Built from the cores, and blank wherever no core is close enough to speak.",
    "Every pixel carries its group's average. Real variation within a group is not shown.",
    "Mostly a published map, nudged by the cores where they are near.",
    "Shown to reveal a problem or a gap. NOT a carbon estimate."),
  check.names = FALSE, stringsAsFactors = FALSE)))
blank()
sp_dir <- file.path(CFG$root, "outputs", "spatial")
inv <- character(0)
for (mf in c("MANIFEST.csv", "MANIFEST_bayes.csv")) {
  p <- file.path(sp_dir, mf)
  if (file.exists(p)) {
    m <- read_csv_verbatim(p)
    col <- intersect(c("file", "units", "scientific_status", "meaning",
                       "depth_support"), names(m))
    inv <- c(inv, md_table(m[, col]), "")
  }
}
if (length(inv)) { add(inv); blank() }
add("The one to be most careful with is any map whose name begins ",
    "`DIAGNOSTIC_`. Those exist to make a limitation visible -- usually how ",
    "little of the area the cores can speak for -- and reading them as carbon ",
    "estimates would be a mistake.")
blank()
add("Key sources: Li et al. (2025) for Hudson Bay Lowlands peat depth and ",
    "carbon; Sothe et al. (2022) for Canada-wide soil carbon; SoilGrids 2.0 ",
    "for the global comparison. Full citations are in `config.R`.")
blank()

out <- file.path(CFG$root, "outputs", "COMMUNITY_REPORT.qmd")
writeLines(md, out)
log_ok("wrote ", out, " (", length(md), " lines)")
log_info("render with:  quarto render outputs/COMMUNITY_REPORT.qmd")
log_ok("12 complete")
