# =============================================================================
# config.R  --  Single shared configuration for the SOC mapping workflow
#
# Fort Severn / Hudson Bay Lowlands community soil carbon workflow.
# Every script sources THIS file and nothing else for paths and parameters.
# Nothing here executes analysis; it only declares.
# =============================================================================

CFG <- local({

  # Locate the project root without depending on the caller's working
  # directory. Tries, in order: an explicit override, the Rscript --file=
  # argument, the path source() is reading, then a walk up from the cwd.
  find_root <- function() {
    has_cfg <- function(d) file.exists(file.path(d, "config.R"))
    walk_up <- function(d, n = 5L) {
      for (i in seq_len(n)) {
        if (has_cfg(d)) return(d)
        p <- dirname(d)
        if (identical(p, d)) break
        d <- p
      }
      NULL
    }

    p <- Sys.getenv("CCNP_ROOT", "")
    if (nzchar(p) && has_cfg(p)) return(normalizePath(p, mustWork = TRUE))

    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(f)) {
      r <- walk_up(normalizePath(dirname(f[[1]]), mustWork = FALSE))
      if (!is.null(r)) return(normalizePath(r, mustWork = TRUE))
    }

    for (i in rev(seq_len(sys.nframe()))) {
      of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(of)) {
        r <- walk_up(normalizePath(dirname(of), mustWork = FALSE))
        if (!is.null(r)) return(normalizePath(r, mustWork = TRUE))
      }
    }

    r <- walk_up(normalizePath(getwd(), mustWork = TRUE))
    if (!is.null(r)) return(normalizePath(r, mustWork = TRUE))

    stop("config.R: cannot locate the project root (no config.R found ",
         "walking up from '", getwd(), "'). Set CCNP_ROOT to the repository ",
         "directory and retry.", call. = FALSE)
  }

  root <- find_root()

  list(

    # ---- provenance -------------------------------------------------------
    project      = "CommunityCarbon_NorthernPeatlands",
    study_area   = "Hudson Bay Lowlands, near Fort Severn, Ontario, Canada",
    seed         = 20260727L,   # every stochastic step re-seeds from this

    # ---- paths ------------------------------------------------------------
    root         = root,
    dir_raw      = file.path(root, "data", "raw"),
    dir_derived  = file.path(root, "data", "derived"),
    dir_tables   = file.path(root, "outputs", "tables"),
    dir_figures  = file.path(root, "outputs", "figures"),
    dir_gee      = file.path(root, "outputs", "gee"),
    # Everything handed to the governing council lives here and nowhere else,
    # so a delivery can be zipped by directory rather than by remembering which
    # files were meant to be shared.
    dir_deliver  = file.path(root, "deliverables"),

    # Raw community core data. Verbatim headers, one row per sampled segment.
    file_cores_raw = file.path(root, "data", "raw", "community_soil_cores.csv"),

    # Output of a PRIOR pipeline, retained only as a comparison target so that
    # 01_ingest_qc.R can demonstrate the errors it silently introduced.
    # NEVER used as an analysis input.
    file_prior_pipeline = file.path(root, "data", "raw",
                                    "prior_pipeline_soil_core_forGEE.csv"),

    # ---- expected column names in the raw file ----------------------------
    # Declared so ingest fails loudly if the supplier renames a column.
    raw_cols = c("year", "Core Id", "Sample Id", "Latitude", "Longitude",
                 "Depth", "Bulk Density", "OM", "SOC",
                 "Organic Carbon Density (g/cm^2)"),

    # ---- geolocation ------------------------------------------------------
    crs_geographic = "EPSG:4326",
    # Equal-area CRS for any distance/area work at this latitude.
    # Canada Atlas Lambert; appropriate for the Hudson Bay Lowlands.
    crs_equal_area = "EPSG:3979",

    # Fort Severn, ON. Used ONLY to test the stated location against the file.
    site_anchor = c(lon = -87.6333, lat = 56.0167),

    # Region in which a Fort Severn core is physically possible. A core outside
    # this box is a geolocation failure, not an outlier. Deliberately generous
    # (~+/- 60 km) so it catches sign/transposition errors, not GPS scatter.
    plausible_bbox = c(xmin = -88.60, ymin = 55.50, xmax = -86.90, ymax = 56.60),

    # Mapping extent for all raster products (GEE side).
    aoi_bbox = c(xmin = -88.05, ymin = 55.88, xmax = -87.45, ymax = 56.22),

    # Longitude sign repair. The supplied file stores western longitudes as
    # POSITIVE, which places the cores in Siberia. 01_ingest_qc.R HARD FAILS on
    # this unless this flag is TRUE, in which case it repairs the sign LOUDLY
    # and stamps every affected row with a provenance flag.
    # Set to FALSE to re-experience the failure.
    allow_longitude_sign_fix = TRUE,

    # ---- depth handling ---------------------------------------------------
    # The `Depth` column is segment THICKNESS in cm, not cumulative bottom
    # depth. 01_ingest_qc.R re-derives this from the data and fails if the
    # evidence contradicts it. See docs in R/depth.R.
    depth_col_semantics = "thickness_cm",

    # Reporting windows. Both are produced; they are NOT interchangeable.
    #  - "common"    : deepest window every core fully covers. The only truly
    #                  like-for-like comparison across all 8 cores. Computed
    #                  from the data, not hard-coded; this is the fallback.
    #  - "reference" : 0-30 cm, chosen to match SoilGrids OCS 0-30 and the
    #                  other published layers. Not every core reaches it, so
    #                  short cores are flagged as LOWER BOUNDS.
    window_reference_cm = c(0, 30),
    window_common_cm    = NULL,   # NULL => derive from data in 02

    # Optional, off by default. If TRUE, 02 additionally reports a variant in
    # which short cores are extended to 30 cm by carrying the deepest observed
    # carbon density downward. This is EXTRAPOLATION and is written to separate
    # columns; it never enters the default product.
    report_extrapolated_variant = TRUE,

    # ---- organic matter -> organic carbon ---------------------------------
    # These are NOT applied to the primary product. The supplied SOC column is
    # used as reported. These constants exist so 01 can identify which factor
    # the data provider used, and so 05 can run a sensitivity analysis.
    factor_van_bemmelen = 0.58,   # 1/1.724. Classic, and known to overestimate
                                  # carbon in true peat.
    factor_peat_lit     = 0.51,   # Typical measured C:OM for northern peat.
                                  # Used for sensitivity on the PM cores only.
    factor_match_tol    = 0.01,   # tolerance when identifying the factor used

    # ---- QC thresholds (flag, never drop) ---------------------------------
    qc = list(
      # Dry bulk density, g/cm3.
      bd_min            = 0.02,   # below this is implausible even for moss
      bd_max            = 2.00,   # above this approaches parent rock density
      bd_peat_max       = 0.40,   # peat (OM > om_peat_threshold) above this is
                                  # internally inconsistent
      om_peat_threshold = 30.0,   # % OM at/above which a sample reads as peat
      om_min            = 0.0,
      om_max            = 100.0,
      soc_max_pct       = 60.0,   # pure organic carbon ceiling for soil
      thickness_min_cm  = 0.1,
      # Relative tolerance when reconciling the supplied Organic Carbon Density
      # column against a first-principles recomputation.
      ocd_rel_tol       = 1e-6
    ),

    # ---- statistics -------------------------------------------------------
    boot_n         = 10000L,   # bootstrap replicates (seeded)
    boot_conf      = 0.95,
    # Cross-validation is ALWAYS at the core level. Segment-level CV would leak
    # (segments within a core share location, campaign, operator and lab batch).
    cv_unit        = "core",

    # ---- local spatial output --------------------------------------------
    # Products written by 09_spatial_products.R without Earth Engine, using the
    # base-R GeoTIFF writer in R/geotiff.R.
    spatial = list(
      cell_m    = 100L,        # local raster cell size, metres
      epsg_out  = 4326L        # geographic; avoids needing a reprojection lib
    ),

    # ---- Google Earth Engine ---------------------------------------------
    # Asset IDs below are taken from the project's `Prior_data` catalogue
    # (Carbon Resources GEE reference library). Access is marked because a
    # PRIVATE asset must be shared with the running account or 04 will fail
    # with a permissions error rather than a data error.
    gee = list(
      # Set to your own Cloud project before running 03/04.
      project = Sys.getenv("EE_PROJECT", "ee-cathalpdoherty2"),

      # --- soil carbon reference layers -----------------------------------
      # SoilGrids 2.0 is built from CONCENTRATION + BULK DENSITY rather than
      # read from the packaged ocs_mean asset. Doing the integration here means
      # the depth intervals and the unit conversion are explicit and auditable,
      # and it yields 0-30 cm and 0-100 cm on the same basis.
      #   soc_mean  is dg/kg  (divide by 10  -> g/kg)
      #   bdod_mean is cg/cm3 (divide by 100 -> g/cm3)
      #   layer OCS (kg C/m2) = (soc/10) * (bdod/100) * thickness_cm / 100
      asset_soilgrids_soc  = "projects/soilgrids-isric/soc_mean",       # PUBLIC
      asset_soilgrids_bdod = "projects/soilgrids-isric/bdod_mean",      # PUBLIC

      # Sothe et al. 2022. TWO depth supports of the SAME product. Never mix.
      asset_sothe_sc_0_30  = paste0("projects/ee-cathalpdoherty2/assets/",
                                    "McMaster_WWFCanada_soil_carbon30cm"), # PRIVATE
      asset_sothe_sc_0_100 = "projects/sat-io/open-datasets/carbon_stocks_ca/sc", # PUBLIC
      asset_sothe_sc_unc   = paste0("projects/carbon-learning-library/assets/",
                                    "McMasterWWFCanadasoilcarbon1muncertainty250mkgm2version3"), # PRIVATE

      # Li et al. 2025, Hudson Bay Lowlands. FULL PEAT COLUMN carbon, kg C/m2.
      # The "30m" in the asset name is PIXEL SIZE, not depth. 04 verifies this
      # against the pixel grid and the published values in CFG$lit rather than
      # trusting either the name or this comment.
      asset_li2025         = paste0("projects/ee-cathalpdoherty2/assets/",
                                    "McMasterCarbon30mkgm2version1"),     # PRIVATE
      asset_li2025_unc     = paste0("projects/ee-cathalpdoherty2/assets/",
                                    "McMasterUncertaintyCarbon30mkgm2"),  # PRIVATE

      # --- comparison context ----------------------------------------------
      # Merged soil/peat point profiles (WOSIS 2023 + CanPeat + Janousek),
      # keyed by a 'dataset' property. Used by 10 to place these 8 cores in a
      # national and global distribution.
      asset_profiles_combined = paste0("projects/north-star-project-470316/",
                                       "assets/combined_profiles"),       # PRIVATE
      asset_profiles_canpeat  = paste0("projects/north-star-project-470316/",
                                       "assets/peat_profiles"),           # PRIVATE
      asset_profiles_wosis    = paste0("projects/north-star-project-470316/",
                                       "assets/wosis_layers_canada"),     # PRIVATE

      # --- covariates. Rationale for each is in 03_gee_covariates.R --------
      asset_dem_arctic    = "UMN/PGC/ArcticDEM/V3/2m_mosaic",
      asset_dem_copernicus= "COPERNICUS/DEM/GLO30",
      asset_s1            = "COPERNICUS/S1_GRD",
      asset_s2            = "COPERNICUS/S2_SR_HARMONIZED",
      asset_worldcover    = "ESA/WorldCover/v200",
      asset_jrc_water     = "JRC/GSW1_4/GlobalSurfaceWater",
      asset_alphaearth    = "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL",

      # GWL_FCS30 (Zhang et al. 2023): a 30 m wetland map with a fine class
      # system. Preferred over WorldCover for defining peatland strata here,
      # because WorldCover cannot separate treed bog on peat from upland forest
      # on mineral soil, which is precisely the distinction this study needs.
      asset_gwl_fcs30     = "projects/sat-io/open-datasets/GWL_FCS30",
      # 180 non-wetland | 181 permanent water | 182 swamp | 183 marsh
      # 184 flooded flat | 185 saline | 186 mangrove | 187 salt marsh
      # 188 tidal flat.  Wetland codes exclude 180 (dry) and 181 (open water).
      gwl_wetland_codes   = c(182L, 183L, 184L, 185L, 186L, 187L, 188L),

      # Growing-season window for optical/radar compositing at 56 N.
      season_start_doy = 182L,   # 1 July
      season_end_doy   = 244L,   # 31 August
      season_years     = 2023:2025,
      s2_cloud_pct     = 40,

      # Output grid and sampling extent.
      # The GEE AOI is expanded around the actual core locations so covariate
      # composites cover the whole sampled neighbourhood rather than only the
      # hand-drawn map box.
      aoi_buffer_km    = 50,
      export_scale_m   = 30L,
      # The 19-band predictor stack is the one export where 30 m is a poor
      # trade: over the AOI envelope it is ~1.2 GB and 64 separate requests,
      # and nothing downstream consumes it while 06 declines to map. 100 m
      # keeps it inspectable at a fraction of the cost. Set to
      # CFG$gee$export_scale_m if you want the full-resolution stack.
      predictor_export_scale_m = 100L,
      # EPSG:4326 keeps every product in this workflow on one grid and avoids
      # the server-side reprojection that a projected CRS forces on large
      # exports. Area calculations that need equal-area are done in EPSG:3979
      # locally, where reprojecting a finished raster is cheap.
      export_crs       = "EPSG:4326"
    ),

    # ---- oriented AOI (script 13) -----------------------------------------
    # An axis-aligned box fits a coastal study area badly. The Hudson Bay
    # shoreline at Fort Severn runs NW-SE and the sampling follows it, so the
    # AOI is a rectangle aligned to that trend with a generous buffer.
    aoi_oriented = list(
      # NULL means "derive the bearing from the cores' principal axis", which
      # for a shoreline-following transect recovers the shoreline trend.
      # Override with a number (degrees from north) to set it by hand.
      bearing_deg = NULL,
      along_buffer_km  = 30,   # beyond the cores at EACH end of the long axis
      across_buffer_km = 12,   # beyond the cores on EACH side
      # Below this share of variance on the first principal component the
      # cores are too blob-like for an oriented rectangle to mean anything,
      # and 13 falls back to an axis-aligned box rather than inventing a trend.
      min_var_explained = 0.60
    ),

    # ---- National Pedon Database (script 14) ------------------------------
    npdb = list(
      dir = "data/raw/National Pedon Database/outputs_npdb",
      file_cores   = "npdb_carbon_cores.csv",
      file_samples = "npdb_carbon_samples.csv",
      # NPDB stocks are Mg C/ha. This workflow reports kg C/m2 throughout.
      # 1 kg/m2 = 10 Mg/ha, so Mg/ha / 10 = kg/m2. Applied explicitly in 14
      # and checked against a plausible range rather than assumed.
      mgha_to_kgm2 = 0.1,
      # Cores whose recorded position is worse than this are flagged: at
      # coarse location a core cannot be attributed to a landscape unit.
      max_location_conf_m = 1000,
      # Comparison rings around Fort Severn, km.
      context_radii_km = c(500, 1000, 2000)
    ),

    # ---- landscape clustering (script 16) ---------------------------------
    # k-means on the AlphaEarth embedding. Unsupervised, so unlike the
    # predictive model in 06 this is NOT constrained by having eight cores:
    # nothing about carbon enters the fit, and the sample size that matters is
    # the pixel count.
    clusters = list(
      k_values = c(3L, 4L, 5L, 6L, 8L, 10L),  # range explored and reported
      k_final  = 6L,                          # k used for the reported product
      n_training_pixels = 5000L
    ),

    # ---- community deliverables (scripts 17-19) ---------------------------
    # The council asked for EPSG:3978 (NAD83 / Canada Atlas Lambert). Note this
    # is NOT the same code as `crs_equal_area` above: 3978 and 3979 share every
    # projection parameter and differ only in datum realisation (NAD83 versus
    # NAD83(CSRS)). The difference is sub-metre in Canada, but the two are not
    # interchangeable in metadata, so both are named explicitly rather than one
    # being used for the other.
    deliver = list(
      epsg   = 3978L,        # deliverable CRS, as specified by the council
      cell_m = 100L          # deliverable raster resolution, projected metres
    ),

    # ---- predictive map (script 18) ---------------------------------------
    # Step 2 of the data-sharing plan. The settings below are STATED CHOICES,
    # fixed before any result was seen, because with eight cores almost any
    # tuning decision can be made to improve a validation statistic after the
    # fact.
    step2 = list(
      num_trees     = 500L,
      # ranger's regression default is 5, which on 8 observations permits at
      # most one split in a tree. That is a threshold, not a forest, so the
      # minimum node size is lowered and the consequence is reported.
      min_node_size = 2L,
      max_depth     = 8L,
      mtry          = NULL,   # NULL -> ceiling(p / 3), the regression convention

      # Permutation importance settings. Expensive by design: a permutation
      # null is the only defensible way to rank covariates at this sample size.
      imp_n_perm = 10L,
      imp_n_null = 50L,
      # Wall-clock budget for the importance step. Comfortable with ranger;
      # without it the base-R forest would run for hours, so the draw count is
      # trimmed to fit and the reduced p-value floor is reported. Raise this for
      # the delivery run.
      imp_budget_sec = 600,

      # Spatial cross-validation. Each held-out core also removes every
      # training core within this radius, so a core is never predicted from a
      # near-neighbour that is effectively itself. The default is the median
      # core spacing derived in 09; NULL means "derive it".
      cv_buffer_km = NULL,

      # PRE-REGISTERED SELECTION RULE. A candidate may be mapped only if it
      # clears all three. Written here, in config, so it cannot be relaxed in
      # the script after the numbers come back.
      require_loco_r2   = 0,    # must beat the observed mean, held out
      require_loso_r2   = 0,    # must transfer between peat and mineral
      require_beat_prior = TRUE # must beat simply handing over the published map
    ),

    # ---- Bayesian map (script 11) -----------------------------------------
    # A published carbon map is the PRIOR; the community cores are the
    # observations. See R/bayes.R for the model. Nothing here is fitted from
    # the data -- every value is a stated choice, so that the weighting is
    # auditable rather than inferred from eight points.
    bayes = list(
      # Distance kernel. The length scale defaults to the median core-to-core
      # spacing, derived from the data in 09 (2.51 km): the scale at which this
      # sampling design resolves anything at all. NULL means "derive it".
      length_scale_km = NULL,
      # Beyond this, a core's weight is set to exactly zero. Without a hard
      # cut-off the map looks informed everywhere, which it is not.
      max_influence_km = 7.5,
      # Sensitivity length scales, reported alongside the default.
      sensitivity_length_scales_km = c(1, 2.5, 5, 10),

      # A mineral core carries no information about a peat pixel. Gating
      # requires a stratum raster; without one the maps are produced per
      # stratum instead, for the analyst to combine against their own layer.
      stratum_gating = TRUE,

      # Share of FULL-COLUMN peat carbon lying in the top 30 cm. Needed only
      # for the Li-prior product, and unavoidable when updating a full-column
      # prior with 30 cm cores. Under uniform density this would be
      # 30/184 = 0.163, but peat compacts: these cores show bulk density rising
      # more than fourfold from the surface to 30 cm, so the top of the column
      # holds LESS than its proportional share. 0.12 is a considered default
      # between the uniform value and the compaction-adjusted end of the range.
      shallow_fraction_of_column = 0.12,
      shallow_fraction_sensitivity = c(0.08, 0.12, 0.163),

      # Organic-matter threshold defining peat, for peat-depth detection.
      om_peat_threshold_pct = 30
    ),

    # ---- published values used as AUDIT TARGETS ---------------------------
    # These are not inputs to any estimate. 04 checks each reference layer's
    # measured statistics against them, so that a mis-specified asset (wrong
    # band, wrong scaling, wrong depth basis) is caught by disagreeing with its
    # own publication rather than by looking plausible.
    lit = list(
      # Li et al. 2025, GRL, doi:10.1029/2024GL110679
      li_peat_depth_mean_cm   = 184,
      li_peat_depth_sd_cm     = 48,
      li_peat_stock_mean_kgm2 = 86,
      li_peat_stock_sd_kgm2   = 35,
      # 1 kg C/m2 = 10 Mg C/ha throughout.
      kgm2_per_Mgha           = 0.1
    ),

    # ---- reference-layer safety ------------------------------------------
    # Physically plausible ranges (kg C / m2) used by 04 to AUDIT rather than
    # assume the depth support and units of each published layer.
    # A 0-30 cm mineral/peat stock cannot plausibly exceed ~60 kg/m2.
    # A full-profile HBL peat stock is typically ~50-400 kg/m2 (Li et al. 2025
    # report a mean of 86 +/- 35).
    ref_plausible = list(
      stock_0_30_kgm2_max  = 60,
      full_profile_kgm2_min = 50
    ),

    # ---- citations --------------------------------------------------------
    citations = c(
      sothe2022 = paste0("Sothe, C., Gonsamo, A., Arabian, J., Kurz, W.A., ",
                         "Finkelstein, S.A., & Snider, J. (2022). Large soil ",
                         "carbon storage in terrestrial ecosystems of Canada. ",
                         "Global Biogeochemical Cycles, 36(2), e2021GB007213. ",
                         "doi:10.1029/2021GB007213"),
      li2025    = paste0("Li, Y., Gonsamo, A., Han, D., Rogers, C.A., ",
                         "Finkelstein, S.A., Hararuk, O., Waddington, J.M., ",
                         "Barreto, C., McLaughlin, J.W., & Snider, J. (2025). ",
                         "Peat Depth and Carbon Storage of the Hudson Bay ",
                         "Lowlands, Canada. Geophysical Research Letters, 52. ",
                         "doi:10.1029/2024GL110679"),
      zhang2023 = paste0("Zhang, X., et al. (2023). GWL_FCS30: a global 30 m ",
                         "wetland map with a fine classification system. ",
                         "Earth System Science Data, 15, 265-293. ",
                         "doi:10.5194/essd-15-265-2023"),
      soilgrids = paste0("Poggio, L., et al. (2021). SoilGrids 2.0: producing ",
                         "soil information for the globe with quantified ",
                         "spatial uncertainty. SOIL, 7, 217-240. ",
                         "doi:10.5194/soil-7-217-2021")
    )
  )
})

# --- derived directories exist ------------------------------------------------
for (.d in c(CFG$dir_derived, CFG$dir_tables, CFG$dir_figures, CFG$dir_gee)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}
rm(.d)

# --- load the pure-function layer ---------------------------------------------
# Sourcing config.R is therefore sufficient to set up a script.
for (.f in sort(list.files(file.path(CFG$root, "R"), pattern = "[.]R$",
                           full.names = TRUE))) {
  source(.f)
}
rm(.f)

invisible(CFG)
