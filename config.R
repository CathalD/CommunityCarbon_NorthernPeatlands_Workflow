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

    # ---- Google Earth Engine ---------------------------------------------
    gee = list(
      # Set to your own Cloud project before running 03/04.
      project = Sys.getenv("EE_PROJECT", "ee-cathalpdoherty2"),

      # Published reference SOC layers.
      asset_soilgrids_ocs = "projects/soilgrids-isric/ocs_mean",
      asset_sothe2022     = paste0("projects/ee-cathalpdoherty2/assets/",
                                   "McMaster_WWFCanada_soil_carbon30cm"),
      asset_li2025        = paste0("projects/ee-cathalpdoherty2/assets/",
                                   "McMasterCarbon30mkgm2version1"),

      # Covariate sources. Rationale for each is in 03_gee_covariates.R.
      asset_dem_arctic    = "UMN/PGC/ArcticDEM/V3/2m_mosaic",
      asset_dem_copernicus= "COPERNICUS/DEM/GLO30",
      asset_s1            = "COPERNICUS/S1_GRD",
      asset_s2            = "COPERNICUS/S2_SR_HARMONIZED",
      asset_worldcover    = "ESA/WorldCover/v200",
      asset_jrc_water     = "JRC/GSW1_4/GlobalSurfaceWater",
      asset_alphaearth    = "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL",

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
      export_crs       = "EPSG:3979"
    ),

    # ---- reference-layer safety ------------------------------------------
    # Physically plausible ranges (kg C / m2) used by 04 to AUDIT rather than
    # assume the depth support and units of each published layer.
    # A 0-30 cm mineral/peat stock cannot plausibly exceed ~60 kg/m2.
    # A full-profile HBL peat stock is typically ~50-400 kg/m2.
    ref_plausible = list(
      stock_0_30_kgm2_max  = 60,
      full_profile_kgm2_min = 50
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
