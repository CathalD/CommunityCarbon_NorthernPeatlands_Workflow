# =============================================================================
# config.R  --  every knob for the MVP pipeline, in one place.
#
# This is the ONLY file the pipeline's behaviour should be tuned from. If a
# step script has a hard-coded number that you find yourself wanting to
# change, it belongs here instead.
#
# Every step script sources R/00_utils.R THEN this file, in that order --
# find_mvp_root() must already be in scope when this file runs.
# =============================================================================

.mvp_root <- find_mvp_root()

CFG <- list(

  # ---- paths --------------------------------------------------------------
  root        = .mvp_root,
  repo_root   = dirname(.mvp_root),
  # Reuses the same raw core CSV as the original workflow -- one source of
  # truth for the field data, not a copy.
  file_cores_raw = file.path(dirname(.mvp_root), "data", "raw",
                             "community_soil_cores.csv"),
  # Hand-drawn, coast-following study area. Preferred over a computed buffer
  # around the cores: it follows the shoreline trend instead of spending most
  # of its area over open water. Step 02 uses this if present and falls back
  # to a buffered convex hull of the cores if it is missing.
  file_aoi = file.path(.mvp_root, "data", "aoi.geojson"),
  dir_current = file.path(.mvp_root, "outputs", "current"),
  dir_versions = file.path(.mvp_root, "outputs", "versions"),
  dir_figures = file.path(.mvp_root, "outputs", "current", "figures"),

  seed = 20260727L,

  # ---- site -----------------------------------------------------------------
  # Fort Severn, ON. The raw file stores longitude as POSITIVE, which places
  # the cores in Siberia; step 1 negates it. This anchor is only used to
  # sanity-check that correction, not to move any data.
  site_anchor = c(lon = -87.6333, lat = 56.0167),
  crs_geographic = "EPSG:4326",
  crs_equal_area = "EPSG:3979",   # Canada Atlas Lambert; used for buffering/distance

  # ---- carbon stock (step 1) ------------------------------------------------
  # One depth window: however deep each core actually was measured. Cores
  # shallower than this are flagged as lower bounds, never rescaled.
  reference_depth_cm = 30,

  qc = list(
    bd_min = 0.02,   # g/cm3, below this is implausible even for moss
    bd_max = 2.00    # g/cm3, above this approaches parent-rock density
  ),

  # ---- covariates + priors (step 2, needs Earth Engine) ----------------------
  gee = list(
    project = Sys.getenv("EE_PROJECT", "ee-cathalpdoherty2"),

    # Covariate stack. Trimmed from the original 6-source stack to the two
    # strongest single-band signals for this landscape (see
    # scripts/03_gee_covariates.R in the original repo for the full
    # rationale) plus one moderate one -- add more back later if the model
    # needs them.
    # TWO DEMs, deliberately. ArcticDEM is 2 m and better where it exists, but
    # it has gaps around Fort Severn -- using it alone returns a masked band
    # and costs you elevation AND slope. Copernicus GLO30 is the base layer and
    # ArcticDEM is mosaicked on top, so ArcticDEM wins where present and
    # Copernicus fills the holes. Step 02 pins the working resolution before
    # any terrain operation; see the comment there, it is load-bearing.
    asset_dem_copernicus = "COPERNICUS/DEM/GLO30",
    asset_dem_arctic     = "UMN/PGC/ArcticDEM/V3/2m_mosaic",

    # Topographic position radii, metres. On a plain with metres of relief over
    # tens of kilometres, absolute elevation says little and POSITION says a
    # lot -- but only if you name the scale. 300 m captures ridge-and-flark
    # microtopography; 2 km captures the broad peat plateaus.
    tpi_radii_m = c(small = 300, large = 2000),
    asset_s1        = "COPERNICUS/S1_GRD",
    asset_s2        = "COPERNICUS/S2_SR_HARMONIZED",
    asset_jrc_water = "JRC/GSW1_4/GlobalSurfaceWater",

    # Published priors. Both private GEE assets already shared with the
    # project's account (see original repo's config.R for provenance).
    asset_li2025        = paste0("projects/ee-cathalpdoherty2/assets/",
                                 "McMasterCarbon30mkgm2version1"),
    asset_sothe_sc_0_30 = paste0("projects/ee-cathalpdoherty2/assets/",
                                 "McMaster_WWFCanada_soil_carbon30cm"),

    # Growing-season window for optical/radar compositing at 56N.
    season_start_doy = 182L,   # 1 July
    season_end_doy   = 244L,   # 31 August
    season_years     = 2023:2025,
    s2_cloud_pct     = 40,

    # Fallback buffer, used ONLY if file_aoi above is missing.
    aoi_buffer_km = 15,

    # 100 m, not 30 m. The supplied AOI spans roughly 106 x 107 km; at 30 m
    # that is ~12 million pixels per band, which is what forced the original
    # workflow into tiled, byte-budgeted batch requests. At 100 m the whole
    # stack is ~1 million pixels per band and exports in one go. This
    # landscape has metres of relief over tens of kilometres, so 30 m detail
    # is not the signal -- drop this to 30L only if you have a reason.
    scale_m       = 100L,
    export_crs    = "EPSG:4326"
  ),

  # ---- random forest (step 4) -----------------------------------------------
  rf = list(
    num_trees     = 500L,
    mtry          = NULL,        # NULL => ranger default
    min_node_size = 5L,
    quantiles     = c(0.1, 0.9)  # for quantile-forest prediction intervals
  ),

  # ---- Bayesian fusion (step 6) ----------------------------------------------
  bayes = list(
    # "external"         : prior = a published map (run 1)
    # "previous_version" : prior = the last version's posterior mean/sd
    prior_source = "external",

    # Which published map to use when prior_source == "external".
    #   "li2025"      : full peat-column carbon (Li et al. 2025). Needs the
    #                   depth split below because the cores measure 0-30 cm,
    #                   not the whole column.
    #   "sothe_0_30"  : already a 0-30 cm product; no depth split needed.
    external_prior = "li2025",

    # Li et al. is full-column; the cores are 0-30 cm. This is the fraction
    # of the full column assumed to sit in the top 30 cm, used ONLY to make
    # the prior comparable to the cores before fusing. The original
    # workflow's own cores implied ~11.6% here (vs. ~16.3% under a naive
    # uniform-density assumption) -- 0.15 splits the difference. Tune this
    # if you have a better local estimate.
    shallow_fraction_of_column = 0.15,

    # Prior uncertainty, spatially constant for the MVP (Li et al.'s
    # published regional uncertainty is +/-35 kg C/m2). A per-pixel
    # uncertainty raster is a natural v2 upgrade -- see README.
    prior_sd_kgm2 = 35,

    # Assumed measurement uncertainty per core. This is a real, if simple,
    # scientific assumption -- tune it down if you trust the core lab
    # values more than this.
    core_sd_kgm2 = 1.0,

    # Distance kernel. NULL => derive as the median nearest-neighbour
    # distance between cores at run time.
    length_scale_km  = NULL,
    max_influence_km = 15
  )
)
