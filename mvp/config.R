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
  dir_current = file.path(.mvp_root, "outputs", "current"),
  dir_versions = file.path(.mvp_root, "outputs", "versions"),

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
    # Copernicus GLO30, NOT ArcticDEM. ArcticDEM covers land north of 60N
    # (plus Alaska, Greenland, Kamchatka); Fort Severn is at 56N and falls
    # OUTSIDE it, so ArcticDEM returns all-NoData here -- silently costing you
    # elevation and slope, which are the strongest predictors in a landscape
    # where centimetres of relief decide peat vs. mineral. GLO30 covers
    # 60S-85N. This is an ImageCollection, so step 02 mosaics it.
    asset_dem       = "COPERNICUS/DEM/GLO30",
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

    # AOI + export. Smaller buffer and coarser scale than the original
    # (50 km / 30 m) -- an MVP-sized covariate stack that downloads in
    # under a minute rather than requiring tiled batch requests.
    aoi_buffer_km = 15,
    scale_m       = 50L,
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
