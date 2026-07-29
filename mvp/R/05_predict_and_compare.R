# =============================================================================
# 05_predict_and_compare.R
#
# Predict the random forest across the whole covariate raster, with a
# quantile-forest uncertainty band. Load the current prior (published map on
# run 1, last version's posterior on every run after) and compute, at each
# core, how far off the prior was -- the residuals step 6 will fuse in.
#
# INPUT   mvp/outputs/current/rf_model.rds        (from 04)
#         mvp/outputs/current/predictors.tif      (from 02)
#         mvp/outputs/current/prior_*.tif          (from 02, run 1 only)
#         mvp/outputs/current/training_data.csv   (from 03)
# OUTPUT  mvp/outputs/current/carbon_prediction.tif      RF surface
#         mvp/outputs/current/carbon_uncertainty.tif     quantile-forest spread
#         mvp/outputs/current/prior_shallow_equiv.tif    prior, made 0-30cm-comparable
#         mvp/outputs/current/prediction_residuals.csv   per-core: observed - prior
#
# REQUIRES: terra, ranger, dplyr
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .this_file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  .this_dir <- if (length(.this_file)) dirname(.this_file[[1]]) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

library(terra)
library(ranger)
library(dplyr)

msg("05  predict landscape + compare to prior")

model_path <- file.path(CFG$dir_current, "rf_model.rds")
pred_path  <- file.path(CFG$dir_current, "predictors.tif")
train_path <- file.path(CFG$dir_current, "training_data.csv")
if (!file.exists(model_path)) stop("Missing ", model_path, " -- run 04 first.")

rf_model   <- readRDS(model_path)
predictors <- rast(pred_path)
training   <- read.csv(train_path, stringsAsFactors = FALSE)

# ---- 1. RF prediction + quantile-forest uncertainty ------------------------

predict_fun <- function(model, data, ...) predict(model, data)$predictions
carbon_prediction <- terra::predict(predictors, rf_model, fun = predict_fun, na.rm = TRUE)
names(carbon_prediction) <- "carbon_prediction_kgm2"

quantile_fun <- function(model, data, ...) {
  q <- predict(model, data, type = "quantiles", quantiles = CFG$rf$quantiles)$predictions
  (q[, 2] - q[, 1]) / 2   # half-width of the 10-90% interval, as an sd-like uncertainty
}
carbon_uncertainty <- terra::predict(predictors, rf_model, fun = quantile_fun, na.rm = TRUE)
names(carbon_uncertainty) <- "carbon_uncertainty_kgm2"

writeRaster(carbon_prediction, file.path(CFG$dir_current, "carbon_prediction.tif"),
           overwrite = TRUE)
writeRaster(carbon_uncertainty, file.path(CFG$dir_current, "carbon_uncertainty.tif"),
           overwrite = TRUE)
msg("wrote carbon_prediction.tif, carbon_uncertainty.tif")

# ---- 2. load the current prior ---------------------------------------------

if (CFG$bayes$prior_source == "external") {

  if (CFG$bayes$external_prior == "li2025") {
    li_path <- file.path(CFG$dir_current, "prior_li2025.tif")
    if (!file.exists(li_path)) stop("Missing ", li_path, " -- run 02 first.")
    prior_full <- rast(li_path)
    # Li et al. is FULL PEAT COLUMN; the cores measure 0-30 cm. Multiplying
    # by shallow_fraction_of_column makes the two comparable. This is a real
    # scientific assumption, not ceremony -- see config.R for its basis.
    prior_shallow <- prior_full * CFG$bayes$shallow_fraction_of_column
    msg("prior: Li et al. 2025 (full column) x shallow_fraction_of_column = ",
       CFG$bayes$shallow_fraction_of_column)
  } else {
    sothe_path <- file.path(CFG$dir_current, "prior_sothe_0_30.tif")
    if (!file.exists(sothe_path)) stop("Missing ", sothe_path, " -- run 02 first.")
    prior_shallow <- rast(sothe_path)
    msg("prior: Sothe et al. 2022, 0-30cm (no depth split needed)")
  }

} else {
  # prior_source == "previous_version": fuse against last run's posterior,
  # not the original published map. This is what makes the pipeline
  # iterative rather than a one-shot rebuild.
  latest <- latest_version_dir(CFG$dir_versions)
  if (is.null(latest)) {
    stop("CFG$bayes$prior_source is 'previous_version' but no versions exist ",
        "yet in ", CFG$dir_versions, ". Run once with prior_source='external' ",
        "first (creates v1), then switch this to 'previous_version'.")
  }
  prior_shallow <- rast(file.path(latest, "carbon_posterior_mean.tif"))
  msg("prior: previous version's posterior mean, from ", latest)
}

# Resample the prior onto the prediction raster's grid -- they may come from
# different sources (a published product vs. this run's own RF output) and
# rarely share an exact grid.
prior_shallow <- resample(prior_shallow, carbon_prediction, method = "bilinear")
names(prior_shallow) <- "prior_shallow_equiv"
writeRaster(prior_shallow, file.path(CFG$dir_current, "prior_shallow_equiv.tif"),
           overwrite = TRUE)
msg("wrote prior_shallow_equiv.tif")

# ---- 3. residuals at the cores ---------------------------------------------

core_pts <- vect(training, geom = c("longitude", "latitude"), crs = CFG$crs_geographic)
prior_at_cores <- terra::extract(prior_shallow, core_pts)$prior_shallow_equiv

residuals <- training %>%
  transmute(core_id, lon = longitude, lat = latitude, campaign,
           observed = stock_kgm2,
           prior_predicted = prior_at_cores,
           residual = observed - prior_predicted)

write.csv(residuals, file.path(CFG$dir_current, "prediction_residuals.csv"), row.names = FALSE)
msg("wrote prediction_residuals.csv")
print(residuals)

msg("05 complete")
