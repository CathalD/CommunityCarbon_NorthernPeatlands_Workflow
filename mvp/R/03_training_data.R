# =============================================================================
# 03_training_data.R
#
# Extract covariate values at each core location and join to its carbon
# stock. This is the modelling frame for step 4.
#
# INPUT   mvp/outputs/current/cores_clean.geojson   (from 01)
#         mvp/outputs/current/predictors.tif        (from 02)
# OUTPUT  mvp/outputs/current/training_data.csv
#
# REQUIRES: sf, terra, dplyr
# =============================================================================

.this_dir <- Sys.getenv("MVP_R_DIR", "")
if (!nzchar(.this_dir)) {
  # Works three ways: Rscript (finds --file=), source() from ANY working
  # directory (finds the sourced file's own path), and run_all.R (sets
  # MVP_R_DIR). Only falls back to getwd() if all three fail.
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(.f)) .f <- unlist(lapply(sys.frames(), function(e) e$ofile))
  .this_dir <- if (length(.f)) dirname(normalizePath(.f[[1]])) else getwd()
}
source(file.path(.this_dir, "00_utils.R"))
source(file.path(.this_dir, "..", "config.R"))

library(sf)
library(terra)
library(dplyr)

msg("03  assemble training data")

cores_path <- file.path(CFG$dir_current, "cores_clean.geojson")
pred_path  <- file.path(CFG$dir_current, "predictors.tif")
if (!file.exists(cores_path)) stop("Missing ", cores_path, " -- run 01 first.")
if (!file.exists(pred_path))  stop("Missing ", pred_path, " -- run 02 first.")

cores <- st_read(cores_path, quiet = TRUE)
predictors <- rast(pred_path)

msg("predictor bands: ", paste(names(predictors), collapse = ", "))

extracted <- terra::extract(predictors, vect(cores))
training <- cores %>%
  st_drop_geometry() %>%
  bind_cols(extracted %>% select(-ID))

# Report NAs PER BAND, not just per core. A band that is NA at every core has
# no coverage over this AOI at all -- a wrong-asset problem, not a data gap --
# and that distinction is invisible in a per-core count.
na_by_band <- sapply(training[names(predictors)], function(v) sum(is.na(v)))
if (any(na_by_band > 0)) {
  msg("covariate bands with missing values at one or more cores:")
  for (b in names(na_by_band)[na_by_band > 0]) {
    flag <- if (na_by_band[[b]] == nrow(training)) "  <-- NA AT EVERY CORE" else ""
    msg(sprintf("   %-18s %d of %d cores%s", b, na_by_band[[b]], nrow(training), flag))
  }
  if (any(na_by_band == nrow(training))) {
    msg("A band that is NA at EVERY core has no coverage over this AOI. Check ",
       "that asset's spatial domain in the Earth Engine catalogue and swap it ",
       "in config.R -- step 4 will drop it, so you would lose that predictor ",
       "silently otherwise.")
  }
}

write.csv(training, file.path(CFG$dir_current, "training_data.csv"), row.names = FALSE)
msg("wrote ", file.path(CFG$dir_current, "training_data.csv"),
   "  (", nrow(training), " cores x ", ncol(training), " columns)")
print(training)

msg("03 complete")
