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
  .args <- commandArgs(trailingOnly = FALSE)
  .this_file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  .this_dir <- if (length(.this_file)) dirname(.this_file[[1]]) else getwd()
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

n_incomplete <- sum(!complete.cases(training %>% select(names(predictors))))
if (n_incomplete > 0) {
  msg("WARNING: ", n_incomplete, " core(s) have at least one NA covariate ",
     "(likely a core sits outside the covariate raster's coverage). ",
     "These rows are kept but will be dropped by ranger at fit time.")
}

write.csv(training, file.path(CFG$dir_current, "training_data.csv"), row.names = FALSE)
msg("wrote ", file.path(CFG$dir_current, "training_data.csv"),
   "  (", nrow(training), " cores x ", ncol(training), " columns)")
print(training)

msg("03 complete")
