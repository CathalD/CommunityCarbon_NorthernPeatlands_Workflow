# =============================================================================
# 04_train_model.R
#
# Fit one random forest on the covariate stack, validate it honestly with
# leave-one-core-out cross-validation (the only CV scheme that doesn't leak
# with this data -- see the README), and report the score plainly.
#
# With 8 cores this score will likely be weak. That is reported, not hidden
# or argued around with a tournament of competing models.
#
# INPUT   mvp/outputs/current/training_data.csv   (from 03)
# OUTPUT  mvp/outputs/current/rf_model.rds
#         mvp/outputs/current/validation_metrics.csv
#         mvp/outputs/current/variable_importance.csv
#
# REQUIRES: ranger, rsample, yardstick, dplyr, purrr, tibble
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

library(ranger)
library(rsample)
library(yardstick)
library(dplyr)
library(purrr)
library(tibble)

msg("04  train & validate random forest")

train_path <- file.path(CFG$dir_current, "training_data.csv")
if (!file.exists(train_path)) stop("Missing ", train_path, " -- run 03 first.")

training <- read.csv(train_path, stringsAsFactors = FALSE)

id_cols <- c("core_id", "campaign", "latitude", "longitude", "lon_repaired",
            "core_depth_cm", "any_bd_flagged", "n_segments",
            "stock_is_lower_bound")
predictor_cols <- setdiff(names(training), c(id_cols, "stock_kgm2", "geometry"))

# Drop all-NA bands FIRST and by name. A band with no coverage over the AOI
# arrives as all-NA; if it stayed in, complete.cases() below would discard
# every core instead of every column, and the run would fail with a confusing
# "0 rows" error rather than telling you which asset is wrong.
all_na <- predictor_cols[sapply(training[predictor_cols], function(v) all(is.na(v)))]
if (length(all_na)) {
  msg("DROPPED ", length(all_na), " predictor(s) that are NA at every core: ",
     paste(all_na, collapse = ", "))
  msg("  Those bands have no coverage over this AOI -- see step 3's output.")
  predictor_cols <- setdiff(predictor_cols, all_na)
}
predictor_cols <- predictor_cols[sapply(training[predictor_cols], is.numeric)]

model_frame <- training[, c("stock_kgm2", predictor_cols)]
n_before <- nrow(model_frame)
model_frame <- model_frame[complete.cases(model_frame), ]
if (nrow(model_frame) < n_before) {
  msg("dropped ", n_before - nrow(model_frame),
     " core(s) with a missing covariate value")
}

# A predictor that doesn't vary across the cores cannot explain anything, and
# with a sample this small every wasted column costs real resolving power.
sdev <- sapply(model_frame[predictor_cols], sd, na.rm = TRUE)
constant <- names(sdev)[!is.finite(sdev) | sdev == 0]
if (length(constant)) {
  msg("dropped ", length(constant), " zero-variance predictor(s): ",
     paste(constant, collapse = ", "))
  predictor_cols <- setdiff(predictor_cols, constant)
  model_frame <- model_frame[, c("stock_kgm2", predictor_cols)]
}

msg(nrow(model_frame), " cores, ", length(predictor_cols), " usable predictors: ",
   paste(predictor_cols, collapse = ", "))
if (!length(predictor_cols)) {
  stop("No usable predictors remain -- nothing to fit. Check steps 2 and 3.")
}

set.seed(CFG$seed)

# ---- 1. leave-one-core-out cross-validation --------------------------------
# The finest split that doesn't leak: every fold trains on all-but-one core
# and predicts the held-out one. Segment-level splitting would leak (every
# segment of a core shares one set of covariate values).

folds <- loo_cv(model_frame)

cv_predictions <- map_dfr(folds$splits, function(split) {
  fit <- ranger(stock_kgm2 ~ ., data = analysis(split),
               num.trees = CFG$rf$num_trees, mtry = CFG$rf$mtry,
               min.node.size = CFG$rf$min_node_size, seed = CFG$seed)
  held_out <- assessment(split)
  tibble(truth = held_out$stock_kgm2,
        estimate = predict(fit, held_out)$predictions)
})

cv_rmse <- rmse_vec(cv_predictions$truth, cv_predictions$estimate)
cv_r2 <- tryCatch(rsq_vec(cv_predictions$truth, cv_predictions$estimate),
                  error = function(e) NA_real_)
null_rmse <- rmse_vec(cv_predictions$truth, rep(mean(model_frame$stock_kgm2),
                                                 nrow(cv_predictions)))

msg(sprintf("leave-one-core-out RMSE: %.3f kg/m2  (mean-only null: %.3f)",
           cv_rmse, null_rmse))
msg(sprintf("leave-one-core-out R2:   %s",
           if (is.na(cv_r2)) "NA (not enough variance in predictions)" else sprintf("%.3f", cv_r2)))
if (cv_rmse >= null_rmse) {
  msg("NOTE: the model did not beat a plain mean-only null on this CV. ",
     "With ", nrow(model_frame), " cores that is an expected result, not a bug -- ",
     "read validation_metrics.csv before trusting the prediction raster in step 5.")
}

# ---- 2. final model, fit on everything -------------------------------------

final_model <- ranger(stock_kgm2 ~ ., data = model_frame,
                      num.trees = CFG$rf$num_trees, mtry = CFG$rf$mtry,
                      min.node.size = CFG$rf$min_node_size,
                      importance = "permutation", quantreg = TRUE,
                      seed = CFG$seed)

saveRDS(final_model, file.path(CFG$dir_current, "rf_model.rds"))

metrics <- tibble(
  metric = c("n_cores", "n_predictors", "loco_rmse_kgm2", "loco_r2",
            "null_rmse_kgm2", "beats_null"),
  value  = c(nrow(model_frame), length(predictor_cols), cv_rmse, cv_r2,
            null_rmse, cv_rmse < null_rmse)
)
write.csv(metrics, file.path(CFG$dir_current, "validation_metrics.csv"), row.names = FALSE)

importance <- tibble(covariate = names(final_model$variable.importance),
                     importance = final_model$variable.importance) %>%
  arrange(desc(importance))
write.csv(importance, file.path(CFG$dir_current, "variable_importance.csv"), row.names = FALSE)
print(importance)

msg("wrote rf_model.rds, validation_metrics.csv, variable_importance.csv")
msg("04 complete")
