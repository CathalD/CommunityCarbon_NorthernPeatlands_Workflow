# =============================================================================
# 06_covariate_model.R
#
# PRODUCT 3  Covariate-based prediction of SOC stock.
#            Status: PREDICTION / EXTRAPOLATION. Unlike Products 1 and 2, this
#            leans on relationships fitted to eight observations and then
#            applied to ground no core ever touched.
#
# This script is built to find out whether such a model earns the right to be
# mapped, and to report the answer either way. The default expectation, stated
# in advance so it cannot be quietly revised afterwards, is that it does not:
# eight points cannot support a covariate model that generalises, and the
# useful output of this script is a defensible demonstration of that, not a
# raster.
#
# CROSS-VALIDATION STRATEGY, AND WHY
#
#   Leave-one-CORE-out is the finest split that is honest here.
#
#   Splitting at the segment level would be catastrophic and is easy to do by
#   accident: there are 22 segments, which looks like a workable sample size.
#   But segments within a core share a location, a campaign, an operator, a
#   laboratory batch and a single set of covariate values -- at 30 m resolution
#   every segment of a core samples the SAME PIXEL. A model trained on
#   segments 1 and 2 of a core and tested on segment 3 is being tested on data
#   it has already seen. That inflates apparent skill dramatically and is the
#   single easiest way to manufacture a fraudulent validation statistic from a
#   dataset this size. This workflow never does it, and 07 records that.
#
#   Leave-one-STRATUM-out is also run, and it is the harder and more relevant
#   test. Mapping means predicting ground unlike the training points. Training
#   on the mineral cores and predicting the peat cores asks whether the model
#   learned anything transferable, or only learned to tell two clusters apart.
#
#   Every model is scored against a mean-only null. A covariate model that
#   cannot beat the mean has no business being painted across a landscape.
#
# INPUT   data/derived/02_core_stocks.csv
#         data/derived/03_covariates_at_cores.csv   (from 03; needs GEE)
# OUTPUT  data/derived/06_model_cv.csv
#         outputs/tables/06_covariate_screening.csv
#         outputs/tables/06_model_verdict.csv
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

log_step("06  COVARIATE MODEL AND HONEST CROSS-VALIDATION")

set.seed(CFG$seed)
log_info("seed: ", CFG$seed)

stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

cov_path <- file.path(CFG$dir_derived, "03_covariates_at_cores.csv")
if (!file.exists(cov_path)) {
  fail_loudly(
    "Covariates have not been extracted",
    paste0("Missing: ", cov_path),
    "Product 3 is the only product that needs Earth Engine covariates.",
    "Products 1 and 2 (script 05) are complete without them.",
    remedy = "Run scripts/03_gee_covariates.R with an authenticated Earth Engine session, then re-run this script."
  )
}
covars <- read_csv_verbatim(cov_path)
log_info("loaded covariates for ", nrow(covars), " cores")

# =============================================================================
# 1. Assemble the modelling frame
# =============================================================================

log_step("06a  MODELLING FRAME")

# The response is the common-support stock: the only window in which all eight
# cores carry a real number rather than a lower bound. Modelling the 0-30 cm
# window would mean either dropping two of eight cores or feeding lower bounds
# to a regression as if they were observations.
y_tbl <- stocks[stocks$scenario == "as_supplied" &
                stocks$window == "common_support", ]
log_info("response: common-support stock (0-",
         sprintf("%.1f", y_tbl$window_hi_cm[1]), " cm), scenario as_supplied")
log_info("chosen because all ", nrow(y_tbl),
         " cores are complete in this window; the 0-30 cm window would ",
         "require feeding two lower bounds to a regression as observations")

dat <- merge(y_tbl[, c("core_id", "campaign", "stock_kgm2")],
             covars, by = c("core_id", "campaign"))
if (nrow(dat) < nrow(y_tbl)) {
  log_warn(nrow(y_tbl) - nrow(dat), " core(s) lost in the covariate join")
}

# Ring-fence what the model may see. AlphaEarth embedding bands and the
# categorical stratum layer are excluded here, for the reasons set out in 03.
embed_cols <- grep("^A[0-9]{2}$", names(dat), value = TRUE)
# Categorical layers are matched by PATTERN, not by a literal name. The Earth
# Engine export was renamed `gwl_fcs30` -> `gwl_class` between runs, and the
# literal exclusion list here silently stopped matching -- which is how a
# wetland CLASS CODE came third in the covariate ranking below, indistinguish-
# able from a real predictor. See R/rf.R:ringfence_categorical().
.ring <- ringfence_categorical(setdiff(names(dat), embed_cols),
                               patterns = CFG$step2$categorical_patterns)
excluded   <- c(.ring$categorical, embed_cols)
if (length(.ring$suspected)) {
  log_warn("layer names that look categorical but are NOT excluded: ",
           paste(.ring$suspected, collapse = ", "))
}
candidates <- setdiff(names(dat),
                      c("core_id", "campaign", "stock_kgm2",
                        "latitude", "longitude", excluded))
candidates <- candidates[vapply(dat[candidates], is.numeric, logical(1))]
# Drop covariates with no variation across the eight points.
keep <- vapply(dat[candidates], function(v) {
  sum(is.finite(v)) >= nrow(dat) - 1L && stats::sd(v, na.rm = TRUE) > 0
}, logical(1))
candidates <- candidates[keep]

log_info(length(candidates), " candidate predictors after screening")
log_info("EXCLUDED from the model by design: ",
         length(embed_cols), " AlphaEarth embedding bands, plus ",
         length(.ring$categorical), " categorical layer(s): ",
         paste(.ring$categorical, collapse = ", "))
log_info("  64 embedding dimensions on 8 observations would fit anything and ",
         "generalise to nothing; the embedding is used for the extrapolation ",
         "diagnostic below instead")

y  <- dat$stock_kgm2
gp <- dat$core_id
st <- dat$campaign

# =============================================================================
# 2. The leakage demonstration
# =============================================================================

log_step("06b  WHY SEGMENT-LEVEL CV WOULD BE FRAUDULENT")

# Show, rather than assert, that segments of a core share covariate values.
seg <- require_artifact(file.path(CFG$dir_derived, "01_segments_qc.csv"),
                        "scripts/01_ingest_qc.R")
n_seg <- table(seg$core_id)
log_info("segments per core: ",
         paste(sprintf("%s=%d", names(n_seg), as.integer(n_seg)), collapse = " "))
log_info("22 segments looks like a usable sample size. It is not:")
log_info("  every segment of a core sits at ONE coordinate pair, so at ",
         CFG$gee$export_scale_m, " m every segment of a core samples the ",
         "SAME PIXEL of every covariate.")
log_info("  the effective sample size for a covariate model is therefore ",
         length(unique(seg$core_id)), ", not ", nrow(seg), ".")
log_warn("Segment-level cross-validation is never performed in this workflow.")

# =============================================================================
# 3. Screening, under leave-one-core-out
# =============================================================================

log_step("06c  SINGLE-COVARIATE SCREENING (leave-one-core-out)")

screen <- screen_single_covariates(dat[candidates], y, gp)
if (is.null(screen) || !nrow(screen)) {
  fail_loudly("No usable covariates survived screening",
              remedy = "Check the covariate extraction in 03.")
}
print_table(screen[, c("rank", "covariate", "rmse", "mae", "bias", "r2_cv",
                       "skill_null")], digits = 3)

n_cand <- screen$n_candidates_screened[1]
best   <- screen[1, ]
log_info(sprintf("best of %d candidates: %s (skill vs null = %+.3f)",
                 n_cand, best$covariate, best$skill_null))
log_warn("SELECTION OPTIMISM: the best of ", n_cand,
         " candidates screened on 8 points is expected to look good by")
log_warn("chance alone. This ranking is NOT a validated result; it is the ",
         "search that a validated result would have to survive.")

# =============================================================================
# 4. Head-to-head against the null, and across strata
# =============================================================================

log_step("06d  CROSS-VALIDATION")

null_cv <- cv_leave_one_group_out(dat[candidates], y, gp, fit_null, predict_null)
log_info("mean-only null model, leave-one-core-out:")
print_table(null_cv$metrics, digits = 3)

best_x  <- dat[, best$covariate, drop = FALSE]
loco    <- cv_leave_one_group_out(best_x, y, gp, fit_ols, predict_ols)
log_info("best single covariate (", best$covariate, "), leave-one-core-out:")
print_table(loco$metrics, digits = 3)

# --- two baselines, because they answer different questions -----------------
# `skill_null` in the tables above compares a CROSS-VALIDATED model against the
# IN-SAMPLE mean. That is the right question for mapping, because the mean you
# would actually paint across the landscape is computed from all the data. But
# it holds the model to a stricter standard than its baseline, so a model can
# fail it while still being better than a like-for-like null.
#
# So both are reported. The verdict rests on r2_cv, the standard test: do the
# cross-validated predictions beat simply using the observed mean? A negative
# r2_cv means they do not.
skill_lfl <- skill_vs_reference(y, loco$predictions$pred, null_cv$predictions$pred)

log_info("the two baselines, side by side:")
print_table(data.frame(
  comparison = c("model (CV) vs in-sample mean   [skill_null]",
                 "model (CV) vs cross-validated null [like for like]",
                 "model (CV) vs observed mean        [r2_cv, the standard test]"),
  value = c(loco$metrics$skill_null, skill_lfl, loco$metrics$r2_cv),
  verdict = c(ifelse(loco$metrics$skill_null > 0, "beats", "does not beat"),
              ifelse(skill_lfl > 0, "beats", "does not beat"),
              ifelse(loco$metrics$r2_cv > 0, "has skill", "NO SKILL")),
  stringsAsFactors = FALSE), digits = 3)

if (skill_lfl > 0 && loco$metrics$r2_cv <= 0) {
  log_warn("Note the split decision: the covariate lowers cross-validated RMSE ",
           sprintf("by %.1f%% against a like-for-like null,", 100 * skill_lfl))
  log_warn("yet its r2_cv is still negative, meaning the predictions remain ",
           "worse than simply using the observed mean.")
  log_warn("Both statements are true. The second is the one that governs ",
           "whether a map should be drawn, so the verdict below is NO.")
}

loso <- cv_leave_one_stratum_out(best_x, y, st, fit_ols, predict_ols)
log_info("best single covariate, leave-one-STRATUM-out:")
print_table(loso$metrics, digits = 3)
log_info("per-fold predictions:")
print_table(loso$predictions, digits = 3)

# Direction and size of the transfer failure, per stratum. A model that has
# learned nothing transferable pulls each held-out stratum toward the mean of
# whatever it was trained on, which shows up as equal and opposite bias.
lo <- loso$predictions
transfer <- do.call(rbind, lapply(unique(lo$stratum), function(s) {
  i <- lo$stratum == s & is.finite(lo$pred)
  data.frame(held_out_stratum = s, n = sum(i),
             obs_mean = mean(lo$obs[i]), pred_mean = mean(lo$pred[i]),
             bias = mean(lo$pred[i] - lo$obs[i]),
             rmse = rmse(lo$obs[i], lo$pred[i]),
             stringsAsFactors = FALSE)
}))
log_info("transfer failure by stratum:")
print_table(transfer, digits = 3)

if (nrow(transfer) == 2L && prod(sign(transfer$bias)) < 0) {
  log_warn(sprintf(paste0("The two biases have OPPOSITE signs (%+.2f and ",
                          "%+.2f kg/m2)."), transfer$bias[1], transfer$bias[2]))
  log_warn("That is the signature of a model with no transferable ",
           "relationship: trained on one stratum it drags the other toward ",
           "its own mean, over-predicting the low stratum and under-")
  log_warn("predicting the high one. The covariate is standing in for the ",
           "peat/mineral split, not for a carbon gradient that generalises.")
}

# In-sample fit, shown only to quantify the gap against cross-validated skill.
insample <- fit_ols(best_x, y)
r2_in <- summary(insample)$r.squared
log_info(sprintf("in-sample R2 for the same model: %.3f", r2_in))
log_warn(sprintf(paste0("in-sample R2 %.3f versus cross-validated R2 %.3f: ",
                        "the gap IS the overfitting"), r2_in, loco$metrics$r2_cv))

# =============================================================================
# 5. Extrapolation diagnostic using the embedding
# =============================================================================

log_step("06e  EXTRAPOLATION DIAGNOSTIC")

if (length(embed_cols) >= 8L) {
  E <- as.matrix(dat[embed_cols])
  D <- as.matrix(stats::dist(E))
  diag(D) <- NA
  nnd <- apply(D, 1, min, na.rm = TRUE)
  ext <- data.frame(core_id = dat$core_id, campaign = dat$campaign,
                    embed_nn_dist = nnd, stringsAsFactors = FALSE)
  log_info("distance in AlphaEarth embedding space to the nearest other core:")
  print_table(ext[order(-ext$embed_nn_dist), ], digits = 3)
  log_info("This is the honest use of the embedding at n=8: not as a")
  log_info("predictor, but as a ruler for how unlike each other these eight")
  log_info("points are, and therefore how far a map would be extrapolating.")
  write_csv_logged(ext, file.path(CFG$dir_tables, "06_embedding_distance.csv"),
                   "how far apart the cores are in embedding space")
} else {
  log_warn("AlphaEarth bands absent; extrapolation diagnostic skipped")
}

# =============================================================================
# 6. Verdict
# =============================================================================

log_step("06f  VERDICT: DOES PRODUCT 3 EARN A MAP?")

has_skill    <- loco$metrics$r2_cv > 0        # the governing test
beats_lfl    <- skill_lfl > 0                 # like-for-like, reported alongside
transfers    <- isTRUE(loso$metrics$r2_cv > 0)
recommend    <- has_skill && transfers

verdict <- data.frame(
  criterion = c("cross-validated predictions beat the observed mean (r2_cv > 0)",
                "lowers RMSE against a like-for-like cross-validated null",
                "transfers across strata (leave-one-stratum-out r2_cv > 0)",
                "n cores",
                "candidates screened",
                "RECOMMENDED FOR MAPPING"),
  value = c(sprintf("%s (r2_cv %+.3f)", ifelse(has_skill, "YES", "NO"),
                    loco$metrics$r2_cv),
            sprintf("%s (skill %+.3f)", ifelse(beats_lfl, "YES", "NO"), skill_lfl),
            sprintf("%s (r2_cv %+.3f)", ifelse(transfers, "YES", "NO"),
                    loso$metrics$r2_cv),
            as.character(nrow(dat)),
            as.character(n_cand),
            ifelse(recommend, "YES", "NO")),
  stringsAsFactors = FALSE
)
print_table(verdict)

log_warn(strrep("-", 68))
if (!recommend) {
  log_warn("PRODUCT 3 IS NOT RECOMMENDED FOR MAPPING.")
  log_warn(sprintf("  Cross-validated r2 is %+.3f: the predictions are worse",
                   loco$metrics$r2_cv))
  log_warn("  than simply using the observed mean. Leave-one-stratum-out ",
           sprintf("r2 is %+.3f.", loso$metrics$r2_cv))
  log_warn("  No prediction raster is written.")
  log_warn("")
  log_warn("  This is the expected outcome and it is a result, not a failure")
  log_warn("  of the pipeline. Eight cores in two clusters cannot identify a")
  log_warn("  covariate relationship that holds across a landscape. Reporting")
  log_warn("  a map anyway would give a spatial pattern the authority of")
  log_warn("  measurement while it carried none.")
  log_warn("")
  log_warn("  Use PRODUCT 1 and PRODUCT 2 from 05, and treat the published")
  log_warn("  layers as the best available spatial pattern until more cores")
  log_warn("  exist. Roughly 30-50 cores, deliberately spread across the")
  log_warn("  wetness and rebound gradients rather than clustered at two")
  log_warn("  accessible sites, would make this question answerable.")
} else {
  log_warn("The model beat the null under both tests. Given n=8 and ", n_cand,
           " candidates screened, treat this as PROVISIONAL: it is a")
  log_warn("hypothesis worth testing with new cores, not a validated map.")
  log_warn("A prediction raster is written, labelled PREDICTION.")
}
log_warn(strrep("-", 68))

# =============================================================================
# 7. Write
# =============================================================================

log_step("06g  OUTPUTS")

cv_out <- rbind(
  cbind(model = "null",                cv = "leave-one-core-out",   null_cv$metrics),
  cbind(model = paste0("ols:", best$covariate), cv = "leave-one-core-out",   loco$metrics),
  cbind(model = paste0("ols:", best$covariate), cv = "leave-one-stratum-out", loso$metrics)
)
# skill_null compares against the in-sample mean; skill_vs_cv_null compares the
# model against the cross-validated null on the same footing. Both are carried
# so a reader is never shown one without the other.
cv_out$skill_vs_cv_null <- c(0, skill_lfl, NA_real_)
write_csv_logged(cv_out, file.path(CFG$dir_derived, "06_model_cv.csv"),
                 "cross-validated skill, core level only, against both baselines")
write_csv_logged(transfer,
                 file.path(CFG$dir_tables, "06_stratum_transfer.csv"),
                 "how badly the model transfers between strata")
write_csv_logged(screen, file.path(CFG$dir_tables, "06_covariate_screening.csv"),
                 "full ranking, so selection optimism stays visible")
write_csv_logged(verdict, file.path(CFG$dir_tables, "06_model_verdict.csv"),
                 "whether Product 3 earns a map")

log_ok("06 complete")
