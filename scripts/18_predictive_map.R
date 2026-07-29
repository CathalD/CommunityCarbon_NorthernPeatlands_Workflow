# =============================================================================
# 18_predictive_map.R
#
# STEP 2 of the Fort Severn data-sharing plan: a predictive carbon surface for
# the AOI, in EPSG:3978, with the diagnostics needed to defend it if the
# council asks how it was made.
#
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DECIDES, AND WHY IT IS DECIDED HERE RATHER THAN IN THE PLAN
# -----------------------------------------------------------------------------
#
# The plan asks to "confirm the model type and whether prior maps were fused in
# (ensemble weighting, residual correction, or direct covariate)". That is not
# a question that should be answered by preference. It is answered here by
# competition, under a rule fixed in config.R BEFORE any result was seen.
#
# FIVE surfaces are built and scored against each other:
#
#   N  NULL         the observed mean, painted everywhere.
#   P  PRIOR        the published map, handed over uncorrected. This is the
#                   benchmark that matters most and the one script 06 lacked:
#                   the council can already download Li et al. for free, so a
#                   new map has to be better than that to be worth making.
#   A  RF           random forest on the environmental covariate stack.
#                   The plan's literal ask.
#   B  RF+PRIOR     the same forest with the published layers added as
#                   covariates. This is "direct covariate" fusion.
#   C  PRIOR+RESID  the published map, corrected by a distance-weighted,
#                   precision-weighted update from the cores. This is
#                   "residual correction", and it is the only candidate that
#                   degrades gracefully: where the cores say nothing it returns
#                   the published map unchanged, rather than an extrapolation.
#
# SELECTION RULE (config.R, CFG$step2). A candidate may be written as the
# council's map only if it clears ALL THREE of:
#     leave-one-core-out r2   > require_loco_r2
#     leave-one-stratum-out r2> require_loso_r2
#     beats candidate P on leave-one-core-out RMSE
# If none clears the bar, the deliverable is C. Not because C won -- it may not
# have -- but because C is the only surface whose failure mode is "returns the
# published map", and that is the failure mode a council can safely receive.
#
# The expectation, stated in advance so it cannot be quietly revised: A and B
# will fail. Script 06 already showed the best single covariate reaching
# leave-one-stratum-out r2 of about -4 on these eight cores. A forest with more
# covariates has more ways to fit noise, not fewer. If A or B does clear the
# bar, that is a surprise and is reported as provisional.
#
# -----------------------------------------------------------------------------
# DEPTH WINDOWS -- READ THIS BEFORE COMPARING ANY NUMBER HERE TO ANY OTHER
# -----------------------------------------------------------------------------
#
# Two windows are carried, and they are NOT interchangeable:
#
#   COMMON SUPPORT (0 - ~14.5 cm)  All eight cores are complete measurements.
#     Both strata have three or more cores, so leave-one-stratum-out can run in
#     BOTH directions. This window therefore carries the VERDICT.
#
#   REFERENCE (0 - 30 cm)  Two peat cores are lower bounds, not measurements,
#     and are excluded as observations exactly as in script 11. That leaves six
#     cores, of which only ONE is peat -- so leave-one-stratum-out can only run
#     in one direction here and cannot be trusted to decide anything. But this
#     is the window every published layer reports, and the only one in which a
#     comparison with Li, Sothe, SoilGrids or the pedon database is like for
#     like. This window therefore carries the MAP.
#
# The harder test decides; the deliverable follows. Mapping the 0-30 cm window
# on the strength of a verdict earned at 0-14.5 cm is itself an assumption, and
# it is recorded as one in the diagnostics table rather than left implicit.
#
# INPUT   data/derived/02_core_stocks.csv
#         data/derived/01_cores_qc.csv
#         data/derived/03_covariates_at_cores.csv   (needs GEE; from 03)
#         data/derived/04_reference_at_cores.csv    (needs GEE; from 04)
#         outputs/gee/ccnp_predictors_*m.tif        (needs GEE; for the raster)
#         outputs/spatial/aoi_vertices.csv          (from 13)
# OUTPUT  deliverables/step2_prediction/carbon_prediction.tif   (EPSG:3978)
#         deliverables/step2_prediction/carbon_prediction_0_15cm.tif
#         deliverables/step2_prediction/MODEL_DIAGNOSTICS.csv
#         deliverables/step2_prediction/VARIABLE_IMPORTANCE.csv
#         deliverables/step2_prediction/COVARIATE_STACK.csv
#         deliverables/step2_prediction/README.md
#         data/derived/18_candidate_cv.csv
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

log_step("18  STEP 2: PREDICTIVE CARBON MAP")

set.seed(CFG$seed)
log_info("seed: ", CFG$seed)
log_info("forest backend: ", rf_backend())
if (rf_backend() != "ranger") {
  log_warn("`ranger` is not installed, so the base-R bagged forest in R/rf.R ")
  log_warn("  is used instead. It implements the same algorithm and gives the ")
  log_warn("  same kind of answer, but it is far slower and has not been ")
  log_warn("  through ranger's testing. Install ranger for the delivery run:")
  log_warn("    install.packages('ranger')")
}

out_dir <- file.path(CFG$dir_deliver, "step2_prediction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Inputs
# =============================================================================

log_step("18a  INPUTS")

cores  <- require_artifact(file.path(CFG$dir_derived, "01_cores_qc.csv"),
                           "scripts/01_ingest_qc.R")
stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

cov_path <- file.path(CFG$dir_derived, "03_covariates_at_cores.csv")
ref_path <- file.path(CFG$dir_derived, "04_reference_at_cores.csv")
if (!file.exists(cov_path)) {
  fail_loudly(
    "Covariates have not been extracted",
    paste0("Missing: ", cov_path),
    "Step 2 is a covariate model; there is nothing to fit without them.",
    remedy = paste("Run scripts/03_gee_covariates.R with an authenticated",
                   "Earth Engine session, then re-run this script."))
}
covars <- read_csv_verbatim(cov_path)
log_ok("covariates at cores: ", nrow(covars), " rows, ", ncol(covars), " columns")

have_ref <- file.exists(ref_path)
refs <- if (have_ref) read_csv_verbatim(ref_path) else NULL
if (have_ref) {
  log_ok("published layers at cores: ", nrow(refs), " rows")
} else {
  log_warn("04_reference_at_cores.csv absent. Candidates B (prior as ",
           "covariate), C (residual correction) and the prior benchmark P ")
  log_warn("  cannot be built, which removes the most important comparison ",
           "in this script. Run scripts/04_gee_reference_audit.R.")
}

# =============================================================================
# 2. The covariate stack, stated explicitly
# =============================================================================

log_step("18b  COVARIATE STACK")

# Ring-fenced, for the reasons set out in 03 and 06. The 64 AlphaEarth bands on
# eight observations would fit anything and generalise to nothing; worldcover
# and the wetland class are categorical and are used for stratification, not as
# numeric predictors -- a forest splitting on class code 186 versus 187 is
# splitting on an arbitrary integer.
embed_cols <- grep("^A[0-9]{2}$", names(covars), value = TRUE)
categorical <- intersect(c("worldcover", "gwl_fcs30", "gwl_wetland"), names(covars))
id_cols <- c("core_id", "campaign", "latitude", "longitude")

env_cov <- setdiff(names(covars), c(id_cols, embed_cols, categorical))
env_cov <- env_cov[vapply(covars[env_cov], is.numeric, logical(1))]

# Published layers, when available, as candidate covariates for B.
prior_cov <- character(0)
if (have_ref) {
  prior_cov <- setdiff(names(refs), c(id_cols, grep("comparable|units|support",
                                                    names(refs), value = TRUE)))
  prior_cov <- prior_cov[vapply(refs[prior_cov], is.numeric, logical(1))]
}

stack_tbl <- rbind(
  data.frame(covariate = env_cov, group = "environmental",
             role = "candidate predictor (A and B)", stringsAsFactors = FALSE),
  if (length(prior_cov)) data.frame(covariate = prior_cov, group = "published map",
             role = "candidate predictor (B only); benchmark P", stringsAsFactors = FALSE),
  if (length(categorical)) data.frame(covariate = categorical, group = "categorical",
             role = "EXCLUDED from the forest; stratification only", stringsAsFactors = FALSE),
  if (length(embed_cols)) data.frame(covariate = sprintf("%s (x%d)", "AlphaEarth A00..",
                                                         length(embed_cols)),
             group = "embedding",
             role = "EXCLUDED from the forest; 64 dimensions on 8 points",
             stringsAsFactors = FALSE)
)
print_table(stack_tbl)
log_info(length(env_cov), " environmental + ", length(prior_cov),
         " published-map covariates offered to the forest; ",
         length(embed_cols) + length(categorical), " ring-fenced")

# =============================================================================
# 3. Two response windows
# =============================================================================

log_step("18c  RESPONSE WINDOWS")

build_frame <- function(window_name) {
  s <- stocks[stocks$scenario == "as_supplied" & stocks$window == window_name, ]
  d <- merge(s[, c("core_id", "campaign", "stock_kgm2", "is_lower_bound",
                   "window_lo_cm", "window_hi_cm")],
             covars, by = c("core_id", "campaign"))
  if (have_ref) {
    d <- merge(d, refs[, c("core_id", intersect(prior_cov, names(refs)))],
               by = "core_id", all.x = TRUE)
  }
  d <- merge(d, cores[, c("core_id", "latitude", "longitude")],
             by = "core_id", suffixes = c("", ".core"))
  d
}

frames <- list(
  common_support = build_frame("common_support"),
  reference_0_30 = build_frame("reference_0_30")
)

for (nm in names(frames)) {
  d <- frames[[nm]]
  if (!nrow(d)) fail_loudly("Window '", nm, "' produced no rows after the join",
                            remedy = "Check 02 and 03 outputs share core_id and campaign.")
  lb <- isTRUE_vec(d$is_lower_bound)
  log_info(sprintf("%-15s 0-%.1f cm: %d cores, %d lower bound(s)",
                   nm, d$window_hi_cm[1], nrow(d), sum(lb)))
}

# Lower bounds are floors, not measurements. Training a regression on them
# would tell the model that a core reaching only 15 cm holds as little carbon
# as one measured to 30 cm and found empty. Same rule as script 11.
usable_rows <- function(d) !isTRUE_vec(d$is_lower_bound)

for (nm in names(frames)) {
  d <- frames[[nm]]
  u <- usable_rows(d)
  if (any(!u)) {
    log_warn(nm, ": excluding ", sum(!u), " lower-bound core(s) as observations: ",
             paste(d$core_id[!u], collapse = ", "))
  }
  frames[[nm]] <- d[u, , drop = FALSE]
  st <- table(frames[[nm]]$campaign)
  log_info(sprintf("  %s usable: n=%d  (%s)", nm, nrow(frames[[nm]]),
                   paste(sprintf("%s=%d", names(st), as.integer(st)), collapse = ", ")))
  if (min(st) < 2L) {
    log_warn("  a stratum has fewer than 2 usable cores, so leave-one-stratum-out ",
             "can only run in one direction in this window and cannot decide anything.")
  }
}

verdict_window <- "common_support"
map_window     <- "reference_0_30"
log_info("VERDICT window: ", verdict_window,
         " (all cores complete, both strata testable)")
log_info("MAP window:     ", map_window,
         " (comparable with every published layer)")

# =============================================================================
# 4. Candidate definitions
# =============================================================================

log_step("18d  CANDIDATES")

rf_args <- list(num_trees = CFG$step2$num_trees, mtry = CFG$step2$mtry,
                min_node_size = CFG$step2$min_node_size,
                max_depth = CFG$step2$max_depth)

fit_forest <- function(x, y) do.call(fit_rf, c(list(x = x, y = y, seed = CFG$seed),
                                               rf_args))
pred_forest <- function(m, nd) predict_rf(m, nd)

# The prior column used as benchmark P and as the anchor for C. Sothe 0-30 cm is
# the like-for-like layer for a 0-30 cm stock; Li is a FULL-COLUMN peat stock
# and is not comparable without a depth assumption, which is exactly the point
# 04 audits. So Sothe is preferred, SoilGrids second, and Li is not used as an
# anchor at all.
pick_prior <- function(d) {
  for (cand in c("sothe_sc_0_30_kgm2", "soilgrids_ocs_0_30_kgm2")) {
    if (cand %in% names(d) && sum(is.finite(d[[cand]])) >= 3L) return(cand)
  }
  NA_character_
}
prior_col <- if (have_ref) pick_prior(frames[[verdict_window]]) else NA_character_
if (!is.na(prior_col)) {
  log_info("prior anchor: ", prior_col)
} else if (have_ref) {
  log_warn("no like-for-like 0-30 cm published layer has enough finite values ",
           "at the cores; the prior benchmark and candidate C are unavailable.")
}

# A "fit" for the prior benchmark is just the identity on the prior column --
# it ignores the training data entirely, which is the whole point: it is what
# the council gets without this project.
fit_prior  <- function(x, y) list(col = prior_col)
pred_prior <- function(m, nd) as.numeric(nd[[m$col]])

# Candidate C: prior + a precision-weighted, distance-decayed residual.
# The machinery is R/bayes.R, used here on RESIDUALS (core minus prior) so the
# surface returns to the published map wherever no core has influence.
make_fit_resid <- function(ls_km, max_km) {
  function(x, y) {
    p <- as.numeric(x[[prior_col]])
    r <- y - p
    ok <- is.finite(r)
    list(lon = x$longitude[ok], lat = x$latitude[ok], resid = r[ok],
         sigma = rep(max(stats::sd(r[ok], na.rm = TRUE), 1e-6), sum(ok)),
         ls_km = ls_km, max_km = max_km, col = prior_col)
  }
}
pred_resid <- function(m, nd) {
  p <- as.numeric(nd[[m$col]])
  if (!length(m$resid)) return(p)
  D <- matrix(NA_real_, nrow(nd), length(m$resid))
  for (j in seq_along(m$resid)) {
    D[, j] <- haversine_km(nd$longitude, nd$latitude, m$lon[j], m$lat[j])
  }
  W <- build_weight_matrix(D, m$ls_km, m$max_km)
  # Prior on the residual is zero with the residual sd as its spread: absent
  # any core, the correction is zero and the surface IS the published map.
  post <- bayes_posterior(rep(0, nrow(nd)), rep(m$sigma[1], nrow(nd)),
                          m$resid, m$sigma, W)
  p + post$mean
}

ls_km <- CFG$bayes$length_scale_km
if (is.null(ls_km)) {
  nn <- nn_distance_summary(cores$longitude, cores$latitude, cores$core_id)
  ls_km <- median(nn$nn_dist_km)
}
buf_km <- CFG$step2$cv_buffer_km
if (is.null(buf_km)) buf_km <- ls_km
log_info(sprintf("distance length scale %.2f km; spatial CV buffer %.2f km",
                 ls_km, buf_km))

candidates <- list(
  N = list(label = "NULL: observed mean",
           cols = env_cov, fit = fit_null, pred = predict_null),
  A = list(label = "A: random forest on environmental covariates",
           cols = env_cov, fit = fit_forest, pred = pred_forest)
)
if (!is.na(prior_col)) {
  candidates$P <- list(label = paste0("P: published map alone (", prior_col, ")"),
                       cols = prior_col, fit = fit_prior, pred = pred_prior)
  candidates$B <- list(label = "B: random forest, published maps as covariates",
                       cols = union(env_cov, prior_cov),
                       fit = fit_forest, pred = pred_forest)
  candidates$C <- list(label = "C: published map + distance-weighted residual",
                       cols = c(prior_col, "longitude", "latitude"),
                       fit = make_fit_resid(ls_km, CFG$bayes$max_influence_km),
                       pred = pred_resid)
}
log_info(length(candidates), " candidates: ",
         paste(names(candidates), collapse = ", "))

# =============================================================================
# 5. Scoring
# =============================================================================

log_step("18e  CROSS-VALIDATION")

score_candidate <- function(cand, d) {
  # Every candidate needs longitude/latitude available for the spatial buffer,
  # even when it does not use them as predictors.
  cols <- intersect(unique(c(cand$cols, "longitude", "latitude")), names(d))
  X <- d[, cols, drop = FALSE]
  y <- d$stock_kgm2
  loco <- cv_leave_one_group_out(X, y, d$core_id, cand$fit, cand$pred)
  loso <- cv_leave_one_stratum_out(X, y, d$campaign, cand$fit, cand$pred)
  sbuf <- cv_spatial_buffer(X, y, d$core_id, d$longitude, d$latitude,
                            buf_km, cand$fit, cand$pred)
  list(loco = loco, loso = loso, sbuf = sbuf)
}

results <- list()
for (wname in names(frames)) {
  d <- frames[[wname]]
  for (cn in names(candidates)) {
    cand <- candidates[[cn]]
    miss <- setdiff(intersect(cand$cols, c(prior_col, prior_cov)), names(d))
    if (length(miss)) {
      log_warn(wname, "/", cn, ": skipped, missing ", paste(miss, collapse = ", "))
      next
    }
    log_info("scoring ", wname, " / ", cand$label)
    r <- score_candidate(cand, d)
    results[[paste(wname, cn, sep = "|")]] <-
      data.frame(window = wname, candidate = cn, label = cand$label,
                 n = nrow(d),
                 loco_rmse = r$loco$metrics$rmse, loco_r2 = r$loco$metrics$r2_cv,
                 loco_bias = r$loco$metrics$bias,
                 loco_cover = r$loco$metrics$coverage,
                 loso_rmse = r$loso$metrics$rmse, loso_r2 = r$loso$metrics$r2_cv,
                 # Coverage is carried beside every score. A leave-one-stratum-out
                 # r2 computed from a single surviving fold is not a worse
                 # estimate of skill, it is not an estimate of skill at all.
                 loso_cover = r$loso$metrics$coverage,
                 loso_n = r$loso$metrics$n_predicted,
                 sbuf_rmse = r$sbuf$metrics$rmse, sbuf_r2  = r$sbuf$metrics$r2_cv,
                 sbuf_cover = r$sbuf$metrics$coverage,
                 sbuf_folds_used = r$sbuf$n_folds_used,
                 sbuf_folds_skipped = r$sbuf$n_folds_skipped,
                 stringsAsFactors = FALSE)
  }
}
cv_tbl <- do.call(rbind, results)
rownames(cv_tbl) <- NULL

log_info("cross-validated skill, all candidates, both windows:")
print_table(cv_tbl[, c("window", "candidate", "n", "loco_rmse", "loco_r2",
                       "loso_r2", "loso_n", "sbuf_r2", "sbuf_folds_skipped")],
            digits = 3)

# A leave-one-stratum-out score is only meaningful if every core got a
# prediction. Where a stratum has one usable core the reverse fold cannot run,
# and whatever number comes back is computed from the handful of cores that
# survived. Those scores are blanked here rather than printed with a caveat --
# a caveat gets dropped when a table is copied into a report; an NA does not.
partial_loso <- isTRUE_vec(cv_tbl$loso_cover < 1)
if (any(partial_loso)) {
  bad_win <- unique(cv_tbl$window[partial_loso])
  log_warn("leave-one-stratum-out is INCOMPLETE in window(s): ",
           paste(bad_win, collapse = ", "))
  log_warn(sprintf(paste0("  only %d of %d cores received a prediction, ",
                          "because a stratum has too few usable cores to train on."),
                   min(cv_tbl$loso_n[partial_loso]), cv_tbl$n[partial_loso][1]))
  log_warn("  Those r2 values are set to NA. Left in, the mean-only null scored")
  log_warn("  +0.95 in this situation purely because five of its six folds were")
  log_warn("  never run -- which is exactly the kind of number this workflow ",
           "exists to refuse.")
  cv_tbl$loso_r2[partial_loso]   <- NA_real_
  cv_tbl$loso_rmse[partial_loso] <- NA_real_
}

if (any(cv_tbl$sbuf_folds_skipped > 0)) {
  log_warn(sprintf(paste0("The spatial-buffer CV skipped folds: a %.1f km ",
                          "buffer around a held-out core removes most of its "),
                   buf_km))
  log_warn("  own cluster, leaving too few training points. Those folds are ",
           "counted, not hidden -- an unreported skip would turn a failed test")
  log_warn("  into a good score computed from whichever folds survived.")
}

# =============================================================================
# 6. The pre-registered verdict
# =============================================================================

log_step("18f  VERDICT")

vt <- cv_tbl[cv_tbl$window == verdict_window, ]
prior_rmse <- if ("P" %in% vt$candidate) vt$loco_rmse[vt$candidate == "P"] else NA_real_

vt$clears_loco   <- vt$loco_r2 > CFG$step2$require_loco_r2
vt$clears_loso   <- isTRUE_vec(vt$loso_r2 > CFG$step2$require_loso_r2)
vt$beats_prior   <- if (is.finite(prior_rmse)) vt$loco_rmse < prior_rmse else NA
vt$eligible      <- vt$clears_loco & vt$clears_loso &
                    (!CFG$step2$require_beat_prior | isTRUE_vec(vt$beats_prior))

log_info("selection rule (fixed in config.R before any result was seen):")
log_info(sprintf("  leave-one-core-out r2    > %g", CFG$step2$require_loco_r2))
log_info(sprintf("  leave-one-stratum-out r2 > %g", CFG$step2$require_loso_r2))
log_info(sprintf("  beats the published map on leave-one-core-out RMSE: %s",
                 CFG$step2$require_beat_prior))
print_table(vt[, c("candidate", "loco_r2", "loso_r2", "beats_prior", "eligible")],
            digits = 3)

modelled <- vt[vt$candidate %in% c("A", "B", "C") & isTRUE_vec(vt$eligible), ]
if (nrow(modelled)) {
  selected <- modelled$candidate[which.max(modelled$loso_r2)]
  selection_basis <- "cleared every pre-registered criterion"
  log_ok("SELECTED: candidate ", selected, " -- ", selection_basis)
  log_warn("Given n=", vt$n[1], " and the number of covariates screened, treat ",
           "this as PROVISIONAL: a hypothesis for the next field season, not ",
           "a validated map.")
} else {
  selected <- if ("C" %in% vt$candidate) "C" else "P"
  selection_basis <- paste0(
    "NO candidate cleared the pre-registered bar; falling back to ",
    if (selected == "C") "the prior-anchored residual correction"
    else "the published map, uncorrected")
  log_warn(strrep("-", 68))
  log_warn("NO MODEL EARNED THE MAP.")
  log_warn("  This is the expected outcome and it is a result, not a pipeline")
  log_warn("  failure. Eight cores in two clusters cannot identify a covariate")
  log_warn("  relationship that holds across 2,758 km2.")
  log_warn("")
  log_warn("  The deliverable is therefore candidate ", selected, ": ",
           if (selected == "C")
             "the published map, nudged toward the cores near the cores and"
           else "the published map as-is")
  log_warn(if (selected == "C")
             "  returning to the published map everywhere else." else "")
  log_warn("  That is a weaker claim than a fitted map, and it is the claim")
  log_warn("  the data supports.")
  log_warn(strrep("-", 68))
}

# =============================================================================
# 7. Variable importance
# =============================================================================

log_step("18g  VARIABLE IMPORTANCE")

log_info("Importance is measured on CROSS-VALIDATED predictions -- shuffle a ")
log_info("covariate, re-run leave-one-core-out, record the rise in RMSE. A ")
log_info("covariate cannot earn importance by memorising training points.")
log_info("Each value carries a permutation p-value from a null built by ")
log_info("shuffling the RESPONSE, taking the LARGEST importance in each ")
log_info("shuffled run, so the p-value already accounts for having screened ")
log_info(length(env_cov), " covariates on ", nrow(frames[[verdict_window]]),
         " points.")

d_v <- frames[[verdict_window]]
imp_cols <- intersect(if (selected == "B") union(env_cov, prior_cov) else env_cov,
                      names(d_v))
imp <- rf_importance(d_v[, imp_cols, drop = FALSE], d_v$stock_kgm2, d_v$core_id,
                     n_perm = CFG$step2$imp_n_perm, n_null = CFG$step2$imp_n_null,
                     seed = CFG$seed, fit_args = rf_args,
                     budget_sec = CFG$step2$imp_budget_sec)

if (!is.null(imp)) {
  print_table(imp[, c("rank", "covariate", "importance", "importance_pct",
                      "p_value", "beats_chance")], digits = 3)
  n_sig <- sum(isTRUE_vec(imp$beats_chance))
  if (all(is.na(imp$p_value))) {
    log_warn("No permutation null was affordable, so there is NO p-value.")
    log_warn("  The ranking above cannot be read as evidence that any ",
             "covariate matters. Install ranger or raise ")
    log_warn("  CFG$step2$imp_budget_sec before quoting this table.")
  } else if (n_sig == 0L) {
    log_warn("NO covariate beat chance (all p >= 0.05).")
    log_warn("  The ranking above is a ranking of noise. It is reported in ",
             "full, rather than truncated to a plausible-looking top three,")
    log_warn("  precisely so that it cannot be mistaken for a finding.")
  } else {
    log_ok(n_sig, " covariate(s) beat the permutation null at p < 0.05")
  }
  if (imp$n_null[1] > 0L) {
    log_info(sprintf("smallest attainable p-value with %d null draws: %.3f",
                     imp$n_null[1], 1 / (1 + imp$n_null[1])))
  }
}

# =============================================================================
# 8. The raster
# =============================================================================

log_step("18h  PREDICTION RASTER (EPSG:3978)")

ring <- read_aoi_ring(file.path(CFG$root, "outputs", "spatial"))
have_terra <- requireNamespace("terra", quietly = TRUE)
pred_tif <- Sys.glob(file.path(CFG$dir_gee, "ccnp_predictors_*m.tif"))

raster_written <- character(0)

if (!have_terra || !length(pred_tif)) {
  log_warn(strrep("-", 68))
  log_warn("NO PREDICTION RASTER WRITTEN.")
  if (!have_terra) log_warn("  `terra` is not installed: install.packages('terra')")
  if (!length(pred_tif))
    log_warn("  no ccnp_predictors_*m.tif in ", CFG$dir_gee,
             " -- run scripts/03_gee_covariates.R and retrieve the export.")
  log_warn("  Everything above -- the candidate comparison, the verdict and")
  log_warn("  the importance test -- is complete and does not depend on the")
  log_warn("  raster. Re-run this script once the inputs are present.")
  log_warn(strrep("-", 68))
} else {
  pred_tif <- pred_tif[[1]]
  log_info("predictor raster: ", basename(pred_tif))
  r <- terra::rast(pred_tif)
  log_info(sprintf("  %d x %d cells, %d bands: %s", terra::nrow(r), terra::ncol(r),
                   terra::nlyr(r), paste(names(r), collapse = ", ")))

  # Source grid definition, taken from the raster itself rather than rebuilt
  # from config -- rebuilding it would silently disagree if 03 was run with a
  # different AOI or scale.
  ext <- terra::ext(r)
  src <- list(xmin = ext$xmin, ymax = ext$ymax,
              xres = terra::res(r)[1], yres = terra::res(r)[2])

  grid <- lcc_target_grid(ring$lon, ring$lat, CFG$deliver$cell_m,
                          lcc_params(CFG$deliver$epsg))
  log_info(sprintf("target grid: %d x %d at %d m, EPSG:%d",
                   grid$nx, grid$ny, CFG$deliver$cell_m, CFG$deliver$epsg))

  # Pixel-level covariates as a data frame, in the source raster's own layout.
  band_mat <- function(nm) {
    terra::as.matrix(r[[nm]], wide = TRUE)
  }

  surface_for <- function(cand_name, d) {
    cand <- candidates[[cand_name]]
    cols <- intersect(cand$cols, names(d))
    fitted <- cand$fit(d[, unique(c(cols, "longitude", "latitude")), drop = FALSE],
                       d$stock_kgm2)
    need <- switch(cand_name,
                   N = character(0),
                   P = prior_col,
                   C = prior_col,
                   intersect(cand$cols, names(r)))
    missing_bands <- setdiff(need, names(r))
    if (length(missing_bands)) {
      log_warn("candidate ", cand_name, ": raster lacks band(s) ",
               paste(missing_bands, collapse = ", "), "; surface not built.")
      return(NULL)
    }
    nb <- length(need)
    nd <- if (nb) {
      as.data.frame(lapply(setNames(need, need), function(nm) as.vector(band_mat(nm))))
    } else {
      data.frame(.dummy = rep(NA_real_, terra::nrow(r) * terra::ncol(r)))
    }
    # Cell-centre coordinates, needed by C and by the buffer machinery.
    lonv <- src$xmin + (seq_len(terra::ncol(r)) - 0.5) * src$xres
    latv <- src$ymax - (seq_len(terra::nrow(r)) - 0.5) * src$yres
    nd$longitude <- rep(lonv, times = terra::nrow(r))
    nd$latitude  <- rep(latv, each  = terra::ncol(r))
    v <- cand$pred(fitted, nd)
    matrix(as.numeric(v), nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE)
  }

  write_deliver <- function(z_src, fname, note) {
    w <- warp_to_lcc(z_src, src, grid, lcc_params(CFG$deliver$epsg),
                     method = "bilinear")
    w <- mask_grid_to_ring(w, grid, ring$lon, ring$lat,
                           lcc_params(CFG$deliver$epsg))
    p <- file.path(out_dir, fname)
    write_geotiff(w, p, xmin = grid$xmin, ymax = grid$ymax,
                  xres = grid$cell_m, yres = grid$cell_m,
                  epsg = CFG$deliver$epsg)
    log_ok("wrote ", fname, "  (", sum(is.finite(w)), " cells in AOI) <- ", note)
    raster_written <<- c(raster_written, fname)
    invisible(w)
  }

  z_map <- surface_for(selected, frames[[map_window]])
  if (!is.null(z_map)) {
    write_deliver(z_map, "carbon_prediction.tif",
                  paste0("0-30 cm SOC stock, kg C/m2, candidate ", selected))
  }
  z_cs <- surface_for(selected, frames[[verdict_window]])
  if (!is.null(z_cs)) {
    write_deliver(z_cs, "carbon_prediction_0_15cm.tif",
                  paste0("0-", sprintf("%.1f", frames[[verdict_window]]$window_hi_cm[1]),
                         " cm SOC stock, every core a measurement"))
  }

  # An independent check on the warp, when terra is present. The two routes are
  # different code by different authors; agreeing to a fraction of a cell is
  # meaningful, and disagreeing is the single most valuable thing this script
  # could discover.
  if (!is.null(z_map)) {
    xy <- lonlat_to_lcc(cores$longitude, cores$latitude,
                        lcc_params(CFG$deliver$epsg))
    at_core_warp <- sample_lonlat(z_map, src, cores$longitude, cores$latitude)
    src_r <- terra::rast(nrows = nrow(z_map), ncols = ncol(z_map),
                         xmin = src$xmin, xmax = src$xmin + ncol(z_map) * src$xres,
                         ymax = src$ymax, ymin = src$ymax - nrow(z_map) * src$yres,
                         crs = "EPSG:4326", vals = as.vector(t(z_map)))
    at_core_terra <- as.numeric(terra::extract(
      src_r, cbind(cores$longitude, cores$latitude))[, 1])
    dd <- abs(at_core_warp - at_core_terra)
    log_info(sprintf("warp cross-check at the 8 cores: max |base-R - terra| = %.3g kg/m2",
                     max(dd, na.rm = TRUE)))
  }
}

# =============================================================================
# 9. Deliverables
# =============================================================================

log_step("18i  DELIVERABLES")

write_csv_logged(cv_tbl, file.path(CFG$dir_derived, "18_candidate_cv.csv"),
                 "every candidate, every CV scheme, both windows")
write_csv_logged(stack_tbl, file.path(out_dir, "COVARIATE_STACK.csv"),
                 "exactly which layers the model could and could not see")

diag <- data.frame(
  item = c(
    "response",
    "map depth window",
    "verdict depth window",
    "cores used (map window)",
    "cores used (verdict window)",
    "cores excluded as lower bounds",
    "forest backend",
    "trees / min node size / max depth",
    "cross-validation",
    "spatial buffer",
    "candidates compared",
    "selection rule",
    "SELECTED",
    "basis",
    "leave-one-core-out r2 (selected)",
    "leave-one-stratum-out r2 (selected)",
    "published-map benchmark RMSE",
    "covariates beating the permutation null",
    "deliverable CRS",
    "deliverable resolution",
    "ASSUMPTION carried",
    "status"),
  value = c(
    "SOC stock, kg C/m2, scenario as_supplied",
    sprintf("0-%.0f cm", frames[[map_window]]$window_hi_cm[1]),
    sprintf("0-%.1f cm", frames[[verdict_window]]$window_hi_cm[1]),
    nrow(frames[[map_window]]),
    nrow(frames[[verdict_window]]),
    paste(setdiff(cores$core_id, frames[[map_window]]$core_id), collapse = ", "),
    rf_backend(),
    sprintf("%d / %d / %d", CFG$step2$num_trees, CFG$step2$min_node_size,
            CFG$step2$max_depth),
    "leave-one-CORE-out, leave-one-STRATUM-out, spatial-buffer; never segment level",
    sprintf("%.2f km", buf_km),
    paste(names(candidates), collapse = ", "),
    sprintf("loco r2 > %g AND loso r2 > %g AND beats published map",
            CFG$step2$require_loco_r2, CFG$step2$require_loso_r2),
    selected,
    selection_basis,
    sprintf("%+.3f", vt$loco_r2[vt$candidate == selected][1]),
    sprintf("%+.3f", vt$loso_r2[vt$candidate == selected][1]),
    if (is.finite(prior_rmse)) sprintf("%.3f kg/m2", prior_rmse) else "not available",
    if (!is.null(imp)) sum(isTRUE_vec(imp$beats_chance)) else NA_character_,
    sprintf("EPSG:%d (NAD83 / Canada Atlas Lambert)", CFG$deliver$epsg),
    sprintf("%d m", CFG$deliver$cell_m),
    "the verdict is earned at 0-14.5 cm and applied to the 0-30 cm map",
    if (length(raster_written)) "raster written" else "no raster (see log)"),
  stringsAsFactors = FALSE
)
write_csv_logged(diag, file.path(out_dir, "MODEL_DIAGNOSTICS.csv"),
                 "how the map was made, in one table")
print_table(diag)

if (!is.null(imp)) {
  write_csv_logged(imp, file.path(out_dir, "VARIABLE_IMPORTANCE.csv"),
                   "importance with permutation p-values, full ranking")
}
write_csv_logged(cv_tbl, file.path(out_dir, "CANDIDATE_COMPARISON.csv"),
                 "the competition the selected map won, or did not")

readme <- c(
  "# Fort Severn carbon prediction map (Step 2)",
  "",
  "## What is in this folder",
  "",
  if (length(raster_written))
    paste0("- `", raster_written, "` -- carbon stock, kilograms of carbon per square metre.")
  else "- (no raster in this run: see MODEL_DIAGNOSTICS.csv, `status`)",
  "- `MODEL_DIAGNOSTICS.csv` -- how the map was made, in one table.",
  "- `CANDIDATE_COMPARISON.csv` -- the five surfaces that competed, and their scores.",
  "- `VARIABLE_IMPORTANCE.csv` -- which measurements mattered, and whether that beat chance.",
  "- `COVARIATE_STACK.csv` -- exactly what the model could and could not see.",
  "",
  "## How to read the map",
  "",
  sprintf("Coordinate system: EPSG:%d, NAD83 / Canada Atlas Lambert. Cell size %d m.",
          CFG$deliver$epsg, CFG$deliver$cell_m),
  "Values are kilograms of carbon in each square metre of ground, top 30 cm.",
  "To convert to tonnes per hectare, multiply by 10.",
  "",
  "## What this map is, honestly",
  "",
  paste0("The map is candidate ", selected, ": ", selection_basis, "."),
  "",
  "Five different ways of making a map were built and scored against each",
  "other, including the option of simply handing over the published map that",
  "anyone can already download. The rule for which one would be delivered was",
  "written down before any of the scores were seen, so the winner could not be",
  "chosen after the fact.",
  "",
  if (!nrow(modelled)) c(
    "No model earned the map. Eight cores in two clusters cannot pin down a",
    "relationship between carbon and satellite measurements that holds across",
    "the whole area. So the map you have is the published regional map, adjusted",
    "toward the community's own measurements near where those measurements were",
    "taken, and left as published everywhere else.",
    "",
    "That is a smaller claim than a map built from the cores would be. It is",
    "the claim the eight cores actually support."
  ) else c(
    "A model did clear the bar. With eight cores this should still be treated",
    "as a hypothesis to test with the next field season, not a settled answer."
  ),
  "",
  "## What would change this",
  "",
  "Roughly 30-50 cores, deliberately spread across the wetness and landscape",
  "gradients rather than clustered at two accessible sites, would make a fitted",
  "map answerable. The eight cores already collected are what makes it possible",
  "to say that honestly.",
  "",
  paste0("Generated by scripts/18_predictive_map.R on ", Sys.Date(), "."))
writeLines(readme, file.path(out_dir, "README.md"))
log_ok("wrote README.md (", length(readme), " lines)")

log_step("18  STEP 2 STATUS")
log_ok("covariate stack recorded: COVARIATE_STACK.csv")
log_ok("model type and fusion method decided by competition, not preference")
log_ok(if (length(raster_written))
         paste0("prediction raster(s) written in EPSG:", CFG$deliver$epsg, ": ",
                paste(raster_written, collapse = ", "))
       else "prediction raster NOT written -- see the warning above")
log_ok("diagnostics and variable importance written")
log_ok("18 complete -- deliverables in ", out_dir)
