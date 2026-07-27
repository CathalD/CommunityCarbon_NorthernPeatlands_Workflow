# =============================================================================
# R/stats.R  --  Estimation, resampling and validation.
#
# Everything here is built around one constraint: there are EIGHT cores, in
# two spatial clusters, and the two clusters are confounded with campaign year,
# landscape type and laboratory batch. That constraint is not something to work
# around; it is the main scientific finding about what this dataset can support.
#
# The functions therefore make the weakness visible rather than smoothing it
# over: cross-validation is always at the core level, the null model is always
# reported alongside, and the variogram diagnostic exists to demonstrate why
# kriging is refused rather than to perform it.
#
# All functions are pure. Every stochastic function takes an explicit seed.
# =============================================================================

# ---- point metrics (pure) ----------------------------------------------------

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
bias <- function(obs, pred) mean(pred - obs, na.rm = TRUE)

#' Nash-Sutcliffe / CV R-squared against the observed mean.
#' Deliberately allowed to go negative: a negative value is the honest signal
#' that the model is worse than predicting the mean, and it must not be hidden.
r2_cv <- function(obs, pred) {
  sse <- sum((obs - pred)^2, na.rm = TRUE)
  sst <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (sst == 0) return(NA_real_)
  1 - sse / sst
}

#' Skill relative to a mean-only null model. > 0 means the covariates helped.
skill_vs_null <- function(obs, pred) {
  1 - rmse(obs, pred) / rmse(obs, rep(mean(obs, na.rm = TRUE), length(obs)))
}

metric_set <- function(obs, pred) {
  data.frame(
    n           = sum(is.finite(obs) & is.finite(pred)),
    rmse        = rmse(obs, pred),
    mae         = mae(obs, pred),
    bias        = bias(obs, pred),
    r2_cv       = r2_cv(obs, pred),
    skill_null  = skill_vs_null(obs, pred),
    stringsAsFactors = FALSE
  )
}

# ---- interval estimation -----------------------------------------------------

#' Percentile bootstrap CI for the mean. Seeded.
#'
#' With n = 3 the bootstrap resamples from three numbers and can only ever
#' produce 10 distinct means; the interval it returns is a description of those
#' three numbers, not a reliable statement about the population. The returned
#' object carries `trustworthy` so callers cannot present it as more than it is.
boot_mean_ci <- function(x, B = 10000L, conf = 0.95, seed = 1L) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) {
    return(data.frame(n = n, mean = if (n == 1L) x else NA_real_,
                      lo = NA_real_, hi = NA_real_, method = "none",
                      trustworthy = FALSE, stringsAsFactors = FALSE))
  }
  set.seed(seed)
  bm <- replicate(B, mean(sample(x, n, replace = TRUE)))
  a <- (1 - conf) / 2
  q <- stats::quantile(bm, c(a, 1 - a), names = FALSE)
  data.frame(n = n, mean = mean(x), lo = q[1], hi = q[2],
             method = sprintf("percentile bootstrap, B=%d", B),
             trustworthy = n >= 8L, stringsAsFactors = FALSE)
}

#' Classical t interval, reported alongside the bootstrap so the reader can see
#' how much the interval depends on the method when n is this small.
t_mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) {
    return(data.frame(n = n, mean = if (n == 1L) x else NA_real_,
                      lo = NA_real_, hi = NA_real_, method = "none",
                      trustworthy = FALSE, stringsAsFactors = FALSE))
  }
  se <- stats::sd(x) / sqrt(n)
  tc <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
  data.frame(n = n, mean = mean(x), lo = mean(x) - tc * se,
             hi = mean(x) + tc * se,
             method = sprintf("t interval, df=%d", n - 1),
             trustworthy = n >= 8L, stringsAsFactors = FALSE)
}

#' Design-based estimate for each stratum, plus a study-area mean formed as a
#' stratum-weighted combination.
#'
#' Weights must come from the AREA each stratum occupies, not from how many
#' cores happen to sit in it. Sample sizes here reflect field logistics, not
#' the landscape, so using them as weights would silently re-weight the map.
#' When areal weights are unavailable the function returns the per-stratum
#' estimates and declines to produce a study-area number.
stratified_estimate <- function(value, stratum, area_weight = NULL,
                                B = 10000L, conf = 0.95, seed = 1L) {
  st <- sort(unique(stratum))
  per <- do.call(rbind, lapply(seq_along(st), function(i) {
    v <- value[stratum == st[i] & is.finite(value)]
    b <- boot_mean_ci(v, B = B, conf = conf, seed = seed + i)
    t <- t_mean_ci(v, conf = conf)
    data.frame(stratum = st[i], n = b$n, mean = b$mean,
               sd = if (b$n > 1) stats::sd(v) else NA_real_,
               se = se_mean(v), cv_pct = cv_pct(v),
               boot_lo = b$lo, boot_hi = b$hi,
               t_lo = t$lo, t_hi = t$hi,
               trustworthy = b$trustworthy,
               stringsAsFactors = FALSE)
  }))

  overall <- NULL
  if (!is.null(area_weight)) {
    w <- area_weight[match(per$stratum, names(area_weight))]
    if (all(is.finite(w)) && sum(w) > 0) {
      w <- w / sum(w)
      mu <- sum(w * per$mean)
      # Stratified variance: sum w_h^2 * SE_h^2
      v  <- sum(w^2 * per$se^2, na.rm = TRUE)
      overall <- data.frame(
        stratum = "STUDY AREA (area-weighted)", n = sum(per$n), mean = mu,
        sd = NA_real_, se = sqrt(v), cv_pct = NA_real_,
        boot_lo = mu - 1.96 * sqrt(v), boot_hi = mu + 1.96 * sqrt(v),
        t_lo = NA_real_, t_hi = NA_real_, trustworthy = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  rbind(per, overall)
}

# ---- cross-validation --------------------------------------------------------

#' Leave-one-out cross-validation at a chosen grouping unit.
#'
#' `group` MUST be the core identifier. Splitting at the segment level would
#' place segments from the same core on both sides of the split; they share a
#' location, a campaign, an operator and a laboratory batch, so the model would
#' be tested on data it had effectively already seen. That inflates apparent
#' skill dramatically and is the single easiest way to produce a fraudulent
#' validation statistic from a dataset this size.
#'
#' Returns per-fold predictions plus the metric set. Pure given fit/predict.
cv_leave_one_group_out <- function(data, y, group, fit_fun, predict_fun) {
  stopifnot(length(y) == nrow(data), length(group) == nrow(data))
  g <- unique(group)
  pred <- rep(NA_real_, nrow(data))
  for (gg in g) {
    te <- group == gg
    tr <- !te
    if (sum(tr) < 2L) next
    m <- try(fit_fun(data[tr, , drop = FALSE], y[tr]), silent = TRUE)
    if (inherits(m, "try-error")) next
    p <- try(predict_fun(m, data[te, , drop = FALSE]), silent = TRUE)
    if (inherits(p, "try-error")) next
    pred[te] <- as.numeric(p)
  }
  list(
    predictions = data.frame(group = group, obs = y, pred = pred,
                             stringsAsFactors = FALSE),
    metrics     = metric_set(y, pred),
    n_folds     = length(g)
  )
}

#' Leave-one-STRATUM-out: train on one stratum, predict the other.
#'
#' This asks a harder and more honest question than leave-one-core-out, namely
#' whether the model has learned anything transferable at all, or has only
#' learned the difference between two clusters of points. For a covariate model
#' intended to map unsampled ground, this is the relevant test, and it is
#' expected to fail here.
cv_leave_one_stratum_out <- function(data, y, stratum, fit_fun, predict_fun) {
  st <- unique(stratum)
  pred <- rep(NA_real_, nrow(data))
  for (s in st) {
    te <- stratum == s
    tr <- !te
    if (sum(tr) < 2L) next
    m <- try(fit_fun(data[tr, , drop = FALSE], y[tr]), silent = TRUE)
    if (inherits(m, "try-error")) next
    p <- try(predict_fun(m, data[te, , drop = FALSE]), silent = TRUE)
    if (inherits(p, "try-error")) next
    pred[te] <- as.numeric(p)
  }
  list(
    predictions = data.frame(stratum = stratum, obs = y, pred = pred,
                             stringsAsFactors = FALSE),
    metrics     = metric_set(y, pred),
    n_folds     = length(st)
  )
}

# ---- spatial structure diagnostics -------------------------------------------

#' Empirical variogram cloud for a small point set.
#'
#' This exists as EVIDENCE, not as a modelling step. With eight points there
#' are only 28 pairs; the function reports how those pairs distribute across
#' lag distance so a reader can see for themselves whether a variogram model
#' could be fitted. Pure.
variogram_cloud <- function(lon, lat, value) {
  n <- length(value)
  d <- pairwise_km(lon, lat)
  iu <- which(upper.tri(d), arr.ind = TRUE)
  data.frame(
    i        = iu[, 1],
    j        = iu[, 2],
    dist_km  = d[iu],
    semivar  = 0.5 * (value[iu[, 1]] - value[iu[, 2]])^2,
    stringsAsFactors = FALSE
  )
}

#' Bin a variogram cloud into lag classes and report pair counts per bin.
#' The pair COUNT is the point of this function: geostatistical practice asks
#' for roughly 30+ pairs per lag before a bin is considered informative.
variogram_binned <- function(cloud, n_bins = 5L) {
  b <- cut(cloud$dist_km, breaks = n_bins)
  do.call(rbind, lapply(levels(b), function(lv) {
    s <- cloud[b == lv & !is.na(b), , drop = FALSE]
    data.frame(lag = lv, n_pairs = nrow(s),
               dist_mean_km = if (nrow(s)) mean(s$dist_km) else NA_real_,
               semivar_mean = if (nrow(s)) mean(s$semivar) else NA_real_,
               informative  = nrow(s) >= 30L,
               stringsAsFactors = FALSE)
  }))
}

#' Nearest-neighbour distance summary: quantifies spatial clustering, which
#' governs how far any interpolation could legitimately reach. Pure.
nn_distance_summary <- function(lon, lat, id) {
  d <- pairwise_km(lon, lat)
  diag(d) <- Inf
  nn <- apply(d, 1, min)
  data.frame(id = id, nn_dist_km = nn,
             nn_of = id[apply(d, 1, which.min)],
             stringsAsFactors = FALSE)
}

# ---- simple, defensible models ----------------------------------------------

#' Fit an ordinary least squares model on a small number of covariates.
#'
#' Deliberately not a random forest or a gradient-boosted ensemble. With eight
#' observations, a flexible learner has enough capacity to reproduce the
#' training data exactly while learning nothing generalisable, and its variable
#' importances would be noise presented with unearned authority. A linear model
#' with one or two terms is something whose failure modes can be read off the
#' coefficients.
fit_ols <- function(x, y) {
  df <- as.data.frame(x)
  df$..y <- y
  stats::lm(..y ~ ., data = df)
}

predict_ols <- function(model, newx) {
  stats::predict(model, newdata = as.data.frame(newx))
}

#' Mean-only null model. The benchmark every covariate model must beat before
#' it earns the right to be mapped.
fit_null <- function(x, y) mean(y, na.rm = TRUE)
predict_null <- function(model, newx) rep(model, nrow(as.data.frame(newx)))

#' Exhaustive single-covariate screening under leave-one-core-out CV.
#'
#' With n = 8 the honest search space is one covariate at a time. Screening
#' many candidates and reporting only the winner is itself a form of
#' overfitting, so the function returns the FULL ranking, including how many
#' candidates were tried, so that selection optimism stays visible.
screen_single_covariates <- function(covars, y, group) {
  nm <- names(covars)
  res <- do.call(rbind, lapply(nm, function(v) {
    x <- covars[, v, drop = FALSE]
    if (!is.numeric(x[[1]]) || stats::sd(x[[1]], na.rm = TRUE) == 0) return(NULL)
    cv <- cv_leave_one_group_out(x, y, group, fit_ols, predict_ols)
    cbind(data.frame(covariate = v, stringsAsFactors = FALSE), cv$metrics)
  }))
  if (is.null(res)) return(res)
  res <- res[order(-res$skill_null), , drop = FALSE]
  res$n_candidates_screened <- length(nm)
  res$rank <- seq_len(nrow(res))
  rownames(res) <- NULL
  res
}
