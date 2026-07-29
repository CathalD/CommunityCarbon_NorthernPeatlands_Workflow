# =============================================================================
# R/rf.R  --  Random forest for Step 2, with a backend that can be checked.
#
# WHY THIS IS SHAPED THIS WAY
#
# The plan asks for `ranger`, and ranger is used when it is installed. But two
# things argue against calling ranger directly from the script:
#
#   1. VARIABLE IMPORTANCE MUST MEAN ONE THING. ranger's "impurity" importance
#      is a sum of split criteria; its "permutation" importance is an OOB error
#      increase; neither is comparable to what another backend would report,
#      and impurity importance is known to favour covariates with many distinct
#      values -- which, on eight points, is every continuous covariate. So
#      importance here is computed OUTSIDE the backend, by permuting a covariate
#      and re-scoring the CROSS-VALIDATED predictions. That definition is the
#      same whichever backend fitted the trees, and it is the definition that
#      matches the question the council will ask: how much worse is the map if
#      we did not know this?
#
#   2. IMPORTANCE ON EIGHT POINTS IS MOSTLY NOISE. A bare ranking would present
#      noise as a finding. Every importance value here is therefore accompanied
#      by a permutation p-value from a null distribution built by shuffling the
#      RESPONSE, which is the only honest way to say whether a covariate beat
#      chance in a sample this small.
#
# The fallback backend is a real bagged regression forest (bootstrap resampling,
# mtry candidate sampling at each split, variance-reduction splits). It exists
# so the cross-validation harness, the importance test and the whole of script
# 18 can be exercised on a machine with no CRAN access, and so ranger's output
# has something independent to be compared against.
#
# EVERY function here is seeded explicitly. A forest is stochastic and an
# unseeded one would give the council a different map on every run.
# =============================================================================

#' Which forest backend is available. Pure with respect to the search path.
rf_backend <- function() {
  if (requireNamespace("ranger", quietly = TRUE)) "ranger" else "base_r_bagged"
}

# ---- base-R backend ----------------------------------------------------------

#' Grow one regression tree on a bootstrap sample, sampling mtry covariates at
#' each split. Internal; returns a nested list of nodes. Pure given `seed`.
#'
#' Splits minimise the within-node sum of squares, the standard CART criterion
#' for regression. Candidate split points are the midpoints between consecutive
#' distinct values, so a split never sits exactly on an observation.
.rf_grow_tree <- function(X, y, mtry, min_node, max_depth) {
  grow <- function(idx, depth) {
    yy <- y[idx]
    if (length(idx) <= min_node || depth >= max_depth ||
        length(unique(yy)) == 1L) {
      return(list(leaf = TRUE, value = mean(yy)))
    }
    p <- ncol(X)
    vars <- sample.int(p, size = min(mtry, p))
    best <- list(sse = Inf)
    for (v in vars) {
      xv <- X[idx, v]
      u  <- sort(unique(xv[is.finite(xv)]))
      if (length(u) < 2L) next
      cuts <- (u[-length(u)] + u[-1L]) / 2
      for (cut in cuts) {
        l <- xv <= cut
        if (sum(l, na.rm = TRUE) < 1L || sum(!l, na.rm = TRUE) < 1L) next
        yl <- yy[which(l)]; yr <- yy[which(!l)]
        sse <- sum((yl - mean(yl))^2) + sum((yr - mean(yr))^2)
        if (is.finite(sse) && sse < best$sse) {
          best <- list(sse = sse, var = v, cut = cut)
        }
      }
    }
    if (!is.finite(best$sse)) return(list(leaf = TRUE, value = mean(yy)))
    l <- X[idx, best$var] <= best$cut
    l[is.na(l)] <- TRUE                     # NA follows the left branch
    list(leaf = FALSE, var = best$var, cut = best$cut,
         left  = grow(idx[which(l)],  depth + 1L),
         right = grow(idx[which(!l)], depth + 1L))
  }
  grow(seq_len(nrow(X)), 0L)
}

#' Predict from one tree. Internal, vectorised over rows. Pure.
.rf_predict_tree <- function(node, X) {
  out <- numeric(nrow(X))
  walk <- function(nd, idx) {
    if (length(idx) == 0L) return(invisible(NULL))
    if (isTRUE(nd$leaf)) { out[idx] <<- nd$value; return(invisible(NULL)) }
    v <- X[idx, nd$var]
    l <- v <= nd$cut
    l[is.na(l)] <- TRUE
    walk(nd$left,  idx[which(l)])
    walk(nd$right, idx[which(!l)])
  }
  walk(node, seq_len(nrow(X)))
  out
}

# ---- public interface --------------------------------------------------------

#' Fit a random forest. Uses ranger when available, the base-R forest otherwise.
#'
#' @param x data frame of numeric covariates.
#' @param y numeric response.
#' @param num_trees,mtry,min_node_size standard forest settings. `mtry` defaults
#'   to ceiling(p/3), the regression convention.
#' @param seed integer; the forest is stochastic and must be reproducible.
#' @return an object carrying its backend, so predict_rf() cannot be handed the
#'   wrong kind of model silently.
#'
#' `min_node_size` defaults to 2 rather than ranger's regression default of 5.
#' With eight observations a minimum node of 5 would allow at most one split in
#' the whole tree, which is not a forest, it is a threshold. This is a real
#' consequence of n=8 and is reported by 18 rather than buried here.
fit_rf <- function(x, y, num_trees = 500L, mtry = NULL, min_node_size = 2L,
                   max_depth = 8L, seed = 1L) {
  x <- as.data.frame(x)
  keep <- vapply(x, function(v) is.numeric(v) && sum(is.finite(v)) > 0L &&
                                stats::sd(v, na.rm = TRUE) > 0, logical(1))
  # A covariate that is constant within this (possibly held-out) fold carries no
  # information and would break ranger; dropping it here keeps CV folds running
  # rather than silently returning NA for the fold.
  x <- x[, keep, drop = FALSE]
  if (!ncol(x)) {
    return(structure(list(backend = "degenerate", value = mean(y, na.rm = TRUE)),
                     class = "ccnp_rf"))
  }
  if (is.null(mtry)) mtry <- max(1L, ceiling(ncol(x) / 3))
  mtry <- min(as.integer(mtry), ncol(x))

  if (rf_backend() == "ranger") {
    fit <- ranger::ranger(
      x = x, y = y, num.trees = num_trees, mtry = mtry,
      min.node.size = min_node_size, max.depth = max_depth,
      keep.inbag = TRUE, seed = seed, num.threads = 1L,
      respect.unordered.factors = "order"
    )
    return(structure(list(backend = "ranger", fit = fit, vars = names(x)),
                     class = "ccnp_rf"))
  }

  set.seed(seed)
  X <- as.matrix(x)
  n <- nrow(X)
  trees <- lapply(seq_len(num_trees), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    .rf_grow_tree(X[idx, , drop = FALSE], y[idx], mtry, min_node_size, max_depth)
  })
  structure(list(backend = "base_r_bagged", trees = trees, vars = colnames(X)),
            class = "ccnp_rf")
}

#' Predict from a fitted forest.
#'
#' @param per_tree if TRUE, also return the matrix of individual tree
#'   predictions. That matrix is what Step 3's uncertainty raster is built from,
#'   so it is produced by the same call that produces the map -- the two can
#'   never come from different fits.
predict_rf <- function(model, newx, per_tree = FALSE) {
  stopifnot(inherits(model, "ccnp_rf"))
  newx <- as.data.frame(newx)

  if (identical(model$backend, "degenerate")) {
    mu <- rep(model$value, nrow(newx))
    return(if (per_tree) list(pred = mu, trees = matrix(mu, ncol = 1L)) else mu)
  }

  missing_vars <- setdiff(model$vars, names(newx))
  if (length(missing_vars)) {
    stop("predict_rf(): covariate(s) absent from newx: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  newx <- newx[, model$vars, drop = FALSE]

  if (identical(model$backend, "ranger")) {
    if (!per_tree) {
      return(as.numeric(stats::predict(model$fit, data = newx)$predictions))
    }
    tp <- stats::predict(model$fit, data = newx, predict.all = TRUE)$predictions
    tp <- matrix(as.numeric(tp), nrow = nrow(newx))
    return(list(pred = rowMeans(tp, na.rm = TRUE), trees = tp))
  }

  X <- as.matrix(newx)
  # A lambda, not `vapply(..., X = X)`: the latter collides with vapply's own
  # `X` argument and passes the trees in as the data.
  tp <- vapply(model$trees, function(tr) .rf_predict_tree(tr, X),
               numeric(nrow(X)))
  tp <- matrix(tp, nrow = nrow(X))
  if (!per_tree) rowMeans(tp, na.rm = TRUE) else
    list(pred = rowMeans(tp, na.rm = TRUE), trees = tp)
}

#' Spread across trees at each prediction point: the model-uncertainty layer.
#'
#' This is the standard deviation of the ensemble members, NOT a predictive
#' interval. It captures disagreement between trees and says nothing about
#' whether the training data covered the ground being predicted -- that is what
#' the Area of Applicability is for. Reporting one as the other is the most
#' common way a machine-learning map understates its own uncertainty, so 18 and
#' 19 always report both.
rf_predict_sd <- function(model, newx) {
  pt <- predict_rf(model, newx, per_tree = TRUE)
  if (ncol(pt$trees) < 2L) return(rep(NA_real_, length(pt$pred)))
  apply(pt$trees, 1L, stats::sd, na.rm = TRUE)
}

# ---- importance, defined outside the backend ---------------------------------

#' Cross-validated permutation importance, with a permutation null.
#'
#' For each covariate: shuffle it, re-run the SAME leave-one-group-out
#' cross-validation, and record how much RMSE rose. Shuffling inside the CV
#' rather than after a single fit means importance is measured on predictions
#' for held-out cores, so a covariate cannot earn importance by memorising the
#' training points.
#'
#' The null distribution comes from shuffling the RESPONSE instead, `n_null`
#' times, and taking the largest importance observed in each shuffled run. Using
#' the MAXIMUM across covariates makes the resulting p-value account for the
#' fact that many covariates were screened -- otherwise, with 20 covariates on 8
#' points, something is nearly certain to look important by chance.
#'
#' @return data frame: covariate, importance (rise in RMSE), importance_pct,
#'   p_value, and `beats_chance`. Pure given `seed`.
#' @param budget_sec wall-clock budget. The permutation null is quadratic in the
#'   covariate count and linear in the draw count, so on the base-R backend a
#'   full-size run can take hours. Rather than either hanging or silently
#'   producing a weaker test, the draw count is TRIMMED to fit the budget and
#'   the trim is reported, along with the smallest p-value the reduced null can
#'   attain -- because that floor, not the nominal 0.05, is what the result must
#'   be read against. `Inf` disables trimming.
rf_importance <- function(x, y, group, n_perm = 10L, n_null = 50L,
                          seed = 1L, fit_args = list(), verbose = TRUE,
                          budget_sec = Inf) {
  x <- as.data.frame(x)
  vars <- names(x)[vapply(x, is.numeric, logical(1))]
  if (!length(vars)) return(NULL)

  cv_rmse <- function(dat, yy, s) {
    fitf <- function(a, b) do.call(fit_rf, c(list(x = a, y = b, seed = s), fit_args))
    cv <- cv_leave_one_group_out(dat, yy, group, fitf,
                                 function(m, nd) predict_rf(m, nd))
    rmse(yy, cv$predictions$pred)
  }

  set.seed(seed)
  t0 <- Sys.time()
  base_rmse <- cv_rmse(x, y, seed)
  per_run <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # This test is intentionally expensive -- a permutation null is the only
  # defensible way to rank covariates on eight points -- so its cost is
  # announced rather than discovered. A silent half-hour looks like a hang.
  runs_for <- function(np, nn) 1L + length(vars) * np + nn * (1L + length(vars))
  if (isTRUE(verbose)) {
    log_info(sprintf(
      "permutation importance: %d covariates x %d shuffles, plus a %d-draw null",
      length(vars), n_perm, n_null))
    log_info(sprintf(
      "  %d cross-validation runs on the '%s' backend, ~%.2f s each: about %.0f min",
      runs_for(n_perm, n_null), rf_backend(), per_run,
      runs_for(n_perm, n_null) * per_run / 60))
  }

  if (is.finite(budget_sec) && runs_for(n_perm, n_null) * per_run > budget_sec) {
    afford <- budget_sec / per_run
    n_perm_new <- max(3L, min(n_perm, as.integer(floor(
      (afford - 1L) / (2 * length(vars))))))
    n_null_new <- max(0L, as.integer(floor(
      (afford - 1L - length(vars) * n_perm_new) / (1L + length(vars)))))
    if (isTRUE(verbose)) {
      log_warn(sprintf("that exceeds the %.0f s budget for this step.", budget_sec))
      log_warn(sprintf("  trimming to %d shuffles and a %d-draw null.",
                       n_perm_new, n_null_new))
      if (n_null_new >= 1L) {
        log_warn(sprintf(paste0("  the smallest p-value a %d-draw null can ",
                                "produce is %.3f, so read the p-values against ",
                                "THAT floor, not against 0.05."),
                         n_null_new, 1 / (1 + n_null_new)))
      } else {
        log_warn("  the budget cannot afford a null at all: importance is ",
                 "reported WITHOUT a p-value, which makes the ranking")
        log_warn("  uninterpretable as evidence. Install ranger, or raise the ",
                 "budget, before using this table for anything.")
      }
    }
    n_perm <- n_perm_new; n_null <- n_null_new
  }

  # Observed importance: mean RMSE rise over n_perm shuffles of each covariate.
  obs <- vapply(vars, function(v) {
    rises <- vapply(seq_len(n_perm), function(k) {
      xp <- x
      xp[[v]] <- sample(xp[[v]])
      cv_rmse(xp, y, seed) - base_rmse
    }, numeric(1))
    mean(rises)
  }, numeric(1))

  # Null: shuffle the RESPONSE, recompute every covariate's importance, keep the
  # largest. This is the distribution of "best covariate found by chance".
  null_max <- vapply(seq_len(n_null), function(b) {
    yb <- sample(y)
    rb <- cv_rmse(x, yb, seed)
    max(vapply(vars, function(v) {
      xp <- x; xp[[v]] <- sample(xp[[v]])
      cv_rmse(xp, yb, seed) - rb
    }, numeric(1)))
  }, numeric(1))

  # +1 in numerator and denominator: the observed value is itself one draw from
  # the null under the null hypothesis, so a p-value of exactly 0 is not
  # attainable and should not be reported.
  #
  # With no null draws at all there is no p-value. NA is returned rather than 1
  # or 0: either number would be a claim about evidence that was never gathered.
  p <- if (!length(null_max)) rep(NA_real_, length(obs)) else
    vapply(obs, function(o) (1 + sum(null_max >= o)) / (1 + length(null_max)),
           numeric(1))

  out <- data.frame(
    covariate      = vars,
    importance     = as.numeric(obs),
    importance_pct = 100 * as.numeric(obs) / base_rmse,
    p_value        = as.numeric(p),
    beats_chance   = isTRUE_vec(as.numeric(p) < 0.05),
    base_rmse      = base_rmse,
    n_null         = length(null_max),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$importance), ]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

# ---- spatial cross-validation -------------------------------------------------

#' Leave-one-core-out with a SPATIAL EXCLUSION BUFFER.
#'
#' Plain leave-one-out is optimistic when points are clustered: a held-out core
#' sitting 200 m from a training core is predicted from what is effectively
#' itself. This variant removes, from each training fold, every core within
#' `buffer_km` of the held-out one. That is the standard spatial-CV correction
#' (Roberts et al. 2017, Ecography; Meyer et al. 2018, Env. Model. Softw.), and
#' with two tight clusters here it amounts to asking the harder question the
#' council's map actually poses: predict this place from somewhere else.
#'
#' Folds whose training set falls below `min_train` are SKIPPED and counted, not
#' silently dropped -- a buffer wide enough to empty the training set would
#' otherwise report an excellent score computed from two surviving folds.
#'
#' @return predictions, metrics, and the fold accounting. Pure given fit/predict.
cv_spatial_buffer <- function(data, y, group, lon, lat, buffer_km,
                              fit_fun, predict_fun, min_train = 3L) {
  stopifnot(length(y) == nrow(data), length(group) == nrow(data),
            length(lon) == nrow(data), length(lat) == nrow(data))
  d <- pairwise_km(lon, lat)
  g <- unique(group)
  pred <- rep(NA_real_, nrow(data))
  acc  <- vector("list", length(g))

  for (k in seq_along(g)) {
    te <- group == g[k]
    # Excluded: anything within the buffer of ANY held-out row.
    near <- apply(d[, which(te), drop = FALSE], 1L, min) <= buffer_km
    tr <- !te & !near
    acc[[k]] <- data.frame(group = g[k], n_train = sum(tr),
                           n_excluded_by_buffer = sum(near & !te),
                           used = sum(tr) >= min_train,
                           stringsAsFactors = FALSE)
    if (sum(tr) < min_train) next
    m <- try(fit_fun(data[tr, , drop = FALSE], y[tr]), silent = TRUE)
    if (inherits(m, "try-error")) next
    p <- try(predict_fun(m, data[te, , drop = FALSE]), silent = TRUE)
    if (inherits(p, "try-error")) next
    pred[te] <- as.numeric(p)
  }

  folds <- do.call(rbind, acc)
  list(
    predictions = data.frame(group = group, obs = y, pred = pred,
                             stringsAsFactors = FALSE),
    metrics     = metric_set(y, pred),
    folds       = folds,
    buffer_km   = buffer_km,
    n_folds_used = sum(folds$used),
    n_folds_skipped = sum(!folds$used)
  )
}
