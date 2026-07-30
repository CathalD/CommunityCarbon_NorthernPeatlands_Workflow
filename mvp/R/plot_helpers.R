# =============================================================================
# plot_helpers.R  --  shared plotting utilities for steps 10, 11 and 12.
#
# Base R graphics throughout, deliberately: these figures need no ggplot2, and
# keeping the dependency list short is worth more here than terse plotting code.
# =============================================================================

# ---- palettes ---------------------------------------------------------------

# Organic warm, mineral cool, our cores black so they always read as the
# subject rather than as one series among thousands.
COL <- c(organic = "#C1660A", mineral = "#1F5FA8", unknown = "#8A8A8A",
        ours = "#111111")

DS_COL <- c("CanPeat"      = "#C1660A",
           "NPDB"          = "#1F5FA8",
           "Janousek"      = "#2E8B6F",
           "WOSIS Canada"  = "#8B5FA8",
           "Fort Severn"   = "#111111")

#' Same colour at a given opacity. Defined here rather than relying on base R's
#' colour-adjust helper, so there is nothing to get wrong about its name or
#' argument signature across R versions.
fade <- function(col, alpha) {
  v <- grDevices::col2rgb(col) / 255
  grDevices::rgb(v[1], v[2], v[3], alpha = alpha)
}

# ---- figure device ----------------------------------------------------------

#' Open a PNG, run the plotting code, close the device even if it errors.
fig <- function(name, expr, width = 1600, height = 1150, res = 185) {
  png(file.path(CFG$dir_figures, paste0(name, ".png")),
      width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  msg("wrote figures/", name, ".png")
}

# ---- depth-profile machinery ------------------------------------------------

GRID_CM <- 0:300

#' Step-interpolate one profile's per-cm values onto the common 1 cm grid.
#' Returns NA below the profile's own bottom -- a profile is never extended
#' past what was actually measured.
densify <- function(upper, lower, value, grid = GRID_CM) {
  out <- rep(NA_real_, length(grid))
  ok <- is.finite(upper) & is.finite(lower) & is.finite(value) & lower > upper
  for (i in which(ok)) out[grid >= upper[i] & grid < lower[i]] <- value[i]
  out
}

#' Per-group percentile envelope across profiles, on the common grid.
#'
#' Uses split() rather than filtering inside a loop: NPDB alone is 9,017
#' profiles over 48,372 layers, and a per-profile subset scan would be
#' quadratic.
#'
#' @param value_col column holding the per-cm quantity to summarise
envelope <- function(layer_df, group_col, value_col = "dens_kgm2_per_cm") {
  groups <- unique(layer_df[[group_col]])
  groups <- groups[!is.na(groups)]
  lapply(setNames(groups, groups), function(g) {
    sub <- layer_df[!is.na(layer_df[[group_col]]) & layer_df[[group_col]] == g, ]
    by_prof <- split(sub, sub$profile_id)
    if (length(by_prof) < 3) return(NULL)
    m <- vapply(by_prof, function(s) {
      densify(s$upper_depth, s$lower_depth, s[[value_col]])
    }, numeric(length(GRID_CM)))
    qs <- function(p) apply(m, 1, quantile, probs = p, na.rm = TRUE, names = FALSE)
    data.frame(depth = GRID_CM,
              n = apply(m, 1, function(v) sum(is.finite(v))),
              p25 = qs(0.25), p50 = qs(0.50), p75 = qs(0.75),
              p10 = qs(0.10), p90 = qs(0.90))
  })
}

#' Draw one envelope. Only where at least `min_n` profiles still contribute, so
#' the band does not narrow to a spurious line at depth as profiles drop out.
draw_env <- function(e, col, min_n = 3, lo = "p10", hi = "p90", lwd = 2.6) {
  if (is.null(e)) return(invisible())
  e <- e[is.finite(e$p50) & e$n >= min_n, ]
  if (!nrow(e)) return(invisible())
  polygon(c(e[[lo]], rev(e[[hi]])), c(e$depth, rev(e$depth)),
          col = fade(col, 0.18), border = NA)
  lines(e$p50, e$depth, col = col, lwd = lwd)
}

#' Overlay the eight Fort Severn cores as individual stepped lines.
draw_our_cores <- function(oseg, col = COL[["ours"]], lwd = 2.2,
                          value_col = "carbon_density_kgm2_per_cm") {
  for (cc in unique(oseg$core_id)) {
    s <- oseg[oseg$core_id == cc, ]
    s <- s[order(s$depth_top_cm), ]
    for (i in seq_len(nrow(s))) {
      lines(rep(s[[value_col]][i], 2),
            c(s$depth_top_cm[i], s$depth_bottom_cm[i]), col = col, lwd = lwd)
    }
  }
}

#' Median cumulative carbon down the profile, per group, plus the share of the
#' full column sitting above `ref` cm.
#'
#' The share is computed PER PROFILE and then medianed. A ratio of medians is
#' not the median of ratios and would quietly answer a different question.
cumulative_by_group <- function(layer_df, group_col, ref = 30) {
  groups <- unique(layer_df[[group_col]])
  groups <- groups[!is.na(groups)]
  lapply(setNames(groups, groups), function(g) {
    sub <- layer_df[!is.na(layer_df[[group_col]]) & layer_df[[group_col]] == g, ]
    by_prof <- split(sub, sub$profile_id)
    if (length(by_prof) < 3) return(NULL)
    m <- vapply(by_prof, function(s) {
      d <- densify(s$upper_depth, s$lower_depth, s$dens_kgm2_per_cm)
      d[is.na(d)] <- 0
      cs <- cumsum(d)
      # Blank below the profile's own bottom so a shallow core does not read as
      # a deep one that stopped accumulating carbon.
      bottom <- suppressWarnings(max(s$lower_depth, na.rm = TRUE))
      cs[GRID_CM > bottom] <- NA_real_
      cs
    }, numeric(length(GRID_CM)))
    share <- vapply(seq_len(ncol(m)), function(j) {
      v <- m[, j]
      tot <- suppressWarnings(max(v, na.rm = TRUE))
      at_ref <- v[GRID_CM == ref]
      if (!is.finite(tot) || tot <= 0 || !length(at_ref) || !is.finite(at_ref)) return(NA_real_)
      at_ref / tot
    }, numeric(1))
    list(
      curve = data.frame(depth = GRID_CM,
                        n = apply(m, 1, function(v) sum(is.finite(v))),
                        p50 = apply(m, 1, quantile, 0.5, na.rm = TRUE, names = FALSE)),
      share_at_ref = median(share, na.rm = TRUE),
      n_profiles = ncol(m)
    )
  })
}
