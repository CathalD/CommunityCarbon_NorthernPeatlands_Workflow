# =============================================================================
# R/depth.R  --  Depth semantics, segment geometry, and like-for-like
#                integration of carbon stock over a depth window.
#
# The single most dangerous assumption in this dataset is what the `Depth`
# column means. Everything here is written so that the meaning is DERIVED from
# the data and checked, never assumed.
#
# Unit conventions used throughout:
#   thickness, top, bottom    cm
#   bulk density              g/cm3
#   SOC                       % by dry mass
#   carbon density            g/cm3   (= bulk density * SOC/100)
#   areal stock               g/cm2   (= carbon density integrated over cm)
#   reported stock            kg/m2   (= g/cm2 * 10)
# =============================================================================

# ---- unit conversions (pure) -------------------------------------------------

gcm2_to_kgm2 <- function(x) x * 10
kgm2_to_tha  <- function(x) x * 10
tha_to_kgm2  <- function(x) x / 10

#' Carbon density from bulk density and carbon percentage. Pure.
#'
#' This is the first-principles definition:
#'   rho_C [g/cm3] = BD [g/cm3] * SOC [%] / 100
#'
#' Note the units: this is a CONCENTRATION (mass of carbon per unit volume of
#' soil), not an areal density. It only becomes an areal stock (g/cm2) after
#' multiplication by a layer thickness in cm.
carbon_density_gcm3 <- function(bulk_density_gcm3, soc_pct) {
  bulk_density_gcm3 * soc_pct / 100
}

# ---- depth column semantics --------------------------------------------------

#' Decide, from evidence, whether `Depth` holds segment THICKNESS or a
#' cumulative BOTTOM depth.
#'
#' The discriminating test: if the column were a cumulative bottom depth it
#' would have to increase monotonically down every core. A single decrease
#' anywhere falsifies that reading.
#'
#' Returns a list with the verdict and the evidence behind it, so the caller
#' can log it rather than take it on faith. Pure.
infer_depth_semantics <- function(core_id, seg_index, depth) {
  ord <- order(core_id, seg_index)
  ci  <- core_id[ord]; dv <- depth[ord]

  violations <- data.frame(core_id = character(0), from = numeric(0),
                           to = numeric(0), kind = character(0),
                           stringsAsFactors = FALSE)
  for (cc in unique(ci)) {
    d <- dv[ci == cc]
    if (length(d) < 2L) next
    dd <- diff(d)
    for (k in which(dd <= 0)) {
      violations <- rbind(violations, data.frame(
        core_id = cc, from = d[k], to = d[k + 1L],
        # A strict decrease is impossible for a cumulative bottom depth.
        # A tie is equally impossible, but by implying a zero-thickness
        # segment rather than an inverted one. Report them distinctly.
        kind = if (dd[k] < 0) "decrease" else "tie",
        stringsAsFactors = FALSE))
    }
  }

  n_dec  <- sum(violations$kind == "decrease")
  n_tie  <- sum(violations$kind == "tie")
  n_multi <- sum(tapply(seq_along(ci), ci, length) > 1L)
  verdict <- if (nrow(violations) > 0L) "thickness_cm" else "ambiguous"

  reasoning <- if (n_dec > 0L) {
    paste0("`Depth` strictly DECREASES at ", n_dec, " step(s) within ",
           length(unique(violations$core_id[violations$kind == "decrease"])),
           " core(s) when read in segment order. A cumulative bottom depth ",
           "cannot decrease, so the column is segment thickness",
           if (n_tie > 0L) paste0(" (a further ", n_tie, " step(s) repeat the ",
             "same value, which a cumulative depth could only do via a ",
             "zero-thickness segment)") else "", ".")
  } else if (n_tie > 0L) {
    paste0("`Depth` repeats the same value at ", n_tie, " step(s) within a ",
           "core. Read as a cumulative bottom depth this would require ",
           "zero-thickness segments, so the column is segment thickness.")
  } else {
    paste0("`Depth` never decreases or repeats within a core, so thickness ",
           "and cumulative-bottom readings are not distinguishable from ",
           "monotonicity alone. Corroborate against total core lengths.")
  }

  list(
    verdict              = verdict,
    n_cores_multi_seg    = n_multi,
    n_decreases          = n_dec,
    n_ties               = n_tie,
    monotonic_violations = violations,
    reasoning            = reasoning
  )
}

#' Attach top/bottom depths to segments by accumulating thickness down each
#' core, in segment order. Pure: returns a new data frame.
#'
#' Requires columns: core_id, seg_index, thickness_cm.
add_depth_bounds <- function(df) {
  stopifnot(all(c("core_id", "seg_index", "thickness_cm") %in% names(df)))
  df <- df[order(df$core_id, df$seg_index), , drop = FALSE]
  sp <- split(seq_len(nrow(df)), df$core_id)
  df$depth_bottom_cm <- NA_real_
  df$depth_top_cm    <- NA_real_
  for (idx in sp) {
    b <- cumsum(df$thickness_cm[idx])
    df$depth_bottom_cm[idx] <- b
    df$depth_top_cm[idx]    <- b - df$thickness_cm[idx]
  }
  rownames(df) <- NULL
  df
}

# ---- integration over a depth window ----------------------------------------

#' Integrate carbon density over [z_lo, z_hi] for one core.
#'
#' Segments that straddle a window boundary are cut, and only the overlapping
#' fraction contributes. Cutting a segment assumes carbon density is uniform
#' WITHIN that segment -- the only assumption available without finer
#' resolution, and stated explicitly rather than buried.
#'
#' The function reports how much of the requested window the core actually
#' covers. If coverage is short, the returned stock is a LOWER BOUND on the
#' true stock over the window, and is flagged as such. It is never silently
#' rescaled to the full window.
#'
#' Pure. Returns a one-row data frame.
integrate_window <- function(depth_top_cm, depth_bottom_cm, carbon_density_gcm3,
                             z_lo, z_hi, tol = 1e-8) {
  stopifnot(z_hi > z_lo)
  ov <- pmax(0, pmin(depth_bottom_cm, z_hi) - pmax(depth_top_cm, z_lo))
  keep <- ov > 0

  stock_gcm2  <- sum(carbon_density_gcm3[keep] * ov[keep])
  covered_cm  <- sum(ov[keep])
  requested   <- z_hi - z_lo

  data.frame(
    window_lo_cm     = z_lo,
    window_hi_cm     = z_hi,
    requested_cm     = requested,
    covered_cm       = covered_cm,
    coverage_frac    = covered_cm / requested,
    n_segments_used  = sum(keep),
    stock_gcm2       = stock_gcm2,
    stock_kgm2       = gcm2_to_kgm2(stock_gcm2),
    is_lower_bound   = covered_cm < requested - tol,
    stringsAsFactors = FALSE
  )
}

#' Deepest window, starting at the surface, that EVERY core fully covers.
#'
#' This is the only depth interval over which all cores are directly
#' comparable with no extrapolation whatsoever. Pure.
common_support_depth <- function(core_id, depth_bottom_cm) {
  tot <- tapply(depth_bottom_cm, core_id, max, na.rm = TRUE)
  list(
    depth_cm       = min(tot),
    limiting_core  = names(tot)[which.min(tot)],
    core_totals_cm = tot[order(tot)]
  )
}

#' EXTRAPOLATION (opt-in, never a default).
#'
#' Extend a short core to the full window by carrying the carbon density of its
#' deepest observed segment downward. Returns the extrapolated stock alongside
#' the observed part so the two are never conflated.
#'
#' In a peat complex this is a weak assumption in a KNOWN direction: carbon
#' density in these profiles changes sharply at the peat/mineral contact, and a
#' core that stopped early may have stopped ON that contact. Treat the result
#' as a scenario, not an estimate. Pure.
extrapolate_window <- function(depth_top_cm, depth_bottom_cm,
                               carbon_density_gcm3, z_lo, z_hi) {
  obs <- integrate_window(depth_top_cm, depth_bottom_cm,
                          carbon_density_gcm3, z_lo, z_hi)
  gap_cm <- max(0, obs$requested_cm - obs$covered_cm)
  deepest <- carbon_density_gcm3[which.max(depth_bottom_cm)]
  add_gcm2 <- deepest * gap_cm

  data.frame(
    gap_cm                  = gap_cm,
    density_carried_gcm3    = if (gap_cm > 0) deepest else NA_real_,
    stock_extrap_gcm2       = obs$stock_gcm2 + add_gcm2,
    stock_extrap_kgm2       = gcm2_to_kgm2(obs$stock_gcm2 + add_gcm2),
    extrap_share_of_stock   = if (obs$stock_gcm2 + add_gcm2 > 0) {
      add_gcm2 / (obs$stock_gcm2 + add_gcm2)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Apply integrate_window() (and optionally extrapolate_window()) to every core.
#' Pure with respect to its inputs; returns one row per core.
integrate_all_cores <- function(seg, z_lo, z_hi, extrapolate = FALSE) {
  stopifnot(all(c("core_id", "depth_top_cm", "depth_bottom_cm",
                  "carbon_density_gcm3") %in% names(seg)))
  ids <- unique(seg$core_id)
  out <- do.call(rbind, lapply(ids, function(cc) {
    s <- seg[seg$core_id == cc, , drop = FALSE]
    r <- integrate_window(s$depth_top_cm, s$depth_bottom_cm,
                          s$carbon_density_gcm3, z_lo, z_hi)
    r <- cbind(data.frame(core_id = cc, stringsAsFactors = FALSE), r)
    if (extrapolate) {
      r <- cbind(r, extrapolate_window(s$depth_top_cm, s$depth_bottom_cm,
                                       s$carbon_density_gcm3, z_lo, z_hi))
    }
    r
  }))
  rownames(out) <- NULL
  out
}
