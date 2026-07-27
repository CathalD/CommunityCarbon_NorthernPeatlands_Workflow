# =============================================================================
# R/qc.R  --  Quality control.
#
# Two tiers, deliberately separated:
#
#   TIER 1  HARD GATES  -> stop the run.
#           Reserved for conditions that would silently corrupt the science
#           downstream and that no later step could detect: wrong hemisphere,
#           missing columns, non-reconciling units.
#
#   TIER 2  FLAGS       -> annotate and carry forward.
#           Every implausible value is surfaced with a named flag and the row
#           is KEPT. Dropping data is an analyst's decision, not the
#           pipeline's, and a flagged row still carries information.
#
# All predicates are pure and vectorised.
# =============================================================================

# ---- flag registry -----------------------------------------------------------

#' Human-readable definition of every flag this workflow can raise.
#' Written out with the results so a reader never has to guess what a flag means.
qc_flag_dictionary <- function() {
  data.frame(
    flag = c(
      "LON_SIGN_REPAIRED", "COORD_OUT_OF_REGION", "COORD_MISSING",
      "BD_ABOVE_MAX", "BD_BELOW_MIN", "BD_PEAT_INCONSISTENT",
      "OM_OUT_OF_RANGE", "SOC_EXCEEDS_OM", "SOC_ABOVE_CEILING",
      "OCD_RECONCILE_FAIL", "THICKNESS_NONPOSITIVE",
      "SEGMENT_INDEX_ABSENT", "DUPLICATE_LAB_VALUES", "SINGLE_SEGMENT_CORE"
    ),
    meaning = c(
      "Longitude sign was wrong for the stated study area and was repaired. The supplied value would place this core in the eastern hemisphere.",
      "Coordinate falls outside the plausible region around Fort Severn even after any sign repair.",
      "Latitude and/or longitude is missing or non-numeric.",
      "Dry bulk density exceeds the plausible ceiling for soil; approaches parent-rock density. Suspect a sample-volume or compaction error.",
      "Dry bulk density below the plausible floor for any soil or moss layer.",
      "Sample reads as peat by organic-matter content but carries a bulk density far too high for peat. The two measurements disagree.",
      "Organic matter percentage outside 0-100.",
      "Organic carbon percentage exceeds organic matter percentage, which is physically impossible: carbon is a component of organic matter.",
      "Organic carbon percentage above the ceiling for soil organic carbon.",
      "The supplied Organic Carbon Density column does not reconcile with a first-principles recomputation from bulk density and SOC.",
      "Segment thickness is zero, negative or missing.",
      "Sample Id carries no trailing integer, so segment order had to be assigned rather than read.",
      "This sample shares an identical organic matter AND organic carbon pair with a sample from a DIFFERENT core. Possible transcription or copy-paste error.",
      "Core is represented by a single segment, so within-core variability cannot be estimated for it."
    ),
    tier = c("repair", rep("flag", 13)),
    stringsAsFactors = FALSE
  )
}

# ---- Tier 1: hard gates ------------------------------------------------------

#' Verify the supplied file has exactly the columns the workflow expects.
gate_columns <- function(df, expected) {
  missing <- setdiff(expected, names(df))
  extra   <- setdiff(names(df), expected)
  if (length(missing)) {
    fail_loudly(
      "Raw data is missing expected columns",
      paste0("Missing: ", paste(missing, collapse = ", ")),
      paste0("Present: ", paste(names(df), collapse = ", ")),
      remedy = "Update CFG$raw_cols in config.R if the supplier renamed columns, and re-check every downstream reference."
    )
  }
  if (length(extra)) {
    log_warn("Unexpected extra column(s) present, ignored: ",
             paste(extra, collapse = ", "))
  }
  invisible(TRUE)
}

#' Geolocation gate.
#'
#' The supplied file stores western longitudes as positive numbers. At 56 N,
#' longitude +87.7 is in western Siberia; the stated study area is Fort Severn,
#' Ontario at -87.7. Every covariate extraction, every reference-layer lookup
#' and every map would be drawn from the wrong continent, and nothing
#' downstream would notice, because the values would all be perfectly valid
#' numbers for the wrong place.
#'
#' So this is a hard gate. It stops unless the analyst has explicitly consented
#' to the repair in config.R, and when it does repair it says so loudly and
#' records the repair per row.
#'
#' Returns list(longitude, flags, repaired, diagnosis).
gate_geolocation <- function(lon, lat, anchor, bbox, allow_fix) {
  n <- length(lon)
  flags    <- replicate(n, character(0), simplify = FALSE)
  repaired <- rep(FALSE, n)

  bad <- !is.finite(lon) | !is.finite(lat)
  if (any(bad)) {
    for (i in which(bad)) flags[[i]] <- c(flags[[i]], "COORD_MISSING")
  }

  in_box <- function(x, y) {
    is.finite(x) & is.finite(y) &
      x >= bbox[["xmin"]] & x <= bbox[["xmax"]] &
      y >= bbox[["ymin"]] & y <= bbox[["ymax"]]
  }

  ok_asis    <- in_box(lon, lat)
  ok_flipped <- in_box(-lon, lat)

  # Rows that are wrong as supplied but correct with the sign flipped.
  needs_flip <- !ok_asis & ok_flipped & !bad

  if (any(needs_flip)) {
    k <- sum(needs_flip)
    # Measure the displacement rather than estimate it: it is the single most
    # persuasive number in this message.
    disp <- haversine_km(lon[needs_flip], lat[needs_flip],
                         -lon[needs_flip], lat[needs_flip])
    diagnosis <- c(
      sprintf("%d of %d rows carry a longitude that is positive where the", k, n),
      sprintf("study area (%s) requires a negative value.",
              "Fort Severn, Ontario"),
      "",
      sprintf("  supplied longitude range : %+.5f to %+.5f",
              min(lon[needs_flip]), max(lon[needs_flip])),
      sprintf("  expected for study area  : %+.5f to %+.5f",
              bbox[["xmin"]], bbox[["xmax"]]),
      sprintf("  anchor (Fort Severn)     : %+.5f, %+.5f",
              anchor[["lon"]], anchor[["lat"]]),
      "",
      "As supplied, these coordinates fall in the EASTERN hemisphere",
      sprintf("(western Siberia), %.0f km from where they belong.", mean(disp)),
      "Every covariate and reference-layer extraction would silently sample",
      "the wrong continent."
    )

    if (!isTRUE(allow_fix)) {
      fail_loudly(
        "Longitude sign is inconsistent with the stated study area",
        diagnosis,
        remedy = paste0(
          "Confirm with the data provider that these are western longitudes, ",
          "then set CFG$allow_longitude_sign_fix <- TRUE in config.R to apply ",
          "the repair with a per-row provenance flag."
        )
      )
    }

    log_warn(strrep("-", 68))
    log_warn("REPAIRING LONGITUDE SIGN on ", k, " of ", n, " rows")
    for (d in diagnosis) log_warn(d)
    log_warn("Repair is authorised by CFG$allow_longitude_sign_fix = TRUE.")
    log_warn("Every affected row is flagged LON_SIGN_REPAIRED.")
    log_warn(strrep("-", 68))

    lon[needs_flip] <- -lon[needs_flip]
    repaired[needs_flip] <- TRUE
    for (i in which(needs_flip)) flags[[i]] <- c(flags[[i]], "LON_SIGN_REPAIRED")
  }

  # Anything still outside the region after repair is flagged, not repaired.
  still_bad <- !in_box(lon, lat) & !bad
  if (any(still_bad)) {
    for (i in which(still_bad)) {
      flags[[i]] <- c(flags[[i]], "COORD_OUT_OF_REGION")
    }
  }

  list(longitude = lon, flags = flags, repaired = repaired,
       n_repaired = sum(repaired), n_out_of_region = sum(still_bad))
}

#' Reconciliation gate for the supplied derived column.
#'
#' Recomputes carbon density from first principles and compares. A systematic
#' failure here means the supplied column means something other than its label,
#' and the run stops rather than propagating an unknown quantity.
gate_ocd_reconciliation <- function(supplied, recomputed, tol,
                                    label = "Organic Carbon Density (g/cm^2)") {
  ok <- near_rel(supplied, recomputed, tol)
  ok[is.na(supplied) & is.na(recomputed)] <- TRUE
  n_bad <- sum(!ok, na.rm = TRUE)

  if (n_bad > 0 && n_bad >= 0.5 * length(ok)) {
    worst <- order(abs(supplied - recomputed), decreasing = TRUE)[1:min(5, n_bad)]
    fail_loudly(
      "Supplied derived column does not reconcile with first principles",
      paste0("Column: ", label),
      sprintf("%d of %d rows disagree beyond a relative tolerance of %g",
              n_bad, length(ok), tol),
      "",
      "Worst disagreements (supplied vs recomputed = BD * SOC/100):",
      sprintf("  row %3d: supplied %.8g, recomputed %.8g, ratio %.6g",
              worst, supplied[worst], recomputed[worst],
              supplied[worst] / recomputed[worst]),
      remedy = "Establish what the supplied column actually measures before using any of it."
    )
  }
  list(ok = ok, n_bad = n_bad)
}

# ---- Tier 2: value flags -----------------------------------------------------

#' Apply every value-plausibility predicate. Returns a list of character
#' vectors, one per row, to be merged with any flags raised earlier. Pure.
flag_values <- function(df, qc) {
  n <- nrow(df)
  f <- replicate(n, character(0), simplify = FALSE)
  add <- function(f, idx, tag) {
    for (i in which(idx)) f[[i]] <- c(f[[i]], tag)
    f
  }
  na_false <- function(x) { x[is.na(x)] <- FALSE; x }

  f <- add(f, na_false(df$bulk_density_gcm3 > qc$bd_max),  "BD_ABOVE_MAX")
  f <- add(f, na_false(df$bulk_density_gcm3 < qc$bd_min),  "BD_BELOW_MIN")
  f <- add(f, na_false(df$om_pct >= qc$om_peat_threshold &
                       df$bulk_density_gcm3 > qc$bd_peat_max),
           "BD_PEAT_INCONSISTENT")
  f <- add(f, na_false(df$om_pct < qc$om_min | df$om_pct > qc$om_max),
           "OM_OUT_OF_RANGE")
  f <- add(f, na_false(df$soc_pct > df$om_pct), "SOC_EXCEEDS_OM")
  f <- add(f, na_false(df$soc_pct > qc$soc_max_pct), "SOC_ABOVE_CEILING")
  f <- add(f, na_false(df$thickness_cm < qc$thickness_min_cm |
                       is.na(df$thickness_cm)), "THICKNESS_NONPOSITIVE")
  f
}

#' Flag identical (OM, SOC) pairs shared across DIFFERENT cores.
#'
#' Two segments from different cores landing on byte-identical lab values is
#' more consistent with a transcription error than with nature. Pure.
flag_cross_core_duplicates <- function(core_id, om_pct, soc_pct) {
  n <- length(core_id)
  key <- paste(om_pct, soc_pct, sep = "|")
  f <- replicate(n, character(0), simplify = FALSE)
  for (k in unique(key[!is.na(om_pct) & !is.na(soc_pct)])) {
    idx <- which(key == k)
    if (length(idx) > 1L && length(unique(core_id[idx])) > 1L) {
      for (i in idx) f[[i]] <- c(f[[i]], "DUPLICATE_LAB_VALUES")
    }
  }
  f
}

#' Merge several per-row flag lists into one, de-duplicated. Pure.
merge_flags <- function(...) {
  ls <- list(...)
  n <- length(ls[[1]])
  lapply(seq_len(n), function(i) unique(unlist(lapply(ls, `[[`, i))))
}

#' Collapse a per-row flag list to a delimited character column. Pure.
flags_to_string <- function(flags, sep = ";") {
  vapply(flags, function(x) if (!length(x)) "" else paste(sort(x), collapse = sep),
         character(1))
}

#' Tally flags for the run log. Pure.
summarise_flags <- function(flags) {
  all <- unlist(flags)
  if (!length(all)) {
    return(data.frame(flag = character(0), n_rows = integer(0),
                      stringsAsFactors = FALSE))
  }
  t <- sort(table(all), decreasing = TRUE)
  data.frame(flag = names(t), n_rows = as.integer(t),
             stringsAsFactors = FALSE)
}

# ---- OM -> carbon conversion forensics ---------------------------------------

#' Identify the organic-matter-to-carbon factor a data provider used.
#'
#' The supplied SOC column is a derived quantity in most community datasets:
#' organic matter measured by loss-on-ignition, multiplied by a conversion
#' factor. Which factor was used is rarely recorded, and it matters, because
#' the classic van Bemmelen factor of 0.58 is known to overestimate carbon in
#' true peat.
#'
#' Rather than assume, this estimates the factor per group by regression
#' through the origin, which is robust to the heavy rounding in the supplied
#' columns. Pure.
identify_om_carbon_factor <- function(om_pct, soc_pct, group) {
  ok <- is.finite(om_pct) & is.finite(soc_pct) & om_pct > 0
  do.call(rbind, lapply(unique(group[ok]), function(g) {
    i <- ok & group == g
    om <- om_pct[i]; soc <- soc_pct[i]
    # Least squares through the origin: minimises sum (soc - b*om)^2
    b  <- sum(om * soc) / sum(om * om)
    r  <- soc - b * om
    data.frame(
      group            = g,
      n                = sum(i),
      factor_est       = b,
      implied_divisor  = 1 / b,
      ratio_min        = min(soc / om),
      ratio_max        = max(soc / om),
      resid_max_abs    = max(abs(r)),
      om_range         = sprintf("%.2f-%.2f", min(om), max(om)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Name a numeric factor against known conventions. Pure.
name_om_carbon_factor <- function(b, tol = 0.01) {
  known <- c("van Bemmelen (0.58 = 1/1.724)" = 0.58,
             "1/1.8 (0.556)"                  = 1 / 1.8,
             "0.50 (half of OM)"              = 0.50,
             "0.51 (northern peat, typical)"  = 0.51)
  d <- abs(known - b)
  if (min(d) <= tol) names(known)[which.min(d)] else sprintf("unrecognised (%.4f)", b)
}
