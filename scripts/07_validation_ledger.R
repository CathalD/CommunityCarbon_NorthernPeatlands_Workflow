# =============================================================================
# 07_validation_ledger.R
#
# Two jobs.
#
# 1. THE LEDGER. For every product this workflow emits, state what its
#    validation actually is, and whether that validation is genuinely external
#    or contaminated. Written unconditionally, because the honest accounting of
#    what is NOT validated is the part most easily left out.
#
# 2. THE BENCHMARK. Where 04 has run, compare the published layers against the
#    community cores at the eight core locations. Done only for layers whose
#    depth support the audit found comparable.
#
# A NOTE ON WHAT "EXTERNAL" MEANS HERE
#
#   Comparing SoilGrids or Sothe et al. against these eight cores is often
#   called external validation. For these particular cores it is closer to
#   QUASI-EXTERNAL, and the distinction is worth keeping:
#
#     - These cores were collected in 2024 and 2025 by a community programme.
#       SoilGrids 2.0 and Sothe et al. 2022 predate them, so they cannot
#       contain these observations. That part is genuinely external.
#
#     - But both products were trained on national and global soil profile
#       databases that include Canadian peatland profiles, potentially from
#       the Hudson Bay Lowlands. The products have therefore seen the
#       POPULATION these cores were drawn from, even though they have not seen
#       these cores. Agreement is weaker evidence than it looks, and
#       disagreement is stronger.
#
#     - Eight points cannot validate a continental product in any case. What
#       follows is a local consistency check: it can reveal a gross local bias,
#       and cannot certify accuracy.
#
# INPUT   data/derived/02_core_stocks.csv, 05_stratum_estimates.csv
#         data/derived/04_reference_at_cores.csv  (optional; needs GEE)
#         data/derived/06_model_cv.csv            (optional; needs GEE)
# OUTPUT  outputs/tables/07_validation_ledger.csv
#         outputs/tables/07_reference_benchmark.csv  (when 04 has run)
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

log_step("07  VALIDATION LEDGER")

stocks <- require_artifact(file.path(CFG$dir_derived, "02_core_stocks.csv"),
                           "scripts/02_depth_harmonise.R")

# =============================================================================
# 1. The ledger
# =============================================================================

log_step("07a  WHAT IS AND IS NOT VALIDATED")

ledger <- data.frame(
  product = c(
    "PRODUCT 1: stratum mean stock",
    "PRODUCT 2: stratum-mean assignment map",
    "PRODUCT 3: covariate prediction",
    "REFERENCE: SoilGrids 2.0 OCS 0-30",
    "REFERENCE: Sothe et al. 2022",
    "REFERENCE: Li et al. 2025 HBL peat carbon",
    "Core stocks, common-support window",
    "Core stocks, 0-30 cm window"
  ),
  scientific_status = c(
    "SAMPLE-BASED ESTIMATE (design-based)",
    "CLASS-MEAN ASSIGNMENT (design-based, projected onto a categorical layer)",
    "PREDICTION / EXTRAPOLATION",
    "REFERENCE (external product)",
    "REFERENCE (external product)",
    "REFERENCE (external product, DIFFERENT DEPTH SUPPORT)",
    "OBSERVATION (interpolation-ready)",
    "OBSERVATION, two cores are LOWER BOUNDS"
  ),
  validation_performed = c(
    "None possible. n=3 and n=5 per stratum; the reported intervals describe the sampled cores, not the landscape.",
    "None. The class-mean assumption is untested: no independent cores exist within either stratum to hold out.",
    "Leave-one-core-out and leave-one-stratum-out cross-validation, scored against a mean-only null.",
    "Local consistency check against 8 cores at 0-30 cm.",
    "Local consistency check against 8 cores at 0-30 cm.",
    "NONE. Barred from comparison: its depth support is not 0-30 cm.",
    "Internal arithmetic only: every derived value recomputed from first principles and reconciled.",
    "As above, plus explicit lower-bound flags."
  ),
  is_external = c(
    "n/a", "n/a", "NO -- same 8 cores fit and tested",
    "QUASI-EXTERNAL", "QUASI-EXTERNAL", "n/a",
    "n/a", "n/a"
  ),
  leakage_risk = c(
    "None: no model is fitted. The risk here is not leakage but unwarranted precision.",
    "None from fitting. The stratum boundaries were NOT learned from the cores, so the map is not circular; but the two strata coincide exactly with the two campaigns.",
    "CONTROLLED. Splits are at core level, never segment level. Residual optimism remains from screening many candidates on 8 points; the full ranking is published so this stays visible.",
    "LOW but non-zero: the training profile databases behind SoilGrids may include Hudson Bay Lowlands profiles from the same population as these cores.",
    "LOW but non-zero, for the same reason. A Canada-specific product is MORE likely than a global one to have used regional profiles.",
    "n/a -- no comparison is made.",
    "None.",
    "None."
  ),
  honest_summary = c(
    "Two group means from eight cores. Report with the cores visible.",
    "A two-colour map. Defensible as a first approximation; not a measurement of spatial pattern.",
    "Expected to fail its own test. If it fails, no raster is produced; that outcome is the result.",
    "Useful as an independent opinion about magnitude. Cannot be certified by 8 points.",
    "As SoilGrids, with a stronger regional prior and therefore a weaker independence claim.",
    "Retained for context on where peat is deep. Never differenced against a 0-30 cm number.",
    "The most trustworthy numbers in this workflow.",
    "Trustworthy for six cores; a floor, not a value, for the other two."
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(ledger))) {
  log_info(strrep("-", 68))
  log_info(ledger$product[i])
  log_info("  status     : ", ledger$scientific_status[i])
  log_info("  validation : ", ledger$validation_performed[i])
  log_info("  external?  : ", ledger$is_external[i])
  log_info("  leakage    : ", ledger$leakage_risk[i])
}
log_info(strrep("-", 68))

log_warn(strrep("-", 68))
log_warn("THE LEAKAGE THAT WAS AVAILABLE AND WAS NOT TAKEN")
log_warn("  This dataset has 22 segments and 8 cores. Fitting and validating")
log_warn("  at the segment level would have produced a validation sample of")
log_warn("  22 and, because every segment of a core shares one coordinate and")
log_warn("  therefore one covariate pixel, cross-validated skill statistics")
log_warn("  that look strong. They would be measuring the model's ability to")
log_warn("  recall cores it had already been trained on.")
log_warn("  No script in this workflow splits at the segment level. 06")
log_warn("  demonstrates the size of the inflation on a controlled example.")
log_warn(strrep("-", 68))

# =============================================================================
# 2. Benchmark against the published layers
# =============================================================================

ref_path <- file.path(CFG$dir_derived, "04_reference_at_cores.csv")
bench <- NULL

if (!file.exists(ref_path)) {
  log_step("07b  REFERENCE BENCHMARK -- SKIPPED")
  log_warn("No reference extraction found at ", basename(ref_path))
  log_warn("Run 04_gee_reference_audit.R with an authenticated Earth Engine ")
  log_warn("session to produce it. The ledger above is complete without it.")
} else {
  log_step("07b  REFERENCE BENCHMARK")
  ref <- read_csv_verbatim(ref_path)

  # Only the six cores that genuinely reach 30 cm may enter a 0-30 cm
  # comparison. Comparing a lower bound against a reference value would
  # manufacture an apparent negative bias in the reference.
  obs <- stocks[stocks$scenario == "as_supplied" &
                stocks$window == "reference_0_30", ]
  usable <- obs[!obs$is_lower_bound, ]
  log_info(nrow(usable), " of ", nrow(obs),
           " cores are complete at 0-30 cm and may enter the comparison")
  log_warn("Excluded as lower bounds: ",
           paste(obs$core_id[obs$is_lower_bound], collapse = ", "))
  log_warn("Including them would compare a floor against a full value and ",
           "invent a negative bias in the reference layers.")

  m <- merge(usable[, c("core_id", "campaign", "stock_kgm2")], ref,
             by = "core_id", all.x = TRUE)

  comparable <- c("soilgrids_ocs_0_30_kgm2", "sothe2022_soc_0_30_kgm2")
  comparable <- intersect(comparable, names(m))

  if (!isTRUE(m$li2025_comparable_0_30[1])) {
    log_warn("Li et al. 2025 is EXCLUDED from this benchmark by the depth-",
             "support audit in 04. Its depth support is: ",
             m$li2025_depth_support[1])
  }

  bench <- do.call(rbind, lapply(comparable, function(v) {
    p <- m[[v]]
    d <- data.frame(
      layer = v,
      n_cores_compared = sum(is.finite(p) & is.finite(m$stock_kgm2)),
      core_mean_kgm2 = mean(m$stock_kgm2, na.rm = TRUE),
      layer_mean_kgm2 = mean(p, na.rm = TRUE),
      bias_layer_minus_core = bias(m$stock_kgm2, p),
      rmse = rmse(m$stock_kgm2, p),
      pearson_r = suppressWarnings(stats::cor(m$stock_kgm2, p,
                                              use = "complete.obs")),
      validation_status = "QUASI-EXTERNAL local consistency check",
      stringsAsFactors = FALSE)
    d
  }))
  print_table(bench, digits = 3)

  log_warn("With ", nrow(usable), " comparable cores, a correlation ",
           "coefficient here is not a meaningful accuracy statistic.")
  log_warn("Read the BIAS column, which can reveal a gross local offset, and ",
           "ignore the correlation, which cannot be estimated from six ",
           "clustered points.")

  write_csv_logged(bench,
                   file.path(CFG$dir_tables, "07_reference_benchmark.csv"),
                   "published layers vs cores, complete cores only")
}

# =============================================================================
# 3. Model verdict carried into the ledger
# =============================================================================

cv_path <- file.path(CFG$dir_derived, "06_model_cv.csv")
if (file.exists(cv_path)) {
  cv <- read_csv_verbatim(cv_path)
  log_step("07c  MODEL VERDICT")
  print_table(cv, digits = 3)
  loco <- cv[cv$model != "null" & cv$cv == "leave-one-core-out", ]
  if (nrow(loco) && loco$skill_null[1] <= 0) {
    log_warn("The covariate model does not beat a mean-only null under ",
             "core-level cross-validation. PRODUCT 3 is not recommended, and ",
             "the ledger entry above stands as written.")
  }
}

# =============================================================================
# 4. Write
# =============================================================================

log_step("07d  OUTPUTS")
write_csv_logged(ledger, file.path(CFG$dir_tables, "07_validation_ledger.csv"),
                 "what is validated, what is not, and where leakage was possible")
log_ok("07 complete")
