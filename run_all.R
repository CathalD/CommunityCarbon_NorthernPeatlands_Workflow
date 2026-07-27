# =============================================================================
# run_all.R  --  Run the workflow end to end.
#
#   Rscript run_all.R          run everything that can run
#   Rscript run_all.R --gee    include the Earth Engine scripts (03, 04, 06)
#
# Scripts 01, 02, 05, 07 and 08 need base R only. Scripts 03, 04 and 06 need
# rgee and an authenticated Earth Engine session; they are skipped by default
# and the pipeline completes without them, producing Products 1 and 2 plus the
# full validation ledger.
# =============================================================================

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

args    <- commandArgs(trailingOnly = TRUE)
use_gee <- "--gee" %in% args

base_scripts <- c("01_ingest_qc.R", "02_depth_harmonise.R",
                  "05_stratified_estimate.R")
gee_scripts  <- c("03_gee_covariates.R", "04_gee_reference_audit.R",
                  "06_covariate_model.R")
tail_scripts <- c("07_validation_ledger.R", "08_report.R")

plan <- if (use_gee) c(base_scripts, gee_scripts, tail_scripts) else
                     c(base_scripts, tail_scripts)

log_step("RUN ALL  --  ", length(plan), " scripts")
if (!use_gee) {
  log_info("Earth Engine scripts are SKIPPED. Pass --gee to include them.")
  log_info("Products 1 and 2 and the validation ledger do not require them.")
}
log_info("plan: ", paste(plan, collapse = " -> "))

t0 <- Sys.time()
failed <- character(0)

for (s in plan) {
  p <- file.path(CFG$root, "scripts", s)
  log_step("RUN  ", s)
  st <- system2("Rscript", p, env = paste0("CCNP_ROOT=", CFG$root))
  if (!identical(st, 0L)) {
    failed <- c(failed, s)
    if (s %in% gee_scripts) {
      log_warn(s, " failed. This is expected without an authenticated Earth ",
               "Engine session; continuing.")
    } else {
      fail_loudly("A core script failed", paste0("Script: ", s),
                  remedy = "Read the error above; it names the problem and the fix.")
    }
  }
}

log_step("DONE in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
if (length(failed)) {
  log_warn("skipped or failed: ", paste(failed, collapse = ", "))
}
log_info("report : ", file.path(CFG$root, "outputs", "REPORT.md"))
log_info("tables : ", CFG$dir_tables)
log_info("derived: ", CFG$dir_derived)
