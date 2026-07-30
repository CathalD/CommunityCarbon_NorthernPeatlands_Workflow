# =============================================================================
# run_all.R  --  convenience wrapper. Runs all 7 steps in order.
#
# Recommended for your FIRST time through: run each mvp/R/0N_*.R file
# individually with Rscript instead, so a problem in one step is easy to
# isolate. Come back to this once you trust the pipeline end to end.
# =============================================================================

.args <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.this_dir <- if (length(.this_file)) dirname(.this_file[[1]]) else getwd()

# Each step script normally finds its own directory via commandArgs(), which
# only works when Rscript invokes it directly. Sourced from here instead,
# commandArgs() would still point at THIS file -- so tell each step where it
# really lives via an env var, which the step's bootstrap checks first.
Sys.setenv(MVP_R_DIR = file.path(.this_dir, "R"))

# Steps 01-07 are the mapping pipeline. 08-10 are the external-comparison
# ADD-ON: they read the mapping pipeline's outputs but nothing in 01-07 depends
# on them, so they can be skipped entirely without affecting the map.
steps <- c("01_clean_and_stocks.R", "01b_plot_profiles.R", "02_covariates.R",
          "03_training_data.R", "04_train_model.R", "05_predict_and_compare.R",
          "06_bayesian_update.R", "07_version_and_export.R",
          "08_external_ingest.R", "09_external_ecosystem.R",
          "10_comparison_outputs.R")

for (s in steps) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("RUNNING ", s, "\n")
  cat(strrep("=", 78), "\n")
  source(file.path(.this_dir, "R", s))
}
