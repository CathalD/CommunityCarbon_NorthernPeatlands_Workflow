# Run once: install everything the MVP pipeline needs.
pkgs <- c("readr", "dplyr", "sf", "terra", "rgee", "ranger",
         "rsample", "yardstick", "purrr", "tibble")

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
} else {
  message("All required packages are already installed.")
}

message("\nEarth Engine (needed only for step 02) is separate:")
message('  rgee::ee_install()')
message('  rgee::ee_Initialize(project = "your-gcp-project-id")')
