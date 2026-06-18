# Test coverage

library(covr)

# ─────────────────────────────────────────────────────────────────────────────
# Load test helpers BEFORE file_coverage() runs the tests.
# covr::file_coverage() uses sys.source() internally, which does not
# auto-load helper-*.R files like testthat::test_dir() would. So we source
# them manually into the global env so tests can see them.
# ─────────────────────────────────────────────────────────────────────────────
helper_files <- list.files("tests/testthat", pattern = "^helper-.*\\.R$", full.names = TRUE)
for (hf in helper_files) {
  source(hf)
}

# Exclude regression tests
test_files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
test_files <- test_files[!grepl("regression", test_files)]

cov <- file_coverage(
  source_files = c(
    "Function/Read_2DNMR_spectrum.R",
    "Function/Vizualisation.R",
    "Function/Peak_picking.R",
    "Function/Peak_fitting.R",
    # CNN module (split into sub-files)
    "Function/CNN_shiny.R",
    "Function/CNN_model.R",
    "Function/CNN_detection.R",
    "Function/CNN_filtering.R",
    "Function/CNN_clustering.R",
    "Function/CNN_main.R",
    "R/utils.R"
  ),
  test_files = test_files
)

percent_coverage(cov)
report(cov)