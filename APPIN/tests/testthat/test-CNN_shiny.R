# test-CNN_shiny.R ----
# Tests for Function/CNN_shiny.R
# This file just loads dependencies and sources other CNN modules.
# The main thing to verify is that after sourcing it, all public functions
# from every sub-module are available in the environment.

test_that("CNN_shiny.R exposes all expected functions after sourcing", {
  # Functions from CNN_model.R
  expect_true(exists("build_peak_predictor", mode = "function"))
  expect_true(exists("get_cnn_model", mode = "function"))
  expect_true(exists("check_cnn_models", mode = "function"))
  expect_true(exists("clear_cnn_cache", mode = "function"))
  
  # Functions from CNN_detection.R
  expect_true(exists("pad_sequence", mode = "function"))
  expect_true(exists("get_detected_peaks_with_intensity", mode = "function"))
  expect_true(exists("predict_peaks_1D_batch", mode = "function"))
  
  # Functions from CNN_filtering.R
  expect_true(exists("filter_peaks_by_proportion", mode = "function"))
  expect_true(exists("filter_noisy_columns", mode = "function"))
  expect_true(exists("clean_peak_clusters_dbscan", mode = "function"))
  expect_true(exists("remove_peaks_ppm_range", mode = "function"))
  
  # Functions from CNN_clustering.R
  expect_true(exists("process_peaks_with_dbscan", mode = "function"))
  
  # Functions from CNN_main.R
  expect_true(exists("run_cnn_peak_picking", mode = "function"))
})

test_that("CNN_MODEL_PATHS global is accessible after sourcing", {
  expect_true(exists("CNN_MODEL_PATHS"))
  expect_type(CNN_MODEL_PATHS, "list")
})

test_that(".cnn_models_cache environment exists", {
  expect_true(exists(".cnn_models_cache"))
  expect_true(is.environment(.cnn_models_cache))
})

test_that("key dependency libraries are loaded", {
  # These are loaded at the top of CNN_shiny.R
  required_pkgs <- c("purrr", "ggplot2", "dplyr", "dbscan", "plotly")
  for (pkg in required_pkgs) {
    expect_true(
      paste0("package:", pkg) %in% search() ||
      requireNamespace(pkg, quietly = TRUE),
      info = paste("Package not loaded:", pkg)
    )
  }
})
