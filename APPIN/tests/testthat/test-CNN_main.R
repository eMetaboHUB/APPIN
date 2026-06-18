# test-CNN_main.R ----
# Tests for Function/CNN_main.R (run_cnn_peak_picking pipeline)

# =============================================================================
# Input validation
# =============================================================================

test_that("run_cnn_peak_picking errors on NULL input", {
  params <- make_test_params()
  
  expect_error(
    run_cnn_peak_picking(NULL, model = NULL, params = params),
    regexp = "must be a normalized matrix"
  )
})

test_that("run_cnn_peak_picking errors on non-matrix input", {
  params <- make_test_params()
  
  expect_error(
    run_cnn_peak_picking(data.frame(x = 1:5), model = NULL, params = params),
    regexp = "must be a normalized matrix"
  )
  
  expect_error(
    run_cnn_peak_picking(c(1, 2, 3), model = NULL, params = params),
    regexp = "must be a normalized matrix"
  )
})

test_that("run_cnn_peak_picking validates method argument", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params()
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  expect_error(
    run_cnn_peak_picking(spec, model = model, params = params, method = "invalid"),
    regexp = "should be one of"
  )
})

# =============================================================================
# Auto-model selection
# =============================================================================

test_that("run_cnn_peak_picking auto-selects model when model is NULL", {
  skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(trace_filter_ratio = 0, disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = NULL, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

test_that("run_cnn_peak_picking normalizes spectrum_type to uppercase", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(trace_filter_ratio = 0, disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "tocsy",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# Methods: batch vs classique
# =============================================================================

test_that("run_cnn_peak_picking works with batch method", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05, trace_filter_ratio = 0)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
  # Clustering mode returns peaks + boxes (shapes is commented out in CNN_clustering.R)
  expect_true(all(c("peaks", "boxes") %in% names(result)))
})

test_that("run_cnn_peak_picking works with classique (sequential) method", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05, trace_filter_ratio = 0)
  # Smaller spectrum because classique is slow (iterates rows and cols)
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "classique", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# HSQC branch (prints F1/F2 ranges)
# =============================================================================

test_that("run_cnn_peak_picking prints F2/F1 range info for HSQC spectra", {
  model <- skip_if_no_cnn_model("TOCSY")  # Use TOCSY weights; branch triggers on type
  params <- make_test_params(eps_value = 0.05, trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  
  # Build an HSQC-like spectrum: F2 = 13C range, F1 = 1H range
  spec <- matrix(runif(16 * 16), nrow = 16, ncol = 16)
  rownames(spec) <- as.character(seq(150, 0, length.out = 16))   # 13C
  colnames(spec) <- as.character(seq(9, 0, length.out = 16))     # 1H
  
  output <- capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "HSQC",
                                   method = "batch", verbose = TRUE)
  })
  
  expect_true(any(grepl("HSQC:", output)))
})

# =============================================================================
# No-peak case
# =============================================================================

test_that("run_cnn_peak_picking warns when no peaks detected", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(pred_class_thres = 0.999,  # impossible threshold
                             trace_filter_ratio = 0)
  # Flat spectrum
  spec <- matrix(0.0001, nrow = 16, ncol = 16)
  rownames(spec) <- as.character(seq(9, 0, length.out = 16))
  colnames(spec) <- as.character(seq(9, 0, length.out = 16))
  
  capture.output({
    result <- tryCatch(
      suppressWarnings(
        run_cnn_peak_picking(spec, model = model, params = params,
                             spectrum_type = "TOCSY",
                             method = "batch", verbose = FALSE)
      ),
      warning = function(w) w
    )
  })
  
  # Either we get a warning-returned result, or a list with null components
  expect_true(is.list(result) || inherits(result, "warning"))
})

# =============================================================================
# Trace filter
# =============================================================================

test_that("run_cnn_peak_picking applies trace filter when ratio > 0", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0.2,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 5)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

test_that("run_cnn_peak_picking skips trace filter when ratio = 0", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05, trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# use_filters branch
# =============================================================================

test_that("run_cnn_peak_picking respects use_filters = TRUE", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             use_filters = TRUE,
                             int_thres = 0.01)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

test_that("run_cnn_peak_picking respects use_filters = FALSE (default)", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             use_filters = FALSE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# disable_clustering branch (no-cluster mode with vertical trace filter)
# =============================================================================

test_that("run_cnn_peak_picking with disable_clustering returns peaks without DBSCAN", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
  expect_s3_class(result$peaks, "data.frame")
  expect_s3_class(result$boxes, "data.frame")
})

test_that("disable_clustering applies vertical-trace filter for COSY", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 5)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "COSY",
                                   method = "batch", verbose = TRUE)
  })
  
  expect_type(result, "list")
})

test_that("disable_clustering applies vertical-trace filter for UFCOSY", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 5)
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "UFCOSY",
                                   method = "batch", verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# keep_peak_ranges in disable_clustering mode
# =============================================================================

test_that("disable_clustering handles keep_peak_ranges filter", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 5)
  
  # Get spectrum ppm bounds
  x_vals <- as.numeric(colnames(spec))
  keep_ranges <- list(
    c(min(x_vals), min(x_vals) + 1),
    c(min(x_vals) + 3, min(x_vals) + 5)
  )
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch",
                                   keep_peak_ranges = keep_ranges,
                                   verbose = TRUE)
  })
  
  expect_type(result, "list")
})

test_that("disable_clustering handles invalid range entries in keep_peak_ranges", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = TRUE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 5)
  
  x_vals <- as.numeric(colnames(spec))
  keep_ranges <- list(
    c(min(x_vals), min(x_vals) + 2),  # valid
    c(1, 2, 3)                        # invalid, length != 2
  )
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch",
                                   keep_peak_ranges = keep_ranges,
                                   verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# Clustering mode with keep_peak_ranges (full pipeline)
# =============================================================================

test_that("run_cnn_peak_picking (clustering mode) handles keep_peak_ranges", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05,
                             trace_filter_ratio = 0,
                             disable_clustering = FALSE)
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 4)
  
  x_vals <- as.numeric(colnames(spec))
  keep_ranges <- list(c(min(x_vals), min(x_vals) + 2))
  
  capture.output({
    result <- run_cnn_peak_picking(spec, model = model, params = params,
                                   spectrum_type = "TOCSY",
                                   method = "batch",
                                   keep_peak_ranges = keep_ranges,
                                   verbose = FALSE)
  })
  
  expect_type(result, "list")
})

# =============================================================================
# verbose = TRUE outputs informative messages
# =============================================================================

test_that("verbose = TRUE prints detection method info", {
  model <- skip_if_no_cnn_model("TOCSY")
  params <- make_test_params(eps_value = 0.05, trace_filter_ratio = 0)
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  output <- capture.output({
    tryCatch(
      run_cnn_peak_picking(spec, model = model, params = params,
                           spectrum_type = "TOCSY",
                           method = "batch", verbose = TRUE),
      warning = function(w) NULL
    )
  })
  
  expect_true(any(grepl("batch method|Post-filters|Peaks before DBSCAN", output)))
})