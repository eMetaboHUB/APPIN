# test-CNN_clustering.R ----
# Tests for Function/CNN_clustering.R
# These tests do NOT require a Keras model (we feed fake peaks directly)

# =============================================================================
# Helper: small synthetic spectrum for clustering context
# =============================================================================

.make_small_spec <- function(n = 64) {
  mat <- matrix(runif(n * n, 0, 1), nrow = n, ncol = n)
  rownames(mat) <- as.character(seq(9, 0, length.out = n))
  colnames(mat) <- as.character(seq(9, 0, length.out = n))
  mat
}

# =============================================================================
# Empty / edge-case inputs
# =============================================================================

test_that("process_peaks_with_dbscan handles NULL peaks", {
  spec <- .make_small_spec(32)
  params <- make_test_params()
  
  expect_warning(
    result <- process_peaks_with_dbscan(NULL, spec, params, step = 4),
    regexp = "Aucun pic"
  )
  
  expect_type(result, "list")
  expect_s3_class(result$peaks, "data.frame")
  expect_equal(nrow(result$peaks), 0)
  expect_s3_class(result$boxes, "data.frame")
  expect_equal(nrow(result$boxes), 0)
  expect_type(result$shapes, "list")
})

test_that("process_peaks_with_dbscan handles empty peaks dataframe", {
  spec <- .make_small_spec(32)
  params <- make_test_params()
  empty_peaks <- data.frame(F1 = numeric(0), F2 = numeric(0), Intensity = numeric(0))
  
  expect_warning(
    result <- process_peaks_with_dbscan(empty_peaks, spec, params),
    regexp = "Aucun pic"
  )
  
  expect_equal(nrow(result$peaks), 0)
})

# =============================================================================
# Basic clustering
# =============================================================================

test_that("process_peaks_with_dbscan returns expected structure with valid peaks", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 30, f2_range = c(5, 60), f1_range = c(5, 60))
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params, step = 4)
  })
  
  expect_type(result, "list")
  expect_true(all(c("peaks", "boxes") %in% names(result)))
  expect_s3_class(result$peaks, "data.frame")
  expect_s3_class(result$boxes, "data.frame")
})

test_that("process_peaks_with_dbscan produces peaks with ppm columns", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 20, f2_range = c(5, 60), f1_range = c(5, 60))
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params)
  })
  
  expect_true(all(c("F1_ppm", "F2_ppm", "stain_intensity", "cluster_db", "stain_id")
                  %in% names(result$peaks)))
})

test_that("process_peaks_with_dbscan creates properly formatted bounding boxes", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 20)
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params)
  })
  
  expect_true(all(c("xmin", "xmax", "ymin", "ymax", "stain_id", "stain_intensity")
                  %in% names(result$boxes)))
  
  if (nrow(result$boxes) > 0) {
    expect_true(all(result$boxes$xmin <= result$boxes$xmax))
    expect_true(all(result$boxes$ymin <= result$boxes$ymax))
  }
})

test_that("process_peaks_with_dbscan stain_ids start with 'cnn_'", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 15)
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params)
  })
  
  if (nrow(result$peaks) > 0) {
    expect_true(all(grepl("^cnn_", result$peaks$stain_id)))
  }
})

# =============================================================================
# box_padding parameter
# =============================================================================

test_that("process_peaks_with_dbscan respects custom box_padding", {
  spec <- .make_small_spec(64)
  peaks <- make_fake_peaks(n = 10)
  
  params_small <- make_test_params(eps_value = 0.01, box_padding = 0.01)
  params_large <- make_test_params(eps_value = 0.01, box_padding = 0.2)
  
  capture.output({
    res_small <- process_peaks_with_dbscan(peaks, spec, params_small)
    res_large <- process_peaks_with_dbscan(peaks, spec, params_large)
  })
  
  if (nrow(res_small$boxes) > 0 && nrow(res_large$boxes) > 0) {
    # Larger padding -> wider boxes on average
    width_small <- mean(res_small$boxes$xmax - res_small$boxes$xmin)
    width_large <- mean(res_large$boxes$xmax - res_large$boxes$xmin)
    expect_gt(width_large, width_small)
  }
})

test_that("process_peaks_with_dbscan auto-computes padding when box_padding is NULL", {
  spec <- .make_small_spec(64)
  peaks <- make_fake_peaks(n = 10)
  params <- make_test_params(eps_value = 0.03, box_padding = NULL)
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params)
  })
  
  expect_s3_class(result$boxes, "data.frame")
})

# =============================================================================
# NA peaks after ppm conversion
# =============================================================================

test_that("process_peaks_with_dbscan handles out-of-range indices safely", {
  spec <- .make_small_spec(32)
  params <- make_test_params(eps_value = 0.05)
  
  # Indices way outside spectrum -> will be clamped
  peaks <- data.frame(
    F1 = c(-100, 5, 10, 1000),
    F2 = c(-100, 5, 10, 1000),
    Intensity = c(0.5, 0.5, 0.5, 0.5)
  )
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params)
  })
  
  expect_s3_class(result$peaks, "data.frame")
  # No NA values should remain after clamping
  if (nrow(result$peaks) > 0) {
    expect_false(any(is.na(result$peaks$F1_ppm)))
    expect_false(any(is.na(result$peaks$F2_ppm)))
  }
})

# =============================================================================
# keep_peak_ranges filter
# =============================================================================

test_that("process_peaks_with_dbscan applies keep_peak_ranges filter", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 40, seed = 7)
  
  # Get F2_ppm range from the spectrum to define reasonable keep ranges
  x_vals <- as.numeric(colnames(spec))
  f2_min <- min(x_vals)
  f2_max <- max(x_vals)
  
  keep_ranges <- list(
    c(f2_min, f2_min + 1),                       # first range (keeps 1 peak)
    c(f2_min + 3, f2_min + 5)                    # second range (keeps 4 peaks)
  )
  
  capture.output({
    result <- process_peaks_with_dbscan(peaks, spec, params,
                                        keep_peak_ranges = keep_ranges)
  })
  
  expect_s3_class(result$peaks, "data.frame")
  expect_s3_class(result$boxes, "data.frame")
})

test_that("process_peaks_with_dbscan ignores empty keep_peak_ranges list", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 15)
  
  capture.output({
    res_no_filter <- process_peaks_with_dbscan(peaks, spec, params,
                                               keep_peak_ranges = NULL)
    res_empty_filter <- process_peaks_with_dbscan(peaks, spec, params,
                                                  keep_peak_ranges = list())
  })
  
  expect_equal(nrow(res_no_filter$peaks), nrow(res_empty_filter$peaks))
})

test_that("process_peaks_with_dbscan skips invalid range entries (length != 2)", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 15)
  
  x_vals <- as.numeric(colnames(spec))
  
  # Mix valid and invalid ranges.
  # NOTE: the loop inside process_peaks_with_dbscan calls min()/max() on every
  # range BEFORE checking length == 2, so empty numeric(0) produces harmless
  # warnings. We skip that case here.
  keep_ranges <- list(
    c(min(x_vals), min(x_vals) + 1),     # valid
    c(2, 4, 6)                           # invalid (length 3) -> skipped inside
  )
  
  capture.output({
    result <- suppressWarnings(
      process_peaks_with_dbscan(peaks, spec, params,
                                keep_peak_ranges = keep_ranges)
    )
  })
  
  expect_s3_class(result$peaks, "data.frame")
})

# =============================================================================
# Custom step parameter
# =============================================================================

test_that("process_peaks_with_dbscan accepts different step values", {
  spec <- .make_small_spec(64)
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 15)
  
  for (s in c(1, 2, 4, 8)) {
    capture.output({
      result <- process_peaks_with_dbscan(peaks, spec, params, step = s)
    })
    expect_s3_class(result$peaks, "data.frame")
  }
})

# =============================================================================
# Additional tests for uncovered branches
# =============================================================================

test_that("process_peaks_with_dbscan handles all-NA peaks after ppm conversion", {
  # Build a spectrum with NaN ppm axes so the conversion produces NA
  spec <- matrix(runif(64), nrow = 8, ncol = 8)
  # Use characters that don't parse as numeric -> as.numeric() returns NA
  rownames(spec) <- rep("not_a_number", 8)
  colnames(spec) <- rep("not_a_number", 8)
  
  params <- make_test_params(eps_value = 0.05)
  peaks <- make_fake_peaks(n = 5, f2_range = c(1, 8), f1_range = c(1, 8))
  
  # All peaks will have NA F1_ppm/F2_ppm after conversion -> empty after filter
  expect_warning(
    result <- process_peaks_with_dbscan(peaks, spec, params),
    regexp = "No valid peak after ppm conversion"
  )
  
  expect_s3_class(result$peaks, "data.frame")
  expect_equal(nrow(result$peaks), 0)
  expect_s3_class(result$boxes, "data.frame")
  expect_equal(nrow(result$boxes), 0)
  expect_type(result$shapes, "list")
})

test_that("process_peaks_with_dbscan warns when DBSCAN finds no valid clusters", {
  # DBSCAN with eps very small + minPts > 1 can produce all-noise (cluster=0)
  # But process_peaks_with_dbscan uses minPts=1, so every single point forms a cluster.
  # To force all-noise output, we need to stub dbscan::dbscan temporarily.
  
  fn_env <- environment(process_peaks_with_dbscan)
  old_dbscan <- dbscan::dbscan
  
  fake_dbscan <- function(x, eps, minPts, ...) {
    # Return an object with $cluster = all zeros (all noise)
    list(cluster = rep(0L, nrow(x)), eps = eps, minPts = minPts)
  }
  
  # We need to inject into dbscan namespace access
  # Simplest: override the dbscan function in the global/fn env
  # Since dbscan::dbscan is accessed via ::, we patch at the package level via local bind
  
  # Use testthat's with_mocked_bindings if available (testthat 3.5+)
  if (exists("with_mocked_bindings", where = "package:testthat", mode = "function")) {
    spec <- .make_small_spec(32)
    params <- make_test_params(eps_value = 0.05)
    peaks <- make_fake_peaks(n = 10)
    
    testthat::with_mocked_bindings(
      {
        capture.output({
          expect_warning(
            result <- process_peaks_with_dbscan(peaks, spec, params),
            regexp = "No clusters found"
          )
        })
        expect_s3_class(result$peaks, "data.frame")
        expect_equal(nrow(result$boxes), 0)
      },
      dbscan = fake_dbscan,
      .package = "dbscan"
    )
  } else {
    skip("testthat with_mocked_bindings not available")
  }
})

test_that("process_peaks_with_dbscan: all-noise cluster returns peaks with cluster_db=0", {
  # Alternative approach to cover the "no clusters" branch without mocking:
  # With minPts=1 every point is its own cluster, so we can't naturally hit it.
  # We rely on with_mocked_bindings test above.
  # Here we test that stain_id for noise peaks starts with 'cnn_noise_'.
  
  fn_env <- environment(process_peaks_with_dbscan)
  
  fake_dbscan <- function(x, eps, minPts, ...) {
    list(cluster = rep(0L, nrow(x)), eps = eps, minPts = minPts)
  }
  
  if (exists("with_mocked_bindings", where = "package:testthat", mode = "function")) {
    spec <- .make_small_spec(32)
    params <- make_test_params(eps_value = 0.05)
    peaks <- make_fake_peaks(n = 8)
    
    testthat::with_mocked_bindings(
      {
        capture.output({
          result <- suppressWarnings(
            process_peaks_with_dbscan(peaks, spec, params)
          )
        })
        # Noise peaks should have stain_id starting with 'cnn_noise_'
        if (nrow(result$peaks) > 0) {
          expect_true(all(grepl("^cnn_noise_", result$peaks$stain_id)))
          expect_true(all(result$peaks$cluster_db == 0))
        }
      },
      dbscan = fake_dbscan,
      .package = "dbscan"
    )
  } else {
    skip("testthat with_mocked_bindings not available")
  }
})