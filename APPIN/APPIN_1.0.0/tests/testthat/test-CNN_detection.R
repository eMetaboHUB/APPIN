# test-CNN_detection.R ----
# Tests for Function/CNN_detection.R

# =============================================================================
# pad_sequence()
# =============================================================================

test_that("pad_sequence pads shorter vectors with zeros", {
  result <- pad_sequence(c(1, 2, 3), 5)
  expect_length(result, 5)
  expect_equal(result, c(1, 2, 3, 0, 0))
})

test_that("pad_sequence returns vector unchanged if already at target length", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(pad_sequence(x, 5), x)
})

test_that("pad_sequence returns vector unchanged if longer than target", {
  x <- c(1, 2, 3, 4, 5, 6, 7)
  expect_equal(pad_sequence(x, 5), x)
})

test_that("pad_sequence handles empty input", {
  result <- pad_sequence(numeric(0), 3)
  expect_length(result, 3)
  expect_equal(result, c(0, 0, 0))
})

test_that("pad_sequence preserves negative values", {
  result <- pad_sequence(c(-1, -2), 4)
  expect_equal(result, c(-1, -2, 0, 0))
})

# =============================================================================
# get_detected_peaks_with_intensity() - Sequential method
# =============================================================================

test_that("get_detected_peaks_with_intensity returns a valid dataframe", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  # Suppress progress bar output
  result <- suppressWarnings(capture.output({
    peaks <- get_detected_peaks_with_intensity(spec, model, target_length = 2048)
  }))
  
  expect_s3_class(peaks, "data.frame")
  expect_true(all(c("F2", "F1", "Intensity") %in% names(peaks)))
})

test_that("get_detected_peaks_with_intensity removes duplicates", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  capture.output({
    peaks <- get_detected_peaks_with_intensity(spec, model)
  })
  
  # Unique should match self
  expect_equal(nrow(peaks), nrow(unique(peaks)))
})

# =============================================================================
# predict_peaks_1D_batch() - Batch method
# =============================================================================

test_that("predict_peaks_1D_batch returns expected dataframe columns (rows axis)", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    peaks <- predict_peaks_1D_batch(spec, model, axis = "rows",
                                    threshold_class = 0.01,
                                    batch_size = 8, verbose = FALSE)
  })
  
  expect_s3_class(peaks, "data.frame")
  expect_true(all(c("F1", "F2", "Intensity", "ppm") %in% names(peaks)))
})

test_that("predict_peaks_1D_batch works on columns axis", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  spec <- make_synthetic_spectrum(n_row = 32, n_col = 32, n_peaks = 3)
  
  capture.output({
    peaks <- predict_peaks_1D_batch(spec, model, axis = "columns",
                                    threshold_class = 0.01,
                                    batch_size = 8, verbose = FALSE)
  })
  
  expect_s3_class(peaks, "data.frame")
  expect_true(all(c("F1", "F2", "Intensity") %in% names(peaks)))
})

test_that("predict_peaks_1D_batch uses sliding window for large spectra", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  # n_col > model_input_length (2048) triggers sliding window
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 2100, n_peaks = 5)
  
  capture.output({
    peaks <- predict_peaks_1D_batch(spec, model, axis = "rows",
                                    threshold_class = 0.01,
                                    batch_size = 4, verbose = FALSE)
  })
  
  expect_s3_class(peaks, "data.frame")
})

test_that("predict_peaks_1D_batch returns empty dataframe when nothing detected", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  # Flat spectrum with no peaks -> no detections with high threshold
  spec <- matrix(0.001, nrow = 16, ncol = 16)
  rownames(spec) <- as.character(seq(9.5, 0.5, length.out = 16))
  colnames(spec) <- as.character(seq(9.5, 0.5, length.out = 16))
  
  capture.output({
    peaks <- predict_peaks_1D_batch(spec, model, axis = "rows",
                                    threshold_class = 0.99,
                                    batch_size = 4, verbose = FALSE)
  })
  
  expect_s3_class(peaks, "data.frame")
  expect_true(all(c("F1", "F2", "Intensity", "ppm") %in% names(peaks)))
})

test_that("predict_peaks_1D_batch validates axis argument", {
  model <- skip_if_no_cnn_model("TOCSY")
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  expect_error(
    predict_peaks_1D_batch(spec, model, axis = "bogus", verbose = FALSE),
    regexp = "should be one of"
  )
})

test_that("predict_peaks_1D_batch output is deduplicated", {
  model <- skip_if_no_cnn_model("TOCSY")
  
  spec <- make_synthetic_spectrum(n_row = 16, n_col = 16, n_peaks = 2)
  
  capture.output({
    peaks <- predict_peaks_1D_batch(spec, model, axis = "rows",
                                    threshold_class = 0.01,
                                    batch_size = 4, verbose = FALSE)
  })
  
  expect_equal(nrow(peaks), nrow(unique(peaks)))
})
