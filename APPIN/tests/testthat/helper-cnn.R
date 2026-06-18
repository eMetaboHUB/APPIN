# helper-cnn.R ----
# Shared helpers for CNN tests
# Loaded automatically by testthat before any test-*.R file

# =============================================================================
# LAZY MODEL LOADING
# =============================================================================
# We load the CNN model once and cache it in a helper env to avoid reloading
# it for every test (Keras model loading is expensive).

.cnn_test_env <- new.env(parent = emptyenv())

#' Get a CNN model for testing (cached)
#'
#' Tries to load the real TOCSY model. If unavailable (no weights, no keras),
#' returns NULL so individual tests can skip gracefully.
get_test_cnn_model <- function(spectrum_type = "TOCSY") {
  key <- paste0("model_", spectrum_type)
  
  if (exists(key, envir = .cnn_test_env)) {
    return(get(key, envir = .cnn_test_env))
  }
  
  model <- tryCatch(
    get_cnn_model(spectrum_type),
    error = function(e) NULL
  )
  
  assign(key, model, envir = .cnn_test_env)
  return(model)
}

#' Skip test if CNN model is unavailable
skip_if_no_cnn_model <- function(spectrum_type = "TOCSY") {
  model <- get_test_cnn_model(spectrum_type)
  if (is.null(model)) {
    testthat::skip(sprintf("CNN model (%s) not available in this environment", spectrum_type))
  }
  invisible(model)
}

# =============================================================================
# SYNTHETIC SPECTRUM GENERATION
# =============================================================================

#' Generate a small synthetic 2D NMR spectrum for testing
#'
#' @param n_row Number of F2 points (rows)
#' @param n_col Number of F1 points (cols)
#' @param ppm_range_f2 Range for F2 axis
#' @param ppm_range_f1 Range for F1 axis
#' @param n_peaks Number of Gaussian peaks to inject
#' @param seed Random seed for reproducibility
#' @return A matrix with numeric rownames/colnames (ppm values)
make_synthetic_spectrum <- function(n_row = 64, n_col = 64,
                                    ppm_range_f2 = c(0.5, 9.5),
                                    ppm_range_f1 = c(0.5, 9.5),
                                    n_peaks = 5,
                                    seed = 42) {
  set.seed(seed)
  
  # Descending ppm (NMR convention: left = high ppm)
  f2_ppm <- seq(ppm_range_f2[2], ppm_range_f2[1], length.out = n_row)
  f1_ppm <- seq(ppm_range_f1[2], ppm_range_f1[1], length.out = n_col)
  
  # Background noise
  mat <- matrix(rnorm(n_row * n_col, mean = 0, sd = 0.01), nrow = n_row, ncol = n_col)
  
  # Add synthetic Gaussian peaks along and off the diagonal
  for (p in seq_len(n_peaks)) {
    cr <- sample(seq_len(n_row), 1)
    cc <- sample(seq_len(n_col), 1)
    amp <- runif(1, 0.5, 1.0)
    sigma <- sample(2:4, 1)
    
    for (i in seq_len(n_row)) {
      for (j in seq_len(n_col)) {
        mat[i, j] <- mat[i, j] + amp * exp(-((i - cr)^2 + (j - cc)^2) / (2 * sigma^2))
      }
    }
  }
  
  rownames(mat) <- as.character(f2_ppm)
  colnames(mat) <- as.character(f1_ppm)
  mat
}

#' Default parameter list for pipeline tests
make_test_params <- function(...) {
  defaults <- list(
    eps_value = 0.01,
    pred_class_thres = 0.01,
    int_thres = 0.001,
    trace_filter_ratio = 0.1,
    use_filters = FALSE,
    disable_clustering = FALSE,
    box_padding = NULL
  )
  overrides <- list(...)
  modifyList(defaults, overrides)
}

#' Generate a fake detected_peaks data frame (bypasses CNN)
#' Useful for testing clustering/filtering without running the model
make_fake_peaks <- function(n = 20, seed = 123,
                            f2_range = c(1, 64), f1_range = c(1, 64)) {
  set.seed(seed)
  data.frame(
    F2 = round(runif(n, f2_range[1], f2_range[2])),
    F1 = round(runif(n, f1_range[1], f1_range[2])),
    Intensity = runif(n, 0.1, 1.0)
  )
}
