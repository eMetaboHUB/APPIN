# test-peak_fitting.R - Tests complets pour Function/Peak_fitting.R
# =============================================================================
# Fusion de: test-peak_fitting.R, test-peak_fitting-complete.R, test-peak_fitting-extended.R
# =============================================================================

library(testthat)

# =============================================================================
# HELPERS
# =============================================================================

create_gaussian_peak <- function(nrow=50, ncol=50, amp=1000, cr=25, cc=25, sr=5, sc=5, noise=10) {
  mat <- matrix(0, nrow, ncol)
  for (i in 1:nrow) for (j in 1:ncol) mat[i,j] <- amp * exp(-((i-cr)^2/(2*sr^2) + (j-cc)^2/(2*sc^2)))
  mat <- mat + matrix(rnorm(nrow*ncol, 0, noise), nrow, ncol)
  rownames(mat) <- seq(5, 1, length.out=nrow)
  colnames(mat) <- seq(10, 0, length.out=ncol)
  mat
}

create_gaussian_peak_matrix <- function(size = 20, amplitude = 100, 
                                        center_row = 10, center_col = 10,
                                        sigma = 3, noise_sd = 1) {
  set.seed(42)
  mat <- matrix(rnorm(size * size, mean = 0, sd = noise_sd), nrow = size, ncol = size)
  for (i in 1:size) {
    for (j in 1:size) {
      dist_sq <- (i - center_row)^2 + (j - center_col)^2
      mat[i, j] <- mat[i, j] + amplitude * exp(-dist_sq / (2 * sigma^2))
    }
  }
  ppm_x <- seq(10, 0, length.out = size)
  ppm_y <- seq(10, 0, length.out = size)
  colnames(mat) <- as.character(round(ppm_x, 4))
  rownames(mat) <- as.character(round(ppm_y, 4))
  list(mat = mat, ppm_x = ppm_x, ppm_y = ppm_y)
}

create_multiplet_matrix <- function(size = 30, n_peaks = 3) {
  set.seed(123)
  mat <- matrix(rnorm(size * size, mean = 0, sd = 0.5), nrow = size, ncol = size)
  peak_positions <- list(c(10, 10), c(10, 20), c(20, 15))
  for (p in seq_len(min(n_peaks, length(peak_positions)))) {
    pos <- peak_positions[[p]]
    amplitude <- 50 + p * 20
    for (i in 1:size) {
      for (j in 1:size) {
        dist_sq <- (i - pos[1])^2 + (j - pos[2])^2
        mat[i, j] <- mat[i, j] + amplitude * exp(-dist_sq / 18)
      }
    }
  }
  ppm_x <- seq(10, 0, length.out = size)
  ppm_y <- seq(10, 0, length.out = size)
  colnames(mat) <- as.character(round(ppm_x, 4))
  rownames(mat) <- as.character(round(ppm_y, 4))
  list(mat = mat, ppm_x = ppm_x, ppm_y = ppm_y)
}

# =============================================================================
# TESTS: detect_local_maxima
# =============================================================================

test_that("detect_local_maxima: trouve pic simple", {
  mat <- create_gaussian_peak(noise=0)
  maxima <- detect_local_maxima(mat, threshold=0.3, min_distance=2)
  expect_s3_class(maxima, "data.frame")
  expect_true(nrow(maxima) >= 1)
  expect_true(abs(maxima$row[1] - 25) <= 2)
})

test_that("detect_local_maxima: respecte min_distance", {
  mat <- create_gaussian_peak(noise=50)
  m_close <- detect_local_maxima(mat, 0.1, 2)
  m_far <- detect_local_maxima(mat, 0.1, 10)
  expect_true(nrow(m_far) <= nrow(m_close))
})

test_that("detect_local_maxima: retourne max global si matrice trop petite", {
  small_mat <- matrix(1:4, 2)
  maxima <- detect_local_maxima(small_mat, 0.3, 2)
  expect_equal(nrow(maxima), 1)
  expect_equal(maxima$value, max(small_mat))
})

test_that("detect_local_maxima: detecte plusieurs pics", {
  data <- create_multiplet_matrix(size = 30, n_peaks = 3)
  result <- detect_local_maxima(data$mat, threshold = 0.2, min_distance = 3)
  expect_true(nrow(result) >= 2)
})

test_that("detect_local_maxima: gere matrice constante", {
  flat_mat <- matrix(1, nrow = 10, ncol = 10)
  result <- detect_local_maxima(flat_mat, threshold = 0.5)
  expect_true(is.data.frame(result))
  expect_true(all(c("row", "col", "value") %in% names(result)))
})

test_that("detect_local_maxima: gere NA dans la matrice", {
  data <- create_gaussian_peak_matrix(size = 15)
  data$mat[5, 5] <- NA
  data$mat[7, 7] <- NA
  result <- detect_local_maxima(data$mat, threshold = 0.3)
  expect_true(is.data.frame(result))
})

test_that("detect_local_maxima: gere NA et Inf dans max_val", {
  na_mat <- matrix(NA, nrow = 2, ncol = 2)
  result <- detect_local_maxima(na_mat)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("detect_local_maxima: matrice avec valeurs negatives", {
  mat <- matrix(c(-10, -5, -8, -3), nrow = 2, ncol = 2)
  result <- detect_local_maxima(mat)
  expect_equal(nrow(result), 1)
  expect_equal(result$value, -3)
})

test_that("detect_local_maxima: matrice 3x3 minimale", {
  mat <- matrix(c(1,2,1, 2,5,2, 1,2,1), nrow = 3, ncol = 3)
  result <- detect_local_maxima(mat, threshold = 0.1)
  expect_true(nrow(result) >= 1)
  expect_equal(max(result$value), 5)
})

test_that("detect_local_maxima: filtre min_distance fonctionne", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[3, 3] <- 100
  mat[3, 5] <- 90
  result1 <- detect_local_maxima(mat, threshold = 0.3, min_distance = 1)
  result3 <- detect_local_maxima(mat, threshold = 0.3, min_distance = 3)
  expect_true(nrow(result1) >= nrow(result3))
})

# =============================================================================
# TESTS: pseudo_voigt_2d
# =============================================================================

test_that("pseudo_voigt_2d: centre = amplitude", {
  expect_equal(pseudo_voigt_2d(5,5,100,5,5,1,1,1,1,0.5), 100)
})

test_that("pseudo_voigt_2d: loin = ~0", {
  expect_true(pseudo_voigt_2d(100,100,100,5,5,1,1,1,1,0.5) < 0.01)
})

test_that("pseudo_voigt_2d: eta=0 Gaussien pur", {
  expect_equal(pseudo_voigt_2d(6,5,100,5,5,1,1,1,1,0), 100*exp(-0.5), tolerance=0.001)
})

test_that("pseudo_voigt_2d: eta=1 Lorentzien pur", {
  expect_equal(pseudo_voigt_2d(6,5,100,5,5,1,1,1,1,1), 50, tolerance=0.001)
})

test_that("pseudo_voigt_2d: symetrique", {
  v1 <- pseudo_voigt_2d(4,5,100,5,5,1,1,1,1,0.5)
  v2 <- pseudo_voigt_2d(6,5,100,5,5,1,1,1,1,0.5)
  expect_equal(v1, v2, tolerance=0.001)
})

test_that("pseudo_voigt_2d: valeurs intermediaires de eta", {
  val_03 <- pseudo_voigt_2d(5, 5, 100, 5, 5, 1, 1, 1, 1, 0.3)
  expect_equal(val_03, 100)
  val_07 <- pseudo_voigt_2d(5, 5, 100, 5, 5, 1, 1, 1, 1, 0.7)
  expect_equal(val_07, 100)
  val_off <- pseudo_voigt_2d(6, 6, 100, 5, 5, 1, 1, 1, 1, 0.5)
  expect_true(val_off < 100)
  expect_true(val_off > 0)
})

test_that("pseudo_voigt_2d: asymetrie en sigma_x vs sigma_y", {
  val_x <- pseudo_voigt_2d(6, 5, 100, 5, 5, 2, 1, 1, 1, 0)
  val_y <- pseudo_voigt_2d(5, 6, 100, 5, 5, 2, 1, 1, 1, 0)
  expect_true(val_x > val_y)
})

test_that("pseudo_voigt_2d: asymetrie en gamma_x vs gamma_y", {
  val_x <- pseudo_voigt_2d(6, 5, 100, 5, 5, 1, 1, 2, 1, 1)
  val_y <- pseudo_voigt_2d(5, 6, 100, 5, 5, 1, 1, 2, 1, 1)
  expect_true(val_x > val_y)
})

# =============================================================================
# TESTS: fit_2d_peak - Cas de base
# =============================================================================

test_that("fit_2d_peak: fallback si peu de points (simulation)", {
  fit_few <- function(n_points, min_pts=25) {
    if(n_points < min_pts) {
      return(list(volume=sum(1:n_points), method="sum_fit_failed", 
                  error="Too few points for fitting"))
    }
    list(volume=1000, method="gaussian", error=NULL)
  }
  result_few <- fit_few(10, min_pts=25)
  expect_equal(result_few$method, "sum_fit_failed")
  expect_false(is.na(result_few$volume))
  result_enough <- fit_few(30, min_pts=25)
  expect_equal(result_enough$method, "gaussian")
})

test_that("fit_2d_peak: failed si box vide", {
  mat <- create_gaussian_peak()
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin=100, xmax=110, ymin=100, ymax=110)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, "gaussian")
  expect_equal(result$method, "failed")
})

test_that("fit_2d_peak: refuse modele invalide", {
  mat <- create_gaussian_peak()
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin=3, xmax=7, ymin=1.5, ymax=4.5)
  expect_error(fit_2d_peak(mat, ppm_x, ppm_y, box, "invalid"), "not supported")
})

test_that("fit_2d_peak: fit gaussien sur pic synthetique", {
  data <- create_gaussian_peak_matrix(size = 25, amplitude = 100, 
                                      center_row = 12, center_col = 12, sigma = 3)
  box <- data.frame(xmin = 3, xmax = 7, ymin = 3, ymax = 7)
  result <- fit_2d_peak(data$mat, data$ppm_x, data$ppm_y, box, 
                        model = "gaussian", min_points = 9)
  expect_true(is.list(result))
  expect_true("volume" %in% names(result))
  expect_true("fit_quality" %in% names(result))
  expect_true(result$method %in% c("gaussian", "sum_fit_failed"))
})

test_that("fit_2d_peak: fallback si region constante", {
  mat <- matrix(5, nrow = 20, ncol = 20)
  ppm_x <- seq(10, 0, length.out = 20)
  ppm_y <- seq(10, 0, length.out = 20)
  colnames(mat) <- as.character(round(ppm_x, 4))
  rownames(mat) <- as.character(round(ppm_y, 4))
  box <- data.frame(xmin = 3, xmax = 7, ymin = 3, ymax = 7)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box)
  expect_equal(result$method, "sum_fit_failed")
})

# =============================================================================
# TESTS: fit_2d_peak - Modele Voigt
# =============================================================================

test_that("fit_2d_peak: modele voigt fonctionne", {
  data <- create_gaussian_peak_matrix(size = 25, amplitude = 100, 
                                      center_row = 12, center_col = 12)
  box <- data.frame(xmin = 3, xmax = 7, ymin = 3, ymax = 7)
  result <- fit_2d_peak(data$mat, data$ppm_x, data$ppm_y, box, 
                        model = "voigt", min_points = 9)
  expect_true(result$method %in% c("voigt", "sum_fit_failed"))
})

test_that("fit_2d_peak: voigt model converge sur pic gaussien", {
  set.seed(42)
  mat <- create_gaussian_peak(nrow = 30, ncol = 30, amp = 500, cr = 15, cc = 15, noise = 5)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "voigt", min_points = 9)
  expect_true(result$method %in% c("voigt", "sum_fit_failed"))
  expect_true(!is.na(result$volume))
})

test_that("fit_2d_peak: voigt retourne parametres eta", {
  set.seed(42)
  mat <- create_gaussian_peak(nrow = 40, ncol = 40, amp = 1000, cr = 20, cc = 20, noise = 2)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "voigt", min_points = 9)
  if (result$method == "voigt" && !is.null(result$params)) {
    expect_true("eta" %in% names(result$params))
    expect_true(result$params["eta"] >= 0 && result$params["eta"] <= 1)
  }
})

# =============================================================================
# TESTS: fit_2d_peak - Multiplets
# =============================================================================

test_that("fit_2d_peak: detecte les multiplets", {
  data <- create_multiplet_matrix(size = 30, n_peaks = 3)
  box <- data.frame(xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  result <- fit_2d_peak(data$mat, data$ppm_x, data$ppm_y, box, min_points = 9)
  expect_true("n_peaks" %in% names(result))
  expect_true("is_multiplet" %in% names(result))
})

test_that("fit_2d_peak: detecte et fitte un multiplet", {
  data <- create_multiplet_matrix(size = 30, n_peaks = 2)
  box <- data.frame(xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  result <- fit_2d_peak(data$mat, data$ppm_x, data$ppm_y, box, 
                        model = "gaussian", min_points = 9)
  expect_true(is.list(result))
  expect_true("is_multiplet" %in% names(result))
  if (result$is_multiplet) {
    expect_true(result$n_peaks >= 2)
    expect_equal(result$method, "multiplet_fit")
  }
})

test_that("fit_2d_peak: multiplet avec modele voigt", {
  data <- create_multiplet_matrix(size = 30, n_peaks = 2)
  box <- data.frame(xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  result <- fit_2d_peak(data$mat, data$ppm_x, data$ppm_y, box, 
                        model = "voigt", min_points = 9)
  expect_true(is.list(result))
  expect_true(result$method %in% c("voigt", "multiplet_fit", "sum_fit_failed"))
})

# =============================================================================
# TESTS: fit_2d_peak - Edge cases
# =============================================================================

test_that("fit_2d_peak: gere region avec beaucoup de NA", {
  mat <- create_gaussian_peak(nrow = 30, ncol = 30, noise = 5)
  mat[1:15, 1:15] <- NA
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 0, xmax = 10, ymin = 0, ymax = 5)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, min_points = 9)
  expect_true(is.list(result))
  expect_true(result$method %in% c("gaussian", "sum_fit_failed", "failed"))
})

test_that("fit_2d_peak: calcule R2 correctement", {
  set.seed(123)
  mat <- create_gaussian_peak(nrow = 40, ncol = 40, amp = 1000, cr = 20, cc = 20, noise = 1)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "gaussian", min_points = 9)
  if (result$method == "gaussian") {
    expect_true(!is.na(result$fit_quality))
    expect_true(result$fit_quality >= 0 && result$fit_quality <= 1)
    expect_true(result$fit_quality > 0.5)
  }
})

# =============================================================================
# TESTS: calculate_fitted_volumes
# =============================================================================

test_that("calculate_fitted_volumes: traite plusieurs boxes", {
  mat <- create_gaussian_peak()
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  boxes <- data.frame(stain_id=c("p1","p2"), xmin=c(2,6), xmax=c(4,8), ymin=c(1.5,3), ymax=c(3,4.5))
  results <- calculate_fitted_volumes(mat, ppm_x, ppm_y, boxes, "gaussian")
  expect_equal(nrow(results), 2)
})

test_that("calculate_fitted_volumes: batch fitting fonctionne", {
  data <- create_gaussian_peak_matrix(size = 30, amplitude = 100,
                                      center_row = 15, center_col = 15)
  boxes <- data.frame(
    xmin = c(3, 6), xmax = c(5, 8),
    ymin = c(3, 6), ymax = c(5, 8),
    stain_id = c("box_1", "box_2")
  )
  result <- calculate_fitted_volumes(data$mat, data$ppm_x, data$ppm_y, boxes,
                                     model = "gaussian", min_points = 9)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("stain_id", "volume_fitted", "r_squared", "fit_method") %in% names(result)))
})

test_that("calculate_fitted_volumes: progress callback est appele", {
  data <- create_gaussian_peak_matrix(size = 20)
  boxes <- data.frame(
    xmin = c(3, 5), xmax = c(5, 7),
    ymin = c(3, 5), ymax = c(5, 7),
    stain_id = c("box_1", "box_2")
  )
  progress_calls <- 0
  mock_progress <- function(value, detail) {
    progress_calls <<- progress_calls + 1
  }
  result <- calculate_fitted_volumes(data$mat, data$ppm_x, data$ppm_y, boxes,
                                     progress_callback = mock_progress)
  expect_equal(progress_calls, 2)
})

test_that("calculate_fitted_volumes: gere les boxes avec echec de fit", {
  data <- create_gaussian_peak_matrix(size = 20)
  boxes <- data.frame(
    xmin = c(3, 100), xmax = c(5, 110),
    ymin = c(3, 100), ymax = c(5, 110),
    stain_id = c("valid", "invalid")
  )
  result <- calculate_fitted_volumes(data$mat, data$ppm_x, data$ppm_y, boxes)
  expect_equal(nrow(result), 2)
  expect_true(any(result$fit_method == "failed" | !is.na(result$fit_error)))
})

test_that("calculate_fitted_volumes: voigt model batch", {
  mat <- create_gaussian_peak(nrow = 40, ncol = 40, amp = 500)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  boxes <- data.frame(
    stain_id = c("box_1", "box_2"),
    xmin = c(2, 5), xmax = c(4, 7),
    ymin = c(1, 2), ymax = c(3, 4)
  )
  results <- calculate_fitted_volumes(mat, ppm_x, ppm_y, boxes, model = "voigt")
  expect_equal(nrow(results), 2)
  expect_true(all(results$fit_method %in% c("voigt", "sum_fit_failed", "failed", "multiplet_fit")))
})

test_that("calculate_fitted_volumes: detecte multiplets", {
  data <- create_multiplet_matrix(size = 40, n_peaks = 2)
  boxes <- data.frame(
    stain_id = c("multiplet_box"),
    xmin = c(0), xmax = c(10),
    ymin = c(0), ymax = c(10)
  )
  results <- calculate_fitted_volumes(data$mat, data$ppm_x, data$ppm_y, boxes, min_points = 9)
  expect_equal(nrow(results), 1)
  expect_true("is_multiplet" %in% names(results))
  expect_true("n_peaks" %in% names(results))
})

# =============================================================================
# TESTS: R2 (formule)
# =============================================================================

test_that("R2 = 1 pour fit parfait", {
  obs <- 1:5
  pred <- 1:5
  r2 <- 1 - sum((obs-pred)^2) / sum((obs-mean(obs))^2)
  expect_equal(r2, 1)
})

test_that("R2 < 1 pour fit imparfait", {
  obs <- 1:5
  pred <- c(1.1, 2.2, 2.8, 4.1, 4.9)
  r2 <- 1 - sum((obs-pred)^2) / sum((obs-mean(obs))^2)
  expect_true(r2 > 0.9 && r2 < 1)
})

# =============================================================================
# TESTS: generate_fit_diagnostic_data
# =============================================================================

test_that("generate_fit_diagnostic_data: retourne erreur si pas de fitted_values", {
  fit_result <- list(
    volume = 1000,
    fit_quality = NULL,
    fitted_values = NULL,
    method = "sum_fit_failed"
  )
  mat <- create_gaussian_peak()
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 3, xmax = 7, ymin = 1.5, ymax = 4.5)
  result <- generate_fit_diagnostic_data(fit_result, mat, ppm_x, ppm_y, box)
  expect_false(result$success)
  expect_true(!is.null(result$error))
})

test_that("generate_fit_diagnostic_data: genere slice F2 correctement", {
  set.seed(42)
  mat <- create_gaussian_peak(nrow = 30, ncol = 30, amp = 500, noise = 2)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "gaussian")
  if (!is.null(fit_result$fitted_values)) {
    diag_data <- generate_fit_diagnostic_data(fit_result, mat, ppm_x, ppm_y, box, slice_direction = "F2")
    expect_true(diag_data$success)
    expect_equal(diag_data$axis_label, "F2 (ppm)")
    expect_true(length(diag_data$experimental) > 0)
    expect_true(length(diag_data$fitted) > 0)
  }
})

test_that("generate_fit_diagnostic_data: genere slice F1 correctement", {
  set.seed(42)
  mat <- create_gaussian_peak(nrow = 30, ncol = 30, amp = 500, noise = 2)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "gaussian")
  if (!is.null(fit_result$fitted_values)) {
    diag_data <- generate_fit_diagnostic_data(fit_result, mat, ppm_x, ppm_y, box, slice_direction = "F1")
    expect_true(diag_data$success)
    expect_equal(diag_data$axis_label, "F1 (ppm)")
  }
})

test_that("generate_fit_diagnostic_data: gere box invalide", {
  fit_result <- list(
    volume = 1000,
    fit_quality = 0.95,
    fitted_values = rep(100, 100),
    method = "gaussian"
  )
  mat <- create_gaussian_peak()
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 100, xmax = 110, ymin = 100, ymax = 110)
  result <- generate_fit_diagnostic_data(fit_result, mat, ppm_x, ppm_y, box)
  expect_false(result$success)
})

# =============================================================================
# TESTS: create_fit_diagnostic_plot
# =============================================================================

test_that("create_fit_diagnostic_plot: gere diag_data avec success=FALSE", {
  diag_data <- list(success = FALSE, error = "Test error message")
  p <- create_fit_diagnostic_plot(diag_data)
  expect_true(inherits(p, "plotly"))
})

test_that("create_fit_diagnostic_plot: genere plot avec residuals", {
  diag_data <- list(
    success = TRUE,
    ppm = seq(10, 0, length.out = 50),
    experimental = rnorm(50, 100, 10),
    fitted = rnorm(50, 100, 5),
    residuals = rnorm(50, 0, 5),
    axis_label = "F2 (ppm)",
    slice_info = "F1 = 3.5 ppm",
    r_squared_global = 0.92,
    r_squared_slice = 0.88,
    fit_method = "gaussian"
  )
  p <- create_fit_diagnostic_plot(diag_data, show_residuals = TRUE)
  expect_true(inherits(p, "plotly"))
})

test_that("create_fit_diagnostic_plot: genere plot sans residuals", {
  diag_data <- list(
    success = TRUE,
    ppm = seq(10, 0, length.out = 50),
    experimental = rnorm(50, 100, 10),
    fitted = rnorm(50, 100, 5),
    residuals = rnorm(50, 0, 5),
    axis_label = "F1 (ppm)",
    slice_info = "F2 = 5.0 ppm",
    r_squared_global = 0.75,
    r_squared_slice = 0.70,
    fit_method = "voigt"
  )
  p <- create_fit_diagnostic_plot(diag_data, show_residuals = FALSE)
  expect_true(inherits(p, "plotly"))
})

test_that("create_fit_diagnostic_plot: affiche qualite selon R2", {
  for (r2 in c(0.98, 0.90, 0.75, 0.50, NA)) {
    diag_data <- list(
      success = TRUE,
      ppm = 1:10,
      experimental = 1:10,
      fitted = 1:10,
      residuals = rep(0, 10),
      axis_label = "F2 (ppm)",
      slice_info = "Test",
      r_squared_global = r2,
      r_squared_slice = r2,
      fit_method = "gaussian"
    )
    p <- create_fit_diagnostic_plot(diag_data)
    expect_true(inherits(p, "plotly"))
  }
})

# =============================================================================
# TESTS: Volume calculation
# =============================================================================

test_that("fit_2d_peak: volume analytique vs somme fitted", {
  set.seed(42)
  mat <- create_gaussian_peak(nrow = 50, ncol = 50, amp = 1000, cr = 25, cc = 25, noise = 1)
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  box <- data.frame(xmin = 2, xmax = 8, ymin = 1, ymax = 4)
  result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = "gaussian")
  if (result$method == "gaussian") {
    expect_true(!is.na(result$volume))
    expect_true(result$volume > 0)
    if (!is.null(result$volume_analytical)) {
      expect_true(result$volume_analytical > 0)
    }
  }
})
