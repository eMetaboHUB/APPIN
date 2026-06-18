# test-regression-CNN.R ----
# Tests de régression pour le CNN peak picking.
# Compare l'output actuel au snapshot de référence généré par
# generate_test_snapshots.R.
#
# Si ces tests échouent après un changement de code, soit:
#   1. Le changement est une amélioration -> régénérer le snapshot
#      (source("tests/generate_test_snapshots.R"); generate_snapshot_cnn())
#   2. Le changement est une régression non voulue -> debugger

SNAPSHOTS_DIR <- "tests/snapshots"
FIXTURES_DIR <- "tests/fixtures"

# Helper: charger un snapshot avec skip gracieux
.load_cnn_snapshot <- function() {
  snapshot_path <- file.path(SNAPSHOTS_DIR, "cnn_peak_picking_ufcosy.rds")
  if (!file.exists(snapshot_path)) {
    testthat::skip(sprintf("Snapshot non trouvé: %s", snapshot_path))
  }
  readRDS(snapshot_path)
}

# Helper: run le CNN peak picking avec les mêmes params que le snapshot
.run_cnn_on_fixture <- function() {
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  if (!dir.exists(ufcosy_path)) {
    testthat::skip("Fixture UFCOSY non disponible")
  }
  
  model <- tryCatch(get_cnn_model("UFCOSY"),
                    error = function(e) NULL)
  if (is.null(model)) {
    testthat::skip("Modèle CNN UFCOSY non chargeable")
  }
  
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  rr_norm <- mat / max(abs(mat))
  
  params <- list(
    eps_value = 0.01,
    pred_class_thres = 0.01,
    int_thres = 0.001,
    trace_filter_ratio = 0.1,
    use_filters = FALSE,
    disable_clustering = FALSE,
    box_padding = NULL
  )
  
  capture.output({
    cnn_result <- run_cnn_peak_picking(
      rr_norm,
      model = model,
      params = params,
      spectrum_type = "UFCOSY",
      method = "batch",
      verbose = FALSE
    )
  })
  
  cnn_result
}

# =============================================================================
# TESTS DE RÉGRESSION
# =============================================================================

test_that("RÉGRESSION CNN: nombre de peaks détectés reste cohérent", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  expected_n <- snapshot$data$peak_stats$n
  actual_n <- if (!is.null(actual$peaks)) nrow(actual$peaks) else 0
  
  # Tolérance: ±10% (les modèles CNN ont une petite variance non-déterministe)
  tolerance <- max(5, ceiling(expected_n * 0.10))
  
  expect_lte(abs(actual_n - expected_n), tolerance,
             info = sprintf("Expected ~%d peaks (±%d), got %d",
                            expected_n, tolerance, actual_n))
})

test_that("RÉGRESSION CNN: nombre de bounding boxes reste cohérent", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  expected_n <- snapshot$data$box_stats$n
  actual_n <- if (!is.null(actual$boxes)) nrow(actual$boxes) else 0
  
  tolerance <- max(5, ceiling(expected_n * 0.10))
  
  expect_lte(abs(actual_n - expected_n), tolerance,
             info = sprintf("Expected ~%d boxes (±%d), got %d",
                            expected_n, tolerance, actual_n))
})

test_that("RÉGRESSION CNN: plage F2_ppm des peaks reste dans le domaine attendu", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  if (is.null(actual$peaks) || nrow(actual$peaks) == 0) {
    skip("Pas de peaks détectés, impossible de comparer les plages")
  }
  
  expected_range <- snapshot$data$peak_stats$f2_ppm_range
  actual_range <- range(actual$peaks$F2_ppm, na.rm = TRUE)
  
  # Tolérance de 0.2 ppm sur les bornes
  expect_lte(abs(actual_range[1] - expected_range[1]), 0.2)
  expect_lte(abs(actual_range[2] - expected_range[2]), 0.2)
})

test_that("RÉGRESSION CNN: plage F1_ppm des peaks reste dans le domaine attendu", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  if (is.null(actual$peaks) || nrow(actual$peaks) == 0) {
    skip("Pas de peaks détectés, impossible de comparer les plages")
  }
  
  expected_range <- snapshot$data$peak_stats$f1_ppm_range
  actual_range <- range(actual$peaks$F1_ppm, na.rm = TRUE)
  
  expect_lte(abs(actual_range[1] - expected_range[1]), 0.2)
  expect_lte(abs(actual_range[2] - expected_range[2]), 0.2)
})

test_that("RÉGRESSION CNN: intensité moyenne des peaks reste stable", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  if (is.null(actual$peaks) || nrow(actual$peaks) == 0) {
    skip("Pas de peaks détectés")
  }
  
  expected_mean <- snapshot$data$peak_stats$intensity_mean
  actual_mean <- mean(actual$peaks$stain_intensity, na.rm = TRUE)
  
  # Tolérance: ±20% sur la moyenne (bruit CNN)
  relative_diff <- abs(actual_mean - expected_mean) / expected_mean
  expect_lt(relative_diff, 0.20,
            info = sprintf("Expected intensity_mean ~%.4f, got %.4f (%.1f%% diff)",
                           expected_mean, actual_mean, relative_diff * 100))
})

test_that("RÉGRESSION CNN: structure du résultat (peaks + boxes) reste identique", {
  actual <- .run_cnn_on_fixture()
  
  # La structure doit rester stable entre versions
  expect_true("peaks" %in% names(actual))
  expect_true("boxes" %in% names(actual))
  
  if (!is.null(actual$peaks) && nrow(actual$peaks) > 0) {
    expected_cols <- c("F2_ppm", "F1_ppm", "stain_intensity", "cluster_db", "stain_id")
    expect_true(all(expected_cols %in% names(actual$peaks)),
                info = paste("Missing columns:",
                             paste(setdiff(expected_cols, names(actual$peaks)),
                                   collapse = ", ")))
  }
  
  if (!is.null(actual$boxes) && nrow(actual$boxes) > 0) {
    expected_cols <- c("xmin", "xmax", "ymin", "ymax", "stain_id", "stain_intensity")
    expect_true(all(expected_cols %in% names(actual$boxes)))
  }
})

test_that("RÉGRESSION CNN: dimensions cohérentes des bounding boxes", {
  snapshot <- .load_cnn_snapshot()
  actual <- .run_cnn_on_fixture()
  
  if (is.null(actual$boxes) || nrow(actual$boxes) == 0) {
    skip("Pas de boxes détectées")
  }
  
  expected_width_f2 <- snapshot$data$box_stats$mean_width_f2
  actual_width_f2 <- mean(actual$boxes$xmax - actual$boxes$xmin, na.rm = TRUE)
  
  expected_width_f1 <- snapshot$data$box_stats$mean_width_f1
  actual_width_f1 <- mean(actual$boxes$ymax - actual$boxes$ymin, na.rm = TRUE)
  
  # Tolérance: ±30% sur la largeur moyenne (les boxes peuvent varier selon les clusters)
  expect_lt(abs(actual_width_f2 - expected_width_f2) / expected_width_f2, 0.30)
  expect_lt(abs(actual_width_f1 - expected_width_f1) / expected_width_f1, 0.30)
})
