# test-vizualisation.R - Tests pour Vizualisation.R
library(testthat)

# =============================================================================
# TEST: seuil_bruit_multiplicatif()
# =============================================================================
test_that("seuil_bruit_multiplicatif calcule facteur * SD", {
  set.seed(42)
  mat <- matrix(rnorm(1000, 0, 10), 100)
  seuil <- seuil_bruit_multiplicatif(mat, 3)
  expect_true(seuil > 25 && seuil < 35)
})

test_that("seuil_bruit_multiplicatif gère NA", {
  mat <- matrix(c(1,2,NA,4,5,NA,7,8,9), 3)
  expect_true(!is.na(seuil_bruit_multiplicatif(mat, 3)))
})

# =============================================================================
# TEST: seuil_max_pourcentage()
# =============================================================================
test_that("seuil_max_pourcentage calcule % du max", {
  mat <- matrix(c(0, 50, 100, 150, 200), 1)
  expect_equal(seuil_max_pourcentage(mat, 0.05), 10)
  expect_equal(seuil_max_pourcentage(mat, 0.10), 20)
})

test_that("seuil_max_pourcentage gère NA", {
  mat <- matrix(c(10, NA, 30, NA, 50), 1)
  expect_equal(seuil_max_pourcentage(mat, 0.1), 5)
})

# =============================================================================
# TEST: modulate_threshold()
# =============================================================================
test_that("modulate_threshold(0) = 0", {
  expect_equal(modulate_threshold(0), 0)
})

test_that("modulate_threshold croissante", {
  expect_true(modulate_threshold(1000) < modulate_threshold(10000))
  expect_true(modulate_threshold(10000) < modulate_threshold(100000))
})

test_that("modulate_threshold formule a * VI^b", {
  VI <- 50000
  expect_equal(modulate_threshold(VI), 0.0006 * VI^1.2)
})

# =============================================================================
# TEST: make_bbox_outline()
# =============================================================================
test_that("make_bbox_outline rectangle fermé", {
  boxes <- data.frame(stain_id="p1", xmin=1, xmax=2, ymin=3, ymax=4)
  outline <- make_bbox_outline(boxes)
  expect_s3_class(outline, "data.frame")
  expect_true(all(c("x","y","group") %in% names(outline)))
  expect_equal(nrow(outline), 5)
})

test_that("make_bbox_outline plusieurs boxes", {
  boxes <- data.frame(stain_id=c("p1","p2"), xmin=c(1,5), xmax=c(2,6), ymin=c(1,5), ymax=c(2,6))
  outline <- make_bbox_outline(boxes)
  expect_true(nrow(outline) >= 10)
})

test_that("make_bbox_outline NULL si vide", {
  expect_null(make_bbox_outline(NULL))
  expect_null(make_bbox_outline(data.frame()))
})

# =============================================================================
# TEST: get_local_Volume()
# =============================================================================
test_that("get_local_Volume somme niveaux", {
  cd <- data.frame(x=c(1,1.001,1.002,2), y=c(3,3.001,3.002,4), level=c(100,150,200,50))
  vol <- get_local_Volume(1.001, 3.001, cd, 0.01)
  expect_equal(vol, 450)
})

test_that("get_local_Volume NA si aucun point", {
  cd <- data.frame(x=1, y=3, level=100)
  expect_true(is.na(get_local_Volume(10, 10, cd, 0.01)))
})

# =============================================================================
# TEST: find_nmr_peak_centroids_optimized()
# =============================================================================
test_that("find_nmr_peak_centroids_optimized valide input", {
  expect_error(find_nmr_peak_centroids_optimized(NULL), "Invalid")
  expect_error(find_nmr_peak_centroids_optimized("x"), "Invalid")
})

test_that("find_nmr_peak_centroids_optimized valide spectrum_type", {
  mat <- matrix(rnorm(100), 10); rownames(mat) <- 1:10; colnames(mat) <- 1:10
  # Un spectrum_type inconnu ne leve plus d'erreur : la fonction emet un
  # warning et bascule sur des parametres generiques par defaut.
  expect_warning(
    find_nmr_peak_centroids_optimized(mat, "INVALID"),
    regexp = "Unknown spectrum_type"
  )
})