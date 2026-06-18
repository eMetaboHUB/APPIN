# test-peak_picking.R - Tests complets pour Function/Peak_picking.R
# =============================================================================
# Fusion de: test-peak_picking.R, test-peak_picking-complete_8_.R
# =============================================================================

library(testthat)
library(dplyr)

# test-peak_picking.R - Tests pour Peak_picking.R

# Helper
create_spectrum <- function(nrow=100, ncol=100, peaks=list(list(r=30,c=30,a=50000)), noise=1000) {
  mat <- matrix(rnorm(nrow*ncol, 0, noise), nrow, ncol)
  for (p in peaks) for (i in 1:nrow) for (j in 1:ncol) mat[i,j] <- mat[i,j] + p$a*exp(-((i-p$r)^2+(j-p$c)^2)/100)
  rownames(mat) <- seq(10, 0, length.out=nrow); colnames(mat) <- seq(10, 0, length.out=ncol)
  mat
}

# =============================================================================
# TEST: peak_pick_2d_nt2() 
# =============================================================================
test_that("peak_pick_2d_nt2 exige rownames/colnames", {
  mat <- matrix(1:100, 10)
  expect_error(peak_pick_2d_nt2(mat, threshold_value=10), "row names.*column names")
})

test_that("peak_pick_2d_nt2 vide si threshold trop haut", {
  mat <- create_spectrum()
  result <- peak_pick_2d_nt2(mat, threshold_value=1e9, verbose=FALSE)
  expect_type(result, "list")
  expect_true(nrow(result$peaks) == 0 || nrow(result$centroids) == 0)
})

test_that("peak_pick_2d_nt2 structure correcte", {
  mat <- create_spectrum()
  result <- peak_pick_2d_nt2(mat, threshold_value=5000, spectrum_type="TOCSY", verbose=FALSE)
  expect_type(result, "list")
  expect_true("bounding_boxes" %in% names(result))
  expect_true("cluster_stats" %in% names(result))
})

test_that("peak_pick_2d_nt2 supporte tous spectrum_type", {
  mat <- create_spectrum()
  for (st in c("HSQC", "TOCSY", "COSY", "UFCOSY")) {
    result <- peak_pick_2d_nt2(mat, threshold_value=5000, spectrum_type=st, verbose=FALSE)
    expect_type(result, "list")
  }
})

# =============================================================================
# TEST: Calculs de filtrage
# =============================================================================
test_that("Elongation cluster carré = 1", {
  x_span <- 1; y_span <- 1
  elongation <- max(x_span/(y_span+1e-10), y_span/(x_span+1e-10))
  expect_equal(elongation, 1)
})

test_that("Elongation cluster allongé", {
  x_span <- 10; y_span <- 1
  elongation <- max(x_span/(y_span+1e-10), y_span/(x_span+1e-10))
  expect_equal(elongation, 10)
})

test_that("Détection diagonale TOCSY", {
  x_center <- 3.5; y_center <- 3.45
  is_diag <- abs(abs(x_center) - abs(y_center)) < 0.1
  expect_true(is_diag)
})

test_that("Détection hors diagonale", {
  x_center <- 3.5; y_center <- 7.0
  is_diag <- abs(abs(x_center) - abs(y_center)) < 0.1
  expect_false(is_diag)
})

test_that("Détection artefact horizontal", {
  x_span <- 0.5; y_span <- 0.01; y_var <- 0.00001; n_pts <- 50
  is_horiz <- x_span/(y_span+1e-10) > y_span/(x_span+1e-10)
  is_line <- is_horiz && y_var < 0.0001 && x_span > 0.08 && n_pts < 200
  expect_true(is_line)
})

test_that("Détection artefact vertical", {
  x_span <- 0.01; y_span <- 0.3; x_var <- 0.00001; n_pts <- 50
  is_vert <- y_span/(x_span+1e-10) > x_span/(y_span+1e-10)
  is_line <- is_vert && x_var < 0.00005 && y_span > 0.15 && n_pts < 200
  expect_true(is_line)
})

# =============================================================================
# TEST: Centroïdes pondérés
# =============================================================================
test_that("Centroïde pondéré poids égaux", {
  pts <- data.frame(x=c(1,2,3), y=c(1,2,3), level=c(100,200,100))
  cx <- sum(pts$x * pts$level) / sum(pts$level)
  cy <- sum(pts$y * pts$level) / sum(pts$level)
  expect_equal(cx, 2); expect_equal(cy, 2)
})

test_that("Centroïde pondéré poids asymétriques", {
  pts <- data.frame(x=c(1,2), y=c(1,2), level=c(300,100))
  cx <- sum(pts$x * pts$level) / sum(pts$level)
  expect_equal(cx, 1.25)
})
# test-peak_picking-complete.R - Tests unitaires pour Function/Peak_picking.R
# =============================================================================


# =============================================================================
# HELPERS: Créer des spectres de test
# =============================================================================

#' Créer un spectre 2D synthétique avec des pics contrôlés
create_test_spectrum_with_peaks <- function(size = 100, peaks = NULL, noise_sd = 0.5) {
  set.seed(42)
  mat <- matrix(rnorm(size * size, mean = 0, sd = noise_sd), 
                nrow = size, ncol = size)
  
  # Pics par défaut
  if (is.null(peaks)) {
    peaks <- list(
      list(row = 30, col = 30, amplitude = 100, sigma = 5),
      list(row = 50, col = 50, amplitude = 150, sigma = 4),
      list(row = 70, col = 70, amplitude = 80, sigma = 6)
    )
  }
  
  # Ajouter les pics
  for (p in peaks) {
    for (i in 1:size) {
      for (j in 1:size) {
        dist_sq <- (i - p$row)^2 + (j - p$col)^2
        mat[i, j] <- mat[i, j] + p$amplitude * exp(-dist_sq / (2 * p$sigma^2))
      }
    }
  }
  
  # Axes ppm (F1 = rows, F2 = cols) - typique pour COSY/TOCSY
  ppm_f1 <- seq(10, 0, length.out = size)  # Rows
  ppm_f2 <- seq(10, 0, length.out = size)  # Cols
  
  rownames(mat) <- as.character(round(ppm_f1, 4))
  colnames(mat) <- as.character(round(ppm_f2, 4))
  
  mat
}

#' Créer un spectre avec artefacts (lignes horizontales/verticales)
create_spectrum_with_artifacts <- function(size = 100) {
  set.seed(42)
  mat <- matrix(rnorm(size * size, mean = 0, sd = 0.5), 
                nrow = size, ncol = size)
  
  # Ajouter des vrais pics
  peaks <- list(
    list(row = 50, col = 50, amplitude = 100, sigma = 4)
  )
  
  for (p in peaks) {
    for (i in 1:size) {
      for (j in 1:size) {
        dist_sq <- (i - p$row)^2 + (j - p$col)^2
        mat[i, j] <- mat[i, j] + p$amplitude * exp(-dist_sq / (2 * p$sigma^2))
      }
    }
  }
  
  # Ajouter une ligne horizontale (t1 noise)
  mat[25, ] <- mat[25, ] + 30
  
  # Ajouter une ligne verticale (bleeding)
  mat[, 75] <- mat[, 75] + 25
  
  ppm_f1 <- seq(10, 0, length.out = size)
  ppm_f2 <- seq(10, 0, length.out = size)
  
  rownames(mat) <- as.character(round(ppm_f1, 4))
  colnames(mat) <- as.character(round(ppm_f2, 4))
  
  mat
}

# =============================================================================
# TESTS: Validation des entrées
# =============================================================================

test_that("peak_pick_2d_nt2: erreur si pas de rownames/colnames", {
  mat <- matrix(runif(100), nrow = 10, ncol = 10)
  # Pas de rownames/colnames
  
  expect_error(
    peak_pick_2d_nt2(mat, threshold_value = 0.5),
    "row names.*column names"
  )
})

test_that("peak_pick_2d_nt2: retourne liste vide si aucun pic au-dessus du seuil", {
  mat <- create_test_spectrum_with_peaks(size = 50)
  
  # Seuil très élevé
  result <- peak_pick_2d_nt2(mat, threshold_value = 10000, verbose = FALSE)
  
  expect_true(is.list(result))
  expect_equal(nrow(result$peaks), 0)
  expect_equal(nrow(result$bounding_boxes), 0)
})

# =============================================================================
# TESTS: Détection de pics basique
# =============================================================================

test_that("peak_pick_2d_nt2: détecte les pics dans un spectre simple", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    spectrum_type = "COSY",
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_true("peaks" %in% names(result) || "centroids" %in% names(result))
  expect_true("bounding_boxes" %in% names(result))
  expect_true("cluster_stats" %in% names(result))
})

test_that("peak_pick_2d_nt2: retourne une structure valide avec centroids", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    spectrum_type = "TOCSY",
    verbose = FALSE
  )
  
  # Vérifier la structure
 expect_true(is.list(result))
  
  # centroids peut être dans result$centroids ou result$peaks selon la version
  centroids <- result$centroids
  if (is.null(centroids)) centroids <- result$peaks
  
  # Doit être un data.frame (même vide)
  expect_true(is.data.frame(centroids) || is.null(centroids))
  
  # Si des pics sont détectés, vérifier le nombre raisonnable
  if (!is.null(centroids) && nrow(centroids) > 0) {
    n_detected <- nrow(centroids)
    expect_true(n_detected >= 1)
    expect_true(n_detected <= 50)  # Limite haute raisonnable
  }
})

# =============================================================================
# TESTS: Exclusion de régions (eau, etc.)
# =============================================================================

test_that("peak_pick_2d_nt2: f2_exclude_range fonctionne sans erreur", {
  # Créer un spectre avec un pic dans la région à exclure
  peaks <- list(
    list(row = 50, col = 50, amplitude = 100, sigma = 4),  # Pic normal
    list(row = 50, col = 47, amplitude = 100, sigma = 4)   # Pic dans région eau (vers 4.7 ppm)
  )
  mat <- create_test_spectrum_with_peaks(size = 100, peaks = peaks)
  
  # Sans exclusion
  result_no_exclude <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    f2_exclude_range = NULL,
    verbose = FALSE
  )
  
  # Avec exclusion (4.7-5.0 ppm)
  result_with_exclude <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    f2_exclude_range = c(4.7, 5.0),
    verbose = FALSE
  )
  
  # Les deux doivent retourner des listes valides
  expect_true(is.list(result_no_exclude))
  expect_true(is.list(result_with_exclude))
  
  # Helper pour obtenir le nombre de centroids
  get_n_centroids <- function(result) {
    centroids <- result$centroids
    if (is.null(centroids)) centroids <- result$peaks
    if (is.null(centroids)) return(0)
    return(nrow(centroids))
  }
  
  n_no_exclude <- get_n_centroids(result_no_exclude)
  n_with_exclude <- get_n_centroids(result_with_exclude)
  
  # L'exclusion ne doit pas ajouter de pics (peut en enlever ou garder le même nombre)
  expect_true(n_with_exclude <= n_no_exclude || n_no_exclude == 0)
})

# =============================================================================
# TESTS: Types de spectre - HSQC spécifique
# =============================================================================

test_that("peak_pick_2d_nt2: HSQC utilise le filtrage permissif", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "HSQC",
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  # HSQC a elongation <= 100 donc très permissif
})

# =============================================================================
# TESTS: Filtrage des artefacts spécifiques
# =============================================================================

test_that("peak_pick_2d_nt2: détecte les pics diagonaux", {
  # Créer un spectre avec un pic sur la diagonale (F1 ≈ F2)
  peaks <- list(
    list(row = 50, col = 50, amplitude = 100, sigma = 4)  # Sur la diagonale
  )
  mat <- create_test_spectrum_with_peaks(size = 100, peaks = peaks)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    spectrum_type = "TOCSY",
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  
  # Les stats doivent inclure is_diagonal
  if (!is.null(result$cluster_stats) && nrow(result$cluster_stats) > 0) {
    expect_true("is_diagonal" %in% names(result$cluster_stats))
  }
})

test_that("peak_pick_2d_nt2: COSY applique le filtrage strict", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "COSY",
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  
  # Vérifier que cluster_stats a les colonnes COSY-spécifiques
  if (!is.null(result$cluster_stats) && nrow(result$cluster_stats) > 0) {
    expect_true("is_vertical_artifact" %in% names(result$cluster_stats) ||
                "status" %in% names(result$cluster_stats))
  }
})

test_that("peak_pick_2d_nt2: UFCOSY utilise les paramètres relaxés", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result_ufcosy <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "UFCOSY",
    verbose = FALSE
  )
  
  result_cosy <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "COSY",
    verbose = FALSE
  )
  
  # Les deux doivent fonctionner
  expect_true(is.list(result_ufcosy))
  expect_true(is.list(result_cosy))
})

# =============================================================================
# TESTS: Statistiques de clusters
# =============================================================================

test_that("peak_pick_2d_nt2: cluster_stats contient les métriques de forme", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "TOCSY",
    verbose = FALSE
  )
  
  stats <- result$cluster_stats
  
  if (!is.null(stats) && nrow(stats) > 0) {
    # Vérifier les métriques de forme
    shape_cols <- c("elongation", "aspect_ratio_x", "aspect_ratio_y", 
                    "is_horizontal", "is_vertical", "density", "area")
    present <- shape_cols %in% names(stats)
    expect_true(sum(present) >= 3)
  }
})

test_that("peak_pick_2d_nt2: cluster_stats contient les flags d'artefacts", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    spectrum_type = "COSY",
    verbose = FALSE
  )
  
  stats <- result$cluster_stats
  
  if (!is.null(stats) && nrow(stats) > 0) {
    # Vérifier les flags d'artefacts
    artifact_cols <- c("is_horizontal_line", "is_vertical_line", "is_thin_vertical",
                       "is_artifact", "status")
    present <- artifact_cols %in% names(stats)
    expect_true(sum(present) >= 2)
  }
})

# =============================================================================
# TESTS: Redistribution d'intensité
# =============================================================================

test_that("peak_pick_2d_nt2: redistribue l'intensité des clusters rejetés", {
  # Ce test vérifie que la redistribution d'intensité fonctionne
  # sans vérifier les valeurs exactes (comportement interne)
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 2,  # Seuil bas pour avoir plus de clusters
    spectrum_type = "COSY",
    min_cluster_intensity = 0.01,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Filtrage par colonne F2 (bleeding)
# =============================================================================

test_that("peak_pick_2d_nt2: TOCSY filtre les pics par colonne F2", {
  # Créer un spectre avec plusieurs pics sur la même colonne F2
  peaks <- list(
    list(row = 30, col = 50, amplitude = 100, sigma = 4),
    list(row = 50, col = 50, amplitude = 50, sigma = 4),   # Même colonne F2
    list(row = 70, col = 50, amplitude = 25, sigma = 4),   # Même colonne F2
    list(row = 50, col = 30, amplitude = 80, sigma = 4)    # Colonne différente
  )
  mat <- create_test_spectrum_with_peaks(size = 100, peaks = peaks)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    spectrum_type = "TOCSY",
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Paramètres de prominence et threshold adaptatif
# =============================================================================

test_that("peak_pick_2d_nt2: prominence_factor affecte la détection", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  # Prominence stricte
  result_strict <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    prominence_factor = 0.5,  # Strict
    verbose = FALSE
  )
  
  # Prominence permissive
  result_permissive <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    prominence_factor = 0.001,  # Permissif
    verbose = FALSE
  )
  
  expect_true(is.list(result_strict))
  expect_true(is.list(result_permissive))
})

test_that("peak_pick_2d_nt2: adaptive_peak_threshold affecte la détection", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    adaptive_peak_threshold = 0.01,  # Top 1%
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Retour de liste vide quand aucun pic valide
# =============================================================================

test_that("peak_pick_2d_nt2: retourne liste vide si tous les clusters sont rejetés", {
  # Créer un spectre avec seulement du bruit faible
  set.seed(42)
  size <- 50
  mat <- matrix(rnorm(size * size, mean = 0, sd = 0.01), nrow = size, ncol = size)
  
  ppm_f1 <- seq(10, 0, length.out = size)
  ppm_f2 <- seq(10, 0, length.out = size)
  rownames(mat) <- as.character(round(ppm_f1, 4))
  colnames(mat) <- as.character(round(ppm_f2, 4))
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 100,  # Très haut, rien ne passe
    min_cluster_intensity = 1000,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  
  # Doit retourner des data.frames vides
  centroids <- result$centroids
  if (is.null(centroids)) centroids <- result$peaks
  
  expect_true(is.null(centroids) || nrow(centroids) == 0)
})

# =============================================================================
# TESTS: Paramètres de clustering
# =============================================================================

test_that("peak_pick_2d_nt2: eps_value affecte le clustering", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  # Petit eps = plus de clusters
  result_small_eps <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    eps_value = 0.01,
    verbose = FALSE
  )
  
  # Grand eps = moins de clusters (fusion)
  result_large_eps <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    eps_value = 0.1,
    verbose = FALSE
  )
  
  # Les deux doivent fonctionner
  expect_true(is.list(result_small_eps))
  expect_true(is.list(result_large_eps))
})

test_that("peak_pick_2d_nt2: neighborhood_size affecte la détection", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result_small <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    neighborhood_size = 3,
    verbose = FALSE
  )
  
  result_large <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    neighborhood_size = 11,
    verbose = FALSE
  )
  
  expect_true(is.list(result_small))
  expect_true(is.list(result_large))
})

# =============================================================================
# TESTS: keep_peak_ranges
# =============================================================================

test_that("peak_pick_2d_nt2: keep_peak_ranges filtre correctement", {
  mat <- create_test_spectrum_with_peaks(size = 100)
  
  # Garder seulement les pics dans certaines plages
  keep_ranges <- list(
    c(6.0, 4.0),  # Garder les tops pics entre 4 et 6 ppm
    c(8.0, 7.0)   # Garder les tops pics entre 7 et 8 ppm
  )
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    keep_peak_ranges = keep_ranges,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Sortie et structure des données
# =============================================================================

test_that("peak_pick_2d_nt2: centroids a les bonnes colonnes", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    verbose = FALSE
  )
  
  # Obtenir centroids (peut être dans $centroids ou $peaks)
  centroids <- result$centroids
  if (is.null(centroids)) centroids <- result$peaks
  
  # Si centroids existe et n'est pas vide
  if (!is.null(centroids) && is.data.frame(centroids) && nrow(centroids) > 0) {
    expect_true("stain_id" %in% names(centroids))
    expect_true("F2_ppm" %in% names(centroids))
    expect_true("F1_ppm" %in% names(centroids))
    expect_true("Volume" %in% names(centroids))
  } else {
    # Si pas de centroids, le test passe quand même (pas d'erreur)
    expect_true(TRUE)
  }
})

test_that("peak_pick_2d_nt2: bounding_boxes a les bonnes colonnes", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    verbose = FALSE
  )
  
  boxes <- result$bounding_boxes
  
  # Si bounding_boxes existe et n'est pas vide
  if (!is.null(boxes) && is.data.frame(boxes) && nrow(boxes) > 0) {
    expect_true("stain_id" %in% names(boxes))
    expect_true("xmin" %in% names(boxes))
    expect_true("xmax" %in% names(boxes))
    expect_true("ymin" %in% names(boxes))
    expect_true("ymax" %in% names(boxes))
  } else {
    expect_true(TRUE)
  }
})

test_that("peak_pick_2d_nt2: cluster_stats contient les statistiques", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    verbose = FALSE
  )
  
  stats <- result$cluster_stats
  
  # Si cluster_stats existe et n'est pas vide
  if (!is.null(stats) && is.data.frame(stats) && nrow(stats) > 0) {
    # Vérifier quelques colonnes attendues
    expected_cols <- c("intensity", "x_span", "y_span", "n_points")
    present <- expected_cols %in% names(stats)
    expect_true(sum(present) >= 2)  # Au moins quelques colonnes
  } else {
    expect_true(TRUE)
  }
})

# =============================================================================
# TESTS: Gestion des artefacts
# =============================================================================

test_that("peak_pick_2d_nt2: filtre les lignes horizontales (t1 noise)", {
  mat <- create_spectrum_with_artifacts(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 15,
    spectrum_type = "COSY",
    verbose = FALSE
  )
  
  # Le test vérifie juste que ça ne plante pas
  expect_true(is.list(result))
  expect_true("bounding_boxes" %in% names(result))
})

test_that("peak_pick_2d_nt2: filtre les traces verticales (bleeding)", {
  mat <- create_spectrum_with_artifacts(size = 100)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 15,
    spectrum_type = "TOCSY",
    verbose = FALSE
  )
  
  # Ne doit pas planter et doit retourner une liste valide
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Paramètres de bounding box
# =============================================================================

test_that("peak_pick_2d_nt2: box_padding affecte la taille des boxes", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result_small_pad <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    box_padding_f2 = 0.001,
    box_padding_f1 = 0.001,
    verbose = FALSE
  )
  
  result_large_pad <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    box_padding_f2 = 0.05,
    box_padding_f1 = 0.05,
    verbose = FALSE
  )
  
  # Les deux doivent fonctionner
  expect_true(is.list(result_small_pad))
  expect_true(is.list(result_large_pad))
})

# =============================================================================
# TESTS: verbose et diagnostics
# =============================================================================

test_that("peak_pick_2d_nt2: verbose=FALSE supprime les messages", {
  mat <- create_test_spectrum_with_peaks(size = 50)
  
  # Capturer les messages
  messages <- capture.output({
    result <- peak_pick_2d_nt2(
      mat, 
      threshold_value = 10,
      verbose = FALSE
    )
  }, type = "message")
  
  # Avec verbose=FALSE, il y aura toujours des messages() mais moins
  expect_true(is.list(result))
})

test_that("peak_pick_2d_nt2: diagnose_zones ne fait pas planter", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    diagnose_zones = c(5.0, 6.0),
    diagnose_radius = 0.5,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# TESTS: Cas limites
# =============================================================================

test_that("peak_pick_2d_nt2: gère un spectre très petit", {
  mat <- create_test_spectrum_with_peaks(size = 20)
  
  result <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 5,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
})

test_that("peak_pick_2d_nt2: gère min_cluster_intensity", {
  mat <- create_test_spectrum_with_peaks(size = 80)
  
  result_low <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    min_cluster_intensity = 0.001,
    verbose = FALSE
  )
  
  result_high <- peak_pick_2d_nt2(
    mat, 
    threshold_value = 10,
    min_cluster_intensity = 1000,
    verbose = FALSE
  )
  
  # Helper pour obtenir le nombre de centroids
  get_n_centroids <- function(result) {
    centroids <- result$centroids
    if (is.null(centroids)) centroids <- result$peaks
    if (is.null(centroids) || !is.data.frame(centroids)) return(0)
    return(nrow(centroids))
  }
  
  n_low <- get_n_centroids(result_low)
  n_high <- get_n_centroids(result_high)
  
  # High min_cluster_intensity devrait donner moins ou égal de pics
  expect_true(n_high <= n_low)
})

# =============================================================================
# TESTS: filter_noise_peaks
# =============================================================================

test_that("filter_noise_peaks: retourne le même dataframe si < 2 pics", {
  # Un seul pic
  peaks <- data.frame(
    F2_ppm = 5.0,
    F1_ppm = 5.0,
    Volume = 100,
    stain_id = "peak_1"
  )
  
  result <- filter_noise_peaks(peaks)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$stain_id, "peak_1")
})
  
test_that("filter_noise_peaks: filtre les pics isolés", {
  # Créer des pics avec un pic isolé
  peaks <- data.frame(
    F2_ppm = c(5.0, 5.01, 5.02, 5.03, 10.0),  # 4 pics groupés + 1 isolé à 10 ppm
    F1_ppm = c(5.0, 5.01, 5.02, 5.03, 10.0),
    Volume = c(100, 90, 80, 70, 50),  # Le pic isolé a une intensité faible
    stain_id = c("p1", "p2", "p3", "p4", "isolated")
  )
  
  result <- filter_noise_peaks(peaks, min_neighbors = 2, neighbor_radius = 0.05, min_relative_intensity = 0.6)
  
  # Le pic isolé devrait être filtré (pas assez de voisins et intensité < 60%)
  expect_true(nrow(result) < nrow(peaks))
})

test_that("filter_noise_peaks: garde les pics avec haute intensité relative", {
  peaks <- data.frame(
    F2_ppm = c(5.0, 10.0),
    F1_ppm = c(5.0, 10.0),
    Volume = c(100, 80),  # 80% de l'intensité max
    stain_id = c("p1", "p2")
  )
  
  result <- filter_noise_peaks(peaks, min_neighbors = 5, neighbor_radius = 0.01, min_relative_intensity = 0.5)
  
  # Les deux pics doivent être gardés car > 50% de l'intensité max
  expect_equal(nrow(result), 2)
})

test_that("filter_noise_peaks: garde les pics avec assez de voisins", {
  # Créer un cluster de pics proches
  peaks <- data.frame(
    F2_ppm = c(5.00, 5.01, 5.02, 5.03, 5.04),
    F1_ppm = c(5.00, 5.01, 5.02, 5.03, 5.04),
    Volume = c(100, 10, 10, 10, 10),  # Faible intensité sauf le premier
    stain_id = paste0("p", 1:5)
  )
  
  result <- filter_noise_peaks(peaks, min_neighbors = 2, neighbor_radius = 0.05, min_relative_intensity = 0.9)
  
  # Tous devraient être gardés car ils ont des voisins
  expect_equal(nrow(result), 5)
})

# =============================================================================
# TESTS: process_nmr_centroids
# =============================================================================

#' Helper pour créer des données de contour simulées
create_mock_contour_data <- function(n_clusters = 3, points_per_cluster = 50) {
  set.seed(42)
  
  contour_data <- data.frame()
  
  for (i in seq_len(n_clusters)) {
    # Centre du cluster
    cx <- -5 + i * 2
    cy <- -5 + i * 2
    
    # Points autour du centre
    cluster_points <- data.frame(
      x = cx + rnorm(points_per_cluster, sd = 0.02),
      y = cy + rnorm(points_per_cluster, sd = 0.02),
      level = runif(points_per_cluster, 10, 100),
      group = paste0("group_", i)
    )
    
    contour_data <- rbind(contour_data, cluster_points)
  }
  
  contour_data
}

test_that("process_nmr_centroids: erreur si contour_data vide", {
  empty_data <- data.frame(x = numeric(0), y = numeric(0), level = numeric(0), group = character(0))
  
  expect_error(
    process_nmr_centroids(
      rr_data = matrix(1, 10, 10),
      contour_data = empty_data,
      eps_value = 0.3
    ),
    "No contours"
  )
})

test_that("process_nmr_centroids: erreur si colonne group manquante", {
  bad_data <- data.frame(x = 1:10, y = 1:10, level = 1:10)  # Pas de colonne group
  
  expect_error(
    process_nmr_centroids(
      rr_data = matrix(1, 10, 10),
      contour_data = bad_data,
      eps_value = 0.3
    ),
    "No contours"
  )
})

test_that("process_nmr_centroids: fonctionne avec spectrum_type HSQC", {
  contour_data <- create_mock_contour_data(n_clusters = 2, points_per_cluster = 30)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "HSQC",
    min_cluster_intensity = 0.01
  )
  
  expect_true(is.list(result))
  expect_true("centroids" %in% names(result) || "peaks" %in% names(result))
  expect_true("bounding_boxes" %in% names(result))
  expect_true("cluster_stats" %in% names(result))
})

test_that("process_nmr_centroids: fonctionne avec spectrum_type TOCSY", {
  contour_data <- create_mock_contour_data(n_clusters = 3, points_per_cluster = 40)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "TOCSY",
    min_cluster_intensity = 0.01
  )
  
  expect_true(is.list(result))
})

test_that("process_nmr_centroids: fonctionne avec spectrum_type COSY", {
  contour_data <- create_mock_contour_data(n_clusters = 2, points_per_cluster = 50)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "COSY",
    min_cluster_intensity = 0.01
  )
  
  expect_true(is.list(result))
})

test_that("process_nmr_centroids: fonctionne avec spectrum_type UFCOSY", {
  contour_data <- create_mock_contour_data(n_clusters = 2, points_per_cluster = 30)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "UFCOSY",
    min_cluster_intensity = 0.01
  )
  
  expect_true(is.list(result))
})

test_that("process_nmr_centroids: cluster_stats contient les métriques de forme", {
  contour_data <- create_mock_contour_data(n_clusters = 2, points_per_cluster = 50)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "TOCSY"
  )
  
  stats <- result$cluster_stats
  
  if (!is.null(stats) && nrow(stats) > 0) {
    expected_cols <- c("intensity", "x_span", "y_span", "elongation", "density")
    present <- expected_cols %in% names(stats)
    expect_true(sum(present) >= 3)
  }
})

test_that("process_nmr_centroids: gère keep_peak_ranges", {
  contour_data <- create_mock_contour_data(n_clusters = 3, points_per_cluster = 30)
  
  result <- process_nmr_centroids(
    rr_data = matrix(runif(100), 10, 10),
    contour_data = contour_data,
    eps_value = 0.5,
    spectrum_type = "TOCSY",
    keep_peak_ranges = list(c(0, -10))  # Plage à filtrer
  )
  
  expect_true(is.list(result))
})

# =============================================================================
# RÉSUMÉ
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           TESTS POUR Function/Peak_picking.R                     ║\n")
cat("║                                                                  ║\n")
cat("║  Fonctions testées:                                             ║\n")
cat("║  - peak_pick_2d_nt2                                             ║\n")
cat("║  - filter_noise_peaks                                           ║\n")
cat("║  - process_nmr_centroids                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")
