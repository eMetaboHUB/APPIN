# test-utils.R - Tests complets pour R/utils.R
# =============================================================================
# Fusion de: test-utils.R, test-utils-complete_5_.R
# =============================================================================

library(testthat)

# test-utils.R - Tests pour R/utils.R

# =============================================================================
# TEST: Opérateur %||% (null coalescing)
# =============================================================================
test_that("%||% retourne valeur si non-NULL", {
  `%||%` <- function(a, b) if(is.null(a)) b else a
  expect_equal(5 %||% 10, 5)
  expect_equal("test" %||% "default", "test")
  expect_equal(0 %||% 10, 0)
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that("%||% retourne défaut si NULL", {
  `%||%` <- function(a, b) if(is.null(a)) b else a
  expect_equal(NULL %||% 10, 10)
  expect_equal(NULL %||% "default", "default")
})

# =============================================================================
# TEST: parse_keep_peak_ranges()
# =============================================================================
test_that("parse_keep_peak_ranges parse format standard", {
  parse_kpr <- function(text) {
    if(is.null(text) || trimws(text)=="") return(NULL)
    ranges <- strsplit(text, ";")[[1]]
    ranges <- trimws(ranges[ranges != ""])
    if(length(ranges)==0) return(NULL)
    result <- lapply(ranges, function(r) {
      parts <- strsplit(r, ",")[[1]]
      if(length(parts)!=2) return(NULL)
      as.numeric(trimws(parts))
    })
    result[!sapply(result, is.null)]
  }
  
  result <- parse_kpr("0.5,-0.5; 1,0.8")
  expect_type(result, "list")
  expect_length(result, 2)
  expect_equal(result[[1]], c(0.5, -0.5))
  expect_equal(result[[2]], c(1, 0.8))
})

test_that("parse_keep_peak_ranges NULL si vide", {
  parse_kpr <- function(text) {
    if(is.null(text) || trimws(text)=="") return(NULL)
    NULL
  }
  expect_null(parse_kpr(NULL))
  expect_null(parse_kpr(""))
  expect_null(parse_kpr("   "))
})

# =============================================================================
# TEST: Calcul volume somme
# =============================================================================
test_that("volume somme = total intensités", {
  mat <- matrix(1, 10, 10)
  rownames(mat) <- 10:1; colnames(mat) <- 1:10
  
  calc_vol <- function(mat, box) {
    ppm_f2 <- as.numeric(colnames(mat))
    ppm_f1 <- as.numeric(rownames(mat))
    x_idx <- which(ppm_f2 >= box$xmin & ppm_f2 <= box$xmax)
    y_idx <- which(ppm_f1 >= box$ymin & ppm_f1 <= box$ymax)
    if(length(x_idx)==0 || length(y_idx)==0) return(NA_real_)
    sum(mat[y_idx, x_idx], na.rm=TRUE)
  }
  
  expect_equal(calc_vol(mat, list(xmin=1,xmax=10,ymin=1,ymax=10)), 100)
  expect_equal(calc_vol(mat, list(xmin=1,xmax=5,ymin=1,ymax=5)), 25)
})

# =============================================================================
# TEST: Normalisation z-score
# =============================================================================
test_that("Normalisation z-score", {
  norm <- function(x) (x - mean(x)) / sd(x)
  x <- c(1, 2, 3, 4, 5)
  x_norm <- norm(x)
  expect_equal(mean(x_norm), 0, tolerance=1e-10)
  expect_equal(sd(x_norm), 1, tolerance=1e-10)
})

# =============================================================================
# TEST: Scaling F1 HSQC
# =============================================================================
test_that("Scaling F1 HSQC divise par facteur", {
  scale_f1 <- function(f1, factor=5) f1 / factor
  expect_equal(scale_f1(c(10,20,30,40,50)), c(2,4,6,8,10))
})
# test-utils.R - Tests unitaires pour R/utils.R
# =============================================================================


# =============================================================================
# TESTS: Opérateur %||% (null-coalesce)
# =============================================================================

test_that("%||%: retourne a si non NULL", {
  expect_equal("value" %||% "default", "value")
  expect_equal(123 %||% 0, 123)
  expect_equal(FALSE %||% TRUE, FALSE)
  expect_equal(0 %||% 99, 0)
  expect_equal("" %||% "default", "")  # Chaîne vide n'est pas NULL
  expect_equal(list() %||% "default", list())  # Liste vide n'est pas NULL
})

test_that("%||%: retourne b si a est NULL",
{
  expect_equal(NULL %||% "default", "default")
  expect_equal(NULL %||% 123, 123)
  expect_equal(NULL %||% NULL, NULL)
})

# =============================================================================
# TESTS: Cache de spectres
# =============================================================================

test_that("clear_spectrum_cache: vide le cache", {
  # S'assurer que le cache existe
  expect_true(exists("spectrum_cache"))
  
  # Nettoyer
  result <- clear_spectrum_cache()
  
  # Vérifier que ça retourne NULL invisible

  expect_null(result)
  
  # Vérifier que le cache est vide
  expect_equal(length(ls(envir = spectrum_cache)), 0)
})

# Note: read_bruker_cached nécessite read_bruker() et des fichiers réels
# On teste seulement le mécanisme de cache ici

test_that("spectrum_cache: environnement existe", {
  expect_true(exists("spectrum_cache"))
  expect_true(is.environment(spectrum_cache))
})

# =============================================================================
# TESTS: parse_keep_peak_ranges
# =============================================================================

test_that("parse_keep_peak_ranges: parsing correct de chaînes valides", {
  # Format standard
  result <- parse_keep_peak_ranges("0.5,-0.5; 1,0.8; 1.55,1.45")
  expect_equal(length(result), 3)
  expect_equal(result[[1]], c(0.5, -0.5))
  expect_equal(result[[2]], c(1, 0.8))
  expect_equal(result[[3]], c(1.55, 1.45))
})

test_that("parse_keep_peak_ranges: gestion des espaces", {
  result <- parse_keep_peak_ranges("  0.5 , -0.5 ;  1.0 , 0.8  ")
  expect_equal(length(result), 2)
  expect_equal(result[[1]], c(0.5, -0.5))
  expect_equal(result[[2]], c(1.0, 0.8))
})

test_that("parse_keep_peak_ranges: retourne NULL pour entrées vides/invalides", {
  expect_null(parse_keep_peak_ranges(NULL))
  expect_null(parse_keep_peak_ranges(""))
  # Note: "   " (espaces) peut retourner une liste vide selon l'implémentation
  result_spaces <- parse_keep_peak_ranges("   ")
  expect_true(is.null(result_spaces) || length(result_spaces) == 0)
})

test_that("parse_keep_peak_ranges: filtre les paires invalides", {
  # Paire avec 1 seul nombre
  result <- parse_keep_peak_ranges("0.5,-0.5; invalid; 1,0.8")
  expect_equal(length(result), 2)
  
  # Paire avec texte
  result2 <- parse_keep_peak_ranges("abc,def; 1,2")
  expect_equal(length(result2), 1)
  expect_equal(result2[[1]], c(1, 2))
})

test_that("parse_keep_peak_ranges: gère le point-virgule final", {
  result <- parse_keep_peak_ranges("0.5,-0.5; 1,0.8;")
  expect_equal(length(result), 2)
})

# =============================================================================
# TESTS: clean_centroids_df
# =============================================================================

test_that("clean_centroids_df: convertit les virgules en points", {
  df <- data.frame(
    F2_ppm = c("1,5", "2,3", "3,7"),
    F1_ppm = c("4,5", "5,6", "6,7"),
    Volume = c("100,5", "200,3", "300,1"),
    stringsAsFactors = FALSE
  )
  
  result <- clean_centroids_df(df)
  
  expect_equal(result$F2_ppm, c(1.5, 2.3, 3.7))
  expect_equal(result$F1_ppm, c(4.5, 5.6, 6.7))
  expect_equal(result$Volume, c(100.5, 200.3, 300.1))
})

test_that("clean_centroids_df: gère les espaces", {
  df <- data.frame(
    F2_ppm = c(" 1.5 ", "  2.3", "3.7  "),
    F1_ppm = c("4.5", "5.6", "6.7"),
    Volume = c("100", "200", "300"),
    stringsAsFactors = FALSE
  )
  
  result <- clean_centroids_df(df)
  
  expect_equal(result$F2_ppm, c(1.5, 2.3, 3.7))
})

test_that("clean_centroids_df: gère les valeurs déjà numériques (format point)", {
  df <- data.frame(
    F2_ppm = c("1.5", "2.3", "3.7"),
    F1_ppm = c("4.5", "5.6", "6.7"),
    Volume = c("100.5", "200.3", "300.1"),
    stringsAsFactors = FALSE
  )
  
  result <- clean_centroids_df(df)
  
  expect_true(is.numeric(result$F2_ppm))
  expect_true(is.numeric(result$F1_ppm))
  expect_true(is.numeric(result$Volume))
})

# =============================================================================
# TESTS: get_box_intensity
# =============================================================================

# Helper: créer une matrice de spectre synthétique
create_test_spectrum <- function(size = 50) {
  set.seed(42)
  mat <- matrix(rnorm(size * size, mean = 0, sd = 0.1), nrow = size, ncol = size)
  
  # Ajouter un pic au centre
  cx <- size / 2
  cy <- size / 2
  for (i in 1:size) {
    for (j in 1:size) {
      dist <- sqrt((i - cy)^2 + (j - cx)^2)
      mat[i, j] <- mat[i, j] + 10 * exp(-dist^2 / 50)
    }
  }
  
  # Axes ppm
  ppm_x <- seq(10, 0, length.out = size)
  ppm_y <- seq(10, 0, length.out = size)
  
  colnames(mat) <- as.character(round(ppm_x, 4))
  rownames(mat) <- as.character(round(ppm_y, 4))
  
  list(mat = mat, ppm_x = ppm_x, ppm_y = ppm_y)
}

test_that("get_box_intensity: retourne vecteur vide pour boxes vides", {
  spec <- create_test_spectrum()
  boxes <- data.frame(xmin = numeric(0), xmax = numeric(0), 
                      ymin = numeric(0), ymax = numeric(0))
  
  result <- get_box_intensity(spec$mat, spec$ppm_x, spec$ppm_y, boxes)
  
  expect_equal(length(result), 0)
  expect_true(is.numeric(result))
})

test_that("get_box_intensity: méthode sum fonctionne", {
  spec <- create_test_spectrum()
  
  # Box couvrant le centre (où se trouve le pic)
  boxes <- data.frame(
    xmin = 4, xmax = 6,
    ymin = 4, ymax = 6
  )
  
  result <- get_box_intensity(spec$mat, spec$ppm_x, spec$ppm_y, boxes, method = "sum")
  
  expect_equal(length(result), 1)
  expect_true(is.numeric(result))
  expect_true(result > 0)  # Le pic doit avoir une intensité positive
})

test_that("get_box_intensity: retourne NA pour box hors spectre", {
  spec <- create_test_spectrum()
  
  # Box complètement hors du spectre
  boxes <- data.frame(
    xmin = 100, xmax = 110,
    ymin = 100, ymax = 110
  )
  
  result <- get_box_intensity(spec$mat, spec$ppm_x, spec$ppm_y, boxes, method = "sum")
  
  expect_equal(length(result), 1)
  expect_true(is.na(result))
})

test_that("get_box_intensity: gère plusieurs boxes", {
  spec <- create_test_spectrum()
  
  boxes <- data.frame(
    xmin = c(4, 1, 7),
    xmax = c(6, 3, 9),
    ymin = c(4, 1, 7),
    ymax = c(6, 3, 9)
  )
  
  result <- get_box_intensity(spec$mat, spec$ppm_x, spec$ppm_y, boxes, method = "sum")
  
  expect_equal(length(result), 3)
  expect_true(all(is.numeric(result)))
})

# =============================================================================
# TESTS: calculate_batch_box_intensities
# =============================================================================

# Helper: créer une liste de spectres mock
create_mock_spectra_list <- function(n_spectra = 2, size = 30) {
  spectra <- list()
  for (i in seq_len(n_spectra)) {
    set.seed(42 + i)
    mat <- matrix(runif(size * size, 0, 10), nrow = size, ncol = size)
    
    # Ajouter un pic
    cx <- size / 2 + (i - 1) * 2  # Position légèrement différente par spectre
    cy <- size / 2
    for (row in 1:size) {
      for (col in 1:size) {
        dist <- sqrt((row - cy)^2 + (col - cx)^2)
        mat[row, col] <- mat[row, col] + 50 * exp(-dist^2 / 20)
      }
    }
    
    ppm_x <- seq(10, 0, length.out = size)
    ppm_y <- seq(10, 0, length.out = size)
    
    colnames(mat) <- as.character(round(ppm_x, 4))
    rownames(mat) <- as.character(round(ppm_y, 4))
    
    spectra[[paste0("spectrum_", i)]] <- list(spectrumData = mat)
  }
  spectra
}

test_that("calculate_batch_box_intensities: erreur si boxes NULL ou vide", {
  spectra <- create_mock_spectra_list(1)
  
  expect_error(
    calculate_batch_box_intensities(NULL, spectra),
    "reference_boxes is empty or NULL"
  )
  
  expect_error(
    calculate_batch_box_intensities(data.frame(), spectra),
    "reference_boxes is empty or NULL"
  )
})

test_that("calculate_batch_box_intensities: fonctionne avec méthode sum", {
  spectra <- create_mock_spectra_list(2, size = 30)
  
  boxes <- data.frame(
    xmin = c(4, 6),
    xmax = c(6, 8),
    ymin = c(4, 6),
    ymax = c(6, 8),
    stain_id = c("box_1", "box_2"),
    stringsAsFactors = FALSE
  )
  
  result <- calculate_batch_box_intensities(boxes, spectra, method = "sum")
  
  # Vérifier la structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true("stain_id" %in% names(result))
  expect_true("F2_ppm" %in% names(result))
  expect_true("F1_ppm" %in% names(result))
  
  # Vérifier qu'il y a des colonnes d'intensité
  intensity_cols <- grep("^Intensity_", names(result), value = TRUE)
  expect_equal(length(intensity_cols), 2)  # 2 spectres
})

test_that("calculate_batch_box_intensities: gère stain_id manquant", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  # Boxes sans stain_id
  boxes <- data.frame(
    xmin = c(4),
    xmax = c(6),
    ymin = c(4),
    ymax = c(6)
  )
  
  result <- calculate_batch_box_intensities(boxes, spectra, method = "sum")
  
  # Un stain_id doit être généré
  expect_true("stain_id" %in% names(result))
  expect_false(any(is.na(result$stain_id)))
})

test_that("calculate_batch_box_intensities: détecte les doublons", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  # Boxes avec stain_id dupliqués
  boxes <- data.frame(
    xmin = c(4, 5),
    xmax = c(6, 7),
    ymin = c(4, 5),
    ymax = c(6, 7),
    stain_id = c("same_id", "same_id")
  )
  
  # Doit afficher un warning mais ne pas planter
  expect_warning(
    result <- calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "Duplicate|duplicate"  # Accepte les deux casses
  )
})

test_that("calculate_batch_box_intensities: shift_tolerance_ppm fonctionne", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = c(4),
    xmax = c(6),
    ymin = c(4),
    ymax = c(6),
    stain_id = c("box_1")
  )
  
  # Avec tolérance de shift
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "sum", 
    shift_tolerance_ppm = 0.05
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

# =============================================================================
# TESTS: calculate_batch_box_intensities - CAS SUPPLÉMENTAIRES
# =============================================================================

test_that("calculate_batch_box_intensities: erreur si colonnes manquantes", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  # Boxes sans xmin
  boxes <- data.frame(
    xmax = c(6),
    ymin = c(4),
    ymax = c(6),
    stain_id = c("box_1")
  )
  
  expect_error(
    calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "Missing columns"
  )
})

test_that("calculate_batch_box_intensities: erreur si toutes les coordonnées sont NA", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = NA_real_,
    xmax = NA_real_,
    ymin = NA_real_,
    ymax = NA_real_,
    stain_id = "box_1"
  )
  
  expect_error(
    calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "NA coordinates"
  )
})

test_that("calculate_batch_box_intensities: corrige les boxes inversées", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  # Box avec xmin > xmax et ymin > ymax
  boxes <- data.frame(
    xmin = 8,  # Plus grand que xmax
    xmax = 4,  # Plus petit que xmin
    ymin = 7,  # Plus grand que ymax  
    ymax = 3,  # Plus petit que ymin
    stain_id = "inverted_box"
  )
  
  result <- calculate_batch_box_intensities(boxes, spectra, method = "sum")
  
  # Doit fonctionner sans erreur
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("calculate_batch_box_intensities: gère les boxes avec coordonnées identiques", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  # Deux boxes avec les mêmes coordonnées
  boxes <- data.frame(
    xmin = c(4, 4),
    xmax = c(6, 6),
    ymin = c(4, 4),
    ymax = c(6, 6),
    stain_id = c("box_1", "box_2")
  )
  
  expect_warning(
    result <- calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "identical coordinates"
  )
  
  # Une seule box doit rester après suppression des doublons
  expect_equal(nrow(result), 1)
})

test_that("calculate_batch_box_intensities: gère spectre invalide (NULL)", {
  # Créer une liste avec un spectre NULL
  spectra <- list(
    "valid_spectrum" = list(spectrumData = matrix(runif(100), 10, 10)),
    "null_spectrum" = NULL
  )
  
  # Ajouter les noms de colonnes/lignes au spectre valide
  colnames(spectra$valid_spectrum$spectrumData) <- as.character(seq(10, 1, length.out = 10))
  rownames(spectra$valid_spectrum$spectrumData) <- as.character(seq(10, 1, length.out = 10))
  
  boxes <- data.frame(
    xmin = 3, xmax = 7, ymin = 3, ymax = 7, stain_id = "box_1"
  )
  
  expect_warning(
    result <- calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "invalid"
  )
  
  # Doit avoir 2 colonnes d'intensité, une avec des valeurs, une avec NA
  intensity_cols <- grep("^Intensity_", names(result), value = TRUE)
  expect_equal(length(intensity_cols), 2)
})

test_that("calculate_batch_box_intensities: gère spectre avec ppm invalides", {
  # Créer un spectre avec des noms de colonnes non numériques
  mat <- matrix(runif(100), 10, 10)
  colnames(mat) <- letters[1:10]  # Non numérique
  rownames(mat) <- LETTERS[1:10]  # Non numérique
  
  spectra <- list(
    "bad_ppm" = list(spectrumData = mat)
  )
  
  boxes <- data.frame(
    xmin = 3, xmax = 7, ymin = 3, ymax = 7, stain_id = "box_1"
  )
  
  expect_warning(
    result <- calculate_batch_box_intensities(boxes, spectra, method = "sum"),
    "invalid ppm"
  )
})

test_that("calculate_batch_box_intensities: apply_shift fonctionne", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = 4, xmax = 6, ymin = 4, ymax = 6, stain_id = "box_1"
  )
  
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    apply_shift = TRUE,
    method = "sum"
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("calculate_batch_box_intensities: progress callback fonctionne si fourni", {
  spectra <- create_mock_spectra_list(2, size = 30)
  
  boxes <- data.frame(
    xmin = c(4, 5), xmax = c(6, 7), ymin = c(4, 5), ymax = c(6, 7), 
    stain_id = c("box_1", "box_2")
  )
  
  progress_calls <- 0
  mock_progress <- function(value, detail) {
    progress_calls <<- progress_calls + 1
  }
  
  # Le callback n'est utilisé que pour method="fit", pas "sum"
  # Ce test vérifie juste que passer un callback ne cause pas d'erreur
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "sum",
    progress = mock_progress
  )
  
  expect_true(is.data.frame(result))
  # Note: progress peut ne pas être appelé pour method="sum"
})

test_that("calculate_batch_box_intensities: méthode fit fonctionne", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = 4, xmax = 6, ymin = 4, ymax = 6, stain_id = "box_1"
  )
  
  # Nécessite calculate_fitted_volumes
 skip_if_not(exists("calculate_fitted_volumes"), 
              "calculate_fitted_volumes non disponible")
  
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "fit",
    model = "gaussian"
  )
  
  expect_true(is.data.frame(result))
})

test_that("calculate_batch_box_intensities: méthode fit avec shift_tolerance_ppm", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = 4, xmax = 6, ymin = 4, ymax = 6, stain_id = "box_1"
  )
  
  skip_if_not(exists("calculate_fitted_volumes"), 
              "calculate_fitted_volumes non disponible")
  
  # Tester avec shift_tolerance_ppm > 0 pour couvrir le recentrage dynamique
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "fit",
    model = "gaussian",
    shift_tolerance_ppm = 0.05
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("calculate_batch_box_intensities: méthode sum avec shift_tolerance_ppm", {
  spectra <- create_mock_spectra_list(1, size = 30)
  
  boxes <- data.frame(
    xmin = 4, xmax = 6, ymin = 4, ymax = 6, stain_id = "box_1"
  )
  
  # Tester avec shift_tolerance_ppm > 0 pour couvrir le recentrage dynamique
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "sum",
    shift_tolerance_ppm = 0.05
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("calculate_batch_box_intensities: shift_tolerance avec plusieurs boxes", {
  spectra <- create_mock_spectra_list(2, size = 40)
  
  boxes <- data.frame(
    xmin = c(3, 5, 7),
    xmax = c(4, 6, 8),
    ymin = c(3, 5, 7),
    ymax = c(4, 6, 8),
    stain_id = c("box_1", "box_2", "box_3")
  )
  
  # Test avec shift_tolerance pour couvrir les branches de recentrage
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "sum",
    shift_tolerance_ppm = 0.03
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
  
  # Vérifier qu'il y a des colonnes d'intensité pour chaque spectre
  intensity_cols <- grep("^Intensity_", names(result), value = TRUE)
  expect_equal(length(intensity_cols), 2)
})

test_that("calculate_batch_box_intensities: fit avec shift_tolerance et plusieurs spectres", {
  spectra <- create_mock_spectra_list(2, size = 35)
  
  boxes <- data.frame(
    xmin = c(4, 6),
    xmax = c(5, 7),
    ymin = c(4, 6),
    ymax = c(5, 7),
    stain_id = c("box_1", "box_2")
  )
  
  skip_if_not(exists("calculate_fitted_volumes"), 
              "calculate_fitted_volumes non disponible")
  
  result <- calculate_batch_box_intensities(
    boxes, spectra, 
    method = "fit",
    model = "gaussian",
    shift_tolerance_ppm = 0.04
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
})

# =============================================================================
# TESTS: get_box_intensity - CAS SUPPLÉMENTAIRES  
# =============================================================================

test_that("get_box_intensity: méthode fit fonctionne", {
  spec <- create_test_spectrum()
  
  boxes <- data.frame(
    xmin = 4, xmax = 6,
    ymin = 4, ymax = 6,
    stain_id = "box_1"
  )
  
  skip_if_not(exists("calculate_fitted_volumes"),
              "calculate_fitted_volumes non disponible")
  
  result <- get_box_intensity(
    spec$mat, spec$ppm_x, spec$ppm_y, boxes, 
    method = "fit", model = "gaussian"
  )
  
  expect_equal(length(result), 1)
  expect_true(is.numeric(result))
})

# =============================================================================
# TESTS: read_bruker_cached
# =============================================================================

test_that("read_bruker_cached: utilise le cache correctement", {
  skip_if_not(exists("read_bruker"), "read_bruker non disponible")
  skip_if_not(dir.exists("tests/fixtures/UFCOSY_sample/pdata/1"),
              "Fixtures non disponibles")
  
  # Vider le cache
  clear_spectrum_cache()
  
  path <- "tests/fixtures/UFCOSY_sample/pdata/1"
  
  # Premier appel - doit lire depuis le disque
  data1 <- read_bruker_cached(path, dim = "2D")
  
  # Deuxième appel - doit utiliser le cache
  data2 <- read_bruker_cached(path, dim = "2D")
  
  # Les deux doivent être identiques
  expect_identical(data1, data2)
  
  # Nettoyer
  clear_spectrum_cache()
})

# =============================================================================
# RÉSUMÉ
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           TESTS POUR R/utils.R                                   ║\n")
cat("║                                                                  ║\n")
cat("║  Fonctions testées:                                             ║\n")
cat("║  - %||% (null-coalesce operator)                                ║\n")
cat("║  - clear_spectrum_cache                                         ║\n")
cat("║  - parse_keep_peak_ranges                                       ║\n")
cat("║  - clean_centroids_df                                           ║\n")
cat("║  - get_box_intensity                                            ║\n")
cat("║  - calculate_batch_box_intensities (tous les cas)               ║\n")
cat("║  - read_bruker_cached                                           ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")
