# test-regression.R - Tests de régression basés sur les snapshots
# 2DNMR-Analyst v3.0
#
# Ces tests comparent les résultats actuels avec des snapshots de référence
# générés quand le comportement était correct.
#
# Si un test échoue, cela signifie qu'une modification a changé le comportement
# d'une fonction — c'est peut-être une régression (bug) ou une amélioration
# intentionnelle (dans ce cas, régénérer les snapshots).
#
library(testthat)
library(digest)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Chemin vers les snapshots et fixtures
get_snapshots_path <- function() {
  paths <- c(
    file.path("..", "snapshots"),
    file.path("tests", "snapshots"),
    "snapshots"
  )
  for (p in paths) {
    if (dir.exists(p)) return(p)
  }
  NULL
}

get_fixtures_path <- function() {
  paths <- c(
    file.path("..", "fixtures", "UFCOSY_sample", "pdata", "1"),
    file.path("tests", "fixtures", "UFCOSY_sample", "pdata", "1"),
    file.path("fixtures", "UFCOSY_sample", "pdata", "1")
  )
  for (p in paths) {
    if (dir.exists(p) && file.exists(file.path(p, "2rr"))) return(p)
  }
  NULL
}

# Helpers
snapshots_available <- function() {
  path <- get_snapshots_path()
  !is.null(path) && length(list.files(path, pattern = "\\.rds$")) > 0
}

fixtures_available <- function() {
  !is.null(get_fixtures_path())
}

load_snapshot <- function(name) {
  path <- get_snapshots_path()
  filepath <- file.path(path, paste0(name, ".rds"))
  if (!file.exists(filepath)) return(NULL)
  readRDS(filepath)
}

# =============================================================================
# TESTS DE RÉGRESSION: read_bruker
# =============================================================================
test_that("RÉGRESSION: read_bruker produit les mêmes dimensions", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  
  snapshot <- load_snapshot("read_bruker_ufcosy")
  skip_if(is.null(snapshot), "Snapshot read_bruker_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  
  # Comparer les dimensions
  expect_equal(
    dim(result$spectrumData),
    snapshot$data$dimensions,
    info = "Les dimensions du spectre ont changé !"
  )
})

test_that("RÉGRESSION: read_bruker produit les mêmes plages ppm", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  
  snapshot <- load_snapshot("read_bruker_ufcosy")
  skip_if(is.null(snapshot), "Snapshot read_bruker_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  
  ppm_f1_range <- range(as.numeric(rownames(result$spectrumData)))
  ppm_f2_range <- range(as.numeric(colnames(result$spectrumData)))
  
  expect_equal(ppm_f1_range, snapshot$data$ppm_f1_range, tolerance = 1e-6,
               info = "La plage ppm F1 a changé !")
  expect_equal(ppm_f2_range, snapshot$data$ppm_f2_range, tolerance = 1e-6,
               info = "La plage ppm F2 a changé !")
})

test_that("RÉGRESSION: read_bruker produit les mêmes données (checksum)", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  
  snapshot <- load_snapshot("read_bruker_ufcosy")
  skip_if(is.null(snapshot), "Snapshot read_bruker_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  
  current_checksum <- digest::digest(result$spectrumData, algo = "md5")
  
  expect_equal(
    current_checksum,
    snapshot$data$matrix_checksum,
    info = "Les données spectrales ont changé ! Checksum différent."
  )
})

# =============================================================================
# TESTS DE RÉGRESSION: Seuils automatiques
# =============================================================================
test_that("RÉGRESSION: seuil_bruit_multiplicatif produit les mêmes valeurs", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("seuil_bruit_multiplicatif"), "Fonction non chargée")
  
  snapshot <- load_snapshot("thresholds_ufcosy")
  skip_if(is.null(snapshot), "Snapshot thresholds_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  seuil_3sigma <- seuil_bruit_multiplicatif(mat, 3)
  seuil_5sigma <- seuil_bruit_multiplicatif(mat, 5)
  
  expect_equal(seuil_3sigma, snapshot$data$seuil_bruit_3sigma, tolerance = 1e-6,
               info = "seuil_bruit_multiplicatif (3σ) a changé !")
  expect_equal(seuil_5sigma, snapshot$data$seuil_bruit_5sigma, tolerance = 1e-6,
               info = "seuil_bruit_multiplicatif (5σ) a changé !")
})

test_that("RÉGRESSION: seuil_max_pourcentage produit les mêmes valeurs", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("seuil_max_pourcentage"), "Fonction non chargée")
  
  snapshot <- load_snapshot("thresholds_ufcosy")
  skip_if(is.null(snapshot), "Snapshot thresholds_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  seuil_5pct <- seuil_max_pourcentage(mat, 0.05)
  seuil_10pct <- seuil_max_pourcentage(mat, 0.10)
  
  expect_equal(seuil_5pct, snapshot$data$seuil_max_5pct, tolerance = 1e-6,
               info = "seuil_max_pourcentage (5%) a changé !")
  expect_equal(seuil_10pct, snapshot$data$seuil_max_10pct, tolerance = 1e-6,
               info = "seuil_max_pourcentage (10%) a changé !")
})

# =============================================================================
# TESTS DE RÉGRESSION: Peak picking
# =============================================================================
test_that("RÉGRESSION: peak_pick_2d_nt2 détecte le même nombre de pics", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("peak_pick_2d_nt2"), "Fonction non chargée")
  
  snapshot <- load_snapshot("peak_picking_ufcosy")
  skip_if(is.null(snapshot), "Snapshot peak_picking_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  # Utiliser le même seuil que lors de la génération du snapshot
  threshold <- snapshot$data$threshold_used
  
  peaks_result <- peak_pick_2d_nt2(
    mat,
    threshold_value = threshold,
    spectrum_type = "UFCOSY",
    verbose = FALSE
  )
  
  n_peaks <- if (!is.null(peaks_result$centroids)) nrow(peaks_result$centroids) else 0
  n_boxes <- if (!is.null(peaks_result$bounding_boxes)) nrow(peaks_result$bounding_boxes) else 0
  
  expect_equal(n_peaks, snapshot$data$n_peaks,
               info = sprintf("Nombre de pics différent: %d vs %d (référence)", 
                              n_peaks, snapshot$data$n_peaks))
  expect_equal(n_boxes, snapshot$data$n_boxes,
               info = sprintf("Nombre de boxes différent: %d vs %d (référence)",
                              n_boxes, snapshot$data$n_boxes))
})

test_that("RÉGRESSION: peak_pick_2d_nt2 détecte les mêmes positions", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("peak_pick_2d_nt2"), "Fonction non chargée")
  
  snapshot <- load_snapshot("peak_picking_ufcosy")
  skip_if(is.null(snapshot), "Snapshot peak_picking_ufcosy non trouvé")
  skip_if(is.null(snapshot$data$centroids), "Pas de centroïdes dans le snapshot")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  threshold <- snapshot$data$threshold_used
  
  peaks_result <- peak_pick_2d_nt2(
    mat,
    threshold_value = threshold,
    spectrum_type = "UFCOSY",
    verbose = FALSE
  )
  
  # Comparer les centroïdes
  if (!is.null(peaks_result$centroids) && nrow(peaks_result$centroids) > 0) {
    # Trier par position pour comparaison
    current <- peaks_result$centroids[order(peaks_result$centroids$F1_ppm, 
                                            peaks_result$centroids$F2_ppm), ]
    reference <- snapshot$data$centroids[order(snapshot$data$centroids$F1_ppm,
                                               snapshot$data$centroids$F2_ppm), ]
    
    # Comparer les positions (avec tolérance)
    if (nrow(current) == nrow(reference)) {
      expect_equal(current$F1_ppm, reference$F1_ppm, tolerance = 0.001,
                   info = "Positions F1 des centroïdes différentes !")
      expect_equal(current$F2_ppm, reference$F2_ppm, tolerance = 0.001,
                   info = "Positions F2 des centroïdes différentes !")
    }
  }
})

# =============================================================================
# TESTS DE RÉGRESSION: Fitting
# =============================================================================
test_that("RÉGRESSION: fit_2d_peak (gaussian) produit le même volume", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("fit_2d_peak"), "Fonction non chargée")
  
  snapshot <- load_snapshot("fitting_ufcosy")
  skip_if(is.null(snapshot), "Snapshot fitting_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  
  # Utiliser la même box que dans le snapshot
  box <- snapshot$data$box_used
  
  fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, "gaussian")
  
  # Comparer les volumes (avec tolérance de 1%)
  expect_equal(
    fit_result$volume, 
    snapshot$data$fit_gaussian$volume,
    tolerance = 0.01 * abs(snapshot$data$fit_gaussian$volume),
    info = "Volume du fitting gaussien différent !"
  )
})

test_that("RÉGRESSION: fit_2d_peak (voigt) produit le même volume", {
  skip_if_not(snapshots_available(), "Snapshots non disponibles")
  skip_if_not(fixtures_available(), "Fixtures non disponibles")
  skip_if_not(exists("fit_2d_peak"), "Fonction non chargée")
  
  snapshot <- load_snapshot("fitting_ufcosy")
  skip_if(is.null(snapshot), "Snapshot fitting_ufcosy non trouvé")
  
  ufcosy_path <- get_fixtures_path()
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  
  box <- snapshot$data$box_used
  
  fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, "voigt")
  
  expect_equal(
    fit_result$volume,
    snapshot$data$fit_voigt$volume,
    tolerance = 0.01 * abs(snapshot$data$fit_voigt$volume),
    info = "Volume du fitting voigt différent !"
  )
})

# =============================================================================
# INFO
# =============================================================================
#
# Si un test échoue :
#
# 1. RÉGRESSION INVOLONTAIRE (bug) :
#    → Corriger le code et relancer les tests
#
# 2. AMÉLIORATION INTENTIONNELLE :
#    → Vérifier que le nouveau comportement est correct
#    → Régénérer les snapshots : source("tests/generate_test_snapshots.R")
#    → generate_all_snapshots(force = TRUE)
#
# =============================================================================