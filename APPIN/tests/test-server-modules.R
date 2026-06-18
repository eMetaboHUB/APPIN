# =============================================================================
# Tests Server-Side pour APPIN
# Utilise shiny::testServer() - pas besoin de navigateur !
# =============================================================================
#
# USAGE:
#   setwd("C:/Users/juguibert/Documents/APPIN_1.0.0")
#   testthat::test_file("tests/test-server-modules.R")
#
# =============================================================================

library(testthat)
library(shiny)

# Déterminer le chemin racine d'APPIN
# (fonctionne que le test soit lancé depuis tests/ ou depuis la racine)
if (file.exists("../R/mod_load_data.R")) {
  APPIN_ROOT <- normalizePath("..")
} else if (file.exists("R/mod_load_data.R")) {
  APPIN_ROOT <- normalizePath(".")
} else {
  
  APPIN_ROOT <- "C:/Users/juguibert/Documents/APPIN_1.0.0"
}

cat("APPIN_ROOT:", APPIN_ROOT, "\n")

# Charger les dépendances APPIN
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(dbscan)
})

# Sourcer les modules
source(file.path(APPIN_ROOT, "R/mod_load_data.R"))
source(file.path(APPIN_ROOT, "R/mod_peak_picking.R"))
source(file.path(APPIN_ROOT, "R/mod_integration.R"))

# Sourcer les fonctions utilitaires
source(file.path(APPIN_ROOT, "Function/Read_2DNMR_spectrum.R"))
source(file.path(APPIN_ROOT, "Function/Peak_picking.R"))
source(file.path(APPIN_ROOT, "Function/Peak_fitting.R"))

# =============================================================================
# HELPERS POUR LES TESTS
# =============================================================================

#' Créer un spectre 2D synthétique pour les tests
#' @param size Taille de la matrice (size x size)
#' @param n_peaks Nombre de pics à générer
#' @return Liste mimant la structure de bruker_data
create_mock_spectrum <- function(size = 128, n_peaks = 5) {
  set.seed(42)
  
  # Créer une matrice de bruit
  mat <- matrix(rnorm(size * size, mean = 0, sd = 0.1), nrow = size, ncol = size)
  
  # Ajouter des pics gaussiens
  for (i in seq_len(n_peaks)) {
    cx <- sample(20:(size-20), 1)
    cy <- sample(20:(size-20), 1)
    sigma <- runif(1, 3, 8)
    amplitude <- runif(1, 5, 20)
    
    for (row in max(1, cy-15):min(size, cy+15)) {
      for (col in max(1, cx-15):min(size, cx+15)) {
        dist <- sqrt((row - cy)^2 + (col - cx)^2)
        mat[row, col] <- mat[row, col] + amplitude * exp(-dist^2 / (2 * sigma^2))
      }
    }
  }
  
  # Créer les axes ppm (typique pour 1H-1H COSY)
  ppm_f2 <- seq(10, 0, length.out = size)  # 1H direct
  
  ppm_f1 <- seq(10, 0, length.out = size)  # 1H indirect
  
  colnames(mat) <- as.character(round(ppm_f2, 4))
  rownames(mat) <- as.character(round(ppm_f1, 4))
  
  list(
    spectrumData = mat,
    ppm_f2 = ppm_f2,
    ppm_f1 = ppm_f1,
    metadata = list(
      nucleus_f2 = "1H",
      nucleus_f1 = "1H",
      experiment = "COSY"
    )
  )
}

#' Créer des boxes de test
#' @param n Nombre de boxes
#' @return data.frame avec xmin, xmax, ymin, ymax, stain_id
create_mock_boxes <- function(n = 3) {
  set.seed(123)
  data.frame(
    xmin = runif(n, 1, 4),
    xmax = runif(n, 5, 8),
    ymin = runif(n, 1, 4),
    ymax = runif(n, 5, 8),
    stain_id = paste0("box_", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# TESTS: mod_load_data
# =============================================================================

test_that("mod_load_data: initialisation correcte", {
  skip_if_not_installed("shinyFiles")
  library(shinyFiles)
  
  # Créer un status_msg mock
  status_msg <- reactiveVal("Ready")
  
  testServer(mod_load_data_server, args = list(status_msg = status_msg), {
    
    # Vérifier que les valeurs initiales sont vides
    expect_equal(spectra_list(), list())  # Liste vide, pas NULL
    expect_equal(length(spectra_list()), 0)
    expect_null(bruker_data())
    
  })
})

test_that("mod_load_data: set_bruker_data fonctionne", {
  skip_if_not_installed("shinyFiles")
  library(shinyFiles)
  
  status_msg <- reactiveVal("Ready")
  
  testServer(mod_load_data_server, args = list(status_msg = status_msg), {
    
    # Créer un spectre mock
    mock_spectrum <- create_mock_spectrum(size = 64, n_peaks = 3)
    
    # Utiliser le setter
    session$returned$set_bruker_data(mock_spectrum)
    
    # Vérifier que bruker_data est mis à jour
    expect_false(is.null(bruker_data()))
    expect_true("spectrumData" %in% names(bruker_data()))
    expect_equal(dim(bruker_data()$spectrumData), c(64, 64))
    
  })
})

test_that("mod_load_data: set_spectra_list fonctionne", {
  skip_if_not_installed("shinyFiles")
  library(shinyFiles)
  
  status_msg <- reactiveVal("Ready")
  
  testServer(mod_load_data_server, args = list(status_msg = status_msg), {
    
    # Créer une liste de spectres mock
    mock_spectra <- list(
      "spectrum1" = create_mock_spectrum(size = 64),
      "spectrum2" = create_mock_spectrum(size = 64)
    )
    
    # Utiliser le setter
    session$returned$set_spectra_list(mock_spectra)
    
    # Vérifier
    expect_equal(length(spectra_list()), 2)
    expect_true("spectrum1" %in% names(spectra_list()))
    expect_true("spectrum2" %in% names(spectra_list()))
    
  })
})

# =============================================================================
# TESTS: mod_integration
# =============================================================================

test_that("mod_integration: initialisation correcte", {
  
  status_msg <- reactiveVal("Ready")
  
  # Mock de load_data
  load_data <- list(
    bruker_data = reactiveVal(NULL)
  )
  
  # Mock de rv (reactive values partagées)
  rv <- list(
    modifiable_boxes = reactiveVal(NULL),
    fit_results_data = reactiveVal(NULL),
    last_fit_method = reactiveVal(NULL)
  )
  
  testServer(mod_integration_server, 
             args = list(status_msg = status_msg, load_data = load_data, rv = rv), {
               
               # Vérifier les valeurs initiales
               expect_null(integration_results())
               expect_false(integration_done())
               
               # Vérifier la méthode par défaut
               session$setInputs(integration_method = "sum")
               expect_equal(effective_integration_method(), "sum")
               
             })
})

test_that("mod_integration: sélection de méthode fonctionne", {
  
  status_msg <- reactiveVal("Ready")
  load_data <- list(bruker_data = reactiveVal(NULL))
  rv <- list(
    modifiable_boxes = reactiveVal(NULL),
    fit_results_data = reactiveVal(NULL),
    last_fit_method = reactiveVal(NULL)
  )
  
  testServer(mod_integration_server, 
             args = list(status_msg = status_msg, load_data = load_data, rv = rv), {
               
               # Tester méthode sum
               session$setInputs(integration_method = "sum", integration_method_fit = "")
               expect_equal(effective_integration_method(), "sum")
               
               # Tester méthode gaussian
               session$setInputs(integration_method = "", integration_method_fit = "gaussian")
               expect_equal(effective_integration_method(), "gaussian")
               
               # Tester méthode voigt
               session$setInputs(integration_method_fit = "voigt")
               expect_equal(effective_integration_method(), "voigt")
               
             })
})

test_that("mod_integration: intégration sum fonctionne", {
  
  status_msg <- reactiveVal("Ready")
  
  # Créer des données mock
  mock_spectrum <- create_mock_spectrum(size = 128, n_peaks = 5)
  mock_boxes <- data.frame(
    xmin = c(2, 4, 6),
    xmax = c(3, 5, 7),
    ymin = c(2, 4, 6),
    ymax = c(3, 5, 7),
    stain_id = c("box_1", "box_2", "box_3"),
    stringsAsFactors = FALSE
  )
  
  load_data <- list(
    bruker_data = reactiveVal(mock_spectrum)
  )
  
  rv <- list(
    modifiable_boxes = reactiveVal(mock_boxes),
    fit_results_data = reactiveVal(NULL),
    last_fit_method = reactiveVal(NULL)
  )
  
  testServer(mod_integration_server, 
             args = list(status_msg = status_msg, load_data = load_data, rv = rv), {
               
               # Configurer la méthode
               session$setInputs(integration_method = "sum", integration_method_fit = "")
               
               # Lancer l'intégration
               session$setInputs(run_integration = 1)
               
               # Vérifier les résultats
               expect_true(integration_done())
               expect_false(is.null(integration_results()))
               
               results <- integration_results()
               expect_equal(nrow(results), 3)  # 3 boxes
               expect_true("stain_id" %in% names(results))
               expect_true("intensity" %in% names(results))
               expect_true(all(results$method == "sum"))
               
             })
})

# =============================================================================
# TESTS: mod_peak_picking (helper functions)
# =============================================================================

test_that("parse_keep_peak_ranges: parsing correct", {
  
  # Test avec input valide
  result <- parse_keep_peak_ranges("0.5,-0.5; 1,0.8; 1.55,1.45;")
  expect_equal(length(result), 3)
  expect_equal(result[[1]], c(0.5, -0.5))
  expect_equal(result[[2]], c(1, 0.8))
  expect_equal(result[[3]], c(1.55, 1.45))
  
  # Test avec input vide
  expect_null(parse_keep_peak_ranges(""))
  expect_null(parse_keep_peak_ranges(NULL))
  
  # Test avec input invalide
  expect_null(parse_keep_peak_ranges("invalid"))
  
  # Test avec espaces
  result <- parse_keep_peak_ranges("  0.5 , -0.5 ;  1 , 0.8  ")
  expect_equal(length(result), 2)
  
})

test_that("mod_peak_picking: valeurs retournées correctes", {
  
  status_msg <- reactiveVal("Ready")
  
  # Mocks minimaux
  load_data <- list(
    bruker_data = reactiveVal(NULL)
  )
  
  data_reactives <- list(
    result_data_list = reactiveVal(list()),
    spectrum_params = reactiveVal(list(neighborhood_size = 5))
  )
  
  rv <- list(
    centroids_data = reactiveVal(NULL),
    fixed_boxes = reactiveVal(NULL),
    modifiable_boxes = reactiveVal(NULL),
    reference_boxes = reactiveVal(NULL),
    contour_plot_base = reactiveVal(NULL)
  )
  
  refresh_nmr_plot <- function() NULL
  parent_input <- reactiveValues(spectrum_type = "TOCSY", selected_subfolder = NULL, contour_start = 0.1)
  
  testServer(mod_peak_picking_server, 
             args = list(
               status_msg = status_msg,
               load_data = load_data,
               data_reactives = data_reactives,
               rv = rv,
               refresh_nmr_plot = refresh_nmr_plot,
               parent_input = parent_input
             ), {
               
               # Configurer les inputs
               session$setInputs(
                 eps_value = 0.01,
                 disable_clustering = FALSE,
                 keep_peak_ranges_text = "0.5,-0.5;"
               )
               
               # Vérifier les valeurs retournées
               expect_equal(session$returned$eps_value(), 0.01)
               expect_false(session$returned$disable_clustering())
               expect_equal(session$returned$keep_peak_ranges_text(), "0.5,-0.5;")
               
             })
})

# =============================================================================
# TESTS: Intégration entre modules
# =============================================================================

test_that("Workflow: load_data -> integration", {
  
  status_msg <- reactiveVal("Ready")
  
  # 1. Créer load_data avec un spectre
  mock_spectrum <- create_mock_spectrum(size = 128, n_peaks = 5)
  
  load_data <- list(
    bruker_data = reactiveVal(mock_spectrum),
    spectra_list = reactiveVal(list("test_spectrum" = mock_spectrum)),
    set_bruker_data = function(data) { load_data$bruker_data(data) }
  )
  
  # 2. Créer des boxes (simulant le résultat de peak_picking)
  mock_boxes <- data.frame(
    xmin = c(1.5, 3.5, 5.5),
    xmax = c(2.5, 4.5, 6.5),
    ymin = c(1.5, 3.5, 5.5),
    ymax = c(2.5, 4.5, 6.5),
    stain_id = c("peak_1", "peak_2", "peak_3"),
    stringsAsFactors = FALSE
  )
  
  rv <- list(
    modifiable_boxes = reactiveVal(mock_boxes),
    fit_results_data = reactiveVal(NULL),
    last_fit_method = reactiveVal(NULL)
  )
  
  # 3. Tester l'intégration
  testServer(mod_integration_server, 
             args = list(status_msg = status_msg, load_data = load_data, rv = rv), {
               
               # Configurer
               session$setInputs(integration_method = "sum", integration_method_fit = "")
               
               # Lancer
               session$setInputs(run_integration = 1)
               
               # Vérifier le workflow complet
               expect_true(integration_done())
               
               results <- integration_results()
               expect_equal(nrow(results), 3)
               expect_true(all(c("stain_id", "F2_ppm", "F1_ppm", "intensity") %in% names(results)))
               
               # Les intensités devraient être des nombres positifs (on a des pics)
               expect_true(all(is.numeric(results$intensity)))
               
             })
})

# =============================================================================
# RAPPORT DE TESTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           TESTS SERVER-SIDE POUR APPIN                           ║\n")
cat("║                                                                  ║\n")
cat("║  Ces tests utilisent shiny::testServer() et ne nécessitent      ║\n")
cat("║  PAS de navigateur. Ils testent la logique des modules.         ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")