# run_tests.R - 2DNMR-Analyst Tests Unitaires v3.0
# ═══════════════════════════════════════════════════════════════════════════════
# 
# Ce fichier permet d'exécuter tous les tests unitaires de l'application.
# 
# USAGE:
#   source("tests/run_tests.R")
#   run_all_tests()        # Tous les tests
#   run_test("cnn")        # Un fichier spécifique
#   list_tests()           # Voir les tests disponibles
#   check_fixtures()       # Vérifier les fixtures disponibles
#
# ═══════════════════════════════════════════════════════════════════════════════

# Installation automatique de testthat si nécessaire
if (!require("testthat", quietly = TRUE)) {
  
  install.packages("testthat")
}
library(testthat)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Fichiers sources à charger
SOURCE_FILES <- c(
  
  # Fonctions métier principales
  "Function/Read_2DNMR_spectrum.R",
  "Function/Vizualisation.R",
  "Function/Peak_picking.R",
  "Function/Peak_fitting.R",
  
  
  # Modules CNN (optionnels - ne génèrent pas d'erreur si absents)
  "Function/CNN_shiny.R",
  
  # Utilitaires Shiny
  "R/utils.R"
)

# Fichiers de test attendus
TEST_FILES <- c(
  "test-read_bruker.R",
  "test-read_bruker_real.R",
  "test-peak_picking.R",
  "test-peak_fitting.R",
  "test-vizualisation.R",
  "test-utils.R",
  "test-data_structures.R",
  "test-cnn.R",
  "test-regression.R"
)

# Chemins
FIXTURES_DIR <- "tests/fixtures"
SNAPSHOTS_DIR <- "tests/snapshots"

# =============================================================================
# FONCTIONS UTILITAIRES
# =============================================================================

#' Charge les fichiers sources nécessaires aux tests
#' @param verbose Afficher les messages de chargement
#' @param include_cnn Inclure les modules CNN
load_sources <- function(verbose = TRUE, include_cnn = TRUE) {
  if (verbose) {
    cat("\n📂 Chargement des sources...\n")
  }
  
  files_to_load <- SOURCE_FILES
  
  # Filtrer les fichiers CNN si non souhaités
  
  
  if (!include_cnn) {
    files_to_load <- files_to_load[!grepl("CNN", files_to_load)]
  }
  
  loaded <- 0
  failed <- 0
  skipped <- 0
  
  for (f in files_to_load) {
    if (file.exists(f)) {
      result <- tryCatch({
        source(f, local = FALSE)
        if (verbose) cat("  ✅", f, "\n")
        "success"
      }, error = function(e) {
        if (verbose) cat("  ❌", f, "-", conditionMessage(e), "\n")
        "error"
      })
      
      if (result == "success") loaded <- loaded + 1
      else failed <- failed + 1
      
    } else {
      if (verbose) cat("  ⚠️ ", f, "(non trouvé)\n")
      skipped <- skipped + 1
    }
  }
  
  if (verbose) {
    cat(sprintf("\n   Chargés: %d | Échoués: %d | Ignorés: %d\n", 
                loaded, failed, skipped))
  }
  
  invisible(list(loaded = loaded, failed = failed, skipped = skipped))
}

#' Vérifie les fixtures disponibles
check_fixtures <- function() {
  cat("\n📁 Vérification des fixtures et snapshots...\n\n")
  
  # === FIXTURES ===
  cat("  FIXTURES:\n")
  
  if (!dir.exists(FIXTURES_DIR)) {
    cat("  ⚠️  Dossier fixtures non trouvé:", FIXTURES_DIR, "\n")
    cat("     Créez-le avec: dir.create('tests/fixtures', recursive = TRUE)\n\n")
    ufcosy_ok <- FALSE
  } else {
    # Vérifier UFCOSY
    ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
    ufcosy_files <- c("2rr", "procs", "proc2s")
    
    cat("  UFCOSY_sample:\n")
    if (dir.exists(ufcosy_path)) {
      for (f in ufcosy_files) {
        fpath <- file.path(ufcosy_path, f)
        if (file.exists(fpath)) {
          size <- file.size(fpath)
          cat(sprintf("    ✅ %s (%.1f KB)\n", f, size / 1024))
        } else {
          cat(sprintf("    ❌ %s (manquant)\n", f))
        }
      }
    } else {
      cat("    ❌ Dossier non trouvé\n")
      cat("       Chemin attendu:", ufcosy_path, "\n")
    }
    
    ufcosy_ok <- dir.exists(ufcosy_path) && 
      all(file.exists(file.path(ufcosy_path, ufcosy_files)))
  }
  
  cat("\n")
  
  # === SNAPSHOTS ===
  cat("  SNAPSHOTS (tests de régression):\n")
  
  if (!dir.exists(SNAPSHOTS_DIR)) {
    cat("    ⚠️  Aucun snapshot trouvé\n")
    cat("       Pour générer: source('tests/generate_test_snapshots.R')\n")
    cat("                     generate_all_snapshots()\n")
    snapshots_ok <- FALSE
  } else {
    snapshot_files <- list.files(SNAPSHOTS_DIR, pattern = "\\.rds$")
    if (length(snapshot_files) == 0) {
      cat("    ⚠️  Dossier vide\n")
      snapshots_ok <- FALSE
    } else {
      for (f in snapshot_files) {
        fpath <- file.path(SNAPSHOTS_DIR, f)
        snap <- tryCatch(readRDS(fpath), error = function(e) NULL)
        if (!is.null(snap)) {
          cat(sprintf("    ✅ %s (créé: %s)\n", 
                      gsub("\\.rds$", "", f),
                      format(snap$created_at, "%Y-%m-%d")))
        } else {
          cat(sprintf("    ❌ %s (corrompu)\n", f))
        }
      }
      snapshots_ok <- TRUE
    }
  }
  
  cat("\n")
  
  # Résumé
  if (ufcosy_ok && snapshots_ok) {
    cat("  ✅ Tout est prêt - tous les tests seront exécutés\n\n")
  } else if (ufcosy_ok) {
    cat("  ⚠️  Fixtures OK, mais pas de snapshots\n")
    cat("     Les tests de régression seront ignorés.\n")
    cat("     Générez les snapshots avec generate_all_snapshots()\n\n")
  } else if (snapshots_ok) {
    cat("  ⚠️  Snapshots OK, mais fixtures manquantes\n")
    cat("     Certains tests seront ignorés.\n\n")
  } else {
    cat("  ⚠️  Fixtures et snapshots manquants\n")
    cat("     Les tests avec fichiers réels seront ignorés.\n\n")
  }
  
  invisible(list(fixtures = ufcosy_ok, snapshots = snapshots_ok))
}

#' Liste les tests disponibles avec nombre de tests
list_tests <- function() {
  test_dir <- "tests/testthat"
  
  if (!dir.exists(test_dir)) {
    cat("\n❌ Dossier tests/testthat non trouvé!\n\n")
    return(invisible(NULL))
  }
  
  files <- list.files(test_dir, pattern = "^test-.*\\.R$")
  
  if (length(files) == 0) {
    cat("\n⚠️  Aucun fichier de test trouvé dans", test_dir, "\n\n")
    return(invisible(NULL))
  }
  
  cat("\n╔═══════════════════════════════════════════════════════╗\n")
  cat("║              📋 TESTS DISPONIBLES                     ║\n")
  cat("╠═══════════════════════════════════════════════════════╣\n")
  
  total_tests <- 0
  
  for (f in sort(files)) {
    fpath <- file.path(test_dir, f)
    lines <- readLines(fpath, warn = FALSE)
    n_tests <- sum(grepl("^test_that\\(", lines))
    total_tests <- total_tests + n_tests
    
    # Nom sans préfixe/suffixe
    name <- gsub("^test-|\\.R$", "", f)
    
    # Indicateur si fixtures requises
    needs_fixtures <- any(grepl("skip_if_not.*fixtures", lines))
    fixture_icon <- if (needs_fixtures) " 📁" else ""
    
    cat(sprintf("║  • %-22s %3d tests%s\n", name, n_tests, fixture_icon))
  }
  
  cat("╠═══════════════════════════════════════════════════════╣\n")
  cat(sprintf("║  TOTAL: %d tests dans %d fichiers                    ║\n", 
              total_tests, length(files)))
  cat("╚═══════════════════════════════════════════════════════╝\n")
  cat("\n  📁 = nécessite fixtures (fichiers Bruker réels)\n")
  cat("\n  Usage: run_test('cnn') ou run_test('peak_picking')\n\n")
  
  invisible(files)
}

# =============================================================================
# FONCTIONS D'EXÉCUTION DES TESTS
# =============================================================================

#' Exécute tous les tests
#' @param load Charger les sources avant les tests
#' @param include_cnn Inclure les tests CNN
#' @param reporter Type de reporter testthat ("summary", "minimal", "progress")
run_all_tests <- function(load = TRUE, include_cnn = TRUE, reporter = "summary") {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
  cat("║                                                                       ║\n")
  cat("║           🧪  2DNMR-Analyst v3.0 - TESTS UNITAIRES  🧪                ║\n")
  cat("║                                                                       ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
  
  # Vérifier le dossier de tests
  if (!dir.exists("tests/testthat")) {
    stop("❌ Dossier tests/testthat non trouvé! ",
         "Exécutez ce script depuis la racine du projet.")
  }
  
  # Charger les sources
  if (load) {
    load_sources(verbose = TRUE, include_cnn = include_cnn)
  }
  
  # Vérifier les fixtures (informatif seulement)
  cat("\n")
  check_fixtures()
  
  # Compter les tests
  test_files <- list.files("tests/testthat", pattern = "^test-.*\\.R$")
  n_files <- length(test_files)
  
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("                    🚀 EXÉCUTION DE %d FICHIERS DE TESTS\n", n_files))
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # Mesurer le temps
  start_time <- Sys.time()
  
  # Exécuter les tests
  results <- test_dir("tests/testthat", reporter = reporter)
  
  # Temps écoulé
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("                    ⏱️  Terminé en %.2f secondes\n", elapsed))
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  invisible(results)
}

#' Exécute un fichier de test spécifique
#' @param name Nom du test (ex: "cnn", "peak_picking", "read_bruker_real")
#' @param load Charger les sources avant
run_test <- function(name, load = TRUE) {
  # Normaliser le nom
  if (!grepl("^test-", name)) name <- paste0("test-", name)
  if (!grepl("\\.R$", name)) name <- paste0(name, ".R")
  
  path <- file.path("tests", "testthat", name)
  
  if (!file.exists(path)) {
    cat("\n❌ Fichier non trouvé:", path, "\n")
    cat("   Utilisez list_tests() pour voir les tests disponibles.\n\n")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  🧪 Exécution: %-38s ║\n", name))
  cat("╚═══════════════════════════════════════════════════════╝\n")
  
  # Charger les sources
  if (load) {
    # Déterminer si CNN nécessaire
    include_cnn <- grepl("cnn", name, ignore.case = TRUE)
    load_sources(verbose = TRUE, include_cnn = include_cnn)
  }
  
  cat("\n")
  
  # Exécuter
  start_time <- Sys.time()
  results <- test_file(path, reporter = "summary")
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  cat(sprintf("\n⏱️  Terminé en %.2f secondes\n\n", elapsed))
  
  invisible(results)
}

#' Exécute uniquement les tests rapides (sans fixtures)
run_fast_tests <- function() {
  cat("\n⚡ Exécution des tests rapides (sans fixtures)...\n")
  
  load_sources(verbose = TRUE)
  
  # Filtrer les fichiers qui nécessitent des fixtures
  test_files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", 
                           full.names = TRUE)
  
  fast_files <- sapply(test_files, function(f) {
    lines <- readLines(f, warn = FALSE)
    !any(grepl("skip_if_not.*fixtures", lines))
  })
  
  fast_files <- test_files[fast_files]
  
  cat(sprintf("\n🧪 %d fichiers de tests rapides\n\n", length(fast_files)))
  
  for (f in fast_files) {
    test_file(f, reporter = "minimal")
  }
}

# =============================================================================
# AFFICHAGE AU CHARGEMENT
# =============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                       ║\n")
cat("║           🧪  APPIN v1.0 - FRAMEWORK DE TESTS  🧪             ║\n")
cat("║                                                                       ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                       ║\n")
cat("║   COMMANDES DISPONIBLES:                                              ║\n")
cat("║                                                                       ║\n")
cat("║   run_all_tests()     → Exécuter tous les tests                       ║\n")
cat("║   run_fast_tests()    → Tests rapides (sans fichiers Bruker)          ║\n")
cat("║   run_test('nom')     → Un fichier spécifique                         ║\n")
cat("║   list_tests()        → Voir les tests disponibles                    ║\n")
cat("║   check_fixtures()    → Vérifier les fixtures Bruker                  ║\n")
cat("║                                                                       ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                       ║\n")
cat("║   EXEMPLES:                                                           ║\n")
cat("║                                                                       ║\n")
cat("║   run_test('cnn')           → Tests du module CNN                     ║\n")
cat("║   run_test('peak_picking')  → Tests détection de pics                 ║\n")
cat("║   run_test('read_bruker_real') → Tests avec fichier UFCOSY réel       ║\n")
cat("║                                                                       ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")