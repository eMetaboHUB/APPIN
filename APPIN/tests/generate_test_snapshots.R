# generate_test_snapshots.R - Génère les fichiers de référence pour les tests de régression
# 2DNMR-Analyst v3.0
#
# USAGE:
#   source("tests/generate_test_snapshots.R")
#   generate_all_snapshots()
#
# Ce script génère des fichiers de référence ("snapshots") à partir du comportement
# actuel des fonctions. Ces snapshots servent ensuite à détecter les régressions.
#
# ⚠️  IMPORTANT: N'exécuter ce script que quand le comportement actuel est CORRECT !
#     Les snapshots générés deviennent la "vérité" pour les tests futurs.
#
# =============================================================================

library(digest)  # Pour les checksums

# =============================================================================
# CONFIGURATION
# =============================================================================

SNAPSHOTS_DIR <- "tests/snapshots"
FIXTURES_DIR <- "tests/fixtures"

# =============================================================================
# FONCTIONS UTILITAIRES
# =============================================================================

#' Initialise le dossier des snapshots
init_snapshots_dir <- function() {
  if (!dir.exists(SNAPSHOTS_DIR)) {
    dir.create(SNAPSHOTS_DIR, recursive = TRUE)
    cat("✅ Dossier créé:", SNAPSHOTS_DIR, "\n")
  }
}

#' Sauvegarde un snapshot avec métadonnées
save_snapshot <- function(data, name, description = "") {
  init_snapshots_dir()
  
  snapshot <- list(
    name = name,
    description = description,
    created_at = Sys.time(),
    r_version = R.version.string,
    data = data,
    checksum = digest::digest(data, algo = "md5")
  )
  
  filepath <- file.path(SNAPSHOTS_DIR, paste0(name, ".rds"))
  saveRDS(snapshot, filepath)
  
  cat(sprintf("✅ Snapshot sauvegardé: %s\n", filepath))
  cat(sprintf("   Checksum: %s\n", snapshot$checksum))
  
  invisible(snapshot)
}

#' Charge un snapshot existant
load_snapshot <- function(name) {
  filepath <- file.path(SNAPSHOTS_DIR, paste0(name, ".rds"))
  if (!file.exists(filepath)) {
    stop("Snapshot non trouvé: ", filepath)
  }
  readRDS(filepath)
}

#' Liste les snapshots disponibles
list_snapshots <- function() {
  files <- list.files(SNAPSHOTS_DIR, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(files) == 0) {
    cat("Aucun snapshot trouvé.\n")
    return(invisible(NULL))
  }
  
  cat("\n📸 Snapshots disponibles:\n\n")
  
  for (f in files) {
    snap <- readRDS(f)
    cat(sprintf("  • %s\n", snap$name))
    cat(sprintf("    Créé: %s\n", snap$created_at))
    cat(sprintf("    Checksum: %s\n", snap$checksum))
    if (snap$description != "") {
      cat(sprintf("    Description: %s\n", snap$description))
    }
    cat("\n")
  }
}

# =============================================================================
# GÉNÉRATION DES SNAPSHOTS
# =============================================================================

#' Génère le snapshot pour read_bruker sur UFCOSY
generate_snapshot_read_bruker <- function() {
  cat("\n📸 Génération snapshot: read_bruker_ufcosy\n")
  
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  
  if (!dir.exists(ufcosy_path)) {
    cat("⚠️  Fixture UFCOSY non trouvée, skip.\n")
    return(invisible(NULL))
  }
  
  # Charger les sources si nécessaire
  if (!exists("read_bruker")) {
    source("Function/Read_2DNMR_spectrum.R")
  }
  
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  
  # Extraire les infos essentielles (pas la matrice complète, trop volumineuse)
  snapshot_data <- list(
    dimensions = dim(result$spectrumData),
    ppm_f1_range = range(as.numeric(rownames(result$spectrumData))),
    ppm_f2_range = range(as.numeric(colnames(result$spectrumData))),
    intensity_range = range(result$spectrumData),
    intensity_mean = mean(result$spectrumData),
    intensity_sd = sd(as.vector(result$spectrumData)),
    n_positive = sum(result$spectrumData > 0),
    n_negative = sum(result$spectrumData < 0),
    # Checksum de la matrice complète
    matrix_checksum = digest::digest(result$spectrumData, algo = "md5")
  )
  
  save_snapshot(
    snapshot_data, 
    "read_bruker_ufcosy",
    "Lecture du spectre UFCOSY de référence"
  )
}

#' Génère le snapshot pour peak_pick_2d_nt2 sur UFCOSY
generate_snapshot_peak_picking <- function() {
  cat("\n📸 Génération snapshot: peak_picking_ufcosy\n")
  
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  
  if (!dir.exists(ufcosy_path)) {
    cat("⚠️  Fixture UFCOSY non trouvée, skip.\n")
    return(invisible(NULL))
  }
  
  # Charger les sources
  if (!exists("read_bruker")) source("Function/Read_2DNMR_spectrum.R")
  if (!exists("peak_pick_2d_nt2")) source("Function/Peak_picking.R")
  
  # Charger le spectre
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  # Paramètres fixes pour reproductibilité
  threshold <- quantile(abs(mat), 0.95)
  
  # Peak picking
  peaks_result <- peak_pick_2d_nt2(
    mat,
    threshold_value = threshold,
    spectrum_type = "UFCOSY",
    verbose = FALSE
  )
  
  # Extraire les données essentielles
  snapshot_data <- list(
    threshold_used = threshold,
    n_peaks = if (!is.null(peaks_result$centroids)) nrow(peaks_result$centroids) else 0,
    n_boxes = if (!is.null(peaks_result$bounding_boxes)) nrow(peaks_result$bounding_boxes) else 0,
    # Sauvegarder les centroïdes (positions des pics)
    centroids = peaks_result$centroids,
    # Sauvegarder les bounding boxes
    bounding_boxes = peaks_result$bounding_boxes,
    # Statistiques de cluster si disponibles
    cluster_stats = peaks_result$cluster_stats
  )
  
  save_snapshot(
    snapshot_data,
    "peak_picking_ufcosy",
    "Peak picking sur UFCOSY avec seuil au 95ème percentile"
  )
}

#' Génère le snapshot pour les seuils automatiques
generate_snapshot_thresholds <- function() {
  cat("\n📸 Génération snapshot: thresholds_ufcosy\n")
  
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  
  if (!dir.exists(ufcosy_path)) {
    cat("⚠️  Fixture UFCOSY non trouvée, skip.\n")
    return(invisible(NULL))
  }
  
  # Charger les sources
  if (!exists("read_bruker")) source("Function/Read_2DNMR_spectrum.R")
  if (!exists("seuil_bruit_multiplicatif")) source("Function/Vizualisation.R")
  
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  snapshot_data <- list(
    seuil_bruit_3sigma = seuil_bruit_multiplicatif(mat, 3),
    seuil_bruit_5sigma = seuil_bruit_multiplicatif(mat, 5),
    seuil_max_5pct = seuil_max_pourcentage(mat, 0.05),
    seuil_max_10pct = seuil_max_pourcentage(mat, 0.10),
    max_intensity = max(abs(mat)),
    sd_intensity = sd(as.vector(mat))
  )
  
  save_snapshot(
    snapshot_data,
    "thresholds_ufcosy",
    "Calculs de seuils automatiques sur UFCOSY"
  )
}

#' Génère le snapshot pour le fitting
generate_snapshot_fitting <- function() {
  cat("\n📸 Génération snapshot: fitting_ufcosy\n")
  
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  
  if (!dir.exists(ufcosy_path)) {
    cat("⚠️  Fixture UFCOSY non trouvée, skip.\n")
    return(invisible(NULL))
  }
  
  # Charger les sources
  if (!exists("read_bruker")) source("Function/Read_2DNMR_spectrum.R")
  if (!exists("fit_2d_peak")) source("Function/Peak_fitting.R")
  
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  ppm_x <- as.numeric(colnames(mat))
  ppm_y <- as.numeric(rownames(mat))
  
  # Trouver le max et créer une box autour
  max_idx <- which(mat == max(mat), arr.ind = TRUE)[1, ]
  center_y <- ppm_y[max_idx[1]]
  center_x <- ppm_x[max_idx[2]]
  
  box <- data.frame(
    xmin = center_x - 0.1,
    xmax = center_x + 0.1,
    ymin = center_y - 0.1,
    ymax = center_y + 0.1
  )
  
  # Fitting
  fit_gaussian <- fit_2d_peak(mat, ppm_x, ppm_y, box, "gaussian")
  fit_voigt <- fit_2d_peak(mat, ppm_x, ppm_y, box, "voigt")
  
  snapshot_data <- list(
    box_used = box,
    max_position = list(x = center_x, y = center_y),
    fit_gaussian = fit_gaussian,
    fit_voigt = fit_voigt
  )
  
  save_snapshot(
    snapshot_data,
    "fitting_ufcosy",
    "Fitting Gaussien et Voigt sur le pic max du UFCOSY"
  )
}

#' Génère le snapshot pour le CNN peak picking sur UFCOSY
#'
#' Capture l'output complet de run_cnn_peak_picking pour détecter les régressions:
#' - Nombre de peaks et de boxes
#' - Distribution statistique des positions et intensités
#' - Checksums MD5 des tables complètes
#' - Un échantillon des premières lignes pour inspection manuelle
generate_snapshot_cnn <- function() {
  cat("\n📸 Génération snapshot: cnn_peak_picking_ufcosy\n")
  
  ufcosy_path <- file.path(FIXTURES_DIR, "UFCOSY_sample", "pdata", "1")
  
  if (!dir.exists(ufcosy_path)) {
    cat("⚠️  Fixture UFCOSY non trouvée, skip.\n")
    return(invisible(NULL))
  }
  
  # Charger les sources
  if (!exists("read_bruker")) source("Function/Read_2DNMR_spectrum.R")
  if (!exists("run_cnn_peak_picking")) source("Function/CNN_shiny.R")
  
  # Vérifier que le modèle est disponible
  model <- tryCatch(
    get_cnn_model("UFCOSY"),
    error = function(e) {
      cat("⚠️  Modèle CNN UFCOSY non disponible:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(model)) {
    cat("⚠️  Skip: CNN model not loadable.\n")
    return(invisible(NULL))
  }
  
  # Charger le spectre
  result <- read_bruker(dir = ufcosy_path, dim = "2D")
  mat <- result$spectrumData
  
  # Normaliser le spectre (CNN attend un spectre normalisé)
  rr_norm <- mat / max(abs(mat))
  
  # Paramètres fixes et reproductibles
  params <- list(
    eps_value = 0.01,
    pred_class_thres = 0.01,
    int_thres = 0.001,
    trace_filter_ratio = 0.1,
    use_filters = FALSE,
    disable_clustering = FALSE,
    box_padding = NULL
  )
  
  # CNN peak picking
  cnn_result <- run_cnn_peak_picking(
    rr_norm,
    model = model,
    params = params,
    spectrum_type = "UFCOSY",
    method = "batch",
    verbose = FALSE
  )
  
  # Préparer des stats robustes pour la régression
  peaks_df <- cnn_result$peaks
  boxes_df <- cnn_result$boxes
  
  peak_stats <- if (!is.null(peaks_df) && nrow(peaks_df) > 0) {
    list(
      n = nrow(peaks_df),
      f2_ppm_range = range(peaks_df$F2_ppm, na.rm = TRUE),
      f1_ppm_range = range(peaks_df$F1_ppm, na.rm = TRUE),
      intensity_range = range(peaks_df$stain_intensity, na.rm = TRUE),
      intensity_mean = mean(peaks_df$stain_intensity, na.rm = TRUE),
      intensity_sum = sum(peaks_df$stain_intensity, na.rm = TRUE),
      n_clusters = length(unique(peaks_df$cluster_db))
    )
  } else {
    list(n = 0)
  }
  
  box_stats <- if (!is.null(boxes_df) && nrow(boxes_df) > 0) {
    list(
      n = nrow(boxes_df),
      mean_width_f2 = mean(boxes_df$xmax - boxes_df$xmin, na.rm = TRUE),
      mean_width_f1 = mean(boxes_df$ymax - boxes_df$ymin, na.rm = TRUE),
      total_area = sum((boxes_df$xmax - boxes_df$xmin) *
                         (boxes_df$ymax - boxes_df$ymin), na.rm = TRUE)
    )
  } else {
    list(n = 0)
  }
  
  snapshot_data <- list(
    spectrum_type = "UFCOSY",
    params_used = params,
    method = "batch",
    peak_stats = peak_stats,
    box_stats = box_stats,
    # Tables complètes pour checksum
    peaks = peaks_df,
    boxes = boxes_df,
    # Checksums pour détection rapide de régression
    peaks_checksum = if (!is.null(peaks_df)) digest::digest(peaks_df, algo = "md5") else NA,
    boxes_checksum = if (!is.null(boxes_df)) digest::digest(boxes_df, algo = "md5") else NA,
    # Échantillon des 5 premiers peaks triés par intensité (pour inspection)
    top_peaks_sample = if (!is.null(peaks_df) && nrow(peaks_df) > 0) {
      peaks_df[order(-peaks_df$stain_intensity)[1:min(5, nrow(peaks_df))], ]
    } else NULL
  )
  
  save_snapshot(
    snapshot_data,
    "cnn_peak_picking_ufcosy",
    "CNN peak picking sur UFCOSY (méthode batch, params par défaut)"
  )
}

# =============================================================================
# FONCTION PRINCIPALE
# =============================================================================

#' Génère tous les snapshots
#' @param force Régénérer même si les snapshots existent déjà
generate_all_snapshots <- function(force = FALSE) {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
  cat("║           📸 GÉNÉRATION DES SNAPSHOTS DE RÉFÉRENCE 📸                 ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  if (!force && length(list.files(SNAPSHOTS_DIR, pattern = "\\.rds$")) > 0) {
    cat("⚠️  Des snapshots existent déjà.\n")
    cat("   Utilisez generate_all_snapshots(force = TRUE) pour les régénérer.\n")
    cat("   Ou supprimez manuellement les fichiers dans", SNAPSHOTS_DIR, "\n\n")
    
    response <- readline("Continuer quand même ? (o/n) : ")
    if (tolower(response) != "o") {
      cat("Annulé.\n")
      return(invisible(NULL))
    }
  }
  
  cat("⚠️  ATTENTION: Ce script va créer des fichiers de référence basés sur\n")
  cat("   le comportement ACTUEL des fonctions. Assurez-vous que tout fonctionne\n")
  cat("   correctement avant de continuer !\n\n")
  
  # Charger toutes les sources
  cat("📂 Chargement des sources...\n")
  source("Function/Read_2DNMR_spectrum.R")
  source("Function/Vizualisation.R")
  source("Function/Peak_picking.R")
  source("Function/Peak_fitting.R")
  source("Function/CNN_shiny.R")
  cat("✅ Sources chargées\n")
  
  # Générer les snapshots
  generate_snapshot_read_bruker()
  generate_snapshot_thresholds()
  generate_snapshot_peak_picking()
  generate_snapshot_fitting()
  generate_snapshot_cnn()
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("✅ Génération terminée !\n")
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("\nUtilisez list_snapshots() pour voir les snapshots créés.\n")
  cat("Lancez run_test('regression') pour exécuter les tests de régression.\n\n")
}

# =============================================================================
# AFFICHAGE AU CHARGEMENT
# =============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║           📸 GÉNÉRATEUR DE SNAPSHOTS DE RÉFÉRENCE 📸                  ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                       ║\n")
cat("║   generate_all_snapshots()  → Créer tous les snapshots                ║\n")
cat("║   list_snapshots()          → Voir les snapshots existants            ║\n")
cat("║                                                                       ║\n")
cat("║   ⚠️  Exécuter UNIQUEMENT si le comportement actuel est correct !     ║\n")
cat("║                                                                       ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")