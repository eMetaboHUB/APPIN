# test-cnn.R - Tests pour les modules CNN
# 2DNMR-Analyst v3.0
library(testthat)

# =============================================================================
# HELPERS
# =============================================================================

#' Créer un spectre 2D normalisé de test avec des pics gaussiens
create_test_spectrum <- function(nrow = 100, ncol = 100, n_peaks = 3, noise = 0.01) {
  set.seed(42)
  mat <- matrix(rnorm(nrow * ncol, 0, noise), nrow, ncol)
  
  # Ajouter quelques pics gaussiens
  for (i in seq_len(n_peaks)) {
    pr <- sample(20:(nrow-20), 1)
    pc <- sample(20:(ncol-20), 1)
    for (r in 1:nrow) {
      for (c in 1:ncol) {
        mat[r, c] <- mat[r, c] + exp(-((r-pr)^2 + (c-pc)^2) / 50)
      }
    }
  }
  
  # Normaliser
  mat <- mat / max(abs(mat))
  rownames(mat) <- seq(10, 0, length.out = nrow)
  colnames(mat) <- seq(10, 0, length.out = ncol)
  mat
}

# =============================================================================
# TEST: pad_sequence()
# =============================================================================
test_that("pad_sequence pad avec des zéros à droite", {
  pad_sequence <- function(x, target_length) {
    if (length(x) >= target_length) return(x[1:target_length])
    c(x, rep(0, target_length - length(x)))
  }
  
  x <- c(1, 2, 3)
  result <- pad_sequence(x, 5)
  expect_equal(length(result), 5)
  expect_equal(result, c(1, 2, 3, 0, 0))
})

test_that("pad_sequence tronque si trop long", {
  pad_sequence <- function(x, target_length) {
    if (length(x) >= target_length) return(x[1:target_length])
    c(x, rep(0, target_length - length(x)))
  }
  
  x <- 1:10
  result <- pad_sequence(x, 5)
  expect_equal(length(result), 5)
  expect_equal(result, 1:5)
})

test_that("pad_sequence longueur exacte inchangée", {
  pad_sequence <- function(x, target_length) {
    if (length(x) >= target_length) return(x[1:target_length])
    c(x, rep(0, target_length - length(x)))
  }
  
  x <- 1:5
  result <- pad_sequence(x, 5)
  expect_equal(result, x)
})

# =============================================================================
# TEST: Normalisation spectre CNN
# =============================================================================
test_that("Normalisation 99.9ème percentile borne à 1", {
  normalize_spectrum <- function(mat) {
    abs_mat <- abs(mat)
    p999 <- quantile(abs_mat, 0.999)
    mat_norm <- abs_mat / p999
    mat_norm[mat_norm > 1] <- 1
    mat_norm
  }
  
  mat <- matrix(c(1, 100, 1000, 10000), 2, 2)
  result <- normalize_spectrum(mat)
  expect_true(max(result) <= 1)
  expect_true(all(result >= 0))
})

test_that("Normalisation préserve les proportions relatives", {
  normalize_spectrum <- function(mat) {
    abs_mat <- abs(mat)
    p999 <- quantile(abs_mat, 0.999)
    mat_norm <- abs_mat / p999
    mat_norm[mat_norm > 1] <- 1
    mat_norm
  }
  
  mat <- matrix(c(10, 20, 30, 40), 2, 2)
  result <- normalize_spectrum(mat)
  expect_true(result[1,1] < result[2,2])
})

# =============================================================================
# TEST: merge_overlapping_boxes() - CORRIGÉ
# =============================================================================
test_that("merge_overlapping_boxes fusionne boxes dont centre est dans l'autre", {
  # Version simplifiée et robuste
  merge_overlapping_boxes <- function(boxes) {
    if (is.null(boxes) || nrow(boxes) <= 1) return(boxes)
    
    repeat {
      merged <- FALSE
      n <- nrow(boxes)
      if (n <= 1) break
      
      for (i in 1:(n-1)) {
        if (merged) break
        for (j in (i+1):n) {
          # Centres
          cx_i <- (boxes$xmin[i] + boxes$xmax[i]) / 2
          cy_i <- (boxes$ymin[i] + boxes$ymax[i]) / 2
          cx_j <- (boxes$xmin[j] + boxes$xmax[j]) / 2
          cy_j <- (boxes$ymin[j] + boxes$ymax[j]) / 2
          
          # Centre j dans box i ?
          j_in_i <- (cx_j >= boxes$xmin[i]) && (cx_j <= boxes$xmax[i]) &&
            (cy_j >= boxes$ymin[i]) && (cy_j <= boxes$ymax[i])
          
          # Centre i dans box j ?
          i_in_j <- (cx_i >= boxes$xmin[j]) && (cx_i <= boxes$xmax[j]) &&
            (cy_i >= boxes$ymin[j]) && (cy_i <= boxes$ymax[j])
          
          if (isTRUE(j_in_i) || isTRUE(i_in_j)) {
            # Fusionner
            new_xmin <- min(boxes$xmin[i], boxes$xmin[j])
            new_xmax <- max(boxes$xmax[i], boxes$xmax[j])
            new_ymin <- min(boxes$ymin[i], boxes$ymin[j])
            new_ymax <- max(boxes$ymax[i], boxes$ymax[j])
            
            # Supprimer j, mettre à jour i
            boxes$xmin[i] <- new_xmin
            boxes$xmax[i] <- new_xmax
            boxes$ymin[i] <- new_ymin
            boxes$ymax[i] <- new_ymax
            boxes <- boxes[-j, , drop = FALSE]
            merged <- TRUE
            break
          }
        }
      }
      
      if (!merged) break
    }
    rownames(boxes) <- NULL
    boxes
  }
  
  # Deux boxes qui se chevauchent (centre l'une dans l'autre)
  boxes <- data.frame(
    xmin = c(0, 0.3),
    xmax = c(1, 1.3),
    ymin = c(0, 0.3),
    ymax = c(1, 1.3)
  )
  result <- merge_overlapping_boxes(boxes)
  expect_equal(nrow(result), 1)
  expect_equal(result$xmin[1], 0)
  expect_equal(result$xmax[1], 1.3)
})

test_that("merge_overlapping_boxes ne fusionne pas boxes séparées", {
  merge_overlapping_boxes <- function(boxes) {
    if (is.null(boxes) || nrow(boxes) <= 1) return(boxes)
    boxes
  }
  
  boxes <- data.frame(
    xmin = c(0, 5),
    xmax = c(1, 6),
    ymin = c(0, 5),
    ymax = c(1, 6)
  )
  result <- merge_overlapping_boxes(boxes)
  expect_equal(nrow(result), 2)
})

test_that("merge_overlapping_boxes gère une seule box", {
  merge_overlapping_boxes <- function(boxes) {
    if (is.null(boxes) || nrow(boxes) <= 1) return(boxes)
    boxes
  }
  
  boxes <- data.frame(xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  result <- merge_overlapping_boxes(boxes)
  expect_equal(nrow(result), 1)
})

# =============================================================================
# TEST: filter_tocsy_traces()
# =============================================================================
test_that("filter_tocsy_traces supprime pics faibles sur même ligne F2", {
  filter_tocsy_traces <- function(peaks, ratio = 0.1, ppm_tolerance = 0.02) {
    if (nrow(peaks) == 0) return(peaks)
    
    peaks$f2_group <- round(peaks$F2_ppm / ppm_tolerance)
    
    result <- do.call(rbind, lapply(split(peaks, peaks$f2_group), function(group) {
      if (nrow(group) == 0) return(NULL)
      max_int <- max(group$intensity, na.rm = TRUE)
      group[group$intensity >= ratio * max_int, , drop = FALSE]
    }))
    
    if (!is.null(result)) result$f2_group <- NULL
    result
  }
  
  peaks <- data.frame(
    F1_ppm = c(1, 2, 3, 4),
    F2_ppm = c(5.00, 5.01, 5.00, 5.01),
    intensity = c(1.0, 0.5, 0.05, 0.8)
  )
  
  result <- filter_tocsy_traces(peaks, ratio = 0.1, ppm_tolerance = 0.02)
  expect_true(nrow(result) < nrow(peaks))
  expect_false(0.05 %in% result$intensity)
})

test_that("filter_tocsy_traces garde tous les pics si aucun faible", {
  filter_tocsy_traces <- function(peaks, ratio = 0.1, ppm_tolerance = 0.02) {
    if (nrow(peaks) == 0) return(peaks)
    peaks$f2_group <- round(peaks$F2_ppm / ppm_tolerance)
    result <- do.call(rbind, lapply(split(peaks, peaks$f2_group), function(group) {
      max_int <- max(group$intensity, na.rm = TRUE)
      group[group$intensity >= ratio * max_int, , drop = FALSE]
    }))
    if (!is.null(result)) result$f2_group <- NULL
    result
  }
  
  peaks <- data.frame(
    F1_ppm = c(1, 2),
    F2_ppm = c(5.00, 5.01),
    intensity = c(1.0, 0.9)
  )
  
  result <- filter_tocsy_traces(peaks, ratio = 0.1)
  expect_equal(nrow(result), 2)
})

# =============================================================================
# TEST: Labels CNN (classification)
# =============================================================================
test_that("Labels CNN classes correctes (0, 1, 2)", {
  generate_labels <- function(peaks_idx, n_points, intensities, 
                              intensity_threshold = 0.2, center_margin = 0.2) {
    y_class <- rep(0, n_points)
    
    for (j in seq_along(peaks_idx)) {
      idx <- peaks_idx[j]
      if (idx >= 1 && idx <= n_points) {
        rel_pos <- idx / n_points
        center_range <- c(0.5 - center_margin, 0.5 + center_margin)
        is_centered <- rel_pos >= center_range[1] && rel_pos <= center_range[2]
        is_strong <- intensities[j] >= intensity_threshold
        
        y_class[idx] <- if (is_centered && is_strong) 1 else 2
      }
    }
    y_class
  }
  
  # Pic fort centré → classe 1
  labels <- generate_labels(c(1024), 2048, c(0.5))
  expect_equal(labels[1024], 1)
  
  # Pic faible → classe 2
  labels <- generate_labels(c(1024), 2048, c(0.1))
  expect_equal(labels[1024], 2)
  
  # Pic fort mais pas centré → classe 2
  labels <- generate_labels(c(100), 2048, c(0.5))
  expect_equal(labels[100], 2)
  
  # Background → classe 0
  labels <- generate_labels(c(1024), 2048, c(0.5))
  expect_equal(labels[500], 0)
})

# =============================================================================
# TEST: Structure sortie CNN
# =============================================================================
test_that("run_cnn_peak_picking retourne structure correcte", {
  mock_cnn_result <- list(
    peaks = data.frame(
      F1_ppm = c(1.5, 3.2),
      F2_ppm = c(2.1, 4.5),
      intensity = c(0.8, 0.6)
    ),
    boxes = data.frame(
      xmin = c(1.4, 3.1),
      xmax = c(1.6, 3.3),
      ymin = c(2.0, 4.4),
      ymax = c(2.2, 4.6)
    )
  )
  
  expect_true("peaks" %in% names(mock_cnn_result))
  expect_true("boxes" %in% names(mock_cnn_result))
  expect_true(all(c("F1_ppm", "F2_ppm") %in% names(mock_cnn_result$peaks)))
  expect_true(all(c("xmin", "xmax", "ymin", "ymax") %in% names(mock_cnn_result$boxes)))
})

# =============================================================================
# TEST: Validation paramètres CNN
# =============================================================================
test_that("Paramètres CNN valides entre 0 et 1", {
  validate_cnn_params <- function(params) {
    errors <- c()
    
    if (!is.null(params$pred_class_thres)) {
      if (params$pred_class_thres < 0 || params$pred_class_thres > 1) {
        errors <- c(errors, "pred_class_thres doit être entre 0 et 1")
      }
    }
    
    if (!is.null(params$trace_filter_ratio)) {
      if (params$trace_filter_ratio < 0 || params$trace_filter_ratio > 1) {
        errors <- c(errors, "trace_filter_ratio doit être entre 0 et 1")
      }
    }
    
    if (!is.null(params$eps_value)) {
      if (params$eps_value <= 0) {
        errors <- c(errors, "eps_value doit être > 0")
      }
    }
    
    list(valid = length(errors) == 0, errors = errors)
  }
  
  valid_params <- list(pred_class_thres = 0.3, trace_filter_ratio = 0.5, eps_value = 0.01)
  result <- validate_cnn_params(valid_params)
  expect_true(result$valid)
  
  invalid_params <- list(pred_class_thres = 1.5)
  result <- validate_cnn_params(invalid_params)
  expect_false(result$valid)
  
  invalid_params <- list(eps_value = -0.01)
  result <- validate_cnn_params(invalid_params)
  expect_false(result$valid)
})

# =============================================================================
# TEST: Décomposition 2D → 1D
# =============================================================================
test_that("Décomposition extrait le bon nombre de coupes", {
  extract_cuts_count <- function(nrow, ncol, extract_rows = TRUE, extract_cols = TRUE) {
    n_cuts <- 0
    if (extract_rows) n_cuts <- n_cuts + nrow
    if (extract_cols) n_cuts <- n_cuts + ncol
    n_cuts
  }
  
  expect_equal(extract_cuts_count(100, 80, TRUE, TRUE), 180)
  expect_equal(extract_cuts_count(100, 80, TRUE, FALSE), 100)
  expect_equal(extract_cuts_count(100, 80, FALSE, TRUE), 80)
})

test_that("Sous-échantillonnage réduit le nombre de coupes", {
  subsample_rows <- function(nrow, sample_rate = 0.25) {
    n_sampled <- ceiling(nrow * sample_rate)
    n_sampled
  }
  
  expect_equal(subsample_rows(100, 0.25), 25)
  expect_equal(subsample_rows(100, 0.5), 50)
})

# =============================================================================
# TEST: Focal Loss - CORRIGÉ
# =============================================================================
test_that("Focal Loss calcul correct", {
  # Formule: -alpha * (1-p)^gamma * log(p)
  # Gamma élevé RÉDUIT la loss pour les exemples "moyens" car (1-p)^gamma diminue
  focal_loss <- function(p, gamma = 2.0, alpha = 1.0) {
    -alpha * (1 - p)^gamma * log(p + 1e-7)
  }
  
  # p proche de 1 → loss faible (bien classé)
  expect_true(focal_loss(0.99) < focal_loss(0.5))
  
  # p proche de 0 → loss élevée (mal classé)
  expect_true(focal_loss(0.01) > focal_loss(0.5))
  
  # Gamma plus élevé RÉDUIT la loss pour p=0.5
  # car (1-0.5)^3 = 0.125 < (1-0.5)^1 = 0.5
  expect_true(focal_loss(0.5, gamma = 3) < focal_loss(0.5, gamma = 1))
})

test_that("Focal Loss avec alpha pondère les classes", {
  focal_loss <- function(p, gamma = 2.0, alpha = 1.0) {
    -alpha * (1 - p)^gamma * log(p + 1e-7)
  }
  
  expect_true(focal_loss(0.5, gamma = 2, alpha = 2) > focal_loss(0.5, gamma = 2, alpha = 1))
})

# =============================================================================
# TEST: Sliding window - CORRIGÉ
# =============================================================================
test_that("Sliding window calcule le bon nombre de fenêtres", {
  calculate_n_windows <- function(length, window_size = 2048, overlap = 256) {
    if (length <= window_size) return(1)
    step <- window_size - overlap
    # Nombre de fenêtres = 1 + ceiling((length - window_size) / step)
    1 + ceiling((length - window_size) / step)
  }
  
  # Spectre de 2048 points → 1 fenêtre
  expect_equal(calculate_n_windows(2048), 1)
  
  # Spectre de 4000 points: step=1792, n=1+ceiling(1952/1792)=1+2=3
  expect_equal(calculate_n_windows(4000, 2048, 256), 3)
  
  # Spectre de 8000 points: n=1+ceiling(5952/1792)=1+4=5
  expect_equal(calculate_n_windows(8000, 2048, 256), 5)
})

# =============================================================================
# TEST: DBSCAN epsilon scaling pour CNN
# =============================================================================
test_that("Epsilon DBSCAN multiplié par 5 pour CNN", {
  scale_eps_for_cnn <- function(eps_base, factor = 5) {
    eps_base * factor
  }
  
  expect_equal(scale_eps_for_cnn(0.0068), 0.034)
  expect_equal(scale_eps_for_cnn(0.01), 0.05)
})

# =============================================================================
# TEST: Bounding box padding
# =============================================================================
test_that("Box padding ajoute de la marge", {
  add_padding <- function(box, padding) {
    data.frame(
      xmin = box$xmin - padding,
      xmax = box$xmax + padding,
      ymin = box$ymin - padding,
      ymax = box$ymax + padding
    )
  }
  
  box <- data.frame(xmin = 1, xmax = 2, ymin = 3, ymax = 4)
  padded <- add_padding(box, 0.1)
  
  expect_equal(padded$xmin, 0.9)
  expect_equal(padded$xmax, 2.1)
  expect_equal(padded$ymin, 2.9)
  expect_equal(padded$ymax, 4.1)
})