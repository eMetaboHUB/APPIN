# CNN_detection.R - Peak Detection Functions ----
#
# Part of the CNN Peak Detection module for 2D NMR Spectra
# Contains sequential and batch detection methods
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB


# Initialize empty dataframe for detected peaks
detected_peaks <- data.frame(F2 = numeric(0), F1 = numeric(0))

#' Pad Sequence to Target Length
#'
#' Pads a numeric vector with zeros on the right to reach the target length.
#' Used to prepare variable-length spectra for fixed-size CNN input.
#'
#' @param x Numeric vector to pad
#' @param target_length Integer, desired output length
#'
#' @return Numeric vector of length target_length
#'
#' @examples
#' pad_sequence(c(1, 2, 3), 5)
#' # Returns: c(1, 2, 3, 0, 0)
pad_sequence <- function(x, target_length) {
  current_length <- length(x)
  if (current_length < target_length) {
    pad_width <- target_length - current_length
    x_padded <- c(x, rep(0, pad_width))
  } else {
    x_padded <- x
  }
  return(x_padded)
}

#' Detect Peaks with Intensity (Sequential Method)
#'
#' Performs peak detection on a 2D NMR spectrum by scanning each row and column
#' sequentially. This method is slower but may be more accurate for small spectra.
#'
#' @param rr_norm Normalized 2D spectrum matrix with ppm values as row/column names
#' @param model Compiled Keras CNN model
#' @param target_length Integer, CNN input size (default: 2048)
#'
#' @return Data frame with columns: F2 (ppm), F1 (ppm), Intensity
#'
#' @details
#' The function:
#' 1. Iterates over all rows (F2 dimension) and applies CNN prediction
#' 2. Iterates over all columns (F1 dimension) and applies CNN prediction
#' 3. Combines and deduplicates detected peaks
#'
#' Classification labels:
#' - 0: Background (no peak)
#' - 1: Peak edge
#' - 2: Peak center
#'
#' @examples
#' peaks <- get_detected_peaks_with_intensity(spectrum_matrix, model)
#'
#' @export
get_detected_peaks_with_intensity <- function(rr_norm, model, target_length = 2048) {
  n_row <- nrow(rr_norm)
  n_col <- ncol(rr_norm)
  detected_peaks <- data.frame(F2 = numeric(0), F1 = numeric(0), Intensity = numeric(0))
  
  # Safe prediction wrapper using TensorFlow directly
  # Avoids issues with Keras predict() in some environments
  safe_predict <- function(model, x) {
    x_tensor <- tensorflow::tf$convert_to_tensor(x, dtype = tensorflow::tf$float32)
    preds <- model$`__call__`(x_tensor, training = FALSE)
    list(preds[[1]]$numpy(), preds[[2]]$numpy())
  }
  
  # === Row-wise detection (F2 dimension) ===
  pb_row <- txtProgressBar(min = 0, max = n_row, style = 3)
  for (i in 1:n_row) {
    spec1D_row <- rr_norm[i, ]
    
    # Pad to target length for CNN input
    spec1D_row_padded <- pad_sequence(spec1D_row, target_length)
    
    input_tensor <- array(spec1D_row_padded, dim = c(1, target_length, 1))
    pred_row <- safe_predict(model, input_tensor)
    
    # Extract predictions only for original (non-padded) points
    class_labels <- apply(pred_row[[1]][1, 1:length(spec1D_row), ], 1, which.max) - 1
    reg_pred <- pred_row[[2]][1, 1:length(spec1D_row), ]  # [original_length, 3]
    
    # Keep points classified as peak edge (1) or peak center (2)
    idx <- which(class_labels %in% c(1, 2))
    if (length(idx) > 0) {
      peaks <- data.frame(
        F2 = as.numeric(rownames(rr_norm))[i],
        F1 = as.numeric(colnames(rr_norm))[idx],
        Intensity = reg_pred[idx, 2]  # Predicted intensity
      )
      detected_peaks <- rbind(detected_peaks, peaks)
    }
    setTxtProgressBar(pb_row, i)
  }
  close(pb_row)
  
  # === Column-wise detection (F1 dimension) ===
  pb_col <- txtProgressBar(min = 0, max = n_col, style = 3)
  for (j in 1:n_col) {
    spec1D_col <- rr_norm[, j]
    
    # Pad to target length for CNN input
    spec1D_col_padded <- pad_sequence(spec1D_col, target_length)
    
    input_tensor <- array(spec1D_col_padded, dim = c(1, target_length, 1))
    pred_col <- safe_predict(model, input_tensor)
    
    class_labels <- apply(pred_col[[1]][1, 1:length(spec1D_col), ], 1, which.max) - 1
    reg_pred <- pred_col[[2]][1, 1:length(spec1D_col), ]
    
    idx <- which(class_labels %in% c(1, 2))
    if (length(idx) > 0) {
      peaks <- data.frame(
        F2 = as.numeric(rownames(rr_norm))[idx],
        F1 = as.numeric(colnames(rr_norm))[j],
        Intensity = reg_pred[idx, 2]
      )
      detected_peaks <- rbind(detected_peaks, peaks)
    }
    setTxtProgressBar(pb_col, j)
  }
  close(pb_col)
  
  # Remove duplicate peaks (detected from both row and column scans)
  detected_peaks <- unique(detected_peaks)
  return(detected_peaks)
}


## SECTION 3: PEAK DETECTION FOR TOCSY SPECTRA (Batch Method) ----


#' Batch Peak Detection for 1D Slices
#'
#' Performs efficient batch prediction on multiple 1D slices of the spectrum.
#' Uses batch processing for significant speed improvement on large spectra.
#'
#' @param spectrum_mat 2D spectrum matrix
#' @param model Compiled Keras CNN model
#' @param axis Character, "rows" or "columns" - which dimension to scan
#' @param threshold_class Numeric, minimum probability to consider a peak (default: 0.01)
#' @param batch_size Integer, number of slices per batch (default: 64)
#' @param model_input_length Integer, CNN input size (default: 2048)
#' @param verbose Logical, print progress information
#'
#' @return Data frame with columns: F1, F2, Intensity, ppm
#'
#' @details
#' For spectra larger than model_input_length:
#' - Downsamples using evenly spaced indices
#' - Maps predictions back to original indices
#'
#' For smaller spectra:
#' - Zero-pads to model_input_length
#' - Only uses valid (non-padded) predictions
#'
#' @export
predict_peaks_1D_batch <- function(spectrum_mat, model,
                                   axis = c("rows", "columns"),
                                   threshold_class = 0.01, batch_size = 64,
                                   model_input_length = 2048, 
                                   window_overlap = 256,  # Overlap entre fenêtres
                                   verbose = TRUE) {
  axis <- match.arg(axis)
  
  
  # Safe prediction wrapper using TensorFlow directly
  safe_predict <- function(model, x) {
    x_tensor <- tensorflow::tf$convert_to_tensor(x, dtype = tensorflow::tf$float32)
    preds <- model$`__call__`(x_tensor, training = FALSE)
    list(preds[[1]]$numpy(), preds[[2]]$numpy())
  }
  
  # Transpose if scanning columns (work with rows internally)
  mat <- if (axis == "columns") t(spectrum_mat) else spectrum_mat
  n_vectors <- nrow(mat)
  n_points <- ncol(mat)
  
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SLIDING WINDOW: instead of downsampling, we cut into windows
  # ═══════════════════════════════════════════════════════════════════════════
  use_sliding_window <- n_points > model_input_length
  window_step <- model_input_length - window_overlap
  
  if (use_sliding_window) {
    n_windows <- ceiling((n_points - model_input_length) / window_step) + 1
  }
  
  detected_list <- vector("list", length = ceiling(n_vectors / batch_size))
  pb <- txtProgressBar(min = 0, max = n_vectors, style = 3)
  idx_list <- 1
  
  # Process in batches for efficiency
  for (start_idx in seq(1, n_vectors, by = batch_size)) {
    end_idx <- min(start_idx + batch_size - 1, n_vectors)
    batch_raw <- mat[start_idx:end_idx, , drop = FALSE]
    
    # ═══════════════════════════════════════════════════════════════════════════
    # OPTIMIZED SLIDING WINDOW: Batch all windows together
    # ═══════════════════════════════════════════════════════════════════════════
    
    detected_batch <- list()
    
    if (use_sliding_window) {
      # Calculer les positions des fenêtres
      window_starts <- seq(1, n_points - model_input_length + 1, by = window_step)
      if (tail(window_starts, 1) + model_input_length - 1 < n_points) {
        window_starts <- c(window_starts, n_points - model_input_length + 1)
      }
      n_windows <- length(window_starts)
      
      # Collecter TOUTES les fenêtres de ce batch de vecteurs
      all_windows <- list()
      window_info <- list()  # Pour tracker quel vecteur et quelle fenêtre
      
      for (i in seq_len(nrow(batch_raw))) {
        x <- batch_raw[i, ]
        vector_idx <- start_idx + i - 1
        
        for (w in seq_along(window_starts)) {
          win_start <- window_starts[w]
          win_end <- win_start + model_input_length - 1
          all_windows[[length(all_windows) + 1]] <- x[win_start:win_end]
          window_info[[length(window_info) + 1]] <- list(
            vector_idx = vector_idx,
            win_start = win_start
          )
        }
      }
      
      # Traiter par gros batches de fenêtres
      window_batch_size <- 128  # Nombre de fenêtres par appel CNN
      n_total_windows <- length(all_windows)
      
      for (wb_start in seq(1, n_total_windows, by = window_batch_size)) {
        wb_end <- min(wb_start + window_batch_size - 1, n_total_windows)
        
        # Construire le tensor pour ce batch de fenêtres
        batch_windows <- do.call(rbind, all_windows[wb_start:wb_end])
        input_tensor <- array(batch_windows, dim = c(nrow(batch_windows), model_input_length, 1))
        
        # Prédiction CNN en un seul appel
        preds <- safe_predict(model, input_tensor)
        prob_peak <- preds[[1]][, , 2]
        reg_intensity <- preds[[2]][, , 2]
        
        # Extraire les pics pour chaque fenêtre
        for (w_idx in seq_len(nrow(batch_windows))) {
          global_w_idx <- wb_start + w_idx - 1
          info <- window_info[[global_w_idx]]
          
          idxs_in_window <- which(prob_peak[w_idx, ] > threshold_class)
          
          if (length(idxs_in_window) > 0) {
            idxs_orig <- info$win_start + idxs_in_window - 1
            
            df <- data.frame(
              F1 = if (axis == "rows") idxs_orig else info$vector_idx,
              F2 = if (axis == "rows") info$vector_idx else idxs_orig,
              Intensity = reg_intensity[w_idx, idxs_in_window],
              ppm = NA
            )
            detected_batch[[length(detected_batch) + 1]] <- df
          }
        }
      }
      
    } else {
      # ═══════════════════════════════════════════════════════════════════
      # NO SLIDING WINDOW: normal batch processing
      # ═══════════════════════════════════════════════════════════════════
      batch_padded <- matrix(0, nrow = nrow(batch_raw), ncol = model_input_length)
      for (i in seq_len(nrow(batch_raw))) {
        x <- batch_raw[i, ]
        batch_padded[i, 1:length(x)] <- x
      }
      
      input_tensor <- array(batch_padded, dim = c(nrow(batch_padded), model_input_length, 1))
      preds <- safe_predict(model, input_tensor)
      prob_peak <- preds[[1]][, , 2]
      reg_intensity <- preds[[2]][, , 2]
      
      for (i in seq_len(nrow(batch_padded))) {
        vector_idx <- start_idx + i - 1
        idxs <- which(prob_peak[i, 1:n_points] > threshold_class)
        
        if (length(idxs) > 0) {
          df <- data.frame(
            F1 = if (axis == "rows") idxs else vector_idx,
            F2 = if (axis == "rows") vector_idx else idxs,
            Intensity = reg_intensity[i, idxs],
            ppm = NA
          )
          detected_batch[[length(detected_batch) + 1]] <- df
        }
      }
    }
    
    detected_list[[idx_list]] <- if (length(detected_batch) > 0) do.call(rbind, detected_batch) else NULL
    setTxtProgressBar(pb, end_idx)
    idx_list <- idx_list + 1
  }
  
  close(pb)
  detected <- do.call(rbind, detected_list)
  if (is.null(detected)) detected <- data.frame(F1 = numeric(), F2 = numeric(), Intensity = numeric(), ppm = numeric())
  
  # Remove duplicates
  detected <- unique(detected)
  return(detected)
}


