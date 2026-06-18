# CNN_main.R - Main Entry Point Function ----
#
# Part of the CNN Peak Detection module for 2D NMR Spectra
# Contains the main pipeline function run_cnn_peak_picking()
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB


#' Run CNN Peak Picking Pipeline
#'
#' @param rr_norm Normalized 2D spectrum matrix
#' @param model Compiled Keras CNN model (if NULL, auto-selects based on spectrum_type)
#' @param params List of parameters
#' @param spectrum_type Character: "TOCSY", "COSY", "HSQC", "UFCOSY"
#' @param method Character, "batch" or "classique"
#' @param target_length Integer, CNN input size (default: 2048)
#' @param threshold_class Numeric, probability threshold
#' @param batch_size Integer, batch size (default: 64)
#' @param step Integer, downsampling factor (default: 4)
#' @param keep_peak_ranges List of numeric vectors. Specific F2 ranges where only top peaks are kept.
#' @param verbose Logical, print progress
#' @return List containing peaks, boxes, shapes
#' @export
run_cnn_peak_picking <- function(rr_norm, model = NULL, params,
                                 spectrum_type = "TOCSY",
                                 method = c("batch", "classique"),
                                 target_length = 2048,
                                 threshold_class = params$pred_class_thres,
                                 batch_size = 64, step = 4,
                                 keep_peak_ranges = NULL,
                                 verbose = TRUE) {
  
  if (is.null(rr_norm) || !is.matrix(rr_norm)) {
    stop("rr_norm must be a normalized matrix")
  }
  
  spectrum_type <- toupper(spectrum_type)
  method <- match.arg(method)
  
  # === Auto-select model based on spectrum_type ===
  if (is.null(model)) {
    if (verbose) cat(sprintf("\n🧠 Auto-selecting CNN model for: %s\n", spectrum_type))
    model <- get_cnn_model(spectrum_type)
  }
  
  if (is.null(model)) {
    stop("CNN model is not loaded")
  }
  
  n_row <- nrow(rr_norm)
  n_col <- ncol(rr_norm)
  
  if (spectrum_type == "HSQC" && verbose) {
    F2_range <- range(as.numeric(rownames(rr_norm)))
    F1_range <- range(as.numeric(colnames(rr_norm)))
    cat(sprintf("HSQC: F2(13C)=%.1f-%.1f, F1(1H)=%.1f-%.1f ppm\n",
                F2_range[1], F2_range[2], F1_range[1], F1_range[2]))
  }
  
  pad_sequence <- function(x, target_length) {
    current_length <- length(x)
    if (current_length < target_length) c(x, rep(0, target_length - current_length)) else x
  }
  
  safe_predict <- function(model, x) {
    x_tensor <- tensorflow::tf$convert_to_tensor(x, dtype = tensorflow::tf$float32)
    preds <- model$`__call__`(x_tensor, training = FALSE)
    list(preds[[1]]$numpy(), preds[[2]]$numpy())
  }
  
  detected_peaks <- data.frame(F2 = numeric(0), F1 = numeric(0), Intensity = numeric(0))
  
  # === Peak Detection ===
  if (method == "classique") {
    if (verbose) cat("Using sequential method\n")
    
    for (i in seq_len(n_row)) {
      spec1D_row <- rr_norm[i, ]
      spec1D_row_padded <- pad_sequence(spec1D_row, target_length)
      input_tensor <- array(spec1D_row_padded, dim = c(1, target_length, 1))
      pred_row <- safe_predict(model, input_tensor)
      class_labels <- apply(pred_row[[1]][1, 1:length(spec1D_row), ], 1, which.max) - 1
      reg_pred <- pred_row[[2]][1, 1:length(spec1D_row), ]
      idx <- which(class_labels %in% c(1, 2))
      if (length(idx) > 0) {
        peaks <- data.frame(F2 = as.numeric(rownames(rr_norm))[i],
                            F1 = as.numeric(colnames(rr_norm))[idx],
                            Intensity = reg_pred[idx, 2])
        detected_peaks <- rbind(detected_peaks, peaks)
      }
    }
    
    for (j in seq_len(n_col)) {
      spec1D_col <- rr_norm[, j]
      spec1D_col_padded <- pad_sequence(spec1D_col, target_length)
      input_tensor <- array(spec1D_col_padded, dim = c(1, target_length, 1))
      pred_col <- safe_predict(model, input_tensor)
      class_labels <- apply(pred_col[[1]][1, 1:length(spec1D_col), ], 1, which.max) - 1
      reg_pred <- pred_col[[2]][1, 1:length(spec1D_col), ]
      idx <- which(class_labels %in% c(1, 2))
      if (length(idx) > 0) {
        peaks <- data.frame(F2 = as.numeric(rownames(rr_norm))[idx],
                            F1 = as.numeric(colnames(rr_norm))[j],
                            Intensity = reg_pred[idx, 2])
        detected_peaks <- rbind(detected_peaks, peaks)
      }
    }
    
  } else if (method == "batch") {
    if (verbose) cat("Using batch method\n")
    peaks_rows <- predict_peaks_1D_batch(rr_norm, model, axis = "rows",
                                         threshold_class = threshold_class,
                                         batch_size = batch_size, verbose = verbose)
    peaks_cols <- predict_peaks_1D_batch(rr_norm, model, axis = "columns",
                                         threshold_class = threshold_class,
                                         batch_size = batch_size, verbose = verbose)
    detected_peaks <- rbind(peaks_rows, peaks_cols) %>% unique()
  }
  
  if (nrow(detected_peaks) == 0) {
    warning("No peaks detected")
    return(list(peaks = NULL, boxes = NULL, shapes = NULL))
  }
  
  # === Trace Filter ===
  trace_filter_ratio <- if (!is.null(params$trace_filter_ratio)) params$trace_filter_ratio else 0.1
  f2_tolerance_ppm <- 0.02
  
  if (trace_filter_ratio > 0) {
    y_vals <- as.numeric(rownames(rr_norm))
    detected_peaks$real_intensity <- sapply(seq_len(nrow(detected_peaks)), function(i) {
      f2_idx <- max(1, min(round(detected_peaks$F2[i]), nrow(rr_norm)))
      f1_idx <- max(1, min(round(detected_peaks$F1[i]), ncol(rr_norm)))
      rr_norm[f2_idx, f1_idx]
    })
    detected_peaks$F2_ppm_temp <- sapply(detected_peaks$F2, function(f2_idx) {
      y_vals[max(1, min(round(f2_idx), length(y_vals)))]
    })
    detected_peaks$f2_group <- round(detected_peaks$F2_ppm_temp / f2_tolerance_ppm) * f2_tolerance_ppm
    group_max <- tapply(detected_peaks$real_intensity, detected_peaks$f2_group, max)
    detected_peaks$group_max_intensity <- group_max[as.character(detected_peaks$f2_group)]
    detected_peaks$intensity_ratio <- detected_peaks$real_intensity / detected_peaks$group_max_intensity
    detected_peaks <- detected_peaks[detected_peaks$intensity_ratio >= trace_filter_ratio, ]
    detected_peaks$f2_group <- NULL
    detected_peaks$group_max_intensity <- NULL
    detected_peaks$intensity_ratio <- NULL
    detected_peaks$F2_ppm_temp <- NULL
    
    if (nrow(detected_peaks) == 0) {
      warning("No peaks after trace filter")
      return(list(peaks = NULL, boxes = NULL, shapes = NULL))
    }
  }
  
  # === Post-filters ===
  use_filters <- if (!is.null(params$use_filters)) params$use_filters else FALSE
  if (verbose) cat("Post-filters:", if(use_filters) "ENABLED" else "DISABLED", "\n")
  
  if (use_filters) {
    result <- filter_peaks_by_proportion(detected_peaks, threshold = 0.5,
                                         intensity_threshold = params$int_thres)
    if (is.null(result$filtered_peaks) || nrow(result$filtered_peaks) == 0) {
      peaks_clean_filtered <- detected_peaks
    } else {
      peaks_clean_filtered <- filter_noisy_columns(result$filtered_peaks, threshold_ratio = 0.3)
      if (is.null(peaks_clean_filtered) || nrow(peaks_clean_filtered) == 0) {
        peaks_clean_filtered <- result$filtered_peaks
      }
    }
  } else {
    peaks_clean_filtered <- detected_peaks
  }
  
  if (verbose) cat("Peaks before DBSCAN:", nrow(peaks_clean_filtered), "\n")
  
  # === No clustering mode ===
  disable_clustering <- if (!is.null(params$disable_clustering)) params$disable_clustering else FALSE
  
  if (disable_clustering) {
    x_vals <- as.numeric(colnames(rr_norm))
    y_vals <- as.numeric(rownames(rr_norm))
    
    # SWAP F1/F2 to match Peak_picking.R convention (same as CNN_clustering.R)
    # CNN F1 -> our F2_ppm, CNN F2 -> our F1_ppm
    peaks_ppm <- peaks_clean_filtered %>%
      mutate(
        F2_idx = pmin(pmax(round(F1), 1), length(x_vals)),  # CNN F1 -> our F2
        F1_idx = pmin(pmax(round(F2), 1), length(y_vals)),  # CNN F2 -> our F1
        F2_ppm = x_vals[F2_idx],
        F1_ppm = y_vals[F1_idx],
        stain_id = paste0("cnn_", row_number()),
        stain_intensity = Intensity,
        cluster_db = row_number()
      ) %>%
      select(F2_ppm, F1_ppm, stain_intensity, cluster_db, stain_id)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FILTER VERTICAL TRACES (bleeding artifacts) - Same logic as Peak_picking.R
    # Groups peaks by F2 column and keeps only top N per column
    # ═══════════════════════════════════════════════════════════════════════════
    if (nrow(peaks_ppm) > 1 && spectrum_type %in% c("TOCSY", "COSY", "UFCOSY")) {
      
      f2_tolerance <- switch(spectrum_type,
                             "TOCSY" = 0.02,
                             "COSY" = 0.02,
                             "UFCOSY" = 0.025,
                             0.02
      )
      
      keep_top_n <- switch(spectrum_type,
                           "TOCSY" = 6,
                           "COSY" = 4,
                           "UFCOSY" = 5,
                           5
      )
      
      n_before <- nrow(peaks_ppm)
      
      # Group by F2 column and rank by intensity
      peaks_ppm <- peaks_ppm %>%
        dplyr::mutate(
          f2_group = round(F2_ppm / f2_tolerance) * f2_tolerance
        ) %>%
        dplyr::group_by(f2_group) %>%
        dplyr::mutate(
          intensity_rank = dplyr::row_number(dplyr::desc(stain_intensity))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(intensity_rank <= keep_top_n) %>%
        dplyr::select(-f2_group, -intensity_rank)
      
      # Reassign stain_id after filtering
      peaks_ppm <- peaks_ppm %>%
        dplyr::mutate(
          stain_id = paste0("cnn_", dplyr::row_number()),
          cluster_db = dplyr::row_number()
        )
      
      n_removed <- n_before - nrow(peaks_ppm)
      if (verbose && n_removed > 0) {
        cat(sprintf("🧹 [%s] Filtered %d vertical trace peaks (kept top %d per F2 column)\n",
                    spectrum_type, n_removed, keep_top_n))
      }
    }
    
    box_padding <- params$eps_value * 2
    boxes <- peaks_ppm %>%
      mutate(xmin = F2_ppm - box_padding, xmax = F2_ppm + box_padding,
             ymin = F1_ppm - box_padding, ymax = F1_ppm + box_padding) %>%
      select(xmin, xmax, ymin, ymax, stain_id, stain_intensity)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FILTER BY KEEP_PEAK_RANGES (same logic as clustering mode)
    # Keeps only top N peaks in each specified F2 range
    # ═══════════════════════════════════════════════════════════════════════════
    if (!is.null(keep_peak_ranges) && is.list(keep_peak_ranges) && length(keep_peak_ranges) > 0) {
      if (verbose) cat("Applying keep_peak_ranges filter (no-clustering mode)...\n")
      
      # DEBUG: Show F2_ppm range and the filter ranges
      if (verbose) {
        cat(sprintf("  DEBUG: peaks F2_ppm range = [%.3f, %.3f]\n", 
                    min(peaks_ppm$F2_ppm), max(peaks_ppm$F2_ppm)))
        cat(sprintf("  DEBUG: %d ranges to apply:\n", length(keep_peak_ranges)))
        for (r in seq_along(keep_peak_ranges)) {
          rng <- keep_peak_ranges[[r]]
          n_in_range <- sum(peaks_ppm$F2_ppm >= min(rng) & peaks_ppm$F2_ppm <= max(rng))
          cat(sprintf("    Range %d: [%.3f, %.3f] -> %d peaks in range\n", 
                      r, min(rng), max(rng), n_in_range))
        }
      }
      
      filtered_peaks <- data.frame()
      filtered_boxes <- data.frame()
      
      for (i in seq_along(keep_peak_ranges)) {
        range <- keep_peak_ranges[[i]]
        
        if (length(range) == 2) {
          lower_bound <- min(range)
          upper_bound <- max(range)
          
          # Filter peaks in this range
          peaks_in_range <- peaks_ppm %>%
            dplyr::filter(F2_ppm >= lower_bound & F2_ppm <= upper_bound)
          
          # Keep fewer peaks in first range (typically reference peak region)
          num_peaks_to_keep <- if (i <= 1) 1 else 4
          
          top_peaks_in_range <- peaks_in_range %>%
            dplyr::arrange(desc(stain_intensity)) %>%
            dplyr::slice_head(n = num_peaks_to_keep)
          
          filtered_peaks <- dplyr::bind_rows(filtered_peaks, top_peaks_in_range)
          
          # Also filter corresponding boxes
          if (nrow(top_peaks_in_range) > 0) {
            boxes_in_range <- boxes %>%
              dplyr::filter(stain_id %in% top_peaks_in_range$stain_id)
            filtered_boxes <- dplyr::bind_rows(filtered_boxes, boxes_in_range)
          }
        }
      }
      
      # Keep peaks outside any specified range
      peaks_outside_ranges <- peaks_ppm %>%
        dplyr::filter(!(
          sapply(1:nrow(peaks_ppm), function(j) {
            any(sapply(keep_peak_ranges, function(range) {
              lower_bound <- min(range)
              upper_bound <- max(range)
              peaks_ppm$F2_ppm[j] >= lower_bound && peaks_ppm$F2_ppm[j] <= upper_bound
            }))
          })
        ))
      
      boxes_outside_ranges <- boxes %>%
        dplyr::filter(stain_id %in% peaks_outside_ranges$stain_id)
      
      peaks_ppm <- dplyr::bind_rows(peaks_outside_ranges, filtered_peaks) %>%
        dplyr::distinct()
      boxes <- dplyr::bind_rows(boxes_outside_ranges, filtered_boxes) %>%
        dplyr::distinct()
      
      if (verbose) {
        cat(sprintf("After keep_peak_ranges filter: %d peaks, %d boxes\n", 
                    nrow(peaks_ppm), nrow(boxes)))
      }
    }
    
    return(list(peaks = peaks_ppm, boxes = boxes, shapes = list()))
  }
  
  # === DBSCAN Clustering ===
  processed <- process_peaks_with_dbscan(peaks_clean_filtered, rr_norm, params, step,
                                         keep_peak_ranges = keep_peak_ranges)
  return(processed)
}