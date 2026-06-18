# CNN_clustering.R - DBSCAN Clustering and Bounding Box Generation ----
#
# Part of the CNN Peak Detection module for 2D NMR Spectra
# Contains DBSCAN clustering and bounding box functions
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB


#' Process Peaks with DBSCAN and Generate Bounding Boxes
#'
#' Takes filtered peaks and applies DBSCAN clustering to group them,
#' then generates bounding boxes for each cluster suitable for Plotly visualization.
#'
#' @param peaks_clean_filtered Data frame with filtered peaks (F1, F2, Intensity)
#' @param rr_norm 2D spectrum matrix
#' @param params List of parameters including eps_value for DBSCAN
#' @param step Integer, downsampling factor for visualization (default: 4)
#' @param keep_peak_ranges List of numeric vectors. Specific F2 ranges where only top peaks are kept.
#'
#' @return List containing:
#'   - peaks: Data frame with cluster centroids (F1, F2, F1_ppm, F2_ppm, stain_intensity, cluster_db)
#'   - boxes: Data frame with bounding box coordinates in ppm
#'   - shapes: List of Plotly shape objects for overlay
#'
#' @details
#' Processing steps:
#' 1. Z-score normalization of peak coordinates
#' 2. DBSCAN clustering with specified epsilon
#' 3. Bounding box calculation for each cluster
#' 4. Conversion from indices to ppm values
#' 5. Generation of Plotly-compatible shape objects
#' 6. Optional filtering by F2 ranges (keep_peak_ranges)
#'
#' @export
process_peaks_with_dbscan <- function(peaks_clean_filtered, rr_norm, params, step = 4,
                                      keep_peak_ranges = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # VERIFICATION: Ensure there are peaks to process
  # ═══════════════════════════════════════════════════════════════════════════
  if (is.null(peaks_clean_filtered) || nrow(peaks_clean_filtered) == 0) {
    warning("process_peaks_with_dbscan: Aucun pic à traiter")
    return(list(
      peaks = data.frame(F1 = numeric(0), F2 = numeric(0), F1_ppm = numeric(0), 
                         F2_ppm = numeric(0), stain_intensity = numeric(0), cluster_db = integer(0)),
      boxes = data.frame(xmin_ppm = numeric(0), xmax_ppm = numeric(0), 
                         ymin_ppm = numeric(0), ymax_ppm = numeric(0), 
                         stain_intensity = numeric(0), cluster_db = integer(0)),
      shapes = list()
    ))
  }
  
  # Helper function to downsample matrix for visualization
  downsample_matrix <- function(mat, step = 4) {
    mat[seq(1, nrow(mat), by = step),
        seq(1, ncol(mat), by = step)]
  }
  
  # Helper function to downsample axis values
  downsample_axis <- function(axis_vals, step = 4) {
    axis_vals[seq(1, length(axis_vals), by = step)]
  }
  
  # --- Extract matrix and axes ---
  # CORRECTED: Match Peak_picking.R convention
  # In NMR 2D spectra (TOCSY/COSY):
  #   rownames = F2 ppm (direct dimension, typically 0-10 ppm for 1H)
  #   colnames = F1 ppm (indirect dimension)
  # CNN_detection uses F2 for row index, F1 for col index
  # But we need to SWAP to match Peak_picking.R output
  z_matrix <- rr_norm
  x_vals <- as.numeric(colnames(rr_norm))  # For F2_ppm (after swap)
  y_vals <- as.numeric(rownames(rr_norm))  # For F1_ppm (after swap)
  
  # Downsample for faster visualization
  z_small <- downsample_matrix(z_matrix, step = step)
  x_small <- downsample_axis(x_vals, step = step)
  y_small <- downsample_axis(y_vals, step = step)
  
  # --- Step 1: Convert the F1/F2 indices to ppm coordinates ---
  # SWAP: F1 from CNN becomes F2_ppm, F2 from CNN becomes F1_ppm
  # This matches Peak_picking.R convention where F2 is the filtering axis
  peaks_with_ppm <- peaks_clean_filtered %>%
    mutate(
      F2_idx = pmin(pmax(round(F1), 1), length(x_vals)),  # CNN F1 -> our F2
      F1_idx = pmin(pmax(round(F2), 1), length(y_vals)),  # CNN F2 -> our F1
      F2_ppm = x_vals[F2_idx],  # F2 from colnames
      F1_ppm = y_vals[F1_idx]   # F1 from rownames
    )
  
  
  # Delete the lines with NA (just in case)
  peaks_with_ppm <- peaks_with_ppm %>%
    dplyr::filter(!is.na(F1_ppm) & !is.na(F2_ppm))
  
  if (nrow(peaks_with_ppm) == 0) {
    warning("process_peaks_with_dbscan: No valid peak after ppm conversion")
    return(list(
      peaks = data.frame(F2_ppm = numeric(0), F1_ppm = numeric(0), 
                         stain_intensity = numeric(0), cluster_db = integer(0),
                         stain_id = character(0)),
      boxes = data.frame(xmin = numeric(0), xmax = numeric(0), 
                         ymin = numeric(0), ymax = numeric(0), 
                         stain_intensity = numeric(0), stain_id = character(0)),
      shapes = list()
    ))
  }
  
  # --- Step 2: DBSCAN clustering on PPM coordinates (like Local Max) ---
  
  eps_cnn <- params$eps_value * 5  
  
  cat(sprintf("DBSCAN avec eps_cnn = %.4f ppm (base eps = %.4f × 5)\n", 
              eps_cnn, params$eps_value))
  
  db <- dbscan::dbscan(peaks_with_ppm[, c("F1_ppm", "F2_ppm")],
                       eps = eps_cnn, minPts = 1)
  peaks_with_ppm <- peaks_with_ppm %>%
    mutate(cluster_db = db$cluster)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # VERIFICATION: Ensure that there are valid clusters (cluster_db > 0)
  # ═══════════════════════════════════════════════════════════════════════════
  peaks_with_clusters <- peaks_with_ppm %>% dplyr::filter(cluster_db > 0)
  
  if (nrow(peaks_with_clusters) == 0) {
    warning("process_peaks_with_dbscan: No clusters found (all peaks are noise)")
    return(list(
      peaks = peaks_with_ppm %>% 
        dplyr::transmute(F2_ppm, F1_ppm, stain_intensity = Intensity, 
                         cluster_db = 0, stain_id = paste0("cnn_noise_", row_number())),
      boxes = data.frame(xmin = numeric(0), xmax = numeric(0), 
                         ymin = numeric(0), ymax = numeric(0), 
                         stain_intensity = numeric(0), stain_id = character(0)),
      shapes = list()
    ))
  }
  
  # --- Step 3: Calculate bounding boxes for each cluster (en ppm) ---
  # Convention NMR: xmin/xmax = F2 (horizontal), ymin/ymax = F1 (vertical)
  bounding_boxes <- peaks_with_clusters %>%
    group_by(cluster_db) %>%
    summarise(
      xmin = min(F2_ppm, na.rm = TRUE),  # F2 = horizontal 
      xmax = max(F2_ppm, na.rm = TRUE),
      ymin = min(F1_ppm, na.rm = TRUE),  # F1 = vertical 
      ymax = max(F1_ppm, na.rm = TRUE),
      center_F2 = (min(F2_ppm, na.rm = TRUE) + max(F2_ppm, na.rm = TRUE)) / 2,
      center_F1 = (min(F1_ppm, na.rm = TRUE) + max(F1_ppm, na.rm = TRUE)) / 2,
      intensity = sum(Intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Add padding around each box to encompass the entire area
  # Uses the box_padding parameter if provided, otherwise calculates automatically
  # ═══════════════════════════════════════════════════════════════════════════
  padding_ppm <- if (!is.null(params$box_padding)) params$box_padding else (eps_cnn * 1.5)
  
  cat(sprintf("Box padding = %.4f ppm\n", padding_ppm))
  
  bounding_boxes <- bounding_boxes %>%
    mutate(
      xmin = xmin - padding_ppm,
      xmax = xmax + padding_ppm,
      ymin = ymin - padding_ppm,
      ymax = ymax + padding_ppm
    )
  
  # --- Step 4 ---
  bounding_boxes <- bounding_boxes %>%
    mutate(
      stain_intensity = intensity
    )
  
  # --- Step 5: Generate Plotly rectangle shapes ---
  shapes_list <- lapply(seq_len(nrow(bounding_boxes)), function(i) {
    list(
      type = "rect",
      x0 = bounding_boxes$xmin[i], x1 = bounding_boxes$xmax[i],
      y0 = bounding_boxes$ymin[i], y1 = bounding_boxes$ymax[i],
      line = list(color = "red", dash = "dot", width = 2),
      fillcolor = "rgba(0,0,0,0)",  # Transparent fill
      xref = "x", yref = "y",
      layer = "above"
    )
  })
  
  # --- Step 6: Build peaks based on bounding box centers ---
  
  peaks_from_boxes <- bounding_boxes %>%
    dplyr::transmute(
      F2_ppm = center_F2,  # horizontal
      F1_ppm = center_F1,  # vertical
      stain_intensity = stain_intensity,
      cluster_db = cluster_db,
      stain_id = paste0("cnn_", cluster_db)
    )
  
  # Reformat the boxes to be compatible with the app
  # The app expects: xmin, xmax, ymin, ymax, stain_id
  boxes_formatted <- bounding_boxes %>%
    dplyr::transmute(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      stain_id = paste0("cnn_", cluster_db),
      stain_intensity = stain_intensity
    )
  
  # === Step 8: Filter by specific F2 ranges (if specified) ===
  # Keeps only top N peaks in each specified range (same logic as Peak_picking.R)
  if (!is.null(keep_peak_ranges) && is.list(keep_peak_ranges) && length(keep_peak_ranges) > 0) {
    cat("Applying keep_peak_ranges filter...\n")
    
    # DEBUG: Show F2_ppm range and the filter ranges
    cat(sprintf("  DEBUG: peaks F2_ppm range = [%.3f, %.3f]\n", 
                min(peaks_from_boxes$F2_ppm), max(peaks_from_boxes$F2_ppm)))
    cat(sprintf("  DEBUG: %d ranges to apply:\n", length(keep_peak_ranges)))
    for (r in seq_along(keep_peak_ranges)) {
      rng <- keep_peak_ranges[[r]]
      n_in_range <- sum(peaks_from_boxes$F2_ppm >= min(rng) & peaks_from_boxes$F2_ppm <= max(rng))
      cat(sprintf("    Range %d: [%.3f, %.3f] -> %d peaks in range\n", 
                  r, min(rng), max(rng), n_in_range))
    }
    
    filtered_peaks <- data.frame()
    filtered_boxes <- data.frame()
    
    for (i in seq_along(keep_peak_ranges)) {
      range <- keep_peak_ranges[[i]]
      
      if (length(range) == 2) {
        lower_bound <- min(range)
        upper_bound <- max(range)
        
        # Filter peaks in this range
        peaks_in_range <- peaks_from_boxes %>%
          dplyr::filter(F2_ppm >= lower_bound & F2_ppm <= upper_bound)
        
        # Keep fewer peaks in first range (typically reference peak region)
        num_peaks_to_keep <- if (i <= 1) 1 else 4
        
        top_peaks_in_range <- peaks_in_range %>%
          dplyr::arrange(desc(stain_intensity)) %>%
          dplyr::slice_head(n = num_peaks_to_keep)
        
        filtered_peaks <- dplyr::bind_rows(filtered_peaks, top_peaks_in_range)
        
        # Also filter corresponding boxes
        if (nrow(top_peaks_in_range) > 0) {
          boxes_in_range <- boxes_formatted %>%
            dplyr::filter(stain_id %in% top_peaks_in_range$stain_id)
          filtered_boxes <- dplyr::bind_rows(filtered_boxes, boxes_in_range)
        }
      }
    }
    
    # Keep peaks outside any specified range
    peaks_outside_ranges <- peaks_from_boxes %>%
      dplyr::filter(!(
        sapply(1:nrow(peaks_from_boxes), function(j) {
          any(sapply(keep_peak_ranges, function(range) {
            lower_bound <- min(range)
            upper_bound <- max(range)
            peaks_from_boxes$F2_ppm[j] >= lower_bound && peaks_from_boxes$F2_ppm[j] <= upper_bound
          }))
        })
      ))
    
    boxes_outside_ranges <- boxes_formatted %>%
      dplyr::filter(stain_id %in% peaks_outside_ranges$stain_id)
    
    peaks_from_boxes <- dplyr::bind_rows(peaks_outside_ranges, filtered_peaks) %>%
      dplyr::distinct()
    boxes_formatted <- dplyr::bind_rows(boxes_outside_ranges, filtered_boxes) %>%
      dplyr::distinct()
    
    cat(sprintf("After keep_peak_ranges filter: %d peaks, %d boxes\n", 
                nrow(peaks_from_boxes), nrow(boxes_formatted)))
  }
  
  # --- Step 9: Return results ---
  return(list(
    peaks = peaks_from_boxes,
    boxes = boxes_formatted
    # shapes = shapes_list
  ))
}