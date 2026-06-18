# CNN_filtering.R - Peak Filtering Functions ----
#
# Part of the CNN Peak Detection module for 2D NMR Spectra
# Contains filtering functions for noise removal and peak cleanup
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB


#' Filter Peaks by Proportion
#'
#' Removes rows and columns that have too many detected peaks (likely noise/artifacts).
#' A high proportion of peaks in a single row/column often indicates systematic noise.
#'
#' @param peaks_clean Data frame with peaks (F1, F2, Intensity columns)
#' @param threshold Numeric, maximum proportion of peaks allowed in a row/column (default: 0.5)
#' @param intensity_threshold Numeric, minimum intensity to keep (optional)
#'
#' @return List containing:
#'   - filtered_peaks: Cleaned peak data frame
#'   - removed_columns: F1 values of removed columns
#'   - removed_rows: F2 values of removed rows
#'
#' @details
#' The function calculates:
#' - peaks_per_col / total_points_per_col for each F1 value
#' - peaks_per_row / total_points_per_row for each F2 value
#' Rows/columns exceeding the threshold are removed.
#'
#' @export
filter_peaks_by_proportion <- function(peaks_clean, threshold = NULL, intensity_threshold = NULL) {
  # ═══════════════════════════════════════════════════════════════════════════
  # VERIFICATION: Ensure there are peaks to process
  # ═══════════════════════════════════════════════════════════════════════════
  if (is.null(peaks_clean) || nrow(peaks_clean) == 0) {
    return(list(
      filtered_peaks = data.frame(F1 = numeric(0), F2 = numeric(0), Intensity = numeric(0)),
      removed_columns = character(0),
      removed_rows = character(0)
    ))
  }
  
  # Optional: initial intensity filtering
  if (!is.null(intensity_threshold)) {
    peaks_clean <- peaks_clean[peaks_clean$Intensity > intensity_threshold, ]
    
    # Vérifier après filtrage intensité
    if (nrow(peaks_clean) == 0) {
      return(list(
        filtered_peaks = data.frame(F1 = numeric(0), F2 = numeric(0), Intensity = numeric(0)),
        removed_columns = character(0),
        removed_rows = character(0)
      ))
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Additional check before tapply
  # ═══════════════════════════════════════════════════════════════════════════
  if (length(unique(peaks_clean$F1)) == 0 || length(unique(peaks_clean$F2)) == 0) {
    return(list(
      filtered_peaks = data.frame(F1 = numeric(0), F2 = numeric(0), Intensity = numeric(0)),
      removed_columns = character(0),
      removed_rows = character(0)
    ))
  }
  
  # Calculate total points per column (F1)
  total_points_per_col <- tryCatch({
    tapply(peaks_clean$F2, peaks_clean$F1, max)
  }, error = function(e) NULL)
  
  if (is.null(total_points_per_col)) {
    return(list(
      filtered_peaks = peaks_clean,
      removed_columns = character(0),
      removed_rows = character(0)
    ))
  }
  
  # Count detected peaks per column
  peaks_per_col <- table(peaks_clean$F1)
  # Calculate proportion of peaks per column
  prop_col <- peaks_per_col / total_points_per_col[names(peaks_per_col)]
  # Identify columns to remove (proportion > threshold)
  cols_to_remove <- if (!is.null(threshold)) names(prop_col[prop_col > threshold]) else character(0)
  
  # Calculate total points per row (F2)
  total_points_per_row <- tryCatch({
    tapply(peaks_clean$F1, peaks_clean$F2, max)
  }, error = function(e) NULL)
  
  if (is.null(total_points_per_row)) {
    return(list(
      filtered_peaks = peaks_clean,
      removed_columns = cols_to_remove,
      removed_rows = character(0)
    ))
  }
  
  # Count detected peaks per row
  peaks_per_row <- table(peaks_clean$F2)
  # Calculate proportion of peaks per row
  prop_row <- peaks_per_row / total_points_per_row[names(peaks_per_row)]
  # Identify rows to remove
  rows_to_remove <- if (!is.null(threshold)) names(prop_row[prop_row > threshold]) else character(0)
  
  # Remove peaks in problematic rows/columns
  filtered_peaks <- peaks_clean[!(peaks_clean$F1 %in% cols_to_remove) &
                                  !(peaks_clean$F2 %in% rows_to_remove), ]
  
  return(list(
    filtered_peaks = filtered_peaks,
    removed_columns = cols_to_remove,
    removed_rows = rows_to_remove
  ))
}

#' Filter Noisy Columns by Relative Intensity
#'
#' Keeps only peaks that have intensity close to the maximum intensity
#' in their respective column. Removes low-amplitude noise.
#'
#' @param peaks_df Data frame with peaks (F1, F2, Intensity columns)
#' @param threshold_ratio Numeric, minimum ratio of peak intensity to column max (default: 0.9)
#' @param min_col_max Numeric, minimum column maximum intensity (optional)
#'
#' @return Filtered data frame with columns F1, F2, Intensity
#'
#' @details
#' A peak is kept if: |Intensity| >= threshold_ratio * MaxIntensity_of_column
#'
#' @export
filter_noisy_columns <- function(peaks_df, threshold_ratio = 0.9, min_col_max = NULL) {
  # ═══════════════════════════════════════════════════════════════════════════
  # VERIFICATION: Ensure there are peaks to process
  # ═══════════════════════════════════════════════════════════════════════════
  if (is.null(peaks_df) || nrow(peaks_df) == 0) {
    return(data.frame(F1 = numeric(0), F2 = numeric(0), Intensity = numeric(0)))
  }
  
  # Calculate maximum absolute intensity per column (F1)
  max_per_col <- tryCatch({
    aggregate(Intensity ~ F1, data = peaks_df, FUN = function(x) max(abs(x)))
  }, error = function(e) {
    warning("filter_noisy_columns: aggregate failed - ", e$message)
    return(NULL)
  })
  
  if (is.null(max_per_col) || nrow(max_per_col) == 0) {
    return(peaks_df[, c("F1", "F2", "Intensity")])
  }
  
  colnames(max_per_col)[2] <- "MaxIntensity"
  
  # Merge max intensity back to each peak
  merged <- merge(peaks_df, max_per_col, by = "F1")
  
  # Filter by relative amplitude
  filtered <- subset(merged, abs(Intensity) >= threshold_ratio * MaxIntensity)
  
  # Optional: remove columns with weak maximum intensity
  if (!is.null(min_col_max)) {
    filtered <- subset(filtered, MaxIntensity >= min_col_max)
  }
  
  # Return to original format
  filtered <- filtered[, c("F1", "F2", "Intensity")]
  return(filtered)
}

#' Clean Peak Clusters with DBSCAN
#'
#' Applies DBSCAN clustering to group nearby peaks. Useful for identifying
#' peak multiplets and removing isolated noise points.
#'
#' @param peaks_df Data frame with peaks (F1, F2 columns)
#' @param ppm_x Numeric vector, F1 ppm axis
#' @param ppm_y Numeric vector, F2 ppm axis
#' @param eps_ppm Numeric, DBSCAN epsilon parameter (default: from params)
#' @param minPts Integer, minimum points per cluster (default: 2)
#'
#' @return Data frame with added 'cluster' column (0 = noise, 1+ = cluster ID)
#'
#' @export
clean_peak_clusters_dbscan <- function(peaks_df, ppm_x, ppm_y, eps_ppm = params$eps_value, minPts = 2) {
  # Remove NA values
  peaks_df <- peaks_df[!is.na(peaks_df$F1) & !is.na(peaks_df$F2), ]
  
  # Extract normalized coordinates for clustering
  coords <- as.matrix(peaks_df[, c("F1", "F2")])
  
  # Apply DBSCAN clustering
  clustering <- dbscan::dbscan(coords, eps = eps_ppm, minPts = minPts)
  peaks_df$cluster <- clustering$cluster
  
  return(peaks_df)
}

#' Remove Peaks in PPM Range
#'
#' Filters out peaks that fall within a specified ppm range on a given axis.
#' Useful for removing solvent peaks, reference peaks, or known artifacts.
#'
#' @param peaks Data frame with peaks
#' @param rr_norm 2D spectrum matrix (used to get ppm axes)
#' @param axis Character, "F1" or "F2" - which axis to filter
#' @param ppm_min Numeric, lower bound of range to remove
#' @param ppm_max Numeric, upper bound of range to remove
#'
#' @return Filtered peaks data frame with peaks in range removed
#'
#' @examples
#' # Remove water region peaks (4.7-5.0 ppm on F1)
#' peaks_clean <- remove_peaks_ppm_range(peaks, spectrum, "F1", 4.7, 5.0)
#'
#' @export
remove_peaks_ppm_range <- function(peaks, rr_norm, axis = "F1", ppm_min, ppm_max) {
  # Get ppm axes from spectrum matrix
  ppm_y <- as.numeric(rownames(rr_norm))   # F2 axis (rows)
  ppm_x <- as.numeric(colnames(rr_norm))   # F1 axis (columns)
  
  # Helper function to convert indices or ppm values to actual ppm
  map_to_ppm <- function(vals, ppm_axis) {
    vnum <- as.numeric(vals)
    n <- length(ppm_axis)
    axis_min <- min(ppm_axis, na.rm = TRUE)
    axis_max <- max(ppm_axis, na.rm = TRUE)
    
    # Check if values are indices (integers in valid range)
    if (all(!is.na(vnum) & vnum >= 1 & vnum <= n & abs(vnum - round(vnum)) < 1e-6)) {
      return(as.numeric(ppm_axis[round(vnum)]))
    }
    # Check if values are already ppm
    if (all(!is.na(vnum) & vnum >= axis_min & vnum <= axis_max)) {
      return(as.numeric(sapply(vnum, function(v) ppm_axis[which.min(abs(ppm_axis - v))])))
    }
    # Mixed case: handle each value individually
    sapply(vnum, function(v) {
      if (is.na(v)) return(NA_real_)
      if (v >= 1 && v <= n && abs(v - round(v)) < 1e-6) {
        return(ppm_axis[round(v)])
      } else if (v >= axis_min && v <= axis_max) {
        return(ppm_axis[which.min(abs(ppm_axis - v))])
      } else {
        return(ppm_axis[which.min(abs(ppm_axis - v))])
      }
    })
  }
  
  # Convert peak positions to ppm
  ppm_vals <- if (axis == "F1") {
    map_to_ppm(peaks$F1, ppm_x)
  } else if (axis == "F2") {
    map_to_ppm(peaks$F2, ppm_y)
  } else {
    stop("axis must be 'F1' or 'F2'")
  }
  
  # Create mask for peaks in the removal range
  mask <- ppm_vals >= ppm_min & ppm_vals <= ppm_max
  removed <- sum(mask, na.rm = TRUE)
  
  message("Removed ", removed, " peaks between ", ppm_min, " and ", ppm_max, " ppm (", axis, ")")
  
  return(peaks[!mask, , drop = FALSE])
}


