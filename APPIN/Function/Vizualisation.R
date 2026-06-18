# 2D NMR Peak Picking and Analysis ----
#
# This script defines functions for visualizing 2D NMR spectra and extracting peak information:
# 1. find_nmr_peak_centroids_optimized: Generates a contour plot from a 2D NMR matrix and extracts contour data.
# 2. get_local_Volume: Computes local contour intensity around a given point.
# 3. Threshold estimation functions: seuil_bruit_multiplicatif, seuil_max_pourcentage, modulate_threshold.
# 4. make_bbox_outline: Creates polygon outlines for bounding box visualization.


## ---- Required libraries ----
library(ggplot2)    # for contour plotting
library(data.table) # for efficient data manipulation
library(dplyr)      # for data wrangling (group_by, summarise, etc.)
library(dbscan)     # for DBSCAN clustering
library(magrittr)   # for piping (%>%)


## ---- Spectra displaying function ----

#' Generate Optimized Contour Plot for 2D NMR Spectrum
#'
#' Creates a contour plot from a 2D NMR intensity matrix with performance optimizations
#' including downsampling, early intensity filtering, and spectrum-type-specific defaults.
#'
#' @param rr_data Numeric matrix. 2D NMR intensity data with chemical shifts as row/column names.
#' @param spectrum_type Character. Type of NMR experiment: "HSQC", "TOCSY", "COSY", or "UFCOSY".
#'   If provided, default parameters will be loaded for that experiment type.
#' @param contour_start Numeric. Starting intensity level for contour lines. Default depends on spectrum_type.
#' @param intensity_threshold Numeric. Minimum intensity to include in the plot (filters noise).
#' @param contour_num Integer. Number of contour levels to draw.
#' @param contour_factor Numeric. Multiplicative factor between successive contour levels (geometric progression).
#' @param zoom_xlim Numeric vector of length 2. X-axis limits (F2/ppm) for zooming.
#' @param zoom_ylim Numeric vector of length 2. Y-axis limits (F1/ppm) for zooming.
#' @param f2_exclude_range Numeric vector of length 2. F2 ppm range to exclude (e.g., water region).
#' @param downsample_factor Integer. Factor for matrix downsampling (default 2). Higher values = faster but less detail.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{plot} - A ggplot2 contour plot object
#'     \item \code{contour_data} - Data frame with contour line coordinates extracted via ggplot_build()
#'   }
#'
#' @details
#' The function applies several optimizations for handling large NMR matrices:
#' \enumerate{
#'   \item Downsampling: Reduces matrix size by the specified factor
#'   \item Early filtering: Only processes points above intensity_threshold
#'   \item Direct data.table creation: Avoids memory-intensive expand.grid()
#'   \item Contour level limiting: Caps TOCSY spectra at 20 levels maximum
#' }
#'
#' Default parameters by spectrum type:
#' \itemize{
#'   \item HSQC: contour_start=8000, intensity_threshold=200, contour_num=8
#'   \item TOCSY: contour_start=100000, intensity_threshold=4000, contour_num=40, excludes water (4.7-5.0 ppm)
#'   \item COSY: contour_start=1000, intensity_threshold=20000, contour_num=60
#'   \item UFCOSY: contour_start=1000, intensity_threshold=20000, contour_num=40
#' }
#'
#' @examples
#' \dontrun{
#' # Load spectrum and generate contour plot
#' result <- find_nmr_peak_centroids_optimized(
#'   rr_data = spectrum$spectrumData,
#'   spectrum_type = "TOCSY",
#'   downsample_factor = 2
#' )
#'
#' # Display the plot
#' print(result$plot)
#'
#' # Access contour coordinates for further processing
#' contours <- result$contour_data
#' }
#'
#' @export

find_nmr_peak_centroids_optimized <- function(rr_data, spectrum_type = NULL, 
                                              contour_start = NULL, intensity_threshold = NULL, 
                                              contour_num = NULL, contour_factor = NULL, 
                                              zoom_xlim = NULL, zoom_ylim = NULL, 
                                              f2_exclude_range = NULL,
                                              downsample_factor = 2) {  # New parameter
  
  # --- Input validation ---
  if (is.null(rr_data) || !is.matrix(rr_data)) {
    stop("Invalid Bruker data. Ensure rr_data is a matrix with proper intensity values.")
  }
  
  # Record which parameters the caller actually supplied (before defaults fill
  # them in). Used by the HMBC noise-derivation so it only kicks in for
  # parameters the user did NOT set -- otherwise it would silently override the
  # UI and changing the settings would have no effect.
  user_contour_start      <- contour_start
  user_intensity_threshold <- intensity_threshold
  user_contour_num        <- contour_num
  user_contour_factor     <- contour_factor
  
  # --- Default parameters for each spectrum type ---
  # These values are empirically determined for typical NMR experiments
  spectrum_defaults <- list(
    HSQC = list(contour_start = 8000, intensity_threshold = 200, contour_num = 8, contour_factor = 1.3),
    TOCSY = list(contour_start = 100000, intensity_threshold = 4000, contour_num = 40, contour_factor = 1.3, f2_exclude_range = c(4.7, 5.0)),
    UFCOSY = list(contour_start = 1000, intensity_threshold = 20000, contour_num = 40, contour_factor = 1.3),
    COSY = list(contour_start = 1000, intensity_threshold = 20000, contour_num = 60, contour_factor = 1.3),
    # HMBC: heteronuclear, signals much weaker than HSQC (2J/3J correlations).
    # Lower thresholds so real cross-peaks are not all filtered out, and a
    # modest contour count to keep the contour grid affordable on the wide 13C axis.
    HMBC = list(contour_start = 5000, intensity_threshold = 150, contour_num = 12, contour_factor = 1.3),
    # J-RES: homonuclear; F2 = 1H shift, F1 = J coupling (narrow). Thresholds
    # similar to TOCSY/COSY; signals are strong, no special weakening needed.
    JRES = list(contour_start = 80000, intensity_threshold = 4000, contour_num = 40, contour_factor = 1.3, f2_exclude_range = c(4.7, 5.0))
  )
  
  # Load defaults if spectrum_type is specified, allowing individual parameter overrides
  if (!is.null(spectrum_type)) {
    if (!spectrum_type %in% names(spectrum_defaults)) {
      # Unknown type: warn and fall back to TOCSY-like defaults instead of stop().
      # A stop() here is caught by the caller's tryCatch, but an unknown type
      # should degrade gracefully rather than abort processing.
      warning(sprintf("Unknown spectrum_type '%s' - falling back to generic defaults.", spectrum_type))
      defaults <- list(contour_start = 100000, intensity_threshold = 4000,
                       contour_num = 20, contour_factor = 1.3)
    } else {
      defaults <- spectrum_defaults[[spectrum_type]]
    }
    contour_start <- ifelse(is.null(contour_start), defaults$contour_start, contour_start)
    intensity_threshold <- ifelse(is.null(intensity_threshold), defaults$intensity_threshold, intensity_threshold)
    contour_num <- ifelse(is.null(contour_num), defaults$contour_num, contour_num)
    contour_factor <- ifelse(is.null(contour_factor), defaults$contour_factor, contour_factor)
    
    # HMBC: signals are weak and vary a lot between samples, so a fixed
    # contour_start often misses every peak (0 peaks). Derive sensible values
    # from the actual noise level of this spectrum -- but ONLY for parameters
    # the caller did not provide. If the UI passed a value, respect it (this is
    # why changing the settings previously had no effect: the derivation used
    # to overwrite the UI values unconditionally).
    if (identical(spectrum_type, "HMBC")) {
      noise_sd <- sd(as.numeric(rr_data), na.rm = TRUE)
      max_abs  <- max(abs(rr_data), na.rm = TRUE)
      
      # user_* captured at function entry == NULL means the caller (UI) did not
      # supply that parameter, so we derive it from the noise; otherwise respect it.
      if (is.finite(noise_sd) && noise_sd > 0) {
        if (is.null(user_contour_start)) {
          contour_start <- min(noise_sd * 7, max_abs * 0.25)
        }
        if (is.null(user_intensity_threshold)) {
          intensity_threshold <- noise_sd * 4
        }
      }
      if (is.null(user_contour_num))    contour_num    <- max(contour_num, 18)
      if (is.null(user_contour_factor)) contour_factor <- 1.25
    }
  }
  
  # --- Guard against NULL/NA/non-finite numeric parameters ---
  # If any of these are NULL or NA (e.g. an unset reactive on a new spectrum
  # type), downstream calls like which(rr_data >= NULL) or seq()/contour level
  # construction can produce huge or invalid results and crash the R process.
  if (is.null(intensity_threshold) || is.na(intensity_threshold) || !is.finite(intensity_threshold)) {
    intensity_threshold <- stats::quantile(abs(rr_data), 0.999, na.rm = TRUE)
    warning("intensity_threshold was invalid; using 99.9th percentile of |intensity| as fallback.")
  }
  if (is.null(contour_start) || is.na(contour_start) || !is.finite(contour_start)) {
    contour_start <- intensity_threshold
  }
  if (is.null(contour_num) || is.na(contour_num) || contour_num < 1) contour_num <- 12
  if (is.null(contour_factor) || is.na(contour_factor) || contour_factor <= 1) contour_factor <- 1.3
  
  # OPTIMIZATION 1: Downsampling the matrix for display
  # Reduces memory usage and rendering time for large spectra.
  # HMBC: the indirect (13C) dimension is acquired with FEW increments, so
  # downsampling makes contours jagged/stair-stepped. Skip it for HMBC.
  effective_downsample <- if (identical(spectrum_type, "HMBC")) 1 else downsample_factor
  if (effective_downsample > 1) {
    seq_x <- seq(1, nrow(rr_data), by = effective_downsample)
    seq_y <- seq(1, ncol(rr_data), by = effective_downsample)
    rr_data <- rr_data[seq_x, seq_y]
  }
  
  # Extract chemical shift axes from matrix row/column names
  ppm_x <- as.numeric(rownames(rr_data))
  ppm_y <- as.numeric(colnames(rr_data))
  
  # OPTIMIZATION 2: Early filtering before expand.grid
  # Only keep matrix indices where intensity exceeds threshold
  # This dramatically reduces the data size for sparse spectra
  high_intensity_indices <- which(rr_data >= intensity_threshold, arr.ind = TRUE)
  
  if (nrow(high_intensity_indices) == 0) {
    warning("No data points above intensity threshold")
    return(list(plot = ggplot() + theme_void(), contour_data = data.frame()))
  }
  
  # --- ANTI-OOM GUARD ---
  # If the threshold is too low for this spectrum (common on HMBC, whose
  # signals are weak so users/defaults may under-set the threshold), the
  # number of points above threshold and the resulting contour grid can blow
  # up memory and the OS kills the R process (the app "just closes"). Cap the
  # number of retained points by raising the threshold to keep the contour
  # build affordable, rather than crashing.
  max_points <- 400000L
  n_above <- nrow(high_intensity_indices)
  if (n_above > max_points) {
    vals <- rr_data[high_intensity_indices]
    # keep the strongest max_points points -> threshold = corresponding quantile
    keep_frac <- max_points / n_above
    adj_threshold <- stats::quantile(vals, probs = 1 - keep_frac, na.rm = TRUE)
    warning(sprintf(
      "Too many points above threshold (%d > %d) for %s; raising threshold to %.4g to avoid memory overflow.",
      n_above, max_points, ifelse(is.null(spectrum_type), "spectrum", spectrum_type), adj_threshold))
    high_intensity_indices <- which(rr_data >= adj_threshold, arr.ind = TRUE)
    if (contour_start < adj_threshold) contour_start <- adj_threshold
  }
  
  # OPTIMIZATION 3: Direct creation of the data.frame without expand.grid
  # Using data.table for memory efficiency with large datasets.
  #
  # NOTE: geom_contour requires a COMPLETE regular grid to interpolate smooth
  # lines. Feeding it only the above-threshold points (a holey point cloud)
  # makes it interpolate across gaps -> jagged/torn contours. This "sparse"
  # path is fine for dense spectra (TOCSY) but breaks on sparse HMBC.
  # For HMBC we therefore build the FULL grid; the threshold only sets which
  # contour levels are drawn (via breaks = contour_levels), not which points
  # exist. The sparse path is preserved for all other spectrum types.
  if (identical(spectrum_type, "HMBC")) {
    intensity_df <- data.table(
      ppm_x = rep(ppm_x, times = length(ppm_y)),
      ppm_y = rep(ppm_y, each  = length(ppm_x)),
      intensity = as.numeric(rr_data)
    )
    # Zero-out the water region instead of deleting rows, so the grid stays
    # rectangular (deleting rows would re-introduce holes -> jagged contours).
    if (!is.null(f2_exclude_range) && length(f2_exclude_range) == 2) {
      intensity_df[ppm_y >= f2_exclude_range[1] & ppm_y <= f2_exclude_range[2],
                   intensity := 0]
    }
  } else {
    intensity_df <- data.table(
      ppm_x = ppm_x[high_intensity_indices[, 1]],
      ppm_y = ppm_y[high_intensity_indices[, 2]],
      intensity = rr_data[high_intensity_indices]
    )
    
    # Exclusion of the water region (typically 4.7-5.0 ppm in F2)
    # Water signal creates artifacts that interfere with peak detection
    if (!is.null(f2_exclude_range) && length(f2_exclude_range) == 2) {
      intensity_df <- intensity_df[!(ppm_y >= f2_exclude_range[1] & ppm_y <= f2_exclude_range[2])]
    }
  }
  
  # OPTIMIZATION 4: Reduce the number of contour levels for TOCSY
  # TOCSY spectra with too many contours can be slow to render
  if (!is.null(spectrum_type) && spectrum_type == "TOCSY" && contour_num > 20) {
    contour_num <- min(contour_num, 20)  # Limit to a maximum of 20 outlines
  }
  
  # Calculate contour levels as geometric progression
  # Each level is contour_factor times the previous one
  contour_levels <- contour_start * contour_factor^(0:(contour_num - 1))
  
  # OPTIMIZATION 5: Build the ggplot contour plot
  # Axes are reversed to match NMR convention (high ppm on left/top)
  p <- ggplot(intensity_df, aes(x = ppm_y, y = ppm_x, z = intensity)) +
    geom_contour(color = "black", breaks = contour_levels, linewidth = 0.3) +  # linewidth reduced
    scale_x_reverse() +
    scale_y_reverse() +
    labs(y = "F1_ppm", x = "F2_ppm") +
    theme_minimal() +
    theme(
      panel.grid = element_blank()  # Remove the grid to speed up rendering
    )
  
  # Apply zoom limits if specified (improves performance by reducing displayed data)
  if (!is.null(zoom_xlim)) p <- p + coord_cartesian(xlim = zoom_xlim)
  if (!is.null(zoom_ylim)) p <- p + coord_cartesian(ylim = zoom_ylim)
  
  # Extract contour line coordinates from the built plot
  # This data is used for subsequent peak detection and clustering
  contour_data <- ggplot_build(p)$data[[1]]
  
  return(list(plot = p, contour_data = contour_data))
}


# # ===== BONUS FUNCTION: Visualize rejected clusters
# 
# plot_rejected_clusters <- function(rr_data, process_result) {
#   
#   # Récupérer les stats
#   stats <- process_result$cluster_stats
#   valid_ids <- process_result$centroids$stain_id
#   
#   # Identifier les rejetés
#   stats <- stats %>%
#     mutate(status = ifelse(stain_id %in% valid_ids, "Valid", "Rejected"))
#   
#   # Plot
#   p <- ggplot(stats, aes(x = elongation, y = intensity_norm, color = status)) +
#     geom_point(size = 3, alpha = 0.7) +
#     geom_vline(xintercept = c(5, 8, 10, 15), linetype = "dashed", alpha = 0.3) +
#     scale_color_manual(values = c("Valid" = "green", "Rejected" = "red")) +
#     scale_x_log10() +
#     labs(
#       title = "Diagnostic: Clusters valides vs rejetés",
#       x = "Élongation (log scale)",
#       y = "Intensité normalisée",
#       color = "Status"
#     ) +
#     theme_minimal()
#   
#   print(p)
#   
#   # Tableau des rejetés
#   rejected <- stats %>% 
#     filter(status == "Rejected") %>%
#     arrange(desc(intensity))
#   
#   cat("\n=== Top 10 clusters rejetés ===\n")
#   print(rejected %>% 
#           select(stain_id, x_center, y_center, intensity, elongation, aspect_ratio, density) %>%
#           head(10))
#   
#   invisible(stats)
# }


#' Compute Local Contour Volume Around a Point
#'
#' Calculates the sum of contour levels within a small neighborhood around
#' a specified chemical shift coordinate. Used to estimate local signal intensity.
#'
#' @param f2_ppm Numeric. F2 (direct dimension) chemical shift in ppm.
#' @param f1_ppm Numeric. F1 (indirect dimension) chemical shift in ppm.
#' @param contour_data Data frame. Contour data from ggplot_build(), must contain x, y, and level columns.
#' @param eps_ppm Numeric. Half-width of the neighborhood in ppm (default 0.0068 ppm ≈ 4 Hz at 600 MHz).
#'
#' @return Numeric. Sum of contour levels within the neighborhood, or NA if no points found.
#'
#' @details
#' This function defines a square neighborhood of size (2*eps_ppm) x (2*eps_ppm)
#' centered on the specified coordinates and sums all contour level values within.
#' Higher values indicate stronger signals at that position.
#'
#' @examples
#' \dontrun
#' # Get local volume at a specific peak position
#' volume <- get_local_Volume(
#'   f2_ppm = 3.45,
#'   f1_ppm = 1.23,
#'   contour_data = result$contour_data
#' )
#' }
#'
#' @export

get_local_Volume <- function(f2_ppm, f1_ppm, contour_data, eps_ppm = 0.0068) {
  # Filter contour points within a square neighborhood centered on (f2_ppm, f1_ppm)
  # The neighborhood size is 2*eps_ppm in each dimension
  local_points <- contour_data %>%
    filter(
      x >= f2_ppm - eps_ppm & x <= f2_ppm + eps_ppm,
      y >= f1_ppm - eps_ppm & y <= f1_ppm + eps_ppm
    )
  
  # Return NA if no contour points found in the neighborhood
  if (nrow(local_points) == 0) {
    return(NA)
  }
  
  # Sum all contour levels as a proxy for local signal volume/intensity
  return(sum(local_points$level, na.rm = TRUE))
}



## Noise threshold ----

#' Estimate Noise Threshold Using Standard Deviation
#'
#' Calculates a signal threshold based on the standard deviation of the intensity matrix,
#' multiplied by a user-defined factor. Assumes noise follows a Gaussian distribution.
#'
#' @param mat Numeric matrix. 2D NMR intensity data.
#' @param facteur Numeric. Multiplicative factor for the standard deviation (default 3).
#'   A factor of 3 corresponds to ~99.7% confidence for Gaussian noise.
#'
#' @return Numeric. Estimated noise threshold.
#'
#' @examples
#' \dontrun{
#' threshold <- seuil_bruit_multiplicatif(spectrum_matrix, facteur = 3)
#' }
#'
#' @export

seuil_bruit_multiplicatif <- function(mat, facteur = 3) {
  # Estimate noise level as the standard deviation of all intensity values
  bruit_estime <- sd(as.numeric(mat), na.rm = TRUE)
  # Threshold = noise_estimate * factor (typically 3 for 3-sigma rule)
  seuil <- bruit_estime * facteur
  return(seuil)
}


#' Estimate Threshold as Percentage of Maximum Intensity
#'
#' Calculates a signal threshold as a fixed percentage of the maximum intensity
#' in the spectrum. Simple approach that adapts to the overall signal strength.
#'
#' @param mat Numeric matrix. 2D NMR intensity data.
#' @param pourcentage Numeric. Fraction of maximum intensity (default 0.05 = 5%).
#'
#' @return Numeric. Estimated threshold value.
#'
#' @examples
#' \dontrun{
#' threshold <- seuil_max_pourcentage(spectrum_matrix, pourcentage = 0.05)
#' }
#'
#' @export

seuil_max_pourcentage <- function(mat, pourcentage = 0.05) {
  # Find the maximum intensity value in the matrix
  max_val <- max(mat, na.rm = TRUE)
  # Threshold = percentage of maximum
  seuil <- max_val * pourcentage
  return(seuil)
}


#' Modulate Threshold Based on Volume Integral
#'
#' Applies a power-law modulation to adjust thresholds based on peak volume.
#' Used for adaptive thresholding where larger peaks need different cutoffs.
#'
#' @param VI Numeric. Volume integral or intensity value.
#'
#' @return Numeric. Modulated threshold value.
#'
#' @details
#' Uses the formula: threshold = a * VI^b, where a=0.0006 and b=1.2.
#' These empirical constants were determined through testing on typical NMR spectra.
#' The power-law relationship accounts for the non-linear relationship between
#' peak intensity and appropriate threshold values.
#'
#' @examples
#' \dontrun{
#' modulated <- modulate_threshold(peak_volume)
#' }
#'
#' @export

modulate_threshold <- function(VI) {
  # Empirical constants for power-law modulation
  a <- 0.0006
  b <- 1.2
  # Power-law relationship: larger VI values get proportionally higher thresholds
  a * VI^b
}


#' Create Bounding Box Outlines for Plotting
#'
#' Converts a data frame of bounding boxes into a format suitable for
#' plotting with geom_path() in ggplot2. Each box becomes a closed polygon.
#'
#' @param boxes Data frame. Must contain columns: xmin, xmax, ymin, ymax.
#'   Optionally includes stain_id for identification.
#'
#' @return Data frame with columns x, y, group for use with geom_path().
#'   Returns NULL if input is NULL or empty.
#'
#' @details
#' Each bounding box is converted to 5 points (closed rectangle) with a unique
#' group identifier. NA values are inserted between boxes to prevent connecting
#' lines when plotting multiple boxes with a single geom_path() call.
#'
#' @examples
#' \dontrun{
#' # Create box outlines for plotting
#' outlines <- make_bbox_outline(peak_boxes)
#'
#' # Add to existing plot
#' p + geom_path(data = outlines, aes(x = x, y = y, group = group), color = "red")
#' }
#'
#' @export

make_bbox_outline <- function(boxes) {
  # Handle NULL or empty input
  if (is.null(boxes) || nrow(boxes) == 0) return(NULL)
  
  # Ensure that stain_id exists for identification
  if (!"stain_id" %in% names(boxes)) {
    boxes$stain_id <- paste0("box_", seq_len(nrow(boxes)))
  }
  
  # Create outlines with UNIQUE group and NA between each box
  # NA values prevent geom_path from connecting separate boxes
  outline_list <- lapply(seq_len(nrow(boxes)), function(i) {
    box <- boxes[i, ]
    
    # Create a closed rectangle (5 points: 4 corners + return to start)
    rect <- data.frame(
      x = c(box$xmin, box$xmax, box$xmax, box$xmin, box$xmin),
      y = c(box$ymin, box$ymin, box$ymax, box$ymax, box$ymin),
      group = paste0("box_", i),  # Unique group per box
      stringsAsFactors = FALSE
    )
    
    # Add a row with NA to separate boxes (prevents connecting lines)
    if (i < nrow(boxes)) {
      rect <- rbind(rect, data.frame(
        x = NA,
        y = NA,
        group = paste0("box_", i),
        stringsAsFactors = FALSE
      ))
    }
    
    rect
  })
  
  # Combine all box outlines into a single data frame
  do.call(rbind, outline_list)
}