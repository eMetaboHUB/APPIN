# 2D NMR Analyst - Utility Functions ----

# Author: Julien Guibert
# Description: Shared utility functions for the 2D NMR analysis application


#' Null-coalesce operator
#' 
#' Returns the first argument if it is not NULL, otherwise returns the second.
#' Similar to the ?? operator in C# or the // operator in Perl.
#'
#' @param a First value to check
#' @param b Default value if a is NULL
#' @return a if not NULL, otherwise b
#' @examples
#' NULL %||% "default"  # returns "default"
#' "value" %||% "default"  # returns "value"


`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# Spectrum  ----

# Global environment for caching loaded Bruker spectra
# This avoids re-reading large files multiple times
spectrum_cache <- new.env(parent = emptyenv())

#' Read Bruker spectrum with caching
#' 
#' Reads a Bruker NMR spectrum from disk, caching the result for subsequent
#' calls with the same path. This significantly improves performance when
#' the same spectrum is accessed multiple times.
#'
#' @param path Character. Path to the Bruker spectrum directory
#' @param dim Character. Dimension of the spectrum, either "1D" or "2D" (default: "2D")
#' @return List containing the spectrum data (from read_bruker function)
#' @note Requires the read_bruker() function from Read_2DNMR_spectrum.R
read_bruker_cached <- function(path, dim = "2D") {
  key <- normalizePath(path, mustWork = FALSE)
  
  if (exists(key, envir = spectrum_cache, inherits = FALSE)) {
    return(get(key, envir = spectrum_cache, inherits = FALSE))
  }
  
  data <- read_bruker(path, dim = dim)
  assign(key, data, envir = spectrum_cache)
  data
}

#' Clear the spectrum cache
#' 
#' Removes all cached spectra from memory. Useful when memory is low
#' or when spectra files have been modified on disk.
#'
#' @return NULL (invisible)
clear_spectrum_cache <- function() {
  rm(list = ls(envir = spectrum_cache), envir = spectrum_cache)
  invisible(NULL)
}


# Box Intensity Calculation ----

#' Calculate intensity within bounding boxes
#' 
#' Calculates the intensity of NMR signals within specified bounding boxes.
#' Supports two methods: simple sum (fast) or peak fitting (more accurate).
#'
#' @param mat Numeric matrix. The spectrum data matrix
#' @param ppm_x Numeric vector. F2 (x-axis) chemical shift values
#' @param ppm_y Numeric vector. F1 (y-axis) chemical shift values
#' @param boxes Data frame. Bounding boxes with columns: xmin, xmax, ymin, ymax
#' @param method Character. Integration method: "sum" (default) or "fit"
#' @param model Character. Fitting model when method="fit": "gaussian" (default) or "voigt"
#' @return Numeric vector of intensities, one per box
#' @note When method="fit", requires calculate_fitted_volumes() from Peak_fitting.R
get_box_intensity <- function(mat, ppm_x, ppm_y, boxes, method = "sum", model = "gaussian") {
  if (nrow(boxes) == 0) return(numeric(0))
  
  if (method == "sum") {
    # Fast summation method
    xmin_v <- as.numeric(boxes$xmin)
    xmax_v <- as.numeric(boxes$xmax)
    ymin_v <- as.numeric(boxes$ymin)
    ymax_v <- as.numeric(boxes$ymax)
    
    vapply(seq_along(xmin_v), FUN.VALUE = 0.0, FUN = function(i) {
      xi <- which(ppm_x >= xmin_v[i] & ppm_x <= xmax_v[i])
      yi <- which(ppm_y >= ymin_v[i] & ppm_y <= ymax_v[i])
      if (length(xi) == 0 || length(yi) == 0) return(NA_real_)
      sum(mat[yi, xi], na.rm = TRUE)
    })
    
  } else {
    # Peak fitting method (slower but more accurate)
    fit_results <- calculate_fitted_volumes(mat, ppm_x, ppm_y, boxes, model = model)
    fit_results$volume_fitted
  }
}


# Peak Range Parsing ----

#' Parse peak range exclusion string
#' 
#' Parses a semicolon-separated string of coordinate pairs defining
#' exclusion zones for peak detection (e.g., solvent regions).
#'
#' @param input_string Character. Format: "max1,min1; max2,min2; ..."
#' @return List of numeric vectors, each with 2 elements (max, min), or NULL
#' @examples
#' parse_keep_peak_ranges("0.5,-0.5; 1,0.8")
#' # Returns: list(c(0.5, -0.5), c(1, 0.8))
parse_keep_peak_ranges <- function(input_string) {
  if (is.null(input_string) || input_string == "") return(NULL)
  
  pairs <- strsplit(input_string, ";")[[1]]
  parsed <- lapply(pairs, function(p) {
    nums <- as.numeric(trimws(unlist(strsplit(p, ","))))
    if (length(nums) == 2 && all(!is.na(nums))) nums else NULL
  })
  parsed[!sapply(parsed, is.null)]
}


# Centroids Data Cleaning ----

#' Clean imported centroids dataframe
#' 
#' Cleans a dataframe of imported peak centroids by converting string
#' representations to numeric values (handles European decimal comma).
#'
#' @param df Data frame with columns: F2_ppm, F1_ppm, Volume
#' @return Data frame with cleaned numeric columns
clean_centroids_df <- function(df) {
  df$F2_ppm <- as.numeric(gsub(",", ".", trimws(df$F2_ppm)))
  df$F1_ppm <- as.numeric(gsub(",", ".", trimws(df$F1_ppm)))
  df$Volume <- as.numeric(gsub(",", ".", trimws(df$Volume)))
  df
}


# Batch Intensity Calculation ----

#' Calculate box intensities across multiple spectra
#' 
#' Calculates peak intensities for a set of reference bounding boxes
#' across multiple spectra. Handles data validation, duplicate detection,
#' and optional spectrum alignment.
#'
#' @param reference_boxes Data frame. Reference boxes with columns:
#'   xmin, xmax, ymin, ymax, and optionally stain_id
#' @param spectra_list Named list of spectrum data objects
#' @param apply_shift Logical. If TRUE, attempts to align spectra (default: FALSE)
#' @param method Character. Integration method: "sum" (default) or "fit"
#' @param model Character. Fitting model: "gaussian" (default) or "voigt"
#' @param progress Function. Progress callback with signature (value, detail)
#' @param shift_tolerance_ppm Numeric. Tolerance for dynamic box recentering (default: 0).
#'   When > 0, each box is temporarily expanded by this amount in all directions,
#'   the local maximum is found within the expanded region, and the box is 
#'   recentered on that maximum ONLY IF the new position captures more intensity
#'   than the original position. This compensates for small chemical shift 
#'   variations between spectra (e.g., due to pH, temperature) while avoiding
#'   false recentering on noise or artifacts.
#'   Typical values: 0.01-0.05 ppm. Set to 0 to disable.
#' @return Data frame with columns:
#'   - stain_id: Box identifier
#'   - F2_ppm, F1_ppm: Box center coordinates
#'   - xmin, xmax, ymin, ymax: Box boundaries
#'   - Intensity_<spectrum_name>: One column per spectrum with intensities
#' @note Requires get_box_intensity() and optionally calculate_fitted_volumes()

#' Compute the F2-only recentering shift for a box (Option B: bounded shift)
#'
#' The search window IS widened by `tol` on F2 (so a peak that has drifted
#' slightly out of the reference box can still be recovered). Within that
#' widened window the INTENSITY-WEIGHTED F2 centroid is computed, matching the
#' peak-picking convention (weighted.mean(F2_ppm, Volume)) used elsewhere in
#' APPIN, rather than a raw pixel maximum (more stable, sub-pixel, less
#' sensitive to single hot pixels). The resulting displacement of the box
#' centre is then CLAMPED to +/- tol on F2. This keeps chemical-shift
#' compensation while preventing a box from snapping onto a far-away dominant
#' feature such as the diagonal / water ridge in homonuclear experiments
#' (COSY / TOCSY / UF-COSY), since such a feature lies further than `tol` away.
#' F1 (y) is never shifted, consistent with the F2-only UI.
#'
#' Both the search window and the final box centre are additionally confined to
#' the exclusive F2 territory [terr_lo, terr_hi] (midpoints to F1-overlapping
#' neighbours). This stops neighbouring boxes from competing for the same peak
#' or overlapping after recentering, so `tol` can be larger without conflicts.
#'
#' @param terr_lo,terr_hi Exclusive F2 territory bounds (ppm). Default +/-Inf
#'   (no neighbour constraint).
#' @return Numeric F2 shift (ppm) to add to xmin/xmax. 0 if no valid signal.
recenter_box_f2_shift <- function(mat, ppm_x, ppm_y,
                                  xmin, xmax, ymin, ymax,
                                  tol,
                                  terr_lo = -Inf, terr_hi = Inf) {
  if (tol <= 0) return(0)
  
  # Search window widened by tol on F2, but never beyond the box's territory
  lo <- max(xmin - tol, terr_lo)
  hi <- min(xmax + tol, terr_hi)
  if (lo >= hi) return(0)
  
  x_idx <- which(ppm_x >= lo & ppm_x <= hi)
  y_idx <- which(ppm_y >= ymin & ppm_y <= ymax)
  if (length(x_idx) == 0 || length(y_idx) == 0) return(0)
  
  submat <- mat[y_idx, x_idx, drop = FALSE]
  
  # Keep only positive signal as weights (negative HSQC lobes -> 0)
  w <- pmax(submat, 0)
  w[!is.finite(w)] <- 0
  
  # Collapse F1 rows: total positive intensity per F2 column
  col_weights <- colSums(w)
  if (sum(col_weights) <= 0) return(0)
  
  # Intensity-weighted F2 centroid (same convention as the peak picking)
  f2_centroid <- sum(ppm_x[x_idx] * col_weights) / sum(col_weights)
  
  box_center_f2 <- (xmin + xmax) / 2
  shift_f2 <- f2_centroid - box_center_f2
  
  # CLAMP the displacement to +/- tol on F2
  if (abs(shift_f2) > tol) shift_f2 <- sign(shift_f2) * tol
  
  # CONFINE the shifted box centre to its exclusive territory
  new_center <- box_center_f2 + shift_f2
  if (new_center < terr_lo) shift_f2 <- terr_lo - box_center_f2
  if (new_center > terr_hi) shift_f2 <- terr_hi - box_center_f2
  
  shift_f2
}

#' Compute exclusive F2 territory bounds for each box
#'
#' Two boxes compete for the same peak only if they also overlap in F1 (same
#' "row" of the spectrum). For each box, the F2 territory is bounded by the
#' midpoint to the nearest F1-overlapping neighbour on each side. Recentering
#' (both the search window and the final shift) is then confined to this
#' territory, so neighbouring integration boxes can never drift onto the same
#' peak or overlap after recentering — which lets `tol` be larger without
#' creating conflicts. Especially important along the TOCSY/COSY diagonal where
#' boxes are densely packed.
#'
#' @return data.frame with columns terr_lo, terr_hi (F2 ppm bounds) per row of boxes.
compute_f2_territories <- function(boxes) {
  n <- nrow(boxes)
  terr_lo <- rep(-Inf, n)
  terr_hi <- rep(Inf, n)
  if (n < 2) return(data.frame(terr_lo = terr_lo, terr_hi = terr_hi))
  
  cx <- (boxes$xmin + boxes$xmax) / 2  # F2 centres
  
  for (i in seq_len(n)) {
    # F1-overlapping neighbours (share a row); exclude self
    f1_overlap <- boxes$ymin < boxes$ymax[i] & boxes$ymax > boxes$ymin[i]
    f1_overlap[i] <- FALSE
    if (!any(f1_overlap)) next
    
    # Nearest neighbour centre on the low-F2 side and high-F2 side
    lo_side <- which(f1_overlap & cx < cx[i])
    hi_side <- which(f1_overlap & cx > cx[i])
    
    if (length(lo_side) > 0) {
      nearest <- lo_side[which.max(cx[lo_side])]
      terr_lo[i] <- (cx[i] + cx[nearest]) / 2  # midpoint
    }
    if (length(hi_side) > 0) {
      nearest <- hi_side[which.min(cx[hi_side])]
      terr_hi[i] <- (cx[i] + cx[nearest]) / 2  # midpoint
    }
  }
  data.frame(terr_lo = terr_lo, terr_hi = terr_hi)
}

calculate_batch_box_intensities <- function(reference_boxes, 
                                            spectra_list, 
                                            apply_shift = FALSE, 
                                            method = "sum",
                                            model = "gaussian",
                                            progress = NULL,
                                            shift_tolerance_ppm = 0) {
  
  # ========== VALIDATIONS ==========
  if (is.null(reference_boxes) || nrow(reference_boxes) == 0) {
    stop("reference_boxes is empty or NULL")
  }
  
  # Copy to avoid modifying the original
  ref_boxes <- as.data.frame(reference_boxes)
  
  # Check required columns
  required_cols <- c("xmin", "xmax", "ymin", "ymax")
  missing_cols <- setdiff(required_cols, names(ref_boxes))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # ========== CLEANING BOXES ==========
  
  # Add stain_id if missing
  if (!"stain_id" %in% names(ref_boxes)) {
    ref_boxes$stain_id <- paste0("box_", seq_len(nrow(ref_boxes)))
  }
  
  # Remove rows with NA coordinates
  ref_boxes <- ref_boxes[
    !is.na(ref_boxes$xmin) & !is.na(ref_boxes$xmax) & 
      !is.na(ref_boxes$ymin) & !is.na(ref_boxes$ymax), , drop = FALSE
  ]
  
  if (nrow(ref_boxes) == 0) {
    stop("All boxes have NA coordinates")
  }
  
  # Fix inverted boxes (xmin > xmax or ymin > ymax)
  for (i in seq_len(nrow(ref_boxes))) {
    if (ref_boxes$xmin[i] > ref_boxes$xmax[i]) {
      tmp <- ref_boxes$xmin[i]
      ref_boxes$xmin[i] <- ref_boxes$xmax[i]
      ref_boxes$xmax[i] <- tmp
    }
    if (ref_boxes$ymin[i] > ref_boxes$ymax[i]) {
      tmp <- ref_boxes$ymin[i]
      ref_boxes$ymin[i] <- ref_boxes$ymax[i]
      ref_boxes$ymax[i] <- tmp
    }
  }
  
  # ========== DUPLICATE MANAGEMENT ==========
  
  # Check stain_id duplicates
  if (any(duplicated(ref_boxes$stain_id))) {
    warning("Duplicate stain_id detected - automatic renaming")
    dup_ids <- ref_boxes$stain_id[duplicated(ref_boxes$stain_id)]
    for (dup_id in unique(dup_ids)) {
      idx <- which(ref_boxes$stain_id == dup_id)
      if (length(idx) > 1) {
        ref_boxes$stain_id[idx[-1]] <- paste0(dup_id, "_dup", seq_along(idx[-1]))
      }
    }
  }
  
  # Check boxes with identical coordinates
  coord_signature <- paste(ref_boxes$xmin, ref_boxes$xmax, 
                           ref_boxes$ymin, ref_boxes$ymax, sep = "_")
  if (any(duplicated(coord_signature))) {
    warning("Boxes with identical coordinates detected - removing duplicates")
    ref_boxes <- ref_boxes[!duplicated(coord_signature), , drop = FALSE]
  }
  
  # ========== CALCULATION OF CENTERS ==========
  ref_boxes$F2_ppm <- (ref_boxes$xmin + ref_boxes$xmax) / 2
  ref_boxes$F1_ppm <- (ref_boxes$ymin + ref_boxes$ymax) / 2
  
  # ========== BUILD RESULT DATAFRAME ==========
  n_boxes <- nrow(ref_boxes)
  message(sprintf("Processing %d boxes across %d spectra", n_boxes, length(spectra_list)))
  
  result_df <- data.frame(
    stain_id = ref_boxes$stain_id,
    F2_ppm = ref_boxes$F2_ppm,
    F1_ppm = ref_boxes$F1_ppm,
    xmin = ref_boxes$xmin,
    xmax = ref_boxes$xmax,
    ymin = ref_boxes$ymin,
    ymax = ref_boxes$ymax,
    stringsAsFactors = FALSE
  )
  
  # ========== EXCLUSIVE F2 TERRITORIES ==========
  # Computed once on the reference boxes; confines recentering so neighbouring
  # boxes (overlapping in F1) cannot compete for the same peak or overlap.
  territories <- compute_f2_territories(ref_boxes)
  
  # ========== CALCULATE INTENSITIES PER SPECTRUM ==========
  for (spectrum_name in names(spectra_list)) {
    spectrum_data <- spectra_list[[spectrum_name]]
    
    if (is.null(spectrum_data) || is.null(spectrum_data$spectrumData)) {
      warning(paste("Spectrum", spectrum_name, "invalid - column filled with NA"))
      col_name <- paste0("Intensity_", make.names(basename(spectrum_name)))
      result_df[[col_name]] <- rep(NA_real_, n_boxes)
      next
    }
    
    mat <- spectrum_data$spectrumData
    ppm_x <- suppressWarnings(as.numeric(colnames(mat)))
    ppm_y <- suppressWarnings(as.numeric(rownames(mat)))
    
    # Check that ppm values are valid
    if (any(is.na(ppm_x)) || any(is.na(ppm_y))) {
      warning(paste("Spectrum", spectrum_name, "has invalid ppm values"))
      col_name <- paste0("Intensity_", make.names(basename(spectrum_name)))
      result_df[[col_name]] <- rep(NA_real_, n_boxes)
      next
    }
    
    # Calculate shift if requested
    shift_f2 <- 0
    shift_f1 <- 0
    
    if (apply_shift && n_boxes > 0) {
      max_idx <- which(mat == max(mat, na.rm = TRUE), arr.ind = TRUE)
      if (length(max_idx) > 0 && nrow(max_idx) > 0) {
        max_f2 <- ppm_x[max_idx[1, 2]]
        max_f1 <- ppm_y[max_idx[1, 1]]
        shift_f2 <- max_f2 - ref_boxes$F2_ppm[1]
        shift_f1 <- max_f1 - ref_boxes$F1_ppm[1]
        # Limit shift to 0.5 ppm max
        if (abs(shift_f2) > 0.5) shift_f2 <- 0
        if (abs(shift_f1) > 0.5) shift_f1 <- 0
      }
    }
    
    if (method == "sum") {
      # Calculate intensities for each box
      intensities <- numeric(n_boxes)
      
      for (i in seq_len(n_boxes)) {
        # Base box coordinates (with global shift if apply_shift is TRUE)
        xmin_base <- ref_boxes$xmin[i] + shift_f2
        xmax_base <- ref_boxes$xmax[i] + shift_f2
        ymin_base <- ref_boxes$ymin[i] + shift_f1
        ymax_base <- ref_boxes$ymax[i] + shift_f1
        
        # Sum over a box (helper; NA if box falls outside the grid)
        box_sum <- function(x0, x1) {
          xi <- which(ppm_x >= x0 & ppm_x <= x1)
          yi <- which(ppm_y >= ymin_base & ppm_y <= ymax_base)
          if (length(xi) == 0 || length(yi) == 0) return(NA_real_)
          sum(mat[yi, xi, drop = FALSE], na.rm = TRUE)
        }
        
        intensities[i] <- box_sum(xmin_base, xmax_base)
        
        # Dynamic recentering: F2-only shift (clamped to +/- tol, confined to
        # the box territory). Strict guard: keep the recentering only if it does
        # not reduce the integrated volume vs the original position.
        if (shift_tolerance_ppm > 0) {
          local_shift_f2 <- recenter_box_f2_shift(
            mat, ppm_x, ppm_y,
            xmin_base, xmax_base, ymin_base, ymax_base,
            tol = shift_tolerance_ppm,
            terr_lo = territories$terr_lo[i] + shift_f2,
            terr_hi = territories$terr_hi[i] + shift_f2
          )
          if (local_shift_f2 != 0) {
            v_shft <- box_sum(xmin_base + local_shift_f2, xmax_base + local_shift_f2)
            v_orig <- intensities[i]
            # Accept recentering only if it does not degrade the volume
            if (is.finite(v_shft) && (!is.finite(v_orig) || v_shft >= v_orig)) {
              intensities[i] <- v_shft
            }
          }
        }
      }
      
    } else {
      # Fitting method - also apply dynamic recentering if tolerance > 0
      base_cols <- c("xmin", "xmax", "ymin", "ymax", "stain_id")
      
      # Original boxes (with global shift only)
      orig_boxes <- ref_boxes[, base_cols]
      orig_boxes$xmin <- orig_boxes$xmin + shift_f2
      orig_boxes$xmax <- orig_boxes$xmax + shift_f2
      orig_boxes$ymin <- orig_boxes$ymin + shift_f1
      orig_boxes$ymax <- orig_boxes$ymax + shift_f1
      
      if (shift_tolerance_ppm > 0) {
        # Candidate recentered boxes (F2-only, clamped + territory-confined)
        shifted_boxes <- orig_boxes
        for (i in seq_len(nrow(shifted_boxes))) {
          local_shift_f2 <- recenter_box_f2_shift(
            mat, ppm_x, ppm_y,
            orig_boxes$xmin[i], orig_boxes$xmax[i],
            orig_boxes$ymin[i], orig_boxes$ymax[i],
            tol = shift_tolerance_ppm,
            terr_lo = territories$terr_lo[i] + shift_f2,
            terr_hi = territories$terr_hi[i] + shift_f2
          )
          shifted_boxes$xmin[i] <- orig_boxes$xmin[i] + local_shift_f2
          shifted_boxes$xmax[i] <- orig_boxes$xmax[i] + local_shift_f2
        }
        
        # Fit BOTH positions with the same model, then keep the recentering
        # per box ONLY if it does not reduce the fitted volume (strict guard).
        # Bounded by tol + territory, so this can never re-snap onto the diagonal.
        fit_orig    <- calculate_fitted_volumes(mat, ppm_x, ppm_y, orig_boxes,    model = model)
        fit_shifted <- calculate_fitted_volumes(mat, ppm_x, ppm_y, shifted_boxes, model = model)
        
        v_orig <- fit_orig$volume_fitted
        v_shft <- fit_shifted$volume_fitted
        v_orig[!is.finite(v_orig)] <- -Inf  # NA orig -> always prefer shifted
        v_shft[!is.finite(v_shft)] <- -Inf  # NA shifted -> keep orig
        intensities <- ifelse(v_shft >= v_orig, fit_shifted$volume_fitted, fit_orig$volume_fitted)
      } else {
        # No tolerance - original boxes only
        fit_results <- calculate_fitted_volumes(mat, ppm_x, ppm_y, orig_boxes, model = model)
        intensities <- fit_results$volume_fitted
      }
    }
    
    col_name <- paste0("Intensity_", make.names(basename(spectrum_name)))
    result_df[[col_name]] <- intensities
  }
  
  # ========== FINAL VERIFICATION ==========
  if (nrow(result_df) != n_boxes) {
    warning(sprintf("Anomaly: %d rows expected, %d obtained", n_boxes, nrow(result_df)))
  }
  
  message(sprintf("Export: %d boxes, %d intensity columns", 
                  nrow(result_df), 
                  sum(grepl("^Intensity_", names(result_df)))))
  
  return(result_df)
}