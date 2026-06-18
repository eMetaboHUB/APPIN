# Function/Peak_fitting.R

# DETECT LOCAL MAX ----

#' Detect local maxima in a 2D region
#' @param mat Matrix of intensities
#' @param threshold Minimum intensity relative to max (0-1)
#' @param min_distance Minimum distance between peaks (in pixels)
#' @return Data frame with local maxima coordinates


detect_local_maxima <- function(mat, threshold = 0.3, min_distance = 2) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # ========== FIX: Handle small regions ==========
  # For very small regions, just return the global maximum
  if (nr < 3 || nc < 3) {
    max_val <- max(mat, na.rm = TRUE)
    if (is.na(max_val) || is.infinite(max_val)) {
      return(data.frame(row = integer(), col = integer(), value = numeric()))
    }
    max_pos <- which(mat == max_val, arr.ind = TRUE)
    if (nrow(max_pos) > 0) {
      return(data.frame(row = max_pos[1, 1], col = max_pos[1, 2], value = max_val))
    }
    return(data.frame(row = integer(), col = integer(), value = numeric()))
  }
  # ========== END FIX ==========
  
  max_val <- max(mat, na.rm = TRUE)
  min_val <- min(mat, na.rm = TRUE)
  thresh_abs <- min_val + threshold * (max_val - min_val)
  
  maxima <- data.frame(row = integer(), col = integer(), value = numeric())
  
  # ========== FIX: Include border pixels (was 2:(nr-1), now 1:nr) ==========
  for (i in 1:nr) {
    for (j in 1:nc) {
      val <- mat[i, j]
      if (is.na(val) || val < thresh_abs) next
      
      # Check if local maximum (8-connectivity) with boundary handling
      row_range <- max(1, i-1):min(nr, i+1)
      col_range <- max(1, j-1):min(nc, j+1)
      neighborhood <- mat[row_range, col_range, drop = FALSE]
      
      if (val >= max(neighborhood, na.rm = TRUE)) {
        maxima <- rbind(maxima, data.frame(row = i, col = j, value = val))
      }
    }
  }
  # ========== END FIX ==========
  
  if (nrow(maxima) == 0) return(maxima)
  
  # Remove peaks too close to each other (keep the highest)
  maxima <- maxima[order(-maxima$value), ]
  keep <- rep(TRUE, nrow(maxima))
  
  for (i in seq_len(nrow(maxima))) {
    if (!keep[i]) next
    for (j in seq_len(nrow(maxima))) {
      if (i >= j || !keep[j]) next
      dist <- sqrt((maxima$row[i] - maxima$row[j])^2 + (maxima$col[i] - maxima$col[j])^2)
      if (dist < min_distance) {
        keep[j] <- FALSE
      }
    }
  }
  
  maxima[keep, ]
}


# Pseudo-Voigt 2D profile function ----

#' Approximation of the Voigt profile as a linear combination of Gaussian and Lorentzian
#' @param x X coordinate
#' @param y Y coordinate
#' @param A Amplitude
#' @param x0 Center X
#' @param y0 Center Y
#' @param sigma_x Gaussian width in X
#' @param sigma_y Gaussian width in Y
#' @param gamma_x Lorentzian width in X
#' @param gamma_y Lorentzian width in Y
#' @param eta Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)
#' @return Intensity value

pseudo_voigt_2d <- function(x, y, A, x0, y0, sigma_x, sigma_y, gamma_x, gamma_y, eta) {
  # Gaussian component
  gauss <- exp(-((x - x0)^2 / (2 * sigma_x^2) + (y - y0)^2 / (2 * sigma_y^2)))
  
  # Lorentzian component
  lorentz <- 1 / (1 + ((x - x0) / gamma_x)^2 + ((y - y0) / gamma_y)^2)
  
  # Pseudo-Voigt: linear combination
  A * (eta * lorentz + (1 - eta) * gauss)
}


# Fit 2D peak to a spectral region (VERSION ROBUSTE avec support MULTIPLET) ----
#' 
#' @param mat Matrix of spectral intensities
#' @param ppm_x Vector of F2 chemical shifts
#' @param ppm_y Vector of F1 chemical shifts
#' @param box Data frame with xmin, xmax, ymin, ymax
#' @param model Type of peak model: "gaussian", "voigt"
#' @param min_points Minimum number of points required for fitting
#' @return List with fitted parameters and volume

fit_2d_peak <- function(mat, ppm_x, ppm_y, box, model = "gaussian", min_points = 25) {
  
  # Extract region
  x_idx <- which(ppm_x >= box$xmin & ppm_x <= box$xmax)
  y_idx <- which(ppm_y >= box$ymin & ppm_y <= box$ymax)
  
  if (length(x_idx) == 0 || length(y_idx) == 0) {
    return(list(volume = NA, fit_quality = NA, params = NULL, 
                method = "failed", error = "No points in region"))
  }
  
  # Check that we have enough points
  
  if (length(x_idx) * length(y_idx) < min_points) {
    # Fallbac on sum
    region <- mat[y_idx, x_idx, drop = FALSE]
    volume_sum <- sum(region, na.rm = TRUE)
    return(list(
      volume = volume_sum,
      fit_quality = NA,
      params = NULL,
      method = "sum_fit_failed",
      error = "Too few points for fitting"
    ))
  }
  
  region <- mat[y_idx, x_idx, drop = FALSE]
  x_sub <- ppm_x[x_idx]
  y_sub <- ppm_y[y_idx]
  
  # Check that the region is not empty or constant
  if (all(is.na(region)) || sd(as.vector(region), na.rm = TRUE) < 1e-10) {
    return(list(
      volume = sum(region, na.rm = TRUE),
      fit_quality = NA,
      params = NULL,
      method = "sum_fit_failed",
      error = "Region is constant or empty"
    ))
  }
  
  #  MULTIPLET DETECTION 
  # Detect if there are multiple peaks in the region
  # NOTE: threshold=0.5 (50% of max) and min_distance=3 to avoid detecting noise as peaks
  local_max <- detect_local_maxima(region, threshold = 0.5, min_distance = 3)
  n_peaks <- nrow(local_max)
  
  # ========== FIX: Handle case where no local maximum is detected ==========
  # This can happen when:
  
  # 1. Region is too small (< 3x3 pixels)
  # 2. Maximum is on the border (excluded by detect_local_maxima)
  # 3. All values are below the 30% threshold
  
  if (n_peaks == 0) {
    # Fallback: use the global maximum of the region as the peak center
    max_val <- max(region, na.rm = TRUE)
    max_pos <- which(region == max_val, arr.ind = TRUE)
    if (nrow(max_pos) > 0) {
      # Create a synthetic local_max entry
      local_max <- data.frame(
        row = max_pos[1, 1],
        col = max_pos[1, 2],
        value = max_val
      )
      n_peaks <- 1
    } else {
      # Truly empty region - fallback to sum
      return(list(
        volume = sum(region, na.rm = TRUE),
        fit_quality = NA,
        params = NULL,
        method = "sum_fit_failed",
        n_peaks = 0,
        is_multiplet = FALSE,
        error = "No peak detected in region"
      ))
    }
  }
  # ========== END FIX ==========
  
  is_multiplet <- n_peaks > 1
  
  # If multiplet detected, fit each peak separately
  if (is_multiplet) {
    
    # Results for each peak
    peak_volumes <- numeric(n_peaks)
    peak_r_squared <- numeric(n_peaks)
    peak_centers_x <- numeric(n_peaks)
    peak_centers_y <- numeric(n_peaks)
    peak_fitted_values <- vector("list", n_peaks)
    
    total_fitted_vals <- rep(0, length(x_sub) * length(y_sub))
    
    for (p in seq_len(n_peaks)) {
      # Peak coordinates
      peak_row <- local_max$row[p]
      peak_col <- local_max$col[p]
      peak_x <- x_sub[peak_col]
      peak_y <- y_sub[peak_row]
      peak_amplitude <- local_max$value[p]
      
      # Create a subregion around this peak
      # Estimate the width of the peak (distance to the nearest neighboring peak / 2)
      
      if (n_peaks > 1) {
        distances <- sqrt((local_max$col - peak_col)^2 + (local_max$row - peak_row)^2)
        distances[p] <- Inf  # Exclure le pic lui-même
        min_dist <- min(distances)
        half_width_pixels <- max(2, floor(min_dist / 2))
      } else {
        half_width_pixels <- 3
      }
      
      # Indices for the sub-region
      col_start <- max(1, peak_col - half_width_pixels)
      col_end <- min(ncol(region), peak_col + half_width_pixels)
      row_start <- max(1, peak_row - half_width_pixels)
      row_end <- min(nrow(region), peak_row + half_width_pixels)
      
      sub_region <- region[row_start:row_end, col_start:col_end, drop = FALSE]
      sub_x <- x_sub[col_start:col_end]
      sub_y <- y_sub[row_start:row_end]
      
      # Create the grid for this peak
      sub_grid <- expand.grid(x = sub_x, y = sub_y)
      sub_grid$z <- as.vector(t(sub_region))
      sub_grid <- sub_grid[!is.na(sub_grid$z), ]
      
      if (nrow(sub_grid) < 9) {
        # Not enough points, use the sum for this peak
        peak_volumes[p] <- sum(sub_region, na.rm = TRUE)
        peak_r_squared[p] <- NA
        peak_centers_x[p] <- peak_x
        peak_centers_y[p] <- peak_y
        next
      }
      
      # Normalization
      z_scale <- max(abs(sub_grid$z), na.rm = TRUE)
      if (z_scale < 1e-10) z_scale <- 1
      sub_grid$z_norm <- sub_grid$z / z_scale
      
      # Initial parameters for this peak
      baseline_init <- quantile(sub_grid$z, 0.1, na.rm = TRUE) / z_scale
      sigma_x_init <- diff(range(sub_x)) / 4
      sigma_y_init <- diff(range(sub_y)) / 4
      
      # Trying to fit
      fit_single <- tryCatch({
        if (model == "gaussian") {
          fit <- minpack.lm::nlsLM(
            z_norm ~ A * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2))) + b,
            data = sub_grid,
            start = list(A = peak_amplitude / z_scale, x0 = peak_x, y0 = peak_y,
                         sx = sigma_x_init, sy = sigma_y_init, b = baseline_init),
            lower = c(A = 0, x0 = min(sub_x), y0 = min(sub_y),
                      sx = diff(range(sub_x)) / 20, sy = diff(range(sub_y)) / 20, b = -Inf),
            upper = c(A = Inf, x0 = max(sub_x), y0 = max(sub_y),
                      sx = diff(range(sub_x)) * 2, sy = diff(range(sub_y)) * 2, b = Inf),
            control = list(maxiter = 100, gtol = 0)
          )
        } else {
          # Voigt (pseudo-Voigt)
          fit <- minpack.lm::nlsLM(
            z_norm ~ A * (eta / (1 + ((x - x0) / gx)^2 + ((y - y0) / gy)^2) + 
                            (1 - eta) * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2)))) + b,
            data = sub_grid,
            start = list(A = peak_amplitude / z_scale, x0 = peak_x, y0 = peak_y,
                         sx = sigma_x_init, sy = sigma_y_init,
                         gx = sigma_x_init, gy = sigma_y_init,
                         eta = 0.5, b = baseline_init),
            lower = c(A = 0, x0 = min(sub_x), y0 = min(sub_y),
                      sx = diff(range(sub_x)) / 20, sy = diff(range(sub_y)) / 20,
                      gx = diff(range(sub_x)) / 20, gy = diff(range(sub_y)) / 20,
                      eta = 0, b = -Inf),
            upper = c(A = Inf, x0 = max(sub_x), y0 = max(sub_y),
                      sx = diff(range(sub_x)) * 2, sy = diff(range(sub_y)) * 2,
                      gx = diff(range(sub_x)) * 2, gy = diff(range(sub_y)) * 2,
                      eta = 1, b = Inf),
            control = list(maxiter = 100, gtol = 0)
          )
        }
        
        params <- coef(fit)
        fitted_vals <- fitted(fit) * z_scale
        
        # Calculate R²
        residuals <- sub_grid$z - fitted_vals
        ss_res <- sum(residuals^2)
        ss_tot <- sum((sub_grid$z - mean(sub_grid$z))^2)
        r2 <- max(0, 1 - (ss_res / ss_tot))
        
        # Volume = sum of fitted values
        vol <- sum(fitted_vals, na.rm = TRUE)
        
        list(volume = vol, r_squared = r2, 
             center_x = params["x0"], center_y = params["y0"],
             fitted_vals = fitted_vals, success = TRUE)
        
      }, error = function(e) {
        # Fallback: sum for this peak
        list(volume = sum(sub_region, na.rm = TRUE), r_squared = NA,
             center_x = peak_x, center_y = peak_y,
             fitted_vals = NULL, success = FALSE)
      })
      
      peak_volumes[p] <- fit_single$volume
      peak_r_squared[p] <- fit_single$r_squared
      peak_centers_x[p] <- fit_single$center_x
      peak_centers_y[p] <- fit_single$center_y
      peak_fitted_values[[p]] <- fit_single$fitted_vals
    }
    
    # Aggregate the results
    total_volume <- sum(peak_volumes, na.rm = TRUE)
    mean_r_squared <- mean(peak_r_squared, na.rm = TRUE)
    
    # Volume-weighted center
    if (sum(peak_volumes, na.rm = TRUE) > 0) {
      weights <- peak_volumes / sum(peak_volumes, na.rm = TRUE)
      center_x <- sum(peak_centers_x * weights, na.rm = TRUE)
      center_y <- sum(peak_centers_y * weights, na.rm = TRUE)
    } else {
      center_x <- mean(peak_centers_x, na.rm = TRUE)
      center_y <- mean(peak_centers_y, na.rm = TRUE)
    }
    
    # For multiplets, we don't reconstruct fitted_values here
    # The visualization will handle it by fitting a single model to the whole region
    
    return(list(
      volume = total_volume,
      fit_quality = mean_r_squared,
      params = c(x0 = center_x, y0 = center_y),
      fitted_values = NULL,  # Will be computed in visualization if needed
      residuals = NULL,
      method = "multiplet_fit",
      n_peaks = n_peaks,
      is_multiplet = TRUE,
      peak_volumes = peak_volumes,
      peak_r_squared = peak_r_squared,
      peak_centers_x = peak_centers_x,
      peak_centers_y = peak_centers_y,
      error = NULL
    ))
  }
  
  # Create grid
  grid <- expand.grid(x = x_sub, y = y_sub)
  grid$z <- as.vector(t(region))
  grid <- grid[!is.na(grid$z), ]  # Enlever les NA
  
  if (nrow(grid) < min_points) {
    return(list(
      volume = sum(region, na.rm = TRUE),
      fit_quality = NA,
      params = NULL,
      method = "sum_fit_failed",
      error = "Too many NA values"
    ))
  }
  
  #  ROBUST ESTIMATE OF INITIAL PARAMETERS 
  
  # Find the overall maximum for centering
  max_idx <- which.max(grid$z)
  center_x_init <- grid$x[max_idx]
  center_y_init <- grid$y[max_idx]
  amplitude_init <- max(grid$z, na.rm = TRUE)
  
  baseline_init <- quantile(grid$z, 0.1, na.rm = TRUE)  # 10e percentile comme baseline
  
  # Width estimates based on approximate FWHM
  # Find points > 50% of the max
  
  half_max <- (amplitude_init + baseline_init) / 2
  points_half <- grid[grid$z > half_max, ]
  
  if (nrow(points_half) > 3) {
    sigma_x_init <- sd(points_half$x) * 1.5  # Empirical factor
    sigma_y_init <- sd(points_half$y) * 1.5
  } else {
    # Fallback on box size
    sigma_x_init <- (box$xmax - box$xmin) / 4
    sigma_y_init <- (box$ymax - box$ymin) / 4
  }
  
  # Ensure that the sigma values are not too small
  sigma_x_init <- max(sigma_x_init, diff(range(x_sub)) / 10)
  sigma_y_init <- max(sigma_y_init, diff(range(y_sub)) / 10)
  
  # Standardize data to improve convergence
  z_scale <- max(abs(grid$z), na.rm = TRUE)
  if (z_scale < 1e-10) z_scale <- 1
  grid$z_norm <- grid$z / z_scale
  
  # MODEL AND FITTING 
  
  if (model == "gaussian") {
    
    fit_formula <- z_norm ~ A * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2))) + b
    
    start_params <- list(
      A = amplitude_init / z_scale,
      x0 = center_x_init,
      y0 = center_y_init,
      sx = sigma_x_init,
      sy = sigma_y_init,
      b = baseline_init / z_scale
    )
    
    # Constraints to avoid aberrant parameters
    lower_bounds <- c(
      A = 0,  # Positive amplitude
      x0 = min(x_sub),
      y0 = min(y_sub),
      sx = diff(range(x_sub)) / 20,  # Min width
      sy = diff(range(y_sub)) / 20,
      b = -Inf
    )
    
    upper_bounds <- c(
      A = Inf,
      x0 = max(x_sub),
      y0 = max(y_sub),
      sx = diff(range(x_sub)) * 2,  # Max width
      sy = diff(range(y_sub)) * 2,
      b = Inf
    )
    
  } else if (model == "voigt") {
    
    # Pseudo-Voigt 2D: linear combination of Gaussian and Lorentzian
    # eta = mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)
    
    fit_formula <- z_norm ~ A * (eta / (1 + ((x - x0) / gx)^2 + ((y - y0) / gy)^2) + 
                                   (1 - eta) * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2)))) + b
    
    start_params <- list(
      A = amplitude_init / z_scale,
      x0 = center_x_init,
      y0 = center_y_init,
      sx = sigma_x_init,      # Gaussian Width
      sy = sigma_y_init,
      gx = sigma_x_init,      # Lorentzian Width
      gy = sigma_y_init,
      eta = 0.5,              # 50-50 mix by default
      b = baseline_init / z_scale
    )
    
    lower_bounds <- c(
      A = 0,
      x0 = min(x_sub),
      y0 = min(y_sub),
      sx = diff(range(x_sub)) / 20,
      sy = diff(range(y_sub)) / 20,
      gx = diff(range(x_sub)) / 20,
      gy = diff(range(y_sub)) / 20,
      eta = 0,                # Pure Gaussian
      b = -Inf
    )
    
    upper_bounds <- c(
      A = Inf,
      x0 = max(x_sub),
      y0 = max(y_sub),
      sx = diff(range(x_sub)) * 2,
      sy = diff(range(y_sub)) * 2,
      gx = diff(range(x_sub)) * 2,
      gy = diff(range(y_sub)) * 2,
      eta = 1,                # Pure Lorentzian
      b = Inf
    )
    
  } else {
    stop("Model not supported. Use 'gaussian' or 'voigt'")
  }
  
  # ATTEMPTED FITTING 
  
  fit_result <- tryCatch({
    
    fit <- minpack.lm::nlsLM(
      fit_formula,
      data = grid,
      start = start_params,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(
        maxiter = 200,
        ftol = 1e-6,
        ptol = 1e-6,
        gtol = 0  # Disable gradient testing to avoid "singular gradient"
      )
    )
    
    params <- coef(fit)
    
    # Rescale settings
    params["A"] <- params["A"] * z_scale
    params["b"] <- params["b"] * z_scale
    
    # Calculate R² on the original data
    fitted_vals <- fitted(fit) * z_scale
    residuals <- grid$z - fitted_vals
    ss_res <- sum(residuals^2)
    ss_tot <- sum((grid$z - mean(grid$z))^2)
    r_squared <- max(0, 1 - (ss_res / ss_tot)) # Force between 0 and 1
    
    
    #  VOLUME CALCULATION 
    volume_fitted_sum <- sum(fitted_vals, na.rm = TRUE)
    
    # Analytical volume (for reference)
    if (model == "gaussian") {
      volume_analytical <- 2 * pi * params["A"] * abs(params["sx"]) * abs(params["sy"])
    } else if (model == "voigt") {
      eta_val <- if ("eta" %in% names(params)) params["eta"] else 0.5
      vol_gauss <- 2 * pi * params["A"] * abs(params["sx"]) * abs(params["sy"])
      vol_lorentz <- pi^2 * params["A"] * abs(params["gx"]) * abs(params["gy"])
      volume_analytical <- (1 - eta_val) * vol_gauss + eta_val * vol_lorentz
    } else {
      volume_analytical <- NA
    }
    
    list(
      volume = volume_fitted_sum,
      volume_analytical = volume_analytical,
      fit_quality = r_squared,
      params = params,
      fitted_values = fitted_vals,
      residuals = residuals,
      method = model,
      error = NULL
    )
    
  }, error = function(e) {
    
    # Always fallback to SUM on any error
    list(
      volume = sum(region, na.rm = TRUE),
      fit_quality = NA,
      params = NULL,
      fitted_values = NULL,
      residuals = NULL,
      method = "sum_fit_failed",
      error = e$message
    )
  })
  
  
  
  # Add information on the number of peaks (for non-multiplets, n_peaks = 1)
  
  fit_result$n_peaks <- 1
  fit_result$is_multiplet <- FALSE
  
  return(fit_result)
}

# FIT DIAGNOSTIC VISUALIZATION ----

#' Generate fit diagnostic plot for a single peak (TopSpin-style)
#' 
#' Creates a 1D slice visualization showing experimental data vs fitted curve,
#' allowing visual assessment of fit quality (similar to TopSpin display).
#' 
#' @param fit_result Output from fit_2d_peak containing fitted_values and fit_quality
#' @param mat Original spectral matrix
#' @param ppm_x F2 chemical shifts vector
#' @param ppm_y F1 chemical shifts vector
#' @param box Data frame with xmin, xmax, ymin, ymax
#' @param slice_direction Character: "F2" (horizontal) or "F1" (vertical) slice
#' @return A list with plot data for rendering (experimental, fitted, residuals)
#' @export
generate_fit_diagnostic_data <- function(fit_result, mat, ppm_x, ppm_y, box, 
                                         slice_direction = "F2") {
  
  # Check if fit data is available
  if (is.null(fit_result$fitted_values) || is.null(fit_result$fit_quality)) {
    return(list(
      success = FALSE,
      error = "No fitted values available (fit may have failed)"
    ))
  }
  
  # Extract region
  x_idx <- which(ppm_x >= box$xmin & ppm_x <= box$xmax)
  y_idx <- which(ppm_y >= box$ymin & ppm_y <= box$ymax)
  
  if (length(x_idx) == 0 || length(y_idx) == 0) {
    return(list(success = FALSE, error = "Invalid box coordinates"))
  }
  
  region <- mat[y_idx, x_idx, drop = FALSE]
  x_sub <- ppm_x[x_idx]
  y_sub <- ppm_y[y_idx]
  
  # Reshape fitted values to matrix form
  n_y <- length(y_idx)
  n_x <- length(x_idx)
  fitted_matrix <- matrix(fit_result$fitted_values, nrow = n_y, ncol = n_x, byrow = TRUE)
  
  # Find maximum position for slice
  max_pos <- which(region == max(region, na.rm = TRUE), arr.ind = TRUE)
  if (nrow(max_pos) > 1) max_pos <- max_pos[1, , drop = FALSE]
  
  if (slice_direction == "F2") {
    # Horizontal slice (along F2) at the row of maximum
    slice_row <- max_pos[1, 1]
    exp_slice <- region[slice_row, ]
    fit_slice <- fitted_matrix[slice_row, ]
    ppm_axis <- x_sub
    axis_label <- "F2 (ppm)"
    slice_pos <- y_sub[slice_row]
    slice_info <- paste0("F1 = ", round(slice_pos, 3), " ppm")
  } else {
    # Vertical slice (along F1) at the column of maximum
    slice_col <- max_pos[1, 2]
    exp_slice <- region[, slice_col]
    fit_slice <- fitted_matrix[, slice_col]
    ppm_axis <- y_sub
    axis_label <- "F1 (ppm)"
    slice_pos <- x_sub[slice_col]
    slice_info <- paste0("F2 = ", round(slice_pos, 3), " ppm")
  }
  
  # Calculate residuals
  residuals <- exp_slice - fit_slice
  
  # Calculate slice-specific R²
  ss_res <- sum(residuals^2, na.rm = TRUE)
  ss_tot <- sum((exp_slice - mean(exp_slice, na.rm = TRUE))^2, na.rm = TRUE)
  slice_r2 <- if (ss_tot > 0) max(0, 1 - ss_res / ss_tot) else NA
  
  return(list(
    success = TRUE,
    ppm = ppm_axis,
    experimental = as.numeric(exp_slice),
    fitted = as.numeric(fit_slice),
    residuals = as.numeric(residuals),
    axis_label = axis_label,
    slice_info = slice_info,
    r_squared_global = fit_result$fit_quality,
    r_squared_slice = slice_r2,
    fit_method = fit_result$method
  ))
}


#' Create a plotly diagnostic plot from fit data
#' 
#' @param diag_data Output from generate_fit_diagnostic_data
#' @param show_residuals Logical: whether to show residuals panel
#' @return A plotly object
#' @export
create_fit_diagnostic_plot <- function(diag_data, show_residuals = TRUE) {
  
  if (!diag_data$success) {
    # Return empty plot with error message
    return(
      plotly::plot_ly() %>%
        plotly::layout(
          title = list(text = paste("⚠️", diag_data$error)),
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    )
  }
  
  # Color scheme
  col_exp <- "#1f77b4"      # Blue for experimental
  
  col_fit <- "#d62728"      # Red for fitted
  col_res <- "#2ca02c"      # Green for residuals
  
  # Quality indicator
  r2 <- diag_data$r_squared_global
  quality_text <- if (is.na(r2)) {
    "R² = NA"
  } else if (r2 >= 0.95) {
    paste0("✓ Excellent fit (R² = ", round(r2, 3), ")")
  } else if (r2 >= 0.85) {
    paste0("● Good fit (R² = ", round(r2, 3), ")")
  } else if (r2 >= 0.70) {
    paste0("◐ Acceptable fit (R² = ", round(r2, 3), ")")
  } else {
    paste0("✗ Poor fit (R² = ", round(r2, 3), ")")
  }
  
  title_text <- paste0(
    "<b>Fit Diagnostic</b> - ", diag_data$fit_method, "<br>",
    "<span style='font-size:12px'>", quality_text, " | Slice: ", diag_data$slice_info, "</span>"
  )
  
  if (show_residuals) {
    # Two-panel plot: main + residuals
    p <- plotly::subplot(
      # Main panel: experimental vs fitted
      plotly::plot_ly() %>%
        plotly::add_lines(
          x = diag_data$ppm, y = diag_data$experimental,
          name = "Experimental", line = list(color = col_exp, width = 1.5)
        ) %>%
        plotly::add_lines(
          x = diag_data$ppm, y = diag_data$fitted,
          name = "Fitted", line = list(color = col_fit, width = 2, dash = "dash")
        ) %>%
        plotly::layout(
          xaxis = list(autorange = "reversed"),
          yaxis = list(title = "Intensity")
        ),
      # Residuals panel
      plotly::plot_ly() %>%
        plotly::add_lines(
          x = diag_data$ppm, y = diag_data$residuals,
          name = "Residuals", line = list(color = col_res, width = 1),
          showlegend = FALSE
        ) %>%
        plotly::add_lines(
          x = range(diag_data$ppm), y = c(0, 0),
          line = list(color = "gray", width = 0.5, dash = "dot"),
          showlegend = FALSE
        ) %>%
        plotly::layout(
          xaxis = list(title = diag_data$axis_label, autorange = "reversed"),
          yaxis = list(title = "Residuals")
        ),
      nrows = 2, shareX = TRUE, heights = c(0.7, 0.3)
    ) %>%
      plotly::layout(
        title = list(text = title_text, x = 0.05),
        showlegend = TRUE,
        legend = list(orientation = "h", y = 1.12)
      )
  } else {
    # Single panel
    p <- plotly::plot_ly() %>%
      plotly::add_lines(
        x = diag_data$ppm, y = diag_data$experimental,
        name = "Experimental", line = list(color = col_exp, width = 1.5)
      ) %>%
      plotly::add_lines(
        x = diag_data$ppm, y = diag_data$fitted,
        name = "Fitted", line = list(color = col_fit, width = 2, dash = "dash")
      ) %>%
      plotly::layout(
        title = list(text = title_text, x = 0.05),
        xaxis = list(title = diag_data$axis_label, autorange = "reversed"),
        yaxis = list(title = "Intensity"),
        showlegend = TRUE,
        legend = list(orientation = "h", y = 1.1)
      )
  }
  
  return(p)
}


# BATCH Peak Fitting ----

#' Batch peak fitting for all boxes avec gestion robuste des erreurs
#' 
#' @param mat Spectral matrix
#' @param ppm_x F2 chemical shifts
#' @param ppm_y F1 chemical shifts
#' @param boxes Data frame of bounding boxes
#' @param model Peak model type
#' @param progress_callback Optional function for progress updates
#' @param min_points Minimum points required for fitting

calculate_fitted_volumes <- function(mat, ppm_x, ppm_y, boxes, 
                                     model = "gaussian",
                                     progress_callback = NULL,
                                     min_points = 25) {
  
  n_boxes <- nrow(boxes)
  results <- vector("list", n_boxes)
  
  # Meters for diagnostics
  n_success <- 0
  n_fallback <- 0
  n_failed <- 0
  n_multiplets <- 0
  
  for (i in seq_len(n_boxes)) {
    if (!is.null(progress_callback)) {
      progress_callback(i / n_boxes, detail = paste("Fitting box", i, "/", n_boxes))
    }
    
    box <- boxes[i, ]
    fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, 
                              model = model, 
                              min_points = min_points)
    
    # Statistics
    if (fit_result$method == model) {
      n_success <- n_success + 1
    } else if (fit_result$method == "multiplet_fit") {
      n_success <- n_success + 1
      n_multiplets <- n_multiplets + 1
    } else if (fit_result$method %in% c("sum_fit_failed", "multiplet_sum")) {
      n_fallback <- n_fallback + 1
    } else {
      n_failed <- n_failed + 1
    }
    
    results[[i]] <- data.frame(
      stain_id = box$stain_id,
      volume_fitted = fit_result$volume,
      r_squared = fit_result$fit_quality,
      center_x = if (!is.null(fit_result$params)) fit_result$params["x0"] else NA,
      center_y = if (!is.null(fit_result$params)) fit_result$params["y0"] else NA,
      fit_method = fit_result$method,
      n_peaks = ifelse(is.null(fit_result$n_peaks), 1, fit_result$n_peaks),
      is_multiplet = ifelse(is.null(fit_result$is_multiplet), FALSE, fit_result$is_multiplet),
      fit_error = ifelse(is.null(fit_result$error), NA, fit_result$error),
      stringsAsFactors = FALSE
    )
  }
  
  
  do.call(rbind, results)
}