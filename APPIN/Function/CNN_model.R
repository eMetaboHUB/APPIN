# CNN_model.R - CNN Model Architecture and Loading ----
#
# Part of the CNN Peak Detection module for 2D NMR Spectra
# Supports multiple models for different spectrum types (TOCSY, HSQC, etc.)
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB


# =============================================================================
# MODEL PATHS CONFIGURATION
# =============================================================================

CNN_MODEL_PATHS <- list(
  default = "saved_model/weights",
  TOCSY   = "saved_model/weights",
  COSY    = "saved_model/weights",
  UFCOSY  = "saved_model/weights",
  HSQC    = "saved_model/weights_hsqc/weights"
)

# Cache pour stocker les modèles chargés
.cnn_models_cache <- new.env(parent = emptyenv())


# =============================================================================
# MODEL ARCHITECTURE
# =============================================================================

#' Build CNN Peak Predictor Model
#'
#' @param n_points Integer, input size (default: 2048)
#' @return A compiled Keras model
#' @export
build_peak_predictor <- function(n_points = 2048) {
  input <- layer_input(shape = c(n_points, 1), name = "input")
  
  x <- input %>%
    layer_conv_1d(filters = 40, kernel_size = 11, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 20, kernel_size = 1, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 10, kernel_size = 11, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 20, kernel_size = 1, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 10, kernel_size = 1, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 30, kernel_size = 11, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 18, kernel_size = 1, padding = "same", activation = "relu") %>%
    layer_conv_1d(filters = 18, kernel_size = 3, padding = "same", activation = "relu")
  
  output_class <- x %>%
    layer_conv_1d(filters = 3, kernel_size = 1, activation = "softmax", name = "class_output")
  
  output_reg <- x %>%
    layer_conv_1d(filters = 3, kernel_size = 1, activation = "linear", name = "reg_output")
  
  model <- keras_model(inputs = input, outputs = list(output_class, output_reg))
  
  focal_loss <- function(gamma = 3.0, alpha = 1.5) {
    function(y_true, y_pred) {
      y_true <- k_cast(y_true, "int32")
      y_true_one_hot <- k_one_hot(y_true, num_classes = 3)
      epsilon <- k_epsilon()
      y_pred <- k_clip(y_pred, epsilon, 1.0 - epsilon)
      cross_entropy <- -y_true_one_hot * k_log(y_pred)
      weight <- alpha * k_pow(1 - y_pred, gamma)
      loss <- weight * cross_entropy
      return(k_sum(loss, axis = -1))
    }
  }
  
  model %>% compile(
    loss = list(
      class_output = focal_loss(gamma = 3.0, alpha = 1.5),
      reg_output = "mse"
    ),
    loss_weights = list(class_output = 1.0, reg_output = 1.0),
    optimizer = optimizer_adam(learning_rate = 1e-4),
    metrics = list(class_output = "accuracy")
  )
  
  return(model)
}


# =============================================================================
# MODEL LOADING FUNCTIONS
# =============================================================================

#' Get CNN Model for Spectrum Type
#'
#' Automatically selects and loads the appropriate CNN model:
#' - TOCSY, COSY, UFCOSY → saved_model/weights
#' - HSQC → saved_model/weights_hsqc
#'
#' @param spectrum_type Character, one of "TOCSY", "COSY", "HSQC", "UFCOSY"
#' @param force_reload Logical, force reloading even if cached
#' @return Compiled Keras model
#' @export
get_cnn_model <- function(spectrum_type = "TOCSY", force_reload = FALSE) {
  
  spectrum_type <- toupper(spectrum_type)
  
  # Déterminer le chemin des poids
  weights_path <- if (spectrum_type %in% names(CNN_MODEL_PATHS)) {
    CNN_MODEL_PATHS[[spectrum_type]]
  } else {
    CNN_MODEL_PATHS[["default"]]
  }
  
  cache_key <- weights_path
  
  # Vérifier le cache
  if (!force_reload && exists(cache_key, envir = .cnn_models_cache)) {
    cat(sprintf("Using cached CNN model for %s\n", spectrum_type))
    return(get(cache_key, envir = .cnn_models_cache))
  }
  
  # Vérifier que les poids existent
  weights_exist <- dir.exists(weights_path) || file.exists(paste0(weights_path, ".index"))
  
  if (!weights_exist) {
    if (spectrum_type == "HSQC") {
      warning(sprintf(
        "HSQC model not found at '%s'. Using default model (trained on 1H only).",
        weights_path
      ))
      weights_path <- CNN_MODEL_PATHS[["default"]]
      cache_key <- weights_path
      
      if (!dir.exists(weights_path) && !file.exists(paste0(weights_path, ".index"))) {
        stop(sprintf("CNN model weights not found at: %s", weights_path))
      }
    } else {
      stop(sprintf("CNN model weights not found at: %s", weights_path))
    }
  }
  
  # Construire et charger le modèle
  cat(sprintf("Loading CNN model for %s from: %s\n", spectrum_type, weights_path))
  
  model <- build_peak_predictor()
  
  tryCatch({
    load_model_weights_tf(model, weights_path)
    cat(sprintf("✅ CNN model loaded successfully for %s\n", spectrum_type))
  }, error = function(e) {
    stop(sprintf("Failed to load CNN weights from %s: %s", weights_path, e$message))
  })
  
  # Mettre en cache
  assign(cache_key, model, envir = .cnn_models_cache)
  
  return(model)
}


#' Check Available CNN Models
#' @export
check_cnn_models <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              CNN MODELS AVAILABILITY                          ║\n")
  cat("╠═══════════════════════════════════════════════════════════════╣\n")
  
  std_path <- CNN_MODEL_PATHS[["default"]]
  std_exists <- dir.exists(std_path) || file.exists(paste0(std_path, ".index"))
  std_status <- if (std_exists) "✅" else "❌"
  cat(sprintf("║  %s Standard (TOCSY/COSY/UFCOSY): %s\n", std_status, std_path))
  
  hsqc_path <- CNN_MODEL_PATHS[["HSQC"]]
  hsqc_exists <- dir.exists(hsqc_path) || file.exists(paste0(hsqc_path, ".index"))
  hsqc_status <- if (hsqc_exists) "✅" else "❌"
  cat(sprintf("║  %s HSQC:                         %s\n", hsqc_status, hsqc_path))
  
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  invisible(c(standard = std_exists, hsqc = hsqc_exists))
}


#' Clear CNN Model Cache
#' @export
clear_cnn_cache <- function() {
  rm(list = ls(envir = .cnn_models_cache), envir = .cnn_models_cache)
  cat("CNN model cache cleared.\n")
}


# =============================================================================
# INITIALIZATION - Backward compatibility
# =============================================================================

# Variable globale new_model pour rétrocompatibilité avec ancien code
tryCatch({
  new_model <- get_cnn_model("TOCSY")
}, error = function(e) {
  warning("Could not load default CNN model: ", e$message)
  new_model <- NULL
})