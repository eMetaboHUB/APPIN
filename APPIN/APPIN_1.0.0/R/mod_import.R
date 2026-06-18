# =============================================================================
# Module: Import
# Description: Handles CSV import for peaks (centroids) and boxes
# =============================================================================

# MODULE UI ----

#' Import Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing import controls
#' @export
mod_import_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("📥 Import Peaks & Boxes"),
      div(
        fileInput(ns("import_centroids_file"), "Peaks CSV:", accept = ".csv"),
        fileInput(ns("import_boxes_file"), "Boxes CSV:", accept = ".csv")
      )
    )
  )
}


# MODULE SERVER ----

#' Import Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values
#' @param refresh_nmr_plot Function. Function to refresh the NMR plot
#'
#' @return NULL (side effects only)
#' @export
mod_import_server <- function(id, rv, refresh_nmr_plot) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Helper function to clean centroids dataframe
    clean_centroids_df <- function(df) {
      if (is.null(df)) return(NULL)
      
      # Ensure required columns exist
      required <- c("stain_id", "F2_ppm", "F1_ppm", "Volume")
      if (!all(required %in% names(df))) return(NULL)
      
      # Convert types
      df$stain_id <- as.character(df$stain_id)
      df$F2_ppm <- as.numeric(df$F2_ppm)
      df$F1_ppm <- as.numeric(df$F1_ppm)
      df$Volume <- as.numeric(df$Volume)
      
      # Remove invalid rows
      df <- df[!is.na(df$F2_ppm) & !is.na(df$F1_ppm), ]
      
      df
    }
    
    # IMPORT CENTROIDS ----
    observeEvent(input$import_centroids_file, {
      req(input$import_centroids_file)
      
      # Try semicolon first (French CSV), then comma
      imported <- tryCatch({
        df <- read.csv(input$import_centroids_file$datapath, sep = ";", stringsAsFactors = FALSE)
        # If only one column, try comma separator
        if (ncol(df) == 1) {
          df <- read.csv(input$import_centroids_file$datapath, sep = ",", stringsAsFactors = FALSE)
        }
        df
      }, error = function(e) { 
        showNotification(paste("Import error:", e$message), type = "error")
        NULL 
      })
      
      if (is.null(imported)) return()
      
      required_cols <- c("stain_id", "Volume", "F2_ppm", "F1_ppm")
      
      if (!all(required_cols %in% colnames(imported))) {
        showNotification(
          paste("❌ File must contain columns:", paste(required_cols, collapse = ", "),
                "\nFound:", paste(colnames(imported), collapse = ", ")),
          type = "error",
          duration = 10
        )
        return()
      }
      
      cleaned <- clean_centroids_df(imported)
      
      if (is.null(cleaned) || nrow(cleaned) == 0) {
        showNotification("❌ No valid data found in file", type = "error")
        return()
      }
      
      rv$centroids_data(cleaned)
      refresh_nmr_plot()
      
      showNotification(paste("✅", nrow(cleaned), "centroids imported"), type = "message")
    })
    
    # IMPORT BOXES ----
    observeEvent(input$import_boxes_file, {
      req(input$import_boxes_file)
      
      # Try semicolon first, then comma
      imported <- tryCatch({
        df <- read.csv(input$import_boxes_file$datapath, sep = ";", stringsAsFactors = FALSE)
        # If only one column, try comma separator
        if (ncol(df) == 1) {
          df <- read.csv(input$import_boxes_file$datapath, sep = ",", stringsAsFactors = FALSE)
        }
        df
      }, error = function(e) {
        showNotification(paste("Import error:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(imported)) return()
      
      required_cols <- c("stain_id", "xmin", "xmax", "ymin", "ymax")
      
      # Validate imported structure
      if (!all(required_cols %in% colnames(imported))) {
        showNotification(
          paste(
            "❌ File must contain columns:",
            paste(required_cols, collapse = ", "),
            "\nFound:", paste(colnames(imported), collapse = ", ")
          ),
          type = "error",
          duration = 10
        )
        return()
      }
      
      # Keep only required columns (+ Volume if present)
      cols_to_keep <- intersect(c(required_cols, "Volume"), colnames(imported))
      imported <- imported[, cols_to_keep, drop = FALSE]
      
      # Convert types - stain_id remains character!
      imported$stain_id <- as.character(imported$stain_id)
      imported$xmin <- as.numeric(imported$xmin)
      imported$xmax <- as.numeric(imported$xmax)
      imported$ymin <- as.numeric(imported$ymin)
      imported$ymax <- as.numeric(imported$ymax)
      
      # Check for valid data
      valid_rows <- !is.na(imported$xmin) & !is.na(imported$xmax) & 
        !is.na(imported$ymin) & !is.na(imported$ymax)
      
      if (sum(valid_rows) == 0) {
        showNotification("❌ No valid box coordinates found", type = "error")
        return()
      }
      
      imported <- imported[valid_rows, , drop = FALSE]
      
      # Update reactive values
      rv$fixed_boxes(imported)
      rv$modifiable_boxes(imported)
      rv$reference_boxes(imported)
      
      # Force plot refresh
      rv$box_intensity_cache(list())
      refresh_nmr_plot(force_recalc = TRUE)
      
      showNotification(paste("✅", nrow(imported), "bounding boxes imported"), type = "message")
    })
    
    invisible(NULL)
  })
}
