# =============================================================================
# Module: Pending Changes
# Description: Handles Apply/Discard workflow and all pending operations
#              including table-based deletions and discards
# =============================================================================

# MODULE UI ----

#' Pending Changes Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing apply/discard buttons
#' @export
mod_pending_changes_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(6, actionButton(ns("apply_changes"), "✅ Apply", class = "btn-success btn-sm btn-block")),
      column(6, actionButton(ns("discard_changes"), "❌ Discard", class = "btn-secondary btn-sm btn-block"))
    )
  )
}


# MODULE SERVER ----

#' Pending Changes Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values
#' @param data_reactives List. Reactive expressions
#' @param load_data List. Return value from mod_load_data_server (for bruker_data)
#' @param refresh_nmr_plot Function. Function to refresh the NMR plot
#' @param parent_input Shiny input. Parent input object for table selections
#' @param parent_session Shiny session. Parent session for plotlyProxy
#'
#' @return NULL (side effects only)
#' @export
mod_pending_changes_server <- function(id, rv, data_reactives, load_data, 
                                       refresh_nmr_plot, parent_input, parent_session) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Helper function to get box intensity
    get_box_intensity <- function(mat, ppm_x, ppm_y, box) {
      if (is.null(mat)) return(NA_real_)
      x_idx <- which(ppm_x >= box$xmin & ppm_x <= box$xmax)
      y_idx <- which(ppm_y >= box$ymin & ppm_y <= box$ymax)
      if (length(x_idx) == 0 || length(y_idx) == 0) return(NA_real_)
      sum(mat[y_idx, x_idx], na.rm = TRUE)
    }
    
    # Helper to reset pending data frames
    reset_pending <- function() {
      rv$pending_centroids(data.frame(
        F2_ppm = numeric(0), F1_ppm = numeric(0),
        Volume = numeric(0), stain_id = character(0),
        stringsAsFactors = FALSE
      ))
      rv$pending_boxes(data.frame(
        xmin = numeric(0), xmax = numeric(0),
        ymin = numeric(0), ymax = numeric(0),
        stain_id = character(0), Volume = numeric(0), status = character(0),
        stringsAsFactors = FALSE
      ))
      rv$pending_fusions(data.frame(
        stain_id = character(0), F2_ppm = numeric(0),
        F1_ppm = numeric(0), Volume = numeric(0),
        stringsAsFactors = FALSE
      ))
      rv$pending_deletions(data.frame(
        stain_id = character(0), F2_ppm = numeric(0),
        F1_ppm = numeric(0), Volume = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Helper to clean up edit state
    cleanup_edit_state <- function() {
      if (isTRUE(rv$preview_trace_added())) {
        tryCatch({
          plotly::plotlyProxy("interactivePlot", parent_session) %>%
            plotly::plotlyProxyInvoke("deleteTraces", -1L)
        }, error = function(e) NULL)
        rv$preview_trace_added(FALSE)
      }
      
      rv$selected_box_for_edit(NULL)
      rv$selected_box_index(NULL)
      rv$original_box_coords(NULL)
      rv$box_has_been_modified(FALSE)
    }
    
    # =========================================================================
    # APPLY CHANGES
    # =========================================================================
    
    observeEvent(input$apply_changes, {
      
      # Show spinner with Apply message
      shinyjs::runjs('
        document.getElementById("spinner_message").innerHTML = "✅ Applying changes...";
        document.getElementById("plot_spinner").style.display = "flex";
      ')
      
      current_centroids <- rv$centroids_data()
      pending_cents <- rv$pending_centroids()
      
      # --- Handle pending centroids ---
      if (!is.null(pending_cents) && nrow(pending_cents) > 0) {
        if (is.null(current_centroids)) current_centroids <- data.frame()
        
        if ("status" %in% names(pending_cents)) {
          cents_to_add <- pending_cents[is.na(pending_cents$status) | 
                                          pending_cents$status != "delete", , drop = FALSE]
          cents_to_delete <- pending_cents[!is.na(pending_cents$status) & 
                                             pending_cents$status == "delete", , drop = FALSE]
          
          # Process deletions
          if (nrow(cents_to_delete) > 0 && nrow(current_centroids) > 0) {
            if ("stain_id" %in% names(cents_to_delete) && "stain_id" %in% names(current_centroids)) {
              ids_to_delete <- cents_to_delete$stain_id
              current_centroids <- current_centroids[!current_centroids$stain_id %in% ids_to_delete, , drop = FALSE]
            }
          }
          
          # Process additions
          if (nrow(cents_to_add) > 0) {
            cents_to_add$status <- NULL
            current_centroids <- dplyr::bind_rows(current_centroids, cents_to_add)
          }
        } else {
          current_centroids <- dplyr::bind_rows(current_centroids, pending_cents)
        }
        
        rv$centroids_data(current_centroids)
      }
      
      # --- Handle pending boxes ---
      current_boxes <- rv$modifiable_boxes()
      pending_bxs <- rv$pending_boxes()
      
      if (!is.null(pending_bxs) && nrow(pending_bxs) > 0) {
        if (!"status" %in% names(pending_bxs)) pending_bxs$status <- "add"
        pending_bxs$status[is.na(pending_bxs$status)] <- "add"
        
        boxes_to_add <- pending_bxs[pending_bxs$status == "add", , drop = FALSE]
        boxes_to_delete <- pending_bxs[pending_bxs$status == "delete", , drop = FALSE]
        boxes_to_edit <- pending_bxs[pending_bxs$status == "edit", , drop = FALSE]
        
        if (is.null(current_boxes)) {
          current_boxes <- data.frame(
            xmin = numeric(0), xmax = numeric(0),
            ymin = numeric(0), ymax = numeric(0),
            stain_id = character(0), Volume = numeric(0),
            stringsAsFactors = FALSE
          )
        }
        
        # Process deletions
        if (nrow(boxes_to_delete) > 0 && nrow(current_boxes) > 0) {
          ids_to_delete <- boxes_to_delete$stain_id
          current_boxes <- current_boxes[!current_boxes$stain_id %in% ids_to_delete, , drop = FALSE]
        }
        
        # Process edits
        if (nrow(boxes_to_edit) > 0 && nrow(current_boxes) > 0) {
          mat <- load_data$bruker_data()$spectrumData
          ppm_x <- if (!is.null(mat)) suppressWarnings(as.numeric(colnames(mat))) else NULL
          ppm_y <- if (!is.null(mat)) suppressWarnings(as.numeric(rownames(mat))) else NULL
          
          for (i in seq_len(nrow(boxes_to_edit))) {
            edit_row <- boxes_to_edit[i, ]
            original_id <- if ("original_stain_id" %in% names(edit_row) && !is.na(edit_row$original_stain_id)) {
              edit_row$original_stain_id
            } else {
              edit_row$stain_id
            }
            
            box_idx <- which(current_boxes$stain_id == original_id)
            
            if (length(box_idx) > 0) {
              current_boxes[box_idx, "xmin"] <- edit_row$xmin
              current_boxes[box_idx, "xmax"] <- edit_row$xmax
              current_boxes[box_idx, "ymin"] <- edit_row$ymin
              current_boxes[box_idx, "ymax"] <- edit_row$ymax
              
              # Recalculate volume
              if (!is.null(mat)) {
                current_boxes[box_idx, "Volume"] <- get_box_intensity(
                  mat, ppm_x, ppm_y, current_boxes[box_idx, , drop = FALSE]
                )
              }
            }
          }
        }
        
        # Process additions
        if (nrow(boxes_to_add) > 0) {
          if (!"stain_id" %in% names(boxes_to_add) || any(is.na(boxes_to_add$stain_id))) {
            boxes_to_add$stain_id <- paste0("box_", seq_len(nrow(boxes_to_add)))
          }
          
          # Calculate volumes for new boxes
          mat <- load_data$bruker_data()$spectrumData
          if (!is.null(mat)) {
            ppm_x <- suppressWarnings(as.numeric(colnames(mat)))
            ppm_y <- suppressWarnings(as.numeric(rownames(mat)))
            boxes_to_add$Volume <- sapply(seq_len(nrow(boxes_to_add)), function(i) {
              get_box_intensity(mat, ppm_x, ppm_y, boxes_to_add[i, , drop = FALSE])
            })
          }
          
          # Remove status columns
          cols_to_remove <- c("status", "original_stain_id")
          boxes_to_add <- boxes_to_add[, !names(boxes_to_add) %in% cols_to_remove, drop = FALSE]
          
          # Align columns
          all_cols <- unique(c(names(current_boxes), names(boxes_to_add)))
          for (col in all_cols) {
            if (!col %in% names(current_boxes)) current_boxes[[col]] <- NA
            if (!col %in% names(boxes_to_add)) boxes_to_add[[col]] <- NA
          }
          boxes_to_add <- boxes_to_add[, names(current_boxes), drop = FALSE]
          current_boxes <- rbind(current_boxes, boxes_to_add)
        }
        
        # Clean up status columns from final data
        cols_to_clean <- c("status", "original_stain_id")
        for (col in cols_to_clean) {
          if (col %in% names(current_boxes)) current_boxes[[col]] <- NULL
        }
        
        rv$modifiable_boxes(current_boxes)
        rv$fixed_boxes(current_boxes)
        rv$reference_boxes(current_boxes)
      }
      
      # Reset pending and refresh
      reset_pending()
      rv$box_intensity_cache(list())
      refresh_nmr_plot(force_recalc = TRUE)
      
      # Update spinner message - will be hidden when plot is rendered
      shinyjs::runjs('document.getElementById("spinner_message").innerHTML = "📊 Updating plot...";')
      
      showNotification("✅ Changes applied", type = "message")
    })
    
    # =========================================================================
    # DISCARD CHANGES
    # =========================================================================
    
    observeEvent(input$discard_changes, {
      reset_pending()
      cleanup_edit_state()
      showNotification("❌ All pending changes discarded", type = "warning")
    })
    
    # =========================================================================
    # DELETE SELECTED PEAKS (from Data tab table)
    # =========================================================================
    
    observeEvent(parent_input$delete_selected_peaks, {
      selected_rows <- parent_input$centroid_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No peaks selected. Use Ctrl+Click to select multiple.", type = "warning")
        return()
      }
      
      df <- rv$centroids_data()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      to_delete <- df[selected_rows, , drop = FALSE]
      to_delete$status <- "delete"
      
      rv$pending_centroids(dplyr::bind_rows(rv$pending_centroids(), to_delete))
      
      showNotification(
        paste("🗑️", n_to_delete, "peak(s) marked for deletion. Click 'Apply' to confirm."),
        type = "message"
      )
    })
    
    # =========================================================================
    # DELETE SELECTED BOXES (from Data tab table)
    # =========================================================================
    
    observeEvent(parent_input$delete_selected_boxes, {
      selected_rows <- parent_input$bbox_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No boxes selected. Use Ctrl+Click to select multiple.", type = "warning")
        return()
      }
      
      df <- data_reactives$bounding_boxes_data()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      to_delete <- df[selected_rows, , drop = FALSE]
      to_delete$status <- "delete"
      
      rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), to_delete))
      
      showNotification(
        paste("🗑️", n_to_delete, "box(es) marked for deletion. Click 'Apply' to confirm."),
        type = "message"
      )
    })
    
    # =========================================================================
    # DISCARD SELECTED PENDING CENTROIDS
    # =========================================================================
    
    observeEvent(parent_input$discard_selected_centroid, {
      selected_rows <- parent_input$pending_centroids_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No pending peaks selected", type = "warning")
        return()
      }
      
      df <- rv$pending_centroids()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      peaks_to_discard <- df[selected_rows, , drop = FALSE]
      df_remaining <- df[-selected_rows, , drop = FALSE]
      rv$pending_centroids(df_remaining)
      
      # Build informative message
      if ("status" %in% names(peaks_to_discard)) {
        n_delete <- sum(peaks_to_discard$status == "delete", na.rm = TRUE)
        n_add <- sum(is.na(peaks_to_discard$status) | peaks_to_discard$status != "delete")
        msg_parts <- c()
        if (n_delete > 0) msg_parts <- c(msg_parts, paste(n_delete, "deletion(s) cancelled"))
        if (n_add > 0) msg_parts <- c(msg_parts, paste(n_add, "addition(s) cancelled"))
        showNotification(paste("↩️", paste(msg_parts, collapse = ", ")), type = "message")
      } else {
        showNotification(paste("🗑️ Removed", n_to_delete, "pending peak(s)"), type = "message")
      }
    })
    
    # =========================================================================
    # DISCARD SELECTED PENDING BOXES
    # =========================================================================
    
    observeEvent(parent_input$discard_selected_box, {
      selected_rows <- parent_input$pending_boxes_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No pending boxes selected", type = "warning")
        return()
      }
      
      df <- rv$pending_boxes()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      boxes_to_discard <- df[selected_rows, , drop = FALSE]
      df_remaining <- df[-selected_rows, , drop = FALSE]
      rv$pending_boxes(df_remaining)
      
      cleanup_edit_state()
      
      # Build informative message
      if ("status" %in% names(boxes_to_discard)) {
        n_delete <- sum(boxes_to_discard$status == "delete", na.rm = TRUE)
        n_add <- sum(boxes_to_discard$status == "add", na.rm = TRUE)
        n_edit <- sum(boxes_to_discard$status == "edit", na.rm = TRUE)
        msg_parts <- c()
        if (n_delete > 0) msg_parts <- c(msg_parts, paste(n_delete, "deletion(s) cancelled"))
        if (n_add > 0) msg_parts <- c(msg_parts, paste(n_add, "addition(s) cancelled"))
        if (n_edit > 0) msg_parts <- c(msg_parts, paste(n_edit, "edit(s) cancelled"))
        showNotification(paste("↩️", paste(msg_parts, collapse = ", ")), type = "message")
      } else {
        showNotification(paste("🗑️ Removed", n_to_delete, "pending box(es)"), type = "message")
      }
    })
    
    # =========================================================================
    # DISCARD SELECTED PENDING FUSIONS
    # =========================================================================
    
    observeEvent(parent_input$discard_selected_fusion, {
      selected_rows <- parent_input$pending_fusions_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No pending fusions selected", type = "warning")
        return()
      }
      
      df <- rv$pending_fusions()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      df_remaining <- df[-selected_rows, , drop = FALSE]
      rv$pending_fusions(df_remaining)
      
      showNotification(paste("🗑️ Removed", n_to_delete, "pending fusion(s)"), type = "message")
    })
    
    # =========================================================================
    # DISCARD SELECTED PENDING DELETIONS
    # =========================================================================
    
    observeEvent(parent_input$discard_selected_deletion, {
      selected_rows <- parent_input$pending_deletions_table_rows_selected
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        showNotification("⚠️ No pending deletions selected", type = "warning")
        return()
      }
      
      df <- rv$pending_deletions()
      if (is.null(df) || nrow(df) == 0) return()
      
      n_to_delete <- length(selected_rows)
      df_remaining <- df[-selected_rows, , drop = FALSE]
      rv$pending_deletions(df_remaining)
      
      showNotification(paste("↩️ Restored", n_to_delete, "deleted point(s)"), type = "message")
    })
    
    # =========================================================================
    # DELETE CENTROID (legacy single selection)
    # =========================================================================
    
    observeEvent(parent_input$delete_centroid, {
      selected <- parent_input$centroid_table_rows_selected
      if (length(selected) > 0) {
        current <- rv$centroids_data()
        to_delete <- current[selected, , drop = FALSE]
        to_delete$status <- "delete"
        rv$pending_centroids(dplyr::bind_rows(rv$pending_centroids(), to_delete))
        rv$centroids_data(current[-selected, , drop = FALSE])
        showNotification("🗑️ Centroid marked for deletion", type = "message")
      } else {
        showNotification("⚠️ Select a centroid first", type = "warning")
      }
    })
    
    # =========================================================================
    # DELETE BOX (legacy single selection)
    # =========================================================================
    
    observeEvent(parent_input$delete_bbox, {
      selected <- parent_input$bbox_table_rows_selected
      if (length(selected) > 0) {
        current <- rv$modifiable_boxes()
        to_delete <- current[selected, , drop = FALSE]
        to_delete$status <- "delete"
        rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), to_delete))
        showNotification(paste("🗑️ Box marked for deletion:", to_delete$stain_id[1]), type = "message")
      } else {
        showNotification("⚠️ Select a box first", type = "warning")
      }
    })
    
    # No return value - side effects only
    invisible(NULL)
  })
}