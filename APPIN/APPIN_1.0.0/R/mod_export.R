# =============================================================================
# Module: Export
# Description: Handles CSV export for peaks, boxes, and batch export
# =============================================================================

# MODULE UI ----

#' Export Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing export controls
#' @export
mod_export_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("📤 Export Data"),
      div(
        fluidRow(
          column(6, downloadButton(ns("export_centroids"), "Peaks", class = "btn-sm btn-block")),
          column(6, downloadButton(ns("export_boxes"), "Boxes", class = "btn-sm btn-block"))
        ),
        br(),
        # Shift tolerance slider for batch export
        div(
          style = "background: #fff3e0; border-radius: 8px; padding: 10px; margin-bottom: 10px; border-left: 4px solid #ff9800;",
          tags$b("⚙️ Batch Export Options", style = "color: #e65100;"),
          sliderInput(
            ns("shift_tolerance_ppm"),
            "Shift tolerance F2 (¹H, ppm):",
            min = 0, max = 0.03, value = 0, step = 0.0005
          ),
          tags$small(
            "Each box is recentered on the local maximum within ±tolerance along F2 only.",
            tags$br(),
            "Typical: 0.002-0.01 ppm. Set to 0 to disable.",
            style = "color: #666;"
          ),
          # --- Conflict preview ---
          # Shows how many boxes can collide given the chosen tolerance.
          # Each box's effective position can shift by ±tol on F2 and F1, so two
          # boxes whose coordinates differ by less than 2×tol on both axes may
          # end up overlapping during the batch export.
          # Conflicting boxes are drawn in red on the main spectrum plot.
          hr(style = "margin: 10px 0;"),
          div(
            style = "margin-top: 8px;",
            tags$b("🔍 Conflict preview", style = "color: #e65100;"),
            uiOutput(ns("conflict_summary"))
          )
        ),
        downloadButton(ns("export_batch_box_intensities"), "📤 Batch Export (all spectra)", 
                       class = "btn-primary btn-sm btn-block")
      )
    )
  )
}


# MODULE SERVER ----

#' Export Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param rv List. Shared reactive values
#' @param load_data List. Return value from mod_load_data_server
#' @param data_reactives List. Reactive expressions
#'
#' @return A list of reactives exposing:
#'   \itemize{
#'     \item \code{shift_tolerance_ppm()} — current tolerance slider value
#'     \item \code{conflict_ids()} — character vector of stain_id flagged as
#'       potentially overlapping. Use this in the main plot to color those
#'       boxes red (others green).
#'   }
#' @export
mod_export_server <- function(id, status_msg, rv, load_data, data_reactives) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # =========================================================================
    # CONFLICT DETECTION ----
    # Helper + reactives that flag boxes whose effective positions may overlap
    # after re-centering within ±shift_tolerance_ppm during batch export.
    # =========================================================================
    
    #' Detect pairs of boxes that may overlap after re-centering
    #'
    #' Two boxes can collide if, after each one is shifted by up to ±tol on
    #' each axis, their rectangles overlap. Mathematically, that means each
    #' axis-projection of box B falls within ±tol of box A's projection.
    #'
    #' @param boxes data.frame with stain_id, xmin, xmax, ymin, ymax
    #' @param tol_f2 numeric, tolerance on the F2 (1H) axis in ppm
    #' @param tol_f1 numeric, tolerance on the F1 (13C) axis in ppm
    #' @return list with `pairs` (data.frame of conflicting stain_id pairs),
    #'   `conflict_ids` (vector), `n_pairs`, `n_boxes_in_conflict`
    detect_box_overlaps <- function(boxes, tol_f2, tol_f1) {
      empty_out <- list(
        pairs = data.frame(id_a = character(0), id_b = character(0),
                           dist_F2 = numeric(0), dist_F1 = numeric(0)),
        conflict_ids = character(0),
        n_pairs = 0L,
        n_boxes_in_conflict = 0L
      )
      if (is.null(boxes) || nrow(boxes) < 2) return(empty_out)
      if (!all(c("xmin", "xmax", "ymin", "ymax", "stain_id") %in% names(boxes))) return(empty_out)
      
      # Expanded bounds: each box can move by ±tol, so its "reachable" envelope
      # extends by tol on each side. Two envelopes overlap iff their bounds
      # interleave on both axes.
      xmin_e <- boxes$xmin - tol_f2
      xmax_e <- boxes$xmax + tol_f2
      ymin_e <- boxes$ymin - tol_f1
      ymax_e <- boxes$ymax + tol_f1
      
      n <- nrow(boxes)
      pairs_list <- vector("list", 0)
      for (i in seq_len(n - 1)) {
        # Vectorised over j > i — much faster than nested loop for n ~ 200
        j_idx <- (i + 1):n
        overlap <- !(xmax_e[i] < xmin_e[j_idx] | xmin_e[i] > xmax_e[j_idx] |
                       ymax_e[i] < ymin_e[j_idx] | ymin_e[i] > ymax_e[j_idx])
        if (any(overlap)) {
          j_hits <- j_idx[overlap]
          # Centre-to-centre distance on each axis (informative, signed-magnitude)
          ci_x <- (boxes$xmin[i] + boxes$xmax[i]) / 2
          ci_y <- (boxes$ymin[i] + boxes$ymax[i]) / 2
          cj_x <- (boxes$xmin[j_hits] + boxes$xmax[j_hits]) / 2
          cj_y <- (boxes$ymin[j_hits] + boxes$ymax[j_hits]) / 2
          pairs_list[[length(pairs_list) + 1]] <- data.frame(
            id_a = boxes$stain_id[i],
            id_b = boxes$stain_id[j_hits],
            dist_F2 = abs(ci_x - cj_x),
            dist_F1 = abs(ci_y - cj_y),
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(pairs_list) == 0) return(empty_out)
      pairs_df <- do.call(rbind, pairs_list)
      conflict_ids <- unique(c(pairs_df$id_a, pairs_df$id_b))
      list(
        pairs = pairs_df,
        conflict_ids = conflict_ids,
        n_pairs = nrow(pairs_df),
        n_boxes_in_conflict = length(conflict_ids)
      )
    }
    
    # Reactive: detect overlaps whenever boxes or tolerance changes
    overlap_info <- reactive({
      boxes <- tryCatch(rv$modifiable_boxes(), error = function(e) NULL)
      tol   <- input$shift_tolerance_ppm %||% 0
      # Tolerance applies on F2 (¹H) only. F1 (¹³C) gets 0 because ¹³C peaks
      # are physically wider in ppm and far less prone to chemical shift drift
      # across spectra, so re-centering along F1 would do more harm than good.
      detect_box_overlaps(boxes, tol_f2 = tol, tol_f1 = 0)
    })
    
    # UI feedback under the slider
    output$conflict_summary <- renderUI({
      info <- overlap_info()
      tol  <- input$shift_tolerance_ppm %||% 0
      n_total <- tryCatch(nrow(rv$modifiable_boxes()), error = function(e) 0L)
      if (is.null(n_total) || is.na(n_total)) n_total <- 0L
      
      if (tol == 0) {
        return(div(style = "font-size: 11px; color: #2e7d32; margin-top: 4px;",
                   "✓ Tolerance = 0: boxes stay fixed, no overlap possible."))
      }
      if (n_total == 0) {
        return(div(style = "font-size: 11px; color: #888; margin-top: 4px;",
                   "(No boxes loaded yet)"))
      }
      if (info$n_pairs == 0) {
        return(div(style = "font-size: 11px; color: #2e7d32; margin-top: 4px;",
                   sprintf("✓ %d boxes OK at ±%.3f ppm — no potential overlap.",
                           n_total, tol)))
      }
      div(
        style = "font-size: 11px; color: #c62828; margin-top: 4px; font-weight: bold;",
        sprintf("⚠ %d box(es) involved in %d overlap pair(s) at ±%.3f ppm",
                info$n_boxes_in_conflict, info$n_pairs, tol),
        tags$br(),
        tags$small(
          style = "font-weight: normal; color: #555;",
          paste("stain_id: ",
                paste(utils::head(info$conflict_ids, 10), collapse = ", "),
                if (length(info$conflict_ids) > 10) sprintf(" … (+%d more)",
                                                            length(info$conflict_ids) - 10)
                else "")
        )
      )
    })
    
    # SHOW CONFLICT MAP MODAL — REMOVED.
    # Conflicting boxes are now shown directly on the main NMR plot in red,
    # so the modal/scatter/table machinery is unnecessary.
    
    # EXPORT CENTROIDS ----
    output$export_centroids <- downloadHandler(
      filename = function() paste0("centroids_", Sys.Date(), ".csv"),
      content = function(file) {
        df <- rv$centroids_data()
        # Use write.csv2 for ";" separator (French Excel compatible)
        if (!is.null(df) && nrow(df) > 0) {
          # Exclude Volume column from peaks export (Volume is for boxes/integration only)
          export_cols <- setdiff(names(df), c("Volume", "stain_intensity", "intensity_plot"))
          df_export <- df[, export_cols[export_cols %in% names(df)], drop = FALSE]
          write.csv2(df_export, file, row.names = FALSE)
        } else {
          write.csv2(data.frame(), file)
        }
      }
    )
    
    # EXPORT BOXES ----
    # Now uses the same calculate_batch_box_intensities() function as batch export
    # to ensure consistent intensity values across both export methods
    output$export_boxes <- downloadHandler(
      filename = function() {
        method <- data_reactives$effective_integration_method()
        method_suffix <- if (method == "sum") "" else paste0("_", method)
        paste0("box_intensities", method_suffix, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(data_reactives$bounding_boxes_data(), load_data$spectra_list())
        
        # Get chosen method (same as batch export)
        method <- data_reactives$effective_integration_method()
        model <- if (method %in% c("gaussian", "voigt")) method else "gaussian"
        
        status_msg(paste0("🔄 Calculating box intensities (", method, " method)..."))
        
        tryCatch({
          ref_boxes <- data_reactives$bounding_boxes_data()
          
          if (is.null(ref_boxes) || nrow(ref_boxes) == 0) {
            showNotification("⚠️ No boxes found", type = "warning")
            write.csv2(data.frame(message = "No boxes found."), file, row.names = FALSE)
            return(invisible(NULL))
          }
          
          # Use the same function as batch export for consistency
          box_intensities <- calculate_batch_box_intensities(
            reference_boxes = ref_boxes,
            spectra_list = load_data$spectra_list(),
            apply_shift = FALSE,
            method = method,
            model = model,
            progress = NULL  # No progress bar for single export
          )
          
          # Replace negative values with 0 (same as batch export)
          intensity_cols <- grep("^Intensity_", names(box_intensities), value = TRUE)
          for (col in intensity_cols) {
            box_intensities[[col]] <- pmax(box_intensities[[col]], 0, na.rm = TRUE)
          }
          
          # Use write.csv2 for ";" separator (French Excel compatible)
          write.csv2(box_intensities, file, row.names = FALSE)
          
          status_msg("✅ Box export complete")
          showNotification(paste("✅ Exported", nrow(ref_boxes), "boxes"), type = "message")
          
        }, error = function(e) {
          showNotification(paste("❌ Export error:", e$message), type = "error")
          status_msg(paste("❌ Error:", e$message))
        })
      }
    )
    
    # EXPORT BATCH BOX INTENSITIES ----
    output$export_batch_box_intensities <- downloadHandler(
      filename = function() {
        method <- data_reactives$effective_integration_method()
        method_suffix <- if (method == "sum") "" else paste0("_", method)
        paste0("batch_box_intensities", method_suffix, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(rv$reference_boxes(), load_data$spectra_list())
        
        # Get chosen method
        method <- data_reactives$effective_integration_method()
        model <- if (method %in% c("gaussian", "voigt")) method else "gaussian"
        
        status_msg(paste0("🔄 Calculating batch intensities (", method, " method)..."))
        
        # Progress bar
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Processing spectra", value = 0)
        
        tryCatch({
          ref_boxes <- rv$reference_boxes()
          
          if (is.null(ref_boxes) || nrow(ref_boxes) == 0) {
            showNotification("⚠️ No reference boxes found", type = "warning")
            return()
          }
          
          # Store used method
          rv$last_fit_method(method)
          
          # Get shift tolerance from UI
          shift_tol <- input$shift_tolerance_ppm %||% 0
          
          # Call batch calculation
          batch_intensities <- calculate_batch_box_intensities(
            reference_boxes = ref_boxes,
            spectra_list = load_data$spectra_list(),
            apply_shift = FALSE,
            method = method,
            model = model,
            progress = function(value, detail) {
              progress$set(value = value, detail = detail)
            },
            shift_tolerance_ppm = shift_tol
          )
          
          # Replace negative values with 0
          intensity_cols <- grep("^Intensity_", names(batch_intensities), value = TRUE)
          for (col in intensity_cols) {
            batch_intensities[[col]] <- pmax(batch_intensities[[col]], 0, na.rm = TRUE)
          }
          
          # Use write.csv2 for ";" separator (French Excel compatible)
          write.csv2(batch_intensities, file, row.names = FALSE)
          
          status_msg("✅ Batch export complete")
          showNotification(paste("✅ Exported", nrow(ref_boxes), "boxes,", 
                                 length(load_data$spectra_list()), "spectra"), type = "message")
          
        }, error = function(e) {
          showNotification(paste("❌ Export error:", e$message), type = "error")
          status_msg(paste("❌ Error:", e$message))
        })
      }
    )
    
    # Expose tolerance and conflict info so the caller (Shine.R) can color
    # conflicting boxes red on the main NMR plot.
    return(list(
      shift_tolerance_ppm = reactive({ input$shift_tolerance_ppm %||% 0 }),
      conflict_ids        = reactive({ overlap_info()$conflict_ids })
    ))
  })
}