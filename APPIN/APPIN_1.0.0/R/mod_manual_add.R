# =============================================================================
# Module: Manual Add
# Description: Handles manual addition of peaks and boxes via form inputs
# =============================================================================

# MODULE UI ----

#' Manual Add Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing manual add controls
#' @export
mod_manual_add_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("➕ Add manually"),
      div(
        # Manual Peak Addition
        tags$b("Peak:"),
        fluidRow(
          column(5, numericInput(ns("manual_f2"), "F2:", value = 4.0, step = 0.01)),
          column(5, numericInput(ns("manual_f1"), "F1:", value = 3.5, step = 0.01))
        ),
        actionButton(ns("add_manual_centroid"), "Add Peak", class = "btn-info btn-sm btn-block"),
        
        hr(),
        
        # Manual Box Addition
        tags$b("Box:"),
        fluidRow(
          column(6, numericInput(ns("manual_xmin"), "xmin:", value = 3.5, step = 0.01)),
          column(6, numericInput(ns("manual_xmax"), "xmax:", value = 4.0, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput(ns("manual_ymin"), "ymin:", value = 2.0, step = 0.01)),
          column(6, numericInput(ns("manual_ymax"), "ymax:", value = 3.0, step = 0.01))
        ),
        actionButton(ns("add_manual_bbox"), "Add Box", class = "btn-info btn-sm btn-block")
      )
    )
  )
}


# MODULE SERVER ----

#' Manual Add Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values (centroids_data, modifiable_boxes, pending_*)
#' @param data_reactives List. Reactive expressions (result_data for contour_data)
#' @param peak_picking List. Return value from mod_peak_picking_server (for eps_value)
#'
#' @return NULL (side effects only)
#' @export
mod_manual_add_server <- function(id, rv, data_reactives, peak_picking) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ADD MANUAL CENTROID ----
    observeEvent(input$add_manual_centroid, {
      req(input$manual_f2, input$manual_f1)
      
      current <- rv$centroids_data()
      if (is.null(current)) {
        current <- data.frame(
          F2_ppm = numeric(0), F1_ppm = numeric(0),
          Volume = numeric(0), stain_id = character(0)
        )
      }
      
      # Generate unique manual peak ID
      existing_ids <- current$stain_id[grepl("^man", current$stain_id)]
      man_number <- if (length(existing_ids) == 0) {
        1
      } else {
        max(as.integer(sub("man", "", existing_ids)), na.rm = TRUE) + 1
      }
      
      # Estimate intensity from contour data
      contour_dat <- data_reactives$result_data()$contour_data
      eps <- peak_picking$eps_value()
      estimated_intensity <- 0
      
      if (!is.null(contour_dat) && nrow(contour_dat) > 0) {
        local_points <- contour_dat[
          abs(contour_dat$x + input$manual_f2) <= eps * 16 &
            abs(contour_dat$y + input$manual_f1) <= eps * 16, , drop = FALSE
        ]
        estimated_intensity <- sum(local_points$level, na.rm = TRUE)
      }
      
      # Create new peak entry
      new_point <- data.frame(
        F2_ppm = as.numeric(input$manual_f2),
        F1_ppm = as.numeric(input$manual_f1),
        Volume = as.numeric(estimated_intensity),
        stain_id = paste0("man", man_number),
        status = "add",
        stringsAsFactors = FALSE
      )
      
      # Ensure all columns from current are present
      missing_cols <- setdiff(colnames(current), colnames(new_point))
      for (mc in missing_cols) new_point[[mc]] <- NA
      
      rv$pending_centroids(dplyr::bind_rows(rv$pending_centroids(), new_point))
      showNotification(paste("✅ Peak added:", new_point$stain_id), type = "message")
    })
    
    # ADD MANUAL BOX ----
    observeEvent(input$add_manual_bbox, {
      req(input$manual_xmin, input$manual_xmax, input$manual_ymin, input$manual_ymax)
      
      # Validate coordinates
      if (input$manual_xmin >= input$manual_xmax) {
        showNotification("❌ xmin must be less than xmax", type = "error")
        return()
      }
      if (input$manual_ymin >= input$manual_ymax) {
        showNotification("❌ ymin must be less than ymax", type = "error")
        return()
      }
      
      current_boxes <- rv$modifiable_boxes()
      
      # Generate unique manual box ID
      existing_manual_ids <- if (!is.null(current_boxes) && nrow(current_boxes) > 0 &&
                                   "stain_id" %in% names(current_boxes)) {
        current_boxes$stain_id[grepl("^manual_box", current_boxes$stain_id)]
      } else {
        character(0)
      }
      
      manual_number <- if (length(existing_manual_ids) == 0) {
        1
      } else {
        max(as.integer(sub("manual_box", "", existing_manual_ids)), na.rm = TRUE) + 1
      }
      
      # Create new box entry
      new_box <- data.frame(
        xmin = input$manual_xmin,
        xmax = input$manual_xmax,
        ymin = input$manual_ymin,
        ymax = input$manual_ymax,
        stain_id = paste0("manual_box", manual_number),
        Volume = NA_real_,
        status = "add",
        stringsAsFactors = FALSE
      )
      
      rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), new_box))
      showNotification(paste("🟦 Box added:", new_box$stain_id), type = "message")
    })
    
    # No return value - side effects only
    invisible(NULL)
  })
}
