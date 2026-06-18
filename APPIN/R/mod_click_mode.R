# =============================================================================
# Module: Click Mode
# Description: Handles click modes (add peak, box selection, delete box)
#              and click coordinates display
# =============================================================================

# MODULE UI ----

#' Click Mode Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing click mode controls
#' @export
mod_click_mode_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("🖱️ Click mode"),
      div(
        radioButtons(
          ns("box_click_mode"), NULL,
          choices = c(
            "Off" = "disabled",
            "Add peak (1 click)" = "add_peak",
            "Add box (selection)" = "box_select",
            "Delete box on click" = "delete_click"
          ),
          selected = "disabled",
          inline = FALSE
        ),
        
        # Add peak mode indicator
        conditionalPanel(
          condition = sprintf("input['%s'] == 'add_peak'", ns("box_click_mode")),
          div(
            class = "info-box",
            style = "background-color: #e8f5e9; padding: 8px; margin: 5px 0; border-radius: 4px;",
            icon("crosshairs"),
            tags$b(" Add peak mode active"), br(),
            tags$small("Click on the spectrum to add a centroid")
          )
        ),
        
        # Box select mode indicator
        conditionalPanel(
          condition = sprintf("input['%s'] == 'box_select'", ns("box_click_mode")),
          div(
            class = "info-box",
            style = "background-color: #e3f2fd; padding: 8px; margin: 5px 0; border-radius: 4px;",
            icon("vector-square"),
            tags$b(" Box selection mode active"), br(),
            tags$small("Use the box select tool "), 
            icon("square", class = "fa-regular"),
            tags$small(" in the plot toolbar, then drag to create a box")
          )
        ),
        
        # Delete mode warning
        conditionalPanel(
          condition = sprintf("input['%s'] == 'delete_click'", ns("box_click_mode")),
          div(
            class = "warning-box",
            style = "padding: 8px; margin: 5px 0;",
            icon("exclamation-triangle"),
            tags$b(" Delete mode active"), br(),
            tags$small("Click inside a box to mark it for deletion")
          )
        ),
        
        # Click coordinates display
        div(class = "click-coords", textOutput(ns("click_coords_display")))
      )
    )
  )
}


# MODULE SERVER ----

#' Click Mode Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values (first_click_for_box, last_click_coords, 
#'           pending_boxes, pending_centroids, centroids_data)
#' @param data_reactives List. Reactive expressions (bounding_boxes_data, result_data)
#' @param peak_picking List. Return value from mod_peak_picking_server (for eps_value)
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{box_click_mode}: Reactive returning current click mode
#'   }
#' @export
mod_click_mode_server <- function(id, rv, data_reactives, peak_picking = NULL) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Track last processed selection to avoid duplicates
    last_selection_id <- reactiveVal(NULL)
    
    # CLICK COORDINATES DISPLAY ----
    output$click_coords_display <- renderText({
      click_data <- plotly::event_data("plotly_click", source = "nmr_plot")
      if (is.null(click_data)) return("Click on the spectrum")
      sprintf("F2 = %.4f ppm, F1 = %.4f ppm", -click_data$x, -click_data$y)
    })
    
    # BOX SELECTION HANDLING via plotly_relayout ----
    # This captures the selection rectangle coordinates directly,
    # independent of any data points or invisible grid
    observeEvent(plotly::event_data("plotly_relayout", source = "nmr_plot"), {
      click_mode <- input$box_click_mode
      
      # Only process if in box_select mode
      if (is.null(click_mode) || click_mode != "box_select") return()
      
      relayout_data <- plotly::event_data("plotly_relayout", source = "nmr_plot")
      
      if (is.null(relayout_data)) return()
      
      # Check if selections exists and is a data.frame
      if (!"selections" %in% names(relayout_data)) return()
      
      sel <- relayout_data$selections
      if (!is.data.frame(sel) || nrow(sel) == 0) return()
      
      # Extract coordinates from the first selection
      x0 <- sel$x0[1]
      x1 <- sel$x1[1]
      y0 <- sel$y0[1]
      y1 <- sel$y1[1]
      
      # Validate we have all coordinates
      if (is.null(x0) || is.null(x1) || is.null(y0) || is.null(y1)) return()
      if (is.na(x0) || is.na(x1) || is.na(y0) || is.na(y1)) return()
      
      # Create a unique ID for this selection to avoid duplicates
      selection_id <- paste(round(x0, 4), round(x1, 4), round(y0, 4), round(y1, 4), sep = "_")
      if (identical(selection_id, last_selection_id())) return()
      last_selection_id(selection_id)
      
      # Convert to positive ppm values (axes are inverted)
      f2_min <- -max(x0, x1)
      f2_max <- -min(x0, x1)
      f1_min <- -max(y0, y1)
      f1_max <- -min(y0, y1)
      
      # Validate that we have a real box
      if (abs(f2_max - f2_min) < 0.001 || abs(f1_max - f1_min) < 0.001) {
        showNotification("⚠️ Selection too small, please draw a larger box", type = "warning")
        return()
      }
      
      # Create new box
      new_box <- data.frame(
        xmin = f2_min,
        xmax = f2_max,
        ymin = f1_min,
        ymax = f1_max,
        stain_id = paste0("sel_box_", format(Sys.time(), "%H%M%S")),
        Volume = NA_real_,
        status = "add",
        stringsAsFactors = FALSE
      )
      
      rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), new_box))
      
      showNotification(
        sprintf("🟦 Box created: F2=[%.3f, %.3f], F1=[%.3f, %.3f]",
                f2_min, f2_max, f1_min, f1_max),
        type = "message", duration = 3
      )
    }, ignoreInit = TRUE, ignoreNULL = FALSE)
    
    # CLICK HANDLING ----
    observeEvent(plotly::event_data("plotly_click", source = "nmr_plot"), {
      click_data <- plotly::event_data("plotly_click", source = "nmr_plot")
      
      # Store last click coordinates
      if (!is.null(click_data) && !is.null(click_data$x) && !is.null(click_data$y) &&
          !is.na(click_data$x) && !is.na(click_data$y)) {
        rv$last_click_coords(list(F2_ppm = click_data$x, F1_ppm = click_data$y))
      }
      
      if (is.null(click_data$x) || is.null(click_data$y)) return()
      if (is.na(click_data$x) || is.na(click_data$y)) return()
      
      f2_ppm <- -click_data$x
      f1_ppm <- -click_data$y
      click_mode <- input$box_click_mode
      
      # ADD_PEAK mode
      if (!is.null(click_mode) && click_mode == "add_peak") {
        
        current <- rv$centroids_data()
        if (is.null(current)) {
          current <- data.frame(
            F2_ppm = numeric(0), F1_ppm = numeric(0),
            Volume = numeric(0), stain_id = character(0)
          )
        }
        
        # Generate unique click peak ID
        existing_ids <- current$stain_id[grepl("^click", current$stain_id)]
        pending <- rv$pending_centroids()
        if (!is.null(pending) && nrow(pending) > 0 && "stain_id" %in% names(pending)) {
          existing_ids <- c(existing_ids, pending$stain_id[grepl("^click", pending$stain_id)])
        }
        
        click_number <- if (length(existing_ids) == 0) {
          1
        } else {
          max(as.integer(sub("click", "", existing_ids)), na.rm = TRUE) + 1
        }
        
        # Estimate intensity from contour data
        estimated_intensity <- 0
        eps <- if (!is.null(peak_picking) && !is.null(peak_picking$eps_value)) {
          peak_picking$eps_value()
        } else {
          0.04  # default fallback
        }
        
        if (!is.null(data_reactives$result_data)) {
          contour_dat <- data_reactives$result_data()$contour_data
          if (!is.null(contour_dat) && nrow(contour_dat) > 0) {
            local_points <- contour_dat[
              abs(contour_dat$x + f2_ppm) <= eps * 16 &
                abs(contour_dat$y + f1_ppm) <= eps * 16, , drop = FALSE
            ]
            estimated_intensity <- sum(local_points$level, na.rm = TRUE)
          }
        }
        
        # Create new peak entry
        new_point <- data.frame(
          F2_ppm = as.numeric(f2_ppm),
          F1_ppm = as.numeric(f1_ppm),
          Volume = as.numeric(estimated_intensity),
          stain_id = paste0("click", click_number),
          status = "add",
          stringsAsFactors = FALSE
        )
        
        # Ensure all columns from current are present
        missing_cols <- setdiff(colnames(current), colnames(new_point))
        for (mc in missing_cols) new_point[[mc]] <- NA
        
        rv$pending_centroids(dplyr::bind_rows(rv$pending_centroids(), new_point))
        showNotification(
          sprintf("📍 Peak added: %s (F2=%.3f, F1=%.3f)", 
                  new_point$stain_id, f2_ppm, f1_ppm),
          type = "message", duration = 3
        )
        return()
      }
      
      # DELETE mode
      if (!is.null(click_mode) && click_mode == "delete_click") {
        boxes <- data_reactives$bounding_boxes_data()
        
        if (is.null(boxes) || nrow(boxes) == 0) {
          showNotification("⚠️ No boxes to delete", type = "warning")
          return()
        }
        
        # Find clicked box
        clicked_box_idx <- which(
          boxes$xmin <= f2_ppm & boxes$xmax >= f2_ppm &
            boxes$ymin <= f1_ppm & boxes$ymax >= f1_ppm
        )
        
        if (length(clicked_box_idx) == 0) {
          showNotification("⚠️ No box at this location", type = "warning")
          return()
        }
        
        # If multiple boxes overlap, select the smallest
        if (length(clicked_box_idx) > 1) {
          box_areas <- (boxes$xmax[clicked_box_idx] - boxes$xmin[clicked_box_idx]) *
            (boxes$ymax[clicked_box_idx] - boxes$ymin[clicked_box_idx])
          clicked_box_idx <- clicked_box_idx[which.min(box_areas)]
        }
        
        box_to_delete <- boxes[clicked_box_idx, , drop = FALSE]
        box_to_delete$status <- "delete"
        rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), box_to_delete))
        
        showNotification(
          paste("🗑️ Box", box_to_delete$stain_id, "marked for deletion. Click 'Apply' to confirm."),
          type = "message"
        )
        return()
      }
    }, priority = 10)
    
    # RETURN VALUES ----
    return(list(
      box_click_mode = reactive({ input$box_click_mode })
    ))
  })
}