# =============================================================================
# Module: Box Editor
# Description: Handles box selection, editing, moving, resizing, and preview
# =============================================================================

# MODULE UI ----

#' Box Editor Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing box editing controls
#' @export
mod_box_editor_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("📦 Edit selected box"),
      div(
        # Selected box info
        verbatimTextOutput(ns("selected_box_info")),
        
        # Coordinates in 2 compact lines
        div(
          class = "edit-box-inputs",
          fluidRow(
            column(6, numericInput(ns("edit_box_xmin"), "xmin:", value = NA, step = 0.01)),
            column(6, numericInput(ns("edit_box_xmax"), "xmax:", value = NA, step = 0.01))
          ),
          fluidRow(
            column(6, numericInput(ns("edit_box_ymin"), "ymin:", value = NA, step = 0.01)),
            column(6, numericInput(ns("edit_box_ymax"), "ymax:", value = NA, step = 0.01))
          )
        ),
        
        # Step and Move buttons on the same line
        div(
          style = "display: flex; align-items: flex-end; gap: 10px; margin-top: 10px;",
          
          # Step input compact
          div(
            class = "step-input-compact",
            style = "width: 70px;",
            numericInput(ns("move_box_step"), "Step:", value = 0.01, min = 0.001, step = 0.005)
          ),
          
          # Move buttons grid
          div(
            class = "move-btn-grid",
            # Line 1
            div(),
            actionButton(ns("move_box_up"), "↑", class = "btn-default btn-xs"),
            div(),
            # Line 2
            actionButton(ns("move_box_left"), "←", class = "btn-default btn-xs"),
            div(
              style = "display: flex; gap: 1px;",
              actionButton(ns("shrink_box"), "−", class = "btn-warning btn-xs"),
              actionButton(ns("expand_box"), "+", class = "btn-success btn-xs")
            ),
            actionButton(ns("move_box_right"), "→", class = "btn-default btn-xs"),
            # Line 3
            div(),
            actionButton(ns("move_box_down"), "↓", class = "btn-default btn-xs"),
            div()
          )
        ),
        
        br(),
        actionButton(ns("apply_box_edit"), "Apply Edit", class = "btn-primary btn-sm btn-block")
      )
    )
  )
}


# MODULE SERVER ----

#' Box Editor Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values
#' @param data_reactives List. Reactive expressions
#' @param parent_input Shiny input. Parent input object for table selections
#' @param parent_session Shiny session. Parent session for plotlyProxy
#'
#' @return NULL (side effects only)
#' @export
mod_box_editor_server <- function(id, rv, data_reactives, parent_input, parent_session) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Null-coalescing operator
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    # SELECTED BOX INFO OUTPUT ----
    output$selected_box_info <- renderText({
      box <- rv$selected_box_for_edit()
      modified <- rv$box_has_been_modified()
      if (is.null(box)) return("No box selected")
      status <- if (modified) " (modified)" else ""
      sprintf("Editing: %s%s", box$stain_id, status)
    })
    
    # BOX SELECTION HANDLER (from table) ----
    observeEvent(parent_input$bbox_table_rows_selected, {
      selected <- parent_input$bbox_table_rows_selected
      
      if (length(selected) > 0 && !is.null(selected)) {
        first_selected <- selected[1]
        boxes <- data_reactives$bounding_boxes_data()
        
        if (!is.null(boxes) && nrow(boxes) >= first_selected) {
          box <- boxes[first_selected, , drop = FALSE]
          rv$selected_box_for_edit(box)
          rv$selected_box_index(first_selected)
          rv$box_has_been_modified(FALSE)
          
          # Store original coordinates
          rv$original_box_coords(list(
            xmin = box$xmin, xmax = box$xmax,
            ymin = box$ymin, ymax = box$ymax
          ))
          
          # Update input fields
          updateNumericInput(session, "edit_box_xmin", value = round(box$xmin, 4))
          updateNumericInput(session, "edit_box_xmax", value = round(box$xmax, 4))
          updateNumericInput(session, "edit_box_ymin", value = round(box$ymin, 4))
          updateNumericInput(session, "edit_box_ymax", value = round(box$ymax, 4))
          
          # Add preview trace (dashed green box)
          x0 <- -box$xmin; x1 <- -box$xmax
          y0 <- -box$ymin; y1 <- -box$ymax
          if (x0 > x1) { tmp <- x0; x0 <- x1; x1 <- tmp }
          if (y0 > y1) { tmp <- y0; y0 <- y1; y1 <- tmp }
          
          plotly::plotlyProxy("interactivePlot", parent_session) %>%
            plotly::plotlyProxyInvoke(
              "addTraces",
              list(
                x = c(x0, x1, x1, x0, x0),
                y = c(y0, y0, y1, y1, y0),
                type = "scatter", mode = "lines",
                line = list(color = "lime", width = 3, dash = "dash"),
                hoverinfo = "text",
                text = paste("Preview:", box$stain_id),
                showlegend = FALSE,
                name = "preview_box"
              )
            )
          
          rv$preview_trace_added(TRUE)
        }
      } else {
        # Deselection - remove preview
        if (isTRUE(rv$preview_trace_added())) {
          plotly::plotlyProxy("interactivePlot", parent_session) %>%
            plotly::plotlyProxyInvoke("deleteTraces", -1L)
          rv$preview_trace_added(FALSE)
        }
        
        rv$selected_box_for_edit(NULL)
        rv$selected_box_index(NULL)
        rv$original_box_coords(NULL)
        rv$box_has_been_modified(FALSE)
      }
    })
    
    # MOVE BOX HANDLERS ----
    observeEvent(input$move_box_up, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_ymin", value = input$edit_box_ymin - delta)
      updateNumericInput(session, "edit_box_ymax", value = input$edit_box_ymax - delta)
      rv$box_has_been_modified(TRUE)
    })
    
    observeEvent(input$move_box_down, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_ymin", value = input$edit_box_ymin + delta)
      updateNumericInput(session, "edit_box_ymax", value = input$edit_box_ymax + delta)
      rv$box_has_been_modified(TRUE)
    })
    
    observeEvent(input$move_box_left, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_xmin", value = input$edit_box_xmin + delta)
      updateNumericInput(session, "edit_box_xmax", value = input$edit_box_xmax + delta)
      rv$box_has_been_modified(TRUE)
    })
    
    observeEvent(input$move_box_right, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_xmin", value = input$edit_box_xmin - delta)
      updateNumericInput(session, "edit_box_xmax", value = input$edit_box_xmax - delta)
      rv$box_has_been_modified(TRUE)
    })
    
    # RESIZE BOX HANDLERS ----
    observeEvent(input$expand_box, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_xmin", value = input$edit_box_xmin - delta)
      updateNumericInput(session, "edit_box_xmax", value = input$edit_box_xmax + delta)
      updateNumericInput(session, "edit_box_ymin", value = input$edit_box_ymin - delta)
      updateNumericInput(session, "edit_box_ymax", value = input$edit_box_ymax + delta)
      rv$box_has_been_modified(TRUE)
    })
    
    observeEvent(input$shrink_box, {
      req(rv$selected_box_for_edit())
      delta <- input$move_box_step %||% 0.01
      updateNumericInput(session, "edit_box_xmin", value = input$edit_box_xmin + delta)
      updateNumericInput(session, "edit_box_xmax", value = input$edit_box_xmax - delta)
      updateNumericInput(session, "edit_box_ymin", value = input$edit_box_ymin + delta)
      updateNumericInput(session, "edit_box_ymax", value = input$edit_box_ymax - delta)
      rv$box_has_been_modified(TRUE)
    })
    
    # PREVIEW UPDATE ON INPUT CHANGE ----
    observeEvent(c(input$edit_box_xmin, input$edit_box_xmax, 
                   input$edit_box_ymin, input$edit_box_ymax), {
      req(rv$selected_box_index())
      req(isTRUE(rv$preview_trace_added()))
      
      x0 <- -input$edit_box_xmin; x1 <- -input$edit_box_xmax
      y0 <- -input$edit_box_ymin; y1 <- -input$edit_box_ymax
      if (x0 > x1) { tmp <- x0; x0 <- x1; x1 <- tmp }
      if (y0 > y1) { tmp <- y0; y0 <- y1; y1 <- tmp }
      
      # Update preview trace
      plotly::plotlyProxy("interactivePlot", parent_session) %>%
        plotly::plotlyProxyInvoke("deleteTraces", -1L) %>%
        plotly::plotlyProxyInvoke(
          "addTraces",
          list(
            x = c(x0, x1, x1, x0, x0),
            y = c(y0, y0, y1, y1, y0),
            type = "scatter", mode = "lines",
            line = list(color = "lime", width = 3, dash = "dash"),
            hoverinfo = "text", text = "Preview (modified)",
            showlegend = FALSE, name = "preview_box"
          )
        )
      
      # Check if modified
      original <- rv$original_box_coords()
      if (!is.null(original)) {
        if (abs(input$edit_box_xmin - original$xmin) > 1e-6 ||
            abs(input$edit_box_xmax - original$xmax) > 1e-6 ||
            abs(input$edit_box_ymin - original$ymin) > 1e-6 ||
            abs(input$edit_box_ymax - original$ymax) > 1e-6) {
          rv$box_has_been_modified(TRUE)
        }
      }
    }, ignoreInit = TRUE)
    
    # APPLY BOX EDIT ----
    observeEvent(input$apply_box_edit, {
      box_to_edit <- rv$selected_box_for_edit()
      if (is.null(box_to_edit)) {
        showNotification("⚠️ Select a box first", type = "warning")
        return()
      }
      
      new_xmin <- input$edit_box_xmin
      new_xmax <- input$edit_box_xmax
      new_ymin <- input$edit_box_ymin
      new_ymax <- input$edit_box_ymax
      
      # Validation
      coords_valid <- TRUE
      error_msg <- ""
      
      if (is.na(new_xmin) || is.na(new_xmax) || is.na(new_ymin) || is.na(new_ymax)) {
        coords_valid <- FALSE
        error_msg <- "Coordinates cannot be empty (NA)"
      } else if (new_xmin == 0 && new_xmax == 0 && new_ymin == 0 && new_ymax == 0) {
        coords_valid <- FALSE
        error_msg <- "All coordinates cannot be zero. Use 'Delete' to remove a box."
      } else if (new_xmin >= new_xmax) {
        coords_valid <- FALSE
        error_msg <- "xmin must be less than xmax"
      } else if (new_ymin >= new_ymax) {
        coords_valid <- FALSE
        error_msg <- "ymin must be less than ymax"
      } else if ((new_xmax - new_xmin) < 0.001 || (new_ymax - new_ymin) < 0.001) {
        coords_valid <- FALSE
        error_msg <- "Box is too small (min size: 0.001 ppm)"
      }
      
      if (!coords_valid) {
        showNotification(paste("❌ Invalid coordinates:", error_msg), type = "error", duration = 5)
        return()
      }
      
      if (!rv$box_has_been_modified()) {
        showNotification("ℹ️ No changes to apply", type = "message")
        # Clean up preview
        if (rv$preview_trace_added()) {
          plotly::plotlyProxy("interactivePlot", parent_session) %>%
            plotly::plotlyProxyInvoke("deleteTraces", list(-1))
          rv$preview_trace_added(FALSE)
        }
        rv$selected_box_for_edit(NULL)
        rv$selected_box_index(NULL)
        rv$original_box_coords(NULL)
        rv$box_has_been_modified(FALSE)
        return()
      }
      
      # Create edited box entry
      edited_box <- data.frame(
        xmin = new_xmin, xmax = new_xmax,
        ymin = new_ymin, ymax = new_ymax,
        stain_id = box_to_edit$stain_id,
        Volume = NA_real_,
        status = "edit",
        original_stain_id = box_to_edit$stain_id,
        stringsAsFactors = FALSE
      )
      
      rv$pending_boxes(dplyr::bind_rows(rv$pending_boxes(), edited_box))
      
      # Clean up preview
      if (isTRUE(rv$preview_trace_added())) {
        plotly::plotlyProxy("interactivePlot", parent_session) %>%
          plotly::plotlyProxyInvoke("deleteTraces", -1L)
        rv$preview_trace_added(FALSE)
      }
      
      # Reset state
      rv$selected_box_for_edit(NULL)
      rv$selected_box_index(NULL)
      rv$original_box_coords(NULL)
      rv$box_has_been_modified(FALSE)
      
      showNotification(paste("✏️ Box edit pending:", box_to_edit$stain_id), type = "message")
    })
    
    # No return value - side effects only
    invisible(NULL)
  })
}
