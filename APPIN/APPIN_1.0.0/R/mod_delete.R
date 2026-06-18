# =============================================================================
# Module: Delete Selected
# Description: Handles deletion of selected peaks and their associated boxes
# =============================================================================

# MODULE UI ----

#' Delete Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing delete controls
#' @export
mod_delete_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("🗑️ Deleting Peaks and Boxes"),
      div(
        tags$p(
          class = "help-text",
          style = "font-size: 0.85em; color: #666;",
          "Use box/lasso selection on the plot to select peaks, then click Delete."
        ),
        actionButton(ns("delete_btn"), "🗑️ Delete Selected", class = "btn-danger btn-sm btn-block")
      )
    )
  )
}


# MODULE SERVER ----

#' Delete Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values (centroids_data, modifiable_boxes, fixed_boxes, pending_deletions)
#'
#' @return NULL (side effects only)
#' @export
mod_delete_server <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # DELETE POINTS ----
    observeEvent(input$delete_btn, {
      req(rv$centroids_data())
      
      # Get selected points from plotly
      sel <- plotly::event_data("plotly_selected", source = "nmr_plot")
      
      if (is.null(sel) || nrow(sel) < 1) {
        showNotification("\u26a0\ufe0f Select at least 1 point using box or lasso selection", type = "error")
        return()
      }
      
      # Convert plotly coordinates back to ppm
      sel$x <- -sel$x
      sel$y <- -sel$y
      
      # Find matching peaks in centroids data
      brushed <- dplyr::semi_join(
        rv$centroids_data(), 
        sel, 
        by = c("F2_ppm" = "x", "F1_ppm" = "y")
      )
      
      if (nrow(brushed) < 1) {
        showNotification("\u26a0\ufe0f Selection did not match any points", type = "error")
        return()
      }
      
      # Remove selected peaks from centroids
      remaining <- dplyr::anti_join(rv$centroids_data(), brushed, by = c("F2_ppm", "F1_ppm"))
      rv$centroids_data(remaining)
      
      # Also remove associated boxes
      if (!is.null(rv$modifiable_boxes()) && nrow(rv$modifiable_boxes()) > 0) {
        boxes <- rv$modifiable_boxes()
        
        # Find boxes that contain the deleted peaks
        boxes_to_remove <- which(
          boxes$xmin <= max(brushed$F2_ppm) & boxes$xmax >= min(brushed$F2_ppm) &
            boxes$ymin <= max(brushed$F1_ppm) & boxes$ymax >= min(brushed$F1_ppm)
        )
        
        removed_boxes <- if (length(boxes_to_remove) > 0) {
          boxes[boxes_to_remove, , drop = FALSE]
        } else {
          boxes[0, ]
        }
        
        # Keep only non-matching boxes
        boxes <- if (length(boxes_to_remove) > 0) {
          boxes[-boxes_to_remove, , drop = FALSE]
        } else {
          boxes
        }
        
        rv$modifiable_boxes(boxes)
        rv$fixed_boxes(boxes)
      }
      
      # Track deletions in pending_deletions
      rv$pending_deletions(dplyr::bind_rows(rv$pending_deletions(), brushed))
      
      showNotification(
        sprintf("\u2705 %d point(s) deleted", nrow(brushed)),
        type = "message"
      )
    })
    
    # No return value - side effects only
    invisible(NULL)
  })
}
