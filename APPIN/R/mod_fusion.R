# =============================================================================
# Module: Fusion
# Description: Handles fusion of selected peaks and their associated boxes
# =============================================================================

# MODULE UI ----

#' Fusion Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing fusion controls
#' @export
mod_fusion_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("🔗 Fusing Peaks and Boxes"),
      div(
        tags$p(
          class = "help-text",
          style = "font-size: 0.85em; color: #666;",
          "Use box/lasso selection on the plot to select peaks, then click Fuse."
        ),
        actionButton(ns("fuse_btn"), "🔗 Fuse Selected", class = "btn-warning btn-sm btn-block")
      )
    )
  )
}


# MODULE SERVER ----

#' Fusion Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param rv List. Shared reactive values (centroids_data, modifiable_boxes, fixed_boxes, pending_fusions)
#'
#' @return NULL (side effects only)
#' @export
mod_fusion_server <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # FUSE POINTS ----
    observeEvent(input$fuse_btn, {
      req(rv$centroids_data())
      
      # Get selected points from plotly
      sel <- plotly::event_data("plotly_selected", source = "nmr_plot")
      
      if (is.null(sel) || nrow(sel) < 2) {
        showNotification("⚠️ Select at least 2 points using box or lasso selection", type = "error")
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
      
      if (nrow(brushed) < 2) {
        showNotification("⚠️ Selection did not match enough points", type = "error")
        return()
      }
      
      # Create fused peak (weighted centroid)
      first_peak_id <- brushed$stain_id[1]
      peak_number <- gsub("[^0-9]", "", first_peak_id)
      
      # Calculate centroid weighted by volume
      total_volume <- sum(as.numeric(brushed$Volume), na.rm = TRUE)
      
      if (total_volume > 0) {
        # Weighted average
        fused_f2 <- sum(brushed$F2_ppm * as.numeric(brushed$Volume), na.rm = TRUE) / total_volume
        fused_f1 <- sum(brushed$F1_ppm * as.numeric(brushed$Volume), na.rm = TRUE) / total_volume
      } else {
        # Simple average if no volume info
        fused_f2 <- mean(brushed$F2_ppm)
        fused_f1 <- mean(brushed$F1_ppm)
      }
      
      fused_point <- data.frame(
        stain_id = paste0("fused_point", peak_number),
        F2_ppm = fused_f2,
        F1_ppm = fused_f1,
        Volume = total_volume,
        stringsAsFactors = FALSE
      )
      
      # Remove fused peaks from centroids and add new fused peak
      remaining <- dplyr::anti_join(rv$centroids_data(), brushed, by = c("F2_ppm", "F1_ppm"))
      
      # Ensure column compatibility
      missing_cols <- setdiff(names(remaining), names(fused_point))
      for (mc in missing_cols) fused_point[[mc]] <- NA
      fused_point <- fused_point[, names(remaining), drop = FALSE]
      
      rv$centroids_data(rbind(remaining, fused_point))
      
      # Also fuse associated boxes
      if (!is.null(rv$modifiable_boxes()) && nrow(rv$modifiable_boxes()) > 0) {
        boxes <- rv$modifiable_boxes()
        
        # Find boxes that overlap with the fused peaks region
        selected_boxes <- which(
          boxes$xmin <= max(brushed$F2_ppm) & boxes$xmax >= min(brushed$F2_ppm) &
            boxes$ymin <= max(brushed$F1_ppm) & boxes$ymax >= min(brushed$F1_ppm)
        )
        
        removed_boxes <- if (length(selected_boxes) > 0) {
          boxes[selected_boxes, , drop = FALSE]
        } else {
          boxes[0, ]
        }
        
        # Remove selected boxes
        boxes <- if (length(selected_boxes) > 0) {
          boxes[-selected_boxes, , drop = FALSE]
        } else {
          boxes
        }
        
        # Create new merged box
        if (nrow(removed_boxes) > 0) {
          new_box <- data.frame(
            xmin = min(removed_boxes$xmin),
            xmax = max(removed_boxes$xmax),
            ymin = min(removed_boxes$ymin),
            ymax = max(removed_boxes$ymax),
            stain_id = paste0("bbox_fused_point", peak_number)
          )
          
          # Ensure column compatibility
          missing_box_cols <- setdiff(names(boxes), names(new_box))
          for (mc in missing_box_cols) new_box[[mc]] <- NA
          new_box <- new_box[, names(boxes), drop = FALSE]
          
          boxes <- rbind(boxes, new_box)
        }
        
        rv$modifiable_boxes(boxes)
        rv$fixed_boxes(boxes)
      }
      
      # Track fusion in pending_fusions
      rv$pending_fusions(dplyr::bind_rows(rv$pending_fusions(), fused_point))
      
      showNotification(
        sprintf("✅ %d points fused into %s", nrow(brushed), fused_point$stain_id),
        type = "message"
      )
    })
    
    # No return value - side effects only
    invisible(NULL)
  })
}
