# =============================================================================
# Module: Reset
# Description: Handles reset of all data and state
# =============================================================================

# MODULE UI ----

#' Reset Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing reset button
#' @export
mod_reset_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    hr(),
    actionButton(ns("reset_all"), "🔄 Reset All", class = "btn-outline-danger btn-sm btn-block")
  )
}


# MODULE SERVER ----

#' Reset Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param rv List. Shared reactive values
#' @param parent_session Shiny session. Parent session for updating inputs
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{reset_triggered}: ReactiveVal that increments when reset is triggered
#'   }
#' @export
mod_reset_server <- function(id, status_msg, rv, parent_session) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Counter to notify parent when reset is triggered
    reset_triggered <- reactiveVal(0)
    
    # RESET ALL ----
    observeEvent(input$reset_all, {
      
      # Reset plots
      rv$nmr_plot(NULL)
      rv$contour_plot_base(NULL)
      
      # Reset centroids/peaks
      if (!is.null(rv$imported_centroids)) rv$imported_centroids(NULL)
      rv$centroids_data(NULL)
      if (!is.null(rv$centroids)) rv$centroids(NULL)
      
      # Reset boxes - ALL box-related variables
      rv$fixed_boxes(data.frame(xmin = numeric(), xmax = numeric(), 
                                ymin = numeric(), ymax = numeric()))
      rv$modifiable_boxes(data.frame())
      rv$reference_boxes(NULL)
      
      # Reset pending changes
      rv$pending_centroids(data.frame(
        F2_ppm = numeric(0), F1_ppm = numeric(0),
        Volume = numeric(0), stain_id = character(0),
        stringsAsFactors = FALSE
      ))
      rv$pending_boxes(data.frame(
        xmin = numeric(0), xmax = numeric(0),
        ymin = numeric(0), ymax = numeric(0)
      ))
      rv$pending_fusions(data.frame(
        stain_id = character(), F2_ppm = numeric(),
        F1_ppm = numeric(), Volume = numeric(),
        stringsAsFactors = FALSE
      ))
      
      # Reset clicks and selections
      if (!is.null(rv$first_click_for_box)) rv$first_click_for_box(NULL)
      if (!is.null(rv$last_click_coords)) rv$last_click_coords(NULL)
      if (!is.null(rv$selected_box_for_edit)) rv$selected_box_for_edit(NULL)
      if (!is.null(rv$selected_box_index)) rv$selected_box_index(NULL)
      if (!is.null(rv$original_box_coords)) rv$original_box_coords(NULL)
      if (!is.null(rv$box_has_been_modified)) rv$box_has_been_modified(FALSE)
      
      # Reset fit results
      rv$fit_results_data(NULL)
      rv$last_fit_method("sum")
      
      # Reset caches
      rv$plot_cache(list())
      rv$contour_cache(list())
      rv$box_intensity_cache(list())
      
      # Reset UI
      updateSelectInput(parent_session, "selected_subfolder", selected = "")
      
      # Remove preview if present
      if (!is.null(rv$preview_trace_added) && isTRUE(rv$preview_trace_added())) {
        tryCatch({
          plotly::plotlyProxy("interactivePlot", parent_session) %>%
            plotly::plotlyProxyInvoke("deleteTraces", -1L)
        }, error = function(e) NULL)
        rv$preview_trace_added(FALSE)
      }
      
      # Notify
      status_msg("🔁 All data reset")
      showNotification("🔁 All data has been reset", type = "message")
      
      # Increment counter to notify parent
      reset_triggered(reset_triggered() + 1)
    })
    
    # Return reset trigger for parent to observe
    return(list(
      reset_triggered = reset_triggered
    ))
  })
}
