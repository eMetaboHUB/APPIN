# =============================================================================
# Module: Session
# Description: Handles session save/load functionality (RDS format)
# =============================================================================

# MODULE UI ----

#' Session Module - UI
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing session save/load controls
#' @export
mod_session_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$details(
      tags$summary("💼 Session"),
      div(
        fluidRow(
          column(6, 
                 div(style = "margin-top: 8px;",
                     downloadButton(ns("save_session"), "💾 Save", 
                                    class = "btn-success btn-sm btn-block")
                 )
          ),
          column(6, 
                 fileInput(ns("load_session_file"), NULL, accept = ".rds", 
                           buttonLabel = "📂 Load", width = "100%")
          )
        ),
        tags$small("Save/load your complete work (peaks, boxes, parameters)", 
                   style = "color: #666;")
      )
    )
  )
}


# MODULE SERVER ----

#' Session Module - Server
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param rv List. Shared reactive values
#' @param load_data List. Return value from mod_load_data_server
#' @param refresh_nmr_plot Function. Function to refresh the NMR plot
#' @param parent_session Shiny session. Parent session for updating inputs
#' @param parent_input Shiny input. Parent input object for reading UI values
#'
#' @return NULL (side effects only)
#' @export
mod_session_server <- function(id, status_msg, rv, load_data, refresh_nmr_plot,
                                parent_session, parent_input) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # SESSION SAVE ----
    output$save_session <- downloadHandler(
      filename = function() {
        paste0("nmr_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        # Collect all session data
        session_data <- list(
          # Version for future compatibility
          version = "1.0",
          timestamp = Sys.time(),
          
          # Main data
          centroids = rv$centroids_data(),
          boxes = rv$modifiable_boxes(),
          fixed_boxes = rv$fixed_boxes(),
          
          # Pending changes
          pending_centroids = rv$pending_centroids(),
          pending_boxes = rv$pending_boxes(),
          pending_fusions = rv$pending_fusions(),
          
          # Fit results
          fit_results = rv$fit_results_data(),
          last_fit_method = rv$last_fit_method(),
          
          # UI parameters (from parent input)
          params = list(
            spectrum_type = parent_input$spectrum_type,
            contour_start = parent_input$contour_start,
            eps_value = parent_input$eps_value,
            neighborhood_size = parent_input$neighborhood_size,
            exclusion_zones = parent_input$exclusion_zones,
            filter_artifacts = parent_input$filter_artifacts,
            filter_diagonal = parent_input$filter_diagonal,
            diagonal_tolerance = parent_input$diagonal_tolerance,
            integration_method = parent_input$integration_method,
            integration_method_fit = parent_input$integration_method_fit
          ),
          
          # Folder path (for reference)
          data_path = if (!is.null(parent_input$directory)) {
            tryCatch(
              parseDirPath(c(Home = normalizePath("~"), getwd = getwd()), parent_input$directory),
              error = function(e) NULL
            )
          } else NULL
        )
        
        # Save as RDS
        saveRDS(session_data, file)
        showNotification("💾 Session saved successfully!", type = "message")
      }
    )
    
    # SESSION LOAD ----
    observeEvent(input$load_session_file, {
      req(input$load_session_file)
      
      tryCatch({
        # Load RDS
        session_data <- readRDS(input$load_session_file$datapath)
        
        # Check version
        if (is.null(session_data$version)) {
          showNotification("⚠️ Old session format, some data may not load correctly", type = "warning")
        }
        
        # Clear the caches to force a recalculation
        rv$plot_cache(list())
        rv$contour_cache(list())
        rv$box_intensity_cache(list())
        
        # Restore main data
        if (!is.null(session_data$centroids) && nrow(session_data$centroids) > 0) {
          rv$centroids_data(session_data$centroids)
        }
        
        if (!is.null(session_data$boxes) && nrow(session_data$boxes) > 0) {
          rv$modifiable_boxes(session_data$boxes)
          rv$fixed_boxes(session_data$boxes)
        }
        
        if (!is.null(session_data$fixed_boxes) && nrow(session_data$fixed_boxes) > 0) {
          rv$fixed_boxes(session_data$fixed_boxes)
        }
        
        # Restore pending changes
        if (!is.null(session_data$pending_centroids) && nrow(session_data$pending_centroids) > 0) {
          rv$pending_centroids(session_data$pending_centroids)
        }
        
        if (!is.null(session_data$pending_boxes) && nrow(session_data$pending_boxes) > 0) {
          rv$pending_boxes(session_data$pending_boxes)
        }
        
        if (!is.null(session_data$pending_fusions) && nrow(session_data$pending_fusions) > 0) {
          rv$pending_fusions(session_data$pending_fusions)
        }
        
        # Restore fit results
        if (!is.null(session_data$fit_results)) {
          rv$fit_results_data(session_data$fit_results)
        }
        
        if (!is.null(session_data$last_fit_method)) {
          rv$last_fit_method(session_data$last_fit_method)
        }
        
        # Restore UI parameters
        params <- session_data$params
        if (!is.null(params)) {
          if (!is.null(params$spectrum_type)) {
            updateSelectInput(parent_session, "spectrum_type", selected = params$spectrum_type)
          }
          if (!is.null(params$contour_start)) {
            updateNumericInput(parent_session, "contour_start", value = params$contour_start)
          }
          if (!is.null(params$eps_value)) {
            updateNumericInput(parent_session, "eps_value", value = params$eps_value)
          }
          if (!is.null(params$neighborhood_size)) {
            updateNumericInput(parent_session, "neighborhood_size", value = params$neighborhood_size)
          }
          if (!is.null(params$exclusion_zones)) {
            updateTextInput(parent_session, "exclusion_zones", value = params$exclusion_zones)
          }
          if (!is.null(params$filter_artifacts)) {
            updateCheckboxInput(parent_session, "filter_artifacts", value = params$filter_artifacts)
          }
          if (!is.null(params$filter_diagonal)) {
            updateCheckboxInput(parent_session, "filter_diagonal", value = params$filter_diagonal)
          }
          if (!is.null(params$diagonal_tolerance)) {
            updateNumericInput(parent_session, "diagonal_tolerance", value = params$diagonal_tolerance)
          }
        }
        
        # Success message with summary
        n_peaks <- if (!is.null(session_data$centroids)) nrow(session_data$centroids) else 0
        n_boxes <- if (!is.null(session_data$boxes)) nrow(session_data$boxes) else 0
        
        # Refresh plot if spectral data is loaded
        if (!is.null(load_data$bruker_data()) && !is.null(rv$contour_plot_base())) {
          refresh_nmr_plot(force_recalc = TRUE)
          
          showNotification(
            paste0("✅ Session loaded! ", n_peaks, " peaks, ", n_boxes, " boxes"),
            type = "message",
            duration = 5
          )
        } else {
          showNotification(
            paste0("✅ Session data loaded! ", n_peaks, " peaks, ", n_boxes, " boxes. ",
                   "Load spectrum data to see them on the plot."),
            type = "message",
            duration = 8
          )
        }
        
        status_msg(paste0("📂 Session loaded: ", n_peaks, " peaks, ", n_boxes, " boxes"))
        
        # Show original path if available
        if (!is.null(session_data$data_path)) {
          showNotification(
            paste0("ℹ️ Original data path: ", session_data$data_path),
            type = "message",
            duration = 8
          )
        }
        
      }, error = function(e) {
        showNotification(paste("❌ Error loading session:", e$message), type = "error")
      })
    })
    
    invisible(NULL)
  })
}
