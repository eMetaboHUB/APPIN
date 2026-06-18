# =============================================================================
# Module: Save & Export (Wrapper)
# Description: Combines all save/export sub-modules into a single interface
#
# Sub-modules:
#   - mod_session.R : Session save/load (RDS)
#   - mod_import.R  : CSV import (peaks, boxes)
#   - mod_export.R  : CSV export + batch export
#   - mod_reset.R   : Reset all data
# =============================================================================

# MODULE UI ----

#' Save & Export Module - UI (Wrapper)
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing all save/export UI components
#' @export
mod_save_export_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Session save/load
    mod_session_ui(ns("session")),
    
    # Import section
    mod_import_ui(ns("import")),
    
    # Export section
    mod_export_ui(ns("export")),
    
    # Reset button
    mod_reset_ui(ns("reset"))
  )
}


# MODULE SERVER ----

#' Save & Export Module - Server (Wrapper)
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param rv List. Named list of reactive values
#' @param load_data List. Return value from mod_load_data_server
#' @param data_reactives List. Named list of reactive expressions
#' @param refresh_nmr_plot Function. Function to refresh the NMR plot
#' @param parent_session Shiny session. Parent session for updating inputs
#' @param parent_input Shiny input. Parent input object for reading UI values
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{reset_triggered}: Reactive that increments when reset is triggered
#'   }
#' @export
mod_save_export_server <- function(id, 
                                   status_msg, 
                                   rv,
                                   load_data,
                                   data_reactives,
                                   refresh_nmr_plot,
                                   parent_session,
                                   parent_input) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize sub-modules
    
    # 1. Session (save/load)
    mod_session_server(
      id = "session",
      status_msg = status_msg,
      rv = rv,
      load_data = load_data,
      refresh_nmr_plot = refresh_nmr_plot,
      parent_session = parent_session,
      parent_input = parent_input
    )
    
    # 2. Import
    mod_import_server(
      id = "import",
      rv = rv,
      refresh_nmr_plot = refresh_nmr_plot
    )
    
    # 3. Export
    export_handle <- mod_export_server(
      id = "export",
      status_msg = status_msg,
      rv = rv,
      load_data = load_data,
      data_reactives = data_reactives
    )
    
    # 4. Reset
    reset_result <- mod_reset_server(
      id = "reset",
      status_msg = status_msg,
      rv = rv,
      parent_session = parent_session
    )
    
    # Return values from sub-modules
    return(list(
      reset_triggered     = reset_result$reset_triggered,
      shift_tolerance_ppm = export_handle$shift_tolerance_ppm,
      conflict_ids        = export_handle$conflict_ids
    ))
  })
}