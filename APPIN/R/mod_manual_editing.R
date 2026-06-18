# =============================================================================
# Module: Manual Editing (Wrapper)
# Description: Combines all manual editing sub-modules into a single interface
# 
# Sub-modules:
#   - mod_click_mode.R      : Click modes (add peak, two-click box, delete-click)
#   - mod_box_editor.R      : Box editing (move, resize, preview)
#   - mod_manual_add.R      : Manual addition of peaks and boxes
#   - mod_fusion.R          : Peak and box fusion
#   - mod_pending_changes.R : Apply/Discard workflow
# =============================================================================

# MODULE UI ----

#' Manual Editing Module - UI (Wrapper)
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing all manual editing UI components
#' @export
mod_manual_editing_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Click mode section
    mod_click_mode_ui(ns("click_mode")),
    
    # Fusion section
    mod_fusion_ui(ns("fusion")),
    
    mod_delete_ui(ns("delete")), 
    
    # Box editor section
    mod_box_editor_ui(ns("box_editor")),
    
    # Manual add section
    mod_manual_add_ui(ns("manual_add")),
    
    # Apply/Discard buttons
    mod_pending_changes_ui(ns("pending"))
  )
}


# MODULE SERVER ----

#' Manual Editing Module - Server (Wrapper)
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param load_data List. Return value from mod_load_data_server
#' @param rv List. Named list of reactive values
#' @param data_reactives List. Named list of reactive expressions
#' @param refresh_nmr_plot Function. Function to refresh the NMR plot
#' @param peak_picking List. Return value from mod_peak_picking_server
#' @param parent_input Shiny input. Parent input object for table selections
#' @param parent_session Shiny session. Parent session for plotlyProxy
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{box_click_mode}: Reactive returning the current click mode
#'   }
#' @export
mod_manual_editing_server <- function(id,
                                      status_msg,
                                      load_data,
                                      rv,
                                      data_reactives,
                                      refresh_nmr_plot,
                                      peak_picking,
                                      parent_input,
                                      parent_session) {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize sub-modules
    
    # 1. Click Mode
    click_mode_result <- mod_click_mode_server(
      id = "click_mode",
      rv = rv,
      data_reactives = data_reactives,
      peak_picking = peak_picking
    )
    
    # 2. Fusion
    mod_fusion_server(
      id = "fusion",
      rv = rv
    )
    
    # 2b. Delete Selected
    mod_delete_server(
      id = "delete",
      rv = rv
    )
    
    # 3. Box Editor
    mod_box_editor_server(
      id = "box_editor",
      rv = rv,
      data_reactives = data_reactives,
      parent_input = parent_input,
      parent_session = parent_session
    )
    
    # 4. Manual Add
    mod_manual_add_server(
      id = "manual_add",
      rv = rv,
      data_reactives = data_reactives,
      peak_picking = peak_picking
    )
    
    # 5. Pending Changes
    mod_pending_changes_server(
      id = "pending",
      rv = rv,
      data_reactives = data_reactives,
      load_data = load_data,
      refresh_nmr_plot = refresh_nmr_plot,
      parent_input = parent_input,
      parent_session = parent_session
    )
    
    # Return values from sub-modules
    return(list(
      box_click_mode = click_mode_result$box_click_mode
    ))
  })
}