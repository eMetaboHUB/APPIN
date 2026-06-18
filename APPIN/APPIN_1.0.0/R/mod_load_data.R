
# 2D NMR Analyst - Module: Load Data ----

# Author: Julien Guibert
# Description: Shiny module for loading Bruker NMR spectra



## Module UI ----


#' Load Data Module - UI
#'
#' Creates the UI components for the data loading section.
#' This includes a directory picker, list of available spectra with
#' checkboxes, and load button.
#'
#' @param id Character. The module's namespace ID
#' @return A tagList containing the module's UI elements
#' @export
#'
#' @examples
#' # In UI definition:
#' mod_load_data_ui("load_data")
mod_load_data_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Directory selection button
    shinyFiles::shinyDirButton(
      ns("directory"), 
      "Select Folder", 
      "Choose directory", 
      class = "btn-primary btn-sm btn-block"
    ),
    
    # Display selected path
    verbatimTextOutput(ns("selected_dir")),
    
    # Dynamic UI for available spectra
    uiOutput(ns("available_spectra_ui"))
  )
}


## Module Server ----

#' Load Data Module - Server
#'
#' Server logic for the data loading module. Handles directory selection,
#' spectrum detection, and loading of Bruker NMR data.
#'
#' @param id Character. The module's namespace ID
#' @param status_msg ReactiveVal. Shared status message reactive value
#' @param trigger_subfolder_update ReactiveVal. Trigger to notify parent 
#'   when spectra list changes (for updating subfolder selector)
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{spectra_list}: Reactive containing named list of loaded spectra
#'     \item \code{bruker_data}: Reactive containing currently selected spectrum
#'     \item \code{main_directory}: Reactive containing the selected main directory path
#'   }
#' @export
#'
#' @examples
#' # In server:
#' load_data <- mod_load_data_server("load_data", status_msg = status_msg)
#' # Access loaded spectra:
#' load_data$spectra_list()
mod_load_data_server <- function(id, status_msg, trigger_subfolder_update = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    

    # REACTIVE VALUES ----

    
    #' @description Internal reactive value storing the list of loaded spectra
    spectra_list <- reactiveVal(list())
    
    #' @description Internal reactive value storing the currently displayed spectrum
    bruker_data <- reactiveVal(NULL)
    

    # DIRECTORY SELECTION ----

    
    # Define root directories for the file browser
    roots <- c(Home = normalizePath("~"), Root = "/")
    
    # Initialize shinyFiles directory chooser
    shinyFiles::shinyDirChoose(input, "directory", roots = roots, session = session)
    
    #' Parse selected directory path
    #'
    #' Reactive that processes the shinyDirChoose selection and returns
    
    #' a normalized, validated directory path.
    #'
    #' @return Character. Normalized path to selected directory, or NULL if invalid
    main_directory <- reactive({
      req(input$directory)
      
      dir_path <- tryCatch({
        selection <- input$directory
        
        # Check that selection is valid
        if (is.null(selection) || length(selection) == 0) {
          return(NULL)
        }
        
        # Get selected root
        selected_root <- selection$root
        
        if (is.null(selected_root) || !selected_root %in% names(roots)) {
          return(NULL)
        }
        
        # Get root base path
        base_path <- roots[[selected_root]]
        
        # Get relative path (folder list)
        path_parts <- selection$path
        
        if (is.null(path_parts) || length(path_parts) == 0) {
          final_path <- base_path
        } else {
          # Filter out empty sections
          path_parts <- unlist(path_parts)
          path_parts <- path_parts[path_parts != ""]
          
          if (length(path_parts) == 0) {
            final_path <- base_path
          } else {
            # Build the complete path
            if (selected_root == "Root") {
              # For Root (/), construct the absolute path directly
              final_path <- paste0("/", paste(path_parts, collapse = "/"))
            } else {
              # For other roots, join with the base_path
              final_path <- file.path(base_path, paste(path_parts, collapse = "/"))
            }
          }
        }
        
        # Normalize the path
        norm_path <- normalizePath(final_path, mustWork = FALSE)
        
        # Check that folder exists
        if (!dir.exists(norm_path)) {
          warning(paste("Directory does not exist:", norm_path))
          return(NULL)
        }
        
        norm_path
        
      }, error = function(e) {
        warning(paste("Error parsing directory:", e$message))
        NULL
      })
      
      dir_path
    })
    

    # DISPLAY SELECTED DIRECTORY ----
    
    #' Render selected directory path
    output$selected_dir <- renderPrint({ 
      main_directory() 
    })
    

    # SUBFOLDER DETECTION ----

    
    #' Detect Bruker spectrum subfolders
    #'
    #' Reactive that scans the selected directory for valid Bruker NMR
    #' spectrum folders (containing 'acqus' and either 'ser' or 'fid' files).
    #'
    #' @return Character vector of paths to valid spectrum folders
    subfolders <- reactive({
      req(main_directory())
      
      all_subfolders <- list.dirs(main_directory(), recursive = TRUE, full.names = TRUE)
      
      # Filter to keep only valid Bruker folders
      all_subfolders[sapply(all_subfolders, function(folder) {
        file.exists(file.path(folder, "acqus")) &&
          (file.exists(file.path(folder, "ser")) || file.exists(file.path(folder, "fid")))
      })]
    })
    

    # AVAILABLE SPECTRA UI ----

    
    #' Render UI for available spectra selection
    #'
    #' Creates a dynamic UI with checkboxes for each detected spectrum,
    #' along with "Select All" / "Deselect All" buttons and a load button.
    output$available_spectra_ui <- renderUI({
      
      # Check if a main folder is selected
      if (is.null(input$directory) || length(input$directory) == 0) {
        return(
          div(
            style = "color: #6c757d; font-style: italic; padding: 10px 0;",
            "Select a folder to see available spectra."
          )
        )
      }
      
      folders <- tryCatch({
        subfolders()
      }, error = function(e) {
        return(character(0))
      })
      
      if (is.null(folders) || length(folders) == 0) {
        return(
          tags$p(
            style = "color: #6c757d; font-style: italic;", 
            "No Bruker spectra found in selected directory"
          )
        )
      }
      
      # Create short display names
      display_names <- basename(folders)
      
      # If multiple have same name, add parent folder
      if (any(duplicated(display_names))) {
        display_names <- sapply(folders, function(f) {
          paste0(basename(dirname(f)), "/", basename(f))
        })
      }
      
      tagList(
        tags$div(
          style = "margin-bottom: 10px;",
          tags$strong(paste(length(folders), "spectra found"))
        ),
        
        # Select All / Deselect All buttons
        fluidRow(
          column(6, actionButton(ns("select_all_spectra"), "✅ All", 
                                 class = "btn-sm btn-outline-success")),
          column(6, actionButton(ns("deselect_all_spectra"), "❌ None", 
                                 class = "btn-sm btn-outline-warning"))
        ),
        
        tags$hr(),
        
        # Scrollable list with checkboxes
        tags$div(
          style = "max-height: 200px; overflow-y: auto; border: 1px solid #ddd; border-radius: 5px; padding: 10px; background-color: #fff;",
          checkboxGroupInput(
            ns("spectra_to_load"),
            label = NULL,
            choices = setNames(folders, display_names),
            selected = folders  # All selected by default
          )
        ),
        
        tags$hr(),
        
        # Load button
        actionButton(ns("load_selected_spectra"), "📥 Load Selected", 
                     class = "btn-primary btn-block")
      )
    })
    

    # SELECT ALL / DESELECT ALL ----

    
    #' Handle "Select All" button click
    observeEvent(input$select_all_spectra, {
      folders <- subfolders()
      updateCheckboxGroupInput(session, "spectra_to_load", selected = folders)
    })
    
    #' Handle "Deselect All" button click
    observeEvent(input$deselect_all_spectra, {
      updateCheckboxGroupInput(session, "spectra_to_load", selected = character(0))
    })
    

    # LOAD SELECTED SPECTRA ----

    
    #' Handle spectrum loading
    #'
    #' Loads selected Bruker spectra when the load button is clicked.
    #' Shows progress bar and notifications for success/failure.
    observeEvent(input$load_selected_spectra, {
      req(input$spectra_to_load)
      
      folders <- input$spectra_to_load
      
      if (length(folders) == 0) {
        showNotification("⚠️ No spectra selected", type = "warning")
        return()
      }
      
      status_msg(paste("🔄 Loading", length(folders), "spectra..."))
      
      # Progress bar
      progress <- shiny::Progress$new()
      progress$set(message = "Loading spectra", value = 0)
      on.exit(progress$close(), add = TRUE)
      
      all_data <- list()
      
      for (i in seq_along(folders)) {
        sub <- folders[[i]]
        data_path <- file.path(sub, "pdata", "1")
        progress$inc(1 / length(folders), 
                     detail = paste0(basename(sub), " (", i, "/", length(folders), ")"))
        
        if (!dir.exists(data_path)) next
        
        data <- tryCatch({
          read_bruker_cached(data_path, dim = "2D")
        }, error = function(e) {
          showNotification(paste("❌ Error reading", basename(sub)), type = "error")
          NULL
        })
        
        if (!is.null(data)) all_data[[sub]] <- data
      }
      
      # Update reactive values
      spectra_list(all_data)
      
      # Set first spectrum as current
      if (length(all_data) > 0) {
        bruker_data(all_data[[1]])
        showNotification(paste("✅", length(all_data), "spectra loaded"), type = "message")
        status_msg(paste("✅", length(all_data), "spectra loaded"))
      } else {
        status_msg("⚠️ No valid spectra found")
      }
      
      # Notify parent that spectra list has changed
      if (!is.null(trigger_subfolder_update)) {
        trigger_subfolder_update(trigger_subfolder_update() + 1)
      }
    })
    

    # RETURN VALUES ----

    
    # Return reactive values and setters for use by other modules
    return(list(
      # Reactive getters
      spectra_list = spectra_list,
      bruker_data = bruker_data,
      main_directory = main_directory,
      subfolders = subfolders,
      
      # Setters (for parent module to update values)
      set_bruker_data = function(data) { bruker_data(data) },
      set_spectra_list = function(data) { spectra_list(data) }
    ))
  })
}