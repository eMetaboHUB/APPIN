#
#                                 APPIN ----
#                          SHINY APPLICATION ----
#

# Application: APPIN
# Author:      Julien Guibert
# Institution: INRAe Toxalim / MetaboHub
# GitHub:      https://github.com/JulienGuibertTlse3/2DNMR-Analyst
# Date : 14/04/2026
# 
# Description:
#   Interactive Shiny application for analyzing 2D NMR spectra from Bruker 
#   instruments. Features include automated peak detection, interactive 
#   bounding box manipulation, peak integration, and batch processing.
#

# FILE STRUCTURE:
#   Section 0: Dependencies (Packages & Source Files)
#   Section 1: User Interface
#   Section 2: Server Logic



## SECTION 0: DEPENDENCIES ----                             


### --- 0.1 Packages ----

# Core Shiny
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyBS) 
library(shinycssloaders)
library(shinydashboard)
library(shinydashboardPlus)

# Data visualization
library(plotly)
library(ggplot2)
library(viridis)
library(DT)

# Data manipulation
library(dplyr)
library(data.table)
library(reshape2)
library(readr)
library(abind)

# Signal processing
library(sp)
library(dbscan)
library(zoo)
library(matrixStats)
library(minpack.lm)
library(imager)

# Optional: TensorFlow/Keras for CNN peak detection
library(tensorflow)
library(keras)
library(reticulate)

### --- 0.2 Source Files ----

# Core functions (Function/ folder)
source("Function/Read_2DNMR_spectrum.R")
source("Function/Vizualisation.R")
source("Function/Peak_picking.R")
source("Function/Peak_fitting.R")
source("Function/CNN_shiny.R")

# Shiny modules (R/ folder)
source("R/utils.R")
source("R/mod_load_data.R")
source("R/mod_session.R")
source("R/mod_import.R")
source("R/mod_export.R")
source("R/mod_reset.R")
source("R/mod_save_export.R")
source("R/mod_peak_picking.R")
source("R/mod_integration.R")
source("R/mod_click_mode.R")
source("R/mod_box_editor.R")
source("R/mod_manual_add.R")
source("R/mod_fusion.R")
source("R/mod_pending_changes.R")
source("R/mod_manual_editing.R") 
source("R/mod_delete.R")  


## SECTION 1: USER INTERFACE ----                        


ui <- fluidPage(
  useShinyjs(),
  
  # External CSS and JavaScript
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$script(src = "plotly_ticks.js"),
    # Spinner animation CSS
    tags$style(HTML("
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
      .manual-spinner {
        border: 4px solid #f3f3f3;
        border-top: 4px solid #3c8dbc;
        border-radius: 50%;
        width: 50px;
        height: 50px;
        animation: spin 1s linear infinite;
        margin: 0 auto;
      }
      #plot_spinner {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        z-index: 1000;
        background: rgba(255,255,255,0.9);
        border-radius: 4px;
      }
      #spinner_message {
        margin-top: 15px;
        color: #3c8dbc;
        font-weight: 500;
        font-size: 14px;
      }
      .plot-container {
        position: relative;
        min-height: 550px;
      }
    "))
  ),
  
  ## DASHBOARD ----
  dashboardPage(
    
    # Header
    dashboardHeader(title = "APPIN"),
    
    # Sidebar
    dashboardSidebar(
      width = 280,
      sidebarMenu(
        id = "main_menu",
        menuItem("📖 Guide", tabName = "guide", icon = icon("book")),
        menuItem("📊 Analysis", tabName = "analysis", icon = icon("chart-area"))
        
      )
      
    ),
    
    # Body
    dashboardBody(
      tabItems(
        
        ### TAB: GUIDE ----
        tabItem(
          tabName = "guide",
          fluidRow(
            column(12,
                   div(
                     class = "guide-container",
                     style = "padding: 20px;",
                     
                     # Styled header
                     div(
                       style = "display: flex; align-items: center; margin-bottom: 20px; padding-bottom: 15px; border-bottom: 2px solid #3c8dbc;",
                       icon("book-open", class = "fa-2x", style = "color: #3c8dbc; margin-right: 15px;"),
                       h2("User Guide", style = "margin: 0; color: #333; font-weight: 600;")
                       
                     ),
                     
                     # Content with improved style
                     div(
                       class = "guide-content",
                       style = "background: #fff; border-radius: 8px; padding: 25px; box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                       
                       uiOutput("toolDescription")
                       
                     )
                     
                   )
                   
            )
            
          )
          
        ),
        
        ### TAB: ANALYSIS ----
        tabItem(
          tabName = "analysis",
          fluidRow(
            
            #### LEFT PANEL - Controls (width = 3) ----
            column(3,
                   
                   # ACCORDION - Only one section open
                   bsCollapse(
                     id = "main_accordion",
                     open = "panel_load",
                     multiple = FALSE,
                     
                     ##### ===== SECTION 1: LOAD DATA =====
                     bsCollapsePanel(
                       title = "📂 1. Load Data",
                       value = "panel_load",
                       style = "primary",
                       
                       # Module UI
                       mod_load_data_ui("load_data")
                     ),
                     
                     ##### ===== SECTION 2: PLOT SETTINGS =====
                     bsCollapsePanel(
                       title = "📈 2. Plot Settings",
                       value = "panel_plot",
                       style = "primary",
                       
                       uiOutput("subfolder_selector"),
                       
                       # === SPECTRUM TYPE SELECTOR (HIGHLIGHTED) ===
                       div(
                         style = "background: linear-gradient(135deg, #fff3cd 0%, #ffeeba 100%); 
                                  border: 2px solid #ffc107; 
                                  border-radius: 8px; 
                                  padding: 10px 12px; 
                                  margin-bottom: 15px;
                                  box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
                         div(
                           style = "display: flex; align-items: center; margin-bottom: 8px;",
                           tags$span(style = "font-size: 18px; margin-right: 8px;", "⚠️"),
                           tags$b("Spectrum Type", style = "color: #856404; font-size: 14px;")
                         ),
                         selectInput("spectrum_type", NULL, 
                                     choices = c("TOCSY", "HSQC", "HMBC", "COSY", "UFCOSY", "JRES"),
                                     selected = "TOCSY",
                                     width = "100%"),
                         tags$small(
                           style = "color: #856404; font-style: italic;",
                           "Please verify the correct spectrum type before analysis"
                         )
                       ),
                       
                       fluidRow(
                         column(8, 
                                numericInput("contour_start", "Threshold:", value = 80000, min = 0, step = 1000)
                                
                         ),
                         column(4,
                                div(style = "padding-top: 25px;",
                                    actionButton("calculate_contour", "Auto", class = "btn-info btn-sm")
                                    
                                )
                                
                         )
                         
                       ),
                       tags$details(
                         tags$summary("⚙️ Advanced"),
                         div(
                           selectInput("seuil_method", NULL, 
                                       choices = c("% of max" = "max_pct", "Noise ×" = "bruit_mult")),
                           conditionalPanel("input.seuil_method == 'max_pct'",
                                            numericInput("pct_val", "Percentage:", value = 0.0001, min = 0.001, max = 1, step = 0.001),
                                            # Warning spécifique HSQC
                                            conditionalPanel(
                                              "input.spectrum_type == 'HSQC'",
                                              div(
                                                style = "background: #f8d7da; 
                                                         border: 1px solid #f5c6cb; 
                                                         border-radius: 5px; 
                                                         padding: 8px 10px; 
                                                         margin-top: 8px;",
                                                tags$span(style = "font-size: 14px;", "⚠️"),
                                                tags$small(
                                                  style = "color: #721c24; margin-left: 5px;",
                                                  HTML("<b>HSQC:</b> '% of max' may give inconsistent thresholds due to varying 
                                                       <sup>13</sup>C intensity ranges. Consider using 'Noise ×' method instead.")
                                                )
                                              )
                                            )
                           ),
                           conditionalPanel("input.seuil_method == 'bruit_mult'",
                                            numericInput("bruit_mult", "Multiplier:", value = 1, min = 0.5, max = 10, step = 0.5)
                                            
                           )
                           
                         )
                         
                       ),
                       br(),
                       actionButton("generate_plot", "📊 Generate Plot", class = "btn-primary btn-block")
                       
                     ),
                     
                     ##### ===== SECTION 3: PEAK PICKING =====
                     bsCollapsePanel(
                       title = "🎯 3. Peak Picking",
                       value = "panel_peaks",
                       style = "primary",
                       
                       # Module UI
                       mod_peak_picking_ui("peak_picking")
                     ),
                     
                     ##### ===== SECTION 4: MANUAL EDITING =====
                     bsCollapsePanel(
                       title = "✏️ 4. Manual Editing",
                       value = "panel_edit",
                       style = "primary",
                       # Module UI
                       mod_manual_editing_ui("manual_editing")
                     ),
                     
                     ##### ===== SECTION 5: INTEGRATION =====
                     bsCollapsePanel(
                       title = "📐 5. Integration",
                       value = "panel_integration",
                       style = "primary",
                       
                       # Module UI
                       mod_integration_ui("integration")
                     ),
                     
                     ##### ===== SECTION 6: SAVE & EXPORT =====
                     bsCollapsePanel(
                       title = "💾 6. Save & Export",
                       value = "panel_export",
                       style = "primary",
                       
                       # Module UI
                       mod_save_export_ui("save_export")
                     )
                     
                   ) # End bsCollapse
                   
            ),
            
            #### RIGHT PANEL - Visualization (width = 9) ----
            column(9,
                   
                   # Status message
                   div(class = "info-box-custom",
                       textOutput("status_message")
                       
                   ),
                   
                   # Main tabs
                   tabBox(
                     width = 12,
                     id = "main_tabs",
                     
                     ##### Tab 1: Spectrum plot ----
                     tabPanel(
                       title = tagList(icon("chart-area"), "Spectrum"),
                       value = "spectrum_tab",
                       div(
                         class = "plot-container",
                         # Hidden spinner - shown only during loading
                         shinyjs::hidden(
                           div(
                             id = "plot_spinner",
                             div(class = "manual-spinner"),
                             div(id = "spinner_message", "Loading spectrum...")
                           )
                         ),
                         # Plot always visible
                         plotlyOutput("interactivePlot", height = "550px", width = "100%")
                       ),
                       conditionalPanel(
                         "output.has_pending_changes",
                         div(class = "warning-box",
                             tags$b("⏳ Pending: "),
                             textOutput("pending_summary", inline = TRUE)
                             
                         )
                         
                       )
                       
                     ),
                     
                     ##### Tab 2: Data tables ----
                     tabPanel(
                       title = tagList(icon("table"), "Data"),
                       value = "data_tab",
                       fluidRow(
                         column(6,
                                h4("🔴 Peaks"),
                                DTOutput("centroid_table"),
                                div(style = "margin-top: 8px;",
                                    actionButton("delete_selected_peaks", "🗑️ Delete Selected", 
                                                 class = "btn-sm btn-danger"),
                                    tags$small(" (Ctrl+Click for multiple)", style = "color: #666; margin-left: 10px;")
                                    
                                )
                                
                         ),
                         column(6,
                                h4("🟦 Boxes"),
                                DTOutput("bbox_table"),
                                div(style = "margin-top: 8px;",
                                    actionButton("delete_selected_boxes", "🗑️ Delete Selected", 
                                                 class = "btn-sm btn-danger"),
                                    tags$small(" (Ctrl+Click for multiple)", style = "color: #666; margin-left: 10px;")
                                    
                                )
                                
                         )
                         
                       )
                       
                     ),
                     
                     ##### Tab 3: Pending changes ----
                     tabPanel(
                       title = tagList(icon("clock"), "Pending"),
                       value = "pending_tab",
                       
                       # Peaks - first row
                       div(style = "margin-bottom: 20px;",
                           h4("🔴 Pending Peaks"),
                           DTOutput("pending_centroids_table"),
                           actionButton("discard_selected_centroid", "Remove Selected", class = "btn-sm btn-danger", style = "margin-top: 5px;")
                           
                       ),
                       hr(),
                       
                       # Boxes - second row
                       div(style = "margin-bottom: 20px;",
                           h4("🟦 Pending Boxes"),
                           DTOutput("pending_boxes_table"),
                           actionButton("discard_selected_box", "Remove Selected", class = "btn-sm btn-danger", style = "margin-top: 5px;")
                           
                       ),
                       hr(),
                       
                       # Fusions - third row
                       div(style = "margin-bottom: 20px;",
                           h4("🔗 Pending Fusions"),
                           DTOutput("pending_fusions_table"),
                           actionButton("discard_selected_fusion", "Remove Selected", class = "btn-sm btn-danger", style = "margin-top: 5px;")
                           
                       ),
                       hr(),
                       
                       # Deletions - fourth row
                       div(style = "margin-bottom: 20px;",
                           h4("🗑️ Pending Deletions"),
                           DTOutput("pending_deletions_table"),
                           actionButton("discard_selected_deletion", "Remove Selected",
                                        class = "btn-sm btn-danger", style = "margin-top: 5px;")
                       )
                       
                     ),
                     
                     # Dans tabBox, ajouter un nouvel onglet
                     tabPanel(
                       title = tagList(icon("chart-line"), "Fit Quality"),
                       value = "fit_quality_tab",
                       
                       # Info box en haut
                       div(class = "info-box-custom",
                           style = "margin-bottom: 20px;",
                           icon("info-circle"),
                           " This tab displays fit quality metrics when using Gaussian or Voigt integration methods. ",
                           "Select a box in the ", tags$b("Data"), " tab to see detailed fit visualization."
                           
                       ),
                       
                       # Global summary
                       fluidRow(
                         column(6,
                                h4("📊 Fit Summary by Method"),
                                DTOutput("fit_summary_table")
                                
                         ),
                         column(6,
                                h4("📋 Fitted Boxes Details"),
                                DTOutput("fit_boxes_detail_table")
                                
                         )
                         
                       ),
                       br(),
                       
                       # Distribution des R²
                       fluidRow(
                         column(12,
                                h4("📈 R² Distribution"),
                                plotlyOutput("fit_quality_plot", height = "400px")
                                
                         )
                         
                       ),
                       br(),
                       
                       # Detailed visualization of a selected box
                       fluidRow(
                         column(12,
                                h4("🔍 Selected Box - 2D Fit"),
                                div(style = "margin-bottom: 10px;",
                                    sliderInput("fit_zoom_level", "Zoom level:", 
                                                min = 0, max = 100, value = 80, step = 5, post = "%",
                                                width = "300px"),
                                    tags$small("0% = max context (zoomed out) | 100% = tight on box (zoomed in)", style = "color: #666;")
                                ),
                                div(style = "border: 1px solid #ddd; border-radius: 8px; padding: 15px; background: #fafafa;",
                                    plotOutput("example_fit_2d", height = "550px")
                                )
                         )
                       ),
                       br(),
                       fluidRow(
                         column(12,
                                h4("📉 Residuals Distribution"),
                                div(style = "border: 1px solid #ddd; border-radius: 8px; padding: 10px; background: #fafafa;",
                                    plotOutput("residuals_plot", height = "350px")
                                    
                                )
                                
                         )
                         
                       ),
                       br(),
                       
                       # Tips
                       div(style = "background: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;",
                           h5("💡 Interpretation Tips"),
                           tags$ul(
                             tags$li(tags$b("R² > 0.9:"), " Excellent fit - peak is well-defined"),
                             tags$li(tags$b("R² 0.7-0.9:"), " Good fit - acceptable quantification"),
                             tags$li(tags$b("R² < 0.7:"), " Poor fit - consider manual inspection or sum method"),
                             tags$li(tags$b("Residuals:"), " Should be randomly distributed around zero")
                             
                           )
                           
                       )
                       
                     )
                     
                   )
                   
            )
            
          )
          
        )
        
      )
      
    )
    
  )
  
)


## SECTION 2: SERVER LOGIC ----                               


server <- function(input, output, session) {
  
  
  ### 2.1 CONFIGURATION ----
  
  
  options(future.globals.maxSize = 10000 * 1024^2)  # Allow up to 10 GB
  future_available <- requireNamespace("future", quietly = TRUE) && 
    requireNamespace("promises", quietly = TRUE)
  
  
  ### 2.2 MODULE: LOAD DATA ----
  
  
  trigger_subfolder_update <- reactiveVal(0)
  load_data <- mod_load_data_server(
    "load_data", 
    status_msg = status_msg,
    trigger_subfolder_update = trigger_subfolder_update
  )
  
  
  ### 2.3 REACTIVE VALUES ----
  
  
  # --- General state ---
  status_msg <- reactiveVal("")
  bruker_data <- reactiveVal(NULL)
  spectra_list <- reactiveVal(list())
  spectra_plots <- reactiveVal(list())
  
  # --- Plots and cache ---
  plot_cache <- reactiveVal(list())
  contour_cache <- reactiveVal(list())
  box_intensity_cache <- reactiveVal(list())
  nmr_plot <- reactiveVal(NULL)
  contour_plot_base <- reactiveVal(NULL)
  result_data <- reactiveVal(NULL)
  result_data_list <- reactiveVal(list())
  
  # --- Centroids ---
  centroids <- reactiveVal(NULL)
  centroids_data <- reactiveVal(NULL)
  imported_centroids <- reactiveVal(NULL)
  
  # --- Boxes ---
  fixed_boxes <- reactiveVal(data.frame(xmin = numeric(), xmax = numeric(),
                                        ymin = numeric(), ymax = numeric()))
  modifiable_boxes <- reactiveVal(data.frame())
  reference_boxes <- reactiveVal()
  
  # Tolerance preview — populated by mod_export via save_export$conflict_ids().
  # Boxes whose stain_id is in this vector are drawn in red on the NMR plot.
  conflict_ids_rv <- reactiveVal(character(0))
  
  # Current shift tolerance (ppm) — populated from mod_export's slider. Used by
  # refresh_nmr_plot to draw a dashed "search envelope" around each box that
  # represents the ±tol zone within which it may be re-centered during batch.
  shift_tol_rv <- reactiveVal(0)
  
  # --- Pending changes ---
  pending_centroids <- reactiveVal(data.frame(
    F2_ppm = numeric(0), F1_ppm = numeric(0),
    Volume = numeric(0), stain_id = character(0),
    stringsAsFactors = FALSE
  ))
  pending_boxes <- reactiveVal(data.frame(
    xmin = numeric(0), xmax = numeric(0),
    ymin = numeric(0), ymax = numeric(0)
  ))
  pending_fusions <- reactiveVal(data.frame(
    stain_id = character(), F2_ppm = numeric(),
    F1_ppm = numeric(), Volume = numeric(),
    stringsAsFactors = FALSE
  ))
  pending_deletions <- reactiveVal(data.frame(
    stain_id = character(), F2_ppm = numeric(),
    F1_ppm = numeric(), Volume = numeric(),
    stringsAsFactors = FALSE
  ))
  
  # --- Click handling ---
  last_click_coords <- reactiveVal(NULL)
  first_click_for_box <- reactiveVal(NULL)
  
  # --- Other state ---
  calculated_contour_value <- reactiveVal(NULL)
  progress_bar <- reactiveVal(NULL)
  data_cc <- reactiveVal(NULL)
  plot_list <- reactiveVal(list())
  
  # --- Fit results ---
  fit_results_data <- reactiveVal(NULL)
  last_fit_method <- reactiveVal("sum")
  
  # --- Hide spinner when plot is updated ---
  # This observeEvent triggers when nmr_plot() changes and hides the spinner
  # after the plot is actually rendered (using session$onFlushed)
  observeEvent(nmr_plot(), {
    session$onFlushed(function() {
      shinyjs::runjs('document.getElementById("plot_spinner").style.display = "none";')
    }, once = TRUE)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #### Subfolder selector UI ----
  
  output$subfolder_selector <- renderUI({
    trigger_subfolder_update()
    spectra <- load_data$spectra_list()
    
    if (is.null(spectra) || length(spectra) == 0) {
      return(div(style = "color: #6c757d; font-style: italic; padding: 10px 0;",
                 "No spectra loaded yet. Load data first."))
    }
    
    subfolder_names <- names(spectra)
    if (is.null(subfolder_names) || length(subfolder_names) == 0) {
      return(div(style = "color: #dc3545;", "Error: No valid spectra names found."))
    }
    
    selectInput("selected_subfolder", "Chosen spectrum:", 
                choices = setNames(subfolder_names, basename(subfolder_names)))
  })
  
  #### Subfolder change handler ----
  
  observeEvent(input$selected_subfolder, {
    spectra <- load_data$spectra_list()
    req(spectra, spectra_plots())
    selected <- input$selected_subfolder
    if (is.null(selected) || !selected %in% names(spectra)) return()
    
    # Show spinner with message for spectrum change
    shinyjs::html("spinner_message", "Loading spectrum... This may take a moment for heavy files.")
    shinyjs::show("plot_spinner")
    
    # Update bruker_data in module with selected spectrum
    load_data$set_bruker_data(spectra[[selected]])
    
    # Update contour plot base with the NEW spectrum
    if (!is.null(spectra_plots()[[selected]])) {
      contour_plot_base(spectra_plots()[[selected]] + ggplot2::labs(title = ""))
      
      # Force rebuild the plot with NEW spectrum + EXISTING boxes/centroids
      refresh_nmr_plot(force_recalc = TRUE)
      
      status_msg(paste0("✅ Spectrum changed: ", basename(selected)))
    }
    
    # Hide spinner after plot is rendered
    session$onFlushed(function() {
      shinyjs::hide("plot_spinner")
    }, once = TRUE)
  }, ignoreInit = TRUE)
  
  
  ### 2.4 SPECTRUM PARAMETERS ----
  
  
  spectrum_params <- reactive({
    switch(input$spectrum_type,
           "TOCSY"  = list(intensity_threshold = 80000, contour_num = 40, contour_factor = 1.5, 
                           eps_value = 0.0068, neighborhood_size = 3),
           "HSQC"   = list(intensity_threshold = 20000, contour_num = 30, contour_factor = 1.3, 
                           eps_value = 0.068, neighborhood_size = 3),
           "COSY"   = list(intensity_threshold = 60000, contour_num = 30, contour_factor = 1.3, 
                           eps_value = 0.068, neighborhood_size = 9),
           "UFCOSY" = list(intensity_threshold = 50000, contour_num = 70, contour_factor = 1.3, 
                           eps_value = 0.014, neighborhood_size = 2),
           "HMBC"   = list(intensity_threshold = 20000, contour_num = 30, contour_factor = 1.3, 
                           eps_value = 0.009, neighborhood_size = 3),
           "JRES"   = list(intensity_threshold = 60000, contour_num = 40, contour_factor = 1.3, 
                           eps_value = 0.06, neighborhood_size = 3),
           # Default fallback: never return NULL (would break downstream params$...)
           list(intensity_threshold = 50000, contour_num = 30, contour_factor = 1.3,
                eps_value = 0.0068, neighborhood_size = 3)
    )
  })
  
  spectrum_params_CNN <- reactive({
    switch(input$spectrum_type,
           "TOCSY"  = list(int_thres = 0.01, int_prop = 0.001, eps_value = 0.0068,
                           pred_class_thres = 0.00001, batch_size = 64, step = 4),
           "UFCOSY" = list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.02,
                           pred_class_thres = 0.001, batch_size = 64, step = 4),
           "HSQC"   = list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.014,
                           pred_class_thres = 0.001, batch_size = 64, step = 4),
           "HMBC"   = list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.009,
                           pred_class_thres = 0.001, batch_size = 64, step = 4),
           "COSY"   = list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.068,
                           pred_class_thres = 0.001, batch_size = 64, step = 4),
           "JRES"   = list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.0068,
                           pred_class_thres = 0.001, batch_size = 64, step = 4),
           # Default fallback instead of stop(): a stop() raised inside a reactive
           # crashes the whole Shiny session (the app closes). Return safe defaults.
           list(int_thres = 0.001, int_prop = 0.5, eps_value = 0.0068,
                pred_class_thres = 0.001, batch_size = 64, step = 4)
    )
  })
  
  
  ### 2.5 PLOT GENERATION ----
  
  
  # Bounding boxes data (reactive with caching)
  
  bounding_boxes_data <- reactive({
    req(modifiable_boxes(), load_data$bruker_data())
    boxes <- modifiable_boxes()
    
    if (is.null(boxes) || nrow(boxes) == 0) {
      return(data.frame(xmin = numeric(0), xmax = numeric(0),
                        ymin = numeric(0), ymax = numeric(0),
                        stain_id = character(0), Volume = numeric(0)))
    }
    
    # Check cache
    cache_key <- paste0(input$selected_subfolder, "_", 
                        digest::digest(boxes[, c("xmin", "xmax", "ymin", "ymax")]))
    cached <- box_intensity_cache()
    if (!is.null(cached[[cache_key]])) return(cached[[cache_key]])
    
    # Calculate intensities
    mat <- load_data$bruker_data()$spectrumData
    if (is.null(mat)) return(boxes)
    
    ppm_x <- suppressWarnings(as.numeric(colnames(mat)))
    ppm_y <- suppressWarnings(as.numeric(rownames(mat)))
    
    if (!"stain_id" %in% names(boxes)) {
      boxes$stain_id <- paste0("box_", seq_len(nrow(boxes)))
    }
    
    boxes$Volume <- get_box_intensity(mat, ppm_x, ppm_y, boxes)
    # Clip negative volumes to 0 (same logic as mod_integration / mod_export):
    # keep the box on the plot, but don't display a negative Volume in the
    # data table or propagate it downstream.
    boxes$Volume <- pmax(boxes$Volume, 0, na.rm = TRUE)
    boxes$Volume[is.na(boxes$Volume)] <- 0
    
    # Update cache
    cached[[cache_key]] <- boxes
    box_intensity_cache(cached)
    boxes
  })
  
  #### 2.5.1 Refresh NMR plot ----
  
  refresh_nmr_plot <- function(force_recalc = FALSE) {
    req(contour_plot_base(), load_data$bruker_data())
    cache_key <- paste0(input$selected_subfolder,
                        "_boxes:", nrow(modifiable_boxes() %||% data.frame()),
                        "_centroids:", nrow(centroids_data() %||% data.frame()))
    cached_plots <- plot_cache()
    if (!force_recalc && !is.null(cached_plots[[cache_key]])) {
      # Tag with a fresh timestamp so nmr_plot() is ALWAYS a new value -- the
      # spinner-hiding observeEvent fires on value change, and two runs with the
      # same box/centroid COUNT would otherwise produce an identical value and
      # leave the spinner stuck / the new result not shown (HMBC/J-RES symptom).
      p <- cached_plots[[cache_key]]
      attr(p, "refresh_stamp") <- as.numeric(Sys.time())
      nmr_plot(p)
      return(invisible(NULL))
    }
    plot_base <- contour_plot_base()
    
    # Add bounding boxes (always red — conflict highlighting and tolerance
    # envelope are now drawn on top via plotlyProxy in an overlay observer
    # below, so they appear/disappear without rebuilding the whole plot).
    boxes <- tryCatch(bounding_boxes_data(), error = function(e) NULL)
    if (!is.null(boxes) && nrow(boxes) > 0) {
      bbox_path_df <- tryCatch(make_bbox_outline(boxes), error = function(e) NULL)
      if (!is.null(bbox_path_df) && nrow(bbox_path_df) > 0) {
        plot_base <- plot_base +
          geom_path(data = bbox_path_df, 
                    aes(x = x, y = y, group = group), 
                    color = "red", linewidth = 0.5, inherit.aes = FALSE)
      }
    }
    
    # Add centroids
    centrs <- imported_centroids() %||% centroids_data()
    if (!is.null(centrs) && nrow(centrs) > 0) {
      centrs_clean <- centrs
      centrs_clean$F2_ppm <- as.numeric(centrs_clean$F2_ppm)
      centrs_clean$F1_ppm <- as.numeric(centrs_clean$F1_ppm)
      if ("Volume" %in% names(centrs_clean)) {
        centrs_clean$intensity_plot <- as.numeric(centrs_clean$Volume)
        centrs_clean$intensity_plot[is.na(centrs_clean$intensity_plot)] <- 0
        centrs_clean$intensity_plot[is.infinite(centrs_clean$intensity_plot)] <- 0
        centrs_clean <- centrs_clean[!is.na(centrs_clean$F2_ppm) & !is.na(centrs_clean$F1_ppm), ]
        if (nrow(centrs_clean) > 0) {
          plot_base <- plot_base +
            geom_point(data = centrs_clean, 
                       aes(x = F2_ppm, y = F1_ppm, color = intensity_plot),
                       size = 1.2, inherit.aes = FALSE) +
            scale_color_gradient(low = "blue", high = "red", name = "Intensity")
        }
      } else {
        centrs_clean <- centrs_clean[!is.na(centrs_clean$F2_ppm) & !is.na(centrs_clean$F1_ppm), ]
        if (nrow(centrs_clean) > 0) {
          plot_base <- plot_base +
            geom_point(data = centrs_clean, aes(x = F2_ppm, y = F1_ppm),
                       color = "red", size = 1.2, inherit.aes = FALSE)
        }
      }
    }
    cached_plots[[cache_key]] <- plot_base
    plot_cache(cached_plots)
    attr(plot_base, "refresh_stamp") <- as.numeric(Sys.time())
    nmr_plot(plot_base)
  }
  
  #### 2.5.2 Generate plots button ----
  
  observeEvent(input$generate_plot, {
    
    # Show spinner with message during plot generation
    shinyjs::html("spinner_message", "Generating 2D spectrum... Please wait ☕")
    shinyjs::show("plot_spinner")
    
    # === LOADING NOTIFICATION ===
    loading_id <- showNotification(
      "☕ Generating 2D spectrum... Please wait",
      type = "message",
      duration = NULL,
      closeButton = FALSE
    )
    
    # Initial checks
    spectra <- load_data$spectra_list()
    if (is.null(spectra) || length(spectra) == 0) {
      shinyjs::hide("plot_spinner")
      removeNotification(loading_id)
      showNotification("⚠️ No spectra loaded. Please load data first.", type = "warning")
      return()
    }
    spectra_names <- names(spectra)
    if (is.null(spectra_names) || length(spectra_names) == 0) {
      showNotification("⚠️ Spectra have no names.", type = "warning")
      
      # Create default names
      spectra_names <- paste0("spectrum_", seq_along(spectra))
      names(spectra) <- spectra_names
      spectra_list(spectra)
    }
    status_msg("🔄 Generating plots...")
    params <- spectrum_params()
    n <- length(spectra)
    if (n == 0) {
      shinyjs::hide("plot_spinner")
      removeNotification(loading_id)
      showNotification("⚠️ No spectra to process.", type = "warning")
      return()
    }
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Processing spectra", value = 0)
    start_time <- Sys.time()
    all_results <- vector("list", n)  # Pre-allocate list with correct size
    for (i in seq_len(n)) {
      data <- spectra[[i]]
      spectrum_name <- spectra_names[i]
      
      # Check that data is valid
      if (is.null(data) || is.null(data$spectrumData)) {
        showNotification(paste("⚠️ Skipping invalid spectrum:", spectrum_name), type = "warning")
        all_results[[i]] <- NULL
        next
      }
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time <- if (i > 1) elapsed / (i - 1) else 0
      remaining <- round(avg_time * (n - i))
      time_msg <- if (i > 1) paste("⏱️ ~", remaining, "sec remaining") else ""
      progress$inc(1/n, detail = paste("Processing", basename(spectrum_name), time_msg))
      result <- tryCatch({
        find_nmr_peak_centroids_optimized(
          data$spectrumData,
          spectrum_type = input$spectrum_type,
          intensity_threshold = modulate_threshold(input$contour_start) %||% 
            modulate_threshold(calculated_contour_value()),
          contour_start = input$contour_start %||% calculated_contour_value(),
          contour_num = params$contour_num,
          contour_factor = params$contour_factor,
          f2_exclude_range = c(4.7, 5.0)
          
        )
      }, error = function(e) {
        showNotification(paste("❌ Error:", basename(spectrum_name), "-", e$message), type = "error")
        NULL
      })
      all_results[[i]] <- result
    }
    
    # Extract plots (keep NULL for failures)
    all_plots <- lapply(all_results, function(res) {
      if (!is.null(res) && !is.null(res$plot)) res$plot else NULL
    })
    
    # Assign names only if the lengths match
    if (length(all_plots) == length(spectra_names)) {
      names(all_plots) <- spectra_names
    } else {
      showNotification("⚠️ Mismatch in plots/names length", type = "warning")
    }
    spectra_plots(all_plots)
    
    # Find first valid result
    first_valid_idx <- which(sapply(all_results, function(x) !is.null(x) && !is.null(x$plot)))[1]
    if (!is.na(first_valid_idx) && length(all_results) > 0) {
      if (length(all_results) == length(spectra_names)) {
        names(all_results) <- spectra_names
      }
      result_data_list(all_results)
      result_data(all_results[[first_valid_idx]])
      
      # Check that plot exists before using it
      if (!is.null(all_results[[first_valid_idx]]$plot)) {
        contour_plot_base(all_results[[first_valid_idx]]$plot + labs(title = ""))
        refresh_nmr_plot(force_recalc = TRUE)  # Force rebuild with new spectrum
      }
      n_success <- sum(sapply(all_results, function(x) !is.null(x)))
      removeNotification(loading_id)
      showNotification(paste("✅", n_success, "/", n, "plots generated"), type = "message")
      status_msg(paste("✅", n_success, "plots generated"))
    } else {
      removeNotification(loading_id)
      showNotification("❌ No valid plots could be generated", type = "error")
      status_msg("❌ Plot generation failed")
    }
    # Hide spinner after plot is fully rendered
    session$onFlushed(function() {
      shinyjs::hide("plot_spinner")
    }, once = TRUE)
    shinyjs::hide("loading_message")
  })
  
  #### 2.5.3 Calculate contour threshold ----
  
  observeEvent(input$calculate_contour, {
    req(load_data$bruker_data())
    mat <- load_data$bruker_data()$spectrumData
    seuil <- switch(input$seuil_method,
                    "max_pct" = seuil_max_pourcentage(mat, pourcentage = input$pct_val),
                    "bruit_mult" = seuil_bruit_multiplicatif(mat, facteur = input$bruit_mult),
                    { showNotification("❌ Unknown method", type = "error"); return(NULL) }
                    
    )
    calculated_contour_value(seuil)
    showNotification(paste0("✅ Threshold: ", round(seuil, 2)), type = "message")
  })
  
  
  ### 2.6 MODULE INITIALIZATION ----
  
  # All major functionality is handled by Shiny modules.
  # Order matters due to dependencies between modules.
  
  # --- Edit state variables (shared with modules via rv) 
  selected_box_for_edit <- reactiveVal(NULL)
  selected_box_index <- reactiveVal(NULL)
  original_box_coords <- reactiveVal(NULL)
  box_has_been_modified <- reactiveVal(FALSE)
  preview_trace_added <- reactiveVal(FALSE)
  
  #### 2.6.1 Reactive values collection for modules ----
  rv <- list(
    # Data
    centroids_data = centroids_data,
    modifiable_boxes = modifiable_boxes,
    fixed_boxes = fixed_boxes,
    reference_boxes = reference_boxes,
    imported_centroids = imported_centroids,
    centroids = centroids,
    # Pending changes
    pending_centroids = pending_centroids,
    pending_boxes = pending_boxes,
    pending_fusions = pending_fusions,
    pending_deletions = pending_deletions,
    # Fit results
    fit_results_data = fit_results_data,
    last_fit_method = last_fit_method,
    # Plot state
    nmr_plot = nmr_plot,
    contour_plot_base = contour_plot_base,
    plot_cache = plot_cache,
    contour_cache = contour_cache,
    box_intensity_cache = box_intensity_cache,
    # Click state
    first_click_for_box = first_click_for_box,
    last_click_coords = last_click_coords,
    # Edit state
    selected_box_for_edit = selected_box_for_edit,
    selected_box_index = selected_box_index,
    original_box_coords = original_box_coords,
    box_has_been_modified = box_has_been_modified,
    preview_trace_added = preview_trace_added,
    # Contour threshold for integration (noise filter)
    contour_start = reactive({ input$contour_start })
  )
  
  
  #### 2.6.2 Module: INTEGRATION (initialized first - others depend on it) ----
  
  integration <- mod_integration_server(
    id = "integration",
    status_msg = status_msg,
    load_data = load_data,
    rv = rv
  )
  
  # --- Data reactives for modules ---
  data_reactives <- list(
    result_data_list = reactive({ result_data_list() }),
    bounding_boxes_data = reactive({ bounding_boxes_data() }),
    effective_integration_method = reactive({ integration$effective_integration_method() }),
    spectrum_params = reactive({ spectrum_params() }),
    calculated_contour_value = reactive({ calculated_contour_value() }),
    result_data = reactive({ result_data() })
  )
  
  
  #### 2.6.3 Module: SAVE & EXPORT ----
  
  save_export <- mod_save_export_server(
    id = "save_export",
    status_msg = status_msg,
    rv = rv,
    load_data = load_data,
    data_reactives = data_reactives,
    refresh_nmr_plot = refresh_nmr_plot,
    parent_session = session,
    parent_input = input
  )
  
  # Internal state: how many overlay traces are currently on the plot, and
  # whether we've injected annotations (arrows pointing at conflicting boxes).
  # Stored in a local environment (NOT a reactiveVal) because reading +
  # writing inside the same observe() would create a reactive loop.
  overlay_state <- new.env(parent = emptyenv())
  overlay_state$n_traces <- 0L
  overlay_state$has_annotations <- FALSE
  
  # Coordinate sign for plotly overlay. Inside ggplot the box coords are
  # stored as NEGATIVE values (NMR convention: -ppm for axis values, ticks
  # are re-labelled with sprintf("%.0f", -x_tickvals)). But when ggplotly()
  # converts to plotly, the displayed axis values match the (positive) ppm,
  # so plotlyProxy traces must be sent in the POSITIVE-ppm system — hence
  # we flip the sign of make_bbox_outline()'s output.
  OVERLAY_SIGN <- -1
  
  # Helper: turn a make_bbox_outline() data.frame (x, y, group) into a single
  # plotly trace where groups are separated by NA gaps. This way one trace
  # covers many rectangles, which is dramatically faster than one-trace-per-box
  # both for addTraces and for plotly's internal rendering.
  outline_to_trace <- function(path_df, color, width, dash = "solid", name = "overlay") {
    if (is.null(path_df) || nrow(path_df) == 0) return(NULL)
    groups <- split(path_df, path_df$group)
    xs <- unlist(lapply(groups, function(g) c(g$x * OVERLAY_SIGN, NA)), use.names = FALSE)
    ys <- unlist(lapply(groups, function(g) c(g$y * OVERLAY_SIGN, NA)), use.names = FALSE)
    list(
      x = xs, y = ys,
      type = "scatter", mode = "lines",
      line = list(color = color, width = width, dash = dash),
      name = name,
      hoverinfo = "skip",
      showlegend = FALSE,
      connectgaps = FALSE
    )
  }
  
  # Bridge: keep conflict_ids_rv AND shift_tol_rv in sync with the export
  # module's reactives, then update the tolerance/conflict overlay on the
  # main plot via plotlyProxy — NO full rebuild, no flicker.
  observe({
    new_ids <- tryCatch(save_export$conflict_ids() %||% character(0),
                        error = function(e) character(0))
    new_tol <- tryCatch(save_export$shift_tolerance_ppm() %||% 0,
                        error = function(e) 0)
    conflict_ids_rv(new_ids)
    shift_tol_rv(new_tol)
    
    proxy <- plotly::plotlyProxy("interactivePlot", session)
    
    # 1) Remove any previous overlay traces. Negative indices = "from the end"
    #    so we don't need to know how many other traces exist on the plot.
    n_prev <- overlay_state$n_traces
    if (!is.null(proxy) && n_prev > 0) {
      plotly::plotlyProxyInvoke(
        proxy, "deleteTraces",
        as.list(seq.int(-1, by = -1, length.out = n_prev))
      )
      overlay_state$n_traces <- 0L
    }
    # 1b) Also clear any previous arrow annotations (added below for conflicts).
    if (!is.null(proxy) && isTRUE(overlay_state$has_annotations)) {
      plotly::plotlyProxyInvoke(proxy, "relayout", list(annotations = list()))
      overlay_state$has_annotations <- FALSE
    }
    
    # 2) Early exit: nothing to overlay when tol=0 AND no conflicts. This is
    #    the default state — the original ggplot-drawn red boxes are enough
    #    and we don't want to paint anything on top. *** This is what makes
    #    the overlay disappear when the user moves the slider back to 0. ***
    if (new_tol <= 0 && length(new_ids) == 0) {
      return()
    }
    
    boxes <- tryCatch(modifiable_boxes(), error = function(e) NULL)
    if (is.null(boxes) || nrow(boxes) == 0 ||
        !all(c("xmin","xmax","ymin","ymax","stain_id") %in% names(boxes))) {
      return()
    }
    
    # Split boxes into OK and conflicting
    is_conf  <- boxes$stain_id %in% new_ids
    boxes_ok   <- boxes[!is_conf, , drop = FALSE]
    boxes_conf <- boxes[ is_conf, , drop = FALSE]
    
    traces_to_add <- list()
    
    # 3) Solid box outlines — only when something needs highlighting:
    #    - we have at least one conflict → repaint everything to show OK/conflict colors
    #    Otherwise (tol>0 but no conflicts), we leave the original red boxes
    #    untouched and just draw the dashed envelopes around them.
    if (length(new_ids) > 0) {
      if (nrow(boxes_ok) > 0) {
        path_ok <- tryCatch(make_bbox_outline(boxes_ok), error = function(e) NULL)
        tr <- outline_to_trace(path_ok, color = "#2e7d32", width = 1.5,
                               name = "boxes_ok_overlay")
        if (!is.null(tr)) traces_to_add[[length(traces_to_add) + 1]] <- tr
      }
      if (nrow(boxes_conf) > 0) {
        path_conf <- tryCatch(make_bbox_outline(boxes_conf), error = function(e) NULL)
        tr <- outline_to_trace(path_conf, color = "#c62828", width = 2.2,
                               name = "boxes_conflict_overlay")
        if (!is.null(tr)) traces_to_add[[length(traces_to_add) + 1]] <- tr
      }
    }
    
    # 4) Tolerance envelopes (dashed) — F2 expansion only.
    if (new_tol > 0) {
      expanded <- boxes
      expanded$xmin <- boxes$xmin - new_tol
      expanded$xmax <- boxes$xmax + new_tol
      # ymin/ymax untouched: tolerance applies on F2 (1H) only
      
      exp_ok   <- expanded[!is_conf, , drop = FALSE]
      exp_conf <- expanded[ is_conf, , drop = FALSE]
      
      if (nrow(exp_ok) > 0) {
        path_eok <- tryCatch(make_bbox_outline(exp_ok), error = function(e) NULL)
        env_color <- if (length(new_ids) > 0) "rgba(46,125,50,0.55)" else "rgba(230,81,0,0.6)"
        tr <- outline_to_trace(path_eok, color = env_color, width = 1, dash = "dash",
                               name = "env_ok_overlay")
        if (!is.null(tr)) traces_to_add[[length(traces_to_add) + 1]] <- tr
      }
      if (nrow(exp_conf) > 0) {
        path_econf <- tryCatch(make_bbox_outline(exp_conf), error = function(e) NULL)
        tr <- outline_to_trace(path_econf, color = "rgba(198,40,40,0.75)",
                               width = 1.2, dash = "dash",
                               name = "env_conflict_overlay")
        if (!is.null(tr)) traces_to_add[[length(traces_to_add) + 1]] <- tr
      }
    }
    
    # 5) Apply the new overlay
    if (length(traces_to_add) > 0 && !is.null(proxy)) {
      plotly::plotlyProxyInvoke(proxy, "addTraces", traces_to_add)
      overlay_state$n_traces <- length(traces_to_add)
    }
    
    # 6) Arrow annotations on conflicting boxes — one arrow per box, pointing
    #    at its center from a short offset. Pure plotly layout annotations
    #    (added via relayout), so they don't take up trace slots and can be
    #    cleared with a single relayout({annotations: []}).
    #
    #    Coordinate system: same as the traces — coords are flipped via
    #    OVERLAY_SIGN to match the ggplotly axis convention.
    #
    #    The arrow is placed up and to the left of the box (in display-ppm
    #    terms: higher-ppm side). Since the axis is reversed, "up-left" of
    #    the box maps to (xmax + small_offset, ymax + small_offset) BEFORE
    #    the OVERLAY_SIGN flip — keeping the label away from typical peak
    #    crowding which sits at lower ppm.
    if (length(new_ids) > 0 && nrow(boxes_conf) > 0 && !is.null(proxy)) {
      # Offsets in ppm — small enough not to leave the box's neighborhood
      # but large enough that the arrow tip + label don't sit on the box.
      off_x <- 0.05    # F2 offset (~0.05 ppm in 1H)
      off_y <- 1.2     # F1 offset (~1.2 ppm in 13C, since F1 spans ~150 ppm)
      
      annots <- lapply(seq_len(nrow(boxes_conf)), function(i) {
        r <- boxes_conf[i, ]
        # Box center (in storage coords, then flipped to display coords)
        cx <- ((r$xmin + r$xmax) / 2) * OVERLAY_SIGN
        cy <- ((r$ymin + r$ymax) / 2) * OVERLAY_SIGN
        # Label anchor — offset from center, on the high-ppm side
        # (which in display coords means cx + off_x, cy + off_y)
        list(
          x = cx, y = cy,                # arrow tip at box center
          ax = cx + off_x, ay = cy + off_y,  # label/tail position
          xref = "x", yref = "y",
          axref = "x", ayref = "y",
          text = paste0("<b>", r$stain_id, "</b>"),
          showarrow = TRUE,
          arrowhead = 3, arrowsize = 1.2, arrowwidth = 1.5,
          arrowcolor = "#c62828",
          font = list(color = "#c62828", size = 11, family = "Arial"),
          bgcolor = "rgba(255,255,255,0.85)",
          bordercolor = "#c62828",
          borderwidth = 1, borderpad = 2,
          opacity = 0.95
        )
      })
      plotly::plotlyProxyInvoke(proxy, "relayout", list(annotations = annots))
      overlay_state$has_annotations <- TRUE
    }
  })
  
  
  #### 2.6.4 Module: PEAK PICKING ----
  
  peak_picking <- mod_peak_picking_server(
    id = "peak_picking",
    status_msg = status_msg,
    load_data = load_data,
    data_reactives = data_reactives,
    rv = rv,
    refresh_nmr_plot = refresh_nmr_plot,
    parent_input = input
  )
  
  
  #### 2.6.5 Module: MANUAL EDITING ----
  
  manual_editing <- mod_manual_editing_server(
    id = "manual_editing",
    status_msg = status_msg,
    load_data = load_data,
    rv = rv,
    data_reactives = data_reactives,
    refresh_nmr_plot = refresh_nmr_plot,
    peak_picking = peak_picking,
    parent_input = input,
    parent_session = session
  )
  
  
  ### 2.7 REMAINING HANDLERS ----
  
  
  
  #### 2.7.1 Export projected centroids ----
  
  
  output$download_projected_centroids <- downloadHandler(
    filename = function() paste0("projected_centroids_", Sys.Date(), ".zip"),
    content = function(zipfile) {
      req(centroids_data(), result_data_list())
      tmp_dir <- tempdir()
      csv_files <- character(0)
      eps_val <- peak_picking$eps_value() %||% 0.04
      reference_centroids <- centroids_data()
      volumes_list <- list()
      eps_factors <- list(HSQC = 10, TOCSY = 8, COSY = 10, UFCOSY = 4, HMBC = 10, JRES = 8)
      for (name in names(result_data_list())) {
        result <- result_data_list()[[name]]
        if (is.null(result$contour_data)) next
        contour_data <- result$contour_data
        spec_type <- result$spectrum_type %||% "COSY"
        eps_factor <- eps_factors[[spec_type]] %||% 1
        delta_F2 <- 0; delta_F1 <- 0
        try({
          ref_hist <- with(reference_centroids, MASS::kde2d(F2_ppm, F1_ppm, n = 200))
          spec_hist <- with(contour_data, MASS::kde2d(-x, -y, n = 200))
          corr <- stats::convolve2d(ref_hist$z, spec_hist$z, type = "open")
          max_idx <- which(corr == max(corr, na.rm = TRUE), arr.ind = TRUE)
          delta_F2 <- (max_idx[1] - nrow(ref_hist$z)) * mean(diff(ref_hist$x))
          delta_F1 <- (max_idx[2] - ncol(ref_hist$z)) * mean(diff(ref_hist$y))
        }, silent = TRUE)
        shifted_centroids <- reference_centroids %>% 
          dplyr::mutate(F2_ppm = F2_ppm + delta_F2, F1_ppm = F1_ppm + delta_F1)
        projected_centroids <- shifted_centroids %>% 
          dplyr::rowwise() %>% 
          dplyr::mutate(Volume = {
            local_points <- contour_data %>% 
              dplyr::filter(sqrt((-x - F2_ppm)^2 + (-y - F1_ppm)^2) <= eps_val * eps_factor)
            sum(local_points$level, na.rm = TRUE)
          }) %>% 
          dplyr::ungroup()
        safe_name <- make.names(basename(name))
        
        output_csv <- file.path(tmp_dir, paste0(safe_name, "_projected_centroids.csv"))
        
        # Use write.csv2 for ";" separator (French Excel compatible)
        write.csv2(projected_centroids, output_csv, row.names = FALSE)
        csv_files <- c(csv_files, output_csv)
        vol_df <- projected_centroids %>% 
          dplyr::select(stain_id, F2_ppm, F1_ppm, Volume) %>% 
          dplyr::rename(!!paste0("Volume_", safe_name) := Volume)
        volumes_list[[safe_name]] <- vol_df
      }
      if (length(volumes_list) > 0) {
        merged_data <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("stain_id", "F2_ppm", "F1_ppm")), volumes_list)
        volume_cols <- grep("^Volume_", names(merged_data), value = TRUE)
        merged_data <- merged_data %>% dplyr::select(stain_id, F2_ppm, F1_ppm, dplyr::all_of(volume_cols))
        merged_data$status <- ifelse(rowSums(is.na(merged_data[volume_cols])) == 0, "ok", "missing")
        summary_csv <- file.path(tmp_dir, "summary_volumes.csv")
        
        # Use write.csv2 for ";" separator (French Excel compatible)
        write.csv2(merged_data, summary_csv, row.names = FALSE)
        csv_files <- c(csv_files, summary_csv)
      }
      zip(zipfile, files = csv_files, flags = "-j")
    }
    
  )
  
  #### 2.7.2 Save directory ----
  save_roots <- c(Home = normalizePath("~"), Root = "/")
  shinyDirChoose(input, "save_directory", roots = save_roots, session = session)
  save_directory <- reactive({ 
    req(input$save_directory)
    tryCatch({
      selection <- input$save_directory
      if (is.null(selection) || length(selection) == 0) {
        return(NULL)
      }
      selected_root <- selection$root
      if (is.null(selected_root) || !selected_root %in% names(save_roots)) {
        return(NULL)
      }
      base_path <- save_roots[[selected_root]]
      path_parts <- selection$path
      if (is.null(path_parts) || length(path_parts) == 0) {
        final_path <- base_path
      } else {
        path_parts <- unlist(path_parts)
        path_parts <- path_parts[path_parts != ""]
        if (length(path_parts) == 0) {
          final_path <- base_path
        } else {
          if (selected_root == "Root") {
            final_path <- paste0("/", paste(path_parts, collapse = "/"))
          } else {
            final_path <- file.path(base_path, paste(path_parts, collapse = "/"))
          }
        }
      }
      norm_path <- normalizePath(final_path, mustWork = FALSE)
      if (!dir.exists(norm_path)) return(NULL)
      norm_path
    }, error = function(e) {
      NULL
    })
  })
  
  output$save_dir_display <- renderPrint({ save_directory() })
  
  
  ### 2.8 UI OUTPUTS ----
  
  
  
  #### 2.8.1 Tool Description (Guide tab) ----
  
  output$toolDescription <- renderUI({
    tags$div(
      
      #'       
      
      #'       # Import de la font Orbitron
      
      #'       tags$style(HTML("
      
      #'   @import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@400;700;900&display=swap');
      
      #' ")),
      
      #'       
      
      #'       div(
      
      #'         style = "
      
      #'     background: rgba(8, 8, 16, 1);
      
      #'     border-radius: 20px;
      
      #'     padding: 25px 35px;
      
      #'     display: flex;
      
      #'     align-items: center;
      
      #'     gap: 30px;
      
      #'     margin-bottom: 20px;
      
      #'     box-shadow: 0 0 40px rgba(196,77,255,0.15), inset 0 0 60px rgba(196,77,255,0.05);
      
      #'     border: 1px solid rgba(196,77,255,0.15);
      
      #'   ",
      
      #'         
      
      #'         img(src = "spin.png", height = "130px", 
      
      #'             style = "filter: drop-shadow(0 0 25px rgba(255,107,107,0.4));"),
      
      #'         
      
      #'         div(
      
      #'           style = "flex: 1;",
      
      #'           h2(
      
      #'             style = "
      
      #'         font-family: 'Orbitron', sans-serif;
      
      #'         font-size: 36px;
      
      #'         font-weight: 900;
      
      #'         letter-spacing: 10px;
      
      #'         margin: 0;
      
      #'         background: linear-gradient(135deg, #ff6b6b 0%, #ffd93d 25%, #6bcb77 50%, #4d96ff 75%, #9b59b6 100%);
      
      #'         -webkit-background-clip: text;
      
      #'         -webkit-text-fill-color: transparent;
      
      #'       ",
      
      #'             "SPIN"
      
      #'           ),
      
      #'           p(
      
      #'             style = "color: rgba(255,255,255,0.6); font-size: 13px; margin-top: 8px;",
      
      #'             HTML("<span style='color: #c44dff; font-weight: 600;'>S</span>harp <span style='color: #c44dff; font-weight: 600;'>P</span>eak <span style='color: #c44dff; font-weight: 600;'>I</span>dentification for 2D <span style='color: #c44dff; font-weight: 600;'>N</span>MR")
      
      #'           )
      
      #'         )
      
      #'       ),
      
      # Quick Start
      div(style = "background: #e8f5e9; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #4caf50;",
          h4(style = "color: #2e7d32; margin-top: 0;", "🚀 Quick Start"),
          tags$ol(style = "margin-bottom: 0;",
                  tags$li(tags$b("Load"), " → Select your Bruker data folder"),
                  tags$li(tags$b("Plot"), " → Generate contour plots"),
                  tags$li(tags$b("Pick"), " → Detect peaks automatically"),
                  tags$li(tags$b("Edit"), " → Refine boxes manually if needed"),
                  tags$li(tags$b("Integrate"), " → Calculate volumes (Sum or Fitting)"),
                  tags$li(tags$b("Export"), " → Save results to CSV or session")
                  
          )
          
      ),
      
      # Tips
      h4("💡 Tips"),
      div(style = "background: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800; margin-bottom: 20px;",
          tags$ul(style = "margin-bottom: 0; font-size: 13px;",
                  tags$li("Start with a QC sample or the most intense spectrum to optimize parameters"),
                  tags$li("Use 'No clustering' option if you do not want to group multiplets"),
                  tags$li("Increase epsilon value to decrease size of cluster and get more solo peaks (e.g., TOCSY)"),
                  tags$li("Click 'Apply' to confirm changes before exporting"),
                  tags$li("Use the 'Data' tab to review and select boxes for editing"),
                  tags$li("If you want to process a batch, you might want to select the folder, only select a QC or the most intense spectrum. Process it (Step 1), then reload every spectra and use the 'batch export' with the peaks selected on the first spectrum."),
                  tags$li("For batch treatment, limit the number of spectrum per batch to 25 for TOCSY, 50 for COSY and HSQC.")
                  
          )
          
      ),
      
      # Main Features
      h4("📋 Features"),
      fluidRow(
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "📂 Data Loading"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li("Load Bruker 2D NMR data (ser/fid files)"),
                           tags$li("Batch processing of multiple spectra"),
                           tags$li("Select specific spectra to analyze"),
                           tags$li("Automatic detection of valid datasets")
                           
                   )
                   
               )
               
        ),
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "📈 Visualization"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li("Interactive Plotly contour plots (zoom, pan)"),
                           tags$li("Adaptive tick labels (positive values, NMR convention)"),
                           tags$li("Smooth zoom/dezoom with temporary axis hiding"),
                           tags$li("Real-time display of peaks and boxes")
                           
                   )
                   
               )
               
        )
        
      ),
      fluidRow(
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "🎯 Peak Detection"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li(tags$b("Local Max method:"), " Local maxima + DBSCAN clustering"),
                           tags$li(tags$b("CNN method:"), " Deep learning detection"),
                           tags$li("Automatic bounding box generation"),
                           tags$li("Configurable clustering parameters")
                           
                   )
                   
               )
               
        ),
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "✏️ Manual Editing"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li("4 collapsible sections (Click mode, Fusing, Edit, Add)"),
                           tags$li("Add/delete boxes directly on graph"),
                           tags$li("Move and resize with preview (green rectangle)"),
                           tags$li("Fuse multiple peaks into one, multi-selection (Ctrl+Click)")
                           
                   )
                   
               )
               
        )
        
      ),
      fluidRow(
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "📐 Integration & Peak Fitting"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li(tags$b("Direct:"), " Sum intensity"),
                           tags$li(tags$b("Fitting:"), " Gaussian, Voigt models"),
                           tags$li("Fit Quality tab with R² metrics"),
                           tags$li("2D fit visualization for each box")
                           
                   )
                   
               )
               
        ),
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px; margin-bottom: 15px; height: 180px;",
                   h5(style = "color: #1976d2;", "💾 Save & Export"),
                   tags$ul(style = "font-size: 13px;",
                           tags$li(tags$b("Session:"), " Save/Load complete work (.rds)"),
                           tags$li(tags$b("Import:"), " CSV files for peaks and boxes"),
                           tags$li(tags$b("Export:"), " CSV (French format ;), Batch export"),
                           tags$li("Pending system: Apply/Discard changes")
                           
                   )
                   
               )
               
        )
        
      ),
      
      # Supported Spectrum Types
      h4("🧪 Supported Spectrum Types"),
      div(style = "display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 20px;",
          div(style = "background: linear-gradient(135deg, #11998e, #38ef7d); color: white; padding: 10px 20px; border-radius: 20px;",
              tags$b("TOCSY"), " - Total Correlation Spectroscopy"
              
          ),
          div(style = "background: linear-gradient(135deg, #667eea, #764ba2); color: white; padding: 10px 20px; border-radius: 20px;",
              tags$b("HSQC"), " - Heteronuclear Single Quantum Coherence"
              
          ),
          div(style = "background: linear-gradient(135deg, #f093fb, #f5576c); color: white; padding: 10px 20px; border-radius: 20px;",
              tags$b("COSY"), " - Correlation Spectroscopy"
              
          ),
          div(style = "background: linear-gradient(135deg, #4facfe, #00f2fe); color: white; padding: 10px 20px; border-radius: 20px;",
              tags$b("UFCOSY"), " - Ultra-Fast COSY"
              
          )
          
      ),
      
      # Workflow Details
      h4("📖 Detailed Workflow"),
      
      # Step 1
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #11998e; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "1")),
            column(11,
                   h5(style = "margin-top: 5px;", "Load Data"),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     "Select a folder containing Bruker NMR data. The tool will automatically detect valid 2D spectra ",
                     "(folders containing 'acqus' and 'ser' or 'fid' files). You can select which spectra to load using checkboxes.")
                   
            )
            
          )
          
      ),
      
      # Step 2
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #667eea; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "2")),
            column(11,
                   h5(style = "margin-top: 5px;", "Generate Plot"),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     "Choose the spectrum type and adjust the intensity threshold. Click 'Auto' to calculate an optimal threshold ",
                     "based on noise level or maximum intensity. Then click 'Generate Plot' to create the contour visualization.")
                   
            )
            
          )
          
      ),
      
      # Step 3
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #f5576c; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "3")),
            column(11,
                   h5(style = "margin-top: 5px;", "Peak Picking"),
                   p(style = "font-size: 13px;", 
                     tags$b("Local Max method:"), " Uses local maxima detection followed by DBSCAN clustering to group nearby points. ",
                     "Adjust 'epsilon' to control cluster size."),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     tags$b("CNN method:"), " Uses a trained convolutional neural network for peak detection. ",
                     "Better for complex or overlapping peaks.")
                   
            )
            
          )
          
      ),
      
      # Step 4
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #4facfe; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "4")),
            column(11,
                   h5(style = "margin-top: 5px;", "Manual Editing"),
                   p(style = "font-size: 13px; margin-bottom: 8px;", 
                     tags$em("4 collapsible sections for a clean interface:")),
                   p(style = "font-size: 13px;", 
                     tags$b("🖱️ Click Mode:"), " Off / Add box (2 clicks) / Delete box on click"),
                   p(style = "font-size: 13px;", 
                     tags$b("🔗 Fusing Peaks & Boxes:"), " Multi-selection of peaks → merge into a single centroid"),
                   p(style = "font-size: 13px;", 
                     tags$b("📦 Edit Selected Box:"), " Modify coordinates, ↑↓←→ buttons, real-time preview (green)"),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     tags$b("➕ Add Manually:"), " Add a peak (F2, F1) or a box (xmin, xmax, ymin, ymax)")
                   
            )
            
          )
          
      ),
      
      # Step 5 - Integration
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #9b59b6; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "5")),
            column(11,
                   h5(style = "margin-top: 5px;", "Integration & Peak Fitting"),
                   p(style = "font-size: 13px;", 
                     tags$b("Direct Integration:"), "Sum (sum of intensities)"),
                   p(style = "font-size: 13px;", 
                     tags$b("Peak Fitting:"), " Gaussian ou Voigt (convolution Gauss-Lorentz)"),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     tags$b("Fit Quality:"), " Dedicated tab with R² metrics and 2D fit visualization")
                   
            )
            
          )
          
      ),
      
      # Step 6 - Export
      div(style = "background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 15px; margin-bottom: 10px;",
          fluidRow(
            column(1, div(style = "background: #fa709a; color: white; width: 30px; height: 30px; border-radius: 50%; text-align: center; line-height: 30px; font-weight: bold;", "6")),
            column(11,
                   h5(style = "margin-top: 5px;", "Save & Export"),
                   p(style = "font-size: 13px; margin-bottom: 8px;", 
                     tags$em("3 collapsible sections:")),
                   p(style = "font-size: 13px;", 
                     tags$b("💼 Session:"), " Complete Save/Load in .rds (peaks, boxes, parameters)"),
                   p(style = "font-size: 13px;", 
                     tags$b("📥 Import:"), " CSV files for peaks and boxes"),
                   p(style = "margin-bottom: 0; font-size: 13px;", 
                     tags$b("📤 Export:"), " CSV (semicolon separator), Batch Export (multiple spectra)")
                   
            )
            
          )
          
      ),
      
      # Output Format
      h4("📄 Output Format"),
      fluidRow(
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px;",
                   h5("Peaks CSV"),
                   tags$code(style = "font-size: 11px;",
                             "stain_id, F2_ppm, F1_ppm"
                             
                   )
                   
               )
               
        ),
        column(6,
               div(style = "background: #f5f5f5; padding: 15px; border-radius: 8px;",
                   h5("Boxes CSV"),
                   tags$code(style = "font-size: 11px;",
                             "stain_id, xmin, xmax, ymin, ymax, Volume"
                             
                   )
                   
               )
               
        )
        
      ),
      br(),
      
      # Recent Improvements
      h4("🆕 What's New (v2.0)"),
      div(style = "background: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3; margin-bottom: 20px;",
          fluidRow(
            column(6,
                   h5(style = "color: #1565c0; margin-top: 0;", "User Interface"),
                   tags$ul(style = "font-size: 13px; margin-bottom: 0;",
                           tags$li("Collapsible sections (Click mode, Fusing, Edit, Add)"),
                           tags$li("Colored accordion with icons"),
                           tags$li("Tables with multi-selection (Ctrl+Click)"),
                           tags$li("Save & Export in 3 sub-sections")
                           
                   )
                   
            ),
            column(6,
                   h5(style = "color: #1565c0; margin-top: 0;", "Axes & Data"),
                   tags$ul(style = "font-size: 13px; margin-bottom: 0;",
                           tags$li("Adaptive ticks (0.1, 0.5, 1, 2, 5...)"),
                           tags$li("Smooth zoom/unzoom (client-side JavaScript)"),
                           tags$li("Pending system (Apply/Discard)")
                           
                   )
                   
            )
            
          )
          
      ),
      
      # Footer
      div(style = "text-align: center; color: #9e9e9e; font-size: 12px; margin-top: 20px;",
          p("Developed for metabolomics research"),
          p("For questions or bug reports, please contact the development team")
          
      )
      
    )
  })
  
  #### 2.8.2 Status messages ----
  
  output$status_message <- renderText({ status_msg() })
  
  output$matrix_dim <- renderPrint({ req(load_data$bruker_data()); dim(load_data$bruker_data()$spectrumData) })
  
  output$seuil_text <- renderText({
    val <- calculated_contour_value()
    if (is.null(val)) return("No threshold calculated.")
    round(val, 5)
  })
  
  #### 2.8.3 Interactive plot ----
  
  output$interactivePlot <- renderPlotly({
    plot_obj <- nmr_plot()
    if (is.null(plot_obj)) {
      return(ggplotly(ggplot() + theme_void() + ggtitle("No spectrum displayed")))
    }
    tryCatch({
      plot_obj <- plot_obj + 
        theme(legend.text = element_text(size = 8), 
              legend.title = element_text(size = 9), 
              legend.key.size = unit(0.4, "cm"))
      
      # Ticks fixes simples : tous les 1 ppm pour F2
      
      # For F1: 1 ppm for TOCSY/COSY/UFCOSY (homonuclear 1H), 10 ppm for HSQC (13C)
      
      # Axis values are negative, invert for display (e.g. -3 becomes 3)
      
      # EXCEPT for true negative ppm values (e.g. -0.5 ppm stays -0.5)
      x_tickvals <- seq(-14, 2, by = 1)  # Covers -2 to 14 ppm in display
      x_ticktext <- sprintf("%.0f", -x_tickvals)  # Inverse le signe : -7 -> "7", 1 -> "-1"
      
      # Adapt F1 according to the spectrum type
      # HSQC and HMBC are heteronuclear: F1 = 13C (wide range, 10 ppm ticks)
      is_heteronuclear <- !is.null(input$spectrum_type) &&
        input$spectrum_type %in% c("HSQC", "HMBC")
      if (is_heteronuclear) {
        
        # HSQC/HMBC: F1 = 13C, ticks every 10 ppm (-10 to 230 ppm in display)
        y_tickvals <- seq(-240, 20, by = 10)
        y_ticktext <- sprintf("%.0f", -y_tickvals)
      } else {
        
        # TOCSY/COSY/UFCOSY: F1 = 1H, ticks every 1 ppm (-2 to 14 ppm in display)
        y_tickvals <- seq(-14, 2, by = 1)
        y_ticktext <- sprintf("%.0f", -y_tickvals)
      }
      p <- suppressWarnings({
        ggplotly(plot_obj, source = "nmr_plot") %>%
          layout(dragmode = "zoom",
                 xaxis = list(
                   showticklabels = TRUE, 
                   ticks = "outside",
                   tickmode = "array",
                   tickvals = x_tickvals,
                   ticktext = x_ticktext,
                   title = list(text = "F2 (ppm)", standoff = 10),
                   gridcolor = "rgba(200,200,200,0.3)"
                   
                 ), 
                 yaxis = list(
                   showticklabels = TRUE, 
                   ticks = "outside",
                   tickmode = "array",
                   tickvals = y_tickvals,
                   ticktext = y_ticktext,
                   title = list(text = "F1 (ppm)", standoff = 10),
                   gridcolor = "rgba(200,200,200,0.3)"
                   
                 )) %>%
          config(modeBarButtonsToAdd = list("select2d", "lasso2d"), displayModeBar = TRUE) %>%
          event_register("plotly_click") %>%
          event_register("plotly_selected") %>%
          event_register("plotly_relayout")
      })
      
      # Invisible grid to capture clicks
      
      # Denser grid means more precise clicks
      x_range <- layer_scales(plot_obj)$x$range$range
      y_range <- layer_scales(plot_obj)$y$range$range
      if (!is.null(x_range) && !is.null(y_range)) {
        
        # Very dense grid: 100x100 = 10000 points for maximum precision
        n_points <- 100
        grid <- expand.grid(x = seq(x_range[1], x_range[2], length.out = n_points),
                            y = seq(y_range[1], y_range[2], length.out = n_points))
        p <- p %>%
          add_trace(x = grid$x, y = grid$y, type = "scatter", mode = "markers",
                    marker = list(size = 6, opacity = 0),
                    hoverinfo = "none", showlegend = FALSE, name = "click_capture")
      }
      p
    }, error = function(e) {
      showNotification(paste("Plot error:", e$message), type = "error")
      ggplotly(ggplot() + theme_void() + ggtitle("Error"))
    })
  })
  
  #### 2.8.4 Click indicators ----
  
  output$clickedCoords <- renderPrint({
    coords <- last_click_coords()
    if (is.null(coords)) "Click on the spectrum" 
    else paste0("F2=", round(-coords$F2_ppm, 4), ", F1=", round(-coords$F1_ppm, 4))
  })
  
  #### 2.8.5 Tables ----
  
  output$centroid_table <- renderDT({
    df <- centroids_data()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame()))
    # Exclude Volume column from peaks display
    display_cols <- setdiff(names(df), c("Volume", "stain_intensity", "intensity_plot"))
    df_display <- df[, display_cols[display_cols %in% names(df)], drop = FALSE]
    datatable(df_display[, seq_len(min(4, ncol(df_display))), drop = FALSE], 
              selection = "multiple",
              rownames = FALSE,
              options = list(pageLength = 10))
  })
  
  output$full_centroid_table <- renderDT({
    df <- centroids_data() %||% data.frame()
    if (nrow(df) == 0) return(datatable(data.frame()))
    # Exclude Volume column from peaks display
    display_cols <- setdiff(names(df), c("Volume", "stain_intensity", "intensity_plot"))
    df_display <- df[, display_cols[display_cols %in% names(df)], drop = FALSE]
    datatable(df_display, selection = "multiple", rownames = FALSE, options = list(pageLength = 10))
  })
  
  output$bbox_table <- renderDT({
    df <- bounding_boxes_data() %||% data.frame()
    if (nrow(df) == 0) return(datatable(data.frame()))
    
    # Reorder columns to put stain_id first
    if ("stain_id" %in% names(df)) {
      other_cols <- setdiff(names(df), "stain_id")
      df <- df[, c("stain_id", other_cols), drop = FALSE]
    }
    
    datatable(df, 
              selection = "multiple",
              rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$pending_centroids_table <- renderDT({ 
    datatable(pending_centroids(), selection = "multiple", rownames = FALSE, options = list(pageLength = 10)) 
  })
  
  output$pending_boxes_table <- renderDT({ 
    df <- pending_boxes()
    if (nrow(df) == 0) return(datatable(data.frame()))
    
    # Reorder columns to put stain_id first if present
    if ("stain_id" %in% names(df)) {
      other_cols <- setdiff(names(df), "stain_id")
      df <- df[, c("stain_id", other_cols), drop = FALSE]
    }
    
    datatable(df, selection = "multiple", rownames = FALSE, options = list(pageLength = 10)) 
  })
  
  output$pending_fusions_table <- renderDT({ 
    req(pending_fusions())
    datatable(pending_fusions(), selection = "multiple", rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10)) 
  })
  
  output$pending_deletions_table <- renderDT({
    req(pending_deletions())
    datatable(pending_deletions(), selection = "multiple", rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 10))
  })
  
  #### 2.8.6 Spectrum type update ----
  
  observeEvent(input$spectrum_type, {
    params <- switch(input$spectrum_type,
                     "TOCSY" = list(contour_start = 80000),
                     "HSQC" = list(contour_start = 20000),
                     "COSY" = list(contour_start = 80000),
                     "UFCOSY" = list(contour_start = 30000),
                     "HMBC" = list(contour_start = 20000),
                     "JRES" = list(contour_start = 60000),
                     list(contour_start = 80000))  # default, avoid NULL
    req(params)
    
    cstart <- params$contour_start
    
    # HMBC: derive an initial contour_start from the actual noise of the loaded
    # spectrum (~7 sigma) instead of a fixed guess. This only PRE-FILLS the UI
    # field; the user can then change it and the change WILL be applied, because
    # find_nmr_peak_centroids_optimized respects any user-supplied value.
    if (identical(input$spectrum_type, "HMBC") && !is.null(load_data$bruker_data())) {
      mat <- load_data$bruker_data()$spectrumData
      if (!is.null(mat)) {
        noise_sd <- sd(as.numeric(mat), na.rm = TRUE)
        max_abs  <- max(abs(mat), na.rm = TRUE)
        if (is.finite(noise_sd) && noise_sd > 0) {
          cstart <- round(min(noise_sd * 7, max_abs * 0.25))
        }
      }
    }
    
    updateNumericInput(session, "contour_start", value = cstart)
    
    # Note: eps_value is now updated by the peak_picking module
  })
  
  #### 2.8.7 Dragmode update ----
  
  observeEvent(input$plotly_dragmode, {
    plotlyProxy("interactivePlot", session) %>%
      plotlyProxyInvoke("relayout", list(dragmode = input$plotly_dragmode))
  })
  
  #### 2.8.8 Dynamic axis ticks on zoom ----
  
  # Store the complete plot ranges for autoscale
  full_plot_ranges <- reactiveVal(list(x = c(-10, 0), y = c(-200, 0)))
  
  # Observer to capture initial ranges when plot is created
  observe({
    plot_obj <- nmr_plot()
    if (!is.null(plot_obj)) {
      x_range <- tryCatch({
        layer_scales(plot_obj)$x$range$range
      }, error = function(e) c(-10, 0))
      y_range <- tryCatch({
        layer_scales(plot_obj)$y$range$range
      }, error = function(e) c(-200, 0))
      if (!is.null(x_range) && !is.null(y_range)) {
        full_plot_ranges(list(x = x_range, y = y_range))
      }
    }
  })
  
  # Function to generate ticks with correct displayed ppm values
  
  # Axis data is opposite of real ppm values (e.g. ppm=3 -> axis=-3)
  generate_positive_ticks <- function(range_vals, target_nticks = 10, decimals = 2) {
    if (any(is.null(range_vals)) || any(is.na(range_vals)) || length(range_vals) < 2) return(NULL)
    
    # Axis range (negative ppm values, e.g. -7.5 to 0.5 if ppm goes from 7.5 to -0.5)
    axis_min <- min(range_vals)
    axis_max <- max(range_vals)
    
    # Valeurs ppm correspondantes (on inverse le signe)
    ppm_max <- -axis_min  # ex: 7.5
    ppm_min <- -axis_max  # ex: -0.5
    span <- ppm_max - ppm_min
    if (span <= 0 || !is.finite(span)) return(NULL)
    
    # Calculate a "nice" step
    rough_step <- span / target_nticks
    if (rough_step <= 0 || !is.finite(rough_step)) return(NULL)
    magnitude <- 10^floor(log10(rough_step))
    nice_steps <- c(0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10)
    step <- magnitude * nice_steps[which.min(abs(nice_steps - rough_step/magnitude))]
    if (step <= 0 || !is.finite(step)) return(NULL)
    
    # Generate tick values in ppm
    tick_start <- floor(ppm_min / step) * step
    tick_end <- ceiling(ppm_max / step) * step
    tick_values_ppm <- seq(tick_start, tick_end, by = step)
    
    # Filter to keep only the ticks within the visible range (with a small margin)
    margin <- step * 0.1
    tick_values_ppm <- tick_values_ppm[
      tick_values_ppm >= (ppm_min - margin) & 
        tick_values_ppm <= (ppm_max + margin)
    ]
    if (length(tick_values_ppm) == 0) return(NULL)
    
    # Axis values are opposite of ppm
    tick_values_axis <- -tick_values_ppm
    
    # Format according to decimal places
    fmt <- paste0("%.", decimals, "f")
    list(
      tickvals = tick_values_axis,
      ticktext = sprintf(fmt, tick_values_ppm)  # Affiche les vraies valeurs ppm
      
    )
  }
  
  # Variable to prevent relayout loops (kept for compatibility)
  last_tick_update <- reactiveVal(0)
  
  # Note: Tick update during zoom is handled client-side (JavaScript)
  
  # for better smoothness. This observer is kept as backup for cases
  
  # where JavaScript doesn't trigger correctly.
  
  observeEvent(event_data("plotly_relayout", source = "nmr_plot"), {
    
    # Ticks are now managed by client-side JavaScript
    
    # for synchronous update with zoom
    NULL
  }, ignoreNULL = FALSE)
  
  #### 2.8.9 Pending changes indicators ----
  
  output$has_pending_changes <- reactive({
    n_centroids <- nrow(pending_centroids() %||% data.frame())
    n_boxes <- nrow(pending_boxes() %||% data.frame())
    n_fusions <- nrow(pending_fusions() %||% data.frame())
    n_deletions <- nrow(pending_deletions() %||% data.frame())    # <-- AJOUTER
    (n_centroids + n_boxes + n_fusions + n_deletions) > 0         # <-- MODIFIER
  })
  
  outputOptions(output, "has_pending_changes", suspendWhenHidden = FALSE)
  
  output$pending_summary <- renderText({
    n_centroids <- nrow(pending_centroids() %||% data.frame())
    n_boxes <- nrow(pending_boxes() %||% data.frame())
    n_fusions <- nrow(pending_fusions() %||% data.frame())
    n_deletions <- nrow(pending_deletions() %||% data.frame())    # <-- AJOUTER
    parts <- c()
    if (n_centroids > 0) parts <- c(parts, paste(n_centroids, "peaks"))
    if (n_boxes > 0) parts <- c(parts, paste(n_boxes, "boxes"))
    if (n_fusions > 0) parts <- c(parts, paste(n_fusions, "fusions"))
    if (n_deletions > 0) parts <- c(parts, paste(n_deletions, "deletions"))  # <-- AJOUTER
    if (length(parts) == 0) return("")
    paste(parts, collapse = ", ")
  })
  
  observeEvent(load_data$spectra_list(), { centroids(NULL) })
  
  #### 2.8.10 Fit Quality Visualizations ----
  
  # Combined data: boxes + fit results
  boxes_with_fit <- reactive({
    boxes <- modifiable_boxes()
    fit_data <- fit_results_data()
    if (is.null(boxes) || nrow(boxes) == 0) return(NULL)
    if (is.null(fit_data)) return(boxes)  # Retourner boxes sans fit info
    
    # Join fit data with boxes
    boxes_merged <- boxes %>%
      left_join(fit_data, by = "stain_id")
    boxes_merged
  })
  
  # Summary table of fitting results
  fit_summary_data <- reactive({
    boxes <- boxes_with_fit()
    if (is.null(boxes) || !"fit_method" %in% names(boxes)) {
      return(NULL)
    }
    
    # Summary by method
    summary_df <- boxes %>%
      filter(!is.na(fit_method)) %>%
      group_by(fit_method) %>%
      summarise(
        n_boxes = n(),
        mean_r2 = mean(r_squared, na.rm = TRUE),
        median_r2 = median(r_squared, na.rm = TRUE),
        min_r2 = min(r_squared, na.rm = TRUE),
        max_r2 = max(r_squared, na.rm = TRUE),
        .groups = "drop"
        
      )
    summary_df
  })
  
  # R² distribution plot
  
  output$fit_quality_plot <- renderPlotly({
    boxes <- boxes_with_fit()
    
    # Check if we have fit data
    if (is.null(boxes) || !"r_squared" %in% names(boxes) || all(is.na(boxes$r_squared))) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No fit quality data available.\n\nSteps to generate fit data:\n1. Select 'Gaussian fit' or 'Voigt fit' in Export section\n2. Click 'Batch Export' button\n3. Return to this tab to see results",
                 size = 4, color = "gray50", hjust = 0.5, vjust = 0.5) +
        theme_void() +
        ggtitle("Fit Quality Distribution")
      return(ggplotly(p))
    }
    
    # Filter NAs
    boxes_with_r2 <- boxes %>% filter(!is.na(r_squared))
    if (nrow(boxes_with_r2) == 0) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No R² values computed yet.\nRun 'Batch Export' with fitting method.",
                 size = 5, color = "gray50") +
        theme_void()
      return(ggplotly(p))
    }
    
    # R² histogram with colors by method
    p <- ggplot(boxes_with_r2, aes(x = r_squared, fill = fit_method)) +
      geom_histogram(bins = 30, color = "white", alpha = 0.8, position = "stack") +
      geom_vline(aes(xintercept = median(r_squared, na.rm = TRUE)), 
                 color = "red", linetype = "dashed", size = 1) +
      scale_fill_manual(
        values = c("gaussian" = "#667eea", "voigt" = "#f5576c", "multiplet_fit" = "#38ef7d", "multiplet_sum" = "#ff9800", "sum_fit_failed" = "#ffd93d", "sum_r2_below" = "#ffcc00"),
        name = "Fit Method"
        
      ) +
      labs(
        title = "Distribution of Fit Quality (R²)",
        subtitle = paste0("Median R² = ", round(median(boxes_with_r2$r_squared, na.rm = TRUE), 3),
                          " | Total boxes = ", nrow(boxes_with_r2)),
        x = "R² (Coefficient of Determination)",
        y = "Number of Peaks"
        
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray30")
        
      )
    ggplotly(p) %>%
      layout(hovermode = "x unified")
  })
  
  # Reactive to get selected box (from fit_boxes_detail_table)
  selected_fit_box <- reactive({
    boxes <- boxes_with_fit()
    if (is.null(boxes) || !"fit_method" %in% names(boxes)) return(NULL)
    
    # Get selection in fit_boxes_detail_table
    selected_row <- input$fit_boxes_detail_table_rows_selected
    if (is.null(selected_row) || length(selected_row) == 0) return(NULL)
    
    # Rebuild same order as in table (sorted by R² descending)
    detail_df <- boxes %>%
      filter(!is.na(fit_method)) %>%
      arrange(desc(r_squared))
    if (selected_row > nrow(detail_df)) return(NULL)
    
    # Get selected stain_id
    selected_stain_id <- detail_df$stain_id[selected_row]
    
    # Return corresponding box from boxes_with_fit (with all columns)
    boxes %>% filter(stain_id == selected_stain_id)
  })
  
  # 2D fit example for a selected box - CONTOUR OVERLAY (TopSpin-style)
  
  output$example_fit_2d <- renderPlot({
    req(load_data$bruker_data())
    box <- selected_fit_box()
    if (is.null(box) || nrow(box) == 0) {
      par(bg = "white", mar = c(1, 1, 1, 1))
      plot.new()
      text(0.5, 0.5, "Select a box in the\n'Fitted Boxes Details' table\nto visualize its fit", cex = 1.4, col = "gray40", font = 2)
      return()
    }
    box <- box[1, ]
    
    # Check that fit exists
    if (!"fit_method" %in% names(box) || is.na(box$fit_method) || grepl("sum", box$fit_method)) {
      par(bg = "white", mar = c(1, 1, 1, 1))
      plot.new()
      text(0.5, 0.5, paste0("'", box$stain_id, "'\n\nNo fit available\n(sum method used)"), 
           cex = 1.2, col = "#e65100", font = 2)
      return()
    }
    
    # Get spectral data
    mat <- load_data$bruker_data()$spectrumData
    ppm_x <- suppressWarnings(as.numeric(colnames(mat)))
    ppm_y <- suppressWarnings(as.numeric(rownames(mat)))
    
    # Get contour threshold - same as main plot
    # NOTE: We do NOT filter the data here (no setting values to 0)
    # The main plot in Vizualisation.R uses raw data with contour breaks starting at contour_start
    # This gives the same visual result without modifying the underlying data
    contour_threshold <- input$contour_start
    if (is.null(contour_threshold) || contour_threshold <= 0) {
      contour_threshold <- max(mat, na.rm = TRUE) * 0.1  # Fallback
    }
    
    # Get zoom level from slider (default 80%)
    # 0% = very zoomed out (lots of context around the box)
    # 100% = tight on the box (no padding)
    zoom_pct <- input$fit_zoom_level
    if (is.null(zoom_pct)) zoom_pct <- 80
    
    # Calculate padding based on zoom level
    # At 0%: padding = 100% of box size (see 2x the box area)
    # At 50%: padding = 50% of box size
    # At 100%: padding = 0% (exactly the box, no extra)
    box_width <- box$xmax - box$xmin
    box_height <- box$ymax - box$ymin
    
    # Linear interpolation: zoom 0% -> 1.0 padding, zoom 100% -> 0.0 padding
    padding_factor <- 1.0 - (zoom_pct / 100)
    padding_x_box <- box_width * padding_factor
    padding_y_box <- box_height * padding_factor
    
    # FIX: Guarantee a MINIMUM absolute padding for small boxes 
    # For very small boxes (e.g. narrow UFCOSY peaks), padding proportional to box size
    # can be smaller than a single ppm step -> no extra context shown when zooming out.
    # We compute a minimum padding based on the spectrum's actual ppm resolution and
    # take the max with the proportional padding. The minimum scales with padding_factor
    # so that zoom 100% (tight on box) is still respected.
    
    # Estimate ppm step (resolution) on each axis
    # Use median of absolute diffs to be robust to non-uniform spacing
    ppm_x_step <- if (length(ppm_x) > 1) median(abs(diff(ppm_x))) else 0.001
    ppm_y_step <- if (length(ppm_y) > 1) median(abs(diff(ppm_y))) else 0.001
    
    # Minimum padding: ~20 points of resolution, capped between 0.02 and 0.15 ppm
    # This guarantees enough visual context regardless of box size
    min_pad_x <- max(min(20 * ppm_x_step, 0.15), 0.02)
    min_pad_y <- max(min(20 * ppm_y_step, 0.15), 0.02)
    
    # Scale minimum by padding_factor so zoom=100% still gives 0 extra padding
    padding_x <- max(padding_x_box, min_pad_x * padding_factor)
    padding_y <- max(padding_y_box, min_pad_y * padding_factor)
    
    
    x_idx <- which(ppm_x >= (box$xmin - padding_x) & ppm_x <= (box$xmax + padding_x))
    y_idx <- which(ppm_y >= (box$ymin - padding_y) & ppm_y <= (box$ymax + padding_y))
    
    # Fallback to exact box if padding gives empty result
    if (length(x_idx) == 0) x_idx <- which(ppm_x >= box$xmin & ppm_x <= box$xmax)
    if (length(y_idx) == 0) y_idx <- which(ppm_y >= box$ymin & ppm_y <= box$ymax)
    
    if (length(x_idx) == 0 || length(y_idx) == 0) {
      plot.new()
      text(0.5, 0.5, "Region out of bounds", cex = 1.5, col = "red")
      return()
    }
    
    # Use RAW data for display (same as Vizualisation.R)
    # The contour levels starting at contour_start will filter out noise
    region <- mat[y_idx, x_idx, drop = FALSE]
    x_sub <- ppm_x[x_idx]
    y_sub <- ppm_y[y_idx]
    
    # Determine which model to use based on fit_method
    use_voigt <- grepl("voigt", box$fit_method, ignore.case = TRUE)
    model_name <- if (use_voigt) "Voigt" else "Gaussian"
    
    # FIT FOR VISUALIZATION (MULTIPLET SUPPORT, GAUSSIAN OR VOIGT)
    # Detect peaks and fit each one separately, then combine for display
    # IMPORTANT: We fit on the BOX region (not the padded view region)
    # but we generate the model surface over the FULL displayed region
    # This ensures the model is always visible regardless of zoom level
    
    # Extract the box region for fitting (without padding) - USE RAW DATA
    box_x_idx <- which(ppm_x >= box$xmin & ppm_x <= box$xmax)
    box_y_idx <- which(ppm_y >= box$ymin & ppm_y <= box$ymax)
    box_region <- mat[box_y_idx, box_x_idx, drop = FALSE]
    box_x_sub <- ppm_x[box_x_idx]
    box_y_sub <- ppm_y[box_y_idx]
    
    fit_result <- tryCatch({
      
      # Detect local maxima on the BOX region (not the padded view)
      local_max <- detect_local_maxima(box_region, threshold = 0.5, min_distance = 3)
      n_peaks <- nrow(local_max)
      
      # Fallback if no peaks detected
      if (n_peaks == 0) {
        max_val <- max(box_region, na.rm = TRUE)
        max_pos <- which(box_region == max_val, arr.ind = TRUE)
        if (nrow(max_pos) > 0) {
          local_max <- data.frame(row = max_pos[1,1], col = max_pos[1,2], value = max_val)
          n_peaks <- 1
        } else {
          stop("No peaks found")
        }
      }
      
      # Initialize combined fitted matrix for the DISPLAYED region (with padding)
      combined_fitted <- matrix(0, nrow = length(y_sub), ncol = length(x_sub))
      
      # Store individual peak fits for potential separate display
      peak_fits <- list()
      
      # Global baseline estimate from box region
      baseline_global <- quantile(as.vector(box_region), 0.1, na.rm = TRUE)
      
      # Fit each peak separately (using box coordinates)
      for (p in seq_len(n_peaks)) {
        peak_row <- local_max$row[p]
        peak_col <- local_max$col[p]
        # Convert box indices to ppm coordinates
        peak_x <- box_x_sub[peak_col]
        peak_y <- box_y_sub[peak_row]
        peak_amplitude <- local_max$value[p]
        
        # Define sub-region around this peak (in box coordinates)
        if (n_peaks > 1) {
          distances <- sqrt((local_max$col - peak_col)^2 + (local_max$row - peak_row)^2)
          distances[p] <- Inf
          min_dist <- min(distances)
          half_width <- max(2, floor(min_dist / 2))
        } else {
          half_width <- max(2, floor(min(nrow(box_region), ncol(box_region)) / 2))
        }
        
        col_start <- max(1, peak_col - half_width)
        col_end <- min(ncol(box_region), peak_col + half_width)
        row_start <- max(1, peak_row - half_width)
        row_end <- min(nrow(box_region), peak_row + half_width)
        
        sub_region <- box_region[row_start:row_end, col_start:col_end, drop = FALSE]
        sub_x <- box_x_sub[col_start:col_end]
        sub_y <- box_y_sub[row_start:row_end]
        
        # Create grid for this sub-region
        sub_grid <- expand.grid(x = sub_x, y = sub_y)
        sub_grid$z <- as.vector(t(sub_region))
        sub_grid <- sub_grid[!is.na(sub_grid$z), ]
        
        if (nrow(sub_grid) < 6) next
        
        # Normalize
        z_scale <- max(abs(sub_grid$z), na.rm = TRUE)
        if (z_scale < 1e-10) z_scale <- 1
        sub_grid$z_norm <- sub_grid$z / z_scale
        
        # Initial parameters for this peak
        sigma_x_init <- diff(range(sub_x)) / 3
        sigma_y_init <- diff(range(sub_y)) / 3
        baseline_init <- quantile(sub_grid$z, 0.1, na.rm = TRUE)
        
        # Fit this peak (Gaussian or Voigt depending on method)
        # Use suppressWarnings to silence "nls.lm: info = -1" (max iterations reached);
        # nlsLM still returns a usable result, the warning is purely informational.
        peak_fit <- tryCatch({
          if (use_voigt) {
            # Pseudo-Voigt fit
            fit <- suppressWarnings(minpack.lm::nlsLM(
              z_norm ~ A * (eta / (1 + ((x - x0) / gx)^2 + ((y - y0) / gy)^2) + 
                              (1 - eta) * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2)))) + b,
              data = sub_grid,
              start = list(A = peak_amplitude / z_scale, x0 = peak_x, y0 = peak_y,
                           sx = sigma_x_init, sy = sigma_y_init,
                           gx = sigma_x_init, gy = sigma_y_init,
                           eta = 0.5, b = baseline_init / z_scale),
              lower = c(A = 0, x0 = min(sub_x), y0 = min(sub_y),
                        sx = sigma_x_init / 10, sy = sigma_y_init / 10,
                        gx = sigma_x_init / 10, gy = sigma_y_init / 10,
                        eta = 0, b = -Inf),
              upper = c(A = Inf, x0 = max(sub_x), y0 = max(sub_y),
                        sx = sigma_x_init * 5, sy = sigma_y_init * 5,
                        gx = sigma_x_init * 5, gy = sigma_y_init * 5,
                        eta = 1, b = Inf),
              control = list(maxiter = 100, gtol = 0)
            ))
            params <- coef(fit)
            list(params = params, model = "voigt")
          } else {
            # Gaussian fit
            fit <- suppressWarnings(minpack.lm::nlsLM(
              z_norm ~ A * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2))) + b,
              data = sub_grid,
              start = list(A = peak_amplitude / z_scale, x0 = peak_x, y0 = peak_y,
                           sx = sigma_x_init, sy = sigma_y_init, b = baseline_init / z_scale),
              lower = c(A = 0, x0 = min(sub_x), y0 = min(sub_y),
                        sx = sigma_x_init / 10, sy = sigma_y_init / 10, b = -Inf),
              upper = c(A = Inf, x0 = max(sub_x), y0 = max(sub_y),
                        sx = sigma_x_init * 5, sy = sigma_y_init * 5, b = Inf),
              control = list(maxiter = 100, gtol = 0)
            ))
            params <- coef(fit)
            list(params = params, model = "gaussian")
          }
        }, error = function(e) {
          # Fallback: use rough Gaussian estimates
          list(
            params = c(A = peak_amplitude / z_scale, x0 = peak_x, y0 = peak_y,
                       sx = sigma_x_init, sy = sigma_y_init, b = baseline_init / z_scale),
            model = "gaussian_fallback"
          )
        })
        
        # Store this peak's parameters
        peak_fits[[p]] <- list(
          params = peak_fit$params,
          model = peak_fit$model,
          z_scale = z_scale,
          center = c(peak_x, peak_y)
        )
        
        # Generate fitted surface for this peak over the ENTIRE region
        params <- peak_fit$params
        for (i in seq_along(y_sub)) {
          for (j in seq_along(x_sub)) {
            x_val <- x_sub[j]
            y_val <- y_sub[i]
            
            if (peak_fit$model == "voigt") {
              # Pseudo-Voigt: eta * Lorentzian + (1-eta) * Gaussian
              eta <- params["eta"]
              lorentz <- 1 / (1 + ((x_val - params["x0"]) / params["gx"])^2 + 
                                ((y_val - params["y0"]) / params["gy"])^2)
              gauss <- exp(-((x_val - params["x0"])^2 / (2 * params["sx"]^2) + 
                               (y_val - params["y0"])^2 / (2 * params["sy"]^2)))
              peak_val <- params["A"] * z_scale * (eta * lorentz + (1 - eta) * gauss)
            } else {
              # Gaussian
              peak_val <- params["A"] * z_scale * exp(
                -((x_val - params["x0"])^2 / (2 * params["sx"]^2) + 
                    (y_val - params["y0"])^2 / (2 * params["sy"]^2))
              )
            }
            
            combined_fitted[i, j] <- combined_fitted[i, j] + peak_val
          }
        }
      }
      
      # Add baseline (only once, not per peak)
      combined_fitted <- combined_fitted + baseline_global
      
      # FIX: Use FITTED centers (x0, y0) from each peak's nlsLM result
      # Previously used local_max (raw maxima detected BEFORE fitting), which made the
      # displayed centers incoherent with the red model contours. Now we extract the
      # actual fitted centers from peak_fits[[p]]$params["x0"/"y0"] so that the green
      # crosses fall exactly on the centers of the red model lobes.
      fitted_centers_x <- numeric(0)
      fitted_centers_y <- numeric(0)
      for (p in seq_along(peak_fits)) {
        if (!is.null(peak_fits[[p]]) && !is.null(peak_fits[[p]]$params)) {
          pp <- peak_fits[[p]]$params
          if (!is.na(pp["x0"]) && !is.na(pp["y0"])) {
            fitted_centers_x <- c(fitted_centers_x, unname(pp["x0"]))
            fitted_centers_y <- c(fitted_centers_y, unname(pp["y0"]))
          }
        }
      }
      
      # Fallback: if no fitted centers (all fits failed), use local maxima as proxy
      if (length(fitted_centers_x) == 0) {
        fitted_centers_x <- box_x_sub[local_max$col]
        fitted_centers_y <- box_y_sub[local_max$row]
      }
      
      list(success = TRUE, 
           fitted_matrix = combined_fitted, 
           n_peaks = n_peaks,
           peak_fits = peak_fits,
           model_used = model_name,
           peak_centers = data.frame(
             x = fitted_centers_x,
             y = fitted_centers_y
           ))
      
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })
    
    # Ensure ascending order for contour()
    x_reorder <- NULL
    y_reorder <- NULL
    
    if (is.unsorted(x_sub)) {
      x_reorder <- order(x_sub)
      x_sub <- x_sub[x_reorder]
      region <- region[, x_reorder, drop = FALSE]
      # FIX: Also reorder the fitted matrix so red contours stay aligned
      # Without this, combined_fitted was built using the ORIGINAL (unsorted) x_sub indices,
      # then plotted against the sorted x_sub -> the model appeared mirrored on the F2 axis.
      if (fit_result$success && !is.null(fit_result$fitted_matrix)) {
        fit_result$fitted_matrix <- fit_result$fitted_matrix[, x_reorder, drop = FALSE]
      }
    }
    if (is.unsorted(y_sub)) {
      y_reorder <- order(y_sub)
      y_sub <- y_sub[y_reorder]
      region <- region[y_reorder, , drop = FALSE]
      # FIX: Same fix on F1 axis 
      if (fit_result$success && !is.null(fit_result$fitted_matrix)) {
        fit_result$fitted_matrix <- fit_result$fitted_matrix[y_reorder, , drop = FALSE]
      }
    }
    
    # TOPSPIN-STYLE VISUALIZATION
    
    # Get contour parameters from main plot settings (same as base plot)
    contour_start <- input$contour_start
    if (is.null(contour_start) || contour_start <= 0) {
      # Fallback: use reasonable default based on spectrum type
      contour_start <- max(region, na.rm = TRUE) * 0.1
    }
    
    # Get contour_factor and contour_num based on spectrum type (same as Vizualisation.R)
    spectrum_type <- input$spectrum_type
    contour_params <- switch(spectrum_type,
                             "HSQC" = list(contour_num = 8, contour_factor = 1.3),
                             "TOCSY" = list(contour_num = 20, contour_factor = 1.3),
                             "COSY" = list(contour_num = 30, contour_factor = 1.3),
                             "UFCOSY" = list(contour_num = 20, contour_factor = 1.3),
                             "HMBC" = list(contour_num = 18, contour_factor = 1.25),
                             "JRES" = list(contour_num = 20, contour_factor = 1.3),
                             list(contour_num = 10, contour_factor = 1.3)  # default
    )
    
    # Calculate contour levels as GEOMETRIC PROGRESSION (same as Vizualisation.R)
    # Each level is contour_factor times the previous one
    contour_levels <- contour_start * contour_params$contour_factor^(0:(contour_params$contour_num - 1))
    
    #FIX: Defensive contour level computation to avoid warnings 
    # Previously: warnings appeared on zoom changes because:
    #   - z_max could be -Inf if region had no valid values
    #   - contour_levels could be empty if z_max < contour_start (small peak in zoomed view)
    #   - seq(contour_start, z_max, ...) could be decreasing if z_max < contour_start
    
    # Safely compute z_max
    valid_region <- region[is.finite(region)]
    z_max <- if (length(valid_region) > 0) max(valid_region) else NA_real_
    
    # Strict validation: only keep levels strictly within data range AND finite
    if (is.finite(z_max) && z_max > contour_start) {
      contour_levels <- contour_levels[
        is.finite(contour_levels) & 
          contour_levels >= contour_start & 
          contour_levels <= z_max
      ]
    } else {
      # No data above contour_start in this region -> no black contours possible
      contour_levels <- numeric(0)
    }
    
    # GGPLOT2 VISUALIZATION (same as base plot)
    
    # Prepare data for ggplot - same approach as Vizualisation.R
    # region has rows = F1 (y_sub), cols = F2 (x_sub)
    # We need to create a long-format data.frame
    
    # Create grid of coordinates
    grid_coords <- expand.grid(F1_idx = seq_along(y_sub), F2_idx = seq_along(x_sub))
    plot_df <- data.frame(
      F2 = x_sub[grid_coords$F2_idx],
      F1 = y_sub[grid_coords$F1_idx],
      z = as.vector(region)  # region is [F1, F2], as.vector goes column by column (F2 varies slower)
    )
    
    # Prepare fitted data if available
    fitted_df <- NULL
    fitted_contour_levels <- NULL
    
    if (fit_result$success && !is.null(fit_result$fitted_matrix)) {
      fitted_matrix <- fit_result$fitted_matrix
      fitted_df <- data.frame(
        F2 = x_sub[grid_coords$F2_idx],
        F1 = y_sub[grid_coords$F1_idx],
        z = as.vector(fitted_matrix)
      )
      
      # FIX: Defensive computation of fitted contour levels
      # Calculate contour levels adapted to the FITTED data
      valid_fitted <- fitted_matrix[is.finite(fitted_matrix)]
      z_max_fitted <- if (length(valid_fitted) > 0) max(valid_fitted) else NA_real_
      
      if (is.finite(z_max_fitted) && z_max_fitted > 0) {
        
        # FIX: Align red contours minimum with the experimental noise floor
        # Previously: fitted_start = z_max_fitted * 0.08 -> red contours could descend
        # WAY below contour_start (the spectrum's noise floor used for black contours).
        # Result: a Gaussian with even slightly too-large sigma would show huge red rings
        # at intensities where the experimental data shows nothing -> the "monstrous red
        # circle around a tiny black peak" problem (visible at low fit_method R² >= 0.99
        # because R² is computed only over the box, not the displayed padded region).
        # 
        # Now: fitted_start = max(contour_start, z_max_fitted * 0.05). At identical
        # intensity, red and black use the same threshold -> a well-fitted peak shows
        # red on top of black; queues of an over-wide Gaussian below noise floor are
        # automatically clipped (exactly what happens visually with the experimental
        # noise floor anyway).
        fitted_start <- max(contour_start, z_max_fitted * 0.05)
        
        fitted_contour_levels <- fitted_start * contour_params$contour_factor^(0:(contour_params$contour_num - 1))
        
        # Strict validation: levels must be finite AND strictly in data range
        fitted_contour_levels <- fitted_contour_levels[
          is.finite(fitted_contour_levels) &
            fitted_contour_levels >= fitted_start &
            fitted_contour_levels <= z_max_fitted
        ]
        
        # If we lost too many levels (e.g. fitted_start very close to z_max_fitted),
        # regenerate with linear spacing
        if (length(fitted_contour_levels) < 4 && fitted_start < z_max_fitted) {
          fitted_contour_levels <- seq(fitted_start, z_max_fitted * 0.99, length.out = 6)
        }
      } else {
        # No valid fitted data -> don't draw red contours
        fitted_df <- NULL
        fitted_contour_levels <- NULL
      }
    }
    
    # Build ggplot with EXPERIMENTAL contours in BLACK
    p <- ggplot(plot_df, aes(x = F2, y = F1, z = z)) +
      scale_x_reverse() +  # NMR convention
      scale_y_reverse() +  # NMR convention
      labs(x = "F2 (ppm)", y = "F1 (ppm)") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    
    # FIX: Add EXPERIMENTAL contours in BLACK only if levels are valid
    # Avoids "Zero contours generated" warning when zoom region has no data above contour_start.
    # Use aes(color = ...) with a literal label so ggplot generates a legend automatically.
    # Slightly thicker than before (0.5 -> 0.65) for better visibility against red contours.
    if (length(contour_levels) >= 2 && all(is.finite(contour_levels))) {
      p <- p + geom_contour(aes(color = "Experimental"), linewidth = 0.65, breaks = contour_levels)
    }
    
    # Add MODEL contours in RED if fit succeeded
    if (fit_result$success && !is.null(fitted_df) && 
        !is.null(fitted_contour_levels) && 
        length(fitted_contour_levels) >= 2 && 
        all(is.finite(fitted_contour_levels))) {
      # Use fitted_contour_levels which are adapted to the model's intensity range
      p <- p + geom_contour(data = fitted_df, aes(x = F2, y = F1, z = z, color = "Model fit"),
                            linewidth = 0.6, breaks = fitted_contour_levels)
      
      # FIX: Display ONLY fitted centers (coherent with red contours)
      # peak_centers now contains the actual fitted x0/y0 from each peak's nlsLM result,
      # so green crosses fall exactly on the centers of the red model lobes.
      # We removed the blue cross (box$center_x/center_y) because it came from
      # mod_integration's fit and could be a volume-weighted barycenter for multiplets,
      # which made it appear off-center relative to the red contours.
      if (!is.null(fit_result$peak_centers) && nrow(fit_result$peak_centers) > 0) {
        p <- p + geom_point(data = fit_result$peak_centers, 
                            aes(x = x, y = y, z = NULL, color = "Fitted center"),
                            shape = 3, size = 4, stroke = 1.8)
      }
    }
    
    # FIX: Color legend mapping labels -> colors
    # Maps the literal labels used in aes(color=...) to actual colors. The legend
    # appears automatically in the plot. We build the legend dynamically based on
    # which layers are actually drawn -> no override.aes mismatch when fit failed.
    legend_labels   <- "Experimental"
    legend_colors   <- c("Experimental" = "black")
    legend_linetype <- "solid"
    legend_shape    <- NA
    legend_lwidth   <- 0.65
    legend_size     <- NA
    
    has_red <- fit_result$success && !is.null(fitted_df) && 
      !is.null(fitted_contour_levels) && 
      length(fitted_contour_levels) >= 2 && 
      all(is.finite(fitted_contour_levels))
    has_green <- has_red && !is.null(fit_result$peak_centers) && 
      nrow(fit_result$peak_centers) > 0
    
    if (has_red) {
      legend_labels   <- c(legend_labels, "Model fit")
      legend_colors   <- c(legend_colors, "Model fit" = "#E31A1C")
      legend_linetype <- c(legend_linetype, "solid")
      legend_shape    <- c(legend_shape, NA)
      legend_lwidth   <- c(legend_lwidth, 0.6)
      legend_size     <- c(legend_size, NA)
    }
    if (has_green) {
      legend_labels   <- c(legend_labels, "Fitted center")
      legend_colors   <- c(legend_colors, "Fitted center" = "#2e7d32")
      legend_linetype <- c(legend_linetype, "blank")
      legend_shape    <- c(legend_shape, 3)
      legend_lwidth   <- c(legend_lwidth, NA)
      legend_size     <- c(legend_size, 4)
    }
    
    p <- p + scale_color_manual(
      name = NULL,
      values = legend_colors,
      breaks = legend_labels,
      guide = guide_legend(
        override.aes = list(
          linetype  = legend_linetype,
          shape     = legend_shape,
          linewidth = legend_lwidth,
          size      = legend_size
        )
      )
    ) + theme(
      legend.position = "top",
      legend.box.margin = margin(0, 0, 5, 0),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "lines")
    )
    
    # Title with peak info
    n_peaks_text <- if (fit_result$success && fit_result$n_peaks > 1) {
      paste0(box$stain_id, " - ", fit_result$model_used, " (", fit_result$n_peaks, " peaks)")
    } else if (fit_result$success) {
      paste0(box$stain_id, " - ", fit_result$model_used)
    } else {
      box$stain_id
    }
    
    # R² info for subtitle
    r2_text <- ifelse(is.na(box$r_squared), "N/A", sprintf("%.3f", box$r_squared))
    r2_col <- if (is.na(box$r_squared)) "gray50" else 
      if (box$r_squared >= 0.90) "#2e7d32" else 
        if (box$r_squared >= 0.80) "#1565c0" else 
          if (box$r_squared >= 0.70) "#f57c00" else "#c62828"
    
    p <- p + 
      ggtitle(n_peaks_text) +
      labs(subtitle = paste0("Method: ", box$fit_method, " | R² = ", r2_text))
    
    # Print the plot
    # suppressWarnings: filet de sécurité pour étouffer les warnings cosmétiques de
    # stat_contour ("Zero contours generated") et min/max sur vecteurs vides qui
    # peuvent émerger sur des configurations de zoom limites. Ces warnings n'indiquent
    # pas un bug : juste que la zone affichée n'a pas de données dans la plage des
    # niveaux de contour. Le rendu reste correct (plot vide ou partiel selon les cas).
    suppressWarnings(print(p))
    
  }, bg = "white")
  
  # Residuals plot
  
  output$residuals_plot <- renderPlot({
    req(load_data$bruker_data())
    box <- selected_fit_box()
    if (is.null(box) || nrow(box) == 0) {
      plot.new()
      text(0.5, 0.5, "Select a box in the\n'Fitted Boxes Details' table\nto see residuals", cex = 1.3, col = "gray50")
      return()
    }
    box <- box[1, ]  # Take the first row if multiple
    if (!"fit_method" %in% names(box) || is.na(box$fit_method) || box$fit_method %in% c("sum_fit_failed", "multiplet_sum") || grepl("sum_r2_below", box$fit_method)) {
      plot.new()
      text(0.5, 0.5, "No fit residuals available\n(sum method used or fit failed)", cex = 1.2, col = "orange")
      return()
    }
    
    # Re-fit to get residuals
    mat <- load_data$bruker_data()$spectrumData
    ppm_x <- suppressWarnings(as.numeric(colnames(mat)))
    ppm_y <- suppressWarnings(as.numeric(rownames(mat)))
    fit_result <- fit_2d_peak(mat, ppm_x, ppm_y, box, model = box$fit_method)
    if (is.null(fit_result$residuals) || length(fit_result$residuals) == 0) {
      plot.new()
      text(0.5, 0.5, "Could not compute residuals", cex = 1.5, col = "red")
      return()
    }
    
    # Residuals histogram
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
    hist(fit_result$residuals, breaks = 30, 
         col = "#f5576c", border = "white",
         main = paste0("Fit Residuals: ", box$stain_id),
         xlab = "Residual (Observed - Fitted)",
         ylab = "Frequency")
    abline(v = 0, col = "blue", lwd = 2, lty = 2)
    abline(v = mean(fit_result$residuals), col = "red", lwd = 2, lty = 2)
    
    # Add stats
    legend("topright", 
           legend = c(
             paste("Mean:", round(mean(fit_result$residuals), 2)),
             paste("SD:", round(sd(fit_result$residuals), 2)),
             paste("R²:", round(box$r_squared, 3)),
             "Blue = 0",
             "Red = Mean"
             
           ),
           bg = "white", cex = 0.9)
  })
  
  # Summary table
  
  output$fit_summary_table <- renderDT({
    summary <- fit_summary_data()
    if (is.null(summary) || nrow(summary) == 0) {
      return(datatable(data.frame(Message = "No fit data available. Run 'Batch Export' with Gaussian or Voigt method.")))
    }
    datatable(summary, 
              options = list(pageLength = 10, dom = 't'),
              rownames = FALSE) %>%
      formatRound(columns = c('mean_r2', 'median_r2', 'min_r2', 'max_r2'), digits = 3)
  })
  
  # Detailed table of fitted boxes
  
  output$fit_boxes_detail_table <- renderDT({
    boxes <- boxes_with_fit()
    if (is.null(boxes) || !"fit_method" %in% names(boxes)) {
      return(datatable(data.frame(Message = "No fit data available.")))
    }
    
    # Select relevant columns
    detail_cols <- c("stain_id", "r_squared", "fit_method", "center_x", "center_y")
    available_cols <- intersect(detail_cols, names(boxes))
    if (length(available_cols) == 0) {
      return(datatable(data.frame(Message = "No fit columns available.")))
    }
    detail_df <- boxes %>%
      select(all_of(available_cols)) %>%
      filter(!is.na(fit_method)) %>%
      arrange(desc(r_squared))
    
    # Rename for display
    names(detail_df) <- gsub("stain_id", "Box Name", names(detail_df))
    names(detail_df) <- gsub("r_squared", "R²", names(detail_df))
    names(detail_df) <- gsub("fit_method", "Method", names(detail_df))
    names(detail_df) <- gsub("center_x", "Center F2", names(detail_df))
    names(detail_df) <- gsub("center_y", "Center F1", names(detail_df))
    datatable(detail_df, 
              options = list(
                pageLength = 10, 
                scrollY = "300px",
                scrollCollapse = TRUE,
                dom = 'ftp'
                
              ),
              rownames = FALSE,
              selection = 'single') %>%
      formatRound(columns = c('R²', 'Center F2', 'Center F1'), digits = 3) %>%
      formatStyle('R²',
                  backgroundColor = styleInterval(c(0.7, 0.9), c('#ffcccc', '#ffffcc', '#ccffcc')))
  })
  
} # End server



## 3 RUN APPLICATION ----                                  


shinyApp(ui = ui, server = server)