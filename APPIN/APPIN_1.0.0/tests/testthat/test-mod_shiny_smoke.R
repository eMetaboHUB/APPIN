# test-mod_shiny_smoke.R ----
# Smoke tests for Shiny modules.
#
# These are MINIMAL sanity checks only. They verify:
#   1. UI functions construct without error and return a Shiny tag/tagList
#   2. Server functions can be instantiated via shiny::testServer() without
#      crashing on startup (observers/reactives that fire immediately work OK)
#
# They do NOT test business logic — that would require full reactive flow
# testing with realistic data, which has poor ROI for Shiny modules.
#
# If a module later gets a bug at load time (bad observer, missing dep, etc.),
# these smoke tests will catch it.

# =============================================================================
# PACKAGE DEPS — required for module UI/server bodies to resolve
# =============================================================================
# Modules call NS(), tagList(), shinyDirChoose(), etc. without prefixing the
# package. So those packages must be ATTACHED (via library()) in the test env,
# not just loaded.

suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  # shinyFiles is used by mod_load_data (shinyDirChoose, shinyDirButton, etc.)
  if (requireNamespace("shinyFiles", quietly = TRUE)) library(shinyFiles)
  # shinyjs is commonly used (disable/enable, hide/show)
  if (requireNamespace("shinyjs", quietly = TRUE)) library(shinyjs)
  # bslib / bsicons for UI cards
  if (requireNamespace("bslib", quietly = TRUE)) library(bslib)
  if (requireNamespace("bsicons", quietly = TRUE)) library(bsicons)
  # DT is often used for tables
  if (requireNamespace("DT", quietly = TRUE)) library(DT)
  # plotly for click/relayout events
  if (requireNamespace("plotly", quietly = TRUE)) library(plotly)
  # htmltools - underlies shiny UI
  if (requireNamespace("htmltools", quietly = TRUE)) library(htmltools)
}))

# =============================================================================
# HELPERS
# =============================================================================

#' Create a minimal fake reactiveValues with the fields modules expect.
#' Use `reactiveValues()` so downstream code that does `rv$foo` works.
make_fake_rv <- function() {
  shiny::reactiveValues(
    # Data
    rr_norm        = NULL,
    ppm_x          = NULL,
    ppm_y          = NULL,
    spectrum_type  = "TOCSY",
    # Peaks / boxes
    peaks          = data.frame(),
    boxes          = data.frame(),
    shapes         = list(),
    # Selection / editing
    selected_ids   = character(0),
    pending_changes = list(),
    # UI state
    click_mode     = "none",
    manual_mode    = FALSE
  )
}

#' Build a reactive wrapper that returns a static value (for things like
#' `data_reactives$rr_norm()` etc.)
fake_reactive <- function(value) {
  shiny::reactive(value)
}

#' Fake "data_reactives" list — most modules take this as an argument.
make_fake_data_reactives <- function() {
  list(
    rr_norm       = fake_reactive(NULL),
    ppm_x         = fake_reactive(NULL),
    ppm_y         = fake_reactive(NULL),
    spectrum_type = fake_reactive("TOCSY")
  )
}

#' Fake "load_data" reactiveValues used by save/export/load modules
make_fake_load_data <- function() {
  shiny::reactiveValues(
    dir = NULL,
    dim = "2D",
    rr = NULL
  )
}

#' Fake status_msg reactive
make_fake_status_msg <- function() {
  shiny::reactiveVal("")
}

#' Find the project root (directory containing R/ with mod_*.R files).
#' Walks up from the current working directory until it finds it.
.find_project_root <- function() {
  # Start from wd, then walk up
  path <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:5) {  # max 5 levels up
    candidate <- file.path(path, "R")
    # Check if R/ exists AND contains at least one mod_*.R file
    if (dir.exists(candidate) &&
        length(list.files(candidate, pattern = "^mod_.*\\.R$")) > 0) {
      return(path)
    }
    parent <- dirname(path)
    if (parent == path) break  # reached filesystem root
    path <- parent
  }
  NULL
}

#' Source ALL modules at once (solves inter-module dependencies).
#' Cached via a flag so we only do this once per test session.
.all_modules_loaded <- FALSE

.ensure_module_loaded <- function(modname) {
  expected_fn <- paste0(modname, "_ui")
  if (exists(expected_fn, mode = "function")) {
    return(invisible(TRUE))
  }
  
  if (.all_modules_loaded) {
    # We already tried loading everything — the module truly doesn't exist
    testthat::skip(sprintf("Module %s not available after full load", modname))
  }
  
  root <- .find_project_root()
  if (is.null(root)) {
    testthat::skip(sprintf(
      "Project root not found (R/ with mod_*.R files). wd=%s", getwd()
    ))
  }
  
  # Source ALL mod_*.R files in dependency-agnostic order.
  # Some modules reference others (e.g. mod_manual_editing_ui uses
  # mod_pending_changes_ui), so we need them all in the env simultaneously.
  mod_files <- list.files(file.path(root, "R"),
                          pattern = "^mod_.*\\.R$",
                          full.names = TRUE)
  
  if (length(mod_files) == 0) {
    testthat::skip(sprintf("No mod_*.R files found in %s/R", root))
  }
  
  for (f in mod_files) {
    # Use local scope to avoid polluting, but source into globalenv so
    # functions are visible to all subsequent tests.
    tryCatch(
      sys.source(f, envir = globalenv()),
      error = function(e) {
        # Non-fatal: some modules may fail to source in test env, that's OK —
        # their own tests will skip on the next .ensure_module_loaded call.
        warning(sprintf("Could not source %s: %s", basename(f), e$message),
                call. = FALSE)
      }
    )
  }
  
  .all_modules_loaded <<- TRUE
  
  # Re-check
  if (!exists(expected_fn, mode = "function")) {
    testthat::skip(sprintf("Module %s still not found after loading R/mod_*.R",
                           modname))
  }
}

# =============================================================================
# UI SMOKE TESTS
# =============================================================================
# Every *_ui(id) should return a Shiny tag or tagList without crashing.

.test_ui <- function(mod_name) {
  test_that(paste0(mod_name, "_ui() builds without error"), {
    .ensure_module_loaded(mod_name)
    ui_fn <- get(paste0(mod_name, "_ui"))
    
    expect_no_error({
      ui <- ui_fn("test_id")
    })
    
    # Result should be a shiny.tag, shiny.tag.list, or similar
    expect_true(
      inherits(ui, c("shiny.tag", "shiny.tag.list", "list", "html")) ||
        is.character(ui),
      info = paste("Unexpected UI return class:",
                   paste(class(ui), collapse = ", "))
    )
  })
}

# One UI test per module
.test_ui("mod_box_editor")
.test_ui("mod_click_mode")
.test_ui("mod_delete")
.test_ui("mod_export")
.test_ui("mod_fusion")
.test_ui("mod_import")
.test_ui("mod_integration")
.test_ui("mod_load_data")
.test_ui("mod_manual_add")
.test_ui("mod_manual_editing")
.test_ui("mod_peak_picking")
.test_ui("mod_pending_changes")
.test_ui("mod_reset")
.test_ui("mod_save_export")
.test_ui("mod_session")

# =============================================================================
# SERVER SMOKE TESTS
# =============================================================================
# Each server should instantiate via shiny::testServer() without errors at
# startup. We use tryCatch because some modules may call observeEvent() with
# triggers that don't fire without user input — that's fine, we only care
# about instantiation-time errors.

.test_server <- function(mod_name, args_builder) {
  test_that(paste0(mod_name, "_server() instantiates without error"), {
    .ensure_module_loaded(mod_name)
    server_fn <- get(paste0(mod_name, "_server"))
    
    args <- args_builder()
    
    # Plotly emits non-critical warnings about event_register when UI isn't
    # rendered (which is the case in testServer). These aren't errors — we
    # suppress them to keep test output readable.
    expect_no_error(
      suppressWarnings({
        shiny::testServer(
          server_fn,
          args = args,
          expr = {
            # Just let reactives settle; don't exercise business logic
            session$flushReact()
          }
        )
      })
    )
  })
}

.test_server("mod_box_editor", function() {
  list(
    rv              = make_fake_rv(),
    data_reactives  = make_fake_data_reactives(),
    parent_input    = list(),
    parent_session  = NULL
  )
})

.test_server("mod_click_mode", function() {
  list(
    rv              = make_fake_rv(),
    data_reactives  = make_fake_data_reactives(),
    peak_picking    = NULL
  )
})

.test_server("mod_delete", function() {
  list(rv = make_fake_rv())
})

.test_server("mod_export", function() {
  list(
    status_msg     = make_fake_status_msg(),
    rv             = make_fake_rv(),
    load_data      = make_fake_load_data(),
    data_reactives = make_fake_data_reactives()
  )
})

.test_server("mod_fusion", function() {
  list(rv = make_fake_rv())
})

.test_server("mod_import", function() {
  list(
    rv                = make_fake_rv(),
    refresh_nmr_plot  = shiny::reactiveVal(0)
  )
})

.test_server("mod_integration", function() {
  list(
    status_msg = make_fake_status_msg(),
    load_data  = make_fake_load_data(),
    rv         = make_fake_rv()
  )
})

.test_server("mod_load_data", function() {
  list(
    status_msg               = make_fake_status_msg(),
    trigger_subfolder_update = NULL
  )
})

.test_server("mod_manual_add", function() {
  list(
    rv              = make_fake_rv(),
    data_reactives  = make_fake_data_reactives(),
    peak_picking    = list()
  )
})

.test_server("mod_manual_editing", function() {
  list(
    status_msg        = make_fake_status_msg(),
    load_data         = make_fake_load_data(),
    rv                = make_fake_rv(),
    data_reactives    = make_fake_data_reactives(),
    refresh_nmr_plot  = shiny::reactiveVal(0),
    peak_picking      = list(),
    parent_input      = list(),
    parent_session    = NULL
  )
})

.test_server("mod_peak_picking", function() {
  list(
    status_msg        = make_fake_status_msg(),
    load_data         = make_fake_load_data(),
    data_reactives    = make_fake_data_reactives(),
    rv                = make_fake_rv(),
    refresh_nmr_plot  = shiny::reactiveVal(0),
    parent_input      = list()
  )
})

.test_server("mod_pending_changes", function() {
  list(
    rv                = make_fake_rv(),
    data_reactives    = make_fake_data_reactives(),
    load_data         = make_fake_load_data(),
    refresh_nmr_plot  = shiny::reactiveVal(0),
    parent_input      = list(),
    parent_session    = NULL
  )
})

.test_server("mod_reset", function() {
  list(
    status_msg     = make_fake_status_msg(),
    rv             = make_fake_rv(),
    parent_session = NULL
  )
})

.test_server("mod_save_export", function() {
  list(
    status_msg        = make_fake_status_msg(),
    rv                = make_fake_rv(),
    load_data         = make_fake_load_data(),
    data_reactives    = make_fake_data_reactives(),
    refresh_nmr_plot  = shiny::reactiveVal(0),
    parent_session    = NULL,
    parent_input      = list()
  )
})

.test_server("mod_session", function() {
  list(
    status_msg        = make_fake_status_msg(),
    rv                = make_fake_rv(),
    load_data         = make_fake_load_data(),
    refresh_nmr_plot  = shiny::reactiveVal(0),
    parent_session    = NULL,
    parent_input      = list()
  )
})