# =============================================================================
# 2D NMR Spectra Analysis - Launch Script
# =============================================================================
#
# This script serves as the main entry point for the APPIN application.
# It handles automatic dependency installation, environment setup, and app launch.
#
# FEATURES:
# - Automatic detection and installation of missing R packages
# - Version checking for critical packages (ensures compatibility)
# - Working directory auto-configuration (when run from RStudio)
# - Source file validation before launch
# - User-friendly console feedback with status indicators
#
# USAGE:
#   Option 1: Open this file in RStudio and click "Source" (or Ctrl+Shift+Enter)
#   Option 2: From R console: source("path/to/run_app.R")
#
# REQUIRED DIRECTORY STRUCTURE:
#   APPIN_1.0.0/
#   ├── run_app.R              <- THIS FILE (entry point)
#   ├── Shine.R                <- Main Shiny application UI/Server
#   └── Function/
#       ├── Read_2DNMR_spectrum.R   <- Bruker file parser
#       ├── Peak_fitting.R          <- Gaussian/Lorentzian fitting
#       ├── Vizualisation.R         <- Contour plot generation
#       ├── Peak_picking.R          <- Peak detection algorithms
#       └── CNN_shiny.R             <- Deep learning peak detection
#
# AUTHOR: Julien Guibert
# REPOSITORY: https://github.com/JulienGuibertTlse3/2DNMR-Analyst
# =============================================================================

cat("
╔══════════════════════════════════════════════════════════════════╗
║                    APPIN - Initialisation                        ║
╚══════════════════════════════════════════════════════════════════╝
\n")

# -----------------------------------------------------------------------------
# STEP 1: CONFIGURE WORKING DIRECTORY
# -----------------------------------------------------------------------------
# When running from RStudio, automatically set the working directory to the
# folder containing this script. This ensures all relative paths work correctly.
# If running from command line or non-interactive mode, user must ensure they
# are in the correct directory.

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  # RStudio environment: extract path from the active source editor
  script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
  if (nchar(script_path) > 0) {
    setwd(script_path)
    cat("📁 Working directory :", getwd(), "\n\n")
  }
} else {
  # Non-RStudio environment: display current directory and warn user
  cat("📁 Current working directory :", getwd(), "\n")
  cat("   (Make sure you are in the 2DNMR-Analyst folder)\n\n")
}

# -----------------------------------------------------------------------------
# STEP 2: DEFINE REQUIRED PACKAGES
# -----------------------------------------------------------------------------
# Complete list of all R packages needed by the application.
# Packages are grouped by functionality for easier maintenance.

packages_required <- c(
  # --- Shiny UI Framework ---
  # Core packages for building the interactive web interface
  "shiny",              # Base Shiny framework
  "shinyFiles",         # File/folder browser dialogs
  "shinydashboard",     # Dashboard layout components
  "shinydashboardPlus", # Extended dashboard features (boxes, cards)
  "shinyBS",            # Bootstrap components (tooltips, modals)
  "shinyjs",            # JavaScript operations from R
  "shinycssloaders",    # Loading spinners for async operations
  
  # --- Data Visualization ---
  # Packages for creating plots and interactive displays
  "plotly",             # Interactive plots (zoom, pan, hover)
  "ggplot2",            # Grammar of graphics plotting
  "DT",                 # Interactive data tables
  
  # --- Data Manipulation ---
  # Packages for efficient data processing and transformation
  "dplyr",              # Data wrangling verbs (filter, mutate, etc.)
  "data.table",         # Fast data operations for large datasets
  "magrittr",           # Pipe operator (%>%)
  "zoo",                # Rolling window functions
  "reshape2",           # Data reshaping (melt/cast)
  
  # --- Statistical Analysis ---
  # Packages for peak detection and signal processing
  "dbscan",             # Density-based clustering (DBSCAN algorithm)
  "sp",                 # Spatial data operations
  "matrixStats",        # Fast matrix statistics (row/col operations)
  "pracma",             # Numerical analysis functions
  "minpack.lm",         # Levenberg-Marquardt nonlinear least squares
  
  # --- Deep Learning ---
  # Packages for CNN-based peak detection (optional feature)
  "tensorflow",         # TensorFlow backend
  
  "keras",              # High-level deep learning API
  "imager",             # Image processing for CNN input
  
  # --- Utilities ---
  # Additional helper packages
  "Rcpp",               # C++ integration (required by some packages)
  "readr",              # Fast file reading
  "viridis",            # Color palettes for visualization
  "abind"               # Array binding operations
)


# -----------------------------------------------------------------------------
# STEP 3: CHECK AND INSTALL PACKAGES WITH VERSION REQUIREMENTS
# -----------------------------------------------------------------------------
# Some packages require minimum versions for compatibility.
# This section checks installed versions and updates if necessary.

# Minimum required versions for critical packages
# Format: package_name = "minimum.version.number"
required_versions <- list(
  # Uncomment these if specific TensorFlow/Keras versions are needed:
  # reticulate = "1.41.0",
  # tensorflow = "2.9.0",
  # keras = "2.9.0",
  shiny = "1.7.0"  # Minimum Shiny version for used features
)

cat("🔍 Checking critical packages...\n\n")

# Iterate through packages with version requirements
for (pkg_name in names(required_versions)) {
  required_version <- required_versions[[pkg_name]]
  needs_install <- FALSE
  
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    # Package not installed at all
    cat("   📦", pkg_name, "not installed\n")
    needs_install <- TRUE
  } else {
    # Package installed - check version
    current_version <- tryCatch(
      as.character(packageVersion(pkg_name)),
      error = function(e) "0.0.0"  # Fallback if version check fails
    )
    
    if (package_version(current_version) < package_version(required_version)) {
      # Installed version is too old
      cat("   ⚠️ ", pkg_name, current_version, "< required version", required_version, "\n")
      needs_install <- TRUE
    } else {
      # Version is sufficient
      cat("   ✅", pkg_name, current_version, "\n")
    }
  }
  
  # Install or update if needed
  if (needs_install) {
    cat("      → Installation/update of", pkg_name, "...\n")
    install.packages(pkg_name, dependencies = TRUE)
  }
}

# Check remaining packages (no specific version requirements)
cat("\n🔍 Checking for other required packages...\n\n")

# Find packages that are not yet installed
missing_packages <- packages_required[!sapply(packages_required, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("📦 Installing missing packages :", paste(missing_packages, collapse = ", "), "\n\n")
  install.packages(missing_packages, dependencies = TRUE)
}


# -----------------------------------------------------------------------------
# STEP 4: LOAD ALL PACKAGES INTO MEMORY
# -----------------------------------------------------------------------------
# Load all required packages, suppressing startup messages to keep console clean.

cat("📚 Loading packages...\n")

for (pkg in packages_required) {
  # suppressWarnings: ignore deprecation warnings during load
  # suppressPackageStartupMessages: hide package welcome messages
  suppressWarnings(suppressPackageStartupMessages(library(pkg, character.only = TRUE)))
}

cat("   ✅ All packages loaded\n")

# -----------------------------------------------------------------------------
# STEP 5: VERIFY SOURCE FILES EXIST
# -----------------------------------------------------------------------------
# Check that all required R source files are present before attempting to launch.
# This prevents cryptic errors during app startup.

cat("\n🔍 Source file verification...\n")

# List of function files that must exist in the Function/ subdirectory
source_files <- c(
  "Function/Read_2DNMR_spectrum.R",  # Bruker file reading
  "Function/Vizualisation.R",         # Contour plot generation
  "Function/Peak_picking.R",          # Peak detection
  "Function/Peak_fitting.R",          # Gaussian/Lorentzian fitting
  "Function/CNN_shiny.R"              # CNN-based detection
)

all_files_ok <- TRUE

# Check each function file
for (f in source_files) {
  if (file.exists(f)) {
    cat("   ✅", f, "\n")
  } else {
    cat("   ❌", f, "- MISSING!\n")
    all_files_ok <- FALSE
  }
}

# Check main application file
if (!file.exists("Shine.R")) {
  cat("   ❌ Shine.R - MISSING!\n")
  all_files_ok <- FALSE
} else {
  cat("   ✅ Shine.R\n")
}

# -----------------------------------------------------------------------------
# STEP 6: LAUNCH THE APPLICATION
# -----------------------------------------------------------------------------
# If all files are present, start the Shiny application.
# Otherwise, display an error message with instructions.

if (all_files_ok) {
  # All checks passed - launch the app
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════╗\n")
  cat("║                 🚀 App launch                                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("The application will open in your browser...\n")
  cat("To stop: click on STOP in RStudio or press Esc\n\n")
  
  # Launch the Shiny app from Shine.R
  # This will block until the app is closed
  shiny::runApp("Shine.R") 
  
} else {
  # Missing files - show error and instructions
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════╗\n")
  cat("║            ❌ ERROR: Missing files                               ║\n")
  cat("╚══════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Please check that you have downloaded all files from:\n")
  cat("https://github.com/JulienGuibertTlse3/2DNMR-Analyst\n\n")
  cat("Ensure the directory structure matches the expected layout shown\n")
  cat("at the top of this file.\n\n")
}