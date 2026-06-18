# CNN_shiny.R - CNN-based Peak Detection for 2D NMR Spectra ----

#
# This module provides CNN (Convolutional Neural Network) based peak detection
# for 2D NMR spectra analysis. It includes functions for model building,
# peak detection using different strategies (batch/sequential), filtering,
# and DBSCAN clustering for bounding box generation.
#
# Author: Julien Guibert
# Institution: INRAe Toxalim / MetaboHUB
#
# Structure:
#   CNN_shiny.R       - Main file (this file) - sources all modules
#   CNN_model.R       - Model architecture and loading
#   CNN_detection.R   - Peak detection (sequential + batch methods)
#   CNN_filtering.R   - Peak filtering functions
#   CNN_clustering.R  - DBSCAN clustering and bounding boxes
#   CNN_main.R        - Main entry point function


## DEPENDENCIES ----


library(purrr)
library(ggplot2)
library(keras)
library(viridis)
library(plotly)
library(reshape2)
library(abind)
library(tibble)
library(dplyr)
library(gridExtra)
library(dbscan)
library(stats)
library(readr)


## SOURCE MODULES ----


# Get the directory of this script
.cnn_script_dir <- tryCatch({
  if (exists("ofile", envir = sys.frame(1))) {
    dirname(sys.frame(1)$ofile)
  } else {
    "."
  }
}, error = function(e) ".")

# Source all CNN modules
source(file.path(.cnn_script_dir, "Function/CNN_model.R"))
source(file.path(.cnn_script_dir, "Function/CNN_detection.R"))
source(file.path(.cnn_script_dir, "Function/CNN_filtering.R"))
source(file.path(.cnn_script_dir, "Function/CNN_clustering.R"))
source(file.path(.cnn_script_dir, "Function/CNN_main.R"))

# Clean up
rm(.cnn_script_dir)
