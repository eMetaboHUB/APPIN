# tests/testthat/setup.R
# =============================================================================
# Charge automatiquement par testthat AVANT tous les test-*.R et helper-*.R.
#
# Role :
#  1) exposer APPIN_ROOT de facon fiable (poste Windows local OU runner CI) ;
#  2) sourcer tout le code metier (Function/*.R, R/*.R) dans l'env global,
#     car les test-*.R appellent ces fonctions sans les sourcer eux-memes.
#
# Regle keras : en CI, keras/tensorflow ne sont pas installes. Les fichiers
# qui en dependent (CNN detection/model, et CNN_shiny qui charge keras en
# tete) echouent au chargement -> on TOLERE leur echec (les tests associes
# skippent via skip_if_not_installed / skip_if_no_cnn_model). En revanche un
# fichier metier NON-ML qui echoue est une vraie erreur -> on s'arrete.
# =============================================================================

# ---- 1. Resolution de APPIN_ROOT -------------------------------------------
.find_appin_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (dir.exists(file.path(path, "Function")) &&
        dir.exists(file.path(path, "R"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  normalizePath("..", winslash = "/", mustWork = FALSE)
}

if (!exists("APPIN_ROOT", envir = globalenv())) {
  assign("APPIN_ROOT", .find_appin_root(), envir = globalenv())
}
.appin_root <- get("APPIN_ROOT", envir = globalenv())

# ---- 2. Fixtures ------------------------------------------------------------
if (!exists("FIXTURES_DIR", envir = globalenv())) {
  assign("FIXTURES_DIR",
         file.path(.appin_root, "tests", "fixtures"),
         envir = globalenv())
}

skip_if_no_fixture <- function(relative_path) {
  full <- file.path(get("FIXTURES_DIR", envir = globalenv()), relative_path)
  if (!file.exists(full) && !dir.exists(full)) {
    testthat::skip(sprintf("Fixture absente : %s", relative_path))
  }
  invisible(full)
}

# ---- 2bis. Attache les packages Shiny necessaires au chargement des modules
# Les R/mod_*.R utilisent des fonctions (ex. shinyDirButton de shinyFiles)
# sans prefixe pkg::, donc il faut attacher ces packages avant de les sourcer.
# Tolerant : si un package manque, on n'echoue pas ici.
for (pkg in c("shiny", "shinyFiles", "shinyjs", "shinyBS",
              "shinycssloaders", "shinydashboard", "shinydashboardPlus",
              "bslib", "bsicons", "DT", "plotly")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

# ---- 3. Chargement du code metier ------------------------------------------
# Fichiers REQUIS : leur echec de chargement doit faire echouer la suite
# (un package manquant ici = bug a corriger dans le workflow/DESCRIPTION).
.required_files <- c(
  "Function/Read_2DNMR_spectrum.R",
  "Function/Vizualisation.R",
  "Function/Peak_picking.R",
  "Function/Peak_fitting.R",
  "Function/CNN_clustering.R",   # DBSCAN pur, pas de keras -> doit charger
  "Function/CNN_filtering.R",    # filtrage pur, pas de keras -> doit charger
  "R/utils.R"
)

# Fichiers OPTIONNELS : peuvent tirer keras/tensorflow (absents en CI).
# Leur echec est tolere ; les tests correspondants skippent d'eux-memes.
.optional_files <- c(
  "Function/CNN_model.R",
  "Function/CNN_detection.R",
  "Function/CNN_main.R",
  "Function/CNN_shiny.R"
)

# Modules Shiny (R/mod_*.R) : sourced en optionnel. Ils ne sont necessaires
# qu'aux tests de modules ; s'ils tirent une dependance Shiny lourde au
# chargement, on ne veut pas casser toute la suite.
.module_files <- list.files(
  file.path(.appin_root, "R"),
  pattern = "^mod_.*\\.R$",
  full.names = FALSE
)
.module_files <- file.path("R", .module_files)

.source_one <- function(relpath, required) {
  full <- file.path(.appin_root, relpath)
  if (!file.exists(full)) {
    if (required) {
      stop(sprintf("setup.R : fichier REQUIS introuvable : %s", relpath),
           call. = FALSE)
    }
    message("setup.R : fichier absent, ignore -> ", relpath)
    return(invisible(FALSE))
  }
  ok <- tryCatch({
    sys.source(full, envir = globalenv())
    TRUE
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (required) {
      stop(sprintf(
        "setup.R : echec du source du fichier REQUIS '%s' : %s\n  -> dependance manquante ? Ajoute le package dans le workflow CI.",
        relpath, msg), call. = FALSE)
    }
    message("setup.R : [optionnel] echec du source de ", relpath, " -> ", msg)
    FALSE
  })
  if (isTRUE(ok)) message("setup.R : source -> ", relpath)
  invisible(ok)
}

# Ordre : requis d'abord, puis optionnels (CNN keras), puis modules Shiny.
for (f in .required_files) .source_one(f, required = TRUE)
for (f in .optional_files) .source_one(f, required = FALSE)
for (f in .module_files)   .source_one(f, required = FALSE)

message("setup.R : APPIN_ROOT = ", .appin_root)
