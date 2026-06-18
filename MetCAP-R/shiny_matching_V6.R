# This file is part of [MetCAP-R]
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, 
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
# 
# See the LICENSE file for more details.

## This script run the shuiny app to match XCMS data with the database files created with the script "Databas Creation.R"
## Author: jfmartin


# shiny apps that match a peaklist (variable metadata) with the different AXIOM DB and CHEBI and HMDB
# FCtype define the subtitle of the plot
# The name of the metabolite if selected by selectInput and must be unique !!!!!!!!!!!!!!!!!!!!
# V00 : define the inputs files and parameters
# v01 : process matching (to finish adapt the GUI submit_match function to shiny)
# v02 : all message to check validity of the user choices input file, inhouse exist and contains at least mz column. mz and rt in peaklist...
# v03 : detect in inhouse DB mz or mass(monoistotopic) and adapt typeMatch (mzVSmz or mzVSmass)
# v04 : correction bug tolannot undefined if only annotcontaminant
# v05 : Modification of mz and rt columns for QTOF corresponding to database update and modification on HMDB matching due to bug tolannot and match2annot
# v06 : Suppression of modules for ORBI matching as OBRI DB has been suppressed

library(shiny)

## source fonctions for match_mass 
source("toolBox.R")
source("match_mass_script_new.R")

##################################  call matching functions ##########################################################

submit_match <- function(VMini, 
                         annotorbi, annotqtof, annotcontaminant, annotCHEBI, annotHMDB, annotinhouse,
                         mzcol, rtcol, useRT = TRUE, device,
                         inhousemass = 0 , inhousert = 0, inhouseMatch = "mzVSmass",
                         ioni, tolmzorbi, tolRTorbi, tolmzqtof,tolRTqtof,IsoAdduct
)
{
  
  #######################    define the tolerance for matching  #######################################
  ######################## for orbitrap (tolorbi) and QTOF (tolqtof)   ####################################
  tolorbi <- c(tolmzorbi,tolRTorbi)
  tolqtof <- c(tolmzqtof,tolRTqtof) 
  
  ############################## matching with the different DB ##########################################
  ########################## no modification are required in the next lines  #############################
  
  VM <- VMini
  
  ## column of m/z 
  mzcol <- which(colnames(VM) == mzcol)
  rtcol <- which(colnames(VM) == rtcol)
  
  if (!useRT) {
    rtcol <- 0
    cat("RT will not be used for matching only mz \n")
  }
  if (length(rtcol) > 1) rtcol <- rtcol[1]
  
  ## match columns of DB Orbi and qtof and contaminant "grave dans le marbre"
  mzdbcol <- 4
  rtdbcol <- 3
  mzdbcontaminant <- 25
  
  ### start matching
  ## QTOF
  if (annotqtof==TRUE) {
    tolannot <- tolqtof
    if (ioni == "negative") DBQTOF <- DBQTOFNEG else DBQTOF <- DBQTOFPOS
    cat("processing qtof \n")
    VM <- match2annot(x=VM,mzx=mzcol,rtx=rtcol,db=DBQTOF,mzdb=mzdbcol,rtdb=rtdbcol,ioni=ioni,tol=tolqtof[1],rtwin=tolqtof[2],typeMatch='mzVSmz')
  }
  
  ## ORBI
  if (annotorbi==TRUE) {
    tolannot <- tolorbi
    if (ioni == "negative") DBORBI <- DBORBINEG else DBORBI <- DBORBIPOS
    cat("processing orbi \n")
    VM <- match2annot(x=VM,mzx=mzcol,rtx=rtcol,db=DBORBI,mzdb=mzdbcol,rtdb=rtdbcol,ioni=ioni,tol=tolorbi[1],rtwin=tolorbi[2],typeMatch='mzVSmz')
  }
  
  ### Pesticides new contaminant EJA DB Screening_Contaminant.xlsx match sur mz envipat
  if (annotcontaminant==TRUE) {
    cat("processing contaminant IsoAdduct=",IsoAdduct,"\n")
    if (device == "msOrbi") tolannot <- tolorbi else  tolannot <- tolqtof
    if (ioni == "negative") {
      DBCP <- dbcontaminantsNEG
      if (IsoAdduct == 0) DBCP <- dbcontaminantsNEG[dbcontaminantsNEG$adduct.j.=="M-H" & dbcontaminantsNEG$currentIso=="[M]",] 
    }
    else {
      DBCP <- dbcontaminantsPOS
      if (IsoAdduct == 0) DBCP <- dbcontaminantsPOS[dbcontaminantsPOS$adduct.j.=="M+H" & dbcontaminantsPOS$currentIso=="[M]",] 
    }
    VM <- match2annot(x = VM,mzx = mzcol, rtx = 0, db = DBCP,mzdb=mzdbcontaminant, rtdb = 0,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmz')
  }
  
  ### Match CHEBI
  if (annotCHEBI==TRUE) VM <- match2annot(x=VM,mzx=mzcol,rtx=0, db=CHEBI,mzdb=5,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmass')
  ### Match HMDB
  if (annotHMDB == TRUE) {
    # Vérification et initialisation de tolannot si nécessaire
    if (!exists("tolannot") || is.null(tolannot)) {
      tolannot <- if (device == "msOrbi") {
        c(tolmzorbi, tolRTorbi)  # Tolérance pour Orbitrap
      } else {
        c(tolmzqtof, tolRTqtof)  # Tolérance pour QTOF
      }
      cat("tolannot initialisé pour HMDB : ", tolannot, "\n")
    }
    
    # Vérification finale avant matching
    if (is.null(tolannot)) {
      stop("Erreur : 'tolannot' n'est pas défini avant le matching HMDB.")
    }
    
    # Matching HMDB
    cat("Démarrage du matching HMDB...\n")
    VM <- match2annot(
      x = VM, 
      mzx = mzcol, 
      rtx = 0, 
      db = xHMDB, 
      mzdb = 4, 
      ioni = ioni, 
      tol = tolannot[1], 
      rtwin = tolannot[2], 
      typeMatch = 'mzVSmass'
    )
    cat("Matching HMDB terminé.\n")
  }
  
  
  #### matching with an inhouse standards TXT file (inhousert == 0 means no use of RT for matching)
  if (annotinhouse==TRUE) {
    
    if (inhousert == 0) rtcol <- 0 
    
    cat("processing inhouse mass:",inhousemass," rt ",inhousert," device ",device,"\n")
    
    if (device == "msOrbi") tolihdb <- tolorbi else tolihdb <- tolqtof
    cat("tolerance used=",tolihdb," ioni:",ioni," mzcol=",mzcol," rtcol=",rtcol,"\n")
    head(VM)
    
    VM <- match2annot(x = VM, mzx = mzcol, rtx = rtcol,
                      db = inhouseDB, mzdb = inhousemass, rtdb = inhousert,
                      ioni = ioni,
                      tol = tolihdb[1], rtwin = tolihdb[2],
                      typeMatch = inhouseMatch)
  }
  #outfil <- paste(infileVM,sep="")
  
  #write.table(x = VM, file=infileVM, sep="\t", quote = FALSE, row.names = FALSE, na="")
  
  return(VM)
  
}


###############################################  S H I N Y ############################################################
# Define UI for application that draws a barplot of factor effect 
ui <- fluidPage(
  titlePanel("Matching peaklist with AXIOM DB"),
  
  fluidRow( 
    column(3,
           h4("input peaklist, mz rt column and acquisition mode"),
           fileInput("vm" , "Choose a variable metadata file", accept = ".txt", multiple = FALSE),
           
           selectInput("mz", "Choose the m/z column", choices=c()),
           selectInput("rt", "Choose the retention time column", choices=c()),
           
           checkboxInput(inputId = "useRT", label =  "Click if you want to used retention time (default)",
                         value = TRUE),
           
           
           radioButtons("ioni", "ionisation",
                        c("positive" = "positive",
                          "negative" = "negative"))
           
    ),
    
    
    column(3, 
           h4("Choose AXIOM DB (Rdata files)"),
           fileInput("DBqtof" , "Choose the Acquity DB", accept = ".Rdata", multiple = FALSE),
           fileInput("DBcont" , "Choose the contaminants DB ", accept = ".Rdata", multiple = FALSE),
           
           checkboxInput(inputId = "IsoAdduc", label =  "Click if you want to match also with isotopes and adducts",
                         value = FALSE),
           
           fileInput("DBCHEBI", "Choose the CHEBI DB", accept = ".Rdata", multiple = FALSE),
           fileInput("DBHMDB" , "Choose the HMDB DB", accept = ".Rdata", multiple = FALSE)
    ),
    
    column(3,
           h4("Matching parameters"),
           sliderInput("tolmzqtof", "mz tolerance for Acquity (ppm)", value = 20,  min = 5, max = 50),
           sliderInput("tolRTqtof", "RT tolerance for Acquity (min)", value = 0.5,  min = 0, max = 2, step = 0.1)
    ),
    
    column(3, 
           h4("inhouse DB (TXT file)"),
           
           helpText("file must contains an ID (or name) and a monoisotopic mass column. if a rt column is present rt will be used for matchinng"),
           
           fileInput("DBinhouse" , "Choose in house DB (TXT file)", accept = ".Rdata", multiple = FALSE),
           
           actionButton(inputId = "button", label = "Click to submit matching"),
    
           downloadButton("annotfile", "Download annotation file"),
           
           tableOutput("tabVM"),
    ),
  )
)


server <- function(session, input, output) {
  
  cat("hello world\n")
  options(shiny.maxRequestSize = 30*1024^2)
  
  data1 <- reactive({
    validate(need(input$vm != "", "Please select a variableMetadata"))
    infile <<- input$vm
    #setwd(infile$datapath)
    ids <- read.table(file = infile$datapath, header = TRUE, sep="\t")
    return(ids)
  })
  
  data2 <- reactive({
    df3 <- data1()
    updateSelectInput(session,"mz" , choices = colnames(df3)) 
    updateSelectInput(session,"rt", choices = colnames(df3)) 
    return(df3)
    
  })
  
  output$tabVM <- renderTable({
    
    dfx <<- data2()
    head(dfx) 
    
  })  
  
  
  observeEvent(input$button, {
    
    # flags that define which DB will be used for matching
    annotorbi <- FALSE
    annotqtof <- FALSE
    annotcontaminant <- FALSE
    annotCHEBI <- FALSE
    annotHMDB <- FALSE
    
    ## init values for inhouse matching DB
    annotinhouse <- FALSE
    ihmass <- 0
    ihrt <- 0
    
    ## user selected parameters problems managment 
    err <- 0
    
    if (!is.null(input$DBorbi)) {
      inOrbi <- input$DBorbi
      load(file = inOrbi$datapath, .GlobalEnv)
      #str(DBORBINEG)
      annotorbi <- TRUE
    }
    if (!is.null(input$DBqtof)) {
      inQtof <- input$DBqtof
      load(file = inQtof$datapath, .GlobalEnv)
      #str(DBQTOFNEG)
      annotqtof <- TRUE
    }
    if (!is.null(input$DBcont)) {
      inCont <- input$DBcont
      load(file = inCont$datapath, .GlobalEnv)
      #str(dbcontaminantsNEG)
      annotcontaminant <- TRUE
    }
    if (!is.null(input$DBCHEBI)) {
      inCHEBI <- input$DBCHEBI
      load(file = inCHEBI$datapath, .GlobalEnv)
      #str(CHEBI)
      annotCHEBI <- TRUE
    }
    if (!is.null(input$DBHMDB)) {
      inHMDB <- input$DBHMDB
      load(file = inHMDB$datapath, .GlobalEnv)
      #str(xHMDB)
      annotHMDB <- TRUE
    }
    if (!is.null(input$DBinhouse)) {
      inhouse <- input$DBinhouse
      inhouseDB <<- read.table(file = inhouse$datapath, header=TRUE, sep="\t", stringsAsFactors = FALSE)
      annotinhouse <- TRUE
      ihmass <- which(colnames(inhouseDB) == "MASS" | colnames(inhouseDB) == "mass")
      
      ## check if monoIsotopic mass or m/z are in the inhouse DB and set the parameter match type
      if (length(ihmass) > 0) {
        ihMatch <- "mzVSmass" 
      } else 
      {
        ihmass <- which(colnames(inhouseDB) == "MZ" | colnames(inhouseDB) == "mz" | colnames(inhouseDB) == "m.z")
        if (length(ihmass) > 0)  ihMatch <- "mzVSmz" 
      }
      
      ihrt <- which(colnames(inhouseDB) == "RT" | colnames(inhouseDB) == "rt")
      ## if no column ret time is detected then ihrt is set to 0
      if (length(ihrt) == 0) ihrt <- 0
      ## if no mass or MASS or mz column in inhouseDB -> error
      if (length(ihmass) == 0) err <- 5
    }
    
    ## no DB selected 
    if (!annotorbi & !annotqtof & !annotcontaminant & !annotCHEBI & !annotHMDB & !annotinhouse) err <- 2
    
    ## no input vm selected
    if (is.null(input$vm)) {err <- 1 }
    else {
      ## mz or rt non numeric in the input vm selected
      imz <- which(colnames(dfx) == input$mz)
      if (!is.numeric(dfx[[imz]])) err <- 3 
      
      irt <- which(colnames(dfx) == input$rt)
      if (!is.numeric(dfx[[irt]]) & input$useRT) err <- 4 
    }
    
    
    ## execution if no error detected
    if (err == 0) {
      cat("device ",input$msdevice,"\n")
      
      ## progression bar during processing
      withProgress(message = 'Matching in process', value = 0, {
        
        resMatch <- submit_match(VMini = dfx,
                                 annotorbi, annotqtof, annotcontaminant, annotCHEBI, annotHMDB, annotinhouse,
                                 mzcol = input$mz, rtcol = input$rt,
                                 useRT = input$useRT,
                                 device = input$msdevice,
                                 inhousemass = ihmass, inhousert = ihrt, inhouseMatch = ihMatch,
                                 ioni = input$ioni, 
                                 tolmzorbi = input$tolmzorbi, 
                                 tolRTorbi = input$tolRTorbi, 
                                 tolmzqtof = input$tolmzqtof,
                                 tolRTqtof = input$tolRTqtof,
                                 IsoAdduc  = input$IsoAdduc
        )
        
      } 
      ) ## end withProgress
      
      
      ##write.table(x = resMatch, file = "annot.txt", row.names = FALSE, quote = FALSE, na = "", sep = "\t")
      ## define download nam of file and action write.table in link with he downloadButton in ui 
      output$annotfile <- downloadHandler(
        filename = function() { "annot.txt"},
        content = function(file) { write.table(resMatch, file, sep="\t",row.names = FALSE) }
      )
    } else 
    {
      if (err == 1) cat("No input VM selected")
      if (err == 2) cat("No DB selected")
      if (err == 3) cat(" mz selected is not a numeric columns")
      if (err == 4) cat(" rt selected is not a numeric columns")
      if (err == 5) cat(" no mass column in inhouse DB")
    }
    
    cat(".... Finished\n")
    
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
