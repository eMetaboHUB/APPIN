#' Read Bruker NMR Spectrum Files
#'
#' Reads and parses Bruker NMR spectrum files (1D or 2D) including metadata
#' and spectral data from processed data directories.
#'
#' @param dir Character. Path to the processed spectrum directory (typically pdata/1/).
#' @param dim Character. Dimension of the spectrum: "1D" or "2D".
#' @param onlyTitles Logical. If TRUE, only reads metadata without loading spectrum data.
#'   Default is FALSE.
#' @param useAsNames Character. Naming convention for the spectrum:
#'   \itemize{
#'     \item "Spectrum titles" - Uses the title from the title file
#'     \item "dir names" - Uses the parent folder name
#'     \item "dir names and EXPNO" - Uses folder name with experiment number (e.g., "sample_0001")
#'   }
#'   Default is "Spectrum titles".
#' @param checkFiles Logical. If TRUE, only verifies that required files exist
#'   without reading data. Default is FALSE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{spectrumData} - Matrix (2D) or vector (1D) of spectral intensities,
#'       with chemical shifts as row/column names
#'     \item \code{spectrumDataName} - Name based on \code{useAsNames} parameter
#'     \item \code{spectrumDataTitle} - Title from the Bruker title file
#'     \item \code{spectrumDataFolderName} - Parent folder name
#'     \item \code{spectrumDataEXPNO} - Experiment number
#'     \item \code{spectrumDataFolderName_EXPNO} - Combined folder name and EXPNO
#'     \item \code{AcquPars} - List of acquisition parameters (NS, BF1, P1, RG, PULPROG, SOLVENT)
#'   }
#'   Returns NULL if \code{dir} is NULL.
#'
#' @details
#' This function reads Bruker TopSpin processed NMR data. The directory should point
#' to a pdata/X/ folder containing procs (and proc2s for 2D) parameter files and
#' the binary spectrum file (1r for 1D, 2rr for 2D).
#'
#' The function handles Bruker's block storage format for 2D spectra, reorganizing
#' the data into a properly ordered matrix with chemical shift axes.
#'
#' @examples
#' \dontrun{
#' # Read a 2D spectrum
#' spectrum <- read_bruker(
#'   dir = "/path/to/experiment/pdata/1",
#'   dim = "2D"
#' )
#'
#' # Only check if files exist
#' read_bruker(dir = "/path/to/experiment/pdata/1", dim = "2D", checkFiles = TRUE)
#'
#' # Read only metadata
#' metadata <- read_bruker(
#'   dir = "/path/to/experiment/pdata/1",
#'   dim = "2D",
#'   onlyTitles = TRUE
#' )
#' }
#'
#' @export


### Read Metadata from Bruker files ----
read_bruker <- function (dir = NULL, dim = NULL, onlyTitles = FALSE, 
                         useAsNames = "Spectrum titles", checkFiles = FALSE) 
{
  # Initialize default acquisition parameters
  AcquPars <- list(NS = 0, BF1 = 0, P1 = 0, RG = 0, PULPROG = "", 
                   SOLVENT = "")
  
  # Map dimension labels to Bruker binary file names
  # 1D spectra use "1r", 2D spectra use "2rr"
  dnDim <- c("1r", "2rr")
  names(dnDim) <- c("1D", "2D")
  spectrum_proc_path <- dir
  datanameTmp <- dnDim[dim]
  
  if (!is.null(spectrum_proc_path)) {
    # Bruker byte order: 0 = little endian, 1 = big endian
    BYTORDP_Dict <- c("little", "big")
    names(BYTORDP_Dict) <- c(0, 1)
    TITLE <- ""
    
    # Try to read spectrum title from the title file
    try(TITLE <- scan(file = paste(spectrum_proc_path, "/title", 
                                   sep = ""), what = "character", sep = "\n", quiet = TRUE)[1], 
        silent = TRUE)
    
    # --- Check files mode: only verify required files exist ---
    if (checkFiles) {
      spectrumData <- NULL
      titleFinal <- NULL
      spectrumDataTitle <- NULL
      spectrumDataFolderName <- NULL
      spectrumDataEXPNO <- NULL
      spectrumDataFolderName_EXPNO <- NULL
      list.filesTMP <- list.files(spectrum_proc_path)
      if (!"procs" %in% list.filesTMP) {
        stop(paste("Could not open spectrum", spectrum_proc_path))
      }
    }
    else {
      if (!onlyTitles) {
        # --- Read acquisition parameters from acqu file ---
        # Navigate from pdata/X/ back to experiment folder to find acqu file
        spectrum_acqu_path <- NULL
        spectrum_acqu_pathTMP <- strsplit(spectrum_proc_path, 
                                          "/")[[1]]
        spectrum_acqu_path <- paste(spectrum_acqu_pathTMP[1:(length(spectrum_acqu_pathTMP) - 
                                                               2)], sep = "/", collapse = "/")
        acqusTMP <- NULL
        try(acqusTMP <- scan(file = paste(spectrum_acqu_path, 
                                          "/acqu", sep = ""), what = "character", sep = "\n", 
                             quiet = TRUE), silent = TRUE)
        
        if (!(is.null(acqusTMP) | (length(acqusTMP) == 0))) {
          # Clean Bruker parameter file formatting (remove # and $ characters)
          acqusTMP <- gsub("#", "", acqusTMP)
          acqusTMP <- gsub("\\$", " ", acqusTMP)
          
          # Extract acquisition parameters: RG (receiver gain)
          RGpos <- grep(" RG=", acqusTMP)
          if (length(RGpos) > 0) 
            AcquPars$RG <- as.numeric(strsplit(acqusTMP[RGpos], 
                                               split = "= ")[[1]][2])
          
          # NS: number of scans
          NSpos <- grep(" NS=", acqusTMP)
          if (length(NSpos) > 0) 
            AcquPars$NS <- as.numeric(strsplit(acqusTMP[NSpos], 
                                               split = "= ")[[1]][2])
          
          # BF1: base frequency of channel 1 (MHz)
          BF1pos <- grep(" BF1=", acqusTMP)
          if (length(BF1pos) > 0) 
            AcquPars$BF1 <- as.numeric(strsplit(acqusTMP[BF1pos], 
                                                split = "= ")[[1]][2])
          
          # PULPROG: pulse program name
          PULPROGpos <- grep(" PULPROG=", acqusTMP)
          if (length(PULPROGpos) > 0) 
            AcquPars$PULPROG <- strsplit(acqusTMP[PULPROGpos], 
                                         split = "= ")[[1]][2]
          
          # SOLVENT: solvent used
          SOLVENTpos <- grep(" SOLVENT=", acqusTMP)
          if (length(SOLVENTpos) > 0) 
            AcquPars$SOLVENT <- strsplit(acqusTMP[SOLVENTpos], 
                                         split = "= ")[[1]][2]
          
          # P1: pulse width (value is on the line AFTER "P=")
          P1pos <- grep(" P=", acqusTMP)
          if (length(P1pos) > 0) 
            AcquPars$P1 <- as.numeric(strsplit(acqusTMP[P1pos + 
                                                          1], , split = " ")[[1]][2])
        }
        
        # --- Read processing parameters from procs file (F2/direct dimension) ---
        proc <- scan(file = paste(spectrum_proc_path, 
                                  "/procs", sep = ""), what = "character", sep = "\n", 
                     quiet = TRUE)
        proc <- gsub("#", "", proc)
        proc <- gsub("\\$", " ", proc)
        
        # SI: size (number of points) in F2 dimension
        SI2 <- as.numeric(strsplit(proc[grep(" SI=", 
                                             proc)], split = "= ")[[1]][2])
        
        # BYTORDP: byte order for binary data
        BYTORDP <- BYTORDP_Dict[strsplit(proc[grep(" BYTORDP=", 
                                                   proc)], split = "= ")[[1]][2]]
        
        # NC_proc: intensity scaling factor (data = raw * 2^NC_proc)
        NC_proc <- as.numeric(strsplit(proc[grep(" NC_proc=", 
                                                 proc)], split = "= ")[[1]][2])
        
        # XDIM: submatrix size for 2D block storage
        XDIM2 <- as.numeric(strsplit(proc[grep(" XDIM=", 
                                               proc)], split = "= ")[[1]][2])
        
        # OFFSET: chemical shift of first point (left edge, ppm)
        OFFSET2 <- as.numeric(strsplit(proc[grep(" OFFSET=", 
                                                 proc)], split = "= ")[[1]][2])
        
        # SF: spectrometer frequency (MHz)
        SF2 <- as.numeric(strsplit(proc[grep(" SF=", 
                                             proc)], split = "= ")[[1]][2])
        
        # SW_p: spectral width in Hz
        SW_p2 <- as.numeric(strsplit(proc[grep(" SW_p=", 
                                               proc)], split = "= ")[[1]][2])
        
        # Calculate chemical shift axis for F2 dimension
        rightlimit2 <- OFFSET2 - SW_p2/SF2
        leftlimit2 <- OFFSET2
        frequencynames2 <- OFFSET2 - (0:(SI2 - 1)) * SW_p2/SF2/SI2
        n <- SI2
        
        # --- For 2D spectra: read proc2s file (F1/indirect dimension) ---
        if (datanameTmp == "2rr") {
          proc2 <- scan(file = paste(spectrum_proc_path, 
                                     "/proc2s", sep = ""), what = "character", 
                        sep = "\n", quiet = TRUE)
          proc2 <- gsub("#", "", proc2)
          proc2 <- gsub("\\$", " ", proc2)
          
          # Extract F1 dimension parameters (same as F2 but for indirect dimension)
          SI1 <- as.numeric(strsplit(proc2[grep(" SI=", 
                                                proc2)], split = "= ")[[1]][2])
          XDIM1 <- as.numeric(strsplit(proc2[grep(" XDIM=", 
                                                  proc2)], split = "= ")[[1]][2])
          OFFSET1 <- as.numeric(strsplit(proc2[grep(" OFFSET=", 
                                                    proc2)], split = "= ")[[1]][2])
          SF1 <- as.numeric(strsplit(proc2[grep(" SF=", 
                                                proc2)], split = "= ")[[1]][2])
          SW_p1 <- as.numeric(strsplit(proc2[grep(" SW_p=", 
                                                  proc2)], split = "= ")[[1]][2])
          
          # Total number of points for 2D = F1 x F2
          n <- SI2 * SI1
          
          # Calculate chemical shift axis for F1 dimension
          rightlimit1 <- OFFSET1 - SW_p1/SF1
          leftlimit1 <- OFFSET1
          frequencynames1 <- OFFSET1 - (0:(SI1 - 1)) * SW_p1/SF1/SI1
        }
        
        # --- Read binary spectrum data ---
        # Bruker stores intensities as 32-bit integers (size = 4 bytes)
        spectrumData <- readBin(paste(spectrum_proc_path, 
                                      "/", datanameTmp, sep = ""), what = "integer", 
                                size = 4, n = n, endian = BYTORDP)
        
        # Apply intensity scaling factor
        spectrumData <- spectrumData * 2^NC_proc
        
        # --- Reorganize 2D data from Bruker block storage format ---
        # Bruker stores 2D data in submatrices (blocks) of size XDIM1 x XDIM2
        # These blocks need to be reassembled into the correct matrix layout
        if (datanameTmp == "2rr") {
          spectrumDataTMP <- spectrumData
          spectrumData <- matrix(spectrumData, 
                                 ncol = SI2, byrow = TRUE)
          counter <- 0
          
          # Loop through all blocks and place them in correct position
          for (j in 1:(nrow(spectrumData)/XDIM1)) {
            for (i in 1:(ncol(spectrumData)/XDIM2)) {
              spectrumData[(j - 1) * XDIM1 + (1:XDIM1), 
                           (i - 1) * XDIM2 + (1:XDIM2)] <- matrix(spectrumDataTMP[counter * 
                                                                                    XDIM2 * XDIM1 + (1:(XDIM2 * XDIM1))], 
                                                                  ncol = XDIM2, byrow = TRUE)
              counter <- counter + 1
            }
          }
          
          # Assign chemical shift values as row/column names
          rownames(spectrumData) <- frequencynames1
          colnames(spectrumData) <- frequencynames2
        }
        else {
          # For 1D: assign chemical shifts as names
          names(spectrumData) <- frequencynames2
        }
      }
      else {
        # onlyTitles mode: skip spectrum data loading
        spectrumData <- NULL
      }
      
      # --- Extract folder and experiment information from path ---
      spectrumDataTitle <- TITLE
      # Folder name is 4 levels up from pdata/X/ (e.g., "sample_name")
      spectrumDataFolderName <- rev(strsplit(spectrum_proc_path, 
                                             "/")[[1]])[4]
      # EXPNO is 3 levels up (e.g., "10")
      spectrumDataEXPNO <- rev(strsplit(spectrum_proc_path, 
                                        "/")[[1]])[3]
      
      # Create formatted name with zero-padded EXPNO (e.g., "sample_name_0010")
      spectrumDataFolderName_EXPNO <- paste(spectrumDataFolderName, 
                                            paste(c("_", "0", "0", "0", "0")[1:max(1, 5 - 
                                                                                     nchar(spectrumDataEXPNO))], sep = "", collapse = ""), 
                                            spectrumDataEXPNO, sep = "")
      
      # Set final name based on user preference
      if (useAsNames == "Spectrum titles") 
        titleFinal <- spectrumDataTitle
      if (useAsNames == "dir names") 
        titleFinal <- spectrumDataFolderName
      if (useAsNames == "dir names and EXPNO") 
        titleFinal <- spectrumDataFolderName_EXPNO
    }
    
    # Return all extracted data as a list
    invisible(list(spectrumData = spectrumData, spectrumDataName = titleFinal, 
                   spectrumDataTitle = spectrumDataTitle, spectrumDataFolderName = spectrumDataFolderName, 
                   spectrumDataEXPNO = spectrumDataEXPNO, spectrumDataFolderName_EXPNO = spectrumDataFolderName_EXPNO, 
                   AcquPars = AcquPars))
  }
}

