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

## This script contains several functions used for the matching Shiny app.
## Author: jfmartin


################################################################################################
## set of R function for annotation of peaklists
################################################################################################


match2annot <- function (x,mzx,rtx=0,db,mzdb,rtdb=0,ioni,spectro=NULL,tol,rtwin,typeMatch)
{
   # function able to match any file containing a mz or Mass(monoIsotopic mass) and rt (not mandatory)
   # ie from xcms diffreport or CAMERA getpeaklist or MSDIAl or a list of monoistopic mass 
   # with any file contaning the same information : mz or Mass (and rt) like inhouse DB standards or HMDB or CHEBI ... 
   # the type of comparison is define with typeMatch <- mzVSmz, mzVSmass, massVSmz, massVSmass
   # mz if you compare a ratio m/z, mass if you compare a monoisotopic mass
   # for peaklist compare to DB with m/z -> typeMatch="mzVSmz"
   # for peaklist compared to HMDB (monoisotopic mass) -> typeMatch="mzVSmass"
   # if rtx=0 retention time is not used for matching 
   # with an inhouse list of mass with a tolrance of mass expressed in ppm.
   # The inhouse file must have with at least 3 columns : 
   #   ENTRY : unique identification of the metabolites
   #   NAME : name of the metabolite
   #   mz: m/z corresponding to positiveor negative ionisation
   # Parameters of the function:
   # input file to be match
   # x     :dataframe of the peaklist
   # mzx   :subscript of the m+h or m-h value in the peaklist x depending on the ionisation mode
   # rtx : :subscript of the retention time value in the peaklist x (must be the same unit as in the DB: minutes or seconds)
   
   # database to use to match
   # db    :dataframe of the inhouse metabolite db
   # mzdb  :subscript of the m+h or m-h depending on the inhouse metabolite db pos or neg
   # rtdb  :subscript of the retention time value in of the inhouse metabolite db
   # ioni  :define the ionisation mode corresponding to the peaklist "positive" or "negative"
   # tol   :define the molecular weigth tolerance expressed in ppm.
   # The function return the x dataframe with the four columns ENTRY,NAME,FORMULA and MASS 
   # If several metabolites match a peaklist mass, all annotations are concatenated separted with ";" 
   # Proton value for matching
   # It is possible to checked the retention time between peaklist and inhouse DB with a RT windows in seconds
   # If user don't want to use RT matching, he has to set a rtx=0 
   # rtwin : windows of rt to match for retention time (must be the same unit as in the peaklist and DB: minutes or seconds)
   
   massH <- 1.0078250321
   electron <- 0.0005486
   # depending on the type of match and ionisation, the value to be match. If a mz must be matched with a 
   # monoisotopic mass. mz is transformed in Mass using H and electron mass depending on the ionisation.
   if (typeMatch == 'mzVSmz' | typeMatch == 'massVSmass') deltaM <- 0 else 
      if (ioni == "positive") deltaM <- - massH + electron else deltaM <- massH - electron
   mzm <- x[[mzx]] + deltaM
   
   x <- data.frame(mzm,x, stringsAsFactors = FALSE)
   x <- x[order(x$mzm),]
   if (rtx >0) db <- db[order(db[[rtdb]]),] else db <- db[order(db[[mzdb]]),] 

   # initialisation of the db information to add in the x peaklist
   result  <- matrix("",nrow=dim(x)[1],ncol=ncol(db))
   colnames(result) <- colnames(db)
   
   for (i in 1:dim(x)[1]) {
      # matching with mz tolerance in ppm = tol using which function for each mass of the x peaklist
      # use rtention time if rtx (in the input file) > 0
      if (rtx >0) {
         mmz <- which((1e6*(abs(x[i,1]- db[[mzdb]]))/db[[mzdb]] <= tol) & 
                         (abs(x[i,rtx+1]- db[[rtdb]]) <= (rtwin+(x[i,rtx+1]/10)) | is.na(db[[rtdb]])))
      } 
      else {
         mmz <- which(1e6*(abs(x[i,1]- db[[mzdb]]))/db[[mzdb]] <= tol)
      }
      
      #if (length(mmz)>0) cat(i," mzx =",x[i,1]," RTx=",x[i,rtx+1],"\n")
      
      if (length(mmz) >0) {
         iondb <- db[mmz,]
         
         for (j in 1:ncol(iondb)) {
            result[i,j]   <- paste(iondb[,j],collapse="|")
         }
      }                    
   }
   xdb <- data.frame(x[,-1],result, stringsAsFactors = FALSE)
   return(xdb)   
}

seekIsoAdduct <- function(inName, stros="13C") {
   
   ## function tag an isotope or adduct(stros) in annotation column (inName type dataframe)
   nli <- length(inName)
   tag <- array(data=0, dim = c(nli,1))
   stag <- grep(pattern = stros, x=inName)
   
   tag[stag] <- 1 
   return(tag)
}

filterLevel1 <- function(VM) {
   
   ## function for filtering the different annotations
   ## present in at least 2/3 of total number of pools
   ## annotated with axiom database (VM$ENTRY non null) 
   ## with Retention time (column RTfil)
   ## first column contains id ions
   thrPools <- round(max(VM$pool)*2/3,digits = 0)
   
   selectL1 <- rep(0,nrow(VM))
   SelOk <- which(VM$pool>= thrPools & VM$ENTRY != "" & !is.na(VM$RT))
   selectL1[SelOk] <- 1
   VMlev1 <- cbind(VM,selectL1)
   return(VMlev1)
}

checkSpectra <- function(VM) {
   ## this function is used after spliting the multi annotation in several lines.
   ## for each compound annotated, counts the number of non duplicated ions extracted in the peaklist and
   ## compare (ratio perfannot) to the number of known expected ions in the DB spectra  
   ## 
   
   annot <- VM[VM$compound != "",]
   
   ## recherche de duplicated
   listIonUniq <- table(annot$compound)
   listIonUniq <- data.frame(listIonUniq, stringsAsFactors = FALSE)
   listIonUniq[[1]] <- as.character(listIonUniq[[1]])
   colnames(listIonUniq) <- c("compound","nobsiis")
   
   # for each compound count the number of No duplicated ions in peaklist that correspond to the DB spectra 
   nbiisNoDup <- rep(NA,nrow(listIonUniq))
   # for each compound ratio between observed ions of spectra in peaklist / expected ions in the DB spectra
   perfannot <- rep(NA,nrow(listIonUniq))
   # fror each non duplicated compound annotated an index to filter those ions 
   SPgroup <- rep(NA,nrow(listIonUniq))
   
   for (i in 1:nrow(listIonUniq)) {

      selection <- annot[annot$compound==listIonUniq[i,1],]
      ## elimination of the duplicate annot
      selection <- selection[!duplicated(selection$ENTRY),]
      #cat(listIonUniq$compound[i]," obs =",nrow(selection),"\n")
      nbiisNoDup[i] <- nrow(selection)
      SPgroup[i] <- i 
      if (selection$nbiis[1]>0) perfannot[i] <-  nbiisNoDup[i]/selection$nbiis[1]  
   }
   
   resuni <- data.frame(listIonUniq, nbiisNoDup, SPgroup, perfannot, stringsAsFactors = FALSE)
   VMnew <- merge(x = VM, y=resuni, by.x="compound", by.y="compound", all.x=TRUE, all.y = TRUE)
   iENTRY <- which(colnames(VMnew)=="ENTRY")
   VMnew <- VMnew[,c(2:iENTRY,1,(iENTRY+1):ncol(VMnew))]
   
   return(VMnew)
}

annotAxiomNew <- function(VM,msdev,ioni,annotPest=FALSE,mzx,rtx,tolannot,
                          DBCP,
                          DBdevice, matchHMDB=TRUE, matchCHEBI=TRUE,
                          annot2woRT=FALSE,splitannot=FALSE) {
   
   ## annotation of varmetadata (VM) file with :
   ##  - db axiom depending on device (msdev) and ionisation (ioni= "positive" or "negative") DBCP contaminant, DBdevice DB qtof or orbi
   ##  - HMDB 
   ##  - CHEBI
   ##  - if annotpest=TRUE match also with base contaminant(DBCP) with or without all the adduct fragment isotopes compute with Envipat 
   ##    and return input VM with annotations
   ##  - if annot2woRT=TRUE a second annotation is performed without taking into account the retention time (rtx=0) 
   
   ## subscripts of mz and db columns in the DB orbi or qtof 
   mzdbcol <- 7
   rtdbcol <- 6
   ## subscript of the 1st column to be splitted before concatenantion of annotations. Useless if no split
   indSplit <- ncol(VM)+1
   
   VM <- match2annot(x=VM,mzx,rtx,db=DBdevice,mzdb=mzdbcol,rtdb=rtdbcol,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmz')
   if (annot2woRT==TRUE) VM <- match2annot(x=VM,mzx,rtx=0,db=DBdevice,mzdb=mzdbcol,rtdb=rtdbcol,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmz')
   
   ### Pesticides new contaminant EJA DB Screening_Contaminant.xlsx match sur mz envipat
   if (annotPest==TRUE) VM <- match2annot(x=VM,mzx,rtx=0, db=DBCP,mzdb=25,rtdb = 0,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmz')
   ### Match CHEBI
   if (matchCHEBI==TRUE) VM <- match2annot(x=VM,mzx,rtx=0, db=CHEBI,mzdb=5,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmass')
   ### Match HMDB
   if (matchHMDB==TRUE) VM <- match2annot(x=VM,mzx,rtx=0, db=xHMDB,mzdb=4,ioni=ioni,tol=tolannot[1],rtwin=tolannot[2],typeMatch='mzVSmass')
   
   if (splitannot){
      beginSplit <- colnames(VM[indSplit])
      ## new subscript for DB mz and rt after concatenation of DB annotations with peaklist
      nindmz <- indSplit+mzdbcol-1
      nindrt <- indSplit+rtdbcol-1
      
      ## split multi annot
      VM <- split.multi.annot(pkl=VM, beginSplit=beginSplit,sepMatch="|")
      VM[[nindmz]] <- as.numeric(VM[[nindmz]])
      VM[[nindrt]] <- as.numeric(VM[[nindrt]])
      VM$nbiis <- as.numeric(VM$nbiis)
      
      ## computes the rt accuracy (delta in seconds between peaklist and db
      VM$rtdelta <- round(60*(abs(VM[[rtx]]-VM[[nindrt]])),0)
      ## computes the mz accuracy and sort for each ion by increasing ppm accuracy
      VM$ppmannot <- (1e6*(abs(VM[[mzx]]-VM[[nindmz]])/VM[[nindmz]]))
      ## VM <- VM[order(VM[[1]],VM$ppmannot),]
      
      ## add some information to assess annotation accuracy
      C13 <- seekIsoAdduct(VM$attribution ,"13C")
      VM <- cbind(VM,C13)
      VM <- filterLevel1(VM)
      
      ## for each annotated compound computes the number of expected ion of the spectra observed in the peaklist
      VM <- checkSpectra(VM)
      VM <- VM[order(VM$compound,VM$perfannot),]
   }
   return(VM)
}


split.multi.annot <- function(pkl, beginSplit, sepMatch="|") {
   #############################################################################################################
   ## split multi annotation separated by "sepMatch" in different lines in order to check 
   ## every possible annotation. 
   ## pkl : peaklist dataframe typically W4M variable metadata 
   ## beginSplit : colname of the column starting the split
   ## sepMatch : separator used to concatenate multiple match of the same peaklist ion
   #############################################################################################################
   
   ## split columns starting in column named beginSplit converted in subscript colbegin
   ## so indfix=columns subscripts fixed and indspl=columns to be splitted in different rows 
   colbegin <- which(colnames(pkl) == beginSplit)
   indfix <- c(1:(colbegin-1))
   indspli <- c(colbegin:ncol(pkl));
   
   ## determine the total number of matches including the "no match" defining the new number of lines of the resulting dataframe
   testN <- strsplit(as.character(pkl[,indspli[1]]),split = sepMatch , fixed = TRUE)
   maxMatch <- 0
   for (i in 1:length(testN)) {
      nM <- length(testN[[i]])
      if (nM ==0) maxMatch <- maxMatch + 1 else maxMatch <- maxMatch + nM
   }
   
   ## initialisation of the new dataframe with all matches and ...
   ## ...a delta rt and an annot ppm between peaklist and db and...
   ## ... a tag if it is a C13
   pklAllMatch <- array("",dim = c(maxMatch,(ncol(pkl)+2)))
   pklAllMatch <- data.frame(pklAllMatch, stringsAsFactors = FALSE)
   for (i in 1:ncol(pkl)) class(pklAllMatch[[i]]) <- class(pkl[[i]])
   colnames(pklAllMatch) <- colnames(pkl)
   lsn <- ncol(pkl)+1
   colnames(pklAllMatch)[lsn:(lsn+1)] <- c("rtdelta","ppmannot")
   
   cat("prep file with all matches ",maxMatch)
   ## processing : for each lines of the initial peaklist...
   nli <- 0 ## indice of lines in the pklAllMatch dataframe for test nli <- 60
   # for (i in 61:63){
   
   for (i in 1:nrow(pkl)) {
      mzl <- unlist(strsplit(as.character(pkl[i,indspli[1]]), split = sepMatch, fixed = TRUE))
      #cat(length(mzl))
      nbm <- length(mzl) 
      
      if (nbm <= 1) {
         nli <- nli+1
         #idIons[nli] <- i
         pklAllMatch[nli,c(1:ncol(pkl))] <- pkl[i,]
      } else 
      {  ## plusieurs matches pour la ligne i en cours
         inli <- nli+1
         nli <- inli+nbm-1
         pklAllMatch[c(inli:nli),indfix] <- pkl[i,indfix]
         #idIons[c(inli:nli)] <- i
         ##for (j in 1:nbm) { ##pour chaque match  
         
         for (n in indspli) {  
            ## for each column to split
            subch <- array(NA,nbm)
            ## change the sep to " "sep" " in order to split in the right number or substring
            tspl <- gsub(pattern=paste("[",sepMatch,"]",sep=""),replacement=paste(" ",sepMatch," ",sep=""),as.character(pkl[i,n]))
            subch <- unlist(strsplit( tspl, split = sepMatch , fixed = TRUE))
            
            if(length(subch)<=nbm & length(subch)!=0) {
               lastnli <- nli-(nbm-length(subch))
               pklAllMatch[c(inli:lastnli),n] <- subch
            }
         }
      }
   }
   for(c in 1:ncol(pklAllMatch)) if (is.character(pklAllMatch[[c]])) pklAllMatch[[c]] <- trimws(pklAllMatch[[c]])
   cat(" completed","\n")
   return(pklAllMatch)
}



