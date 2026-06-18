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

## This script contains several functions used for the database creation and matching Shiny app.
## Author: jfmartin

###  Tools dedicated to DB creation for matching with peaklist using package enviPat

checkFormulas <- function(formulas, adduct=NULL){
   # test formulas
   if(length(formulas) == 0) stop('no formulas given')
   pass <- data.frame(formula=formulas, pass=TRUE, mess=NA, stringsAsFactors=FALSE)
   test <- grepl('[\\+-]', pass$formula)
   if(any(test)) pass[which(pass[, 2]), ][which(test), c(2, 3)] <- list(
      rep(FALSE, times=length(which(test))), rep('no sign in formulas!', times=length(which(test))))
   test <- grepl('\\[[[:digit:]]\\]', pass[which(pass[, 2]), 'formula'])
   if(any(test)) pass[which(pass[, 2]), ][which(test), c(2, 3)] <- list(
      rep(FALSE, times=length(which(test))), rep('no isotopes in formulas!', times=length(which(test))))
   
   res <- check_chemform(isotopes, pass[which(pass[, 2]), 'formula'])
   if(any(res$warning)) pass[which(pass[, 2]), ][which(res$warning), c(2, 3)] <- list(
      rep(FALSE, times=length(which(res$warning))), rep('no isotopes in formulas!', times=length(which(res$warning))))
   pass$new_formula <- NA
   pass[which(pass[, 2]), 'new_formula'] <- res[which(!res$warning), 'new_formula']
   
   # make new formulas with adduct's formula (if it's not null)
   if(!is.null(adduct)){
      if(!adduct %in% adducts$Name) stop('incorrect adduct')
      adduct <- adducts[which(adducts$Name == adduct), ]
      if(adduct$Multi > 1) pass[which(pass[, 2]), 'new_formula'] <- sapply(
         pass[which(pass[, 2]), 'new_formula'], function(formula) mergeform(formula, adduct$Multi))
      if(adduct$Formula_add != 'FALSE') pass[which(pass[, 2]), 'new_formula'] <- sapply(
         pass[which(pass[, 2]), 'new_formula'], function(formula) mergeform(formula, adduct$Formula_add))
      if(adduct$Formula_ded != 'FALSE'){
         test <- sapply(pass[which(pass[, 2]), 'new_formula'], function(formula) check_ded(formula, adduct$Formula_ded)) 
         if(any(test == "TRUE")) pass[which(pass[, 2]), ][which(test == "TRUE"), c(2, 3)] <- list(
            rep(FALSE, times=length(which(test == "TRUE"))), 
            rep('cannot substract adduct formula', times=length(which(test == "TRUE"))))
         pass[which(pass[, 2]), 'new_formula'] <- sapply(pass[which(pass[, 2]), 'new_formula'], 
                                                         function(formula) subform(formula, adduct$Formula_ded))
      }
   }
   pass
}

getIsotopicPattern <- function(formulas, adduct=NULL, charge=0, instrument=NULL, resolution=NULL, 
                               nknots=6, spar=0.2, threshold=0.1, emass=0.00054858, algo=2, ppm=FALSE, dmz = "get", frac = 1/4, 
                               env = "Gaussian", detect = "centroid", plotit = FALSE, verbose = FALSE){
   tryCatch({
      if(!is.numeric(charge)) stop('charge is not numeric')
      
      checked <- checkFormulas(formulas, adduct)
      if(any(!checked$pass)){
         print('################## ERROR ####################')
         print(checked[which(!checked$pass), c('formula', 'mess')])
         return(NULL)
      }
      else checked <- check_chemform(isotopes, checked$new_formula)
      
      if(!is.null(adduct)) charge <- charge + adducts[which(adducts$Name == adduct), 'Charge']
      if(!is.null(instrument) & !is.null(resolution)) stop('choose between resolution or instrument, not both')
      else if(is.null(instrument) & is.null(resolution)) return(
         isopattern(isotopes, checked[, 2], threshold=threshold, charge=charge, emass=emass, 
                    plotit=plotit, algo=algo, verbose=verbose))
      else if(!is.null(instrument)){
         if(!instrument %in% names(resolution_list)) stop('instrument incorrect')
         resmassIndex <- which(names(resolution_list) == instrument)
         return(isowrap(isotopes, checked, resmass=resolution_list[[resmassIndex]], 
                        charge=charge, nknots=nknots, spar=spar, threshold=threshold, emass=emass,
                        algo=algo, ppm=FALSE, dmz=dmz, frac=frac, env=env, detect=detect, plotit=plotit, verbose=verbose))
      }
      else if(!is.null(resolution)){
         if(!is.numeric(resolution)) stop('resolution is not numeric')
         else if(resolution <= 0) stop('resolution must be > 0')
         return(isowrap(isotopes, checked, resolution=resolution, charge=charge, 
                        nknots=nknots, spar=spar, threshold=threshold, emass=emass,
                        algo=algo, ppm=FALSE, dmz=dmz, frac=frac, env=env, detect=detect, plotit=plotit, verbose=verbose))
      }
   }, error=function(e){
      print('########## ERROR ################')
      cat(e$message, sep='\n', '\n')
      return(list())
   })
}

isowrap <- function(isotopes, checked, nknots, spar, threshold, charge, emass, 
                    algo, ppm, dmz, frac, env, detect, plotit, verbose, resmass=FALSE, resolution=FALSE){
   if (length(resmass) > 1){
      resolution <- getR(checked, resmass = resmass, nknots = nknots, 
                         spar = spar, plotit = plotit)
   }
   pattern <- isopattern(isotopes, checked[, 2], threshold = threshold, 
                         charge = charge, emass = emass, plotit = plotit, algo = algo, verbose = verbose)
   profiles <- envelope(pattern, ppm = ppm, dmz = dmz, frac = frac, 
                        env = env, resolution = resolution, plotit = plotit, verbose = verbose)
   centro <- vdetect(profiles, detect = detect, plotit = plotit, verbose = verbose)
   return(centro)
}






creaDBtoMatch <- function(iniDBl, colFormula, colID, ioni, adduct, labIso, spectro, reso) {
   ## function to create a database of mass + adduct + Isoptopes + relative intenisty
   ## based on enviPat package
   ## input : 
   #  - iniDBl : dataframe with initial list of mass with at least
   #  - ID in column N� colID
   #  - Formula in column N� colFormula
   #  - ioni : type of ionisation "pos" or "neg"
   #  - adduct vector character with the list of adduct to compute :c('M+H','M+Na')
   #  - labIso : vector of labels of isotopes :  c("[M]","[M+1]","[M+2]","[M+3]","[M+4]","[M+5]","[M+6]","[M+7]","[M+8]","[M+9]","[M+10]",...
   #  - spectro : enviPat name of spectro in resolution_list
   #  - reso : resolution value if your spectro is unknown in resolution_list
   ## output : dataframe iniDBl + adduct + isotopes labels 
   
   
   nbi <- nrow(iniDBl)
   nbad <- length(adduct)
   
   ## Blank detection and replacement by "" and lower to upper transform for C O N H S
   iniDBl[,colFormula] <- gsub(" ", "", as.character(iniDBl[,colFormula]), fixed = TRUE)
   iniDBl[,colFormula] <- trim(as.character(iniDBl[,colFormula]))
   
   ## Replace lowercase by uppercase
   iniDBl[,colFormula] <- gsub("c", "C", as.character(iniDBl[,colFormula]), fixed = TRUE)
   iniDBl[,colFormula] <- gsub("o", "O", as.character(iniDBl[,colFormula]), fixed = TRUE)
   iniDBl[,colFormula] <- gsub("n", "N", as.character(iniDBl[,colFormula]), fixed = TRUE)
   iniDBl[,colFormula] <- gsub("h", "H", as.character(iniDBl[,colFormula]), fixed = TRUE)
   iniDBl[,colFormula] <- gsub("s", "S", as.character(iniDBl[,colFormula]), fixed = TRUE)
   
   ## Detection of "+" in formula, suppress "+" and substract 1 H
   chk1 <- grepl("\\+",iniDBl[,colFormula])
   if (sum(chk1)>0) {
      iniDBl[chk1,colFormula] <- gsub("+", "", as.character(iniDBl[chk1,colFormula]), fixed = TRUE)
      iniDBl[chk1,colFormula] <- subform(iniDBl[chk1,colFormula],"H1")
   }
   
   ## Detection of "-" in formula and substract H
   chk2 <- grepl("\\-",iniDBl[,colFormula])
   if (sum(chk2)>0) {
      iniDBl[chk2,colFormula] <- gsub("-", "", as.character(iniDBl[chk2,colFormula]), fixed = TRUE)
      iniDBl[chk2,colFormula] <- mergeform(iniDBl[chk2,colFormula],"H1")
   }
   
   ## Calcul isotopologues with theoretical relative abundance
   ###################################################################################

   for (i in 1:nbi)  {
      ## for each ion 
      currentForm <- iniDBl[i,colFormula]
      
      cat(i,currentForm,"\n")
      currentId <- iniDBl[i,colID]

      
      for (j in 1:nbad) {
         ## for each adduct, NbIso contains the number or isotopes computed. Could be different for the 
         ## molecules depending on the atoms and adducts
         res <- getIsotopicPattern(currentForm,adduct[j],charge=0, instrument=spectro, resolution=reso)
         if (!is.null(res)){
            NbIso <- dim(res[[1]])[1]
            currentIso <- labIso[1:NbIso]
            #cat(" Iso",NbIso,"\n")         
            currentRes <- data.frame(currentId,currentForm,adduct[j],currentIso,res[[1]])
            if (j==1) currentAdductRes <- currentRes else 
               currentAdductRes <- rbind(currentAdductRes,currentRes)
         }
      }
      
      if (i==1) resAll <- currentAdductRes else resAll <- rbind(resAll,currentAdductRes)
   }
   
   final <- merge(x=iniDBl, y=resAll, by.x=colID, by.y=1)
   final <- data.frame(ioni,final)
   
   return(final)
}

dateVersion <- function() {
   ## return a character string for versionning database axiom with the format _YYMMDD
   hoy <- Sys.time()
   YY <- as.character(substr(hoy,3,4))
   MM <- as.character(substr(hoy,6,7))
   DD <- as.character(substr(hoy,9,10))
   verlabel <- paste("_",YY,MM,DD,sep ="")
   return(verlabel)
}

