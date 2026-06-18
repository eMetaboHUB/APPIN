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

## This script allows the creation of a database usable for the matching shiny app.
## Use a database at the same format as the one provided in the Git project

source("toolbox.R")
source("DB_box.R")

library(openxlsx)
###############################################################
INdb <- "Database_Test.xlsx"  #### Indicate the name of your Database at format corresponding to the template provided in the Git project
dateUp <- dateVersion()        #### Allows the update of the date version       

std <- read.xlsx(xlsxFile= INdb) 

 std$monoistopic_molecular_weigth <- as.numeric(std$monoistopic_molecular_weigth)
 std$rt <- as.numeric(std$rt)
 std$mz <- as.numeric(std$mz)

std <- std[std$spectro =="QTOF-C18",] #### Define the spectrometer used corresponding to the one indicated in your Database


#################################### P O S I T I V E  I O N  M O D E #######################################


DBQTOFPOS <- std[(std$ionisation == "pos" & !is.na(std$mz)),c(1,4,10,11,8,29,18,14,15,16,5)]
colnames(DBQTOFPOS) <- c("ENTRY","Name","RT","mz","Formula","Subclass","CHEBI","Inchi","InchiKey","Smiles","attribution")   ####Name the columns


#################################### N E G A T I V E  I O N  M O D E #######################################

DBQTOFNEG <- std[(std$ionisation == "neg" & !is.na(std$mz)),c(1,4,10,11,8,29,18,14,15,16,5)]
colnames(DBQTOFNEG) <- c("ENTRY","Name","RT","mz","Formula","Subclass","CHEBI","Inchi","InchiKey","Smiles","attribution")

setwd(paste(repPar,"DB/",sep=""))
filposname <- paste("DBQTOFPOS",dateUp,sep="") #### Change the name of the DB in function of the spectrometer / method
filnegname <- paste("DBQTOFNEG",dateUp,sep="") #### Change the name of the DB in function of the spectrometer / method
save(DBQTOFPOS, file=paste(filposname,".Rdata",sep=""))
save(DBQTOFNEG, file=paste(filnegname,".Rdata",sep=""))
## write.table(DBQTOFPOS,file=paste(filposname,".txt",sep=""),sep="\t", row.names=F,quote=F)
## write.table(DBQTOFNEG,file=paste(filnegname,".txt",sep=""),sep="\t", row.names=F,quote=F)



