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

## This script is a tool to automatize the annotation propagation of cured data matched with internal database in POSITIVE ion mode.
## It automatize the following process :
## - Adding of annotation level columns
## - Suppression of HMDB data matched with non parent ions
## - Suppression of aberrant annotations
## - Automatic annotation pattern-dependant
## - Global propagation based on isotope
## - Specific RT propagation based on PCGROUP followed by adduct/isotope loop propagation
## - Specific RT propagation based on compound name followed by adduct/isotope loop propagation
## - Specific annotation on sodium adduct signals
## - Global file reorganization (relocalisation of columns and suppression of non-interest columns)

# Please define the RT to use for different annotation, and keep a security margin under the RT MAX (C18 RT = 2.5; HILIC NEG RT = 14.5; SCREENING RT = )
RTUsed = 2.5

###################################################################### Preparation of the data and loading the functions ######################################################################
#################################### Load the necessary packages for the script
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyselect)
library(tidyverse)
library(stringr)
library(data.table)
library(tidyr)

#################################### Call the different functions for the script
source("Functions_Scripts_MetCAP-R.R")


#################################### Create the dataframe from the output document of the shiny_matching
data <- read()

#################################### Convert data to data.frame (necessary for several functions)
data <- as.data.frame(data)

#################################### Formatting the file
cols <- c("mz", "mzmin", "mzmax", "rt", "rtmax", "rtmin", "adduct", "RT", "mz.1")
data[cols] <- lapply(data[cols], function(x) {          
  if (is.character(x)) {
    x <- gsub(",", ".", x)
    return(x)                   
  } else {                                       
    return(x)                                    
  }
})

data[] <- lapply(data, function(x) {          
  if (is.character(x)) {
    x <- gsub("\"", "", x)
    return(x)                   
  } else {                                       
    return(x)                                    
  }
})

#########Suppress empty rows and columns
data <- clean_empty_rows(data)
data <- remove_empty_columns_dt(data)

#################################### Adding new columns for automatic annotation
data$MSMS <- ""
data$lvl <- ""
data$annotation_confidence <- ""

#################################### Relocalisation of the new columns
data <- data %>% relocate (MSMS, lvl, annotation_confidence, .before = NULL, .after = pcgroup) 

#################################### Suppress the HMDB data matching with [M+1], [M+2] ... isotopes
data <- data %>%
  mutate(entry = if_else(grepl("M\\+[0-9]", isotopes), NA_character_, as.character(entry)),
         name.1 = if_else(grepl("M\\+[0-9]", isotopes), NA_character_, as.character(name.1)),
         formula = if_else(grepl("M\\+[0-9]", isotopes), NA_character_, as.character(formula)),
         mass = if_else(grepl("M\\+[0-9]", isotopes), NA_character_, as.character(mass)))


###################################################################### Automatic annotation ######################################################################
#################################### Split columns for comparisons necessary to propagate annotation
columns_to_split <- c("ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles" ,"attribution")

for (col in columns_to_split) {
  print(col)
  # Split the data in different columns for each sparator "|"
  data[[col]] <- strsplit(as.character(data[[col]]), "\\|")
  # Create new colulmns for each separated element
  max_length <- max(sapply(data[[col]], length))
  for (i in 1:max_length) {
    new_col_name <- paste0(col, "_part", i)  # Create a new name for each new column as "oldName_partX"
    data[[new_col_name]] <- sapply(data[[col]], function(x) ifelse(length(x) >= i, x[i], NA))
  }
  data[[col]] <- NULL
}

#################################### Suppress the attribution [M+H]+ if the isotope is [M+X]
data_corrected_motif <- remove_motif_isotope_POS(data)

#################################### Suppress the attribution 13C in lines with an [M] isotope
data_corrected_motif13C <- Remove_13C(data_corrected_motif)

#################################### Suppress the attribution different from 13CX in lines with [M+X] isotopes
data_corrected_isotopes <- Remove_no_13C(data_corrected_motif13C)

#################################### Suppress the 13CX that not correspond with the [M+X] isotopes (13C2 for isotopes [M+2], 13C1 for isotopes [M+1]...)
### It is possible to add more lines and change n value 
data_corrected_attribution <- Remove_diff_13Cn(data_corrected_isotopes, 1)
data_corrected_attribution <- Remove_diff_13Cn(data_corrected_attribution, 2)
data_corrected_attribution <- Remove_diff_13Cn(data_corrected_attribution, 3)
data_corrected_attribution <- Remove_diff_13Cn(data_corrected_attribution, 4)
#data_corrected_attribution <- Remove_diff_13Cn(data_corrected_isotopes, n)

#################################### Fill the annotation columns for lines with [M+H]+ attribution
data_cleaned_filled <- data_corrected_attribution %>%
  mutate(
    MSMS = ifelse(Verify_motif_POS(data_corrected_attribution), "CID", MSMS),
    lvl = ifelse(Verify_motif_POS(data_corrected_attribution), "2", lvl),
    annotation_confidence = ifelse(Verify_motif_POS(data_corrected_attribution), "a,b,c", annotation_confidence)
  )

#################################### Fill the annotation columns for lines with [M] isotopes and [M+H]+ attribution
data_cleaned_filled <- data_cleaned_filled %>%
  mutate(
    MSMS = ifelse(verify_isotope_POS(data_cleaned_filled) & Verify_motif_POS(data_cleaned_filled), "CID", MSMS),
    lvl = ifelse(verify_isotope_POS(data_cleaned_filled) & Verify_motif_POS(data_cleaned_filled), "2", lvl),
    annotation_confidence = ifelse(verify_isotope_POS(data_cleaned_filled) & Verify_motif_POS(data_cleaned_filled), "a,b,c", annotation_confidence)
  )

#################################### Merge the columns previously splited
data_cleaned_merged <- merge_part_columns(data_cleaned_filled)

#################################### Propagate annotation to "superior" isotopes [M+(X+1)]
data_cleaned_merged <- copy_compound_formula(data_cleaned_merged)

#################################### Create sub data for the different retention time
data_sup_2.5 <- data_cleaned_merged %>%
  filter(rt> RTUsed)

data_inf_2.5 <- data_cleaned_merged %>%
  filter(rt<= RTUsed)

###################################################################### Annotation propagation for RT>2.5 min ######################################################################
#################################### Propagate annotation to the other members of the PCGROUP under the condition there is no other annotated compound in this PCGROUP
data_sup_2.5 <- propagate_info(data_sup_2.5)

data_sup_2.5$motif_found_bis <- NULL

#################################### Propagate annotation to lines with the same compound name
#################################### Split columns for comparisons necessary to propagate annotation
for (col in columns_to_split) {
  print(col) 
  data_sup_2.5[[col]] <- strsplit(as.character(data_sup_2.5[[col]]), "\\|")
  max_length <- max(sapply(data_sup_2.5[[col]], length))
  for (i in 1:max_length) {
    new_col_name <- paste0(col, "_part", i)
    data_sup_2.5[[new_col_name]] <- sapply(data_sup_2.5[[col]], function(x) ifelse(length(x) >= i, x[i], NA))
  }
  data_sup_2.5[[col]] <- NULL
}

#################################### Verify empty cells are NA 
data_sup_2.5 <- data_sup_2.5 %>%
  mutate(across(c(MSMS, lvl, annotation_confidence), ~na_if(.x, "")))

data_sup_2.5 <- propagation_compound(data_sup_2.5)

#################################### Merge columns
data_sup_2.5 <- merge_part_columns(data_sup_2.5)

#################################### Split of adduct columns
columns_to_splitB <- c("adduct")
for (col in columns_to_splitB) {
  print(col)
  # Split data of the column with space separator
  data_sup_2.5[[col]] <- strsplit(as.character(data_sup_2.5[[col]]), split = " ")
  max_length <- max(sapply(data_sup_2.5[[col]], length)) 
  for (i in 1:max_length) {
    new_col_name <- paste0(col, "_part", i) 
    data_sup_2.5[[new_col_name]] <- sapply(data_sup_2.5[[col]], function(x) ifelse(length(x) >= i, x[i], NA))
  }
  data_sup_2.5[[col]] <- NULL
}

#################################### Propagate annotation to the adducts in the same PCGROUP
data_sup_2.5 <- propagate_adduct_pcgroup(data_sup_2.5)

#################################### Propagate annotation to "superior" isotopes [M+(X+1)]
data_sup_2.5 <- copy_compound_formula_B(data_sup_2.5)

#################################### Annotation propagation loop 
fonctions <- list(propagate_adduct_B, copy_compound_formula_B)

repeat {
  data_sup_2.5_repeated <- data_sup_2.5
  
  for (f in fonctions) {
    data_sup_2.5 <- f(data_sup_2.5)
  }
  
  if (identical(data_sup_2.5, data_sup_2.5_repeated)) {
    break
  }
}

#################################### Adduct columns merge
data_sup_2.5 <- merge_part_columns_B(data_sup_2.5)


###################################################################### Annotation propagation for RT<2.5 min ######################################################################
#################################### Split columns for comparisons necessary to propagate annotation
for (col in columns_to_split) {
  print(col) 
  data_inf_2.5[[col]] <- strsplit(as.character(data_inf_2.5[[col]]), "\\|")
  max_length <- max(sapply(data_inf_2.5[[col]], length))
  for (i in 1:max_length) {
    new_col_name <- paste0(col, "_part", i)  # Créer un nom pour chaque nouvelle colonne
    data_inf_2.5[[new_col_name]] <- sapply(data_inf_2.5[[col]], function(x) ifelse(length(x) >= i, x[i], NA))
  }
  data_inf_2.5[[col]] <- NULL
}

#################################### Suppress data not linked to attribution [M+H]+ on annotate lines
data_inf_2.5 <- Filter_attribution_POS(data_inf_2.5)

#################################### Propagate annotation to lines with the same compound name
#################################### Verify empty cells are NA
data_inf_2.5 <- data_inf_2.5 %>%
  mutate(across(c(MSMS, lvl, annotation_confidence), ~na_if(.x, "")))

data_inf_2.5 <- propagation_compound(data_inf_2.5)

#################################### Merge the columns previously splited
data_inf_2.5 <- merge_part_columns(data_inf_2.5)

#################################### Split of adduct columns
columns_to_splitB <- c("adduct")
for (col in columns_to_splitB) {
  print(col) 
  data_inf_2.5[[col]] <- strsplit(as.character(data_inf_2.5[[col]]), split = " ")
  max_length <- max(sapply(data_inf_2.5[[col]], length)) 
  for (i in 1:max_length) {
    new_col_name <- paste0(col, "_part", i) 
    data_inf_2.5[[new_col_name]] <- sapply(data_inf_2.5[[col]], function(x) ifelse(length(x) >= i, x[i], NA))
  }
  data_inf_2.5[[col]] <- NULL
}

############################################## propagation adducts 
data_inf_2.5 <- propagate_adduct(data_inf_2.5)

#################################### Propagate annotation to "superior" isotopes [M+(X+1)]
data_inf_2.5 <- copy_compound_formula_B(data_inf_2.5)

#################################### Annotation propagation loop
fonctions <- list(propagate_adduct, copy_compound_formula_B)

repeat {
  data_inf_2.5_repeated <- data_inf_2.5
  
  for (f in fonctions) {
    data_inf_2.5 <- f(data_inf_2.5)
  }
  
  if (identical(data_inf_2.5, data_inf_2.5_repeated)) {
    break
  }
}

#################################### Adduct columns merge
data_inf_2.5 <- merge_part_columns_B(data_inf_2.5)

#################################### Relocalisation of the columns
data_inf_2.5$PCGROUP_Conflict <- ""
data_inf_2.5 <- data_inf_2.5 %>% relocate (PCGROUP_Conflict, .before = NULL, .after = pcgroup)
data_sup_2.5 <- data_sup_2.5 %>% relocate (PCGROUP_Conflict, .before = NULL, .after = pcgroup)
data_sup_2.5 <- data_sup_2.5 %>% relocate (adduct, .before = NULL, .after = isotopes)
data_sup_2.5 <- data_sup_2.5 %>% relocate (ENTRY, compound, RT, mz.1, Formula, Subclass, CHEBI, Inchi, InchiKey, Smiles, attribution, .before = NULL, .after = annotation_confidence) 
data_inf_2.5 <- data_inf_2.5 %>% relocate (adduct, .before = NULL, .after = isotopes)
data_inf_2.5 <- data_inf_2.5 %>% relocate (ENTRY, compound, RT, mz.1, Formula, Subclass, CHEBI, Inchi, InchiKey, Smiles, attribution, .before = NULL, .after = annotation_confidence)

#################################### Rename columns for correspondace between the sub data
colnames(data_sup_2.5) <- colnames(data_inf_2.5)
#################################### Merge sub data for RT>2.5 min and RT<2.5 min
data_annoted <- rbind(data_sup_2.5, data_inf_2.5)    

#################################### Annotate with CIDm/z-22 if already annotate and attribution is [M+Na]+ 
data_annoted <- data_annoted %>%
  mutate(
    MSMS = ifelse(Verify_motif_Na(data_annoted), "CIDm/z-22", MSMS))

#################################### Add 4 in lvl column if empty
data_annoted <- data_annoted %>%
  mutate(lvl = ifelse(lvl == "" | is.na(lvl), 4, lvl))

#################################### Remove non-necessary colulns and reorganize the file
data_annoted <- data_annoted %>% select(-ENTRY,-mz.1,-Subclass,-CHEBI,-Formula,-InchiKey,-Smiles,-formula,-mass)
data_annoted <- data_annoted %>% relocate (entry, name.1, .before = NULL, .after = attribution)

#################################### Export the final file at excel format
write.xlsx(data_annoted, file = "data_annoted.xlsx", #Name the file with the extension
           sheetName = "Sheet1", #Name the sheet
           colnames = TRUE, rownames = TRUE, append = FALSE)
