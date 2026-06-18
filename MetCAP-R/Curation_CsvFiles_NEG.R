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

## This script is a tool to facilitate the curation of data obtained by LC-MS/MS in NEGATIVE ion mode and processed by XCMS and CAMERA on W4M.
## It automatize the following process :
## - Suppression of signals acquired before the solvent front
## - Suppression of NaCl clusters
## - PCGROUP correction
## - Suppression of multicharged ions and associated PCGROUP

############## Upload the packages ##############
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyselect)
library(tidyverse)
library(stringr)
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

############## Suppress lines with mz with .9, .8, .7, .6 et .5 for condition RT<=1.5 ##############
data_mz_filtered <- data %>%
  filter(!(rt <= 1.5 &                                      
             str_detect(mz, "\\.9|\\.8|\\.7|\\.6|\\.5")))   

############## Suppress data with RT<0.8 ##############
data_filtered_rt <- filter(data_mz_filtered, rt>= 0.8)

############## PCGROUP correction ##############
### Create separate columns for Isotope data and the associated Code
df_filtration_isotope <- data_filtered_rt %>%
  mutate(                                                   
    Code = sub("\\[([0-9]+)\\]\\[.*", "\\1", isotopes),      
    Iso = sub("\\[[0-9]+\\]\\[(.*?)\\].*", "\\1", isotopes)  
  )

### Find the PCGROUP value for each 'Code' with Iso == "M" and replace it for M+X with divergent PCGROUP
data_filtred_pcgroup <- df_filtration_isotope %>%
  group_by(Code) %>%                                            
  mutate(
    pcgroup_updated = ifelse(Iso == "M", pcgroup, NA)  
  ) %>%
  ungroup() %>%                                         
  group_by(Code) %>%
  mutate(
    pcgroup = coalesce(pcgroup_updated[Iso == "M"][1], pcgroup)     
  ) %>%
  ungroup() %>%
  select(-pcgroup_updated)      

### Suppress lines M+X if there is no associated M with same Code
data_filtred_pcgroupB <- data_filtred_pcgroup %>%
  group_by(Code) %>%
  mutate(
    has_M = any(Iso == "M")
  ) %>%
  ungroup() %>%
  filter(
    !(Iso %in% c("M+1", "M+2", "M+3","M+4") & !has_M)
  ) %>%
  select(-has_M,-Code, -Iso)

############## Create two sub-data with RT>1.5 and with RT<1.5 ##############
filtered_data_rt_ion <- data_filtred_pcgroupB %>%
  filter(rt >= 1.5)

data_filtered_rt_oob <- data_filtred_pcgroupB %>%
  filter(rt <= 1.5)

############## Remove PCGROUP with multicharged isotopes in data with RT>1.5 with exceptions ##############
###List of exceptions
List_of_exceptions <- c("^GMP$", "^UMP$", "^AMP$", "^UDP$", "^ADP$", "^GDP$", "^Leucine$", "^Isoleucine$", "^Inosine$", "^DL-beta-Leucine$", "^Cyclic ADP-ribose$", "^dADP$", "^NADPH$", "^Nicotinamide adenine dinucleotide phosphate$", "^dGMP$", "^dUMP$", "^3'-AMP$", "^5’-AMP$", "^Cyclic AMP$", "^Uridine 5'-diphospho-N-acetylglucosamine$", "^dAMP$", "^Phenylalanine$", "^DL-Tryptophan$")

###### Step 1: Create a list of PCGROUP to not exclude (with the exceptions)
# 1) Copy column compound for modifications
filtered_data_rt_ion <- filtered_data_rt_ion %>%
  mutate(compound_clean = compound)

# 2) Suppress the parentheses (actual code bug for exceptions if there is presence of parentheses)
filtered_data_rt_ion <- filtered_data_rt_ion %>%
  mutate(compound_clean = str_replace_all(compound_clean, "\\([^)]*\\)", "")) %>%
  mutate(compound_clean = str_squish(compound_clean))

# 3) Split column compound_clean_versX
split_compound <- strsplit(filtered_data_rt_ion$compound_clean, "\\|")
max_len <- max(lengths(split_compound))

for (i in seq_len(max_len)) {
  filtered_data_rt_ion[[paste0("compound_clean_vers", i)]] <-
    sapply(split_compound, function(x) ifelse(length(x) >= i, str_trim(x[i]), NA))
}

# 4) Indicate column to investigate
cols_clean <- grep("^compound_clean_vers", names(filtered_data_rt_ion), value = TRUE)

# 5) Build a simple pattern for exceptions (technically optionnal)
pattern_exceptions <- str_c(paste0("\\b", List_of_exceptions, "\\b"), collapse = "|")

# 6) Filter the pcgroup to keep
pcgroup_to_keep <- filtered_data_rt_ion %>%
  rowwise() %>%
  filter(any(str_detect(c_across(all_of(cols_clean)), pattern_exceptions))) %>%
  ungroup() %>%
  pull(pcgroup)

###### Step 2: Create the list of PCGROUP to exclude (with 2- or 3- isotopes and without exceptions)
pcgroup_with_isotopes_to_remove <- filtered_data_rt_ion %>%
  filter(pcgroup %in% setdiff(unique(filtered_data_rt_ion$pcgroup), pcgroup_to_keep)) %>%
  filter(str_detect(isotopes, "2\\-|3\\-")) %>%
  pull(pcgroup)

###### Step 3 : Combine PCGROUP to keep and to exclude
pcgroup_to_remove <- unique(c(!pcgroup_to_keep, pcgroup_with_isotopes_to_remove))

###### Step 4 : Apply to the data
data_filtered <- filtered_data_rt_ion %>%
  filter(!(pcgroup %in% pcgroup_to_remove))  # Exclude lines for wich PCGROUP is in pcgroup_to_remove

###### Step 5: Clean the file of the new columns
### Merge compound columns
merge_compound_columns <- function(df) {
  # Find prefix of "_versX" columns
  prefixes <- unique(sub("_vers\\d+$", "", names(df)))
  # For each prefix, merge associated columns
  for (prefix in prefixes) {
    # Select all columns with same prefix
    cols_to_merge <- grep(paste0("^", prefix, "_vers"), names(df), value = TRUE)
    if (length(cols_to_merge) > 0) {
      # Merge columns for each lines with "|" as separation
      df[[prefix]] <- apply(df[, cols_to_merge, drop = FALSE], 1, function(row) {
        paste(na.omit(row), collapse = "|")
      })
      # Suppress columns with "_versX" suffix after merging
      df <- df[, !(names(df) %in% cols_to_merge)]
    }
  }
  return(df)
}

data_filtered <- merge_compound_columns(data_filtered)

###Delete compound_clean column
data_filtered <- data_filtered %>% select(-compound_clean)

############## Merge the two sub-data ##############
datamerged <- rbind(data_filtered, data_filtered_rt_oob)

############## Export the final dataframe at xlsx format ##############
write.xlsx(datamerged, file = "FileName.xlsx",          
           sheetName = "Sheet1",                               
           colnames = TRUE, rownames = TRUE, append = FALSE)   
