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

## This file contains all functions for the use of the "Data_Curation" POS & NEG Scripts and the "Automatic_Annotation" POS & NEG Scripts

##############################Version 1.0 - Compilation of the functions from the Automatic annotation POS and NEG
##############################Version 1.1 - Translation of all indications
##############################Version 1.2 - Adding propagation at RT > 2.5 by corresponding compound name
##############################Version 1.2.2 - Modification of the different scores and adding a possibility to modify SD in the score function samplemean/blankemean
##############################Version 1.2.3 - adding detection in POS annotation of pattern ([M])+ in addition to [M+H]+
##############################Version 1.2.4 - Creating a succession of linear functions to score the sample_value
##############################Version 1.3 - Modifying the score on sample_value with a new formula
##############################Version 1.3.2 - Modifying the score on sample_value with a new formula
##############################Version 1.3.3 - Modification of Remove_13C to avoid conflicts with HMDB name.1 conlumn
##############################Version 1.3.4 - Duplication of propagate_adduct to add a RT condition for loop at RT >2.5


###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################

#################################### READ function to open file
read <- function() {
  file_path <- file.choose()     # Open the window to choose the file
  
      # Load the necessary packages to read the different formats
  require(readr)
  require(readxl)
  require(haven)
  require(jsonlite)
  require(data.table)
  
  ext <- tolower(tools::file_ext(file_path))     # Extract the extension of the file
  
      # Read the file depending of the format
  data <- switch(ext,
                 csv = readr::read_csv(file_path),
                 tsv = readr::read_tsv(file_path),
                 tab = readr::read_tsv(file_path),
                 tabular = readr::read_tsv(file_path),
                 txt = data.table::fread(file_path),
                 xls = readxl::read_excel(file_path),
                 xlsx = readxl::read_excel(file_path),
                 sav = haven::read_sav(file_path),
                 dta = haven::read_dta(file_path),
                 sas7bdat = haven::read_sas(file_path),
                 json = jsonlite::fromJSON(file_path),
                 stop(paste("Format non supporté :", ext))
  )
  return(data)
}

#################################### Clean row function
clean_empty_rows <- function(df) {
      # Convert in data.frame
  if (!is.data.frame(df)) {
    df <- tryCatch(as.data.frame(df), error = function(e) return(df))
  }
      # Replace NA by empty cells
  df[is.na(df)] <- ""
      # Suppress lines where cells are all empty or blanks
  df <- df[!apply(df, 1, function(row) all(trimws(as.character(row)) == "")), , drop = FALSE]
  return(df)
}

#################################### Clean column function
remove_empty_columns_dt <- function(dt) {
      # Test column by column
  empty_cols <- sapply(dt, function(col) {
    all(is.na(col) | (col == "") | (trimws(col) == ""))
  })
  cols_to_keep <- names(dt)[!empty_cols]
  dt[, cols_to_keep]
}

#################################### Function to verify pattern in POS dataframe in attribution column
Verify_motif_POS <- function(df) {
      # Select columns beginning by "attribution_part"
  attribution_cols <- grep("^attribution_part", names(df), value = TRUE)
      # Verify if the pattern "[M+H]+" is present in any "attribution" column
  motif_found <- apply(df[, attribution_cols], 1, function(row) {
    any(grepl("\\[M\\+H\\]\\+$", row) | grepl("\\[\\(M\\)\\]\\+$", row))      # Verify if "[M+H]+" or "([M])+" is in any value of the line
  })
      # Return logical vector (TRUE if the pattern is found in at least one column, FALSE else)
  return(motif_found)
}

#################################### function to verify pattern in NEG dataframe in attribution column
Verify_motif_NEG <- function(df) {
      # Select columns beginning by "attribution_part"
  attribution_cols <- grep("^attribution_part", names(df), value = TRUE)
      # Verify if the pattern "[M-H]-" is present in any "attribution" column
  motif_found <- apply(df[, attribution_cols], 1, function(row) {
    any(grepl("\\[M\\-H\\]\\-$", row))    # Verify if "[M-H]-" is in any value of the line
  })
      # Return logical vector (TRUE if the pattern is found in at least one column, FALSE else)
  return(motif_found)
}

#################################### function to find code value and isotope in the column "isotopes"
extract_isotope_info <- function(isotope_str) {
        # Use a regular expression to extract the code number and the isotope number if present
  match <- regexec("\\[(\\d+)\\]\\[M(\\+(\\d+))?\\]", isotope_str)
  result <- regmatches(isotope_str, match)
  if (length(result[[1]]) > 0) {
        # Extract the code number (always in position 2)
    number <- as.numeric(result[[1]][2])
        # If the isotope number is not present, X_value is defined as 0
    X_value <- ifelse(length(result[[1]]) > 3 && !is.na(result[[1]][4]), 
                      as.numeric(result[[1]][4]), 0)
        # If the string is in type of [Number][M] (without +X), define explicitly X_value =0
    if (grepl("\\[M\\]", isotope_str) && !grepl("\\+\\d+", isotope_str)) {
      X_value <- 0
    }
    cat("Extracted Number: ", number, ", Extracted X_value: ", X_value, "\n")
    return(c(number, X_value))
  } else {
    return(c(NA, NA))  # Return NA if no corresponding value
  }
}

#################################### function to suppress [M+H]+ in attribution column for isotopes M+X
remove_motif_isotope_POS <- function(df) {
                # Apply pattern verification on the data.frame
  motif_found <- Verify_motif_POS(df)
                # Apply on the lines where the pattern is found
  for (i in 1:nrow(df)) {
    if (motif_found[i]) {
                # Find the "attribution_partX" columns with the pattern for this line
      attribution_cols <- grep("^attribution_part", names(df), value = TRUE)
      matching_cols <- attribution_cols[sapply(attribution_cols, function(col) {
        grepl("\\[M\\+H\\]\\+$", df[i, col]) | grepl("\\[\\(M\\)\\]\\+$", df[i, col])
      })]
                # For each column with the pattern, verify the corresponding isotope
      for (col in matching_cols) {
                # Extract the isotope information for the source line
        source_isotope <- extract_isotope_info(df[i, "isotopes"])
        source_number <- source_isotope[1]
        source_X <- source_isotope[2]
                # Verify if X_value is valid (no NA) and superior to 0
        if (!is.na(source_X) && source_X > 0) {
          suffix <- sub("^attribution_part", "", col)        # Extract the suffix "X" from the "attribution_partX" column
                # Identify columns to suppress (which don't have the same _partX suffix)
          cols_to_remove <- grep(paste0("_part", suffix), names(df))
          # Suppress these columns (put as NA)
          df[i, cols_to_remove] <- NA
        }
      }
    }
  }
  return(df)
}

#################################### function to suppress [M-H]- in attribution column for isotopes M+X 
remove_motif_isotope_NEG <- function(df) {
              # Apply pattern verification on the data.frame
  motif_found <- Verify_motif_NEG(df)
              # Apply on the lines where the pattern is found
  for (i in 1:nrow(df)) {
    if (motif_found[i]) {
              # Find the "attribution_partX" columns with the pattern for this line
      attribution_cols <- grep("^attribution_part", names(df), value = TRUE)
      matching_cols <- attribution_cols[sapply(attribution_cols, function(col) {
        grepl("\\[M\\-H\\]\\-$", df[i, col])
      })]
              # For each column with the pattern, verify the corresponding isotope
      for (col in matching_cols) {
              # Extract the isotope information for the source line
        source_isotope <- extract_isotope_info(df[i, "isotopes"]) 
        source_number <- source_isotope[1]
        source_X <- source_isotope[2]
              # Verify if X_value is valid (no NA) and superior to 0
        if (!is.na(source_X) && source_X > 0) {
          suffix <- sub("^attribution_part", "", col)        # Extract the suffix "X" from the "attribution_partX" column
              # Identify columns to suppress (which don't have the same _partX suffix)
          cols_to_remove <- grep(paste0("_part", suffix), names(df))
              # Suppress these columns (put as NA)
          df[i, cols_to_remove] <- NA
        }
      }
    }
  }
  return(df)
}

#################################### function to suppress 13C in lines with isotopes M in attribution column
Remove_13C <- function(df) {
        # Browse all colums starting by "attribution_part"
  for (col in grep("_part", names(df), value = TRUE)) {
        # Select lines where "isotope" column contain "[M]" and "attribution" column contain "13C"
    lignes_motif_M_et_13C <- grepl("\\[M\\]", df$isotopes) & grepl("13C", df[[col]])
        # If pattern "13C" is found in a column, suppresse informations in this columns and all orthers with same "_partX"
    if (any(lignes_motif_M_et_13C)) {
      suffixe_partX <- sub("^attribution_part", "", col)      # Extract the suffix "X" from colum "attribution_partX"
      colonnes_a_supprimer <- grep(paste0(".*_part", suffixe_partX, "$"), names(df), value = TRUE)
        # Replace values by NA in all corresponding columns
      df[lignes_motif_M_et_13C, colonnes_a_supprimer] <- NA
    }
  }
  return(df)
}

#################################### function to suppress non 13C in lines with isotopes M+X in attribution column
Remove_no_13C <- function(df) {
        # Identify lines where pattern [M(+n)] is present without being [M]
  lines_iso <- grepl("\\[(\\d+)\\]\\[M(\\+(\\d+))?\\]", df$isotopes) & !grepl("\\[M\\]", df$isotopes)
  for (col in grep("^attribution_part", names(df), value = TRUE)) {
        # Verify if the column contain "13C" for lines corresponding to lines_iso
    no13C_in_col <- !grepl("13C", df[[col]])
        # Replace by NA in column for lines where condition is fulfilled
    df[[col]][lines_iso & no13C_in_col] <- NA
        # Extract the suffix X from the column "attribution_part"
    suffixe <- sub("^attribution_part", "", col)
        # Find columns associated to this suffix
    colonnes_a_supprimer <- grep(paste0("_part", suffixe), names(df), value = TRUE)
        # Replace values by NA in all corresponding columns
    for (assoc_col in colonnes_a_supprimer) {
      df[[assoc_col]][lines_iso & no13C_in_col] <- NA
    }
  }
  return(df)
}
 
#################################### function to suppress non 13CX in lines with isotopes M+Y in attribution column
Remove_diff_13Cn <- function(df, n) {
        # Identify lines where pattern [M(+n)] is present
  lines <- grepl(paste0("\\[M\\+", n, "\\]"), df$isotopes)
        # Iterate on all columns "attribution_partX"
  for (col in grep("^attribution_part", names(df), value = TRUE)) {
        # Verify if the column contain "13Cn" for corresponding lines
    no13C_in_col <- !grepl(paste0("13C", n), df[[col]])
        # Replace by NA in columns where condition is fulfilled
    df[[col]][lines & no13C_in_col] <- NA
        # Extract the suffix X from the column "attribution_part"
    suffixe <- sub("^attribution_part", "", col)
        # Find columns associated to this suffix
    colonnes_a_supprimer <- grep(paste0("_part", suffixe), names(df), value = TRUE)
        # Replace values by NA in all corresponding columns
    for (assoc_col in colonnes_a_supprimer) {
      df[[assoc_col]][lines & no13C_in_col] <- NA
    }
  }
  return(df)
}

#################################### function to verify if there is [M]+ in column isotopes
verify_isotope_POS <- function(df) {
  isotope_col <- grep("isotopes", names(df), value = TRUE)
  if (length(isotope_col) == 0) {
    warning("Aucune colonne contenant 'isotopes' trouvée.")
    return(rep(FALSE, nrow(df)))  # retourne un vecteur FALSE pour chaque ligne
  }
  isotope_motif <- apply(df[, isotope_col, drop = FALSE], 1, function(row) {
    any(grepl("\\[M\\]\\+", row))
  })
  return(isotope_motif)
}

#################################### function to verify if there is [M]- in column isotopes
verify_isotope_NEG <- function(df) {
  isotope_col <- grep("isotopes", names(df), value = TRUE)
  if (length(isotope_col) == 0) {
    warning("Aucune colonne contenant 'isotopes' trouvée.")
    return(rep(FALSE, nrow(df)))  # retourne un vecteur FALSE pour chaque ligne
  }
  isotope_motif <- apply(df[, isotope_col, drop = FALSE], 1, function(row) {
    any(grepl("\\[M\\]\\-", row))
  })
  return(isotope_motif)
}

#################################### function to merge column splited previously (_partX)
merge_part_columns <- function(df) {
          # Find prefix of "_partX" columns
  prefixes <- unique(sub("_part\\d+$", "", names(df)))
          # For each prefix, merge associated columns
  for (prefix in prefixes) {
          # Select all columns with same prefix
    cols_to_merge <- grep(paste0("^", prefix, "_part"), names(df), value = TRUE)
    if (length(cols_to_merge) > 0) {
          # Merge columns for each lines with "|" as separation
      df[[prefix]] <- apply(df[, cols_to_merge, drop = FALSE], 1, function(row) {
        paste(na.omit(row), collapse = "|")
      })
          # Suppress columns with "_partX" suffix after merging
      df <- df[, !(names(df) %in% cols_to_merge)]
    }
  }
  return(df)
}

#################################### function to propagate data to same isotopes
copy_compound_formula <- function(df) {
                  # Create a vector to indicate lines with MSMS, lvl and annotation confidence filled
  filled_rows <- which(!is.na(df$MSMS) & !is.na(df$lvl) & !is.na(df$annotation_confidence))
  for (source_row in filled_rows) {
                  # Extract the isotope information for the source line
    source_isotope <- extract_isotope_info(df$isotopes[source_row])
    source_number <- source_isotope[1]
    source_X <- source_isotope[2]
                  # If isotope informations are valid
    if (!is.na(source_number) && !is.na(source_X)) {
                  # Search for the lines which have the same code number and an isotope number superior
      target_rows <- which(df$isotopes != "" & !is.na(df$isotopes))       # Ignore empty lines
      for (target_row in target_rows) {
        target_isotope <- extract_isotope_info(df$isotopes[target_row])
        target_number <- target_isotope[1]
        target_X <- target_isotope[2]
                  # Verify if the code number is the same and if the isotope number is superior to source line
        if (!is.na(target_number) && !is.na(target_X) && target_number == source_number && target_X > source_X) {
                  # Copy of information from source line to target line
          df$lvl[target_row] <- df$lvl[source_row]
          df$annotation_confidence[target_row] <- df$annotation_confidence[source_row]
                  # Copy of the other information if the attribution column is empty
          if (is.na(df$attribution[target_row]) || df$attribution[target_row] == "") {
            df$ENTRY[target_row] <- df$ENTRY[source_row]
            df$compound[target_row] <- df$compound[source_row]
            df$RT[target_row] <- df$RT[source_row]
            df$mz.1[target_row] <- df$mz.1[source_row]
            df$Formula[target_row] <- df$Formula[source_row]
            df$Subclass[target_row] <- df$Subclass[source_row]
            df$CHEBI[target_row] <- df$CHEBI[source_row]
            df$Inchi[target_row] <- df$Inchi[source_row]
            df$InchiKey[target_row] <- df$InchiKey[source_row]
            df$Smiles[target_row] <- df$Smiles[source_row]
          }
        }
      }
    }
  }
  return(df)
}

#################################### function to verify the presence of automatic annotation (MSMS column)
Verify_motif_bis <- function(df) {
        # Find column containing "MSMS"
  attribution_cols <- grep("MSMS", names(df), value = TRUE)
        # Secutity if no columns found
  if (length(attribution_cols) == 0) {
    warning("Aucune colonne contenant 'MSMS' trouvée.")
    return(rep(FALSE, nrow(df)))
  }
        # Apply the search to pattern "CID"
  motif_found_bis <- apply(df[, attribution_cols, drop = FALSE], 1, function(row) {
    any(grepl("CID", row, ignore.case = TRUE))  # case-insensitive si utile
  })
  return(motif_found_bis)
}

#################################### function to propagate info to the pcgroup if there is only a single identification inside it
propagate_info <- function(df) {
              # Add a column 'PCGROUP_Conflict' 
  df$PCGROUP_Conflict <- NA
              # Add a column 'motif_found' by verifying if the pattern "MSMS" is present
  df$motif_found_bis <- Verify_motif_bis(df)
              # Select lines containing the pattern
  lines_with_motif <- df %>%
    filter(motif_found_bis) %>%
    select(pcgroup, ENTRY, compound, RT, mz.1, Formula, Subclass, CHEBI, Inchi, InchiKey, Smiles, annotation_confidence, lvl)
              #  Security if no propagation possible
  if (nrow(lines_with_motif) == 0) {
    message("No lines with motif - No propagation")
    return(df)
  }
              # For each line with the pattern, propagate information to lines with same pcgroup
  for (i in 1:nrow(lines_with_motif)) {
              # extract values for source line
    pcgroup_val <- lines_with_motif$pcgroup[i]
    ENTRY_val <- lines_with_motif$ENTRY[i]
    compound_val <- lines_with_motif$compound[i]
    RT_val <- lines_with_motif$RT[i]
    mz.1_val <- lines_with_motif$mz.1[i]
    Formula_val <- lines_with_motif$Formula[i]
    Subclass_val <- lines_with_motif$Subclass[i]
    CHEBI_val <- lines_with_motif$CHEBI[i]
    Inchi_val <- lines_with_motif$Inchi[i]
    InchiKey_val <- lines_with_motif$InchiKey[i]
    Smiles_val <- lines_with_motif$Smiles[i]
    annot_val <- lines_with_motif$annotation_confidence[i]
    lvl_val <- lines_with_motif$lvl[i]
              # Verify if there is lines in same pcgroup without pattern
    lines_in_group <- df %>% filter(pcgroup == pcgroup_val)
              # Verify how many lines in the pcgroup contains the pattern "[M+H]+"
    motif_count <- sum(lines_in_group$motif_found_bis)
    # If more than one motif OR zero motif, mark conflict on target lines
    if (motif_count != 1) {
      df$PCGROUP_Conflict[df$pcgroup == pcgroup_val] <- "Conflict"
      next
    }
              # If there is more than one line in the pcgroup with the pattern, no propagation
    if (motif_count <= 1) {  
              # Filtrate lines to update
      lines_to_update <- lines_in_group %>%
        filter(!motif_found_bis)
              # Update lines found without pattern
      if (nrow(lines_to_update) > 0) {
              #Propagate information to lines without pattern in the same pcgroup overwritting ols values
        df <- df %>%
          mutate(
            ENTRY = ifelse(pcgroup == pcgroup_val & !motif_found_bis, ENTRY_val, ENTRY),
            compound = ifelse(pcgroup == pcgroup_val & !motif_found_bis, compound_val, compound),
            RT = ifelse(pcgroup == pcgroup_val & !motif_found_bis, RT_val, RT),
            mz.1 = ifelse(pcgroup == pcgroup_val & !motif_found_bis, mz.1_val, mz.1),
            Formula = ifelse(pcgroup == pcgroup_val & !motif_found_bis, Formula_val, Formula),
            Subclass = ifelse(pcgroup == pcgroup_val & !motif_found_bis, Subclass_val, Subclass),
            CHEBI = ifelse(pcgroup == pcgroup_val & !motif_found_bis, CHEBI_val, CHEBI),
            Inchi = ifelse(pcgroup == pcgroup_val & !motif_found_bis, Inchi_val, Inchi),
            InchiKey = ifelse(pcgroup == pcgroup_val & !motif_found_bis, InchiKey_val, InchiKey),
            Smiles = ifelse(pcgroup == pcgroup_val & !motif_found_bis, Smiles_val, Smiles),
            annotation_confidence= ifelse(pcgroup == pcgroup_val & !motif_found_bis, annot_val, annotation_confidence),
            lvl = ifelse(pcgroup == pcgroup_val & !motif_found_bis, lvl_val, lvl)
          )
      }
    }
  }
  return(df)  
}

#################################### function to suppress data not linked with [M+H]+ in attribution column
Filter_attribution_POS <- function(df) {
              # Identify columns "attribution_partX"
  attribution_cols <- grep("^attribution_part\\d+$", names(df), value = TRUE)
              # Extract the number "partX" from the column names
  part_indices <- gsub("attribution_part", "", attribution_cols)
              # For each line
  for (i in 1:nrow(df)) {
    if (!is.na(df[i, "MSMS"]) && df[i, "MSMS"] == "CID") {
              # Find the index where the pattern "[M+H]+" is present
      keep_indices <- c()
      for (idx in part_indices) {
        col_name <- paste0("attribution_part", idx)
        val <- as.character(df[i, col_name])
        if (!is.na(val) && grepl("\\[M\\+H\\]\\+$", val) | grepl("\\[\\(M\\)\\]\\+$", val)) {
          keep_indices <- c(keep_indices, idx)
        }
      }
              # Suppress the values of columns "_partX" for the X not kept
      all_part_cols <- grep("_part\\d+$", names(df), value = TRUE)
      for (col in all_part_cols) {
              # Extract the X for each column
        x <- gsub(".*_part", "", col)
        if (!(x %in% keep_indices)) {
          df[i, col] <- NA
        }
      }
    }
  }
  return(df)
}

#################################### function to suppress data not linked with [M-H]- in attribution column
Filter_attribution_NEG <- function(df) {
              # Identify columns "attribution_partX"
  attribution_cols <- grep("^attribution_part\\d+$", names(df), value = TRUE)
              # Extract the number "partX" from the column names
  part_indices <- gsub("attribution_part", "", attribution_cols)
              # For each line
  for (i in 1:nrow(df)) {
    if (!is.na(df[i, "MSMS"]) && df[i, "MSMS"] == "CID") {
              # Find the index where the pattern "[M-H]-" is present
      keep_indices <- c()
      for (idx in part_indices) {
        col_name <- paste0("attribution_part", idx)
        val <- as.character(df[i, col_name])
        if (!is.na(val) && grepl("^\\[M\\-H\\]\\-$", val)) {
          keep_indices <- c(keep_indices, idx)
        }
      }
             # Suppress the values of columns "_partX" for the X not kept
      all_part_cols <- grep("_part\\d+$", names(df), value = TRUE)
      for (col in all_part_cols) {
            # Extract the X for each column
        x <- gsub(".*_part", "", col)
        if (!(x %in% keep_indices)) {
          df[i, col] <- NA
        }
      }
    }
  }
  return(df)
}

#################################### function to propagate to same compound (compound column)
propagation_compound <- function(df) {
  all_part_cols <- grep("_part\\d+$", names(df), value = TRUE)
  compound_cols <- grep("compound_part\\d+$", names(df), value = TRUE)
  df[compound_cols] <- lapply(df[compound_cols], as.character)
  df[all_part_cols] <- lapply(df[all_part_cols], as.character)
  df$rt <- as.numeric(df$rt)
  for (i in 1:nrow(df)) {
    if (is.na(df[i, "MSMS"])) {
      valeurs_i <- as.character(df[i, compound_cols])
      valeurs_i <- valeurs_i[!is.na(valeurs_i)]
      if (length(valeurs_i) == 0) next
      for (j in 1:nrow(df)) {
        if (i == j) next
        if (abs(df[i, "rt"] - df[j, "rt"]) >= 0.1) next
        valeurs_j <- as.character(df[j, compound_cols])
        valeurs_j <- valeurs_j[!is.na(valeurs_j)]
        shared_vals <- intersect(valeurs_i, valeurs_j)
        if (length(shared_vals) > 0 && !is.na(df[j, "MSMS"])) {
          if (is.na(df[i, "lvl"]) && !is.na(df[j, "lvl"])) {
            df[i, "lvl"] <- df[j, "lvl"]
          }
          if (is.na(df[i, "annotation_confidence"]) && !is.na(df[j, "annotation_confidence"])) {
            df[i, "annotation_confidence"] <- df[j, "annotation_confidence"]
          }
          keep_indices <- c()
          for (col in compound_cols) {
            val <- as.character(df[i, col])
            if (!is.na(val) && val %in% shared_vals) {
              x <- gsub(".*_part", "", col)
              keep_indices <- c(keep_indices, x)
            }
          }
          for (col in all_part_cols) {
            x <- gsub(".*_part", "", col)
            if (!(x %in% keep_indices)) {
              df[i, col] <- NA
            }
          }
        }
      }
    }
  }
  return(df)
}

#################################### function to propagate to same adduct number
propagate_adduct <- function(df) {
              # Identify the columns "adduct_partX"
  adduct_cols <- grep("^adduct_part\\d+$", names(df), value = TRUE)
              # Extract the lines number in the adduct columns
  extract_numeric_values <- function(row) {
    values <- row[adduct_cols]
    nums <- suppressWarnings(as.numeric(unlist(values)))
    nums[!is.na(nums)]
  }
              # Add a temporary column with the numerical values
  df$numeric_keys <- apply(df, 1, extract_numeric_values)
              # Define source ("lvl" filled) and target ("lvl" empty) lines
  df_sources <- df %>% filter(!is.na(lvl) & lvl != "")
  df_targets <- df %>% filter(is.na(lvl) | lvl == "")
  for (i in seq_len(nrow(df_targets))) {
    target_row <- df_targets[i, ]
    target_keys <- target_row$numeric_keys[[1]]
             # Find the first source line with a numerical value
    match_found <- FALSE
    for (j in seq_len(nrow(df_sources))) {
      source_keys <- df_sources$numeric_keys[[j]]
      if (length(intersect(target_keys, source_keys)) > 0) {
        df_targets[i, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")] <- df_sources[j, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")]
        match_found <- TRUE
        break
      }
    }
  }
            # Reform the data.frame
  df_final <- bind_rows(
    df_sources,
    df_targets
  ) %>%
    select(-numeric_keys) %>%
    arrange(row_number())  
  return(df_final)
}

#################################### function to propagate to same adduct number with RT<0.02
propagate_adduct_B <- function(df) {
  # Identify the columns "adduct_partX"
  adduct_cols <- grep("^adduct_part\\d+$", names(df), value = TRUE)
  # Extract the lines number in the adduct columns
  extract_numeric_values <- function(row) {
    values <- row[adduct_cols]
    nums <- suppressWarnings(as.numeric(unlist(values)))
    nums[!is.na(nums)]
  }
  # Add a temporary column with the numerical values
  df$numeric_keys <- apply(df, 1, extract_numeric_values)
  # Define source ("lvl" filled) and target ("lvl" empty) lines
  df_sources <- df %>% filter(!is.na(lvl) & lvl != "")
  df_targets <- df %>% filter(is.na(lvl) | lvl == "")
  for (i in seq_len(nrow(df_targets))) {
    target_row <- df_targets[i, ]
    target_keys <- target_row$numeric_keys[[1]]
    # Find the first source line with a numerical value
    match_found <- FALSE
    for (j in seq_len(nrow(df_sources))) {
      source_keys <- df_sources$numeric_keys[[j]]
      if (length(intersect(target_keys, source_keys)) > 0 &&
          !is.na(target_row$rt) &&
          !is.na(df_sources$rt[j]) &&
          abs(target_row$rt - df_sources$rt[j]) < 0.02
      ) {
        df_targets[i, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")] <- df_sources[j, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")]
        match_found <- TRUE
        break
      }
    }
  }
  # Reform the data.frame
  df_final <- bind_rows(
    df_sources,
    df_targets
  ) %>%
    select(-numeric_keys) %>%
    arrange(row_number())  
  return(df_final)
}

#################################### function to propagate to same adduct number if same pcgroup #################################### A tester
propagate_adduct_pcgroup <- function(df) {
  # Identify the columns "adduct_partX"
  adduct_cols <- grep("^adduct_part\\d+$", names(df), value = TRUE)
  # Extract the lines number in the adduct columns
  extract_numeric_values <- function(row) {
    values <- row[adduct_cols]
    nums <- suppressWarnings(as.numeric(unlist(values)))
    nums[!is.na(nums)]
  }
  # Add a temporary column with the numerical values
  df$numeric_keys <- apply(df, 1, extract_numeric_values)
  # Define source ("lvl" filled) and target ("lvl" empty) lines
  df_sources <- df %>% filter(!is.na(lvl) & lvl != "")
  df_targets <- df %>% filter(is.na(lvl) | lvl == "")
  for (i in seq_len(nrow(df_targets))) {
    target_row <- df_targets[i, ]
    target_keys <- target_row$numeric_keys[[1]]
    target_pcgroup <- target_row$pcgroup
    # Find the first source line with a numerical value
    match_found <- FALSE
    for (j in seq_len(nrow(df_sources))) {
      source_keys <- df_sources$numeric_keys[[j]]
      source_pcgroup <- df_sources$pcgroup[j]
      if (length(intersect(target_keys, source_keys)) > 0 && target_pcgroup == source_pcgroup) {
        df_targets[i, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")] <- df_sources[j, c("lvl", "annotation_confidence", "ENTRY", "compound", "RT", "mz.1", "Formula", "Subclass", "CHEBI", "Inchi","InchiKey" ,"Smiles")]
        match_found <- TRUE
        break
      }
    }
  }
  # Reform the data.frame
  df_final <- bind_rows(
    df_sources,
    df_targets
  ) %>%
    select(-numeric_keys) %>%
    arrange(row_number())  
  return(df_final)
}

#################################### function to propagate to same isotope for data with rt<2.5 min
copy_compound_formula_B <- function(df) {
  filled_rows <- which(!is.na(df$lvl) & !is.na(df$annotation_confidence))
  for (source_row in filled_rows) {
    source_isotope <- extract_isotope_info(df$isotopes[source_row])
    source_number <- source_isotope[1]
    source_X <- source_isotope[2]
    if (!is.na(source_number) && !is.na(source_X)) {
      target_rows <- which(df$isotopes != "" & !is.na(df$isotopes))  # Ignorer les lignes vides
      for (target_row in target_rows) {
        target_isotope <- extract_isotope_info(df$isotopes[target_row])
        target_number <- target_isotope[1]
        target_X <- target_isotope[2]
        if (!is.na(target_number) && !is.na(target_X) && target_number == source_number && target_X > source_X) {
          if (is.na(df$lvl[target_row]) || df$lvl[target_row] == "") {
            df$lvl[target_row] <- df$lvl[source_row]
            df$annotation_confidence[target_row] <- df$annotation_confidence[source_row]
            df$ENTRY[target_row] <- df$ENTRY[source_row]
            df$compound[target_row] <- df$compound[source_row]
            df$RT[target_row] <- df$RT[source_row]
            df$mz.1[target_row] <- df$mz.1[source_row]
            df$Formula[target_row] <- df$Formula[source_row]
            df$Subclass[target_row] <- df$Subclass[source_row]
            df$CHEBI[target_row] <- df$CHEBI[source_row]
            df$Inchi[target_row] <- df$Inchi[source_row]
            df$InchiKey[target_row] <- df$InchiKey[source_row]
            df$Smiles[target_row] <- df$Smiles[source_row]
          }
        }
      }
    }
  }
  return(df)
}

#################################### function to merge the adduct_partX columns
merge_part_columns_B <- function(df) {
  prefixes <- unique(sub("_part\\d+$", "", names(df)))
  for (prefix in prefixes) {
    cols_to_merge <- grep(paste0("^", prefix, "_part"), names(df), value = TRUE)
    if (length(cols_to_merge) > 0) {
      df[[prefix]] <- apply(df[, cols_to_merge, drop = FALSE], 1, function(row) {
        paste(na.omit(row), collapse = " ")
      })
      df <- df[, !(names(df) %in% cols_to_merge)]
    }
  }
  return(df)
}

#################################### function to verify the presence of [M+Na]+ in attribution column
Verify_motif_Na <- function(df) {
        # Select column containing "attribution"
  attribution_cols <- grep("attribution", names(df), value = TRUE)
        # Apply line per line
  sapply(seq_len(nrow(df)), function(i) {
    row <- df[i, attribution_cols]
    has_motif <- any(grepl("^\\[M\\+Na\\]\\+$", row))  # Focus on patter "[M+Na]+"
    has_lvl <- !is.na(df$lvl[i]) && df$lvl[i] != ""
    has_motif && has_lvl
  })
}

#################################### function to verify the presence of [M+Cl]- in attribution column
Verify_motif_Cl <- function(df) {
  attribution_cols <- grep("attribution", names(df), value = TRUE)
  sapply(seq_len(nrow(df)), function(i) {
    row <- df[i, attribution_cols]
    has_motif <- any(grepl("^\\[M\\+Cl\\]\\-$", row))  # Focus on patter "[M+Cl]-"
    has_lvl <- !is.na(df$lvl[i]) && df$lvl[i] != ""
    has_motif && has_lvl
  })
}

#################################### function to fill the column ppm_score
Ppm_score <- function(data, mzmin, mzmax) {
  # Vérifier que les colonnes existent
  if (!(mzmin %in% names(data)) || !(mzmax %in% names(data))) {
    stop("Les colonnes spécifiées n'existent pas dans le data.frame.")
  }
  # Calcul du delta
  data <- data %>%
    dplyr::mutate(ppm_diff = .data[[mzmax]] - .data[[mzmin]])
  # Calcul de la valeur max
  z_max <- max(data$ppm_diff, na.rm = TRUE)
  # Ajout du score
  if (z_max == 0) {
    data$ppm_score <- 1
  } else {
    data <- data %>%
      dplyr::mutate(ppm_score = 1 - (ppm_diff / z_max))
  }
  # Nettoyage
  data <- data %>% dplyr::select(-ppm_diff)
  return(data)
}


#################################### function to fill the column rt_score
rt_score <- function(data, rtmin, rtmax) {
  data <- data %>%
    mutate(rt_diff = .data[[rtmax]] - .data[[rtmin]])
  z_max <- max(data$rt_diff, na.rm = TRUE)
  if (z_max == 0) {
    data$rt_score <- 1 
  } else {
    data <- data %>%
      mutate(rt_score = 1 - (rt_diff / z_max))
  }
  data %>% select(-rt_diff)
}

