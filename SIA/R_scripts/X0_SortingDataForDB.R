
########################### load libraries ################################

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

########################### Define path ################################

base_dir <- "D:/LongCOVID"
data_dir <- file.path(base_dir, "Important Tables")
metadata_dir <- file.path(base_dir, "Metadata")
output_dir <- file.path(data_dir, "db_formatted_tables")  # Store output in a 'Results' folder


if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


data_file <- file.path(data_dir, "20250207_longCOVID_AddedFilterationColumns_db.csv")
metadata_file <- file.path(metadata_dir, "rawfiles_Metadata.csv")
rawfiles <- file.path(data_dir, "lipidomics_rawfiles_db.csv")
lipidmaps_file <- file.path(metadata_dir, "LipidMAPS_Classification.csv") 

######################data_dir########################### Read Data ################################

data <- read_csv(data_file)
metadata <- read_csv(metadata_file)
rawfiles <- read_csv(rawfiles)
lipidmaps <- read_csv(lipidmaps_file) 


########################### Biomolecules Table ################################

biomolecules <- data.frame(
  UniqueID = data$UniqueID,
  biomolecule_id = seq(38166, 38166 + nrow(data) - 1), # Assign sequential IDs
  standardized_name = ifelse(data$Annotation == "", 
                            paste0("Unknown_", data$UniqueID), 
                            paste0(data$Annotation, "_RTmz_", data$UniqueID)), 
  omics_id = 2, 
  keep = ifelse(data$KeepQuantitation == TRUE & (is.na(data$FoundIn50perFiles) | data$FoundIn50perFiles == TRUE),
                "1", 
                ifelse(data$KeepQuantitation == FALSE, 
                       paste0("0; ", data$FilterReasons), 
                       ifelse(data$FoundIn50perFiles == FALSE, 
                              "0; Notfound_in_50pecent_samples", 
                              "1")))
)


print(biomolecules)


########################### Biomolecules Metadata ##################################

data <- data %>%
  left_join(lipidmaps, by = "LipidClass")  


biomolecules_metadata <- data %>%
  mutate(
    biomolecule_id = seq(38166, 38166 + nrow(data) - 1),
    standardized_name = ifelse(Annotation == "", paste0("Unknown_", UniqueID), 
                               paste0(Annotation, "_RTmz_", UniqueID)),
    across(c(LipidCategoty, LipidClass, MainClass, SubClass, 
             NumFattyAcylCarbons, NumFattyAcylUnsaturations), 
           ~ifelse(grepl("^Unknown_", standardized_name), NA, .)),
    AverageUnsaturations = ifelse(grepl("^Unknown_", standardized_name), 
                                  NA, NumFattyAcylUnsaturations / NumFattyAcyls),
    UnsaturationLevel = case_when(
      AverageUnsaturations <= 0.5 ~ "Very Low",
      AverageUnsaturations <= 1 ~ "Low",
      AverageUnsaturations <= 2 ~ "Medium",
      AverageUnsaturations <= 3 ~ "High",
      AverageUnsaturations > 3 ~ "Very High",
      TRUE ~ "Unknown"
    ),
    omics_id = 2,
    keep = ifelse(KeepQuantitation & (is.na(FoundIn50perFiles) | FoundIn50perFiles), "1", 
                  ifelse(!KeepQuantitation, paste0("0; ", FilterReasons), 
                         ifelse(!FoundIn50perFiles, "0; Notfound_in_50pecent_samples", "1")))
  ) %>%
  select(biomolecule_id, UniqueID, standardized_name, LipidCategoty, LipidClass, 
         MainClass, SubClass, NumFattyAcylCarbons, NumFattyAcylUnsaturations, 
         AverageUnsaturations, UnsaturationLevel, omics_id, keep)



########################### Lipidomics meansurements Table ###############################

data_subset = data %>% select(UniqueID, matches(".raw$"))

lipidomics_measurements <- data_subset %>%
  pivot_longer(cols = -UniqueID, names_to = "rawfile_name_R", values_to = "raw_abundance") %>%
  mutate(
    measurement_id = row_number(),  # Assign unique measurement IDs
    normalized_abundance = log2(raw_abundance)  # Log2 normalization
  )


lipidomics_measurements <- lipidomics_measurements %>%
  left_join(biomolecules %>% select(UniqueID, standardized_name, biomolecule_id), by = "UniqueID") %>%
  left_join(rawfiles %>% select(rawfile_id, rawfile_name_R), by = "rawfile_name_R")



########################### db formatted Tables ###############################

lipidomics_measurements <- lipidomics_measurements %>%
  select(measurement_id, biomolecule_id, UniqueID, standardized_name,  rawfile_id, rawfile_name_R, raw_abundance, normalized_abundance)

biomolecules <- biomolecules %>%
  select(biomolecule_id, UniqueID, standardized_name, omics_id, keep)

########################### Save Results ################################

biomolecules_file <- file.path(output_dir, "biomolecules_table.csv")
write_csv(biomolecules, biomolecules_file)


lipid_measurements_file <- file.path(output_dir, "lipidomics_measurements_table.csv")
write_csv(lipidomics_measurements, lipid_measurements_file)


biomolecules_metadata_file <- file.path(output_dir, "biomolecules_metadata.csv")
write_csv(biomolecules_metadata, biomolecules_metadata_file)
