#### make a proteomics_measurements table for the long covid database

library(data.table)
library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(rrcovNA)



#### Section Detailing How I got to Protein Quantity from Spectronaut Output ####
#read in report
file <- "20250203_170537_20250123_WFB_JaitovichLongCovid_AllPlates_Samples_dDIA_stringent_htrms_NPANPB_Report.tsv"

names(fread(paste0("data/spectronaut_output/", file), nrows = 0))

spec <- fread(paste0("data/spectronaut_output/", file), 
              sep = "\t",
              select = c("R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Quantity"))

####PGs####
#Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
spec %>%
  dplyr::summarise(n = dplyr::n(), .by = c(R.Condition, PG.ProteinGroups, R.FileName)) %>%
  dplyr::filter(n > 1L)   

spec_PG <- pivot_wider(spec,
                       names_from = R.FileName,
                       values_from = PG.Quantity,
                       values_fill = NA,
                       values_fn = list(PG.Quantity = mean))

colSums(!is.na(spec_PG))

#fixing columns so they all line up
spec_PG <- spec_PG %>%
  group_by(PG.ProteinGroups) %>%
  summarise_all(~ na.omit(.) %>% .[1]) %>%
  ungroup()

plot(colSums(!is.na(spec_PG)))

#save ungrouped pg quant matrix
fwrite(spec_PG, file = "data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined.csv")

#Make sure I have an NPA and NPB file since I searched separately
spec_PG_NPA <- spec_PG

#merge NPA and NPB dataframes
spec_PG <- merge(spec_PG_NPA, spec_PG_NPB, by = "PG.ProteinGroups", all = T)


#make dataframes just of samples, QCs

#QCs
QCs_names <- spec_PG[,grepl(paste(c("PG.ProteinGroups","QC"), collapse = "|"), names(spec_PG))]
spec_PG_QCs <- data.frame(QCs_names)
write.csv(spec_PG_QCs, file = "data/processed/PG_Matrix_AllPlates_QCs.csv")

#samples
samples_names <- spec_PG[,grepl(paste(c("PG.ProteinGroups", "NPA", "NPB"), collapse = "|"), names(spec_PG))]
spec_PG_samples <- data.frame(samples_names)
write.csv(spec_PG_samples, file = "data/processed/PG_Matrix_AllPlates_Samples.csv")


#use this to take either NPA or NPB for all samples for each pg, depending on which set is more complete
#make a filter to choose whichever row has fewer missing values for NPA or NPB. If same, choose NPA.
NPA <- apply(spec_PG_samples[,grep("NPA", names(spec_PG_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPA[is.na(NPA)] <-0

NPB <- apply(spec_PG_samples[,grep("NPB", names(spec_PG_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPB[is.na(NPB)] <-0

filter_NP <- NPA>=NPB
table(is.na(filter_NP))

#separate NPA and NPB into different dataframes
NPA_df <- spec_PG_samples[,c(1,grep("NPA", names(spec_PG_samples)))]
NPB_df <- spec_PG_samples[,c(1,grep("NPB", names(spec_PG_samples)))]

#apply filter to choose correct rows for either dataframe
NPA_df_filter <- NPA_df[filter_NP,]
NPB_df_filter <- NPB_df[!filter_NP,]
names(NPA_df_filter) <- sub(pattern = "_NPA", replacement = "", x = names(NPA_df_filter))
names(NPB_df_filter) <- sub(pattern = "_NPB", replacement = "", x = names(NPB_df_filter))
NPA_df_filter_1 <- NPA_df_filter[,c(1, order(names(NPA_df_filter)[-1])+1)]
NPB_df_filter_1 <- NPB_df_filter[,c(1, order(names(NPB_df_filter)[-1])+1)]

spec_PG_NPs <- rbind(NPA_df_filter, NPB_df_filter)
#save combined quant matrix
write.csv(spec_PG_NPs, file = "data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined_NPs.csv")




#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### Read in the rawfiles_all for mapping ####

rawfiles_all <- dbGetQuery(con, "SELECT rawfile_name, rawfile_id 
                           FROM rawfiles_all")



#### Make NPA and NPB long format protein quant for database ####
# read in matrices
data_matrix <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined.csv")

# remove condition column
data_matrix <- data_matrix[,-2]

# do normalization then export to perseus for imputation and reimport
# Log2 transform the protein groups
data_matrix_tf <- log2(data_matrix[, -1])
data_matrix_tf <- data.frame(ID = data_matrix$PG.ProteinGroups, data_matrix_tf, stringsAsFactors = FALSE)

# Impute with Marcel's function

imputeIntensities_mm <- function(df, seed_value = 1234, left_shift = 1.8, sd_width = 0.3) {
  # Set seed for reproducibility
  set.seed(seed_value)
  # Replace all exact matches of 0 with NA in numeric columns
  df <- df %>%
    mutate(across(where(is.numeric), ~ na_if(.x, 0)))
  # Initialize lists to store new columns and summary information
  imputationColumns_list <- list()
  imputationSummary_list <- list()
  # Loop over each numeric column to perform imputation
  numeric_columns <- names(df)[sapply(df, is.numeric)]
  for (col in numeric_columns) {
    current_data <- df[[col]]
    # Print progress message by default
    message(paste("Processing column:", col))
    # Calculate standard deviation (SD) of non-NA values
    sd_val <- sd(current_data, na.rm = TRUE)
    # Generate imputed values based on a left-centric normal distribution
    imputedValues_vec <- rnorm(
      sum(is.na(current_data)),  # Number of missing values to impute
      mean = mean(current_data, na.rm = TRUE) - left_shift * sd_val,  # Shift the mean to the left
      sd = sd_width * sd_val  # Use a narrower SD (specified by sd_width)
    )
    # Create a binary indicator column for imputed values
    isImputedBoolean_vec <- ifelse(is.na(current_data), TRUE, FALSE)
    # Combine imputed and original values
    combinedValues_vec <- replace(current_data, is.na(current_data), imputedValues_vec)
    # Store the newly created columns in the list
    imputationColumns_list[[paste0(col, "_isImputedBoolean")]] <- isImputedBoolean_vec
    imputationColumns_list[[paste0(col, "_withImputedValues")]] <- combinedValues_vec
    # Store the summary statistics in the list
    validCount_val <- sum(!is.na(current_data))
    missingCount_val <- sum(is.na(current_data))
    totalCount_val <- length(current_data)
    vvPercentage_val <- (validCount_val / totalCount_val) * 100
    mvPercentage_val <- (missingCount_val / totalCount_val) * 100
    imputationSummary_list[[col]] <- data.frame(
      Run = col,
      validValues = validCount_val,
      missingValues = missingCount_val,
      totalValues = totalCount_val,
      vvPercentage = vvPercentage_val,
      mvPercentage = mvPercentage_val
    )
  }
  # Combine the list of summary data frames into a single data frame
  imputationSummary_df <- dplyr::bind_rows(imputationSummary_list)
  # Add the new columns to the original data frame
  imputed_df <- bind_cols(df, imputationColumns_list %>% as.data.frame()) %>%
    select(order(names(.)))
  # Return the list with imputed data and summary data frames
  list(imputed_data = imputed_df, summary_data = imputationSummary_df)
}

functionOutput_list <- imputeIntensities_mm(data_matrix_tf)
data_imp <- functionOutput_list$imputed_data
data_summ <- functionOutput_list$summary_data

# organize wide data frame with imputed values
data_imp <- data_imp %>%
  select(ID, contains("_withImputedValues"))
colnames(data_imp) <- word(colnames(data_imp), 1, 6, "_")
colnames(data_imp)[1] <- "PG.ProteinGroups"
data_matrix_tf_imp <- data_imp


#fix column names to match rawfile_name in rawfiles_all table
colnames(data_matrix_tf_imp) <- gsub("^X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^X", "", colnames(data_matrix))
colnames(data_matrix_tf_imp) <- gsub("^N\\.\\.X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^N\\.\\.X", "", colnames(data_matrix))

#### Make longer and merge abundances with normalized and imputed ####
data_matrix_long <- data_matrix %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "raw_abundance")

data_matrix_tf_imp_long <- data_matrix_tf_imp %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "normalized_abundance")

data_matrix_all_long <- merge(data_matrix_long, data_matrix_tf_imp_long, by = c("PG.ProteinGroups", "rawfile_name"))


#add rawfile_id
#fix hyphens turn into dots
rawfiles_all$rawfile_name <- sub("-", ".", rawfiles_all$rawfile_name, fixed = T)
data_matrix_all_long <- data_matrix_all_long %>%
  inner_join(rawfiles_all, by = "rawfile_name")

#add biomolecule_id
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name 
                           FROM biomolecules")

colnames(data_matrix_all_long)[colnames(data_matrix_all_long) == "PG.ProteinGroups"] <- "standardized_name"

data_matrix_all_long <- data_matrix_all_long %>%
  left_join(biomolecules, by = "standardized_name")

data_matrix_all_long$measurement_id <- seq(1,nrow(data_matrix_all_long))

# how many measurements in the whole experiment
table(!is.na(data_matrix_all_long$raw_abundance))


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "proteomics_measurement", data_matrix_all_long, overwrite = T)
dbDisconnect(con)

