#### make a proteomics_measurements table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Section Detailing How I got to Protein Quantity from Spectronaut Output ####
#read in report
file <- "20240815_154253_20240815_WFB_JaitovichLongCOVID_AllPlates_SamplesQCs_DIA_NPA_Report.tsv"

spec <- read.delim(paste0("data/spectronaut_output/", file), sep = "\t")

####PGs####
#Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
selects <- (is.na(spec$PG.Qvalue) | spec$PG.Qvalue <= 0.01) &
  (is.na(spec$EG.Qvalue) | spec$EG.Qvalue <= 0.01) 

spec_raw <- spec[selects, c("R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Quantity")]

spec_raw %>%
  dplyr::group_by(PG.ProteinGroups, R.Condition) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)  

spec_PG <- pivot_wider(spec_raw,
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
write.csv(spec_PG, file = "data/processed/PG_Matrix_AllPlates_QCsSamples_NPA.csv")

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
write.csv(spec_PG_NPs, file = "data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")




#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### Read in the rawfiles_all for mapping ####

rawfiles_all <- dbGetQuery(con, "SELECT rawfile_name, rawfile_id 
                           FROM rawfiles_all")



#### Make NPA and NPB long format protein quant for database ####
#read in matrices
NPA_matrix <- read.csv("data/processed/PG_Matrix_AllPlates_QCsSamples_NPA.csv")
NPB_matrix <- read.csv("data/processed/PG_Matrix_AllPlates_QCsSamples_NPB.csv")

#remove first column
NPA_matrix <- NPA_matrix[-1]
NPB_matrix <- NPB_matrix[-1]

#do normalization then export to perseus for imputation and reimport
#Log2 transform the protein groups
NPA_matrix_tf <- log2(NPA_matrix[, -c(1, 2)])
NPA_matrix_tf <- data.frame(ID = NPA_matrix$PG.ProteinGroups, NPA_matrix_tf, stringsAsFactors = FALSE)

NPB_matrix_tf <- log2(NPB_matrix[, -c(1, 2)])
NPB_matrix_tf <- data.frame(ID = NPB_matrix$PG.ProteinGroups, NPB_matrix_tf, stringsAsFactors = FALSE)


#write csv files to go impute in a different program
write.csv(NPA_matrix_tf, file = "data/processed/PG_Matrix_AllPlates_Samples_NPA_tf.csv")
write.csv(NPB_matrix_tf, file = "data/processed/PG_Matrix_AllPlates_Samples_NPB_tf.csv")

#read in imputed file from perseus
NPA_matrix_tf_imp <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPA_tf_imp.csv")
NPB_matrix_tf_imp <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPB_tf_imp.csv")


#fix column names to match rawfile_name in rawfiles_all table
colnames(NPA_matrix) <- gsub("^X", "", colnames(NPA_matrix))
colnames(NPB_matrix) <- gsub("^X", "", colnames(NPB_matrix))

colnames(NPA_matrix_tf_imp) <- gsub("^N\\.\\.X", "", colnames(NPA_matrix_tf_imp))
colnames(NPB_matrix_tf_imp) <- gsub("^N\\.\\.X", "", colnames(NPB_matrix_tf_imp))

NPA_matrix <- NPA_matrix[,-2]
NPB_matrix <- NPB_matrix[,-2]


#### Make longer and merge abundances with normalized and imputed ####
NPA_matrix_long <- NPA_matrix %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "raw_abundance")

NPB_matrix_long <- NPB_matrix %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "raw_abundance")


NPA_matrix_tf_imp_long <- NPA_matrix_tf_imp %>%
  pivot_longer(cols = -T..ID, names_to = "rawfile_name", values_to = "normalized_abundance")

NPB_matrix_tf_imp_long <- NPB_matrix_tf_imp %>%
  pivot_longer(cols = -T..ID, names_to = "rawfile_name", values_to = "normalized_abundance")

#change colnames to match
colnames(NPA_matrix_tf_imp_long)[colnames(NPA_matrix_tf_imp_long) == "T..ID"] <- "PG.ProteinGroups"
colnames(NPB_matrix_tf_imp_long)[colnames(NPB_matrix_tf_imp_long) == "T..ID"] <- "PG.ProteinGroups"


NPA_matrix_all_long <- merge(NPA_matrix_long, NPA_matrix_tf_imp_long, by = c("PG.ProteinGroups", "rawfile_name"))
NPB_matrix_all_long <- merge(NPB_matrix_long, NPB_matrix_tf_imp_long, by = c("PG.ProteinGroups", "rawfile_name"))

#add rawfile_id
#fix hyphens turn into dots
rawfiles_all$rawfile_name <- sub("-", ".", rawfiles_all$rawfile_name, fixed = T)
NPA_matrix_all_long <- merge(NPA_matrix_all_long, rawfiles_all, by = "rawfile_name")
NPB_matrix_all_long <- merge(NPB_matrix_all_long, rawfiles_all, by = "rawfile_name")

#add biomolecule_id
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name 
                           FROM biomolecules")

colnames(NPA_matrix_all_long)[colnames(NPA_matrix_all_long) == "PG.ProteinGroups"] <- "standardized_name"
colnames(NPB_matrix_all_long)[colnames(NPB_matrix_all_long) == "PG.ProteinGroups"] <- "standardized_name"

NPA_matrix_all_long <- merge(NPA_matrix_all_long, biomolecules, by = "standardized_name")
NPB_matrix_all_long <- merge(NPB_matrix_all_long, biomolecules, by = "standardized_name")

#add NPB to the end of NPA
PGs_matrix_all_long <- rbind(NPA_matrix_all_long, NPB_matrix_all_long)
PGs_matrix_all_long$measurement_id <- seq(1,nrow(PGs_matrix_all_long))

# how many measurements in the whole experiment
table(!is.na(PGs_matrix_all_long$raw_abundance))


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "proteomics_measurement", PGs_matrix_all_long, overwrite = T)
dbDisconnect(con)

