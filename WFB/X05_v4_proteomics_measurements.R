#### make a proteomics_measurements table for the long covid database

library(data.table)
library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(rrcovNA)



#### Section Detailing How I got to Protein Quantity from Spectronaut Output ####
# read in report
file <- "20250203_170537_20250123_WFB_JaitovichLongCovid_AllPlates_Samples_dDIA_stringent_htrms_NPANPB_Report.tsv"

names(fread(paste0("data/spectronaut_output/", file), nrows = 0))

spec <- fread(paste0("data/spectronaut_output/", file), 
              sep = "\t",
              select = c("R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Quantity"))

#### PGs ####
# Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
spec %>%
  dplyr::summarise(n = dplyr::n(), .by = c(R.Condition, PG.ProteinGroups, R.FileName)) %>%
  dplyr::filter(n > 1L)   

spec_PG <- pivot_wider(spec,
                       names_from = R.FileName,
                       values_from = PG.Quantity,
                       values_fill = NA,
                       values_fn = list(PG.Quantity = mean))

colSums(!is.na(spec_PG))

# fixing columns so they all line up
spec_PG <- spec_PG %>%
  group_by(PG.ProteinGroups) %>%
  summarise_all(~ na.omit(.) %>% .[1]) %>%
  ungroup()

plot(colSums(!is.na(spec_PG)))

# save ungrouped pg quant matrix
fwrite(spec_PG, file = "data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined.csv")






#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

rawfiles_all <- dbGetQuery(con, "SELECT rawfile_name, rawfile_id 
                           FROM rawfiles_all")

dbDisconnect(con)

#### Make NPA and NPB long format protein quant for database ####
# read in matrices
data_matrix <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined.csv")

# remove condition column
data_matrix <- data_matrix[,-2]

# filter for only ids in ids_to_keep (see 02_Pep_PG_Numbers for filtering strategy, 50% completeness in one study group at least)
ids_to_keep <- readLines("data/processed/proteomics_ids_to_keep.txt")
data_matrix1 <- data_matrix %>%
  filter(PG.ProteinGroups %in% ids_to_keep)

# do normalization then export to perseus for imputation and reimport
# Log2 transform the protein groups
data_matrix_tf <- log2(data_matrix1[, -1])
data_matrix_tf <- data.frame(ID = data_matrix1$PG.ProteinGroups, data_matrix_tf, stringsAsFactors = FALSE)


# impute with ImpSeq
data_matrix_tf1 <- data_matrix_tf %>%
  dplyr::select(where(is.double))

start_time <- Sys.time()
  data_imp <- impSeq(data_matrix_tf1)
end_time <- Sys.time()
print(end_time - start_time)

data_matrix_tf_imp <- data_matrix_tf %>%
  dplyr::select(ID) %>%
  cbind(data_imp)



#fix column names to match rawfile_name in rawfiles_all table
colnames(data_matrix_tf_imp) <- gsub("^X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^X", "", colnames(data_matrix))
colnames(data_matrix_tf_imp) <- gsub("^N\\.\\.X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^N\\.\\.X", "", colnames(data_matrix))

#### Make longer and merge abundances with normalized and imputed ####
data_matrix_long <- data_matrix %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "raw_abundance")

data_matrix_tf_imp_long <- data_matrix_tf_imp %>%
  rename(PG.ProteinGroups = ID) %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "normalized_abundance")

data_matrix_all_long <- data_matrix_long %>%
  left_join(data_matrix_tf_imp_long, by = c("PG.ProteinGroups", "rawfile_name"))


#add rawfile_id
#fix hyphens turn into dots
rawfiles_all$rawfile_name <- sub("-", ".", rawfiles_all$rawfile_name, fixed = T)
data_matrix_all_long <- data_matrix_all_long %>%
  inner_join(rawfiles_all, by = "rawfile_name")

#add biomolecule_id
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name 
                           FROM biomolecules")
dbDisconnect(con)

colnames(data_matrix_all_long)[colnames(data_matrix_all_long) == "PG.ProteinGroups"] <- "standardized_name"

data_matrix_all_long <- data_matrix_all_long %>%
  left_join(biomolecules, by = "standardized_name")

data_matrix_all_long$measurement_id <- seq(1,nrow(data_matrix_all_long))

# how many measurements in the whole experiment
table(!is.na(data_matrix_all_long$raw_abundance))


## Add the new biomolecule_id, one for each ProteinGroup_NP combination ---- 
proteomics <- data_matrix_all_long

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
rawfiles <- dbGetQuery(con, 'SELECT rawfile_name, Sample, sample_id, ome_id, keep, rawfile_id, run_type
                           FROM rawfiles_all')
dbDisconnect(con)

rawfiles <- rawfiles %>%
  dplyr::select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl('Sample', run_type))

df <- proteomics %>%
  left_join(rawfiles, by = 'rawfile_id') %>%
  mutate(NP = word(rawfile_name, 6, 6, "_"))

asdf <- df %>%
  dplyr::select(-biomolecule_id) %>%
  group_by(standardized_name, NP) %>%
  summarize() %>%
  ungroup() %>%
  mutate(biomolecule_id = row_number())

df <- df %>%
  dplyr::select(-biomolecule_id) %>%
  left_join(asdf, by = c("standardized_name", "NP"))

df <- df %>%
  arrange(biomolecule_id)

proteomics_measurement <- df %>%
  mutate(measurement_id = row_number()) %>%
  dplyr::select(measurement_id, standardized_name, biomolecule_id, rawfile_name, rawfile_id, raw_abundance, normalized_abundance)


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "proteomics_measurement", proteomics_measurement, overwrite = T)
dbDisconnect(con)













# drop extra table
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbExecute(con, "DROP TABLE IF EXISTS proteomics_measurements")
dbDisconnect(con)





