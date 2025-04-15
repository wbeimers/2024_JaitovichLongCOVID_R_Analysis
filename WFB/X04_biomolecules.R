#### make a biomolecules table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


#### biomolecule s####

#import a protein group  table containing all IDed peptides and proteins, and give them identifiers
all_PGs <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")
all_PGs <- all_PGs[-1]
all_PGs <- tibble::rownames_to_column(all_PGs, "biomolecule_id")
all_PGs <- all_PGs[, c(1, 2)]
colnames(all_PGs)[colnames(all_PGs) == "PG.ProteinGroups"] <- "standardized_name"
all_PGs$omics_id <- 1
all_PGs$keep <- 1


dbWriteTable(con, "biomolecules", all_PGs, append = T)

dbDisconnect(con)


#### Changing the Protein Groups ----
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

db_table <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules")

# remove proteomics data
filtered_data <- db_table %>%
  filter(omics_id != 1) 

all_PGs <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPAandNPBcombined.csv")
all_PGs <- tibble::rownames_to_column(all_PGs, "biomolecule_id")
all_PGs <- all_PGs[, c(1, 2)]
colnames(all_PGs)[colnames(all_PGs) == "PG.ProteinGroups"] <- "standardized_name"
all_PGs$omics_id <- 1
all_PGs$keep <- 1

biomolec <- all_PGs %>%
  rbind(filtered_data)

dbWriteTable(con, "biomolecules", biomolec, overwrite = T)

#### Making a unique biomolecule_id for each ProteinGroup_NP combination. So 2 biomolecule_id per protein group depending on which NP it came from. For filtering ----

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
db_table <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules")
dbDisconnect(con)

# only take proteomics runs and make duplicates of each row for NPA/NPB comparison

db_table_proteins <- db_table %>%
  filter(omics_id == 1) %>%
  dplyr::select(-biomolecule_id)
db_table_proteins <- db_table_proteins[rep(1:nrow(db_table_proteins), each = 2), ]
db_table_other <- db_table %>%
  filter(omics_id != 1) %>%
  dplyr::select(-biomolecule_id)
db_table_all <- db_table_proteins %>%
  bind_rows(db_table_other) %>%
  mutate(biomolecule_id = row_number()) %>% 
  dplyr::select(biomolecule_id, everything())

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "biomolecules", db_table_all, overwrite = T)
dbDisconnect(con)


## Fix keep column ----

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
db_table <- dbGetQuery(con, "SELECT *
                           FROM biomolecules")
dbDisconnect(con)


db_table1 <- db_table %>%
  mutate(keep = case_when(
    keep == "1.0" ~ "1",
    T ~ keep
  ))


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "biomolecules", db_table1, overwrite = T)
dbDisconnect(con)
