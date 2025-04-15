#### add Salma's tables from P:\Projects\WFB_SIA_2024_Jaitovich_LongCOVID\Lipidomics\db_formatted_tables

library(data.table)
library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(rrcovNA)



#### Add lipids to the biomolecules table ####

biomolecules_lipids <- fread("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/biomolecules_table.csv") %>%
  select(-UniqueID) %>%
  mutate(standardized_name = standarized_name) %>%
  select(-standarized_name)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbWriteTable(con, "biomolecules", biomolecules_lipids, append = T)

dbDisconnect(con)



#### Add lipidomics measurements table ####

measurements_lipids <- fread("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/lipidomics_measurements_table.csv")

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbExecute(con, "DROP TABLE IF EXISTS lipidomics_measurments")

dbWriteTable(con, "lipidomics_measurements", measurements_lipids, overwrite = T)

dbDisconnect(con)



#### Add new rawfiles table ####

rawfiles_lipids <- fread("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/lipidomics_rawfiles_db.csv") %>%
  rename(Sample = sample) %>%
  rename(dilution_concentration = DilutionConcentration)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

rawfiles_all <- dbGetQuery(con, "SELECT * 
                       FROM rawfiles_all")

dbDisconnect(con)

rawfiles_all <- rawfiles_all %>%
  select(-date, -time) %>%
  filter(ome_id == 1)

rawfiles_all <- rbind(rawfiles_all, rawfiles_lipids)  


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbExecute(con, "DROP TABLE IF EXISTS rawfiles_all")

dbWriteTable(con, "rawfiles_all", rawfiles_all, overwrite = T)

dbDisconnect(con)



## Put new biomolecule IDs in the lipid table

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
lipidomics <- dbGetQuery(con, 'SELECT rawfile_name, raw_abundance, normalized_abundance, measurement_id, standardized_name
                                    FROM lipidomics_measurements')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                                 FROM biomolecules')
dbDisconnect(con)

lipidomics1 <- lipidomics %>%
  left_join(biomolecules %>% dplyr::select(biomolecule_id, standardized_name), by = "standardized_name")

lipidomics1 <- lipidomics1 %>%
  dplyr::select(measurement_id, standardized_name, biomolecule_id, everything())


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "lipidomics_measurements", lipidomics1, append = T, row.names = F)
dbDisconnect(con)

  

## Lipidomics Metadata ----

metadata_lipids <- fread("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/biomolecules_metadata.csv") %>%
  dplyr::select(-omics_id, -keep)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "lipidomics_metadata", metadata_lipids, row.names = F)
dbDisconnect(con)



#### Remove a rawfile from the samples ####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

rawfiles_all <- dbGetQuery(con, "SELECT * 
                       FROM rawfiles_all")

dbDisconnect(con)

rawfiles_all[rawfiles_all$rawfile_id == 979, "keep"] <- "0; replicate_injection"
  
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbExecute(con, "DROP TABLE IF EXISTS rawfiles_all")

dbWriteTable(con, "rawfiles_all", rawfiles_all, overwrite = T)

dbDisconnect(con)