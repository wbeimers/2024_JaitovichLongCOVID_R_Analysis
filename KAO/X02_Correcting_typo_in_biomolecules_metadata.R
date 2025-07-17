#### Editing the biomolecules_metadata table #### 

library(DBI)
library(RSQLite)

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### List tables in db #####

dbGetInfo(con)
dbListTables(con)

#### Fetch biomolecules table and biomolecules_metadata tables ####

biomolecules <- dbReadTable(con, "biomolecules")
biomolecules_metadata <- dbReadTable(con, "biomolecules_metadata")

dbDisconnect(con)

lipidomics_metadata <- read.csv("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/biomolecules_metadata_v2.csv")

#### Substitute  metadata_type = "Lipd_Class" with "Lipid_Class" #### 

biomolecules_metadata$metadata_type[biomolecules_metadata$metadata_type == "Lipd_Class"] <- "Lipid_Class"
unique(biomolecules_metadata$metadata_type)

#### Replace biomolecules_metadata table in the db #### 
#dim(biomolecules_metadata)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbWriteTable(con, "biomolecules_metadata", biomolecules_metadata, overwrite = T)

## check 

dbListTables(con)

## disconnect
dbDisconnect(con)

## Change made 3/21/2025 at 12:30pm 
