
library(DBI)
library(RSQLite)


#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### List tables in db #####

dbGetInfo(con)
dbListTables(con)
dbReadTable(con, "biomolecules")


#### Fetch biomolecules table and biomolecules metadata tables, to be combined ####

biomolecules <- dbReadTable(con, "biomolecules")
#lipidomics_metadata <- dbReadTable(con, "lipidomics_metadata")
proteomics_metadata <- dbReadTable(con, "proteomics_metadata")
transcriptomics_metadata <- dbReadTable(con, "rnaseq_metadata")
transcriptomics_metadata2 <- dbReadTable(con, "rnaseq_measurements")

dbDisconnect(con)

