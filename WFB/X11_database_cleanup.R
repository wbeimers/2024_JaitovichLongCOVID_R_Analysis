
## Clean Up Database ----
# 1. remove standardized_name and rawfile_name from biomolecule tables
# 2. remove lipidomics_runs and proteomics_runs and
# 2. remove rawfile_name_R from lipidomics_measurements


library(DBI)
library(RSQLite)


# Connect to your .sqlite file
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


## Remove columns from measurements tables ----
# lipidomics
dbExecute(con, "ALTER TABLE lipidomics_measurements DROP COLUMN standardized_name;")
dbExecute(con, "ALTER TABLE lipidomics_measurements DROP COLUMN rawfile_name_R;")

# proteomics
dbExecute(con, "ALTER TABLE proteomics_measurement DROP COLUMN standardized_name;")
dbExecute(con, "ALTER TABLE proteomics_measurement DROP COLUMN rawfile_name;")

# transcriptomics
dbExecute(con, "ALTER TABLE rnaseq_measurements DROP COLUMN standardized_name;")


## Drop lipidomics and proteomics runs tables ----
dbExecute(con, "DROP TABLE IF EXISTS lipidomics_runs;")
dbExecute(con, "DROP TABLE IF EXISTS proteomics_runs;")


dbDisconnect(con)



