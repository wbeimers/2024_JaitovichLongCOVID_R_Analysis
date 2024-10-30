#### make a biomolecules table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


####biomolecules####

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

