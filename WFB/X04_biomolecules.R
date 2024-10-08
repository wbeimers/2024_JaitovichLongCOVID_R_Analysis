#### make a biomolecules table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


####biomolecules####

#import a protein group and peptide table containing all IDed peptides and proteins, and give them identifiers



