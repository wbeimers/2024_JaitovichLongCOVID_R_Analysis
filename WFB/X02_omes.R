#### make an omes table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


####omes table####

omes <- data.frame(omics_id = c(1, 2),
                   omics_name = c("Proteomics", "Lipidomics"))


dbWriteTable(con, "omes", omes, append = T)

dbDisconnect(con)



