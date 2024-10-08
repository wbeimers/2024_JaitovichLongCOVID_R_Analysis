#### make a patient_metadata table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


####sample_id table####
#before importing into RStudio I changed the names of some of the original sample ID and cohort to fit with R. To see the original
#can go to look at the P: projects folder 

patient_metadata <- read.csv("data/metadata/Coon_lab_SEER_sample_list_2024_metadata.csv")

patient_metadata <- tibble::rownames_to_column(patient_metadata, "sample_id")

colnames(patient_metadata)[colnames(patient_metadata) == "Sample.ID"] <- "Sample"

#fix the wrong database numbers
patient_metadata$Sample <- gsub("1132.2", "1032.2", patient_metadata$Sample)
patient_metadata$Sample <- gsub("1134.2", "1034.2", patient_metadata$Sample)
patient_metadata$Sample <- gsub("1139.2", "1039.2", patient_metadata$Sample)
patient_metadata$Sample <- gsub("1140.2", "1040.2", patient_metadata$Sample)
patient_metadata$Sample <- gsub("1150.2", "1050.2", patient_metadata$Sample)


#drop old patient metadata table before adding new one
dbExecute(con, "DROP TABLE IF EXISTS patient_metadata")

# Write the new patient_metadata table to the database
dbWriteTable(con, "patient_metadata", patient_metadata, append = F, overwrite = T)

dbDisconnect(con)



