#### make a patient_metadata table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


#### sample_id table ####
# before importing into RStudio I changed the names of some of the original sample ID and cohort to fit with R. 
# To see the original can go to look at the P: projects folder 

patient_metadata <- read.csv("data/metadata/Coon_lab_SEER_sample_list_2024_metadata.csv")

patient_metadata <- tibble::rownames_to_column(patient_metadata, "sample_id")

colnames(patient_metadata)[colnames(patient_metadata) == "Sample.ID"] <- "Sample"


#drop old patient metadata table before adding new one
dbExecute(con, "DROP TABLE IF EXISTS patient_metadata")

# Write the new patient_metadata table to the database
dbWriteTable(con, "patient_metadata", patient_metadata, append = F, overwrite = T)

dbDisconnect(con)


## Fix wrong Sample numbers in rawfiles_all

# Establish a connection to the DB
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

df <- dbGetQuery(con, "SELECT Sample, sample_id
           FROM rawfiles_all
           ")

df_fixed <- df %>%
  mutate(Sample = case_when(
    Sample == '1032.2' ~ '1132.2',
    Sample == '1034.2' ~ '1134.2',
    Sample == '1039.2' ~ '1139.2',
    Sample == '1040.2' ~ '1140.2',
    Sample == '1050.2' ~ '1150.2',
    T ~ Sample
  )) 

## R Function to replace values 
replaceValues <- function(row) {
  dbExecute(
    con,
    "UPDATE rawfiles_all SET Sample = ? WHERE sample_id = ?",
    params = list(row["Sample"], row["sample_id"])
  )
}

##### Iterate over data frame to update normalized_abundance values in sqlite db ######

apply(df_fixed, 1, replaceValues)

dbDisconnect(con)





## Fix wrong Sample numbers in patient_metadata

# Establish a connection to the DB
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

df <- dbGetQuery(con, "SELECT Sample, sample_id
           FROM patient_metadata
           ")

df_fixed <- df %>%
  mutate(Sample = case_when(
    Sample == '1032.2' ~ '1132.2',
    Sample == '1034.2' ~ '1134.2',
    Sample == '1039.2' ~ '1139.2',
    Sample == '1040.2' ~ '1140.2',
    Sample == '1050.2' ~ '1150.2',
    T ~ Sample
  )) 

## R Function to replace values 
replaceValues <- function(row) {
  dbExecute(
    con,
    "UPDATE patient_metadata SET Sample = ? WHERE sample_id = ?",
    params = list(row["Sample"], row["sample_id"])
  )
}

##### Iterate over data frame to update normalized_abundance values in sqlite db ######

apply(df_fixed, 1, replaceValues)

dbDisconnect(con)
