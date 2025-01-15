install.packages("DBI")
install.packages("RSQLite")

library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


#Readfiles

Lipid_metadata = read.csv("Lipidomics_Metadata_LipidSamplesOnly_Edited.csv", header = T)
raw_all = read.csv("rawfiles_all.csv", header = T)
metadata_database = read.csv("metadata_database.csv", header =T)

Lipid_metadata = merge(metadata_database, Lipid_metadata, by= "Sample" ,all= TRUE)



#Get Timestamp for all rawfiles
rawfiles_dir <- "H:/Coon Lab/Projects/4) LongCOVID/Used for analysis"

rawfiles <- list.files(rawfiles_dir, full.names = TRUE)

rawfile_data <- data.frame(
  rawfile_name = basename(rawfiles),  # Extract the file name only
  timestamp = as.character(file.info(rawfiles)$mtime)  # Get modification time with date and time
)

# Add the timestamp to the Lipid_metadata table
Lipid_metadata <- merge(Lipid_metadata, rawfile_data, by = "rawfile_name", all.x = TRUE)


write.csv(Lipid_metadata, "UseThisMetadata.csv", row.names = F) #Add Date and Time in new columns in excel

#Read it again
Lipid_metadata = read.csv("UseThisMetadata.csv", header = T)



#combine lipid metadata to proteomcis metadata for rawfiles_all
common_columns <- intersect(colnames(Lipid_metadata), colnames(raw_all))

Lipid_metadata_common <- Lipid_metadata[, common_columns]
raw_all_common <- raw_all[, common_columns]


combined_data <- rbind(raw_all_common, Lipid_metadata_common)

View(combined_data)

rawfiles_all = combined_data

write.csv(rawfiles_all, "rawfiles_all.csv", row.names = F)



#Addn to database
dbExecute(con, "DROP TABLE IF EXISTS rawfiles_all")

# Write the new patient_metadata table to the database (that deletes the older version and upload a new one)
dbWriteTable(con, "rawfiles_all", rawfiles_all, append = F, overwrite = T)


#That just adds new rows to whatever was in there
dbWriteTable(con, "rawfiles_all", rawfiles, append = T)



dbDisconnect(con)