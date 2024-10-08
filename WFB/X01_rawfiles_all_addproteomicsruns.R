#### make a rawfiles_all table for the long covid database


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


####raw file info####
#find info on run order for each of the QC runs and make a new dataframe combining them
raw_folders <- list("20240612_Plate_01",
                    "20240615_Plate_02",
                    "20240620_Plate_03",
                    "20240626_Plate_04",
                    "20240628_Plate_05",
                    "20240704_Plate_06",
                    "20240708_Plate_07",
                    "20240726_Plate_08_rerun",
                    "20240716_Plate_09",
                    "20240721_Plate_10"
)
subfolders <- list("Samples", "QCs")
rawfiles <- NULL
for (i in raw_folders) {
  for (j in subfolders) {
    path <- paste0("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Proteomics/LCMS_Runs/", i, "/", j, "/")
    files <- list.files(path = path,
                        pattern = "\\.raw.*",
                        recursive = T,
                        full.names = T)
    rawfiles <- rbind(rawfiles, file.info(files))
  }
}
#turn the rows into a column
rawfiles <- tibble::rownames_to_column(rawfiles, "names")

#make the rawfile_name column of the name of each rawfile
rawfiles$rawfile_name <- unname(sapply(rawfiles$names, function(x) strsplit(x, "/")[[1]][8]))
rawfiles$rawfile_name <- sub(pattern = ".raw", replacement = "", x = rawfiles$rawfile_name)

#rename mtime to timestamp
colnames(rawfiles)[colnames(rawfiles) == "mtime"] <- "timestamp"

#make a rawfile order by runtime
rawfiles <- rawfiles[with(rawfiles, order(as.POSIXct(timestamp))),]
rawfiles$rawfile_id <- seq_along(rawfiles$timestamp)

#remove unnecessary columns
rawfiles <- subset(rawfiles, select = -c(names, size, isdir, mode, ctime, atime, exe))

#add a column for type of run (blank, control, samples)
rawfiles$run_type <- ifelse(grepl("NPA|NPB", rawfiles$rawfile_name), "Sample",
                        ifelse(grepl("QC", rawfiles$rawfile_name), "QC",
                               ifelse(grepl("PC1|DC|MPE|PC2|MSControl", rawfiles$rawfile_name), "Control", 
                                      ifelse(grepl("Blank", rawfiles$rawfile_name), "Blank", NA))))

#add a column to keep or not (1 = keep, 0 = not)
rawfiles$keep <- ifelse(grepl("Blank", rawfiles$run_type), 0, 1)

#add a column for batch number
rawfiles$batch <- sapply(strsplit(rawfiles$rawfile_name, "_"), function(x) paste(x[4], collapse = "_"))

#add a column for the ome_id
rawfiles$ome_id <- 1

#add a column for sample to go back to other samples
rawfiles$Sample <- sapply(strsplit(rawfiles$rawfile_name, "_"), function(x) paste(x[5], collapse = "_"))
rawfiles$Sample <- sub(pattern = "-", replacement = ".", x = rawfiles$Sample)

#now make sample_id to match (or can use the sample to match??)
rawfiles <- merge(rawfiles, patient_metadata[, c("Sample", "sample_id")], by = "Sample", all.x = T)
rawfiles$sample_id[is.na(rawfiles$sample_id)] <- -1


dbWriteTable(con, "rawfiles_all", rawfiles, append = T)

dbDisconnect(con)
