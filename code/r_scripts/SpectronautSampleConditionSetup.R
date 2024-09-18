#Need to specify conditions going back to raw runs



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
details_samples <- NULL
for (i in raw_folders) {
  for (j in subfolders) {
    path <- paste0("P:/Projects/WFB_2024_Jaitovich_LongCOVID/LCMS_Runs/", i, "/", j, "/")
    files <- list.files(path = path,
                        pattern = "\\.raw.*",
                        recursive = T,
                        full.names = T)
    details_samples <- rbind(details_samples, file.info(files))
  }
}
details_samples <- tibble::rownames_to_column(details_samples, "names")
details_samples$Runs <- unname(sapply(details_samples$names, function(x) strsplit(x, "/")[[1]][7]))
details_samples$Runs <- sub(pattern = ".raw", replacement = "", x = details_samples$Runs)
#only looking at QCs and either NPA or NPB 
details_samples <- subset(details_samples, grepl("QC|NPB", details_samples$Runs))
details_samples$Sample <- sapply(strsplit(details_samples$Runs, "_"), function(x) paste(x[5], collapse = "_"))

#combine with groups metadata
groups <- read.csv("data/metadata/sample_groups.csv")

setup <- merge(details_samples, groups, by = "Sample", all.x = TRUE)

write.csv(setup, "data/metadata/spectronaut_setup_NPAandQC.csv")
