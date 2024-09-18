####First set proper working directory under Session tab####
library(devtools)
library(tidyverse)

####SPECTRONAUT####
#read in default report
file <- "20240815_154253_20240815_WFB_JaitovichLongCOVID_AllPlates_SamplesQCs_DIA_NPA_Report.tsv"

spec <- read.delim(paste0("data/spectronaut_output/", file), sep = "\t")

####PGs####
#Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
selects <- (is.na(spec$PG.Qvalue) | spec$PG.Qvalue <= 0.01) &
  (is.na(spec$EG.Qvalue) | spec$EG.Qvalue <= 0.01) 

spec_raw <- spec[selects, c("R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Quantity")]

spec_raw %>%
  dplyr::group_by(PG.ProteinGroups, R.Condition) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)  

spec_PG <- pivot_wider(spec_raw,
                       names_from = R.FileName,
                       values_from = PG.Quantity,
                       values_fill = NA,
                       values_fn = list(PG.Quantity = mean))

colSums(!is.na(spec_PG))

#fixing columns so they all line up
spec_PG <- spec_PG %>%
  group_by(PG.ProteinGroups) %>%
  summarise_all(~ na.omit(.) %>% .[1]) %>%
  ungroup()

plot(colSums(!is.na(spec_PG)))

#save ungrouped pg quant matrix
write.csv(spec_PG, file = "data/processed/PG_Matrix_AllPlates_QCsSamples_NPA.csv")

#Make sure I have an NPA and NPB file since I searched separately
spec_PG_NPA <- spec_PG


#merge NPA and NPB dataframes
spec_PG <- merge(spec_PG_NPA, spec_PG_NPB, by = "PG.ProteinGroups", all = T)



#make dataframes just of samples, QCs, and controls
#controls
#controls_names <- spec_PG[, grepl(paste(c("PG.ProteinGroups", "PC1", "PC2", "MPE", "DC", "MSControl"), collapse = "|"), names(spec_PG))]
#spec_PG_controls <- data.frame(controls_names)
#write.csv(spec_PG_controls, file = "data/processed/PG_Matrix_Plate_01-05_controls.csv")

#QCs
QCs_names <- spec_PG[,grepl(paste(c("PG.ProteinGroups","QC"), collapse = "|"), names(spec_PG))]
spec_PG_QCs <- data.frame(QCs_names)
write.csv(spec_PG_QCs, file = "data/processed/PG_Matrix_AllPlates_QCs.csv")

#samples
samples_names <- spec_PG[,grepl(paste(c("PG.ProteinGroups", "NPA", "NPB"), collapse = "|"), names(spec_PG))]
spec_PG_samples <- data.frame(samples_names)
write.csv(spec_PG_samples, file = "data/processed/PG_Matrix_AllPlates_Samples.csv")


#use this to take either NPA or NPB for all samples for each pg, depending on which set is more complete
#make a filter to choose whichever row has fewer missing values for NPA or NPB. If same, choose NPA.
NPA <- apply(spec_PG_samples[,grep("NPA", names(spec_PG_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPA[is.na(NPA)] <-0

NPB <- apply(spec_PG_samples[,grep("NPB", names(spec_PG_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPB[is.na(NPB)] <-0

filter_NP <- NPA>=NPB
table(is.na(filter_NP))

#separate NPA and NPB into different dataframes
NPA_df <- spec_PG_samples[,c(1,grep("NPA", names(spec_PG_samples)))]
NPB_df <- spec_PG_samples[,c(1,grep("NPB", names(spec_PG_samples)))]

#apply filter to choose correct rows for either dataframe
NPA_df_filter <- NPA_df[filter_NP,]
NPB_df_filter <- NPB_df[!filter_NP,]
names(NPA_df_filter) <- sub(pattern = "_NPA", replacement = "", x = names(NPA_df_filter))
names(NPB_df_filter) <- sub(pattern = "_NPB", replacement = "", x = names(NPB_df_filter))
NPA_df_filter_1 <- NPA_df_filter[,c(1, order(names(NPA_df_filter)[-1])+1)]
NPB_df_filter_1 <- NPB_df_filter[,c(1, order(names(NPB_df_filter)[-1])+1)]

spec_PG_NPs <- rbind(NPA_df_filter, NPB_df_filter)
#save combined quant matrix
write.csv(spec_PG_NPs, file = "data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")

#make sure there are a certain number of values in each row
countrow <- apply(spec_PG_NPs[, -c(1, 2)], 1, function(x) table(unlist(unname(x)) > 0)[1])
non_NA <- countrow>398
spec_PG_NPs_100 <- spec_PG_NPs[non_NA,]
table(non_NA)

#save combined quant matrix
write.csv(spec_PG_NPs_100, file = "data/processed/PG_Matrix_AllPlates_Samples_NPs_100%.csv")

hist(countrow,
     breaks = 40)

colSums(!is.na(spec_PG))
colSums(!is.na(spec_PG_NPs))

plot(colSums(!is.na(spec_PG_NPs)))


####Peps####
#Pivot report into peptide quant values for each sample
selects <- (is.na(spec$PG.Qvalue) | spec$PG.Qvalue <= 0.01) &
  (is.na(spec$EG.Qvalue) | spec$EG.Qvalue <= 0.01) 

#stripped squence collapses precursors into a peptide group
spec_raw <- spec[selects, c("R.FileName","R.Condition", "EG.ModifiedSequence", "EG.TotalQuantity..Settings.")]

spec_raw %>%
  dplyr::group_by(EG.ModifiedSequence, R.Condition) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)  

spec_pep <- pivot_wider(spec_raw,
                       names_from = R.FileName,
                       values_from = EG.TotalQuantity..Settings.,
                       values_fill = NA,
                       values_fn = list(EG.TotalQuantity..Settings. = mean))

colSums(!is.na(spec_pep))

#fixing columns so they all line up
spec_pep <- spec_pep %>%
  group_by(EG.ModifiedSequence) %>%
  summarise_all(~ na.omit(.) %>% .[1]) %>%
  ungroup()

colSums(!is.na(spec_pep))

#save ungrouped pep quant matrix
write.csv(spec_pep, file = "data/processed/pep_Matrix_Plate_01-05_All.csv")


#make dataframes just of samples, QCs, and controls
#controls
controls_names <- spec_pep[, grepl(paste(c("EG.ModifiedSequence", "PC1", "PC2", "MPE", "DC", "MSControl"), collapse = "|"), names(spec_pep))]
spec_pep_controls <- data.frame(controls_names)
write.csv(spec_pep_controls, file = "data/processed/pep_Matrix_Plate_01-05_controls.csv")

#QCs
QCs_names <- spec_pep[,grepl(paste(c("EG.ModifiedSequence","QC"), collapse = "|"), names(spec_pep))]
spec_pep_QCs <- data.frame(QCs_names)
write.csv(spec_pep_QCs, file = "data/processed/pep_Matrix_Plate_01-05_QCs.csv")

#samples dataframe
samples_names <- spec_pep[,grepl(paste(c("EG.ModifiedSequence", "NPA", "NPB"), collapse = "|"), names(spec_pep))]
spec_pep_samples <- data.frame(samples_names)
write.csv(spec_pep_samples, file = "data/processed/pep_Matrix_Plate_01-05_Samples.csv")


#use this to take either NPA or NPB for all samples for each pg, depending on which set is more complete
#make a filter to choose whichever row has fewer missing values for NPA or NPB. If same, choose NPA.
NPA <- apply(spec_pep_samples[,grep("NPA", names(spec_pep_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPA[is.na(NPA)] <-0

NPB <- apply(spec_pep_samples[,grep("NPB", names(spec_pep_samples))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPB[is.na(NPB)] <-0

filter_NP <- NPA>=NPB
table(is.na(filter_NP))

#separate NPA and NPB into different dataframes
NPA_df <- spec_pep_samples[,c(1,grep("NPA", names(spec_pep_samples)))]
NPB_df <- spec_pep_samples[,c(1,grep("NPB", names(spec_pep_samples)))]

#apply filter to choose correct rows for either dataframe
NPA_df_filter <- NPA_df[filter_NP,]
NPB_df_filter <- NPB_df[!filter_NP,]
names(NPA_df_filter) <- sub(pattern = "_NPA", replacement = "", x = names(NPA_df_filter))
names(NPB_df_filter) <- sub(pattern = "_NPB", replacement = "", x = names(NPB_df_filter))
NPA_df_filter_1 <- NPA_df_filter[,c(1, order(names(NPA_df_filter)[-1])+1)]
NPB_df_filter_1 <- NPB_df_filter[,c(1, order(names(NPB_df_filter)[-1])+1)]

spec_pep_NPs <- rbind(NPA_df_filter, NPB_df_filter)
#save combined quant matrix
write.csv(spec_pep_NPs, file = "data/processed/pep_Matrix_Plate_01-05_Samples_NPs.csv")

#make sure there are a certain number of values in each row
countrow <- apply(spec_pep_NPs[, -c(1, 2)], 1, function(x) table(unlist(unname(x)) > 0)[1])
non_NA <- countrow>199
spec_pep_NPs_100 <- spec_pep_NPs[non_NA,]
table(non_NA)

#save combined quant matrix
write.csv(spec_pep_NPs_50, file = "data/processed/pep_Matrix_Plate_01-05_Samples_NPs_50%.csv")


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
#do this for looking at the combination of samples after combining npa and npb
details_samples <- subset(details_samples, grepl("NPA", details_samples$Runs))
details_samples <- details_samples[with(details_samples, order(as.POSIXct(mtime))),]
details_samples$order <- seq_along(details_samples$mtime)
details_samples$Samples <- sapply(strsplit(details_samples$Runs, "_"), function(x) paste(x[5], collapse = "_"))
samples <- details_samples$Samples
samples <- sub(pattern = "_NA", replacement = "", x = samples)
samples <- sub(pattern = "_NPA", replacement = "", x = samples)
samples <- sub(pattern = "-", replacement = ".", x = samples)
details_samples$Samples <- samples
write.csv(details_samples, "data/metadata/AllPlates_Samples_file_info.csv")
