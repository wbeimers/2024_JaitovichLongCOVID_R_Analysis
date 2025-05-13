
############################    *Load Libraries* ####################################

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)


############################    *Read Data* ####################################

data <- read.csv("D:/LongCOVID/Important Tables/20240919_Degreaser_Evaluated_ALL_plusUnknownFoundinMoreThan330_KAO20250122.csv", header= TRUE)
Metadata = read.csv("D:/LongCOVID/Metadata/Metadata_All_lipidomics.csv", header=TRUE)
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)

############################    *Edit Data* ####################################


#Add Columns for Mass Error Filter and Foundin100Files for Unknowns
data_edited <- data %>%
  mutate(MassError = ppmError > 10 | ppmError < -10) %>%  # TRUE if outside -10 to 10 range
  mutate(KeepID = as.logical(KeepID)) %>%  # Convert KeepID to logical
  mutate(FoundIn50perFiles = case_when(
    KeepID == FALSE & FoundInNFiles >= 167  ~ TRUE,   # If KeepID is FALSE and FoundInNFiles >= 85
    KeepID == FALSE & FoundInNFiles < 167   ~ FALSE,  # If KeepID is FALSE and FoundInNFiles < 85
    KeepID == TRUE  ~ NA                   # If KeepID is TRUE, assign NA
  )) %>%
  mutate(UniqueID = paste(RT, mz, sep = "_"))  # Create UniqueID by concatenating RT and mz


write.csv(data_edited, "D:/LongCOVID/Important Tables/20250207_longCOVID_AddedFilterationColumns.csv", row.names = FALSE) #edited the order of the columns in excel

data_edited <- read.csv("D:/LongCOVID/Important Tables/20250207_longCOVID_EditedFilterationColumns.csv", header= TRUE)


############################    *Filter Data* ####################################

#data_filtered <- data_edited %>% filter(MassError == FALSE, FoundIn50perFiles %in% c(TRUE, NA))

data_filtered <- data_edited %>% 
  filter(FoundIn50perFiles %in% c(TRUE, NA))

data_filtered = data_filtered %>%
mutate(
Annotation = if_else(KeepID == FALSE, NA_character_, Annotation),
LipidClass = if_else(KeepID == FALSE, NA_character_, LipidClass),
LipidSumComposition = if_else(KeepID == FALSE, NA_character_, LipidSumComposition)
)  

write.csv(data_filtered, "D:/LongCOVID/Important Tables/20250207_longCOVID_FilteredIDs.csv", row.names=FALSE)


data_quantifiable <- data_filtered %>%
  mutate(KeepQuantitation = as.logical(KeepQuantitation)) %>%
  filter(KeepQuantitation == TRUE)

Data_information = data_quantifiable[, c(1:31)]

write.csv(data_quantifiable,"D:/LongCOVID/Important Tables/20250207_longCOVID_FilteredandQuantifiableIDs.csv", row.names=FALSE)

####################################################