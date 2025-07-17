
############################    *Load Libraries* ####################################

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)


############################    *Read Data* ####################################
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)
data_filtered <- read.csv("D:/LongCOVID/Important Tables/20250207_longCOVID_FilteredIDs.csv", header= TRUE)

Data_information = data_filtered[,c(1:31)]

############################    *QC only Data* ####################################

qc_filenames <- runorder %>%
  filter(runorder$run_type == "QcRep") %>%
  pull(rawfile_name_R)  

QCs <- data_filtered[, qc_filenames, drop = FALSE]
QC_expression_log2 <- log2(QCs)
data_QC = cbind(Data_information, QCs)

write.csv(data_QC, "D:/LongCOVID/Important Tables/DataSubsets/20250207_longCOVID_AllQCs_RunOrder.csv", row.names=FALSE)

############################ *Samples only Data - No Blanks, QCs nor Dilution* ####################################
sample_filenames <- runorder %>%
  filter(Metadata$run_type == "Sample") %>%
  pull(rawfile_name_R)  

study_samples <- data_filtered[, sample_filenames, drop = FALSE]
Expression =study_samples
data_samples = cbind(Data_information, study_samples)

write.csv(data_samples, "D:/LongCOVID/Important Tables/DataSubsets/20250207_longCOVID_AllSamples_RunOrder.csv", row.names=FALSE)

###################################################################