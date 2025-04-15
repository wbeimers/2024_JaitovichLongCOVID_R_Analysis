
############################  *Read Data* ####################################

runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)

######### Keep == 0 for Files Excluded #################

runorder <- read.csv("D:/LongCOVID/Metadata/RunOrderMetadata.csv", header= TRUE)

Samples_excluded = c("X20240808_SIA_Batch2_36.raw", 
                     "X20240729_SIA_Batch6_16.raw", 
                     "X20240729_SIA_Batch6_58.raw",
                     "X20240801_SIA_Batch7_2.raw", 
                     "X20240806_SIA_Batch10_102.raw",  
                     "X20240809_SIA_Batch1_1030.raw" )


runorder_v2 <- runorder %>%
  mutate(
    keep = ifelse(rawfile_name_R %in% Samples_excluded, 0, 1)  
  )


write.csv(runorder_v2, "D:/LongCOVID/Metadata/rawfiles_Metadata.csv", row.names = FALSE)

#---------------------------------------------------------------------