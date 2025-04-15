

############################  *Read Data* ####################################

data_samples <- read.csv("D:/LongCOVID/Important Tables/DataSubsets/20250207_longCOVID_AllSamples_RunOrder.csv", header= TRUE)
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)

############### Median Abundance #####################

filenames <- runorder %>%
  filter(runorder$run_type == "Sample") %>%
  pull(rawfile_name_R)  

data_samples <- data_samples%>%
  select(UniqueID, any_of(filenames))

############################  *Boxplot* ####################################

melted_data <- data_samples %>%
  pivot_longer(
    cols = -UniqueID,  # Keep UniqueID as an identifier, convert others to long format
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Batch = ifelse(grepl("Batch\\d+", variable, ignore.case = TRUE),
                   sub(".*(Batch\\d+).*", "\\1", variable),
                   NA)
  )

median_abundance <- melted_data %>%
  group_by(variable) %>%  # Group by sample (column "variable" contains sample names)
  summarize(median_abundance = median(log2(value), na.rm = TRUE)) %>%  # Compute median abundance
  arrange(median_abundance)  # Sort in ascending order

# View the samples with the lowest median abundance
write.csv(median_abundance, "D:/LongCOVID/Important Tables/MedianAbundanceTable_AllFeaturesFiltered.csv", row.names = FALSE)

#-----------------------------------------------------------------------