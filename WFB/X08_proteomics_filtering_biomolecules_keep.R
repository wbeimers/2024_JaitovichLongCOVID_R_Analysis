



## Figure out how many IDs to keep and what their biomolecule_id is so we can filter easily with the biomolecules table ----
# 1. Combine NPs by greater presence
# 2. look for 50% presence of the protein in one of the cohorts and only keep those that pass
# 3. Go back and only make normalized data only for those. (only do this if the 50% cutoff changes)






# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_name, Sample, sample_id, ome_id, rawfile_id, run_type
                           FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles <- rawfiles %>%
  filter(ome_id == 1) %>%
  filter(grepl('Sample', run_type))

df <- proteomics %>%
  left_join(rawfiles, by = 'rawfile_id') %>%
  left_join(metadata, by = 'sample_id')



# Combine by which NP has more completeness by protein group. One NP for each protein group
df <- df %>%
  mutate(NP = case_when(
    grepl("NPA", rawfile_name) == T ~ "NPA",
    grepl("NPB", rawfile_name) == T ~ "NPB"
  )) %>%
  group_by(standardized_name, NP) %>%
  mutate(na_count = sum(is.na(raw_abundance))) %>%
  ungroup() %>%
  group_by(standardized_name) %>%
  mutate(keep_group = NP[which.min(na_count)]) %>%  
  filter(NP == keep_group) %>%  
  dplyr::select(-na_count, -keep_group, -NP)  



# Also filter for completeness
# Show how many non-NA values there are for each protein group in each study group
na_summary <- df %>%
  group_by(Cohort, standardized_name) %>%
  summarise(na_ratio = mean(!is.na(raw_abundance)), .groups = 'drop')

# Make a list of IDs to keep where there are at least 50% non-NA values in one of the cohorts
ids_to_keep <- na_summary %>%
  group_by(standardized_name) %>%
  summarise(max_na_ratio = max(na_ratio)) %>%
  filter(max_na_ratio >= 0.5) %>% 
  pull(standardized_name)

filtered_df <- df %>%
  filter(standardized_name %in% ids_to_keep) 

biomolecules_to_keep <- unique(filtered_df$biomolecule_id)

# take the biomolecule IDs, and make those 1s in the biomolecules table, and the rest change to 0

biomolecules1 <- biomolecules %>%
  mutate(keep = case_when(
    omics_id == "1" & biomolecule_id %in% biomolecules_to_keep ~ "1",
    omics_id == "1" & !(biomolecule_id %in% biomolecules_to_keep) ~ "0",
    TRUE ~ keep  
  ))


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "biomolecules", biomolecules1, overwrite = T)
dbDisconnect(con)



