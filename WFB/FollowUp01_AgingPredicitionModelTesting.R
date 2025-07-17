
################## organize proteomics data for testing aging predictive model ######################


#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(RSQLite)
library(data.table)


# Colors #
col <- brewer.pal(8, "Dark2") 

pal1 <- c("#66C2A5",
          "#FFD92F",
          "#8DA0CB",
          "#FC8D62",
          "#A6CEE3",
          "#E78AC3",
          "#A6D854",
          "#FDB462",
          "#B3B3B3",
          "#B2DF8A")

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')


pal <- c("Acute" = "#E78AC3", 
         "Acute_fu" = '#AA3377', 
         "Acute_NC" = "#B3B3B3", 
         "Healthy" = '#229100', 
         "PASC" = '#66CCEE', 
         "PASC_fu" = '#4477AA')


col1 <- viridis_pal(option = "rocket")(100)[round(c(0.25, 0.5, 0.75) * 100)]


# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT *
                               FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT *
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                             FROM biomolecules_metadata')

dbDisconnect(con)



# proteomics
# Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles_p <- rawfiles %>%
  filter(ome_id == 1) %>%
  select(-keep) %>%
  filter(grepl('Sample', run_type))

df_p <- proteomics %>%
  inner_join(rawfiles_p, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_p <- biomolecules %>%
  filter(omics_id == 1) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_p <- df_p %>%
  filter(biomolecule_id %in% biomolecules_p) 




## make correct format for python lasso regression check ----
filtered_df_p_python <- filtered_df_p %>%
  filter(PG_change_collection_cutoff == 0) %>% # filter out post-tube change samples for strange behavior
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "gene_name") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id") %>% # add gene name for proper organization
  select(sample_id, Age, Cohort, metadata_value, normalized_abundance) %>%
  pivot_wider(id_cols = c(sample_id, Age, Cohort),
              names_from = metadata_value,
              values_from = normalized_abundance)

write_csv(filtered_df_p_python, "data/processed/ProteinsForAgingPredictionModel.csv")




## make correct format for coefficients ----
coef <- read_csv("data/processed/Lehallier2019_ModelCoefficients.csv") %>%
  mutate(Variable = word(Variable, 1, 1, "\\."))

write_csv(coef, "data/processed/Lehallier2019_ModelCoefficients_fixed.csv")

