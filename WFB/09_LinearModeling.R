#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(missForest)
library(VIM)
library(viridis)
library(RSQLite)
library(ggrepel)
library(limma)
library(broom)


# Colors #
# Make a classic palette
col <- brewer.pal(8, "Set2") 

pal <- c("#66C2A5",
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

# Make a Custom Gradient
col1 <- colorRampPalette(col)(16)

# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


proteomics <- dbGetQuery(con, "SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement")
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules")
rawfiles <- dbGetQuery(con, "SELECT rawfile_name, Sample, sample_id, ome_id, keep , rawfile_id, run_type
                           FROM rawfiles_all")
metadata <- dbGetQuery(con, "SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`
                           FROM patient_metadata")

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample and QC proteomics runs
rawfiles <- rawfiles %>%
  select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl("Sample|QC", run_type))

metadata <- metadata %>%
  select(-Sample) %>%
  mutate(sample_id = as.integer(sample_id))

df <- proteomics %>%
  left_join(rawfiles, by = "rawfile_id") %>%
  left_join(metadata, by = "sample_id")


## LIMIT TO NPA FOR NOW, AFTER RE-SEARCHING COMBINE NPs BY COMPLETENESS
# Also filter for completeness

df <- df %>%
  filter(grepl("NPA", rawfile_name)) #Select only sample runs

# Show how many non-NA values there are for each protein group in each study group
na_summary <- df %>%
  group_by(Cohort, standardized_name) %>%
  summarise(na_ratio = mean(!is.na(raw_abundance)), .groups = "drop")

# Make a list of IDs to keep where there are at least 50% non-NA values in one of the cohorts
ids_to_keep <- na_summary %>%
  group_by(standardized_name) %>%
  summarise(max_na_ratio = max(na_ratio)) %>%
  filter(max_na_ratio >= 0.5) %>% 
  pull(standardized_name)


filtered_df <- df %>%
  filter(standardized_name %in% ids_to_keep)



#### Linear Model ####
#LM results including study group, age, sex, and qol score
results <- filtered_df %>%
  group_by(standardized_name) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(normalized_abundance ~ Cohort + Age + Sex + BMI + `SF.36.QOL.Score`, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(adjusted_p_value = p.adjust(p.value, method = "fdr"))

filtered_results <- results %>%
  filter(adjusted_p_value < 0.05)
  

## limma
abundance_df <- filtered_df %>%
  select(sample_id, standardized_name, normalized_abundance) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  column_to_rownames("standardized_name")

metadata_df <- metadata %>%
  mutate(`SF.36.QOL.Score` = case_when(
    `SF.36.QOL.Score` == "N/A" ~ 0,
    T ~ as.numeric(`SF.36.QOL.Score`))) %>%
  mutate(Sex = factor(Sex, levels = c("M", "F"))) %>%
  mutate(Cohort = factor(Cohort, levels = c("Acute", "Acute_fu", "Acute_NC", "Healthy", "PASC", "PASC_fu"))) %>%
  mutate(Cohort = relevel(Cohort, ref = "Healthy"))

design <- model.matrix(~ Cohort + Age + Sex + BMI + `SF.36.QOL.Score`, data = metadata_df) 

fit <- lmFit(abundance_df, design)
fit<- eBayes(fit)
results <- topTable(fit, coef = "CohortAcute", number = Inf)



