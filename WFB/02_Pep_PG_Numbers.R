#### Libraries/colors/files ####
# libraries 
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(RSQLite)


# colors
# Make a classic palette
col <- brewer.pal(8, 'Set2') 

# Make another palette
pal <- c('#66C2A5',
         '#FFD92F',
         '#8DA0CB',
         '#FC8D62',
         '#A6CEE3',
         '#E78AC3',
         '#A6D854',
         '#FDB462',
         '#B3B3B3',
         '#B2DF8A')

pal <- c('#EE6677', # Acute
         '#AA3377', # Acute_fu
         '#CCBB44', # Acute_NC
         '#228833', # Healthy
         '#66CCEE', # PASC
         '#4477AA') # PASC_fu

# Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),'white', colorRampPalette(col1)(100))

# plot colors
pie(rep(1, length(col)), col = col , main='') 


# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')


proteomics <- dbGetQuery(con, 'SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_name, Sample, sample_id, ome_id, keep , rawfile_id, run_type
                           FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles <- rawfiles %>%
  select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl('Sample', run_type))

metadata <- metadata %>%
  select(-Sample) %>%
  mutate(sample_id = as.integer(sample_id))

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
  select(-na_count, -keep_group, -NP)  
  
length(unique(df$standardized_name))

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
  filter(standardized_name %in% ids_to_keep) %>%
  filter(!is.na(normalized_abundance))



#### plot:missingness heatmap ----
filtered_df_NA <- filtered_df %>%
  mutate(na = ifelse(is.na(raw_abundance), 'NA', 'Non-NA')) %>%
  arrange(sample_id)

ggplot(filtered_df_NA, aes(sample_id, standardized_name, fill = na)) +
  geom_tile() + 
  scale_fill_manual(values = c('NA' = 'gray80', 'Non-NA' = 'steelblue')) +
  theme_minimal() +
  guides(fill = guide_legend(title = 'Missing')) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        plot.title = element_blank(), 
        axis.title = element_blank())
ggsave('reports/figures/AllPlates_Samples_MissingnessHeatmap_postfiltering.png', 
       width = 24, height = 16, units = 'cm')




#### plot:proteinvsmissing ####
#Calc percent of samples each protein was detected in
percent_non_na <- df %>%
  group_by(standardized_name) %>%  
  summarise(percent_non_na = mean(!is.na(raw_abundance)) * 100) %>%
  arrange(desc(percent_non_na)) %>%
  mutate(Order = seq(1,nrow(.)),
         color = "Y")

sum(percent_non_na$percent_non_na == 100) #1288
sum(percent_non_na$percent_non_na >= 50) #6166

ggplot(percent_non_na, aes(Order, percent_non_na)) + 
  #geom_point(size = 3) + 
  #geom_area(fill = col[8], alpha = 0.2) + # Shade the area under the line
  geom_line(color = col[3], linewidth = 0.5) + # Add the line on top
  geom_vline(xintercept = 1288, linetype = 'solid', color = 'black', linewidth = 0.2) +
  geom_vline(xintercept = 6166, linetype = 'solid', color = col[2], linewidth = 0.2) + 
  geom_hline(yintercept = 50, linetype = 'solid', color = col[2], linewidth = 0.2) +   
  #scale_fill_manual(values = col[5]) +
  ggtitle('Protein Group Counts') +
  xlab('Protein Groups') +
  ylab('% of Samples Quantified') +
  #ylim(0, 20) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,1)) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 101)) +
  guides(fill = guide_legend(title='Impairment\nStatus')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) 
ggsave('reports/figures/AllPlates_Samples_PGvsMissing.pdf', width = 16, height = 8, units = 'cm')



#### plot:RunOrdervsProteinNumber ####
run_ids <- df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  inner_join(df %>% ungroup() %>% select(sample_id, rawfile_id, Cohort), by = 'sample_id') %>%
  distinct() %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())
  
ggplot(run_ids, aes(order, count, color = Cohort)) + 
  geom_point() + 
  scale_color_manual(values = pal) +
  xlab('Run Order') +
  ylab('PG IDs') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Study Group')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(1, 0.05), 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
        ) 
ggsave('reports/figures/AllPlates_Sample_PGvsRunOrder_Set_Point.pdf', 
       width = 16, height = 6, units = 'cm')



#### plot:contaminations ####
# contaminant lists
contaminant <- read.csv("data/metadata/contamination_proteins.csv") %>%
  mutate(single_id = word(Protein_IDs, 1, sep = ";")) %>%
  mutate(single_id = word(single_id, 1, sep = "-"))


# check for exact matches to contamination vector
contaminant1 <- contaminant %>%
  mutate(is = single_id %in% filtered_df$standardized_name) %>%
  filter(is == T) %>%
  mutate(standardized_name = single_id) %>%
  dplyr::select(standardized_name, Type)

contaminant_df <- filtered_df %>%
  mutate(
    group = case_when(
      standardized_name %in% 
        (contaminant1 %>% 
           filter(Type == "Erythrocyte") %>% 
           pull(standardized_name)) ~ "erythrocyte",
      standardized_name %in% 
        (contaminant1 %>% 
           filter(Type == "Platelet") %>% 
           pull(standardized_name)) ~ "platelet",
      standardized_name %in% 
        (contaminant1 %>% 
           filter(Type == "Coagulation") %>% 
           pull(standardized_name)) ~ "coagulation",
      TRUE ~ "none" # Default group if no match is found
    )
  ) %>%
  filter(group != "none")

# Separate plots for each facet
p1 <- ggplot(contaminant_df %>% filter(group == "erythrocyte"), 
             aes(sample_id, normalized_abundance, color = standardized_name)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("erythrocyte") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

p2 <- ggplot(contaminant_df %>% filter(group == "platelet"), 
             aes(sample_id, normalized_abundance, color = standardized_name)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("platelet") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

p3 <- ggplot(contaminant_df %>% filter(group == "coagulation"), 
             aes(sample_id, normalized_abundance, color = standardized_name)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("coagulation") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

# Combine the plots with patchwork
p1 / p2 / p3
ggsave('reports/figures/AllPlates_Sample_Contaminants_Set_Point.pdf', 
       width = 16, height = 12, units = 'cm')

## Platelet analysis ----
platelet <- contaminant1 %>%
  filter(Type == "Platelet") %>%
  pull(standardized_name)

platelet_df <- filtered_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% platelet], na.rm = TRUE),
    ratio = sum_subset / sum_all
  ) %>%
  left_join(metadata, by = "sample_id")

platelet_df_mean_sd <- platelet_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(platelet_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = pal) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio + (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio - (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("platelet") +
  xlab('Sample ID') +
  ylab('sum(platelet proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.006, 0.009)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/AllSamples_PlateletContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_platelet <- run_ids %>%
  left_join(platelet_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_platelet, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = pal) +
  xlab('platelet contamination index (higher = more platelets)') +
  ylab('Seer Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/AllPlates_Samples_PlateletContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')


## Erythrocyte analysis ----
erythrocyte <- contaminant1 %>%
  filter(Type == "Erythrocyte") %>%
  pull(standardized_name)

erythrocyte_df <- filtered_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% erythrocyte], na.rm = TRUE),
    ratio = sum_subset / sum_all  
  ) %>%
  left_join(metadata, by = "sample_id")

erythrocyte_df_mean_sd <- erythrocyte_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(erythrocyte_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = pal) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio + (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio - (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("erythrocyte") +
  xlab('Sample ID') +
  ylab('sum(erythrocyte proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.0045, 0.007)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/AllSamples_ErythrocyteContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_erythrocyte <- run_ids %>%
  left_join(erythrocyte_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_erythrocyte, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = pal) +
  xlab('erythrocyte contamination index (higher = more erythrocyte)') +
  ylab('Seer Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/AllPlates_Samples_ErythrocyteContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')


## Coagulation analysis ----
coagulation <- contaminant1 %>%
  filter(Type == "Coagulation") %>%
  pull(standardized_name)

coagulation_df <- filtered_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% coagulation], na.rm = TRUE),
    ratio = sum_subset / sum_all 
  ) %>%
  left_join(metadata, by = "sample_id")

coagulation_df_mean_sd <- coagulation_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(coagulation_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = pal) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio + (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio - (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("coagulation") +
  xlab('Sample ID') +
  ylab('sum(coagulation proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.00225, 0.0032)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/AllSamples_CoagulationContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_coagulation <- run_ids %>%
  left_join(coagulation_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_coagulation, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = pal) +
  xlab('coagulation contamination index (higher = more coagulation proteins)') +
  ylab('Seer Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/AllPlates_Samples_CoagulationContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')




# plot:sample order log2 quant boxplots ----
quant_df <- filtered_df %>%
  select(rawfile_id, normalized_abundance, sample_id, Sample, Cohort) %>%
  arrange(sample_id) %>%
  mutate(order = dense_rank(sample_id)) %>%
  mutate(Cohort = as.factor(Cohort))


ggplot(quant_df, aes(factor(order), normalized_abundance, fill = Cohort, color = Cohort)) + 
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.1) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = NULL,
       y = "Log2 Abundance") +
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  guides(fill = guide_legend(title="Cohort")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size = 10),
        legend.position = 'bottom',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave(paste0('reports/figures/AllPlates_Sample_QuantBoxplots_Cohort.pdf'), 
       width = 16, height = 6, units = "cm")



## plot:single protein boxplots of cohorts ----
for (i in unique(filtered_df$standardized_name)) {
  
  poi <- i
  
  single_prot_df <- filtered_df %>%
    filter(grepl(poi, standardized_name))
  
  ggplot(single_prot_df, aes(Cohort, normalized_abundance, fill = Cohort, color = Cohort)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = NULL,
         y = "Log2 Abundance") +
    scale_y_continuous(expand = c(0,0), limits = c(min(single_prot_df$normalized_abundance) / 1.1, 
                                                   max(single_prot_df$normalized_abundance) * 1.1)) +
    theme_classic() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 7),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none"
    )
  
  ggsave(paste0('reports/figures/SingleProteinPlots/Proteomics_AllPlates_Sample_singleprotein_', poi, '_distribution_Cohort.pdf'), 
         width = 8, height = 6, units = "cm")
}



## plot:abundance rank vs missingness ----

df_ab_miss <- df %>%
  group_by(standardized_name) %>%  
  summarise(percent_non_na = mean(!is.na(raw_abundance)) * 100,
            median_abundance = max(normalized_abundance, na.rm = T)) %>%
  arrange(desc(median_abundance)) %>%
  mutate(Order = seq(1,nrow(.)),
         color = "Y")

ggplot(df_ab_miss, 
       aes(median_abundance, percent_non_na)) + 
  geom_point(shape = 16,
             size = 0.2,
             color = col[3]) +
  labs(x = "Max log2(Abundance)",
       y = "Percent non-NA") +
  scale_y_continuous(expand = c(0,0)) +
  #guides(fill = guide_legend(title="Cohort")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size = 10),
        legend.position = 'bottom',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave(paste0('reports/figures/AllPlates_Sample_MaxAbundancevsPercentComplete.pdf'), 
       width = 8, height = 6, units = "cm")
  
  
  
  #### plot:first/second PASC protein numbers ####
item <- run_ids %>%
  left_join(metadata %>% select(PASC_Cohort, sample_id), by = "sample_id") %>%
  filter(!is.na(PASC_Cohort))

  
ggplot(item, aes(PASC_Cohort, count, alpha = PASC_Cohort)) + 
  geom_boxplot(width = 0.4, 
               outliers = F,
               size = 0.2,
               color = "black",
               fill = pal[5]) +
  geom_jitter(width = 0.1, 
              size = 0.2,
              color = "black") +
  scale_alpha_manual(values = c(1, 0.5)) +
  labs(x = NULL,
       y = "Protein IDs") +
  scale_y_continuous(expand = c(0,0), limits = c(0,7500)) +
  guides(fill = guide_legend(title="Cohort")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size = 10),
        legend.position = 'bottom',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave(paste0('reports/figures/AllPlates_Sample_PASC_Cohort_ProteinIDs_boxplot.pdf'), 
       width = 8, height = 6, units = "cm")


t.test(count ~ PASC_Cohort, data = item)
