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
col <- brewer.pal(6, 'Set2') 

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

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')

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
metadata <- dbGetQuery(con, 'SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample and QC proteomics runs
rawfiles <- rawfiles %>%
  select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl('Sample|QC', run_type))

metadata <- metadata %>%
  select(-Sample) %>%
  mutate(sample_id = as.integer(sample_id))

df <- proteomics %>%
  left_join(rawfiles, by = 'rawfile_id') %>%
  left_join(metadata, by = 'sample_id')


## LIMIT TO NPA FOR NOW, AFTER RE-SEARCHING COMBINE NPs BY COMPLETENESS
# Also filter for completeness

df <- df %>%
  filter(grepl('NPA', rawfile_name)) #Select only sample runs

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
percent_non_na <- rowMeans(!is.na(spec_PG_NPs)) * 100
non_na_counts_ordered <- percent_non_na[order(-percent_non_na)]

sum(percent_non_na == 100) #1215
sum(percent_non_na >= 50) #6095

#find percent between 
p <- rowSums(!is.na(spec_PG_NPs[, -c(1, 2)]))
op <- p[order(-p)]
fiftyp <- op[1:7608]
divided <- 100*fiftyp/40
p.mv <- mean(divided)
print(p.mv)

missing <- data.frame(Order = seq_along(non_na_counts_ordered), Value = non_na_counts_ordered)
missing$color <- 'Y'

ggplot(missing, aes(Order, Value)) + 
  #geom_point(size = 3) + 
  geom_area(fill = col[8], alpha = 0.5) + # Shade the area under the line
  geom_line(color = col[3], linewidth = 2) + # Add the line on top
  geom_vline(xintercept = 1215, linetype = 'solid', color = 'black', linewidth = 1) +
  geom_vline(xintercept = 6095, linetype = 'solid', color = col[2], linewidth = 1) + 
  geom_hline(yintercept = 50, linetype = 'solid', color = col[2], linewidth = 1) +   
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
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave('reports/figures/AllPlates_Samples_PGvsMissing.pdf', width = 24, height = 16, units = 'cm')



#### plot:RunOrdervsProteinNumber ####
run_ids <- filtered_df_NA %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  left_join(filtered_df_NA %>% select(sample_id, rawfile_id, Sample, Cohort), by = 'sample_id') %>%
  distinct() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())
  
ggplot(run_ids, aes(order, count, color = Cohort)) + 
  geom_point() + 
  scale_color_manual(values = pal) +
  xlab('Run Order') +
  ylab('Filtered PG IDs') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7000)) +
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
erythrocyte <- c('P69905', 'P68871', 'P00915', 'P02042', 'P00918')
platelet <- c('P21333', 'Q9Y490', 'P35579', 'P60709', 'P18206')
coagulation <- c('P02675', 'P02679', 'P02671', 'P00488', 'P01008')

contaminant_df <- filtered_df %>%
  mutate(
    group = case_when(
      standardized_name %in% erythrocyte ~ "erythrocyte",
      standardized_name %in% platelet ~ "platelet",
      standardized_name %in% coagulation ~ "coagulation",
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
        legend.position = 'right',
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
        legend.position = 'right',
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
        legend.position = 'right',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

# Combine the plots with patchwork
p1 / p2 / p3
ggsave('reports/figures/AllPlates_Sample_Contaminants_Set_Point.pdf', 
       width = 16, height = 12, units = 'cm')


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
poi <- "Q9BXD5"

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
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
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
ggsave(paste0('reports/figures/AllPlates_Sample_singleprotein_', poi, '_distribution_Cohort.pdf'), 
              width = 8, height = 6, units = "cm")
