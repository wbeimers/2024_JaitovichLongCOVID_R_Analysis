## Figure S1: Contamination Look ----

col1 <- brewer.pal(8, "Set1") 



##* All sample protein PCA ----

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
  filter(biomolecule_id %in% biomolecules_p) %>%
  filter(Cohort != "Acute_NC") %>%
  filter(Cohort != "PASC_fu")



pca_set <- filtered_df_p %>%
  select(standardized_name, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance)

t_pca_set <- as.data.frame(t(pca_set))

# turn row into colnames
colnames(t_pca_set) <- as.character(t_pca_set[1, ])
# Remove the row
t_pca_set <- t_pca_set[-1, ]
#change to numeric
t_pca_set <- data.frame(lapply(t_pca_set, as.numeric), row.names = rownames(t_pca_set))
# re-add sample_id
t_pca_set <- tibble::rownames_to_column(t_pca_set, "sample_id")
ncol(t_pca_set)

t_pca_set <- t_pca_set %>%
  mutate(sample_id = as.integer(sample_id)) 

pca_score <- prcomp(t_pca_set[,c(2:ncol(t_pca_set))],
                    scale. = T)
summary(pca_score)

# scree
explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
variance <- data.frame(proportion = explained_variance,
                       PC = 1:(ncol(pca_set)-1))

# pca scores
scores <- as.data.frame(pca_score$x)
scores <- scores %>%
  mutate(sample_id = t_pca_set$sample_id) %>%
  left_join(filtered_df_p %>% 
              ungroup() %>%
              select(sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, Collection_date, PG_change_collection_cutoff) %>%
              distinct(), 
            by = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score)) %>%
  mutate(Collection_date = as.integer(Collection_date))


ggplot(scores, aes(PC1, PC2, fill = Cohort)) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  #stat_ellipse(aes(color = Cohort), 
  #             geom = "path", 
  #             show.legend = FALSE,
  #             linewidth = 0.2) +
  #geom_text_repel(aes(label = Sample), size = 2) +
  scale_fill_manual(values = pal) +
  #scale_fill_viridis_c(option = "viridis", direction = -1) +
  xlab(paste("PC1", round(variance$proportion[1]*100, 2))) +
  ylab(paste("PC2", round(variance$proportion[2]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "inside", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave("reports/figures/Proteomics_AllPlates_Proteomics_Samples_Cohort_PCA_PC1PC2.pdf", 
       width = 6, height = 4, units = "cm")



##* PG ids vs run order or collection date ----

metadata_PGids <- metadata %>%
  filter(Cohort != "Acute_NC") %>%
  filter(Cohort != "PASC_fu") %>%
  left_join(rawfiles_p %>%
              select(sample_id, rawfile_id),
            by = "sample_id") %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())

ggplot(metadata_PGids, aes(order, PG_IDs, fill = Cohort)) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) + 
  scale_fill_manual(values = pal) +
  xlab('Run Order') +
  ylab('Filtered PG IDs') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Cohort')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "inside", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave('reports/figures/AllPlates_AllSamples_PGvsRunOrder_Cohort_Point.pdf', 
       width = 8, height = 4, units = 'cm')



## Contaminant Analysis ----
contaminant <- read.csv("data/metadata/contamination_proteins.csv") %>%
  mutate(single_id = word(Protein_IDs, 1, sep = ";")) %>%
  mutate(single_id = word(single_id, 1, sep = "-"))

# check for exact matches to contamination vector
contaminant1 <- contaminant %>%
  mutate(is = single_id %in% filtered_df_p$standardized_name) %>%
  mutate(single_id = case_when(
    single_id == "P0CG48" ~ "P62979",
    T ~ single_id
  )) %>% # this changes that id to the same as PG.ProteinGroups so there is a match
  mutate(is = single_id %in% filtered_df_p$standardized_name) %>%
  filter(is == T) %>%
  mutate(standardized_name = single_id) %>%
  dplyr::select(standardized_name, Type)


contaminant_df <- filtered_df_p %>%
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

##* Platelet analysis ----
platelet <- contaminant1 %>%
  filter(Type == "Platelet") %>%
  pull(standardized_name)

platelet_df <- filtered_df_p %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% platelet], na.rm = TRUE),
    ratio = sum_subset / sum_all
  ) %>%
  left_join(metadata_PGids, by = "sample_id")

platelet_df_mean_sd <- platelet_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(platelet_df, 
       aes(order, ratio, fill = factor(PG_change_collection_cutoff))) +
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) + 
  scale_fill_manual(values = col1) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio + (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio - (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("platelet") +
  xlab('Run Order') +
  ylab('sum(platelet proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.006, 0.0085)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_PlateletContamination.pdf', 
       width = 6, height = 3, units = 'cm')


# Correlation of contaminant index to protein ids
ggplot(platelet_df, aes(ratio, PG_IDs, fill = factor(PG_change_collection_cutoff))) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  scale_fill_manual(values = col1) +
  xlab('platelet contamination index (higher = more platelets)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_PlateletContamCount_cohort_Point.pdf', 
       width = 4, height = 3, units = 'cm')


## Erythrocyte analysis ----
erythrocyte <- contaminant1 %>%
  filter(Type == "Erythrocyte") %>%
  pull(standardized_name)

erythrocyte_df <- filtered_df_p %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% erythrocyte], na.rm = TRUE),
    ratio = sum_subset / sum_all
  ) %>%
  left_join(metadata_PGids, by = "sample_id")

erythrocyte_df_mean_sd <- erythrocyte_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(erythrocyte_df, 
       aes(order, ratio, fill = factor(PG_change_collection_cutoff))) +
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) + 
  scale_fill_manual(values = col1) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio + (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio - (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("erythrocyte") +
  xlab('Run Order') +
  ylab('sum(erythrocyte proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.0045, 0.007)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_ErythrocyteContamination.pdf', 
       width = 6, height = 3, units = 'cm')

# Correlation of contaminant index to protein ids
ggplot(erythrocyte_df, aes(ratio, PG_IDs, fill = factor(PG_change_collection_cutoff))) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  scale_fill_manual(values = col1) +
  xlab('erythrocyte contamination index (higher = more erythrocytes)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_ErythrocyteContamCount_cohort_Point.pdf', 
       width = 4, height = 3, units = 'cm')


## Coagulation analysis ----
coagulation <- contaminant1 %>%
  filter(Type == "Coagulation") %>%
  pull(standardized_name)

coagulation_df <- filtered_df_p %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[standardized_name %in% coagulation], na.rm = TRUE),
    ratio = sum_subset / sum_all
  ) %>%
  left_join(metadata_PGids, by = "sample_id")

coagulation_df_mean_sd <- coagulation_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(coagulation_df, 
       aes(order, ratio, fill = factor(PG_change_collection_cutoff))) +
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) + 
  scale_fill_manual(values = col1) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio + (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio - (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("coagulation") +
  xlab('Run Order') +
  ylab('sum(coagulation proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.00225, 0.003)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_CoagulationContamination.pdf', 
       width = 6, height = 3, units = 'cm')

# Correlation of contaminant index to protein ids
ggplot(coagulation_df, aes(ratio, PG_IDs, fill = factor(PG_change_collection_cutoff))) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  scale_fill_manual(values = col1) +
  xlab('coagulation contamination index (higher = more coagulation)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_Protein_CoagulationContamCount_cohort_Point.pdf', 
       width = 4, height = 3, units = 'cm')






# Separate plots for each facet
p1 <- ggplot(contaminant_df %>% filter(group == "erythrocyte"), 
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
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
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
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
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
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
ggsave('reports/figures/NeatRePrep_ContaminantProfiles.pdf', 
       width = 16, height = 12, units = 'cm')



















