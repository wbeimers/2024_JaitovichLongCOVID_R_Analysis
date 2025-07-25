library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(viridis)
library(RSQLite)
library(ggrepel)
library(data.table)
library(ROTS)
library(patchwork)
library(pheatmap)
library(broom)





## Make proteomics/lipidomics/transcriptomics dataframes ----

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


# lipidomics
rawfiles_l <- rawfiles %>%
  filter(ome_id == 2) %>%
  filter(keep == "1") %>%
  filter(run_type == "Sample") %>%
  select(-keep)

df_l <- lipidomics %>%
  inner_join(rawfiles_l, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_l <- biomolecules %>%
  filter(omics_id == 2) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_l <- df_l %>%
  filter(biomolecule_id %in% biomolecules_l)


# combined lipid and protein (no transcript in Acute samples)
df_pl <- df_p %>%
  mutate(ome = "p") %>%
  bind_rows(df_l %>%
              mutate(ome = "l"))

filtered_df_pl <- filtered_df_p %>%
  mutate(ome = "p") %>%
  bind_rows(filtered_df_l %>%
              mutate(ome = "l")) %>%
  filter(batch != 1) %>% # remove samples from Batch 1
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Acute") ~ "Acute_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID",
    Cohort %in% c("PASC") ~ "Long_COVID")) 






#### Acute -> PASC -> Acute_fu -> Healthy Differences ----
##* ANOVA ----
# start with filtered_df_pl for only proteins and lipids (due to acute), also no batch 1

# means
APAH_means <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id, Cohort) %>%
  summarise(mean_value = mean(normalized_abundance), .groups = "drop")

# anova and post-hoc tukey hsd for each biomolecule/comparison
APAH_comparison <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  do({
    model <- aov(normalized_abundance ~ Cohort, data = .)
    
    # Extract the ANOVA result
    tidy_result <- tidy(model)
    anova_pvalue <- tidy_result %>% 
      filter(term == "Cohort") %>% 
      pull(p.value)
    
    # Tukey HSD post-hoc test
    tukey <- TukeyHSD(model)
    tukey_df <- as.data.frame(tukey$Cohort)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$biomolecule_id <- unique(.$biomolecule_id)
    tukey_df$anova_pvalue <- anova_pvalue  # Add unadjusted p-value
    
    tukey_df
  }) %>%
  ungroup()

# anova_results <- filtered_df_pl %>%
#   filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
#   filter(PG_change_collection_cutoff == 0) %>%
#   group_by(biomolecule_id) %>%
#   do(tidy(aov(normalized_abundance ~ Cohort, data = .))) %>%
#   filter(term == "Cohort") 

# find the proportion of significant biomolecules for each comparison

APAH_signif <- APAH_comparison %>%
  group_by(comparison) %>%
  summarise(
    n_significant = sum(`p adj` < 0.05, na.rm = TRUE),
    n_total = n(),
    proportion_significant = n_significant / n_total
  ) %>%
  arrange(desc(proportion_significant))


## plot:padj histograms for each comparison
ggplot(APAH_comparison, aes(`p adj`)) + 
  geom_histogram(size = 0.2,
                 color = "black",
                 fill = "lightgray",
                 bins = 60) +
  xlab('P adjusted') +
  ylab('Count') +
  scale_y_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_text(size = 5),
        strip.background = element_rect(color = "black", size = 0.2),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  ) +
  facet_wrap(. ~ comparison)
ggsave('reports/figures/AllPlates_APAH_comparison_pvalue_hist.pdf', 
       width = 14, height = 8, units = 'cm')


##* Healthy vs Acute_fu ttest----
Acutefu_Healthy_comp <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  summarise(t_test = list(t.test(normalized_abundance ~ Cohort, data = cur_data()))) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result)


ggplot(Acutefu_Healthy_comp, aes(p.value)) + 
  geom_histogram(size = 0.2,
                 color = "black",
                 fill = "lightgray",
                 bins = 60) +
  xlab('P value') +
  ylab('Count') +
  scale_y_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_text(size = 5),
        strip.background = element_rect(color = "black", size = 0.2),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_AcutefuHealthy_comparison_pvalue_hist.pdf', 
       width = 6, height = 4, units = 'cm')



##* Healthy vs Acute ttest----
Acute_Healthy_comp <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  summarise(t_test = list(t.test(normalized_abundance ~ Cohort, data = cur_data()))) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result)


ggplot(Acute_Healthy_comp, aes(p.value)) + 
  geom_histogram(size = 0.2,
                 color = "black",
                 fill = "lightgray",
                 bins = 60) +
  xlab('P value') +
  ylab('Count') +
  scale_y_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_text(size = 5),
        strip.background = element_rect(color = "black", size = 0.2),
        legend.position = "none", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)
  )
ggsave('reports/figures/AllPlates_AcuteHealthy_comparison_pvalue_hist.pdf', 
       width = 6, height = 4, units = 'cm')



## Heatmap Equivalence ----
library(pheatmap)


expression_matrix <- filtered_df_a %>%
  filter(Cohort %in% c("Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  dplyr::select(biomolecule_id, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  tibble::column_to_rownames(var = "biomolecule_id") %>%
  select(where(~ !any(is.na(.))))

sample_annot <- filtered_df_a %>%
  ungroup() %>%
  select(sample_id, 
         Cohort, 
         #Age, 
         #Sex, 
         #BMI, 
         #SF.36.QOL.Score, 
         #PASC_Cohort, 
         #PG_change_collection_cutoff
  ) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") 

dat <- filtered_df_a %>%
  select(ome, biomolecule_id) %>%
  distinct()

row_annot <- data.frame(ome = dat$ome,
                        row.names = dat$biomolecule_id) %>%
  mutate(ome1 = case_when(
    ome == "p" ~ "protein",
    ome == "l" ~ "lipid",
    ome == "t" ~ "transcript"
  )) %>%
  select(-ome) %>%
  rename(Ome = ome1)

k <- 2

pheat <- pheatmap(t(expression_matrix),
                  #color = viridis(24, direction = 1, option = "plasma"),
                  color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(15)),
                  breaks = c(-4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4),
                  cluster_rows = T,
                  #cutree_rows = k, 
                  #gaps_row = T,
                  cluster_cols = T,
                  treeheight_row = 10,
                  treeheight_col = 0,
                  show_rownames = F,
                  show_colnames = F,
                  border_color = NA,
                  scale = "column",
                  annotation_row = sample_annot,
                  annotation_col = row_annot,
                  annotation_legend = F,
                  legend = F,
                  annotation_colors = list(Cohort = c(Acute_fu = pal[2], Healthy = pal[4]),
                                           Ome = c(protein = col[1], lipid = col[2], transcript = col[3])),
                  fontsize = 5,
                  fontsize_col = 5,
                  fontsize_row = 5,
                  filename = paste0("reports/figures/Heatmap_plots/Heatmap_AllOmes_AcutevsHealthy.png"),
                  width = 4.72,
                  height = 1.18)




