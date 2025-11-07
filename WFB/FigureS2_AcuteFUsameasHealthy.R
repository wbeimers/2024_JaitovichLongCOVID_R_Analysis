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
write_csv(APAH_signif, "data/processed/proteinlipid_group_significant_differences_ANOVA.csv")

## find the proportion of |log2(fc)| >= 1.5 for plot below
APAH_fc <- APAH_comparison %>%
  group_by(comparison) %>%
  summarise(
    n_fc = sum(abs(diff) >= 0.263, na.rm = TRUE),
    n_total = n(),
    proportion_fc = n_fc / n_total
  ) %>%
  arrange(desc(proportion_fc))
write_csv(APAH_fc, "data/processed/proteinlipid_group_fc_differences_ANOVA.csv")


## plot:padj histograms for each comparison ----
ggplot(APAH_comparison, aes(`p adj`)) + 
  geom_histogram(size = 0.2,
                 color = "black",
                 fill = "lightgray",
                 bins = 60) +
  geom_vline(xintercept = 0.05,
             color = "indianred",
             linetype = "dashed",
             linewidth = 0.2) +
  geom_text(data = APAH_signif,
            aes(x = 0,
                y = 4000,
                label = n_significant),
            size = 2) +
  geom_text(data = APAH_signif,
            aes(x = 0.1,
                y = 4000,
                label = (8701 - n_significant)),
            size = 2) +
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
        strip.text = element_text(size = 6),
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

## plot:log2 fc distributions for each comparison ----
ggplot(APAH_comparison, aes(`diff`)) + 
  geom_histogram(size = 0.2,
                 color = "black",
                 fill = "lightgray",
                 bins = 60) +
  geom_vline(xintercept = c(-0.263, 0.263),
             color = "indianred",
             linetype = "dashed",
             linewidth = 0.2) +
  geom_text(data = APAH_fc,
            aes(x = 2,
                y = 2000,
                label = n_fc),
            size = 2) +
  geom_text(data = APAH_fc,
            aes(x = 0,
                y = 2000,
                label = (8701 - n_fc)),
            size = 2) +
  xlab('Log2 Fold Change') +
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
        strip.text = element_text(size = 6),
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
ggsave('reports/figures/AllPlates_APAH_comparison_fc_hist.pdf', 
       width = 14, height = 8, units = 'cm')


## plot:umap of protein/lipid driven sample clustering ----
library(uwot)

filtered_df_pl_wide <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  select(biomolecule_id, sample_id, normalized_abundance) %>%
  pivot_wider(names_from = sample_id,
              values_from = normalized_abundance) %>%
  mutate(biomolecule_id = as.character(biomolecule_id)) %>% 
  select(where(~ !any(is.na(.)))) %>%
  t(.) %>%
  as.data.frame(.)
# turn row into colnames
colnames(filtered_df_pl_wide) <- as.character(filtered_df_pl_wide[1, ])
# Remove the row
filtered_df_pl_wide <- filtered_df_pl_wide[-1, ]
#change to numeric
filtered_df_pl_wide <- data.frame(lapply(filtered_df_pl_wide, as.numeric), row.names = rownames(filtered_df_pl_wide))
# re-add sample_id
filtered_df_pl_wide <- tibble::rownames_to_column(filtered_df_pl_wide, "sample_id")  



filtered_df_pl_umap <- as.data.frame(umap(filtered_df_pl_wide,
                                          min_dist = 0.05,
                                          n_neighbors = 10,
                                          dens_scale = 0.2)) %>%
  bind_cols(filtered_df_pl_wide %>%
              select(sample_id)) %>%
  mutate(sample_id = as.integer(sample_id)) %>%
  left_join(patient_metadata,
            by = "sample_id")
  
ggplot(filtered_df_pl_umap, 
       aes(V1, V2,
           fill = Cohort)) +
  geom_point(shape = 21,
             size = 2,
             alpha = 0.9,
             color = "black",
             stroke = 0.1) +
  #stat_ellipse(aes(color = Cohort), 
  #             geom = "path", 
  #             show.legend = FALSE,
  #             linewidth = 0.2) +
  #geom_text_repel(aes(label = Sample), size = 2) +
  scale_fill_manual(values = pal)  +
  #xlab(paste("PC1", round(variance$proportion[1]*100, 2))) +
  #ylab(paste("PC2", round(variance$proportion[2]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave(paste0("reports/figures/UMAP_pl_0p05mindist_10nn_0p2dens.pdf"), 
       width = 8, height = 6, units = "cm")




## Effect size / power analysis ----
library(effectsize)
library(pwr)
library(emmeans)
# anova and extract eta and f
APAH_pairwise_poweranalysis <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  do({
    model <- aov(normalized_abundance ~ Cohort, data = .)
    em <- emmeans(model, pairwise ~ Cohort)
    pw <- as.data.frame(em$contrasts)
    
    group_ns <- table(.$Cohort)
    
    # Extract group names from the contrast label
    pw <- pw %>%
      mutate(
        group1 = sub(" - .*", "", contrast),
        group2 = sub(".*- ", "", contrast)) %>%
      mutate(
        n1 = as.numeric(group_ns[group1]),
        n2 = as.numeric(group_ns[group2]),
        n_eff = 2 / (1 / n1 + 1 / n2),   # harmonic mean sample size
        d = t_to_d(t.ratio, df = df)$d,
      ) %>%
      rowwise() %>%
      mutate(
        power = pwr.t.test(
          d = abs(d),
          n = n_eff,
          sig.level = 0.05,
          type = "two.sample"
        )$power
      ) %>%
      ungroup()
  
    pw
  }) %>%
  ungroup()
write_csv(APAH_pairwise_poweranalysis,
          "data/processed/ANOVA_proteinlipid_poweranalysis.csv")

# get median effect sizes
power_med <- APAH_pairwise_poweranalysis %>%
  group_by(contrast) %>%
  summarize(
    median_d = median(abs(d), na.rm = TRUE),
    iqr_d = IQR(abs(d), na.rm = TRUE)
  )


ggplot(APAH_pairwise_poweranalysis, 
       aes(x = contrast,
           y = abs(d))) +
  geom_violin(width = 1,
              linewidth = 0.2) +
  geom_boxplot(width = 0.1,
               linewidth = 0.2,
               outliers = F) +
  labs(title = NULL, 
       y = "Effect Size (Cohen's d)", 
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7, angle = 60),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_text(size = 7),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave(paste0("reports/figures/effectsizedistribution_anovas_pl.pdf"), 
       width = 6, height = 6, units = "cm")

# check how many samples i would need for sufficient power
pwr.t.test(d = median(abs(APAH_pairwise_poweranalysis$d[APAH_pairwise_poweranalysis$contrast == "Acute_fu - PASC"])), 
           sig.level = 0.05, power = 0.8, type = "two.sample")



## pvalues eta explained/unexplained ----
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

pvalues <- dbGetQuery(con, 'SELECT *
                            FROM pvalues')
formulas <- dbGetQuery(con, 'SELECT *
                             FROM formula_table')
dbDisconnect(con)

# options:
# analysis_group (1, 2, 3, 0)
anal <- 7
# comparison (Age, Sex, QoL, BMI)
comp <- "group7_PASCnoPASC"
# formula (1, 2, 3, etc.)
form <- 57

volc_plot_PASC <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 & effect_size > 0 ~ "UP",
    q_value <= 0.05 & effect_size < 0 ~ "DOWN",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id)


ggplot(volc_plot_PASC, 
       aes(x = eta_squared)) +
  geom_histogram(binwidth = 0.01) +
  labs(title = NULL, 
       y = "Biomolecules", 
       x = "Eta Squared") +
  scale_y_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_text(size = 7),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave(paste0("reports/figures/effectsizedistribution_pl.pdf"), 
       width = 8, height = 6, units = "cm")

quantile(volc_plot_PASC$eta_squared, probs = 0.95)
quantile(volc_plot_PASC$eta_squared)




























##* Healthy vs PASC ttest----
Healthy_PASC_comp <- filtered_df_pl %>%
  filter(Cohort %in% c("PASC", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  summarise(t_test = list(t.test(normalized_abundance ~ Cohort, data = cur_data()))) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result) %>%
  mutate(contrast = "Healthy-PASC")


ggplot(Healthy_PASC_comp, aes(p.value)) + 
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


##* Acute_fu vs PASC ttest----
Acutefu_PASC_comp <- filtered_df_pl %>%
  filter(Cohort %in% c("PASC", "Acute_fu")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  summarise(t_test = list(t.test(normalized_abundance ~ Cohort, data = cur_data()))) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result) %>%
  mutate(contrast = "Acute_fu-PASC")


ggplot(Acutefu_PASC_comp, aes(p.value)) + 
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


##* Combined vs PASC ttest----
Combined_PASC_comp <- filtered_df_pl %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  mutate(group = case_when(
    Cohort %in% c("Healthy", "Acute_fu") ~ "control",
    Cohort == "PASC" ~ "PASC"
  )) %>%
  group_by(biomolecule_id) %>%
  summarise(t_test = list(t.test(normalized_abundance ~ group, data = cur_data()))) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result) %>%
  mutate(contrast = "Combined-PASC")


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


## * Checking equivalence from ttest ----
sample_sizes <- tibble(
  contrast = c("Acute_fu-PASC", "Healthy-PASC", "Combined-PASC"),
  n1 = c(15, 24, 39),  # control group sizes
  n2 = c(139, 139, 139)   # treatment group sizes
)

ttest_results_all <- Healthy_PASC_comp %>%
  bind_rows(Acutefu_PASC_comp) %>%
  bind_rows(Combined_PASC_comp) %>%
  left_join(sample_sizes, by = "contrast") %>%
  mutate(
    cohen_d = statistic / sqrt(1/n1 + 1/n2),
    abs_d = abs(cohen_d)
  )

ttest_results_all_wide <- ttest_results_all %>%
  select(biomolecule_id, contrast, cohen_d, p.value) %>%
  pivot_wider(names_from = contrast, values_from = c(cohen_d, p.value))

# plot similarity
ttest_results_all_wide %>%
  ggplot(aes(x = `cohen_d_Healthy-PASC`,
             y = `cohen_d_Acute_fu-PASC`)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Effect size (Healthy vs PASC)",
    y = "Effect size (Acute follow-up vs PASC)",
    title = "Concordance of Cohen's d between controls"
  )
ggsave('reports/figures/AcuteFUHealthy_comparison_ttest_effectsize_correlation.pdf', 
       width = 6, height = 4, units = 'cm')

# correlation
cor.test(
  ttest_results_all_wide$`cohen_d_Healthy-PASC`,
  ttest_results_all_wide$`cohen_d_Acute_fu-PASC`,
  use = "pairwise.complete.obs"
)

# significance overlap
ttest_results_all_wide %>%
  mutate(
    sig_C1 = `p.value_Healthy-PASC` < 0.05,
    sig_C2 = `p.value_Acute_fu-PASC` < 0.05
  ) %>%
  summarize(
    overlap = mean(sig_C1 & sig_C2),
    prop_sig_C1 = mean(sig_C1),
    prop_sig_C2 = mean(sig_C2)
  )







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




