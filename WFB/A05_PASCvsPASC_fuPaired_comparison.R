#### Overview ####

# 1. Volcano Plots of Effect Size after filtering
# 2. fgsea of the ranked effect size lists
# 3. Correlation of Qol and molecule abondance
# 4. Change in QoL correlated to change in biomolecule between the 2 timepoints


#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(viridis)
library(RSQLite)
library(ggrepel)
library(data.table)
library(fgsea)
library(gridExtra)
library(patchwork)
library(colorspace)


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



# plot colors
pie(rep(1, length(col)), col = col , main="") 


colorRampPalette(col)(15)


# files #
# Import patient metadata and biomolecule measurements
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                               FROM proteomics_measurement')
lipidomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                              FROM lipidomics_measurements')
transcriptomics <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, sample_id, Counts, normalized_counts
                                    FROM rnaseq_measurements")
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_id, rawfile_name, Sample, sample_id, run_type, ome_id, keep
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')
dbDisconnect(con)


# proteomics
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


# transcriptomics
df_t <- transcriptomics %>%
  filter(!is.na(sample_id)) %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_t <- biomolecules %>%
  filter(omics_id == 3) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_t <- df_t %>%
  filter(biomolecule_id %in% biomolecules_t) %>%
  rename(normalized_abundance = normalized_counts)


# combined
df_a <- df_p %>%
  mutate(ome = "p") %>%
  bind_rows(df_l %>%
              mutate(ome = "l")) %>%
  select(-rawfile_id, -rawfile_name, -run_type, -ome_id) %>%
  bind_rows(df_t %>%
              rename(raw_abundance = Counts,
                     normalized_abundance = normalized_counts) %>%
              mutate(ome = "t"))

filtered_df_a <- filtered_df_p %>%
  mutate(ome = "p") %>%
  bind_rows(filtered_df_l %>%
              mutate(ome = "l")) %>%
  select(-rawfile_id, -rawfile_name, -run_type, -ome_id) %>%
  bind_rows(filtered_df_t %>%
              rename(raw_abundance = Counts) %>%
              mutate(ome = "t"))


## Limit df to paired samples in the analysis
filtered_df_a_1 <- filtered_df_a %>%
  filter(analysis_group_2 == 1)



con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

pvalues <- dbGetQuery(con, 'SELECT biomolecule_id, analysis_group, test, comparison, formula, predictor, effect_size, eta_squared, lratio, p_value, q_value
                            FROM pvalues')
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, unique_patient_id, Collection_date, PG_change_collection_cutoff
                             FROM patient_metadata')
formulas <- dbGetQuery(con, 'SELECT *
                             FROM formula_table')

dbDisconnect(con)



## Paired T-Test Analysis ----
# options:
# analysis_group (1, 2, 3, 0)
anal <- 2
# comparison (Age, Sex, QoL, BMI)
comp <- "QoL"
# formula (1, 2, 3, etc.)
form <- 30
# ome
omea <- "transcript"


volc_plot <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 ~ "YES",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id)

counts <- volc_plot %>%
  filter(diffexp %in% c("YES", "NO")) %>%
  count(diffexp)


##* Volcano ----
ggplot(volc_plot,
       aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             stroke = 0.2) +
  scale_size_manual(values = c(0.5, 1.5), guide = "none") +
  scale_alpha_manual(values = c(0.2, 0.8), guide = "none") +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_x_continuous(limits = c(-max(abs(volc_plot$effect_size)), max(abs(volc_plot$effect_size)))) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab(paste("Effect Size", comp)) +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 
#+
#  geom_text(data = counts[counts$diffexp == "YES",], 
#            aes(x = -0.001, y = Inf, 
#                label = paste(n)),
#            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
#  geom_text(data = counts[counts$diffexp == "NO",], 
#            aes(x = 0.001, y = Inf, 
#                label = paste(n)),
#            hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) 
ggsave(paste0("reports/figures/Volcano_group_", anal, "_", comp, "_formula", form, "_LRM.pdf"), 
       width = 8, height = 6, units = "cm")



##* Boxplots of Significant Molecules ----

bm <- "PC 16:0_16:1_RTmz_20.63_732.55307"

single_bm_df <- filtered_df_a_1 %>%
  filter(standardized_name == bm)

ggplot(single_bm_df,
       aes(Cohort, normalized_abundance, color = Cohort, fill = Cohort)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
         scale_fill_manual(values = pal) +
         scale_color_manual(values = pal) +
         ggtitle(paste(bm, "Abundance")) +
         labs(x = NULL,
              y = "Log2 Abundance") +
         scale_y_continuous(expand = c(0,0), limits = c(min(single_bm_df$normalized_abundance) / 1.1, 
                                                        max(single_bm_df$normalized_abundance) * 1.1)) +
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
       ggsave(paste0('reports/figures/SingleProteinPlots/AnalysisGroup2_PASCvsPASC_fu_singlebm_', bm, '_distribution_Cohort.pdf'), 
              width = 8, height = 6, units = "cm")













## Correlation Analysis ----


       
##* Correlation Setup ----
## Calculate deltaQoL between each patient/timepoint
changes_qol <- filtered_df_a_1 %>%
  select(Cohort, SF.36.QOL.Score, unique_patient_id) %>%
  distinct() %>%
  pivot_wider(names_from = Cohort, values_from = SF.36.QOL.Score) %>%
  mutate(delta_qol = PASC_fu - PASC)

## Calculate deltaBiomolecule between each patient/timepoint
changes_biomolecule <- filtered_df_a_1 %>%
  select(biomolecule_id, Cohort, normalized_abundance, unique_patient_id) %>%
  pivot_wider(names_from = Cohort, values_from = normalized_abundance) %>%
  mutate(delta_biomolecule = PASC_fu - PASC)

changes_full <- changes_biomolecule %>%
  left_join(changes_qol,
            by = "unique_patient_id") %>%
  left_join(biomolecules %>% select(-keep),
              by = "biomolecule_id")

## Calculate Correlations ----

cor_results <- changes_full %>%
  group_by(biomolecule_id) %>%
  summarize(
    cor = cor(delta_biomolecule, delta_qol, use = "complete.obs"),
    pval = cor.test(delta_biomolecule, delta_qol)$p.value
  ) %>%
  mutate(qval = p.adjust(pval, method = "BH"))


## Plot:Correlation Scatter Plot ----

bm <- 38504

ggplot(changes_full %>% filter(biomolecule_id == bm), 
       aes(x = delta_qol, y = delta_biomolecule)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Change in QoL vs Change in Biomolecule")








# Plot:QoL vs Abundance with Cohort ----

bm <- 38504

ggplot(filtered_df_a_1 %>% filter(biomolecule_id == bm), 
       aes(x = normalized_abundance, y = SF.36.QOL.Score, color = Cohort, group = unique_patient_id)) +
  geom_line(color = "gray",
            linewidth = 0.2) +
  geom_point(size = 1) +
  #geom_smooth(method = "lm") +
  labs(title = "Change in QoL vs Change in Biomolecule")
ggsave(paste0('reports/figures/SingleProteinPlots/AnalysisGroup2_PASCvsPASC_fu_singlebm_', bm, '_scatter_Cohort.pdf'), 
       width = 8, height = 6, units = "cm")




## Plot:deltaQoL vs Abundance ----

changes_full_diff <- changes_full


bm <- 38504

ggplot(changes_full %>% filter(biomolecule_id == bm), 
       aes(x = delta_qol)) +
  geom_point(aes(y = PASC.x),
             color = "blue") +
  geom_point(aes(y = PASC_fu.x),
             color = "red") +
  #geom_smooth(method = "lm") +
  labs(title = "Change in QoL vs Change in Biomolecule")



