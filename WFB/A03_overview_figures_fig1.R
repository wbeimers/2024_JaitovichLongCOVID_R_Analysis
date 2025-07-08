#### Overview ####

# 1. Plots of metadata for overview figure



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
library(ggpubr)


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

tab <- c("#1f77b4",
         "#ff7f0e",
         "#2ca02c",
         "#d62728",
         "#9467bd",
         "#8c564b",
         "#e377c2",
         "#7f7f7f",
         "#bcbd22",
         "#17becf")



# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')

dbDisconnect(con)



## plot:pie - all samples ----
counts <- metadata %>%
  count(Cohort, name = "count") %>%
  arrange(desc(count)) %>% 
  mutate(Cohort = factor(Cohort, levels = rev(Cohort))) %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Cohort, "\n", count)  
  )
  


ggplot(counts, aes(x = "", y = count, fill = Cohort)) +
  geom_bar(stat = "identity", 
           width = 0.99,
           color = "black",
           linewidth = 0.2,
           linetype = "solid") + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  geom_text(aes(y = ypos, label = label), 
            color = "black", 
            size = 2, 
            nudge_x = 0.3) +
  theme_void() +
  theme(legend.position = "none") +
  expand_limits(y = c(0, counts$count))
ggsave(paste0("reports/figures/pie_allsamples.pdf"), 
       width = 6, height = 6, units = "cm")



## plot:pie - Analysis Group 1 ----
counts <- metadata %>%
  filter(analysis_group_1 == 1) %>%
  count(Cohort, name = "count") %>%
  bind_rows(data.frame(Cohort = "Excluded", count = 400 - sum(counts$count))) %>%
  arrange(desc(count)) %>% 
  mutate(Cohort = factor(Cohort, levels = rev(Cohort))) %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Cohort, "\n", count)  
  )


ggplot(counts, aes(x = "", y = count, fill = Cohort)) +
  geom_bar(stat = "identity", 
           width = 0.99,
           color = "black",
           linewidth = 0.2,
           linetype = "solid") + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  geom_text(aes(y = ypos, label = label), 
            color = "black", 
            size = 2, 
            nudge_x = 0.3) +
  theme_void() +
  theme(legend.position = "none") +
  expand_limits(y = c(0, counts$count))
ggsave(paste0("reports/figures/pie_analysisgroup1.pdf"), 
       width = 6, height = 6, units = "cm")



## plot:pie - Analysis Group 2 ----
counts <- metadata %>%
  filter(analysis_group_2 == 1) %>%
  count(Cohort, name = "count") %>%
  bind_rows(data.frame(Cohort = "Excluded", count = 400 - sum(counts$count))) %>%
  arrange(desc(count)) %>% 
  mutate(Cohort = factor(Cohort, levels = rev(Cohort))) %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Cohort, "\n", count)  
  )


ggplot(counts, aes(x = "", y = count, fill = Cohort)) +
  geom_bar(stat = "identity", 
           width = 0.99,
           color = "black",
           linewidth = 0.2,
           linetype = "solid") + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  geom_text(aes(y = ypos, label = label), 
            color = "black", 
            size = 2, 
            nudge_x = 0.3) +
  theme_void() +
  theme(legend.position = "none") +
  expand_limits(y = c(0, counts$count))
ggsave(paste0("reports/figures/pie_analysisgroup2.pdf"), 
       width = 6, height = 6, units = "cm")



## plot:pie - Analysis Group 3 ----
counts <- metadata %>%
  filter(analysis_group_3 == 1) %>%
  count(Cohort, name = "count") %>%
  bind_rows(data.frame(Cohort = "Excluded", count = 400 - sum(counts$count))) %>%
  arrange(desc(count)) %>% 
  mutate(Cohort = factor(Cohort, levels = rev(Cohort))) %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Cohort, "\n", count)  
  )


ggplot(counts, aes(x = "", y = count, fill = Cohort)) +
  geom_bar(stat = "identity", 
           width = 0.99,
           color = "black",
           linewidth = 0.2,
           linetype = "solid") + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  geom_text(aes(y = ypos, label = label), 
            color = "black", 
            size = 2, 
            nudge_x = 0.3) +
  theme_void() +
  theme(legend.position = "none") +
  expand_limits(y = c(0, counts$count))
ggsave(paste0("reports/figures/pie_analysisgroup3.pdf"), 
       width = 6, height = 6, units = "cm")


## plot:bar - all samples ----
counts <- metadata %>%
  count(Cohort, name = "count") %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Cohort, "\n", count)  
  )

ggplot(counts, aes(x = reorder(Cohort, desc(count)), y = count, fill = Cohort)) +
  geom_bar(stat = "identity", 
           width = 0.9) + 
  #coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  geom_text(aes(y = count, label = count), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  labs(y = "Samples") +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
              panel.grid.major = element_blank(), 
              axis.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_blank(), 
              axis.title = element_text(size = 6),
              legend.title = element_text(size = 6),
              legend.text = element_text(size = 6),
              axis.line = element_line(linewidth = 0.2),
              axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bar_allsamples.pdf"), 
       width = 6, height = 4, units = "cm")


## plot:bar - Analysis Group 1 ----
counts <- metadata %>%
  count(Cohort, analysis_group_1) %>%
  group_by(Cohort) %>%
  mutate(
    total = sum(n), 
    ypos = cumsum(n) - n / 2
  )

ggplot(counts, aes(x = reorder(Cohort, desc(total)), y = n, fill = Cohort, alpha = as.factor(analysis_group_1))) +
  geom_bar(stat = "identity",
           position = "stack", 
           width = 0.9) + 
  #coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  scale_alpha_manual(values = c(0.2, 1)) +
  geom_text(aes(y = ypos, label = n), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  labs(y = "Samples") +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_blank(), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bar_analysisgroup1.pdf"), 
       width = 6, height = 4, units = "cm")


## plot:bar - Analysis Group 2 ----
counts <- metadata %>%
  count(Cohort, analysis_group_4) %>%
  group_by(Cohort) %>%
  mutate(
    total = sum(n), 
    ypos = cumsum(n) - n / 2
  )

ggplot(counts, aes(x = reorder(Cohort, desc(total)), y = n, fill = Cohort, alpha = as.factor(analysis_group_4))) +
  geom_bar(stat = "identity",
           position = "stack", 
           width = 0.9) + 
  #coord_polar(theta = "y") + 
  scale_fill_manual(values = pal) +
  scale_alpha_manual(values = c(0.2, 1)) +
  geom_text(aes(y = ypos, label = n), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  labs(y = "Samples") +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_blank(), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bar_analysisgroup4.pdf"), 
       width = 6, height = 4, units = "cm")



## plot:bmi boxplot - all samples ----
ggplot(metadata, aes(x = "", y = BMI)) +
  geom_violin(alpha = 0.5, 
              width = 0.8, 
              size = 0.2,
              color = "black",
              fill = "gray") +
  geom_boxplot(width = 0.2, 
               alpha = 0.8, 
               outliers = F,
               size = 0.2,
               fill = "gray") +
  #geom_text(aes(y = count, label = count), 
  #          color = "black", 
  #          size = 2, 
  #          nudge_y = -10) +
  labs(y = "BMI") +
  #scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bmi_boxplot_allsamples.pdf"), 
       width = 3, height = 4, units = "cm")



## plot:age boxplot - all samples ----
ggplot(metadata, aes(x = "", y = Age)) +
  geom_violin(alpha = 0.5, 
              width = 0.8, 
              size = 0.2,
              color = "black",
              fill = "gray") +
  geom_boxplot(width = 0.2, 
               alpha = 0.8, 
               outliers = F,
               size = 0.2,
               fill = "gray") +
  #geom_text(aes(y = count, label = count), 
  #          color = "black", 
  #          size = 2, 
  #          nudge_y = -10) +
  labs(y = "Age") +
  #scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/age_boxplot_allsamples.pdf"), 
       width = 3, height = 4, units = "cm")



## plot:sex bar - all samples ----
counts <- metadata %>%
  count(Sex, name = "count") %>%
  mutate(
    cumulative = cumsum(count), 
    ypos = cumulative - count / 2,
    label = paste0(Sex, "\n", count)  
  )

ggplot(counts, aes(x = Sex, y = count)) +
  geom_bar(stat = "identity", 
           width = 0.9,
           fill = "gray") + 
  #coord_polar(theta = "y") + 
  #scale_fill_manual(values = col) +
  geom_text(aes(y = count, label = count), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  labs(y = "Samples",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bar_sex_allsamples.pdf"), 
       width = 4, height = 4, units = "cm")



## plot:qol boxplot - all samples ----
ggplot(metadata, aes(x = "", y = SF.36.QOL.Score)) +
  geom_violin(alpha = 0.5, 
              width = 0.8, 
              size = 0.2,
              color = "black",
              fill = "gray") +
  geom_boxplot(width = 0.2, 
               alpha = 0.8, 
               outliers = F,
               size = 0.2,
               fill = "gray") +
  annotate("text",x = "", y = 910, label = "n = 236", size = 2) +
  #geom_text(aes(y = count, label = count), 
  #          color = "black", 
  #          size = 2, 
  #          nudge_y = -10) +
  labs(y = "QoL Score") +
  scale_y_continuous(expand = c(0,0), limits = c(0,950)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/QoL_boxplot_allsamples.pdf"), 
       width = 3, height = 4, units = "cm")



## Biomolecule Data Overview ----
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



## Plot:bar - biomolecule overview ----
bm_counts <- df_a %>%
  group_by(ome) %>%
  summarise(biomolecules = n_distinct(standardized_name)) %>%
  mutate(filtered = "N")

bm_counts_filtered <- filtered_df_a %>%
  group_by(ome) %>%
  summarise(biomolecules = n_distinct(standardized_name)) %>%
  mutate(filtered = "Y")

bm_counts_all <- bind_rows(bm_counts, bm_counts_filtered) %>%
  mutate(ome = factor(ome, levels = c("p", "l", "t")))


ggplot(bm_counts_all, aes(x = ome, y = biomolecules, alpha = filtered, fill = ome)) +
  geom_bar(stat = "identity", 
           width = 0.8,
           position = position_dodge()) + 
  #coord_polar(theta = "y") + 
  scale_fill_manual(values = col) +
  scale_alpha_manual(values = c(0.4, 1)) +
  guides(fill = "none") +
  geom_text(aes(y = biomolecules, label = biomolecules), 
            color = "black", 
            size = 2, 
            position = position_dodge(width = 0.8),
            vjust = 1.2) +
  labs(y = "# IDs",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20000)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"))
ggsave(paste0("reports/figures/bar_biomoleculeIDS_allsamples.pdf"), 
       width = 6, height = 4, units = "cm")


## Plot:pie - biomolecule filtered ----

bm_counts_filtered <- filtered_df_a %>%
  group_by(ome) %>%
  summarise(biomolecules = n_distinct(standardized_name)) %>%
  mutate(filtered = "Y") %>%
  arrange(desc(biomolecules)) %>%
  mutate(
    proportion = biomolecules / sum(biomolecules),
    ypos = cumsum(proportion) - 0.5 * proportion,
    label = paste0(ome, ": ", biomolecules)
  )

ggplot(bm_counts_filtered, aes(x = "", y = proportion, fill = ome)) +
  geom_bar(stat = "identity", 
           width = 1,
           color = "white") + 
  coord_polar("y") + 
  scale_fill_manual(values = c(col[2], col[1], col[3])) +
  guides(fill = "none") +
  geom_text(aes(y = ypos, label = biomolecules), 
            color = "black", 
            size = 2) +
  labs(y = NULL,
       x = NULL) +
  theme_void() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        plot.title = element_blank(), 
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"))
ggsave(paste0("reports/figures/pie_biomoleculeIDS_filtered.pdf"), 
       width = 3, height = 3, units = "cm")


## Plot:connection plot - QoL Paired ----
paired_metadata <- metadata %>%
  filter(analysis_group_2 == 1)


ggplot(paired_metadata, aes(x = Cohort, y = SF.36.QOL.Score, group = unique_patient_id)) +
  geom_line(color = "gray",
            size = 0.5) + 
  geom_point(size = 1.5,
             alpha = 0.8) + 
  labs(y = "QoL",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 900)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/point_QoL_change_paired.pdf"), 
       width = 6, height = 4, units = "cm")

# all paried samples (including tube change ones)
paired_metadata1 <- metadata %>%
  filter(!is.na(Paired_samples)) %>%
  filter(!grepl("Acute", Cohort)) %>%
  filter(unique_patient_id != "P1134")

ggplot(paired_metadata1, aes(x = Cohort, y = SF.36.QOL.Score, group = unique_patient_id)) +
  geom_line(color = "gray",
            size = 0.5,
            alpha = 0.8) + 
  geom_point(size = 1.5,
             alpha = 0.8) + 
  labs(y = "QoL",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 900)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_blank(), 
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/point_QoL_change_paired_allsamples.pdf"), 
       width = 6, height = 4, units = "cm")



## plot:bmi boxplot - analysis group 1 ----

ggplot(metadata %>%
         filter(analysis_group_1 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("PASC", "PASC_fu") ~ "Long_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
         distinct(unique_patient_id, .keep_all = T), 
       aes(x = PASCnoPASC, y = BMI, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "BMI",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bmi_boxplot_analysiscohort1.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:age boxplot - analysis group 1 ----

ggplot(metadata %>%
         filter(analysis_group_1 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("PASC", "PASC_fu") ~ "Long_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
         distinct(unique_patient_id, .keep_all = T), 
       aes(x = PASCnoPASC, y = Age, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "Age",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/age_boxplot_analysiscohort1.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:sex bar - analysis group 1 ----
counts <- metadata %>%
  filter(analysis_group_1 == 1) %>%
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("PASC", "PASC_fu") ~ "Long_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
  distinct(unique_patient_id, .keep_all = T) %>%
  group_by(PASCnoPASC) %>%
  count(Sex, name = "count") %>%
  mutate(
    total = sum(count), 
    perc = count / total,
    ypos = total - count / 2
  )

ggplot(counts, aes(x = PASCnoPASC, y = count, color = Sex, fill = PASCnoPASC)) +
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7,
           size = 0.2) + 
  #coord_polar(theta = "y") + 
  #scale_fill_manual(values = col) +
  geom_text(aes(y = count, label = paste0(count, "\n", round(perc * 100, 1))), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c(tab[4], tab[1])) +
  labs(y = "Samples",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "inside",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave(paste0("reports/figures/bar_sex_analysisgroup1.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:bmi boxplot - analysis group 2 ----

ggplot(metadata %>%
         filter(analysis_group_6 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("Acute") ~ "Acute_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID")), 
       aes(x = PASCnoPASC, y = BMI, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "BMI",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bmi_boxplot_analysiscohort2.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:age boxplot - analysis group 2 ----

ggplot(metadata %>%
         filter(analysis_group_6 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("Acute") ~ "Acute_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID")), 
       aes(x = PASCnoPASC, y = Age, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "Age",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/age_boxplot_analysiscohort2.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:sex bar - analysis group 2 ----
counts <- metadata %>%
  filter(analysis_group_6 == 1) %>%
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Acute") ~ "Acute_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID")) %>%
  group_by(PASCnoPASC) %>%
  count(Sex, name = "count") %>%
  mutate(
    total = sum(count), 
    perc = count / total,
    ypos = total - count / 2
  )

ggplot(counts, aes(x = PASCnoPASC, y = count, color = Sex, fill = PASCnoPASC)) +
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7,
           size = 0.2) + 
  #coord_polar(theta = "y") + 
  #scale_fill_manual(values = col) +
  geom_text(aes(y = count, label = paste0(count, "\n", round(perc * 100, 1))), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c(tab[4], tab[1])) +
  labs(y = "Samples",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,120)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "inside",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave(paste0("reports/figures/bar_sex_analysisgroup2.pdf"), 
       width = 3.5, height = 3.5, units = "cm")




asdfasdf <- metadata %>%
  filter(analysis_group_6 == 1) %>%
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Acute") ~ "Acute_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID"))
table(asdfasdf$Cohort)


asdfasdfasdf <- metadata %>%
  filter(Cohort == "Acute")
table(asdfasdfasdf$Sex)









## 2020 vs 2024 COVID Biomolecules ----
# 2020
samples_2020 <- 128
p_2020 <- 517
l_2020 <- 646
t_2020 <- 13263

df_2020 <- data.frame(p = p_2020,
                      l = l_2020,
                      t = t_2020) %>%
  pivot_longer(cols = everything(),
               names_to = "ome",
               values_to = "count") %>%
  mutate(year = "2020")

# 2024
samples_2024 <- 380
p_2024 <- 6972
l_2024 <- 1729
t_2024 <- 16167

df_2024 <- data.frame(p = p_2024,
                      l = l_2024,
                      t = t_2024) %>%
  pivot_longer(cols = everything(),
               names_to = "ome",
               values_to = "count") %>%
  mutate(year = "2024")

df_both <- df_2020 %>%
  bind_rows(df_2024) %>%
  group_by(year) %>%
  mutate(total = sum(count),
         proportion = count / total,
         ypos = cumsum(proportion) - 0.5 * proportion) %>%
  arrange(desc(count))


## plot:pie - biomolecule 2020 vs 2024 ----

ggplot(df_both, aes(x = year, y = proportion, fill = ome)) +
  geom_bar(stat = "identity", 
           width = 1,
           color = "white") + 
  coord_polar("y") + 
  scale_fill_manual(values = c(col[2], col[1], col[3])) +
  guides(fill = "none") +
  geom_text(aes(y = proportion, label = count), 
            color = "black", 
            size = 2) +
  labs(y = NULL,
       x = NULL) +
  theme_void() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        plot.title = element_blank(), 
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")) + 
  facet_wrap( ~ year)
ggsave(paste0("reports/figures/pie_biomoleculeIDS_filtered_2020vs2024.pdf"), 
       width = 6, height = 3, units = "cm")


## plot:pie - sample 2020 vs 2024 ----

df_samplenum <- data.frame(
  year = c("2020", "2024"),
  value = c(samples_2020, samples_2024)
)

ggplot(df_samplenum, aes(x = year, y = 1, size = value)) +
  geom_point(shape = 21, fill = "skyblue") +
  scale_size_area(max_size = 30) +  # size area = true proportional scaling
  theme_void() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
ggsave(paste0("reports/figures/pie_samples_filtered_2020vs2024.pdf"), 
       width = 6, height = 3, units = "cm")





## plot:bmi boxplot - analysis group 7 ----

ggplot(metadata %>%
         filter(analysis_group_7 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("PASC") ~ "Long_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
         distinct(unique_patient_id, .keep_all = T), 
       aes(x = PASCnoPASC, y = BMI, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "BMI",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bmi_boxplot_analysiscohort7.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:age boxplot - analysis group 7 ----

ggplot(metadata %>%
         filter(analysis_group_7 == 1) %>%
         mutate(PASCnoPASC = case_when(
           Cohort %in% c("PASC") ~ "Long_COVID",
           Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
         distinct(unique_patient_id, .keep_all = T), 
       aes(x = PASCnoPASC, y = Age, fill = PASCnoPASC, color = PASCnoPASC)) +
  geom_jitter(alpha = 0.5, 
              width = 0.1, 
              size = 0.2) +
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.2) +
  stat_compare_means(
    method = "t.test",      
    label = "p.format"        
  ) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c("darkgray", "lightgray")) +
  labs(y = "Age",
       x = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/age_boxplot_analysiscohort7.pdf"), 
       width = 3.5, height = 3.5, units = "cm")



## plot:sex bar - analysis group 7 ----
counts <- metadata %>%
  filter(analysis_group_7 == 1) %>%
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("PASC") ~ "Long_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_Long_COVID")) %>%
  distinct(unique_patient_id, .keep_all = T) %>%
  group_by(PASCnoPASC) %>%
  count(Sex, name = "count") %>%
  mutate(
    total = sum(count), 
    perc = count / total,
    ypos = total - count / 2
  )

ggplot(counts, aes(x = PASCnoPASC, y = count, color = Sex, fill = PASCnoPASC)) +
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7,
           size = 0.2) + 
  #coord_polar(theta = "y") + 
  #scale_fill_manual(values = col) +
  geom_text(aes(y = count, label = paste0(count, "\n", round(perc * 100, 1))), 
            color = "black", 
            size = 2, 
            nudge_y = -10) +
  scale_fill_manual(values = c("darkgray", "lightgray")) +
  scale_color_manual(values = c(tab[4], tab[1])) +
  labs(y = "Samples",
       x = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = "inside",
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_blank(), 
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave(paste0("reports/figures/bar_sex_analysisgroup7.pdf"), 
       width = 3.5, height = 3.5, units = "cm")
