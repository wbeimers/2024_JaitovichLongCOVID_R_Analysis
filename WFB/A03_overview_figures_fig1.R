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

ggplot(counts, aes(x = Cohort, y = count, fill = Cohort)) +
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
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 7),
              plot.title = element_blank(), 
              axis.title = element_text(size = 7),
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 7),
              axis.line = element_line(linewidth = 0.2),
              axis.ticks = element_line(linewidth = 0.2))
ggsave(paste0("reports/figures/bar_allsamples.pdf"), 
       width = 8, height = 4, units = "cm")


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



## plot:age boxplot - all samples ----
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

# THINK ABOUT FILTERING OUT POST-TUBE CHANGE SAMPLES


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



## Plot:bar - biomolecule overview
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
       width = 8, height = 4, units = "cm")



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
