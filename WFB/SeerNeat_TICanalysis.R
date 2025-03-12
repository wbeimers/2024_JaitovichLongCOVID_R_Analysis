#### analyze the neat re-prep of some of the samples 
library(pheatmap)
library(data.table)
library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(rrcovNA)


pal <- c('#EE6677', # Acute
         '#AA3377', # Acute_fu
         '#CCBB44', # Acute_NC
         '#228833', # Healthy
         '#66CCEE', # PASC
         '#4477AA') # PASC_fu




tic_seer <- fread("data/metadata/summary_tic_results.csv") %>%
  filter(grepl("Plate", RawFilename)) %>%
  filter(grepl("NPA|NPB", RawFilename)) %>%
  filter(!grepl("20240712", RawFilename)) %>%
  mutate(NP = word(RawFilename, 6, 6, "_")) %>%
  mutate(Plate = word(RawFilename, 3, 4, "_")) %>%
  mutate(Sample = word(RawFilename, 5, 5, "_")) %>%
  mutate(Sample = sub("-", ".", Sample)) %>%
  mutate(Sample = case_when(
    Sample == '1032.2' ~ '1132.2',
    Sample == '1034.2' ~ '1134.2',
    Sample == '1039.2' ~ '1139.2',
    Sample == '1040.2' ~ '1140.2',
    Sample == '1050.2' ~ '1150.2',
    T ~ Sample
  )) %>%
  left_join(metadata, by = "Sample")

tic_neat <- fread("data/metadata/summary_tic_results.csv") %>%
  filter(grepl("NeatRerun", RawFilename)) %>%
  filter(!grepl("blank", RawFilename)) %>%
  mutate(Sample = word(RawFilename, 5, 5, "_")) %>%
  mutate(Sample = sub("-", ".", Sample)) %>%
  left_join(metadata, by = "Sample")


## Make different boxplots by different variables ----

ggplot(tic, aes(Plate, Ms1TicSum, color = NP, fill = NP)) + 
  geom_boxplot(width = 0.25, 
               alpha = 0.5, 
               outliers = F,
               size = 0.2) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  labs(x = "Plate",
       y = "MS1 TIC") +
  scale_y_continuous(expand = c(0,0), limits = c((min(tic$Ms1TicSum) / 1.1), 
                                                 (max(tic$Ms1TicSum) * 1.1))) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7), 
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
ggsave(paste0('reports/figures/Proteomics_TICBoxplots_Plate.pdf'), 
       width = 16, height = 6, units = "cm")



## Combine neat and Seer ----
tic_seer_comb <- tic_seer %>%
  dplyr::select(Ms1TicSum, NP, Sample, Cohort)

tic_comb <- tic_neat %>%
  dplyr::select(Ms1TicSum, Sample, Cohort) %>%
  mutate(NP = "neat") %>%
  bind_rows(tic_seer_comb) %>%
  filter(Sample %in% tic_neat$Sample) %>%
  mutate(NP = factor(NP, levels = c("NPA", "NPB", "neat")))


ggplot(tic_comb, aes(Sample, Ms1TicSum, fill = NP)) + 
  #geom_point(size = 1) +
  geom_col(position = position_dodge(),
           width = 0.7) +
  scale_fill_manual(values = col) +
  labs(x = "Sample",
       y = "MS1 TIC") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 
                                                 max(tic_comb$Ms1TicSum) * 1.1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 90), 
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
ggsave(paste0('reports/figures/Proteomics_TICBoxplots_NeatvsSeer.pdf'), 
       width = 24, height = 10, units = "cm")



## Correlation with protein group IDs? ----
# join PG ID count from run_ids dataframe
tic_neat <- tic_neat %>%
  left_join(run_ids_cor %>% dplyr::select(sample_id, count), by = "sample_id") %>%
  mutate(cor = cor(Ms1TicSum, count))

ggplot(tic_neat, aes(Ms1TicSum, count, color = Cohort)) + 
  geom_point() +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('MS1 TIC') +
  ylab('Protein Groups') +
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
ggsave('reports/figures/NeatRePrep_MS1TICNeatCount_cohort_Point.pdf', 
       width = 8, height = 6, units = 'cm')


# do same thing for NPA NPB with seer
tic_seer <- tic_seer %>%
  left_join(run_ids %>% dplyr::select(sample_id, count), by = "sample_id")

ggplot(tic_seer %>% filter(NP == "NPA"), aes(Ms1TicSum, count, color = Cohort)) + 
  geom_point() +
  scale_color_manual(values = pal) +
  xlab('NPA MS1 TIC') +
  ylab('Protein Groups') +
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
ggsave('reports/figures/AllPlates_Samples_NPA_MS1TICNeatCount_cohort_Point.pdf', 
       width = 8, height = 6, units = 'cm')





