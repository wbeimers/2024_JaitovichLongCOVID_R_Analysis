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
library(pheatmap)
library(ROTS)
library(data.table)


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
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')


proteomics <- dbGetQuery(con, 'SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_name, Sample, sample_id, ome_id, keep , rawfile_id, run_type
                           FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples, unique_patient_id, Collection_date, PG_change_collection_cutoff
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles <- rawfiles %>%
  select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl('Sample', run_type))

df <- proteomics %>%
  left_join(rawfiles, by = 'rawfile_id') %>%
  left_join(metadata, by = 'sample_id')

biomolecules_to_keep <- biomolecules %>%
  filter(omics_id == 1) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df <- df %>%
  filter(biomolecule_id %in% biomolecules_to_keep)



#### DiffExp ####
# Get grouping vector and data matrix
# Groups: Acute, Acute_fu, Acute_NC, Healthy, PASC, PASC_fu
# PASC_Cohort: first, second
# PG_change_collection_cutoff: 0, 1
group1 <- 0
group2 <- 1

diffexp_df <- filtered_df %>%
  filter(PG_change_collection_cutoff %in% c(group1, group2)) %>%
  select(standardized_name, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = "sample_id", values_from = "normalized_abundance") %>%
  tibble::column_to_rownames(var = "standardized_name")

grouping_vector <- metadata %>%
  filter(sample_id %in% colnames(diffexp_df)) %>%
  arrange(match(sample_id, colnames(diffexp_df))) %>%
  pull(PG_change_collection_cutoff)


##* ROTS DEA ----
# run ROTS
results <- ROTS(data = diffexp_df, 
                groups = grouping_vector, 
                B = 5000, 
                K = 500, 
                seed = 1234)
  
summary(results, fdr = 0.05)
  
# output df
volc_plot <- data.frame(results$logfc) %>%
  tibble::rownames_to_column(var = "standardized_name") %>%
  mutate(logfc = results.logfc) %>%
  dplyr::select(-results.logfc) %>%
  mutate(pvalue = results$pvalue) %>%
  mutate(qvalue = results$FDR)
  
volc_plot <- volc_plot %>%
  mutate(neglogpvalue = -log10(pvalue)) %>%
  mutate(diffexp = case_when(
    logfc >= 0.263 & qvalue <= 0.05 ~ "UP",
    logfc <= -0.263 & qvalue <= 0.05 ~ "DOWN",
    T ~ "NO"
  ))

fwrite(volc_plot, paste0("data/processed/AllPlates_Samples_Volcano_", group1, "_", group2, ".csv"))

volc_plot <- fread(paste0("data/processed/AllPlates_Samples_Volcano_", group1, "_", group2, ".csv"))

counts <- volc_plot %>%
  filter(diffexp %in% c("DOWN", "UP", "NO")) %>%
  count(diffexp)

ggplot(volc_plot, aes(logfc, neglogpvalue, color = diffexp, size = diffexp)) + 
  geom_point() +
  geom_vline(xintercept=c(-0.263, 0.263), 
             col="black",
             size = 0.2) +
  #geom_hline(yintercept=-log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_color_manual(values = c("#9e1b45", col[8], "#9e1b45")) +
  scale_size_manual(values = c(0.5,0.1,0.5)) +
  scale_x_continuous(limits = c(-max(abs(volc_plot$logfc)), max(abs(volc_plot$logfc))), breaks = seq(-5, 5, by = 2)) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab(paste0("Log2 Fold Change (", group1, "/", group2, ")")) +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom") +
  #facet_wrap(. ~ factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")),
  #           scales = "free",
  #           ncol = 2) +
  #geom_text(data = volc_plot %>% distinct(method),
  #          aes(x = -3.5, y = 8, label = method),
  #          size = 3, 
  #          color = "black", 
  #          hjust = 0) +
  geom_text(data = counts[counts$diffexp == "DOWN",], 
            aes(x = -3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "UP",], 
            aes(x = 3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "NO",], 
            aes(x = 0, y = Inf, 
                label = paste(n)),
            hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) 
ggsave(paste0("reports/figures/AllPlates_Samples_50pconditionMissingImputed_Volcano_small_", group1, "_", group2, ".pdf"), 
       width = 8, height = 6, units = "cm")


