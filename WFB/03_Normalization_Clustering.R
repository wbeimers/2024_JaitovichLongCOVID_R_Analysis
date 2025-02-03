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
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


proteomics <- dbGetQuery(con, "SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement")
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules")
rawfiles <- dbGetQuery(con, "SELECT rawfile_name, Sample, sample_id, ome_id, keep , rawfile_id, run_type
                           FROM rawfiles_all")
metadata <- dbGetQuery(con, "SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`
                           FROM patient_metadata")

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample and QC proteomics runs
rawfiles <- rawfiles %>%
  select(-keep) %>%
  filter(ome_id == 1) %>%
  filter(grepl("Sample|QC", run_type))

metadata <- metadata %>%
  select(-Sample) %>%
  mutate(sample_id = as.integer(sample_id))

df <- proteomics %>%
  left_join(rawfiles, by = "rawfile_id") %>%
  left_join(metadata, by = "sample_id")


## LIMIT TO NPA FOR NOW, AFTER RE-SEARCHING COMBINE NPs BY COMPLETENESS
# Also filter for completeness

df <- df %>%
  filter(grepl("NPA", rawfile_name)) #Select only sample runs

# Show how many non-NA values there are for each protein group in each study group
na_summary <- df %>%
  group_by(Cohort, standardized_name) %>%
  summarise(na_ratio = mean(!is.na(raw_abundance)), .groups = "drop")

# Make a list of IDs to keep where there are at least 50% non-NA values in one of the cohorts
ids_to_keep <- na_summary %>%
  group_by(standardized_name) %>%
  summarise(max_na_ratio = max(na_ratio)) %>%
  filter(max_na_ratio >= 0.5) %>% 
  pull(standardized_name)


filtered_df <- df %>%
  filter(standardized_name %in% ids_to_keep)
  
  


#### PCA ####
pca_set <- filtered_df %>%
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

pca_score <- prcomp(t_pca_set[,c(2:6490)],
                    scale. = T)
summary(pca_score)

# make variables to plot
# scree
explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
variance <- data.frame(proportion = explained_variance,
                       PC = 1:400)

# pca scores
scores <- as.data.frame(pca_score$x)
scores <- scores %>%
  mutate(sample_id = t_pca_set$sample_id) %>%
  left_join(filtered_df %>% 
              select(sample_id, rawfile_id, Sample, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`) %>%
              distinct(), 
            by = "sample_id")


ggplot(scores, aes(PC2, PC3, fill = Cohort)) + 
  geom_point(shape = 21,
           size = 1,
           color = "black",
           stroke = 0.1) +
  stat_ellipse(aes(color = Cohort), 
               geom = "path", 
               show.legend = FALSE,
               linewidth = 0.2) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlab(paste("PC2", round(variance$proportion[2]*100, 2))) +
  ylab(paste("PC3", round(variance$proportion[3]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = c(0.95, 0.05), 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave("reports/figures/AllPlates_Samples_cohort_PCA_PC2PC3.pdf", 
       width = 8, height = 6, units = "cm")


## loadings ----
loadings <- pca_score$rotation %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  tibble::rownames_to_column(var = 'standardized_name')


ggplot(loadings, aes(PC1, PC2)) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  geom_text_repel(aes(label = standardized_name), size = 2) +
  xlab(paste("PC1 Loadings", round(variance$proportion[1]*100, 2))) +
  ylab(paste("PC2 Loadings", round(variance$proportion[2]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = c(0.95, 0.05), 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave("reports/figures/AllPlates_Samples_cohort_PCAloadings.pdf", 
       width = 8, height = 6, units = "cm")



