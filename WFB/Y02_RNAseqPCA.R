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


rnaseq <- dbGetQuery(con, "SELECT standardized_name, SYMBOL, Sample, Counts, biomolecule_id
                         FROM rnaseq_measurements")
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
metadata <- dbGetQuery(con, 'SELECT Sample, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples, unique_patient_id, Collection_date
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample and QC proteomics runs

df <- rnaseq %>%
  inner_join(metadata, by = "Sample")



#### PCA ####
pca_set <- df %>%
  dplyr::select(standardized_name, Counts, Sample) %>%
  group_by(standardized_name, Sample) %>%
  summarise(Counts = mean(Counts)) %>%
  pivot_wider(names_from = Sample, values_from = Counts)

t_pca_set <- as.data.frame(t(pca_set))

# turn row into colnames
colnames(t_pca_set) <- as.character(t_pca_set[1, ])
# Remove the row
t_pca_set <- t_pca_set[-1, ]
#change to numeric
t_pca_set <- data.frame(lapply(t_pca_set, as.numeric), row.names = rownames(t_pca_set))
# re-add sample_id
t_pca_set <- tibble::rownames_to_column(t_pca_set, "Sample")
ncol(t_pca_set)


pca_score <- prcomp(t_pca_set[,c(2:19772)],
                    scale. = T)
summary(pca_score)

# make variables to plot
# scree
explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
variance <- data.frame(proportion = explained_variance,
                       PC = 1:269)

# pca scores
scores <- as.data.frame(pca_score$x)
scores <- scores %>%
  mutate(Sample = t_pca_set$Sample) %>%
  left_join(df %>% 
            ungroup() %>%
            dplyr::select(Sample, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, Collection_date) %>%
            distinct(), 
          by = "Sample") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))


ggplot(scores, aes(PC1, PC2, fill = as.Date(Collection_date, format = "%Y%m%d"))) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  #stat_ellipse(aes(color = Cohort), 
  #             geom = "path", 
  #             show.legend = FALSE,
  #             linewidth = 0.2) +
  #scale_fill_manual(values = c(pal[2], pal[4], pal[5], pal[6])) +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  xlab(paste("PC1", round(variance$proportion[1]*100, 2))) +
  ylab(paste("PC2", round(variance$proportion[2]*100, 2))) +
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
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )
ggsave("reports/figures/RNAseq_Samples_CollectionDate_PCA_PC1PC2.pdf", 
       width = 12, height = 6, units = "cm")


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