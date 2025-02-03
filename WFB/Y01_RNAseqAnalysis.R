#### libraries/colors/files ####
# libraries 
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


# colors
# Make a classic palette
col <- brewer.pal(6, 'Set2') 

# Make another palette
pal <- c('#66C2A5',
         '#FFD92F',
         '#8DA0CB',
         '#FC8D62',
         '#A6CEE3',
         '#E78AC3',
         '#A6D854',
         '#FDB462',
         '#B3B3B3',
         '#B2DF8A')

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')

pal1 <- c('#AA3377',
          '#228833', 
          '#66CCEE', 
          '#4477AA')

# Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),'white', colorRampPalette(col1)(100))

# plot colors
pie(rep(1, length(col)), col = col , main='') 


# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

rnaseq <- dbGetQuery(con, 'SELECT standardized_name, SYMBOL, Sample, Counts, biomolecule_id
                         FROM rnaseq_measurements')

metadata <- dbGetQuery(con, 'SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`
                           FROM patient_metadata')

dbDisconnect(con)


## Merge transcript and metadata
rnaseq <- rnaseq %>%
  left_join(metadata, by = 'Sample') %>%
  filter(!is.na(sample_id))

#### sample order ----
quant_rnaseq <- rnaseq %>%
  select(Counts, sample_id, Sample, Cohort) %>%
  arrange(sample_id) %>%
  mutate(order = dense_rank(sample_id)) %>%
  mutate(Cohort = as.factor(Cohort))


ggplot(quant_rnaseq, aes(factor(order), Counts, fill = Cohort, color = Cohort)) + 
  geom_boxplot(width = 0.4, 
               alpha = 0.25, 
               outliers = F,
               size = 0.1) +
  scale_fill_manual(values = pal1) +
  scale_color_manual(values = pal1) +
  labs(x = NULL,
       y = "Counts") +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(title="Cohort")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), 
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
ggsave(paste0('reports/figures/AllPlates_Sample_rnseqCountBoxplots_Cohort.pdf'), 
       width = 16, height = 6, units = "cm")



#### heatmap ----
expression_matrix <- rnaseq %>%
  mutate(log2Counts = log2(Counts)) %>%
  select(standardized_name, sample_id, Counts) %>%
  pivot_wider(names_from = sample_id, values_from = Counts) %>%
  tibble::column_to_rownames(var = "standardized_name")

#Make annotation dataframe for the sample groups
sample_annot <- rnaseq %>%
  select(sample_id, Cohort) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id")

#make annotation dataframe for the genes
#heatmap_annot <- human_gene_list %>%
#  dplyr::select(c("Entry", "Gene.Ontology..cellular.component."))
#heatmap_annot <- semi_join(heatmap_annot, twentyfive_NPs_pgs_f_imputed_log2,
#                           by = join_by("Entry" == "allgenes"))
#rownames(heatmap_annot) <- heatmap_annot[,1]
#heatmap_annot <- mutate(heatmap_annot, Intracellular = ifelse(grepl("cytosol|cytoplasm", Gene.Ontology..cellular.component.), "Yes", "No"))
#heatmap_annot <- mutate(heatmap_annot, Membrane = ifelse(grepl("membrane", Gene.Ontology..cellular.component.), "Yes", "No"))
#heatmap_annot <- mutate(heatmap_annot, Plasma = ifelse(grepl("plasma", Gene.Ontology..cellular.component.), "Yes", "No"))
#heatmap_annot <- heatmap_annot[,-1]
#heatmap_annot <- heatmap_annot[,-1]


pheatmap(expression_matrix,
         #color = viridis(24, direction = 1, option = "plasma"),
         color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
         breaks = c(-11, -9, -7, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11),
         cluster_rows = T,
         cluster_cols = T,
         treeheight_row = 0,
         treeheight_col = 10,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         #annotation_row = heatmap_annot,
         annotation_col = sample_annot,
         annotation_colors = list(Cohort = c(Acute_fu = pal[2], Healthy = pal[4], PASC = pal[5], PASC_fu = pal[6])),
         filename = "reports/figures/AllPlates_Samples_rnaseq_heatmap.png",
         width = 8,
         height = 8)






