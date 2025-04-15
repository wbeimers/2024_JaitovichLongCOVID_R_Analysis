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
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples, unique_patient_id, Collection_date
                           FROM patient_metadata')

dbDisconnect(con)


## Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample and QC proteomics runs
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


## Can make a combined NP unfiltered data frame here 
df <- df %>%
  mutate(NP = case_when(
    grepl("NPA", rawfile_name) == T ~ "NPA",
    grepl("NPB", rawfile_name) == T ~ "NPB"
  )) %>%
  group_by(standardized_name, NP) %>%
  mutate(na_count = sum(is.na(raw_abundance))) %>%
  ungroup() %>%
  group_by(standardized_name) %>%
  mutate(keep_group = NP[which.min(na_count)]) %>%  
  filter(NP == keep_group) %>%  
  dplyr::select(-na_count, -keep_group, -NP)  

length(unique(df$standardized_name))




#### Heatmap:all ####
expression_matrix <- filtered_df %>%
  select(standardized_name, sample_id, normalized_abundance) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  tibble::column_to_rownames(var = "standardized_name")

#Make annotation dataframe for the sample groups
sample_annot <- filtered_df %>%
  ungroup() %>%
  select(sample_id, Cohort, Age, Sex, BMI, SF.36.QOL.Score) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))

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
         annotation_colors = list(Cohort = c(Acute = pal[1], Acute_fu = pal[2], Acute_NC = pal[3], Healthy = pal[4], PASC = pal[5], PASC_fu = pal[6])),
         filename = "reports/figures/AllPlates_Samples_heatmap_Cohort.png",
         width = 8,
         height = 8)



# Get clustering Information for proteins
protein_hclust <- hclust(dist(expression_matrix), method = "complete")

plot(protein_hclust, labels = F)
k = 16
rect.hclust(protein_hclust, k = k, border = "red")

protein_hclust_clusters <- cutree(as.hclust(protein_hclust), k = k)







#### Heatmap:PASC_Cohort ####
expression_matrix <- filtered_df %>%
  filter(!is.na(PASC_Cohort)) %>%
  select(standardized_name, sample_id, normalized_abundance) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  tibble::column_to_rownames(var = "standardized_name")

#Make annotation dataframe for the sample groups
sample_annot <- filtered_df %>%
  ungroup() %>%
  filter(!is.na(PASC_Cohort)) %>%
  select(sample_id, PASC_Cohort, Age, Sex, BMI, SF.36.QOL.Score) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))

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
         annotation_colors = list(PASC_Cohort = c(first = col[3], second = col[5])),
         filename = "reports/figures/AllPlates_Samples_heatmap_PASC_Cohort.png",
         width = 8,
         height = 8)




## Heatmap Showing NAs ----

expression_matrix <- df %>%
  dplyr::select(standardized_name, raw_abundance, sample_id) %>%
  mutate(raw_abundance = log2(raw_abundance)) %>%
  pivot_wider(names_from = sample_id, values_from = raw_abundance) %>%
  tibble::column_to_rownames(var = "standardized_name")

# Make annotation dataframe for the sample groups
sample_annot <- df %>%
  ungroup() %>%
  dplyr::select(sample_id, Cohort, Age, Sex, BMI, SF.36.QOL.Score) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))


# Make annotation dataframe for the rows of proteins
# remove NAs to cluster
clust <- expression_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) # column mean

# do protein clustering
protein_hclust <- hclust(dist(clust), method = "complete")

# check clustering and choose number of clusters
plot(protein_hclust, labels = F)
k = 8
rect.hclust(protein_hclust, k = k, border = "red")

# make annotation dataframe from clustering info
row_annot <- data.frame(cluster = cutree(as.hclust(protein_hclust), k = k),
                        row.names = rownames(clust))


pheatmap(expression_matrix,
         clustering_distance_rows = dist(clust), 
         clustering_distance_cols = dist(t(clust)),
         #color = viridis(24, direction = 1, option = "plasma"),
         color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
         na_col = "lightgray",
         breaks = c(-11, -9, -7, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11),
         cluster_rows = T,
         #cutree_rows = 8, 
         #gaps_row = T,
         cluster_cols = T,
         cutree_cols = 5,
         gaps_col = T,
         treeheight_row = 0,
         treeheight_col = 10,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         #annotation_row = row_annot,
         annotation_col = sample_annot,
         annotation_colors = list(Cohort = c(Acute = pal[1], Acute_fu = pal[2], Acute_NC = pal[3], Healthy = pal[4], PASC = pal[5], PASC_fu = pal[6])),
         filename = "reports/figures/AllPlates_Samples_heatmap_Cohort_NAs_clusterscols.png",
         width = 8,
         height = 8)

