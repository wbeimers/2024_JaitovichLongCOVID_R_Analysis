#### Overview ####

# 1. PCA of each ome separately, evaluate relationships with sex, age, cohort, QOL
# 2. PCA of combined omics, evaluate relationships with sex, age, cohort, QOL
# 3. PCA of cohorts (PASC first/second,  PASC paried fu samples, and PASC + healthy)
# 4. Heatmap/heirarchical clustering of combined omics, all data
# 5. Heatmap/heirarchical clustering of combined omics, cohorts



#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(viridis)
library(RSQLite)
library(ggrepel)
library(pheatmap)
library(data.table)
library(igraph)
library(ggraph)
library(STRINGdb)


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


pal <- c("#E78AC3", 
         '#AA3377', 
         "#B3B3B3", 
         '#228833', 
         '#66CCEE', 
         '#4477AA')



# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
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
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, unique_patient_id, Collection_date, PG_change_collection_cutoff
                             FROM patient_metadata')

dbDisconnect(con)


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
  select(-rawfile_id, -rawfile_name, -Sample, -run_type, -ome_id) %>%
  bind_rows(df_t %>%
              rename(raw_abundance = Counts,
                     normalized_abundance = normalized_counts) %>%
              mutate(ome = "t"))

filtered_df_a <- filtered_df_p %>%
  mutate(ome = "p") %>%
  bind_rows(filtered_df_l %>%
              mutate(ome = "l")) %>%
  select(-rawfile_id, -rawfile_name, -Sample, -run_type, -ome_id) %>%
  bind_rows(filtered_df_t %>%
              rename(raw_abundance = Counts) %>%
              mutate(ome = "t"))
  



## PCA Setup ----
df_list <- list(p = filtered_df_p,
                l = filtered_df_l,
                t = filtered_df_t,
                a = filtered_df_a)

# this for loop we can modify the subset of data chosen for the PCA plots

# options: 
# exclude post-cutoff samples ("PreTube", "AllTube")
exclude <- "PreTube"
# PASC, PASC+fu, Paired, All
opt <- "All"
# exclude Acute and Acute_NC ("YesAcute", "NoAcute")
acute <- "NoAcute"


for (i in names(df_list)) {
  
  # Filter based on the exclude variable
  if (exclude == "PreTube") {
    x <- df_list[[i]] %>%
      filter(PG_change_collection_cutoff == 0)
  } else {
    x <- df_list[[i]]
  }
  
  # Filter based on the opt variable
  if (opt == "PASC") {
    x <- x %>%
      filter(Cohort == "PASC")
  } else if(opt == "PASC+fu"){
    x <- x %>%
      filter(Cohort %in% c("PASC", "PASC_fu"))
  } else if(opt == "Paired"){
    x <- x %>%
      group_by(unique_patient_id, standardized_name) %>%
      filter(n() > 1) %>%
      ungroup()
  } else {
    x <- x
  }
  
  # Filter based on the acute variable
  if (acute == "NoAcute") {
    x <- x %>%
      filter(!Cohort %in% c("Acute", "Acute_NC"))
  } else {
    x <- x
  } 
  
  # take normalized abundance values for PCA
  pca_set <- x %>%
    select(standardized_name, normalized_abundance, sample_id) %>%
    pivot_wider(names_from = sample_id, values_from = normalized_abundance)
  
  if (i == "a") {
    pca_set <- pca_set %>%
      select(where(~ !any(is.na(.))))
  }
  
  t_pca_set <- as.data.frame(t(pca_set))
  
  # turn row into colnames
  colnames(t_pca_set) <- as.character(t_pca_set[1, ])
  # Remove the row
  t_pca_set <- t_pca_set[-1, ]
  #change to numeric
  t_pca_set <- data.frame(lapply(t_pca_set, as.numeric), row.names = rownames(t_pca_set))
  # re-add sample_id
  t_pca_set <- tibble::rownames_to_column(t_pca_set, "sample_id") %>%
    mutate(sample_id = as.integer(sample_id))
  
  pca_score <- prcomp(t_pca_set[,c(2:ncol(t_pca_set))],
                      scale. = T)
  summary(pca_score)
  
  # make variables to plot
  # scree
  explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
  variance <- data.frame(proportion = explained_variance,
                         PC = 1:(ncol(pca_set)-1))
  
  # pca scores
  scores <- as.data.frame(pca_score$x) %>%
    mutate(sample_id = t_pca_set$sample_id) %>%
    left_join(df_list[[i]] %>% 
                select(sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, unique_patient_id, Collection_date, PG_change_collection_cutoff) %>%
                distinct(), 
              by = "sample_id")
  
  # loadings
  loadings <- pca_score$rotation %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'standardized_name')
  
  # make separate dataframes
  variance_name <- paste0("pca_variance_", i)
  pca_score_name <- paste0("pca_scores_", i)
  loadings_name <- paste0("loadings_", i)
  assign(variance_name, variance)
  assign(pca_score_name, scores)
  assign(loadings_name, loadings)
  
  }


## PCA ----
# options
# choose 2 components to plot
C1 <- "PC1"
C2 <- "PC2"

# choose a grouping variable
# discrete: Cohort, Sex, PASC_Cohort, unique_patient_id, PG_change_collection_cutoff
# continuous: Age, BMI, SF.36.QOL.Score, Collection_date
option <- "Cohort"

pca_list <- list(p = pca_scores_p,
                 l = pca_scores_l,
                 t = pca_scores_t,
                 a = pca_scores_a)

loadings_list <- list(p = loadings_p,
                      l = loadings_l,
                      t = loadings_t,
                      a = loadings_a)

for (i in names(pca_list)) {
  
  # Filter based on the option variable
  if (option == "PASC_Cohort") {
    x <- pca_list[[i]] %>%
      filter(!is.na(PASC_Cohort))
  } else if (option == "unique_patient_id") {
    x <- pca_list[[i]] %>%
      group_by(unique_patient_id) %>%
      filter(n() > 1) %>%
      ungroup()
  } else {
    x <- pca_list[[i]]
  }
  
  x <- x %>%
    mutate(PG_change_collection_cutoff = as.factor(PG_change_collection_cutoff))
  
  # Choose colors for PCA
  if (option == "unique_patient_id") {
    grouping <- "Cohort"
  } else {
    grouping <- option
  }
  
  
  p <- ggplot(x, aes(.data[[C1]], .data[[C2]], fill = .data[[grouping]])) + 
    geom_point(shape = 21,
               size = 1.5,
               alpha = 0.9,
               color = "black",
               stroke = 0.1) +
    #stat_ellipse(aes(color = Cohort), 
    #             geom = "path", 
    #             show.legend = FALSE,
    #             linewidth = 0.2) +
    #geom_text_repel(aes(label = Sample), size = 2) +
    xlab(paste(C1, round(get(paste0("pca_variance_", i))$proportion[1]*100, 2))) +
    ylab(paste(C2, round(get(paste0("pca_variance_", i))$proportion[2]*100, 2))) +
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
  
  # change colors
  pal <- c("Acute" = "#E78AC3", 
           "Acute_fu" = '#AA3377', 
           "Acute_NC" = "#B3B3B3", 
           "Healthy" = '#229100', 
           "PASC" = '#66CCEE', 
           "PASC_fu" = '#4477AA')
  
  # Add appropriate scale
  if(grouping %in% c("Age", "BMI", "SF.36.QOL.Score", "Collection_date")) {
    p <- p + scale_fill_viridis_c(option = "mako", direction = -1)
  } else if(grouping == "PASC_Cohort") {
    p <- p + scale_fill_manual(values = c(col[1], col[2])) 
  } else if(grouping == "Sex") {
    p <- p + scale_fill_manual(values = c(col[1], col[2])) 
  } else if(grouping == "PG_change_collection_cutoff") {
    p <- p + scale_fill_manual(values = c(col[1], col[2])) 
  }  else {
    p <- p + scale_fill_manual(values = pal) 
  }

  p
  ggsave(paste0("reports/figures/PCA_plots/PCA_", i, "_", option, "_", C1, C2, "_", exclude, "_", opt, "_", acute, ".pdf"), 
         width = 12, height = 6, units = "cm")
  
  
  
  # Filter based on the option variable
  y <- loadings_list[[i]] 

  p1 <- ggplot(y, aes(.data[[C1]], .data[[C2]])) + 
    geom_point(shape = 21,
               size = 1,
               color = "black",
               stroke = 0.1) +
    geom_text_repel(aes(label = standardized_name), size = 2) +
    xlab(paste(C1, "loadings", round(get(paste0("pca_variance_", i))$proportion[1]*100, 2))) +
    ylab(paste(C2, "loadings", round(get(paste0("pca_variance_", i))$proportion[2]*100, 2))) +
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
  
  p1
  ggsave(paste0("reports/figures/PCA_plots/loadings_", i, "_", option, "_", C1, C2, "_", exclude, "_", opt, "_", acute, ".pdf"), 
         width = 8, height = 6, units = "cm")
  
}



## Heatmap ----
# Make expression matrices
# options: 
# exclude post-cutoff samples ("PreTube", "AllTube")
exclude <- "PreTube"
# PASC, PASC+fu, Paired, All
opt <- "All"
# exclude Acute and Acute_NC ("YesAcute", "NoAcute")
acute <- "NoAcute"
# with NA ("YesNA", "NoNA")
na <- "NoNA"
# clutser number for rows
clusters <- 8


df_list <- list(p = filtered_df_p,
                l = filtered_df_l,
                t = filtered_df_t)

df_list1 <- list(p = df_p,
                l = df_l,
                t = df_t)

for (i in names(df_list)) {
  
  pal <- c("#E78AC3", 
           '#AA3377', 
           "#B3B3B3", 
           '#228833', 
           '#66CCEE', 
           '#4477AA')
  
  # Filter with or without NA
  if (na == "YesNA") {
    x <- df_list1[[i]]
  } else {
    x <- df_list[[i]]
  } 
  
  # Filter based on the exclude variable
  if (exclude == "PreTube") {
    x <- x %>%
      filter(PG_change_collection_cutoff == 0)
  } else {
    x <- x
  }
  
  # Filter based on the opt variable
  if (opt == "PASC") {
    x <- x %>%
      filter(Cohort == "PASC")
  } else if(opt == "PASC+fu"){
    x <- x %>%
      filter(Cohort %in% c("PASC", "PASC_fu"))
  } else if(opt == "Paired"){
    x <- x %>%
      group_by(unique_patient_id, standardized_name) %>%
      filter(n() > 1) %>%
      ungroup()
  } else {
    x <- x
  }
  
  # Filter based on the acute variable
  if (acute == "NoAcute") {
    x <- x %>%
      filter(!Cohort %in% c("Acute", "Acute_NC"))
  } else {
    x <- x
  } 
  
  # Heatmap with or without NA
  if (na == "YesNA") {
    expression_matrix <- x %>%
      dplyr::select(standardized_name, raw_abundance, sample_id) %>%
      mutate(raw_abundance = log2(raw_abundance)) %>%
      pivot_wider(names_from = sample_id, values_from = raw_abundance) %>%
      tibble::column_to_rownames(var = "standardized_name")
    
    sample_annot <- x %>%
      ungroup() %>%
      select(sample_id, Cohort, Age, Sex, BMI, SF.36.QOL.Score, PASC_Cohort, PG_change_collection_cutoff) %>%
      distinct() %>%
      tibble::column_to_rownames(var = "sample_id") %>%
      mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))
    
    clust <- expression_matrix %>%
      mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
    
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
             filename = paste0("reports/figures/Heatmap_plots/Heatmap_", i, "_", opt, "_", exclude, "_", na, "_", acute, ".png"),
             width = 8,
             height = 8)
    
    hclust1 <- hclust(dist(clust), method = "complete")
    hclust_clusters <- cutree(as.hclust(hclust1), k = clusters)
    fwrite(hclust_clusters, paste0("reports/figures/Heatmap_plots/Clusters_", i, "_", opt, "_", exclude, "_", na, "_", acute, "_", clusters, ".txt"))
  } else {
    expression_matrix <- x %>%
      select(standardized_name, sample_id, normalized_abundance) %>%
      pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
      tibble::column_to_rownames(var = "standardized_name")
    
    sample_annot <- x %>%
      ungroup() %>%
      select(sample_id, Cohort, Age, Sex, BMI, SF.36.QOL.Score, PASC_Cohort, PG_change_collection_cutoff) %>%
      distinct() %>%
      tibble::column_to_rownames(var = "sample_id") %>%
      mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))
    
    pheatmap(expression_matrix,
             #color = viridis(24, direction = 1, option = "plasma"),
             color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
             breaks = c(-11, -9, -7, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11),
             cluster_rows = T,
             cutree_rows = clusters, 
             gaps_row = T,
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
             filename = paste0("reports/figures/Heatmap_plots/Heatmap_", i, "_", opt, "_", exclude, "_", na, "_", acute, ".png"),
             width = 8,
             height = 8)
    
    hclust1 <- hclust(dist(expression_matrix), method = "complete")
    hclust_clusters <- cutree(as.hclust(hclust1), k = clusters)
    fwrite(as.data.frame(hclust_clusters) %>%
             tibble::rownames_to_column(var = "ID"), paste0("reports/figures/Heatmap_plots/Clusters_", i, "_", opt, "_", exclude, "_", na, "_", acute, "_", clusters, ".csv"))
  } 

}



## Heatmap:AnalysisGroup1 Overlap ----
# get information from A02_effectsizevolcano_Action2.R for volc_plot

bmol <- volc_plot %>%
  filter(q_value <= 0.01) %>%
  pull(biomolecule_id)

expression_matrix <- filtered_df_a %>%
  filter(biomolecule_id %in% bmol) %>%
  filter(analysis_group_1 == 1) %>%
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
         SF.36.QOL.Score, 
         #PASC_Cohort, 
         #PG_change_collection_cutoff
         ) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))

dat <- filtered_df_a %>%
  filter(biomolecule_id %in% bmol) %>%
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

pheat <- pheatmap(expression_matrix,
         #color = viridis(24, direction = 1, option = "plasma"),
         color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(15)),
         breaks = c(-4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4),
         cluster_rows = T,
         #cutree_rows = k, 
         #gaps_row = T,
         cluster_cols = T,
         treeheight_row = 0,
         treeheight_col = 10,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         annotation_row = row_annot,
         annotation_col = sample_annot,
         annotation_colors = list(Cohort = c(Acute_fu = pal[2], Healthy = pal[4], PASC = pal[5], PASC_fu = pal[6]),
                                  Ome = c(protein = col[1], lipid = col[2], transcript = col[3]),
                                  SF.36.QOL.Score = rev(colorRampPalette(c("white", "black"))(4))),
         filename = paste0("reports/figures/Heatmap_plots/Heatmap_AllOmes_AnalysisGroup1_LRMdiffexpbmol_0,01cutoff.png"),
         width = 8,
         height = 8)

cluster_assignments <- cutree(pheat$tree_row, k = k)
fwrite(as.data.frame(cluster_assignments) %>%
         tibble::rownames_to_column(var = "biomolecule_id"), paste0("reports/figures/Heatmap_plots/Clusters_AllOmes_AnalysisGroup1_4methodoverlap_2clusters.csv"))



## STRING Network ----
# choose cluster
# 1. load the list of clusters
i <- "p" # "p", "l", "t"

cluster <- fread(paste0("reports/figures/Heatmap_plots/Clusters_", i, "_", opt, "_", exclude, "_", na, "_", acute, "_", clusters, ".csv"))

# 2. choose which cluster to look at
cluster1 <- cluster %>%
  filter(hclust_clusters == 1)

# 3. start string up
string_db <- STRINGdb$new(version = "12.0", 
                          species = 9606, 
                          score_threshold = 200, 
                          network_type = "full", 
                          link_data = "combined_only", 
                          input_directory = "")

# 4. map ids
cluster1_mapped <- string_db$map(cluster1, "ID", removeUnmappedRows = TRUE)

hits <- cluster1_mapped$STRING_id[1:200]

string_db$plot_network(hits)



## Protein Co-Expression Network ----
# use df_p for including NAs, filtered_df_p to exclude












expression_matrix <- filtered_df_p %>%
  filter(PG_change_collection_cutoff == 0) %>% # exclude post- tube change
  #filter(Cohort == "PASC") %>% # filter for specific cohort
  filter(!Cohort %in% c("Acute", "Acute_NC")) %>% # get rid of acute samples 
  dplyr::select(standardized_name, normalized_abundance, sample_id) %>%
  #mutate(raw_abundance = log2(raw_abundance)) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  tibble::column_to_rownames(var = "standardized_name")

# 1. Make correlation matrix
cor_matrix <- cor(t(expression_matrix), method = "pearson")

# 2. Make adjacency matrix
adjacency_matrix <- cor_matrix
adjacency_matrix[abs(cor_matrix) < 0.5] <- 0  # Set threshold to 0.5

# 3. Plot
graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = T)
E(graph)$weight <- abs(E(graph)$weight) # take absolute value of correlations
degree_centrality <- degree(graph)
communities <- cluster_louvain(graph)

ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
  geom_node_point(aes(size = degree_centrality), color = "blue") +
  theme_minimal()





