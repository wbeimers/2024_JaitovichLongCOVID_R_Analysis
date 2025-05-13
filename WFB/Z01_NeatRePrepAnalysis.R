#### analyze the neat re-prep of some of the samples 
library(pheatmap)
library(data.table)
library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(rrcovNA)
library(ROTS)


pal <- c('#EE6677', # Acute
         '#AA3377', # Acute_fu
         '#CCBB44', # Acute_NC
         '#228833', # Healthy
         '#66CCEE', # PASC
         '#4477AA') # PASC_fu


# read in report
file <- "20250228_WFB_JaitovichLongCOVID_NeatRePrep_dDIA_stringent_Report_WFB_Report (Normal).tsv"

names(fread(paste0("data/spectronaut_output/", file), nrows = 0))

spec <- fread(paste0("data/spectronaut_output/", file), 
              sep = "\t",
              select = c("R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Quantity"))

#### PGs ####
# Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
spec %>%
  dplyr::summarise(n = dplyr::n(), .by = c(R.Condition, PG.ProteinGroups, R.FileName)) %>%
  dplyr::filter(n > 1L)   

spec_PG <- pivot_wider(spec,
                       names_from = R.FileName,
                       values_from = PG.Quantity,
                       values_fill = NA,
                       values_fn = list(PG.Quantity = mean))

colSums(!is.na(spec_PG))

# fixing columns so they all line up
spec_PG <- spec_PG %>%
  group_by(PG.ProteinGroups) %>%
  summarise_all(~ na.omit(.) %>% .[1]) %>%
  ungroup()

plot(colSums(!is.na(spec_PG)))

# save ungrouped pg quant matrix
fwrite(spec_PG, file = "data/processed/PG_Matrix_NeatRePrep.csv")



# make one data matrix for raw and one for imputed
data_matrix <- spec_PG %>%
  dplyr::select(-R.Condition)
data_matrix1 <- spec_PG %>%
  dplyr::select(-R.Condition)

# Log2 transform the protein groups
data_matrix_tf <- log2(data_matrix1[, -1])
data_matrix_tf <- data.frame(ID = data_matrix1$PG.ProteinGroups, data_matrix_tf, stringsAsFactors = FALSE)

# impute with ImpSeq
data_matrix_tf1 <- data_matrix_tf %>%
  dplyr::select(where(is.double))

start_time <- Sys.time()
data_imp <- impSeq(data_matrix_tf1)
end_time <- Sys.time()
print(end_time - start_time)

data_matrix_tf_imp <- data_matrix_tf %>%
  dplyr::select(ID) %>%
  cbind(data_imp)



#fix column names to match rawfile_name in rawfiles_all table
colnames(data_matrix_tf_imp) <- gsub("^X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^X", "", colnames(data_matrix))
colnames(data_matrix_tf_imp) <- gsub("^N\\.\\.X", "", colnames(data_matrix_tf_imp))
colnames(data_matrix) <- gsub("^N\\.\\.X", "", colnames(data_matrix))
colnames(data_matrix) <- sub("-", ".", colnames(data_matrix))

#### Make longer and merge abundances with normalized and imputed ####
data_matrix_long <- data_matrix %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "raw_abundance")

data_matrix_tf_imp_long <- data_matrix_tf_imp %>%
  dplyr::rename(PG.ProteinGroups = ID) %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "rawfile_name", values_to = "normalized_abundance")

data_matrix_all_long <- data_matrix_long %>%
  left_join(data_matrix_tf_imp_long, by = c("PG.ProteinGroups", "rawfile_name"))


#### raw file info ####
# find info on run order for each of the QC runs and make a new dataframe combining them
rawfiles <- NULL
path <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Proteomics/20250224_NeatRePrep/runs"
files <- list.files(path = path,
                    pattern = "\\.raw.*",
                    recursive = T,
                    full.names = T)
rawfiles <- rbind(rawfiles, file.info(files))

#turn the rows into a column
rawfiles <- tibble::rownames_to_column(rawfiles, "names") %>%
  filter(!grepl("blank", names))

#make the rawfile_name column of the name of each rawfile
rawfiles$rawfile_name <- unname(sapply(rawfiles$names, function(x) strsplit(x, "/")[[1]][7]))
rawfiles$rawfile_name <- sub(pattern = ".raw", replacement = "", x = rawfiles$rawfile_name)
rawfiles$rawfile_name <- sub(pattern = "-", replacement = ".", x = rawfiles$rawfile_name)

#rename mtime to timestamp
colnames(rawfiles)[colnames(rawfiles) == "mtime"] <- "timestamp"

#make a rawfile order by runtime
rawfiles <- rawfiles[with(rawfiles, order(as.POSIXct(timestamp))),]
rawfiles$rawfile_id <- seq_along(rawfiles$timestamp)

#remove unnecessary columns
rawfiles <- subset(rawfiles, select = -c(names, size, isdir, mode, ctime, atime, exe))

#add a column for sample to go back to other samples
rawfiles$Sample <- sapply(strsplit(rawfiles$rawfile_name, "_"), function(x) paste(x[5], collapse = "_"))
rawfiles$Sample <- sub(pattern = "-", replacement = ".", x = rawfiles$Sample)


# merge back with the long format data frame
data_matrix_all_long <- data_matrix_all_long %>%
  left_join(rawfiles, by = "rawfile_name")

con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

metadata <- dbGetQuery(con, 'SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples, unique_patient_id, Collection_date
                           FROM patient_metadata')

dbDisconnect(con)


neat_data_df <- data_matrix_all_long %>%
  left_join(metadata, by = "Sample")

## get run_ids from the 02_Pep_PG_Numbers.R script and combine here to plot together ----

neat_run_ids <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  inner_join(neat_data_df %>% ungroup() %>% dplyr::select(sample_id, rawfile_id, Cohort), by = 'sample_id') %>%
  distinct() %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())

ggplot(neat_run_ids, aes(order, count, color = Cohort)) + 
  geom_point() + 
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('Run Order') +
  ylab('Filtered PG IDs') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4000)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Study Group')) +
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
ggsave('reports/figures/NeatRePrep_PGvsRunOrder_Set_Point.pdf', 
       width = 16, height = 6, units = 'cm')


# combined plot in sample id order
run_ids_com <- run_ids %>%
  filter(sample_id %in% neat_run_ids$sample_id) %>%
  mutate(method = "seer")

combined_run_ids <- neat_run_ids %>%
  mutate(method = "neat") %>%
  bind_rows(run_ids_com)

ggplot(combined_run_ids, aes(sample_id, count, color = method)) + 
  geom_point() + 
  geom_line() +
  scale_color_manual(values = c(col[2], col[3])) +
  xlab('Sample ID') +
  ylab('Filtered PG IDs') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7000)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Study Group')) +
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
ggsave('reports/figures/NeatRePrep_PGvsSampleID_Method_Point.pdf', 
       width = 16, height = 6, units = 'cm')


## Make correlation of neat protein IDs and seer protein IDs and see if they correlate ----

run_ids_cor <- run_ids_com %>%
  dplyr::select(sample_id, count) %>%
  dplyr::rename(count_seer = count) %>%
  left_join(neat_run_ids %>% dplyr::select(sample_id, count), by = "sample_id") %>%
  mutate(cor = cor(count_seer, count)) %>%
  left_join(metadata, by = "sample_id")

ggplot(run_ids_cor, aes(count, count_seer, color = Cohort)) + 
  geom_point() +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('Neat Count') +
  ylab('Seer Count') +
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
ggsave('reports/figures/NeatRePrep_neatcountseercount_cohort_Point.pdf', 
       width = 8, height = 6, units = 'cm')
  


## plot:single protein boxplots of cohorts ----
for (i in unique(neat_data_df$PG.ProteinGroups)) {
  
  poi <- i
  
  single_prot_df <- neat_data_df %>%
    filter(grepl(poi, PG.ProteinGroups))
  
  ggplot(single_prot_df, aes(Cohort, normalized_abundance, fill = Cohort, color = Cohort)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_manual(values = c(pal[4], pal[6])) +
    scale_color_manual(values = c(pal[4], pal[6])) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = NULL,
         y = "Log2 Abundance") +
    scale_y_continuous(expand = c(0,0), limits = c((min(single_prot_df$normalized_abundance) / 1.1), 
                                                   (max(single_prot_df$normalized_abundance) * 1.1))) +
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
  
  ggsave(paste0('reports/figures/NeatRePrep_SingleProtein/NeatRePrep_singleprotein_', poi, '_distribution_Cohort.pdf'), 
         width = 8, height = 6, units = "cm")
}




## PCA ----

pca_set <- neat_data_df %>%
  dplyr::select(PG.ProteinGroups, normalized_abundance, sample_id) %>%
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

pca_score <- prcomp(t_pca_set[,c(2:3228)],
                    scale. = T)
summary(pca_score)

# make variables to plot
# scree
explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
variance <- data.frame(proportion = explained_variance,
                       PC = 1:70)

# pca scores
scores <- as.data.frame(pca_score$x)
scores <- scores %>%
  mutate(sample_id = t_pca_set$sample_id) %>%
  left_join(neat_data_df %>% 
              ungroup() %>%
              dplyr::select(sample_id, rawfile_id, Sample, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, Collection_date) %>%
              distinct(), 
            by = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score)) %>%
  mutate(Collection_date = as.integer(Collection_date))


ggplot(scores, aes(PC1, PC2, fill = as.Date(Collection_date, format = "%Y%m%d"), shape = Cohort)) + 
  geom_point(size = 1,
             color = "black",
             stroke = 0.1) +
  #stat_ellipse(aes(color = Cohort), 
  #             geom = "path", 
  #             show.legend = FALSE,
  #             linewidth = 0.2) +
  #geom_text_repel(aes(label = Sample), size = 2) +
  #scale_fill_manual(values = c(pal[4], pal[6])) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "Collection_date") +
  scale_shape_manual(values = c(21, 24)) +
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
ggsave("reports/figures/NeatRePrep_collectiondate_PCA_PC1PC2.pdf", 
       width = 12, height = 6, units = "cm") 



## Heatmap Showing NAs ----

expression_matrix <- neat_data_df %>%
  dplyr::select(PG.ProteinGroups, raw_abundance, sample_id) %>%
  mutate(raw_abundance = log2(raw_abundance)) %>%
  pivot_wider(names_from = sample_id, values_from = raw_abundance) %>%
  tibble::column_to_rownames(var = "PG.ProteinGroups")

# Make annotation dataframe for the sample groups
sample_annot <- neat_data_df %>%
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
         na_col = "white",
         breaks = c(-11, -9, -7, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11),
         cluster_rows = T,
         cluster_cols = T,
         cutree_rows = 8,  
         gaps_row = TRUE,
         treeheight_row = 0,
         treeheight_col = 10,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         annotation_row = row_annot,
         annotation_col = sample_annot,
         annotation_colors = list(Cohort = c(Healthy = pal[4], PASC_fu = pal[6])),
         filename = "reports/figures/NeatRePrep_heatmap_Cohort_withNAs_clusters.png",
         width = 8,
         height = 8)



# Get clustering Information for proteins
hm <- pheatmap(expression_matrix)
row_cluster <- data.frame(cluster = cutree(hm$tree_row, k = 5))
row_cluster_3 <- filter(row_cluster, row_cluster$cluster == 3)
row_cluster_3 <- rownames(row_cluster_3)




## DiffExpVolcano 50% Filtered ----
# Show how many non-NA values there are for each protein group in each study group
na_summary <- neat_data_df %>%
  group_by(Cohort, PG.ProteinGroups) %>%
  summarise(na_ratio = mean(!is.na(raw_abundance)), .groups = 'drop')

# Make a list of IDs to keep where there are at least 50% non-NA values in one of the cohorts
ids_to_keep <- na_summary %>%
  group_by(PG.ProteinGroups) %>%
  summarise(max_na_ratio = max(na_ratio)) %>%
  filter(max_na_ratio >= 0.5) %>% 
  pull(PG.ProteinGroups)

neat_data_df_filt <- neat_data_df %>%
  filter(PG.ProteinGroups %in% ids_to_keep) %>%
  filter(!is.na(normalized_abundance))


# Get grouping vector and data matrix
# Groups: Healthy, PASC_fu
group1 <- "Healthy"
group2 <- "PASC_fu"

diffexp_df <- neat_data_df_filt %>%
  filter(Cohort %in% c(group1, group2)) %>%
  dplyr::select(PG.ProteinGroups, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = "sample_id", values_from = "normalized_abundance") %>%
  tibble::column_to_rownames(var = "PG.ProteinGroups")

grouping_vector <- metadata %>%
  filter(sample_id %in% colnames(diffexp_df)) %>%
  arrange(match(sample_id, colnames(diffexp_df))) %>%
  pull(Cohort)


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
  tibble::rownames_to_column(var = "PG.ProteinGroup") %>%
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

fwrite(volc_plot, paste0("data/processed/NeatRePrep_Volcano_", group1, "_", group2, ".csv"))

volc_plot <- fread(paste0("data/processed/NeatRePrep_Volcano_", group1, "_", group2, ".csv"))

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
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
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
ggsave(paste0("reports/figures/NeatRePrep_Volcano_50pconditionMissingImputed_small_", group1, "_", group2, ".pdf"), 
       width = 8, height = 6, units = "cm")


## Enrichment analysis from 50% presence volcano ---- 
## EnrichR
library(enrichR)
setEnrichrServer("https://maayanlab.cloud/Speedrichr/")
library(UniProt.ws)
up <- UniProt.ws(taxId = 9606)

dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
         "GO_Biological_Process_2023")

#subset data by fold change and pvalue


#split protein groups and olink groups into individual rows
volc_plot1 <- volc_plot %>%
  separate_rows(PG.ProteinGroup, sep = ";") 


gene_list <- volc_plot1
mapping <- UniProt.ws::select(up, keys = gene_list$PG.ProteinGroup, keytype = "UniProtKB", columns = "gene_primary")
gene_list <- gene_list %>%
  left_join(mapping, by = c("PG.ProteinGroup" = "From")) 
  
background <- gene_list$Gene.Names..primary.
background <- unlist(strsplit(background, ";"))
write.csv(background, "data/processed/NeatRePrep_50percent_background.csv")
  
gene_list <- gene_list %>%
  filter(diffexp == "DOWN")
  
input <- gene_list$Gene.Names..primary.
input <- unlist(strsplit(input, ";"))
write.csv(input, "data/processed/NeatRePrep_50percent_PASC_fu_Input.csv")
  
enriched <- enrichr(input, dbs, background = background)
  
enrichr <- enriched$GO_Molecular_Function_2023

enrichr <- enrichr %>%
  mutate(`-log10qvalue` = -log10(Adjusted.P.value)) %>%
  mutate(Rank = as.numeric(enrichr$Rank))

ggplot(enrichr, aes(Odds.Ratio, `-log10qvalue`)) +
  geom_point(shape = 16,
             size = 2,
             aes(color = Rank)) +
  geom_text_repel(aes(label = Term), size = 2) +
  scale_color_viridis(option = "plasma") +
  xlab("Odds Ratio") +
  ylab("-log10(q-value)") +
  #xlim(c(-150, 150)) +
  ylim(c(0, 3)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(color = "gray", fill = NA, size = 0.2))
ggsave("reports/figures/NeatRePrep_50percent_ENRICHr_GOCC_PASCfu.pdf", 
       width = 18, height = 12, units = "cm")


## Enrichment analysis from unique detections in PASC_fu (not healthy) ----
# limit to only uniquely detected in PASC_fu and not Healthy
neat_data_df_PASC <- neat_data_df %>%
  group_by(PG.ProteinGroups) %>%
  filter(
    any(!is.na(raw_abundance[Cohort == "PASC_fu"])) &
      all(is.na(raw_abundance[Cohort == "Healthy"]))
  )

neat_data_df_Healthy <- neat_data_df %>%
  group_by(PG.ProteinGroups) %>%
  filter(
    any(!is.na(raw_abundance[Cohort == "Healthy"])) &
      all(is.na(raw_abundance[Cohort == "PASC_fu"]))
  )

# count which proteins have most completeness and are only detected in PASC_fu
neat_data_df_PASC_counts <- neat_data_df_PASC %>%
  group_by(PG.ProteinGroups) %>%
  summarise(non_na_count = sum(!is.na(raw_abundance)))



gene_list <- unique(neat_data_df$PG.ProteinGroups)
mapping <- UniProt.ws::select(up, keys = gene_list, keytype = "UniProtKB", columns = "gene_primary")
background <- unlist(strsplit(mapping$Gene.Names..primary., ";"))
write.csv(background, "data/processed/NeatRePrep_All_background.csv")

gene_list <- neat_data_df_PASC_counts$PG.ProteinGroups
mapping <- UniProt.ws::select(up, keys = gene_list, keytype = "UniProtKB", columns = "gene_primary")
input <- unlist(strsplit(mapping$Gene.Names..primary., ";"))
write.csv(input, "data/processed/NeatRePrep_All_PASCfu_unique_Input.csv")

enriched <- enrichr(input, dbs, background = background)

enrichr <- enriched$GO_Cellular_Component_2023

enrichr <- enrichr %>%
  mutate(`-log10qvalue` = -log10(Adjusted.P.value)) %>%
  mutate(Rank = as.numeric(enrichr$Rank))

ggplot(enrichr, aes(Odds.Ratio, `-log10qvalue`)) +
  geom_point(shape = 16,
             size = 2,
             aes(color = Rank)) +
  geom_text_repel(aes(label = Term), size = 2) +
  scale_color_viridis(option = "plasma") +
  xlab("Odds Ratio") +
  ylab("-log10(q-value)") +
  #xlim(c(-150, 150)) +
  ylim(c(0, 8)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(color = "gray", fill = NA, size = 0.2))
ggsave("reports/figures/NeatRePrep_PASCfuUnique_ENRICHr_GOCC.pdf", 
       width = 18, height = 12, units = "cm")



## Contaminant Analysis ----
contaminant <- read.csv("data/metadata/contamination_proteins.csv") %>%
  mutate(single_id = word(Protein_IDs, 1, sep = ";")) %>%
  mutate(single_id = word(single_id, 1, sep = "-"))


# check for exact matches to contamination vector
contaminant1 <- contaminant %>%
  mutate(is = single_id %in% neat_data_df$PG.ProteinGroups) %>%
  mutate(single_id = case_when(
    single_id == "P0CG48" ~ "P0CG47;P0CG48;P62979;P62987",
    T ~ single_id
  )) %>% # this changes that id to the same as PG.ProteinGroups so there is a match
  mutate(is = single_id %in% spec_PG$PG.ProteinGroups) %>%
  filter(is == T) %>%
  mutate(PG.ProteinGroups = single_id) %>%
  dplyr::select(PG.ProteinGroups, Type)


contaminant_df <- neat_data_df %>%
  mutate(
    group = case_when(
      PG.ProteinGroups %in% 
        (contaminant1 %>% 
        filter(Type == "Erythrocyte") %>% 
        pull(PG.ProteinGroups)) ~ "erythrocyte",
      PG.ProteinGroups %in% 
        (contaminant1 %>% 
        filter(Type == "Platelet") %>% 
        pull(PG.ProteinGroups)) ~ "platelet",
      PG.ProteinGroups %in% 
        (contaminant1 %>% 
        filter(Type == "Coagulation") %>% 
        pull(PG.ProteinGroups)) ~ "coagulation",
      TRUE ~ "none" # Default group if no match is found
    )
  ) %>%
  filter(group != "none")


##* Platelet analysis ----
platelet <- contaminant1 %>%
  filter(Type == "Platelet") %>%
  pull(PG.ProteinGroups)

platelet_df <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[PG.ProteinGroups %in% platelet], na.rm = TRUE),
    ratio = sum_subset / sum_all
  ) %>%
  left_join(metadata, by = "sample_id")

platelet_df_mean_sd <- platelet_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(platelet_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = c(pal[4], pal[6])) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio + (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = platelet_df_mean_sd$median_ratio - (2 * platelet_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("platelet") +
  xlab('Sample ID') +
  ylab('sum(platelet proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.008, 0.02)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/NeatRePrep_PlateletContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_platelet <- run_ids_cor %>%
  left_join(platelet_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_platelet, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('platelet contamination index (higher = more platelets)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_PlateletContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')


## Erythrocyte analysis ----
erythrocyte <- contaminant1 %>%
  filter(Type == "Erythrocyte") %>%
  pull(PG.ProteinGroups)

erythrocyte_df <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[PG.ProteinGroups %in% erythrocyte], na.rm = TRUE),
    ratio = sum_subset / sum_all  
  ) %>%
  left_join(metadata, by = "sample_id")

erythrocyte_df_mean_sd <- erythrocyte_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(erythrocyte_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = c(pal[4], pal[6])) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio + (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = erythrocyte_df_mean_sd$median_ratio - (2 * erythrocyte_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("erythrocyte") +
  xlab('Sample ID') +
  ylab('sum(erythrocyte proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.01, 0.015)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/NeatRePrep_ErythrocyteContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_erythrocyte <- run_ids_cor %>%
  left_join(erythrocyte_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_erythrocyte, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('erythrocyte contamination index (higher = more erythrocyte)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_ErythrocyteContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')


## Coagulation analysis ----
coagulation <- contaminant1 %>%
  filter(Type == "Coagulation") %>%
  pull(PG.ProteinGroups)

coagulation_df <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(
    sum_all = sum(normalized_abundance, na.rm = TRUE),
    sum_subset = sum(normalized_abundance[PG.ProteinGroups %in% coagulation], na.rm = TRUE),
    ratio = sum_subset / sum_all 
  ) %>%
  left_join(metadata, by = "sample_id")

coagulation_df_mean_sd <- coagulation_df %>%
  ungroup() %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE)
  )


ggplot(coagulation_df, 
       aes(sample_id, ratio, color = Cohort)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_manual(values = c(pal[4], pal[6])) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio,
             linewidth = 0.2) +
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio + (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  geom_hline(yintercept = coagulation_df_mean_sd$median_ratio - (2 * coagulation_df_mean_sd$sd_ratio),
             linetype = 3,
             linewidth = 0.2) +  
  ggtitle("coagulation") +
  xlab('Sample ID') +
  ylab('sum(coagulation proteins)/sum(all proteins)') +
  scale_y_continuous(expand = c(0,0), limits = c(0.005, 0.009)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )
ggsave('reports/figures/NeatRePrep_CoagulationContamination.pdf', 
       width = 12, height = 6, units = 'cm')


## Correlation of contaminant index to protein ids

run_ids_coagulation <- run_ids_cor %>%
  left_join(coagulation_df %>%
              dplyr::select(sample_id, ratio),
            by = "sample_id")

ggplot(run_ids_coagulation, aes(ratio, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('coagulation contamination index (higher = more coagulation proteins)') +
  ylab('Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_CoagulationContamCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')







# Separate plots for each facet
p1 <- ggplot(contaminant_df %>% filter(group == "erythrocyte"), 
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("erythrocyte") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

p2 <- ggplot(contaminant_df %>% filter(group == "platelet"), 
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("platelet") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

p3 <- ggplot(contaminant_df %>% filter(group == "coagulation"), 
             aes(sample_id, normalized_abundance, color = PG.ProteinGroups)) +
  geom_point(size = 1,
             alpha = 1,
             stroke = 0) + 
  scale_color_viridis(discrete = T) +
  ggtitle("coagulation") +
  xlab('Sample ID') +
  ylab('log2 Abundance') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Proteins')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none',
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  )

# Combine the plots with patchwork
p1 / p2 / p3
ggsave('reports/figures/NeatRePrep_ContaminantProfiles.pdf', 
       width = 16, height = 12, units = 'cm')



## Check protein content correlation to Neat IDs or Seer IDs ----

protein_bca <- fread("data/metadata/HealthyPASC_fuBCA.csv")

bca_cor <- run_ids_cor %>%
  left_join(protein_bca, by = "Sample")

ggplot(bca_cor, aes(Concentration, count_seer, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('BCA Concentration') +
  ylab('Seer Count') +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_ProteinBCA_SeerCount_cohort_Point.pdf', 
       width = 10, height = 6, units = 'cm')


## plot:CollectionDatevsIDs ----

neat_run_ids <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  left_join(metadata, by = 'sample_id') %>%
  left_join(rawfiles, by = 'Sample') %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number()) %>%
  arrange(Collection_date) %>%
  mutate(collection_order = row_number()) %>%
  mutate(Collection_date = as.character(Collection_date))


ggplot(neat_run_ids, aes(collection_order, count, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('Collection Order') +
  ylab('Neat Count') +
  #scale_x_date(
  #  date_breaks = "1 year",  # Change as needed (day, week, month, year)
  #  minor_breaks = "1 month",
  #  date_labels = "%Y"     # Format: Jan 2021
  #) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_PGvsCollectionOrder_Set_Point.pdf', 
       width = 12, height = 6, units = 'cm')


# plot:BCAvsCollectionDate

neat_run_ids_bca <- neat_data_df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  left_join(metadata, by = 'sample_id') %>%
  left_join(rawfiles, by = 'Sample') %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number()) %>%
  arrange(Collection_date) %>%
  mutate(collection_order = row_number()) %>%
  mutate(Collection_date = as.character(Collection_date)) %>%
  left_join(protein_bca, by = "Sample")


ggplot(neat_run_ids_bca, aes(collection_order, Concentration, color = Cohort)) + 
  geom_point(size = 1) +
  scale_color_manual(values = c(pal[4], pal[6])) +
  xlab('Collection Order') +
  ylab('Neat Count') +
  #scale_x_date(
  #  date_breaks = "1 year",  # Change as needed (day, week, month, year)
  #  minor_breaks = "1 month",
  #  date_labels = "%Y"     # Format: Jan 2021
  #) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/NeatRePrep_CollectionOrdervsBCA_Set_Point.pdf', 
       width = 12, height = 6, units = 'cm')

