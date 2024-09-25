####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(impute)
library(viridis)

#Colors#
#Make a classic palette
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

#Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),"white", colorRampPalette(col1)(100))

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files#
spec_PG_NPs <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")
spec_PG_NPs_50 <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs_50%.csv")
spec_PG_NPs_50_tf_imp <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs_50%_tf_imp.csv")
groups <- read.csv("data/metadata/sample_groups.csv")
file_info <- read.csv("data/metadata/AllPlates_Samples_file_info.csv")
spec_PG_NPs <- spec_PG_NPs[-1]
spec_PG_NPs_50 <- spec_PG_NPs_50[-1]


####Heatmap####
#Axes are proteins on one side and samples with hierarchical clustering on the other
expression_matrix <- spec_PG_NPs
rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
#log2transform if not done yet
expression_matrix <- log2(expression_matrix)
#fix colnames
n <- colnames(expression_matrix)
n <- sapply(strsplit(n, "_"), function(x) paste(x[5], collapse = "_"))
colnames(expression_matrix) <- n
#Make annotation dataframe for the sample groups
sample_annot <- groups
row.names(sample_annot) <- sample_annot$Sample
sample_annot <- sample_annot[2]

#make annotation dataframe for the genes
heatmap_annot <- human_gene_list %>%
  dplyr::select(c("Entry", "Gene.Ontology..cellular.component."))
heatmap_annot <- semi_join(heatmap_annot, twentyfive_NPs_pgs_f_imputed_log2,
                           by = join_by("Entry" == "allgenes"))
rownames(heatmap_annot) <- heatmap_annot[,1]
heatmap_annot <- mutate(heatmap_annot, Intracellular = ifelse(grepl("cytosol|cytoplasm", Gene.Ontology..cellular.component.), "Yes", "No"))
heatmap_annot <- mutate(heatmap_annot, Membrane = ifelse(grepl("membrane", Gene.Ontology..cellular.component.), "Yes", "No"))
heatmap_annot <- mutate(heatmap_annot, Plasma = ifelse(grepl("plasma", Gene.Ontology..cellular.component.), "Yes", "No"))
heatmap_annot <- heatmap_annot[,-1]
heatmap_annot <- heatmap_annot[,-1]

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
         annotation_colors = list(Set = c(Acute = col[1], Acute_fu = col[2], Acute_NC = col[3], Healthy = col[4], PASC = col[5], PASC_fu = col[6])),
         filename = "reports/figures/NoFilter_heatmap_AllPlates.png",
         width = 8,
         height = 6)


#Get clustering Information for proteins
hm <- pheatmap(expression_matrix)
row_cluster <- data.frame(cluster = cutree(hm$tree_row, k = 5))
row_cluster_3 <- filter(row_cluster, row_cluster$cluster == 3)
row_cluster_3 <- rownames(row_cluster_3)


colors <- viridis(100, direction = 1, option = "magma")
image(matrix(1:length(colors), nrow = 1), 
      col = colors, 
      axes = FALSE, 
      main = "Color Visualization")
dev.off()

