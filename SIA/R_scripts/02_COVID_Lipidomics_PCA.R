
############################    *Load Libraries* ####################################

library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)
library(plotly)
library(ggfortify)
library(ggrepel)
library(reshape2)
library(scales)
library(Matrix)


############################    *Read Data* ####################################

data_samples <- read.csv("D:/LongCOVID/Important Tables/Expression_QuantifiableOnly/20250207_longCOVID_SamplesOnly_QuantifiableOnly_RunOrder.csv", header= TRUE)
Metadata = read.csv("D:/LongCOVID/Metadata/Metadata_All_lipidomics.csv", header=TRUE)
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)


############################    *Color Scheme* ####################################

batch_colors <- c("Batch1" = "#FFDDAA", "Batch2" = "#C3E4A6", "Batch3" = "#00CC99",
                  "Batch4" = "#9DAE9F", "Batch5" = "#DDDDCC", "Batch6" = "#A2A9B9",
                  "Batch7" = "#99CCDF", "Batch8" = "#0099CC", "Batch9" = "#415A87",
                  "Batch10" = "#FFDDEE")


cohort_colors <- c(
  "Acute" = "#EE6677",
  "Acute_fu" = "#AA3377",
  "Acute_NC" = "#CCBB44",
  "Healthy" = "#228833",
  "PASC" = "#66CCEE",
  "PASC_fu" = "#4477AA"
)

print(unique(Metadata$Cohort))

############################ *PCA - Check Before Normalization* ####################################


UniqueID = data.frame(UniqueID = data_samples$UniqueID)


filenames <- Metadata %>%
  filter(Metadata$run_type == "Sample", keep == 1) %>%
  pull(rawfile_name_R)  

data_samples <- data_samples%>%
  select(UniqueID,any_of(filenames))

Expression_log2 <-log2(data_samples[,-c(1)]) 

# Perform PCA on the transposed data
transposed_expression <- t(Expression_log2) 
colnames(transposed_expression) <- UniqueID$UniqueID



#PCA####

pca_result <- prcomp(transposed_expression, scale. = TRUE) 

# Extract PCA scores for the first two principal components (PC1 and PC2)
pca_scores <- as.data.frame(pca_result$x[, 1:2])

# #Metadata for PCA
pca_scores$rawfile_name_R <- rownames(pca_scores)

pca_scores <- merge(pca_scores, Metadata, by = "rawfile_name_R")


# Extract the explained variances for PC1 and PC2
var <- pca_result$sdev^2
pc1_var <- var[1]
pc2_var <- var[2]

# Calculate the percentages of variance explained
total_variance <- sum(var)
pc1_percentage <- percent(pc1_var/ total_variance, accuracy = 0.1)
pc2_percentage <- percent(pc2_var / total_variance, accuracy = 0.1)

pca_scores$batch <- factor(pca_scores$batch, levels = unique(pca_scores$batch))

######## Batch #####
pdf("D:/LongCOVID/PCA/PCA_Batch_Quuantifiable.pdf", width = 15, height = 7)
ggplot(pca_scores, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("PCA Plot- LongCOVID") +
  xlab(paste("PC1 (", pc1_percentage, "%)", sep = "")) +
  ylab(paste("PC2 (", pc2_percentage, "%)", sep = "")) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black", linewidth = 0.9),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
    )
dev.off()
                             



######## Cohort #####

pca_scores$Cohort <- factor(pca_scores$Cohort, levels = names(cohort_colors))

pdf("D:/LongCOVID/PCA/PCA_Cohort_Quantifiable.pdf", width = 15, height = 7)
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust point size and transparency
  ggtitle("PCA Plot - long COVID") +
  xlab(paste("PC1 (", pc1_percentage, "%)", sep = "")) +
  ylab(paste("PC2 (", pc2_percentage, "%)", sep = "")) +
  scale_color_manual(values = cohort_colors) +  # Apply cohort colors
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black", linewidth = 0.9),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )
dev.off()


