
library(data.table)
library(readxl)
library(ggplot2)
library(pheatmap)
library(purrr)
library(dplyr)
library(tidyr)
library(stats)
library(plotly)
library(ggfortify)
library(ggrepel)
library(reshape2)
library(scales)
library(Matrix)
library(lme4)
library(lmerTest)
library(broom)




QC_norm_matrix <- data.matrix(QC_expression_norm)


# List to store the subsets
batch_qc_list <- list()

# Loop through batch numbers from 1 to 12
for (batch_number in 1:12) {
  # Create a subset based on the column names containing the batch number
  batch_qc_list[[paste0("Batch", batch_number, "_QC")]] <- QC_expression_norm[, grepl(paste0("Batch", batch_number, "_"), colnames(QC_expression_norm))]
}

# Extract Batch#_QC for 12 batches from the list to separate data frames
for (batch_name in names(batch_qc_list)) {
  # Assign the data frame to a variable with a meaningful name
  assign(batch_name, batch_qc_list[[batch_name]])
}
  
rsd <- function(Name, index = 1:7789){sd(Name[index], na.rm = T)/mean(Name[index], na.rm = T) * 100}

QC_rsd <- apply(QC_expression_norm, 1, function(x) mean(unlist(x))/sd(unlist(x)))

QC_rsd <- data.frame(RSD= QC_rsd) #normalized


hist(QC_rsd$RSD)

# Create a list to store the batch-specific data frames
batch_QC_rsd <- list(Batch1_QC, Batch2_QC, Batch3_QC, Batch4_QC, Batch5_QC, Batch6_QC, Batch7_QC, Batch8_QC, Batch9_QC, Batch10_QC, Batch11_QC, Batch12_QC)


# Loop through each batch-specific data frame
for (i in 1:9) {
  # Calculate the RSD for the current batch and add it as a new column
  batch_QC_rsd[[i]]$rsd <- apply(batch_QC_rsd[[i]], 1, rsd)
  # Create a new data frame for the RSD values of the current batch
  assign(paste0("B", i, "_QC_rsd"), data.frame(RSD = batch_QC_rsd[[i]]$rsd))
}

for (i in 11:12) {
  # Calculate the RSD for the current batch and add it as a new column
  batch_QC_rsd[[i]]$rsd <- apply(batch_QC_rsd[[i]], 1, rsd)
  # Create a new data frame for the RSD values of the current batch
  assign(paste0("B", i, "_QC_rsd"), data.frame(RSD = batch_QC_rsd[[i]]$rsd))
  
}

all_batches_rsd <- cbind(QC_rsd, B1_QC_rsd, B2_QC_rsd, B3_QC_rsd,B4_QC_rsd,B5_QC_rsd,B6_QC_rsd,B7_QC_rsd,B8_QC_rsd,B9_QC_rsd,B11_QC_rsd,B12_QC_rsd)
colnames(all_batches_rsd) <- c("Combined","Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6", "Batch7", "Batch8", "Batch9", "Batch11", "Batch12")
write.csv(all_batches_rsd, "rsd_per_batch.csv", row.names = FALSE)


_________________________________________
##Another 2 codes for rsd calculation per batch to check####

###1) RSD All Batches####  



QC_rsd = apply(QC_data_normalized[, grep("Control|control|1x|1X|Pooled|pooled", names(QC_data_normalized))], 1, function(x) mean(unlist(x))/sd(unlist(x)))
QC_rsd = data.frame(QC_rsd)

data_filtered_normalized_rsd <- cbind (QC_data_normalized, QC_rsd)

rsd <- function(Name, index = 1:7789){sd(Name[index], na.rm = T)/mean(Name[index], na.rm = T) * 100}

data_filtered_normalized_rsd$QC_rsd =  apply(data_filtered_normalized_rsd[, grep("Control|control|1x|1X|Pooled|pooled", names(data_filtered_normalized_rsd))], 1, function(x) mean(unlist(x))/sd(unlist(x)))


QC_RSD1 <- apply(data_filtered_normalized_rsd[, grep("Control|control|1x|1X|Pooled|pooled", names(data_filtered_normalized_rsd))], 1, function(x) mean(unlist(x))/sd(unlist(x)))
QC_RSD1 = data.frame(QC_RSD1)



hist(data_filtered_normalized_rsd$QC_rsd, breaks=20)

median(data_filtered_normalized_rsd$QC_rsd)


Monkey_norm_rsd = data_filtered_normalized_rsd

write.csv(Monkey_norm_rsd, "Monkey_norm_rsd_QConly.csv", row.names = TRUE)


###2) QCs_per_Batch: RSD per Batch####


QC_expression_norm_subset = data_filtered_normalized_rsd[data_filtered_normalized_rsd$DetectedPeaks_in20.QC == TRUE, ]

write.csv(QC_expression_norm_subset, "QC_normalized_filtered_Features_Found_only_inQC.csv", row.names = FALSE)

QC_expression_norm_subset = QC_expression_norm_subset[, -c(1:64)]


# List to store the subsets
batch_qc_list <- list()

# Loop through batch numbers from 1 to 12
for (batch_number in 1:12) {
  # Create a subset based on the column names containing the batch number
  batch_qc_list[[paste0("Batch", batch_number, "_QC")]] <- data_filtered_normalized_rsd[, grepl(paste0("Batch", batch_number, "_"), colnames(data_filtered_normalized_rsd))]
}

# Extract Batch#_QC for 12 batches from the list to separate data frames
for (batch_name in names(batch_qc_list)) {
  # Assign the data frame to a variable with a meaningful name
  assign(batch_name, batch_qc_list[[batch_name]])
}



Batch1_QC$RSD = apply(Batch1_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch1_QC))], 1, rsd)

B1_QC_rsd = data.frame(RSD1 =Batch1_QC$RSD)

Batch2_QC$RSD = apply(Batch2_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch2_QC))], 1, rsd)
B2_QC_rsd = data.frame(RSD2 =Batch2_QC$RSD)

Batch3_QC$RSD = apply(Batch3_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch3_QC))], 1, rsd)
B3_QC_rsd = data.frame(RSD3 =Batch3_QC$RSD)

Batch4_QC$RSD = apply(Batch4_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch4_QC))], 1, rsd)
B4_QC_rsd = data.frame(RSD4 =Batch4_QC$RSD)

Batch5_QC$RSD = apply(Batch5_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch5_QC))], 1, rsd)
B5_QC_rsd = data.frame(RSD5 =Batch5_QC$RSD)

Batch6_QC$RSD = apply(Batch6_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch6_QC))], 1, rsd)
B6_QC_rsd = data.frame(RSD6 =Batch6_QC$RSD)

Batch7_QC$RSD = apply(Batch7_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch7_QC))], 1, rsd)
B7_QC_rsd = data.frame(RSD7 =Batch7_QC$RSD)

Batch8_QC$RSD = apply(Batch8_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch8_QC))], 1, rsd)
B8_QC_rsd = data.frame(RSD8 =Batch8_QC$RSD)

Batch9_QC$RSD = apply(Batch9_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch9_QC))], 1, rsd)
B9_QC_rsd = data.frame(RSD9 =Batch9_QC$RSD)

Batch11_QC$RSD = apply(Batch11_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch11_QC))], 1, rsd)
B11_QC_rsd = data.frame(RSD11 =Batch11_QC$RSD)

Batch12_QC$RSD = apply(Batch12_QC[, grep("Control|control|1x|1X|Pooled|pooled", names(Batch12_QC))], 1, rsd)
B12_QC_rsd = data.frame(RSD12 =Batch12_QC$RSD)



Name_filtered_norm = data.frame(Name=data_filtered_normalized_rsd$All.Features)
Name_index = data.frame(UniqueID = data_filtered_normalized_rsd$UniqueID)
all_batches_rsd <- cbind(Name_index, Name_filtered_norm, QC_rsd, B1_QC_rsd, B2_QC_rsd, B3_QC_rsd,B4_QC_rsd,B5_QC_rsd,B6_QC_rsd,B7_QC_rsd,B8_QC_rsd,B9_QC_rsd,B11_QC_rsd,B12_QC_rsd)
colnames(all_batches_rsd) <- c("UniqueID","Name", "Combined","Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6", "Batch7", "Batch8", "Batch9", "Batch11", "Batch12")
write.csv(all_batches_rsd, "rsd_per_batch.csv", row.names = FALSE)


rsd_per_batch_foundinQC = QC_expression_norm_subset[, grep("RSD|rsd", colnames(QC_expression_norm_subset))]

rsd_per_batch_foundinQC = cbind(data_filtered_normalized_rsd[data_filtered_normalized_rsd$DetectedPeaks_in20.QC == TRUE, ], rsd_per_batch_foundinQC) 

__________________________________________________________________________________________________



-------------------------
  
##3) Violin####

all_batches_rsd = all_batches_rsd[, -c(1,2)]

### RSD per Batch####  #no batch10
# Extract the batch names from the column names
batch_names <- colnames(all_batches_rsd)
all_batches_rsd <- cbind(Name_index, all_batches_rsd)

all_batches_rsd$UniqueID <- gsub("[A-Za-z]","", all_batches_rsd$UniqueID)

all_batches_rsd <- gather(all_batches_rsd, key= "Batch", value= "RSD", -UniqueID) #mutate data frame to a long one for ggplot

# Define a custom color palette for the 12 batches
# Replace "batch1", "batch2", ..., "batch12" with appropriate names for your batches
custom_colors <- c("grey","#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928",
                   "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")



# Map the sample order to the "Sample" column in pca_df
all_batches_rsd$Batch <- factor(all_batches_rsd$Batch , levels = unique(all_batches_rsd$Batch))


# Plot multiple violins for each batch with custom colors

ggplot(all_batches_rsd, aes(x = Batch, y = RSD, fill=Batch)) +
  geom_violin( color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  ggtitle("Distribution of RSDs for Different Batches (Normalized)") +
  xlab("Batch") +
  ylab("RSD (%)") +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
    axis.title = element_text(face = 'bold')
  )


median_values_rsd <- all_batches_rsd %>%
  group_by(Batch) %>%
  summarize(Median_RSD = median(RSD))

# Print the median values table
print(median_values_rsd)

median_values_rsd <- median_values_rsd[-c(1),]
Median_RSD <- median(median_values_rsd$Median_RSD)

print(Median_RSD)

write.csv(median_values_rsd, "median_values_rsd_beforecutoff.csv", row.names = TRUE)



----------------------------------
  
####RSD for all Batches after normalization####

ggplot(rsd_df_norm, aes(x = " ", y = RSD)) +
  geom_violin( fill = "skyblue") +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  ggtitle("Distribution of RSDs for all QCs (normalized)") +
  xlab("Batch") +
  ylab("RSD (%)") +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
    axis.title = element_text(face = 'bold')
  )

rsd_df1 <- data.frame(Monkey_norm_rsd$RSD)

####RSD for all batches/ filtered RSD####
ggplot(rsd_df1, aes(x = " ", y = Monkey_norm_rsd$RSD)) +
  geom_violin( fill = "maroon") +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  ggtitle("Distribution of RSDs for all QCs 85% cutoff (normalized & filtered)") +
  xlab("Batch") +
  ylab("RSD (%)") +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
    axis.title = element_text(face = 'bold')
  )

# .......................................................  ######
---------------------------------  
  #"FILTERED & NORMALIZED DATA SORTING"####
---------------------------------
  
  metabolite_rsd_filtered <- rsd_df1
metabolite_rsd_filtered$row_num <- rownames(metabolite_rsd_filtered)

Monkey_rsd_filtered <- Data_colnames_edited
Monkey_rsd_filtered$row_num <- seq_len(nrow(Data_colnames_edited))
Monkey_rsd_filtered <- merge(Monkey_rsd_filtered,metabolite_rsd_filtered, by= 'row_num') 
Monkey_rsd_filtered <- Monkey_rsd_filtered[,-c(1,653)] #remove row number and RSD columns


#Log2_intensity of filtered data
Monkey_rsd_filtered_log <- log2(Monkey_rsd_filtered[, 27:651])

----------------------------------------------------------------  
  # .......................................................  ######
-------------------------------------------
  #"PLOTTING NORMALIZED $ FILTERED DATA"####
# .......................................................  ######


##1) BoxPlot after RSD Filtering####


Monkey_rsd_filtered_log$Name <- row.names(Monkey_rsd_filtered_log) 
melted_expression_box <- melt(Monkey_rsd_filtered_log, id.vars = "Name")

# Create a new column to indicate if the sample is a QC or not
melted_expression_box$Sample_Type <- ifelse(grepl("Control|control|1x|1X|Pooled|pooled", melted_expression_box$variable), "QC", "Sample")

melted_expression_box$color <- ifelse(grepl("Control|control|1x|1X|Pooled|pooled", melted_expression_box$variable), "#67D2CA", "#556EB5")


#subset melted_expression_box
melted_expression_box1 <- melted_expression[1:859320,] # Batches 1,2
melted_expression_box2 <- melted_expression[859321:1718640,] # Batch4, Control batch7, Batch3 
melted_expression_box3 <- melted_expression[1718641:2546712,] # Batches 5,6
melted_expression_box4 <- melted_expression[2546713:3335724,] # Batches 7,8
melted_expression_box5 <- melted_expression[3335725:4062240,] # Batches 9,10
melted_expression_box6 <- melted_expression[4062241:4882500,] # Batches 11,12

#First Box plot  
ggplot(melted_expression_box1, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10)),
    plot.margin = margin(l =  20)
  )


#2nd Box plot  
ggplot(melted_expression_box2, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot2 of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))
  ) 

#3rd Box plot  
ggplot(melted_expression_box3, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot3 of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))
  ) 
#4th Box plot  
ggplot(melted_expression_box4, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot4 of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))
  ) 

#5th Box plot  
ggplot(melted_expression_box5, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot5 of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))
  ) 
#6th Box plot  
ggplot(melted_expression_box6, aes(x = variable, y = value, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("Boxplot6 of Log2 Intensity for Abundance of Metabolites") +
  ylim(1,40)+
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))
  ) 




