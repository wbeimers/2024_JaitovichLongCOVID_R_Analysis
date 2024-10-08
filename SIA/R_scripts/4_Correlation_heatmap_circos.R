# libraries
library(corrplot)
library(RColorBrewer)
library(circlize)

#read data from volcano plot and metadata:
Linear_model_data_noQC = read.csv("Linear_model_data_noQC_normalized_by_feature.csv", header = TRUE)
Volcano_data = read.csv("20230315_Volcano_diet_simple.csv", header=TRUE)
significant_features = read.csv("20230315_NamedSignificantFeatures_diet_simple.csv", header=TRUE)
significant_features = significant_features[,-c(1)]


metadata_linearModel = Linear_model_data_noQC[,c(1:21)]

# Filter result_simple_model based on significant features
significant_results <- Linear_model_data_noQC[,colnames(Linear_model_data_noQC) %in% significant_features$Intensity_Column ]
significant_results = cbind(metadata_linearModel, significant_results)
test = significant_results[,c("MonkeyID","PullDate", "Diet", "Sex", "SampleName")]
test2 = significant_results[, c(22:75)]



# Calculate correlation matrix
correlation_matrix <- cor(test2, method = "pearson")
filtered_matrix <- correlation_matrix


# Filter out correlations below the threshold value
threshold <- 0.6

filtered_matrix[abs(filtered_matrix) < threshold] <- 0

# Omit repeated values and keep only the lower triangular part
filtered_matrix[upper.tri(filtered_matrix)] <- 0

# Set diagonal elements to zero
diag(filtered_matrix) <- 0

# Perform hierarchical clustering on the filtered correlation matrix
hc_rows <- hclust(as.dist(1 - abs(filtered_matrix)), method = "complete")
row_order <- hc_rows$order

# Define color palette for positive and negative correlations
positive_color <- "steelblue"
negative_color <- "maroon"

# Assign colors based on the sign of the correlation
linkcolors <- ifelse(filtered_matrix > 0, positive_color, negative_color)

column_mapping <- setNames(significant_features$Feature, significant_features$Intensity_Column)

# Rename columns and rows
colnames(filtered_matrix) <- column_mapping[colnames(correlation_matrix)]
rownames(filtered_matrix) <- column_mapping[colnames(correlation_matrix)]

order= c(" Docosahexaenoic acid ethyl ester_peaksplit_avg",
         " Deoxycholic acid isomer_(needs confirmation)",
         " Testosterone glucuronide",
         "LysoPC 20:3  ",
         "LysoPC 20:5  ",
         "LysoPC 20:3  ",
         "LysoPE 16:1  ",
         "LysoPE 18:2  ",
         "LysoPE 20:3  ",
         "LysoPI 20:3  ",
         "PC 30:1  ",
         "PC 32:1  ",
         "PC 32:3  ",
         "PE 34:3  ",
         "PC 35:5  ",
         "PE 16:0_18:2  ",
         "PI 18:1_18:0  ",
         "Plasmenyl-PC P-34:2  ",
         "Plasmenyl-PC P-36:1  ",
         "Plasmenyl-PE P-34:2  ",
         "Plasmenyl-PE P-38:4  ",
         "DG 32:1  ",
         "DG 34:3  ",
         "DG 34:3  ",
         "DG 34:2  ",
         "DG 16:0_18:1  ",
         "DG 16:0_18:1  ",
         "DG 18:3_18:2  ",
         "DG 18:2_18:1  ",
         "DG 18:2_18:1  ",
         "DG 36:2  ",
         "DG 36:2  ",
         "DG 18:1_18:0  ",
         "TG 48:3  ",
         "Alkanyl-TG O-16:0_16:0_18:2  ",
         "TG 50:4  ",
         "TG 50:4  ",
         "TG 50:3  ",
         "TG 50:3  ",
         "TG 16:0_18:1_18:1  ",
         "TG 51:3  ",
         "TG 51:2  ",
         "TG 16:0_18:1_18:0  ",
         "TG 16:0_18:1_18:0  ",
         "TG 52:5  ",
         "TG 53:2  ",
         "TG 54:2  ",
         "TG 54:1  ",
         "TG 54:4  ",
         "TG 58:6  ",
         "TG 56:2  ",
         "TG 56:4  ",
         "TG 56:3  ",
         "SM d40:3  "
)



# Plot the circular plot with colors based on correlation sign and clustering
chordDiagram(filtered_matrix, transparency = 0.5, col = linkcolors, order= unique(order))

write.csv(filtered_matrix,"circos_filtered_matrix.csv", row.names = FALSE)

filtered_matrix = read.csv("circos_filtered_matrix.csv", header=TRUE)
column_order <- order(colnames(filtered_matrix))



############## Heat Map Correlation #############################

corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black")

hc <- hclust(as.dist(1 - correlation_matrix))

# Reorder the correlation matrix based on the hierarchical clustering
reordered_corr <- correlation_matrix[hc$order, hc$order]

# Plot the reordered correlation matrix with clustering
corrplot(reordered_corr, method = "color", type = "upper", tl.col = "black")

