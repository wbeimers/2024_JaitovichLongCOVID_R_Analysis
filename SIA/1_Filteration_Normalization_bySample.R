
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


-------------------------------------
#ReadData###
-------------------------------------
  
###*Read Monkey Data*####  
Monkey <- read.csv("20230223_Anderson_Monkey_IntensityOnly.csv", header= TRUE)

Expression = Monkey[, -c(1:64)]
Metadata = read.csv("Final_metadata1.csv", header=TRUE)
Data_information = Monkey[, c(1:64)]



Name <- Monkey$UniqueID #Features Name stored
metabolites_features <- Monkey$Metabolites.Name..CD33.
lipids_features <- Monkey$Lipid.Annotation..LipiDex2.

#Edit colnames for original data (manually removed Group Area: from excel/same as Filename in Batches Metadata)
Data_colnames_edited <- Monkey
expression_colnames_edited <- Data_colnames_edited[, -c(1:64)]
expression_log2_colnames_edited <- log2(expression_colnames_edited)


expression <- expression_colnames_edited
expression_log2= log2(expression)
expression_log2_matrix= as.matrix(expression_log2)



Final_metadata = read.csv("Final_metadata1.csv", header = TRUE)

Final_metadata = Metadata

#Read Run_order file (added X in front of dates to fix an error)
runorder <- read.csv("RunOrderMetadata - copy.csv", header= TRUE)


ordered_filenames <- runorder$Filename


# Subset and reorder columns in the expression data frame based on runorder

expression_reordered <- expression_log2_colnames_edited %>%
  select(one_of(ordered_filenames))



###*Subset QC Data*####


Monkey_QC <- subset(Monkey, select = grep("UniqueID|Control|control|1x|1X|Pooled|pooled", colnames(Monkey), ignore.case = TRUE))


# Fix QCs in Batches 8 and 10

# Identify the QC sample filenames from Final_metadata
qc_filenames <- Final_metadata$Sample_Type[Final_metadata$Sample_Type == "Pooled Monkey Serum"]

# Subset the Monkey_QC data frame to include only QC columns
Monkey_QC_expression <- expression_colnames_edited
QC_expression_log2 <- expression_log2_colnames_edited


-------------------------------------------------------------------------------------------------------------------------------------------------------

###RawFilteration_based on #Detected Peaks after Gap Filling"


QC_RSD1 <- apply(Monkey[, grep("Control|control|1x|1X|Pooled|pooled", names(Monkey))], 1, function(x) mean(unlist(x))/sd(unlist(x)))
                            
Monkey_Filtered = Monkey$DetectedPeaks_.10Sample == TRUE & !(QC_RSD1 > 85 & Monkey$DetectedPeaks_in20.QC == TRUE)
                            
Monkey_Filtered <- Monkey[Monkey$DetectedPeaks_.10Sample == TRUE & !(QC_RSD1 > 85 & Monkey$DetectedPeaks_in20.QC == TRUE), ]
Expression_Filtered <- Monkey_Filtered[, -c(1:64)]
Data_information_filtered = Monkey_Filtered[, c(1:64)]

                            
########################################################

#Read Run_order file (added X in front of dates to fix an error)
runorder <- read.csv("RunOrderMetadata - copy.csv", header= TRUE)


ordered_filenames <- runorder$Filename


# Subset and reorder columns in the expression data frame based on runorder

expression_reordered <- Monkey_Filtered %>%
  select(one_of(ordered_filenames))

expression_reordered_log2 <- log2(expression_reordered)

expression_reordered_log2 = cbind(Monkey_Filtered$UniqueID, expression_reordered_log2)

names(expression_reordered_log2)[names(expression_reordered_log2) == "Monkey_Filtered$UniqueID"] <- "UniqueID"
###*Subset QC Data*####

Monkey_QC <- subset(Monkey_Filtered, select = grep("UniqueID|Control|control|1x|1X|Pooled|pooled", colnames(Monkey_Filtered), ignore.case = TRUE))

Monkey_QC_runorder <- subset(expression_reordered, select = grep("UniqueID|Control|control|1x|1X|Pooled|pooled", colnames(expression_reordered), ignore.case = TRUE))


# Subset the Monkey_QC data frame to include only QC columns
Monkey_QC_expression <- Monkey_QC

Monkey_QC_expression_runorder <- Monkey_QC_runorder

QC_expression_log2 <- log2(Monkey_QC[,-c(1,2)])
QC_expression_runorder <- log2(Monkey_QC_expression_runorder[,-c(1,2)])  

---------------------------------------------------------------
#Normalization####
#"NORMALIZATION"####
# .......................................................  ###### 

##1) BoxPlot for QCs - Run Order####


QC_expression_log2_runorder <- QC_expression_runorder  %>%
  select(matches("Control|control|1x|1X|Pooled|pooled"))


-----------
  
#create a melted data
  QC_expression_log2_runorder$UniqueID <- Monkey_Filtered$UniqueID 
  
  batch_box_QC <- character(ncol(QC_expression_log2_runorder))
  batch_box_QC[1:12] <- "Batch1"
  batch_box_QC[13:22] <- "Batch2"
  batch_box_QC[23:33] <- "Batch4"
  batch_box_QC[34:43] <- "Batch3"
  batch_box_QC[44:52] <- "Batch5"
  batch_box_QC[53:61] <- "Batch6"
  batch_box_QC[62:70] <- "Batch7"
  batch_box_QC[71:72] <- "Batch8"
  batch_box_QC[73:81] <- "Batch9"
  batch_box_QC[82] <- "Batch10"
  batch_box_QC[83:95 ] <- "Batch11"
  batch_box_QC[96:107] <- "Batch12"
  
  
  
  melted_expression_box_QC_runorder <- melt(QC_expression_log2_runorder, id.vars = "UniqueID")
  
  melted_expression_box_QC_runorder <- melted_expression_box_QC_runorder %>%
    mutate(Batch = ifelse(grepl("Batch\\d+", variable, ignore.case = TRUE),
                          sub(".*(Batch\\d+).*", "\\1", variable),
                          NA))
  
  
  
  QC_Box_colors <- c("Batch1" = "#FFDDAA", "Batch2" = "#C3E4A6", "Batch3" = "#00CC99",
                     "Batch4" = "#9DAE9F", "Batch5" = "#DDDDCC", "Batch6" = "#A2A9B9",
                     "Batch7" = "#99CCDF", "Batch8" = "#0099CC", "Batch9" = "#415A87",
                     "Batch10" = "#FFDDEE", "Batch11" = "#CC0099", "Batch12" = "#0055AA")
  
  
  
  # Assuming 'melted_expression_box_QC_runorder' is your data frame
  
  # Define the order of levels for th
  
  batch_order <- paste0("Batch", 1:12)
  
  # Reorder the levels of the Batch factor
  melted_expression_box_QC_runorder$Batch <- factor(
    melted_expression_box_QC_runorder$Batch,
    levels = batch_order
  )
  
  
  
  
  
###########################################################  
  #Normalize by sample checlk (didn't use)
  
  batch_sample_medians <- melted_expression_box_QC_runorder %>%
    group_by(Batch, variable) %>%
    summarize(sample_median = median(value))
  
  # Calculate the median of sample medians for each batch
  batch_medians <- batch_sample_medians %>%
    group_by(Batch) %>%
    summarize(batch_median = median(sample_median))
  
  # View the resulting batch medians
  print(batch_medians)
 ############################################################## 

  
  
  ggplot(melted_expression_box_QC_runorder, aes(x = variable, y = value, fill= Batch)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    xlab("Samples") +
    ylab("Log2 Intensity") +
    ggtitle("Boxplot of Log2 Intensity for Abundance of Metabolites") +
    ylim(1,40)+
    scale_fill_manual(values = QC_Box_colors) +  # Use the custom color palette
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
  
  
  median_values_QC <- melted_expression_box_QC_runorder %>%
    group_by(Batch) %>%
    summarize(Median_QC = median(value))
  
  
  # Print the median values table
  print(median_values_QC)
  
  batch_order <- c("Batch1", "Batch2", "Batch4", "Batch3", "Batch5", "Batch6", "Batch7", "Batch8", "Batch9", "Batch10", "Batch11", "Batch12")
  
  # Reorder the levels of the Batch factor according to the desired order
  median_values_QC <- median_values_QC %>%
    arrange(match(Batch, batch_order))
  
  
  # Create a plot of median values with a best fitting line
  ggplot(median_values_QC, aes(x = as.numeric(factor(Batch)), y = Median_QC)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "#fdbf6f", linetype = "dashed") +
    ylim(16,20)+
    xlab("Run Order") +
    ylab("Median of Log2 Intensity") +
    ggtitle("Median Log2 Intensity for Each Batch") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold", margin = margin(b = 30),hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          axis.text.x = element_blank(),   # Remove x-axis text
          axis.line = element_line(color = "black", size = 0.8),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  
  ##2) Grand Median for QCs####
  
  Grand_median <- median(median_values_QC$Median_QC)
  print(Grand_median)
  
  
  ## 3) Correction Value####  
  within_Batch_correction_Value <-  median_values_QC$correction_value <- (median_values_QC$Median_QC) - Grand_median #within_batch_median(Median_QC) - Grand Median 
  
  Norm_table = median_values_QC
  
  
  library(stringr)
  
  batch_data_frames <- list()
  
  # Loop through each batch number (assuming 12 batches)
  for (batch_number in 1:12) {
    # Initialize a data frame with the same number of rows as expression_log2
    batch_data <- data.frame(matrix(NA, nrow = nrow(expression_reordered_log2), ncol = 0))
    
    # Loop through all column names
    for (col_name in colnames(expression_reordered_log2)) {
      # Split the column name by underscores
      col_parts <- unlist(strsplit(col_name, "_"))
      
      # Check if the column name contains "Batch" and the batch number
      if (length(col_parts) >= 3 && col_parts[3] == paste0("Batch", batch_number)) {
        # Extract and add the column to the batch_data data frame
        batch_data[[col_name]] <- expression_reordered_log2[[col_name]]
      }
    }
    
    # Store the batch-specific data frame in the list
    batch_data_frames[[paste0("expression_log2_batch", batch_number)]] <- batch_data
  }
  
 
  # Access batch-specific data frames, e.g., expression_log2_batch1, expression_log2_batch2, ...
  expression_log2_batch1 <- batch_data_frames[["expression_log2_batch1"]]
  expression_log2_batch2 <- batch_data_frames[["expression_log2_batch2"]]
  expression_log2_batch3 <- batch_data_frames[["expression_log2_batch3"]]
  expression_log2_batch4 <- batch_data_frames[["expression_log2_batch4"]]
  expression_log2_batch5 <- batch_data_frames[["expression_log2_batch5"]]
  expression_log2_batch6 <- batch_data_frames[["expression_log2_batch6"]]
  expression_log2_batch7 <- batch_data_frames[["expression_log2_batch7"]]
  expression_log2_batch8 <- batch_data_frames[["expression_log2_batch8"]]
  expression_log2_batch9 <- batch_data_frames[["expression_log2_batch9"]]
  expression_log2_batch10 <- batch_data_frames[["expression_log2_batch10"]]
  expression_log2_batch11 <- batch_data_frames[["expression_log2_batch11"]]
  expression_log2_batch12 <- batch_data_frames[["expression_log2_batch12"]]
  
  
  ## 4) Within-batch-normalization####
  
  for (batch_number in 1:12) {
    # Get the correction value for the current batch from median_values_QC
    correction_value <- median_values_QC$correction_value[median_values_QC$Batch == paste0("Batch", batch_number)]
    # Multiply the correction value with the corresponding expression_log2_batch data frame
    expression_log2_batch_name <- paste0("expression_log2_batch", batch_number)
    if (exists(expression_log2_batch_name)) {
      assign(expression_log2_batch_name, get(expression_log2_batch_name) - correction_value)
    }
  }
  
  # gather data frames into one
  batch_data_list <- list()
  
  # Loop through each batch number (assuming 12 batches)
  for (batch_number in 1:12) {
    # Get the name of the current batch data frame
    expression_log2_batch_name <- paste0("expression_log2_batch", batch_number)
    
    # Check if the data frame exists
    if (exists(expression_log2_batch_name)) {
      # Add the current batch data frame to the list
      batch_data_list[[expression_log2_batch_name]] <- get(expression_log2_batch_name)
    }
  }
  
  # Combine the list of batch data frames into a single data frame without changing column names
  expression_log2_normalized_runorder <- dplyr::bind_cols(batch_data_list)
  
  #Reorder after the loop#
  
  expression_log2_normalized_runorder <- expression_log2_normalized_runorder %>%
    select(one_of(ordered_filenames))
  
  
  ## 5) Recalculating within-batch Median and Grand Median for QCs####
  expression_QC_log_norm <- expression_log2_normalized_runorder  %>%
    select(matches("Control|control|1x|1X|Pooled|pooled"))
  
  
  
  

  ##6) Normalized new data frames####
  
  expression_normalized <- 2^expression_log2_normalized_runorder
  data_normalized <- cbind(Data_information_filtered, expression_normalized)


  QC_expression_norm <- 2^expression_QC_log_norm
  QC_data_normalized <- cbind(Data_information_filtered, QC_expression_norm)
  

  QC_data_normalized_log <- cbind(Data_information_filtered, expression_QC_log_norm)
  

  
  write.csv(Norm_table, "normalization_table.csv", row.names = FALSE)
  write.csv(data_normalized, "FullMonkeydata_normalized.csv", row.names = FALSE)
  write.csv(QC_data_normalized_log, "Filtered_Normalized_QC_Log2.csv", row.names = FALSE)  
  write.csv(QC_data_normalized, "Filtered_Normalized_QC_expression.csv", row.names = FALSE)
  write.csv(FilteredMonkey_normalized, "expression_normalized_AllrawData.csv", row.names = FALSE)
  
  
  ## 7) Boxplot after normalization####
  
  #create a melted data
  
  expression_box_QC_norm <- character(ncol(expression_QC_log_norm))
  expression_QC_log_norm$Name <- Monkey_Filtered$UniqueID 
  
  expression_box_QC_norm <- data.frame(expression_box_QC_norm)
  
  melted_expression_box_QC_norm <- melt(expression_QC_log_norm, id.vars = "Name")
  
  melted_expression_box_QC_norm <- melted_expression_box_QC_norm %>%
    mutate(Batch = ifelse(grepl("Batch\\d+", variable, ignore.case = TRUE),
                          sub(".*(Batch\\d+).*", "\\1", variable),
                          NA))
  
  
  
  QC_Box_colors <- c("Batch1" = "#FFDDAA", "Batch2" = "#C3E4A6", "Batch3" = "#00CC99",
                     "Batch4" = "#9DAE9F", "Batch5" = "#DDDDCC", "Batch6" = "#A2A9B9",
                     "Batch7" = "#99CCDF", "Batch8" = "#0099CC", "Batch9" = "#415A87",
                     "Batch10" = "#FFDDEE", "Batch11" = "#CC0099", "Batch12" = "#0055AA")
  
  
  
  # Assuming 'melted_expression_box_QC_norm' is your data frame
  
  # Define the order of levels for th
  
  batch_order <- paste0("Batch", 1:12)
  
  # Reorder the levels of the Batch factor
  melted_expression_box_QC_norm$Batch <- factor(
    melted_expression_box_QC_norm$Batch,
    levels = batch_order
  )
  
  
  ggplot(melted_expression_box_QC_norm, aes(x = variable, y = value, fill= Batch)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    xlab("Samples") +
    ylab("Log2 Intensity") +
    ggtitle("Boxplot of Log2 Intensity for Abundance of Metabolites (Normalized)") +
    ylim(1,40)+
    scale_fill_manual(values = QC_Box_colors) +  # Use the custom color palette
    theme(
      axis.text.x = element_blank(),
      legend.position = "top",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.1, margin = margin(b = 1)),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.title.align = 0.5,  # Center align the legend title
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10)),
      plot.margin = margin(l =  20)
    )
  
  median_QC_after_norm <- melted_expression_box_QC_norm %>%
    group_by(Batch) %>%
    summarize(Median_QC = median(value))
  
  # Print the median values table
  print(median_QC_after_norm)
  
  batch_order <- c("Batch1", "Batch2", "Batch4", "Batch3", "Batch5", "Batch6", "Batch7", "Batch8", "Batch9", "Batch10", "Batch11", "Batch12")
  
  # Reorder the levels of the Batch factor according to the desired order
  median_QC_after_norm <- median_QC_after_norm %>%
    arrange(match(Batch, batch_order))
  
  # Create a plot of median values with a best fitting line
  ggplot(median_QC_after_norm, aes(x = as.numeric(factor(Batch)), y = Median_QC)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "maroon", linetype = "dashed") +
    ylim(16,20)+
    xlab("Run Order") +
    ylab("Median of Log2 Intensity") +
    ggtitle("Median Log2 Intensity for Each Batch after Normalization") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold", margin = margin(b = 30),hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          axis.text.x = element_blank(),   # Remove x-axis text
          axis.line = element_line(color = "black", size = 0.8),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  
  #################################################
  
  
  V1000_beforeNorm <- Monkey_Filtered[Monkey_Filtered$UniqueID == "V100", grep("Control|control|1x|1X|Pooled|pooled", names(Monkey_Filtered), ignore.case = TRUE)]
  V1000_afterNorm <- data_normalized[data_normalized$UniqueID == "V100", grep("Control|control|1x|1X|Pooled|pooled", names(data_normalized), ignore.case = TRUE)]
  
   
  # Combine the intensities into a single data frame or matrix
  all_intensities <- data.frame(
    Original = V1000_beforeNorm,
    Normalized = V1000_afterNorm
  )
  
  # Plot box plot for all intensities
  boxplot(all_intensities, col = c("blue", "red"), main = "Intensities Before and After Normalization",
          xlab = "Normalization", ylab = "Intensity", names = c("Original", "Normalized"))
  