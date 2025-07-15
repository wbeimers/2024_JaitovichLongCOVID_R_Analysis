
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
library(stringr)


############################    *Read Data* ####################################

data_samples <- read.csv("D:/LongCOVID/Important Tables/DataSubsets/20250207_longCOVID_AllSamples_RunOrder.csv", header= TRUE)
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)

############################  *Sorting Tables and Assigning Variables* ####################################

UniqueID = data.frame(UniqueID = data_samples$UniqueID)

filenames <- runorder %>%
  filter(runorder$run_type == "Sample", keep == 1) %>%
  pull(rawfile_name_R)  

data_samples <- data_samples%>%
  select(UniqueID,any_of(filenames))

log2_Expression <-log2(data_samples[,-c(1)]) 


batch_colors <- c("Batch1" = "#FFDDAA", "Batch2" = "#C3E4A6", "Batch3" = "#00CC99",
                  "Batch4" = "#9DAE9F", "Batch5" = "#DDDDCC", "Batch6" = "#A2A9B9",
                  "Batch7" = "#99CCDF", "Batch8" = "#0099CC", "Batch9" = "#415A87",
                  "Batch10" = "#FFDDEE")

cohort_colors <- c(
  "Acute covid" = "#EE6677",
  "Acute f/u (no PASC)" = "#AA3377",
  "Acute non-COVID" = "#CCBB44",
  "Healthy" = "#228833",
  "PASC" = "#66CCEE",
  "PASC f/u" = "#4477AA"
)

############################  *Batch Median - Befotr Norm* #########################

melted_data <- data_samples %>%
  pivot_longer(
    cols = -UniqueID,  # Keep UniqueID as an identifier, convert others to long format
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Batch = ifelse(grepl("Batch\\d+", variable, ignore.case = TRUE),
                   sub(".*(Batch\\d+).*", "\\1", variable),
                   NA)
  )

sample_medians <- melted_data %>%
  group_by(Batch, variable) %>%
  summarize(sample_median = median(log2(value)))

# Calculate the median of sample medians for each batch
batch_medians <- sample_medians %>%
  group_by(Batch) %>%
  summarize(batch_median = median(sample_median))

############################  *Grand Median* #########################

Grand_median <- median(batch_medians$batch_median)

##########################  *Correction Values* #######################

batch_medians$correction_value <- (batch_medians$batch_median) - Grand_median   #within_batch_median(Median_QC) - Grand Median 

##########################  *Splitting data into different Batch data frames  * #######################

split_batches <- function(log2_Expression, num_batches = 10) {
  batch_data_frames <- list()
  
  # Loop through each batch number
  for (batch_number in 1:num_batches) {
    # Initialize an empty data frame with the same number of rows as log2_Expression
    batch_data <- data.frame(matrix(NA, nrow = nrow(log2_Expression), ncol = 0))
    
    # Loop through column names
    for (col_name in colnames(log2_Expression)) {
      # Split the column name by underscores
      col_parts <- unlist(strsplit(col_name, "_"))
      
      # Check if the column belongs to the current batch
      if (length(col_parts) >= 3 && col_parts[3] == paste0("Batch", batch_number)) {
        batch_data[[col_name]] <- log2_Expression[[col_name]]
      }
    }
    
    # Store the batch-specific data frame in the list
    batch_data_frames[[paste0("log2_batch", batch_number)]] <- batch_data
  }
  
  # Assign batch data frames to global environment
  list2env(batch_data_frames, envir = .GlobalEnv)
  
  return(batch_data_frames)  # Return list of batch data frames
}

batch_data_frames <- split_batches(log2_Expression, num_batches = 10)

##########################  *Normalization* #######################

apply_batch_correction <- function(batch_medians, num_batches = 10) {
  for (batch_number in 1:num_batches) {
    # Get the correction value for the current batch
    correction_value <- batch_medians$correction_value[batch_medians$Batch == paste0("Batch", batch_number)]
    
    # Construct the batch data frame name
    log2_batch_name <- paste0("log2_batch", batch_number)
    
    # Check if the data frame exists
    if (exists(log2_batch_name, envir = .GlobalEnv)) {
      # Apply correction
      assign(log2_batch_name, get(log2_batch_name, envir = .GlobalEnv) - correction_value, envir = .GlobalEnv)
    }
  }
}


apply_batch_correction(batch_medians, num_batches = 10)

##########################  *Combining normalized batches into one* #######################

combine_batches <- function(num_batches = 10) {
  batch_data_list <- list()
  
  # Loop through each batch number
  for (batch_number in 1:num_batches) {
    log2_batch_name <- paste0("log2_batch", batch_number)
    
    # Check if the batch data frame exists
    if (exists(log2_batch_name, envir = .GlobalEnv)) {
      batch_data_list[[log2_batch_name]] <- get(log2_batch_name, envir = .GlobalEnv)
    }
  }
  
  # Combine the list of batch data frames into one without changing column names
  expression_log2_normalized <- dplyr::bind_cols(batch_data_list)
  
  return(as.data.frame(expression_log2_normalized))
}

expression_log2_normalized <- combine_batches(num_batches = 10)


############################  *Batch Median - After Norm* #########################

Norm_data = cbind(UniqueID, expression_log2_normalized)
Norm_data_exp = Norm_data[,-c(1)]
expression_normalized <- 2^Norm_data_exp
expression_normalized = cbind(UniqueID, expression_normalized)


melted_norm_data <- expression_normalized %>%
  pivot_longer(
    cols = -UniqueID,  # Keep UniqueID as an identifier, convert others to long format
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Batch = ifelse(grepl("Batch\\d+", variable, ignore.case = TRUE),
                   sub(".*(Batch\\d+).*", "\\1", variable),
                   NA)
  )

sample_medians_Afternorm <- melted_norm_data %>%
  group_by(Batch, variable) %>%
  summarize(sample_median = median(log2(value)))

# Calculate the median of sample medians for each batch
batch_medians_Afternorm <- sample_medians_Afternorm %>%
  group_by(Batch) %>%
  summarize(batch_median = median(sample_median))

##################### *Save Tables* ########################

write.csv(batch_medians, "D:/LongCOVID/Normalization/normalization_table.csv", row.names = FALSE)
write.csv(sample_medians_Afternorm, "D:/LongCOVID/Important Tables/MedianAbundanceTable_AllFeaturesFilteredNormalized.csv", row.names = FALSE)
write.csv(Norm_data, "D:/LongCOVID/Normalization/LongCOVID_Normalized_data_log2.csv", row.names = FALSE)
write.csv(expression_normalized, "D:/LongCOVID/Normalization/LongCOVID_Normalized_data_ActualValues.csv", row.names = FALSE)


##################### *Boxplot Before and After Norm - All Samples* ########################


melted_data <- melted_data %>%
  mutate(variable = factor(variable, levels = runorder$rawfile_name_R[order(runorder$RunOrder)])) 


p1 = ggplot(melted_data, aes(x = variable, y = log2(value), fill= Batch)) +
  geom_boxplot(width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("(Before Norm)") +
  ylim(1,40)+
  scale_fill_manual(values = batch_colors) +  # Use the custom color palette
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5, margin = margin(b = 1)),
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


melted_norm_data <- melted_norm_data %>%
  mutate(variable = factor(variable, levels = runorder$rawfile_name_R[order(runorder$RunOrder)])) 


p2 = ggplot(melted_norm_data, aes(x = variable, y = log2(value), fill= Batch)) +
  geom_boxplot(width = 0.7) +
  xlab("Samples") +
  ylab("Log2 Intensity") +
  ggtitle("(Normalized)") +
  ylim(1,40)+
  scale_fill_manual(values = batch_colors) +  # Use the custom color palette
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5, margin = margin(b = 1)),
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


pdf("D:/LongCOVID/Boxplot/Boxplot_Combined_Cohort_Batch_RunOrder_FilteredFeatures_Before and After Normalization.pdf", width = 40, height = 18)  # Adjust height to fit both plots
grid.arrange(p1, p2, ncol = 1)  
dev.off()




##################### *Boxplot Before and After Norm - per Batch* ########################

melted_data <- melted_data %>%
  mutate(Status = "Before Normalization")

melted_norm_data <- melted_norm_data %>%
  mutate(Status = "After Normalization")

melted_combined <- bind_rows(melted_data, melted_norm_data)

custom_batch_order <- c("Batch2", "Batch3", "Batch4", "Batch5", "Batch6", 
                        "Batch7", "Batch8", "Batch9", "Batch10", "Batch1")


melted_combined$Batch <- factor(melted_combined$Batch, levels = custom_batch_order)
melted_combined$Status <- factor(melted_combined$Status, levels = c("Before Normalization", "After Normalization"))



pdf("D:/LongCOVID/Boxplot/Boxplot_Combined_Cohort_Batch_RunOrder_FilteredFeatures_BatchesBefore and After Normalization.pdf", width = 15, height = 9)  
ggplot(melted_combined, aes(x = Batch, y = log2(value), fill = Status)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, width = 0.6) +  
  scale_fill_manual(values = c("Before Normalization" = "#caf0f8", "After Normalization" = "#415a77")) +
  theme_minimal() +  
  labs(
    title = "Batch Effect Before and After Normalization",
    x = "Batch",
    y = "Log2 (Abundance)",
    fill = "Normalization Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "#1d3557"),  # Centered, bold title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  # Rotated x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),  # Larger y-axis text
    axis.title.x = element_text(face = "bold", size = 14),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold y-axis title
    legend.position = "top",  # Move legend to the top
    legend.title = element_text(size = 13, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 12),  # Increase legend text size
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),  # Dashed grid lines for readability
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
dev.off()

