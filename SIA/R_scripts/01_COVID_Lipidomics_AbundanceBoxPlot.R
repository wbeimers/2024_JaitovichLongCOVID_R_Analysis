
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


############################  *Read Data* ####################################

data_samples <- read.csv("D:/LongCOVID/Important Tables/DataSubsets/20250207_longCOVID_AllSamples_RunOrder.csv", header= TRUE)
runorder <- read.csv("D:/LongCOVID/Metadata/rawfiles_Metadata.csv", header= TRUE)


############################   *Color Scheme* ####################################

batch_colors <- c("Batch1" = "#FFDDAA", "Batch2" = "#C3E4A6", "Batch3" = "#00CC99",
                  "Batch4" = "#9DAE9F", "Batch5" = "#DDDDCC", "Batch6" = "#A2A9B9",
                  "Batch7" = "#99CCDF", "Batch8" = "#0099CC", "Batch9" = "#415A87",
                  "Batch10" = "#FFDDEE")

batch_order <- paste0("Batch", 1:10)


cohort_colors <- setNames(brewer.pal(length(unique(melted_boxplot_data$Cohort)), "Set2"), 
                          unique(melted_boxplot_data$Cohort))
cohort_colors <- c(
  "Acute covid" = "#EE6677",
  "Acute f/u (no PASC)" = "#AA3377",
  "Acute non-COVID" = "#CCBB44",
  "Healthy" = "#228833",
  "PASC" = "#66CCEE",
  "PASC f/u" = "#4477AA"
)

############################   *Expression Data* ####################################

UniqueID = data.frame(UniqueID = data_samples$UniqueID)

sample_filenames <- runorder %>%
  filter(runorder$run_type == "Sample") %>%
  pull(rawfile_name_R)  

Expression <- data_samples%>%
  select(any_of(sample_filenames))

Expression_log2 <- log2(Expression)

Expression_log2_matrix= as.matrix(Expression_log2)

ordered_filenames <- runorder %>%
  filter(runorder$run_type == "Sample") %>%
  arrange(RunOrder) %>%
  pull(rawfile_name_R) 

expression_reordered <- Expression_log2 %>%
  select(one_of(ordered_filenames))


############################  *Boxplot* ####################################

Boxplot <- cbind(UniqueID, Expression_log2)

melted_boxplot_data <- Boxplot %>%
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


melted_boxplot_data <- melted_boxplot_data %>%
  left_join(runorder %>% select(rawfile_name_R, RunOrder, Cohort),  
            by = c("variable" = "rawfile_name_R")) 

melted_boxplot_data <- melted_boxplot_data %>%
  mutate(variable = factor(variable, levels = runorder$rawfile_name_R[order(runorder$RunOrder)])) 



######### Cohort #################

pdf("D:/LongCOVID/Boxplot/Boxplot_Cohort_RunOrder_FilteredFeatures.pdf", width = 20, height = 9)
p1=ggplot(melted_boxplot_data, aes(x = variable, y = value, fill = Cohort)) +
  geom_boxplot() +  # Hide outliers for clarity
  scale_fill_manual(values = cohort_colors) +  
  theme_minimal() +
  labs(title = "LongCOVID - Abundance", x = "Samples", y = "Log2(Relative Abundance)", fill = "Cohort") +
  theme(
    legend.position = "right"
  )
dev.off()

######### Batch #################

pdf("D:/LongCOVID/Boxplot/Boxplot_Batch_RunOrder_FilteredFeatures.pdf", width = 20, height = 9)
p2=ggplot(melted_boxplot_data, aes(x = variable, y = value, fill = Batch)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers for clarity
  scale_fill_manual(values = batch_colors) +  # Assign Set2 colors to batches
  theme_minimal() +
  labs(title = "LongCOVID - Abundance", x = "Samples", y = "Value", fill = "Batch") +
  theme(
    legend.position = "right",
  )
dev.off()

############# Combined Boxplot ##################

pdf("D:/LongCOVID/Boxplot/Boxplot_Combined_Cohort_Batch_RunOrder_FilteredFeatures.pdf", width = 20, height = 18)  # Adjust height to fit both plots
grid.arrange(p1, p2, ncol = 1)  
dev.off()

#-------------------------------------------------