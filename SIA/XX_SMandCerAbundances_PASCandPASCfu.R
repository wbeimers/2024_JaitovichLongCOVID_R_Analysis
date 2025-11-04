
############################    *Load Libraries* ####################################

library(data.table)
library(readxl)
library(ggplot2)
library(tidyr)


########################### Important Functions ################################

connect_db <- function(db_path) {
  dbConnect(RSQLite::SQLite(), dbname = db_path)
}

query_data <- function(con, table) {
  dbGetQuery(con, paste0("SELECT * FROM ", table))
}


pal <- c(
  "Acute" = "#e09cc0",
  "Acute_fu" = "#a1537a",
  "Acute_NC" = "#b9babd",
  "Healthy" = "#488a56",
  "PASC" = "#89bed4",
  "PASC_fu" = "#667e9e"
)


############################    *Read Data* ####################################

db_path <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite"

con <- connect_db(db_path)

metadata <- query_data(con, "patient_metadata")
biomolecules <- query_data(con, "biomolecules")
biomolecules_metadata <- query_data(con, "biomolecules_metadata")
rawfiles <- query_data(con, "rawfiles_all")
formulas <- query_data(con, "formula_table")
pvalues <- query_data(con, "pvalues")

df_lipids <- query_data(con, "lipidomics_measurements")


dbDisconnect(con)


############################  *Sorting Data* ####################################

metadata <- metadata %>%
  mutate(
    SF.36.QOL.Score = ifelse(is.na(SF.36.QOL.Score) | SF.36.QOL.Score == "NULL", 900, as.numeric(SF.36.QOL.Score))
  )


df_lipids <- df_lipids %>%
  inner_join(rawfiles %>%
               dplyr::select(-timestamp, -run_type, -rawfile_name) %>%
               rename(keep_rawfiles = keep), by = "rawfile_id") %>%
  inner_join(metadata %>%
               dplyr::select(-Sample), by = "sample_id") %>%
  inner_join(biomolecules %>%
               dplyr::select(-standardized_name) %>%
               rename(keep_biomolecules = keep), by = "biomolecule_id") %>%
  filter(keep_rawfiles == "1", keep_biomolecules == "1",
         batch != 1)


# Filter for SM and Cer lipids
df_plot <- df_lipids %>%
  filter(grepl("^SM|^Cer", standardized_name))  # starts with SM or Cer


############################  *Plotting Data* ####################################

#Abundance of SMs and Cers Boxplots across patient groups
ggplot(df_plot, aes(x = Cohort, y = normalized_abundance, fill = Cohort)) +
  geom_boxplot() +
  facet_wrap(~standardized_name, scales = "free_y") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  labs(title = "Abundance of SMs and Ceramides Across Patient Groups", 
       y = "Normalized Abundance", x = "Patient Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Filtering only PASC and PASC_fu
df_plot_subset <- df_plot %>%
  filter(Cohort %in% c("PASC", "PASC_fu"))

#QoL relationship with SM and Cer abundance in PASc and PASC_fu
ggplot(df_plot_subset, aes(x = SF.36.QOL.Score, y = normalized_abundance, color = Cohort)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~standardized_name, scales = "free_y") +
  scale_color_manual(values = pal) +
  theme_bw() +
  labs(title = "Relationship Between Lipid Abundance and QoL Score (PASC and PASC Follow-up)",
       y = "Normalized Abundance", x = "SF-36 QoL Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

test = unique(df_plot_subset$standardized_name)

############################  *Plotting Averages Data* ####################################


#Averaged SM and Cer abundances

df_plot <- df_lipids %>%
  filter(grepl("^SM|^Cer", standardized_name)) %>%
  mutate(lipid_class = ifelse(grepl("^SM", standardized_name), "SM", "Cer"))

# average abundance for SM and Cer per sample
df_avg <- df_plot %>%
  group_by(Sample, Cohort, SF.36.QOL.Score, lipid_class) %>%
  summarize(avg_abundance = mean(normalized_abundance, na.rm = TRUE)) %>%
  ungroup()


ggplot(df_avg, aes(x = Cohort, y = avg_abundance, fill = Cohort)) +
  geom_boxplot() +
  facet_wrap(~lipid_class, scales = "free_y") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  labs(title = "Average Abundance of SMs and Ceramides Across Patient Groups", 
       y = "Average Normalized Abundance", x = "Patient Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 

# Filter PASC and PASC_fu only
df_avg_subset <- df_avg %>%
  filter(Cohort %in% c("PASC", "PASC_fu"))

ggplot(df_avg_subset, aes(x = SF.36.QOL.Score, y = avg_abundance, color = Cohort)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~lipid_class, scales = "free_y") +
  scale_color_manual(values = pal) +
  theme_bw() +
  labs(title = "Relationship Between Average Lipid Abundance and QoL Score (PASC and PASC Follow-up)",
       y = "Average Normalized Abundance", x = "SF-36 QoL Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##############  *SM/Cer ratio (not accurate due to more SMs identified than Ceramides)* ####################################


# Label SM and Cer lipids
df_plot <- df_lipids %>%
  filter(grepl("^SM|^Cer", standardized_name)) %>%
  mutate(lipid_class = ifelse(grepl("^SM", standardized_name), "SM", "Cer"))

# Average abundance per Sample and lipid class
df_ratio_sample <- df_plot %>%
  group_by(Sample, Cohort, SF.36.QOL.Score, Age, Sex, BMI, unique_patient_id, lipid_class) %>%
  summarize(total_abundance = mean(normalized_abundance, na.rm = TRUE)) %>%
  pivot_wider(names_from = lipid_class, values_from = total_abundance)

# SM/Cer ratio
df_ratio_sample <- df_ratio_sample %>%
  mutate(SM_Cer_Ratio = SM / Cer)

# Average per patient
df_ratio_patient <- df_ratio_sample %>%
  group_by(unique_patient_id) %>%
  summarize(
    mean_SM = mean(SM, na.rm = TRUE),
    mean_Cer = mean(Cer, na.rm = TRUE),
    mean_SM_Cer_Ratio = mean(SM_Cer_Ratio, na.rm = TRUE),
    Age = first(Age),
    Sex = first(Sex),
    BMI = first(BMI)
  )




############################    *Prepare Data* ####################################

# 1. Select only Ceramides
df_cer <- df_lipids %>%
  filter(grepl("^Cer", standardized_name))

# 2. Calculate **average** normalized abundance per sample
df_cer_avg <- df_cer %>%
  group_by(Sample, unique_patient_id, Cohort, SF.36.QOL.Score, Age, Sex, BMI) %>%
  summarize(avg_ceramide_abundance = mean(normalized_abundance, na.rm = TRUE)) %>%
  ungroup()

# 3. Focus only on PASC and PASC_fu for longitudinal analysis
df_cer_pasc <- df_cer_avg %>%
  filter(Cohort %in% c("PASC", "PASC_fu"))

############################    *Paired Line Plot (Longitudinal)* ####################################

ggplot(df_cer_pasc, aes(x = Cohort, y = avg_ceramide_abundance, group = unique_patient_id)) +
  geom_line(aes(color = unique_patient_id), alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "Paired Changes in Ceramide Abundance (PASC → PASC_fu)",
       x = "Cohort", y = "Average Ceramide Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################    *Paired Wilcoxon Test* ####################################

# Reshape data to wide format for paired testing
df_cer_wide <- df_cer_pasc %>%
  select(unique_patient_id, Cohort, avg_ceramide_abundance) %>%
  pivot_wider(names_from = Cohort, values_from = avg_ceramide_abundance)

# Perform paired Wilcoxon signed-rank test
wilcox_test_result <- wilcox.test(df_cer_wide$PASC, df_cer_wide$PASC_fu, paired = TRUE)
print(wilcox_test_result)

############################    *Link to Clinical Variables (Age, Sex, BMI)* ####################################

# Scatterplot: Ceramide vs Age
ggplot(df_cer_avg, aes(x = Age, y = avg_ceramide_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Ceramide Abundance vs Age", x = "Age", y = "Average Ceramide Abundance")

# Boxplot: Ceramide by Sex
ggplot(df_cer_avg, aes(x = Sex, y = avg_ceramide_abundance, fill = Sex)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Ceramide Abundance by Sex", x = "Sex", y = "Average Ceramide Abundance")

# Scatterplot: Ceramide vs BMI
ggplot(df_cer_avg, aes(x = BMI, y = avg_ceramide_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Ceramide Abundance vs BMI", x = "BMI", y = "Average Ceramide Abundance")








############################    *Prepare SM Data* ####################################

# 1. Select only Sphingomyelins
df_sm <- df_lipids %>%
  filter(grepl("^SM", standardized_name))

# 2. Calculate **average** normalized abundance per sample
df_sm_avg <- df_sm %>%
  group_by(Sample, unique_patient_id, Cohort, SF.36.QOL.Score, Age, Sex, BMI) %>%
  summarize(avg_sm_abundance = mean(normalized_abundance, na.rm = TRUE)) %>%
  ungroup()

# 3. Focus only on PASC and PASC_fu for longitudinal tracking
df_sm_pasc <- df_sm_avg %>%
  filter(Cohort %in% c("PASC", "PASC_fu"))

############################    *Paired Line Plot for SMs* ####################################

ggplot(df_sm_pasc, aes(x = Cohort, y = avg_sm_abundance, group = unique_patient_id)) +
  geom_line(aes(color = unique_patient_id), alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "Paired Changes in Sphingomyelin Abundance (PASC → PASC_fu)",
       x = "Cohort", y = "Average Sphingomyelin Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################    *Paired Wilcoxon Test for SMs* ####################################

# Reshape for paired test
df_sm_wide <- df_sm_pasc %>%
  select(unique_patient_id, Cohort, avg_sm_abundance) %>%
  pivot_wider(names_from = Cohort, values_from = avg_sm_abundance)

# Perform paired Wilcoxon signed-rank test
wilcox_test_result_sm <- wilcox.test(df_sm_wide$PASC, df_sm_wide$PASC_fu, paired = TRUE)
print(wilcox_test_result_sm)


############################    *Focus on PASC and PASC_fu Separately* ####################################

# Subset ceramide data for PASC only
df_cer_pasc_only <- df_cer_avg %>%
  filter(Cohort == "PASC")

# Subset ceramide data for PASC_fu only
df_cer_pascfu_only <- df_cer_avg %>%
  filter(Cohort == "PASC_fu")

# Subset sphingomyelin data for PASC only
df_sm_pasc_only <- df_sm_avg %>%
  filter(Cohort == "PASC")

# Subset sphingomyelin data for PASC_fu only
df_sm_pascfu_only <- df_sm_avg %>%
  filter(Cohort == "PASC_fu")

############################    *Scatterplots: Ceramide and SM vs QoL* ####################################

# Ceramide vs QoL - PASC
ggplot(df_cer_pasc_only, aes(x = SF.36.QOL.Score, y = avg_ceramide_abundance)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Ceramide Abundance vs QoL (PASC only)", x = "SF-36 QoL Score", y = "Average Ceramide Abundance")

# Ceramide vs QoL - PASC_fu
ggplot(df_cer_pascfu_only, aes(x = SF.36.QOL.Score, y = avg_ceramide_abundance)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Ceramide Abundance vs QoL (PASC follow-up only)", x = "SF-36 QoL Score", y = "Average Ceramide Abundance")

# Sphingomyelin vs QoL - PASC
ggplot(df_sm_pasc_only, aes(x = SF.36.QOL.Score, y = avg_sm_abundance)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Sphingomyelin Abundance vs QoL (PASC only)", x = "SF-36 QoL Score", y = "Average Sphingomyelin Abundance")

# Sphingomyelin vs QoL - PASC_fu
ggplot(df_sm_pascfu_only, aes(x = SF.36.QOL.Score, y = avg_sm_abundance)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Sphingomyelin Abundance vs QoL (PASC follow-up only)", x = "SF-36 QoL Score", y = "Average Sphingomyelin Abundance")

############################    *Correlation Calculations (Spearman)* ####################################

# Ceramide vs QoL correlation in PASC
cor.test(df_cer_pasc_only$avg_ceramide_abundance, df_cer_pasc_only$SF.36.QOL.Score, method = "spearman")

# Ceramide vs QoL correlation in PASC_fu
cor.test(df_cer_pascfu_only$avg_ceramide_abundance, df_cer_pascfu_only$SF.36.QOL.Score, method = "spearman")

# Sphingomyelin vs QoL correlation in PASC
cor.test(df_sm_pasc_only$avg_sm_abundance, df_sm_pasc_only$SF.36.QOL.Score, method = "spearman")

# Sphingomyelin vs QoL correlation in PASC_fu
cor.test(df_sm_pascfu_only$avg_sm_abundance, df_sm_pascfu_only$SF.36.QOL.Score, method = "spearman")




# Join Collection_date into df_cer_avg and df_sm_avg
df_cer_avg <- df_cer_avg %>%
  left_join(metadata %>% select(Sample, Collection_date), by = "Sample")

df_sm_avg <- df_sm_avg %>%
  left_join(metadata %>% select(Sample, Collection_date), by = "Sample")

# Ceramides
df_cer_pasc_only <- df_cer_avg %>% filter(Cohort == "PASC")
df_cer_pascfu_only <- df_cer_avg %>% filter(Cohort == "PASC_fu")

# Sphingomyelins
df_sm_pasc_only <- df_sm_avg %>% filter(Cohort == "PASC")
df_sm_pascfu_only <- df_sm_avg %>% filter(Cohort == "PASC_fu")


ggplot(df_cer_pascfu_only, aes(x = SF.36.QOL.Score, y = avg_ceramide_abundance)) +
  geom_point(size = 4, alpha = 0.9, shape = 16) +
  geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8, alpha = 0.4) +
  scale_color_gradient(
    low = "#aae5f9", high = "#08306b", name = "Collection Date",
    guide = guide_colorbar(barwidth = 0.6, barheight = 6)
  ) +
  scale_x_continuous(name = "SF-36 QoL Score") +
  scale_y_continuous(name = "Average Ceramide Abundance") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Ceramide Abundance vs QoL (PASC)")












# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1. Get IDs for paired samples (patients with both PASC and PASC_fu)
paired_ids <- df_cer_avg %>%
  filter(Cohort %in% c("PASC", "PASC_fu")) %>%
  count(unique_patient_id) %>%
  filter(n == 2) %>%
  pull(unique_patient_id)

# 2. Ceramide: Filter to paired samples
df_cer_paired <- df_cer_avg %>%
  filter(unique_patient_id %in% paired_ids, Cohort %in% c("PASC", "PASC_fu")) %>%
  select(unique_patient_id, Cohort, avg_ceramide_abundance, SF.36.QOL.Score)

# 3. Ceramide: Pivot to wide format and calculate deltas
df_cer_deltas <- df_cer_paired %>%
  pivot_wider(names_from = Cohort, values_from = c(avg_ceramide_abundance, SF.36.QOL.Score)) %>%
  mutate(
    delta_ceramide = avg_ceramide_abundance_PASC_fu - avg_ceramide_abundance_PASC,
    delta_qol = SF.36.QOL.Score_PASC_fu - SF.36.QOL.Score_PASC
  )

# 4. SM: Filter to paired samples
df_sm_paired <- df_sm_avg %>%
  filter(unique_patient_id %in% paired_ids, Cohort %in% c("PASC", "PASC_fu")) %>%
  select(unique_patient_id, Cohort, avg_sm_abundance, SF.36.QOL.Score)

# 5. SM: Pivot to wide format and calculate deltas
df_sm_deltas <- df_sm_paired %>%
  pivot_wider(names_from = Cohort, values_from = c(avg_sm_abundance, SF.36.QOL.Score)) %>%
  mutate(
    delta_sm = avg_sm_abundance_PASC_fu - avg_sm_abundance_PASC,
    delta_qol = SF.36.QOL.Score_PASC_fu - SF.36.QOL.Score_PASC
  )


plot_cer <- ggplot(df_cer_deltas, aes(x = delta_qol, y = delta_ceramide)) +
  geom_point(size = 2.5, color = "#a1537a") +  # muted purple points
  geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Δ QoL Score (Follow-up − Baseline)",
    y = "Δ Ceramide Abundance",
    title = "QoL Improvement vs Ceramide Reduction"
  )


# 7. Plot: QoL delta vs SM delta
plot_sm <- ggplot(df_sm_deltas, aes(x = delta_qol, y = delta_sm)) +
  geom_point(size = 2.5, color = "#667e9e") +
  geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Δ QoL Score (Follow-up − Baseline)",
    y = "Δ SM Abundance",
    title = "QoL Improvement vs SM Reduction"
  )

# 8. Combine plots
combined_plot <- plot_cer + plot_sm + plot_layout(ncol = 2)

# 9. View or save
print(combined_plot)
# ggsave("Delta_Lipids_vs_QoL.pdf", combined_plot, width = 10, height = 5)








library(ggplot2)
library(gridExtra)

# PASC Ceramide
plot1 <- ggplot(df_cer_pasc_only, aes(x = SF.36.QOL.Score, y = avg_ceramide_abundance)) +
  geom_point(size = 3, color = "#89bed4") +
  geom_smooth(method = "lm", color = "#89bed4") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Ceramide Abundance vs QoL (PASC only)", x = "SF-36 QoL Score", y = "Average Ceramide Abundance")

# PASC follow-up Ceramide
plot2 <- ggplot(df_cer_pascfu_only, aes(x = SF.36.QOL.Score, y = avg_ceramide_abundance)) +
  geom_point(size = 3, color = "#667e9e") +
  geom_smooth(method = "lm", color = "#667e9e") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Ceramide Abundance vs QoL (PASC follow-up only)", x = "SF-36 QoL Score", y = "Average Ceramide Abundance")

# PASC Sphingomyelin
plot3 <- ggplot(df_sm_pasc_only, aes(x = SF.36.QOL.Score, y = avg_sm_abundance)) +
  geom_point(size = 3, color = "#89bed4") +
  geom_smooth(method = "lm", color = "#89bed4") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Sphingomyelin Abundance vs QoL (PASC only)", x = "SF-36 QoL Score", y = "Average Sphingomyelin Abundance")

# PASC follow-up Sphingomyelin
plot4 <- ggplot(df_sm_pascfu_only, aes(x = SF.36.QOL.Score, y = avg_sm_abundance)) +
  geom_point(size = 3, color = "#667e9e") +
  geom_smooth(method = "lm", color = "#667e9e") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Sphingomyelin Abundance vs QoL (PASC follow-up only)", x = "SF-36 QoL Score", y = "Average Sphingomyelin Abundance")

# Save all plots in a single PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/SM_Cer/Ceramide_Sphingomyelin_vs_QoL.pdf", width = 18, height = 7)
grid.arrange(plot1, plot2, plot3, plot4, ncol = 4)
dev.off()
