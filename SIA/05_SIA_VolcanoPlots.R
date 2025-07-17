
########################### load libraries ################################

library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggrepel)
library(gridExtra)


########################### Important Functions ################################

connect_db <- function(db_path) {
  dbConnect(RSQLite::SQLite(), dbname = db_path)
}

query_data <- function(con, table) {
  dbGetQuery(con, paste0("SELECT * FROM ", table))
}

omics_colors <- c("Protein" = "#1b9e77", "Lipid" = "#d95f02", "Transcript" = "#7570b3", "Unknown" = "beige")

lipid_colors <- c(
  "Fatty Acyls" = "#67C2A5", 
  "Sterol Lipids" = "#A5CF56", 
  "Glycerolipids" = "#F58B63", 
  "Sphingolipids" = "#E38BBB", 
  "Glycerophospholipids" = "#8B9FC9",
  "NA" = "beige"
)

############################ Pull data from DB ###################################

db_path <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite"

con <- connect_db(db_path)
biomolecules <- dbReadTable(con, "biomolecules")
biomolecules_metadata <- query_data(con, "biomolecules_metadata")
pvalues <- dbReadTable(con, "pvalues")
dbDisconnect(con)


############################ Edit data tables ###################################

# === Merge metadata ===

# Pivot biomolecules_metadata to wide format
lipid_metadata <- biomolecules_metadata %>%
  filter(metadata_type %in% c(
    "Lipid_Category", "Lipid_Class", "Main_Class", "Sub_Class", 
    "NumFattyAcylCarbons", "NumFattyAcylUnsaturations", "Unsaturation_Level"
  )) %>%
  pivot_wider(
    id_cols = biomolecule_id,
    names_from = metadata_type,
    values_from = metadata_value
  )


Volcano_data <- pvalues %>%
  filter(formula %in% c(35:38)) %>%
  inner_join(biomolecules %>%
               dplyr::select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  left_join(lipid_metadata, by = "biomolecule_id") %>%
  mutate(
    ome_for_color = case_when(
      omics_id == 2 & grepl("unknown", standardized_name, ignore.case = TRUE) ~ "Unknown",
      omics_id == 1 ~ "Protein",
      omics_id == 2 ~ "Lipid",
      omics_id == 3 ~ "Transcript",
      TRUE ~ "Other"
    ),
    point_alpha = ifelse(q_value < 0.05, 1, 0.2)
  )

############################## Plot Data by omic Type ########################################

plot_volcano_dynamic <- function(data, comparison_name, omics_colors) {
  df <- data %>%
    filter(comparison == comparison_name) 
  
  # Dynamic x-axis limits
  x_min <- min(df$effect_size, na.rm = TRUE)
  x_max <- max(df$effect_size, na.rm = TRUE)
  expand_range <- 0.2 * max(abs(c(x_min, x_max)))
  x_limits <- c(x_min - expand_range, x_max + expand_range)
  
  # Volcano plot
  ggplot(df, aes(
    x = effect_size,
    y = -log10(q_value),
    fill = ome_for_color,
    alpha = point_alpha
  )) +
    geom_point(shape = 21, colour = "grey30", size = 2.5, stroke = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
       geom_text_repel(
      data = subset(df, abs(effect_size) > 0.0005 & q_value < 0.05 & omics_id %in% c("1", "2", "3")),
      aes(label = standardized_name),
      size = 3,
      max.overlaps = 40,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = 'grey50'
    ) +
    scale_fill_manual(values = omics_colors, na.value = "lightgrey", name = "Biomolecule") +
    scale_alpha_identity() +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(0, max(-log10(df$q_value), na.rm = TRUE)),
                       breaks = scales::pretty_breaks(n = 6)) +
    labs(
      x = "Effect Size",
      y = expression(-log[10]("q-value")),
      title = paste(comparison_name)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
}

comparisons <- Volcano_data$comparison

comparisons <- c("QoL", "Age", "Sex", "BMI")

plot_list <- lapply(comparisons, function(comp) {
  plot_volcano_dynamic(Volcano_data, comp, omics_colors)
})

grid.arrange(grobs = plot_list, ncol = 4)


####### Plotting List ########################

plot_keys <- Volcano_data %>%
  dplyr::distinct(analysis_group, test) %>%
  dplyr::arrange(analysis_group, test)


plot_list <- list()

for (i in seq_len(nrow(plot_keys))) {
  ag <- plot_keys$analysis_group[i]
  tst <- plot_keys$test[i]
  
  sub_data <- Volcano_data %>%
    filter(analysis_group == ag, test == tst)
  
  comparisons <- unique(sub_data$comparison)
  
  for (comp in comparisons) {
    p <- plot_volcano_dynamic(sub_data, comp, omics_colors)
    
    plot_name <- paste0("AG", ag, "_", tst, "_", comp)
    plot_list[[plot_name]] <- p
  }
}


grid.arrange(grobs = plot_list[1:16], ncol = 2)



############## Save Plot ####################

output_dir <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/VolcanoPlots"

pdf(file.path(output_dir, "volcano_plots_combined_analysisGroup1_v2Unknown_noBatch1.pdf"), width = 12, height = 10)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()



################################################### Plot Data by Lipid Class #######################################

plot_volcano_lipids <- function(data, comparison_name, lipid_colors) {
  df <- data %>%
    filter(comparison == comparison_name)
  
  # Create separate layers
  df_lipids <- df %>%
    filter(omics_id == 2) %>%
    mutate(
      Lipid_Category = ifelse(q_value < 0.05,
                              ifelse(is.na(Lipid_Category), "NA", Lipid_Category),
                              "NotSignificant")  # override to gray
    )
  
  df_other <- df %>%
    filter(omics_id != 2)
  
  # Axis limits
  x_min <- min(df$effect_size, na.rm = TRUE)
  x_max <- max(df$effect_size, na.rm = TRUE)
  expand_range <- 0.2 * max(abs(c(x_min, x_max)))
  x_limits <- c(x_min - expand_range, x_max + expand_range)
  
  ggplot() +
    # Non-lipid background
    geom_point(
      data = df_other,
      aes(x = effect_size, y = -log10(q_value)),
      shape = 21, fill = "grey80", color = "grey30", size = 2.5, stroke = 0.3, alpha = 0.5
    ) +
    
    # Lipids on top
    geom_point(
      data = df_lipids,
      aes(x = effect_size, y = -log10(q_value), fill = Lipid_Category, alpha = point_alpha),
      shape = 21, color = "grey30", size = 2.5, stroke = 0.3
    ) +
    
    # Significance threshold
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    
    # Label only significant lipids
    geom_text_repel(
      data = df_lipids %>% filter(q_value < 0.05),
      aes(x = effect_size, y = -log10(q_value), label = standardized_name),
      size = 3, max.overlaps = 40,
      box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50'
    ) +
    
    # Fill color scale with "NotSignificant" gray
    scale_fill_manual(
      values = c(lipid_colors, NotSignificant = "grey80"),
      na.value = "beige",
      name = "Lipid Category"
    ) +
    scale_alpha_identity() +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(
      limits = c(0, max(-log10(df$q_value), na.rm = TRUE)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    
    labs(
      x = "Effect Size",
      y = expression(-log[10]("q-value")),
      title = paste(comparison_name)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
}

comparisons <- c("QoL", "Age", "Sex", "BMI")

comparisons <- c("PASC_noPASC", "Age", "Sex", "BMI")

plot_list_lipids <- lapply(comparisons, function(comp) {
  plot_volcano_lipids(Volcano_data, comp, lipid_colors)
})

gridExtra::grid.arrange(grobs = plot_list_lipids, ncol = 2)

pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/VolcanoPlots/VolcanoPlots_Lipids_analysisGroup1_noBatch1_PASC_noPASC.pdf", width = 14, height = 10)
gridExtra::grid.arrange(grobs = plot_list_lipids, ncol = 2)
dev.off()
