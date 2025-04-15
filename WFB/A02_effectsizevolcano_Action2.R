#### Overview ####

# 1. Volcano Plots of Effect Size after filtering
# 2. fgsea of the ranked effect size lists


#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(viridis)
library(RSQLite)
library(ggrepel)
library(data.table)
library(fgsea)
library(gridExtra)
library(patchwork)
library(colorspace)


# Colors #
col <- brewer.pal(8, "Dark2") 

pal1 <- c("#66C2A5",
          "#FFD92F",
          "#8DA0CB",
          "#FC8D62",
          "#A6CEE3",
          "#E78AC3",
          "#A6D854",
          "#FDB462",
          "#B3B3B3",
          "#B2DF8A")

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')


pal <- c("Acute" = "#E78AC3", 
         "Acute_fu" = '#AA3377', 
         "Acute_NC" = "#B3B3B3", 
         "Healthy" = '#229100', 
         "PASC" = '#66CCEE', 
         "PASC_fu" = '#4477AA')



# plot colors
pie(rep(1, length(col)), col = col , main="") 


colorRampPalette(col)(15)


# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

pvalues <- dbGetQuery(con, 'SELECT biomolecule_id, analysis_group, test, comparison, formula, predictor, effect_size, eta_squared, lratio, p_value, q_value
                            FROM pvalues')
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
metadata <- dbGetQuery(con, 'SELECT sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, unique_patient_id, Collection_date, PG_change_collection_cutoff
                             FROM patient_metadata')
formulas <- dbGetQuery(con, 'SELECT *
                             FROM formula_table')

dbDisconnect(con)



## Volcano Plot Preparation ----
# options:
# analysis_group (1, 2, 3, 0)
anal <- 1
# comparison (Age, Sex, QoL, BMI)
comp <- "PASC_noPASC"
# formula (1, 2, 3, etc.)
form <- 35
# ome
omea <- "protein"


volc_plot <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 ~ "YES",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id)

counts <- volc_plot %>%
  filter(diffexp %in% c("YES", "NO")) %>%
  count(diffexp)

# Color Based on GO TERM
GO_colors <- c(
  #"GO:0000786", "GO:0005884" # trascript
  "GO:0060589", "GO:0035173", "GO:0006959", "GO:0030545" # protein
               )

con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                                          FROM biomolecules_metadata')
GO_terms <- dbGetQuery(con, 'SELECT *
                             FROM GO_terms')
dbDisconnect(con)

biomolecule_metadata_GO <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  separate_rows(metadata_value, sep = ";") %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id), .groups = "drop") %>%
  filter(metadata_value %in% GO_colors)

volc_plot <- volc_plot %>%
  left_join(biomolecule_metadata_GO, by = c("biomolecule_id" = "GO_terms")) %>% 
  mutate(metadata_value = ifelse(is.na(metadata_value), "NO", metadata_value))


ggplot(volc_plot %>%
         filter(metadata_value == "NO") %>%
         filter(ome != "protein"), 
       aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.1,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value == "NO") %>%
               filter(ome == "protein"), 
             mapping = aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.6,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value != "NO") %>%
               filter(ome == "protein"), 
             mapping = aes(effect_size, neglogpvalue,
                           fill = metadata_value,
                           size = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2,
             alpha = 1) +
  scale_fill_manual(values = c("#1B9E77", "#57C39C", "#83EAC2", "#FAFFFD")) +
  scale_size_manual(values = c(0.5, 1.5), guide = "none") +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_x_continuous(limits = c(-max(abs(volc_plot$effect_size)), max(abs(volc_plot$effect_size)))) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab(paste("Effect Size", comp)) +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 
#+
#  geom_text(data = counts[counts$diffexp == "YES",], 
#            aes(x = -0.001, y = Inf, 
#                label = paste(n)),
#            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
#  geom_text(data = counts[counts$diffexp == "NO",], 
#            aes(x = 0.001, y = Inf, 
#                label = paste(n)),
#            hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) 
ggsave(paste0("reports/figures/Volcano_group_", anal, "_", comp, "_formula", form, "_ProteinGOterms_col.pdf"), 
       width = 8, height = 6, units = "cm")



## Enrichment Analysis ----
GO_set <- fread("data/metadata/GOtermset_PASCnoPASC_Transcripts_5to100_90coverage_10overlap.csv")

# set up GO term list
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                                          FROM biomolecules_metadata')
GO_terms <- dbGetQuery(con, 'SELECT *
                             FROM GO_terms')
dbDisconnect(con)

biomolecule_metadata_GO <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  separate_rows(metadata_value, sep = ";") %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id), .groups = "drop") 
#%>% filter(metadata_value %in% GO_set$GO_term)

GO_term_list <- split(biomolecule_metadata_GO$GO_terms,
                      biomolecule_metadata_GO$metadata_value)

# Rank of biomolecules by effect size and pvalue
# only do proteins or transcripts one at a time
gene_list <- volc_plot %>%
  filter(ome == omea)

gene_list <- gene_list %>%
  mutate(rank = effect_size * neglogpvalue)

gene_list <- setNames(gene_list$rank,
                      gene_list$biomolecule_id)


# do enrichment
fgsea <- fgsea(pathways = GO_term_list, 
               stats    = gene_list,
               minSize  = 5,
               maxSize  = 1000) %>%
  left_join(GO_terms %>%
              select(GO_term, name),
            by = c("pathway" = "GO_term"))

fwrite(fgsea, paste0("data/processed/fgsea_AnalysisGroup1_PASCnoPASC_", omea, "_limitedGOterms.csv"))


##plot:GSEA ----
# List of pathways you want to plot
selected_pathways <- c("GO:0005788", "GO:0030027", "GO:0031012", "GO:0035173")


#plots <- lapply(selected_pathways, function(pathway) {
#  enrichment_plot <- plotEnrichment(GO_term_list[[pathway]], gene_list) +
#    ggtitle(pathway) +
#    theme_minimal() +
#    theme(plot.title = element_text(size = 12, hjust = 0.5))
#  return(enrichment_plot)
#})
#
# Arrange the plots in a grid (or combine them)
#grid.arrange(grobs = plots, ncol = 1)
#ggsave(paste0("reports/figures/FGSEA_group_top4_", anal, "_", comp, "_formula", form, ".pdf"), 
#       width = 8, height = 6, units = "cm")









# Generate enrichment data for each pathway
pd_list <- setNames(lapply(selected_pathways, function(pathway) {
  plotEnrichmentData(pathway = GO_term_list[[pathway]], stats = gene_list)
}), selected_pathways)

combined_data_curve <- do.call(rbind, lapply(names(pd_list), function(pathway) {
  pd <- pd_list[[pathway]]
  curve_data <- pd$curve  
  curve_data$pathway <- pathway  
  return(curve_data)  
}))

combined_data_ticks <- do.call(rbind, lapply(names(pd_list), function(pathway) {
  pd <- pd_list[[pathway]]
  ticks_data <- pd$ticks  
  ticks_data$pathway <- pathway  
  return(ticks_data)  
}))


# Now you can plot the combined data using ggplot
curve <- ggplot(combined_data_curve, aes(rank, ES, color = pathway)) +
  geom_line(size = 0.4, 
            alpha = 0.8) + 
  geom_hline(yintercept = 0, 
             color = "black",
             size = 0.2) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "gray", linetype = "solid", size = 0.2),
    legend.position = "none",
    axis.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  labs(x = NULL, y = "Enrichment Score")

ticks <- ggplot(combined_data_ticks, aes(color = pathway)) +
  geom_segment(aes(x = rank, y = -0.6/16,
                  xend = rank, yend = 0.6/16),
               size = 0.2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, max(combined_data_ticks$rank), by = 1000)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    panel.spacing = unit(0, "lines"),
    strip.text = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.25, "cm")
  ) +
  labs(x = "Rank", 
       y = NULL) +
  facet_grid(pathway ~ .)

curve / ticks + plot_layout(heights = c(2, 1))

ggsave(paste0("reports/figures/FGSEA_group_top4_", omea, "_", anal, "_", comp, "_formula", form, ".pdf"), 
       width = 16, height = 8, units = "cm")

goterms <- combined_data_curve %>%
  select(pathway) %>%
  distinct() %>%
  left_join(GO_terms, by = c("pathway" = "GO_term"))













#### Analysis Group 1 Overview ####
# Find the overlap of the ML binary/continuous and the LM binary/continuous
# Take the top 1000 features of each analysis and look at the overlap

# ML binary
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

ML_coef <- dbGetQuery(con, 'SELECT *
                            FROM machine_learning_coefficients')

dbDisconnect(con)


# ML binary
coef_bin <- ML_coef %>%
  filter(formula_id == 44) %>%
  filter(mean != 0) %>%
  arrange(desc(`percent nonzero`)) %>%
  mutate(rank = row_number())


# LM binary
lm_bin_top <- pvalues %>%
  filter(analysis_group == 1) %>%
  filter(comparison == "PASC_noPASC") %>%
  filter(formula == 35) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  arrange(p_value) %>%
  mutate(rank = row_number()) %>%
  filter(q_value < 0.05) 


# ML continuous
coef_cont <- ML_coef %>%
  filter(formula_id == 43) %>%
  filter(mean != 0) %>%
  arrange(desc(`percent nonzero`)) %>%
  mutate(rank = row_number()) %>%
  filter(rank >= 1 & rank <= 1000) 


# LM cont

## Combine and look at overlap
library(eulerr)
overlaps <- list(
  #ML_cont = coef_cont$biomolecule_id,
  ML_bin = coef_bin$biomolecule_id,
  #LM_cont = lm_cont_top1000$biomolecule_id,
  LM_bin = lm_bin_top$biomolecule_id
)

euler_ov <- euler(overlaps)

pdf("reports/figures/EulerPlot_AnalysisGroup1_Overlap_top40.pdf", width = 2, height = 2)
plot(euler_ov, 
     fills = c(col[4], col[5], col[6]), 
     edges = T,
     quantities = T)
dev.off()



## filter for only overlapping ones, and combine into one long dataframe
common <- Reduce(intersect, overlaps)

common_df <- coef_bin %>%
  filter(biomolecule_id %in% common) %>%
  select(biomolecule_id, rank) %>%
  rename(rank_ML_bin = rank) %>%
  mutate(biomolecule_id = as.integer(biomolecule_id)) %>%
  left_join(lm_bin_top %>%
              select(biomolecule_id, rank) %>%
              rename(rank_LM_bin = rank),
            by = "biomolecule_id")

common_df_long <- common_df %>%
  pivot_longer(cols = c(-biomolecule_id),
               names_to = "rank_type",
               values_to = "rank")

## plot:overlapping feature ranks ----
ggplot(common_df_long, aes(x = rank, y = rank_type, color = as.factor(omics_id))) +
  geom_point(size = 1.5,
             alpha = 0.8) +
  scale_color_manual(values = c(col[1], col[3])) +
  labs(x = "Rank", y = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave("reports/figures/AnalysisGroup1_Method_Overlap.pdf", 
       width = 8, height = 4, units = "cm")



## Make a composite rank of the ranks from each source, to find the overall most important features

common_df <- common_df %>%
  mutate(composite_rank = rank_ML_cont * rank_LM_cont * rank_LM_bin) %>%
  arrange(composite_rank) %>%
  mutate(final_rank = rank(composite_rank, ties.method = "min")) %>%
  select(-composite_rank) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name),
             by = "biomolecule_id")

## Try the RRHO script
library(RRHO)

asd <- RRHO(
  list1 = common_df %>%
    select(biomolecule_id, rank_ML_cont), 
  list2 = common_df %>%
    select(biomolecule_id, rank_LM_cont),
  labels = c("ML", "LM"),
  alternative = "enrichment",
  plots = T,
  outputdir = "reports/figures",
  BY = FALSE,
  log10.ind = FALSE)



writeLines(as.character(common_df$biomolecule_id), "data/processed/AnalysisGroup1_Top1000_BiomoleculeOverlap.txt")


## Overlapping feature single protein boxplots ----
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
proteomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                               FROM proteomics_measurement')
lipidomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                              FROM lipidomics_measurements')
transcriptomics <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, sample_id, Counts, normalized_counts
                                    FROM rnaseq_measurements")
rawfiles <- dbGetQuery(con, 'SELECT rawfile_id, rawfile_name, Sample, sample_id, run_type, ome_id, keep
                             FROM rawfiles_all')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                             FROM biomolecules_metadata')
dbDisconnect(con)


## Make proteomics/lipidomics/transcriptomics dataframes ----

# proteomics
# Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles_p <- rawfiles %>%
  filter(ome_id == 1) %>%
  select(-keep) %>%
  filter(grepl('Sample', run_type))

df_p <- proteomics %>%
  inner_join(rawfiles_p, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_p <- biomolecules %>%
  filter(omics_id == 1) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_p <- df_p %>%
  filter(biomolecule_id %in% biomolecules_p)

# THINK ABOUT FILTERING OUT POST-TUBE CHANGE SAMPLES


# lipidomics
rawfiles_l <- rawfiles %>%
  filter(ome_id == 2) %>%
  filter(keep == "1") %>%
  filter(run_type == "Sample") %>%
  select(-keep)

df_l <- lipidomics %>%
  inner_join(rawfiles_l, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_l <- biomolecules %>%
  filter(omics_id == 2) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_l <- df_l %>%
  filter(biomolecule_id %in% biomolecules_l)


# transcriptomics
df_t <- transcriptomics %>%
  filter(!is.na(sample_id)) %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_t <- biomolecules %>%
  filter(omics_id == 3) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_t <- df_t %>%
  filter(biomolecule_id %in% biomolecules_t) %>%
  rename(normalized_abundance = normalized_counts)


# combined
df_a <- df_p %>%
  mutate(ome = "p") %>%
  bind_rows(df_l %>%
              mutate(ome = "l")) %>%
  select(-rawfile_id, -rawfile_name, -run_type, -ome_id) %>%
  bind_rows(df_t %>%
              rename(raw_abundance = Counts,
                     normalized_abundance = normalized_counts) %>%
              mutate(ome = "t"))

filtered_df_a <- filtered_df_p %>%
  mutate(ome = "p") %>%
  bind_rows(filtered_df_l %>%
              mutate(ome = "l")) %>%
  select(-rawfile_id, -rawfile_name, -run_type, -ome_id) %>%
  bind_rows(filtered_df_t %>%
              rename(raw_abundance = Counts) %>%
              mutate(ome = "t"))

## Filter for analysis_group_1 samples, and overlapping biomolecules
## plot:boxplots - shared biomolecules ----
comparison_df_a <- filtered_df_a %>%
  filter(analysis_group_1 == 1) %>%
  filter(biomolecule_id %in% common_df$biomolecule_id)

for (i in unique(comparison_df_a$standardized_name)) {
  
  poi <- i
  
  single_feat_df <- comparison_df_a %>%
    filter(standardized_name == poi)
  
  om <- unique(single_feat_df$ome)
  
  ggplot(single_feat_df, aes(Cohort, normalized_abundance, fill = Cohort, color = Cohort)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = NULL,
         y = "Log2 Abundance") +
    scale_y_continuous(expand = c(0,0), limits = c(min(single_feat_df$normalized_abundance) / 1.1, 
                                                   max(single_feat_df$normalized_abundance) * 1.1)) +
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
  
  ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AnalysisGroup1_PreTube_singleprotein_', poi, '_distribution_Cohort.pdf'), 
         width = 8, height = 6, units = "cm")
}

## plot:pointplots - QoL shared biomolecules ----
comparison_df_a_1 <- comparison_df_a %>%
  mutate(SF.36.QOL.Score = ifelse(is.na(SF.36.QOL.Score), 900, SF.36.QOL.Score))


for (i in unique(comparison_df_a_1$standardized_name)) {
  
  poi <- i
  
  single_feat_df <- comparison_df_a_1 %>%
    filter(standardized_name == poi)
  
  om <- unique(single_feat_df$ome)
  
  ggplot(single_feat_df, aes(normalized_abundance, SF.36.QOL.Score, fill = Cohort)) + 
    geom_smooth(method = "lm", 
                se = TRUE, 
                color = "gray", 
                fill = "lightgray",
                size = 0.5) +
    geom_point(alpha = 0.9, 
               shape = 21,
               stroke = 0.2,
               color = "black",
               size = 1.5) +
    scale_fill_manual(values = pal) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = "Log2 Abundance",
         y = "QOL") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 950)) +
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
  
  ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AnalysisGroup1_PreTube_singleprotein_', poi, '_scatter_Cohort.pdf'), 
         width = 8, height = 6, units = "cm")
}


## Check for overlapping Gene Names ----

common_df_gene <- common_df %>%
  left_join(biomolecules_metadata %>%
              filter(metadata_type %in% c("Entry_name", "gene_symbol")) %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id")


## Enrichment of Overlap Heatmap Clusters ----
# set up GO term list
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                                          FROM biomolecules_metadata')
GO_terms <- dbGetQuery(con, 'SELECT *
                             FROM GO_terms')
dbDisconnect(con)

biomolecule_metadata_GO <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  separate_rows(metadata_value, sep = ";") %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id), .groups = "drop")

GO_term_list <- split(biomolecule_metadata_GO$GO_terms,
                      biomolecule_metadata_GO$metadata_value)

# get gene list and background list
clusters <- fread("reports/figures/Heatmap_plots/Clusters_AllOmes_AnalysisGroup1_4methodoverlap_2clusters.csv")

####### switch between clusters here!!!
background_list <- filtered_df_a %>%
  filter(ome != "l") %>%
  select(biomolecule_id, ome) %>%
  distinct() %>%
  pull(biomolecule_id)


## Katie Function for testing enrichment ##
enrichment <- function(set, reference_sets, background){
  # output p_value
  # output fdr
  nset <- length(set)
  set <- as.character(set)
  nbackground <- length(background)
  background <- as.character(background)  
  output <- data.frame(reference = names(reference_sets), 
                       pvalue = rep(1, length(names(reference_sets))), 
                       fdr_pvalue = rep(NA, length(names(reference_sets))),
                       stringsAsFactors = F)
  for (i in 1:nrow(output)){
    hits <- length(intersect(set, reference_sets[[output[i,1]]]))
    hitsBackground <- length(intersect(background, reference_sets[[output[i,1]]]))
    if (hits > 0){
      output$pvalue[i] <- phyper(hits-1, hitsBackground, length(background) - hitsBackground, length(set), lower.tail = F)
    } 
    if (length(reference_sets[[output[i,1]]]) == 1){
      output$pvalue[i] <- NA
    }
  }
  output$fdr_pvalue <- p.adjust(output$pvalue, method = "BH")
  output
} 

gene_list <- common_df_gene %>%
  left_join(clusters, by = "biomolecule_id") %>%
  filter(cluster_assignments == 1) %>%
  pull(biomolecule_id)

enrichment_df_clust1 <- enrichment(gene_list, GO_term_list, background_list) %>%
  left_join(GO_terms, by = c("reference" = "GO_term"))

gene_list <- common_df_gene %>%
  left_join(clusters, by = "biomolecule_id") %>%
  filter(cluster_assignments == 2) %>%
  pull(biomolecule_id)

enrichment_df_clust2 <- enrichment(gene_list, GO_term_list, background_list) %>%
  left_join(GO_terms, by = c("reference" = "GO_term"))



## Plot enrichment barplots
enrichment_df_clust1_plot <- enrichment_df_clust1 %>%
  mutate(neglog10 = -log10(fdr_pvalue)) %>%
  filter(neglog10 >= 1.301) %>%
  mutate(cluster = 1)

enrichment_df_clust2_plot <- enrichment_df_clust2 %>%
  mutate(neglog10 = -log10(fdr_pvalue)) %>%
  filter(neglog10 >= 1.301) %>%
  mutate(cluster = 2)

enrichment_df_plot_combined <- enrichment_df_clust1_plot %>%
  bind_rows(enrichment_df_clust2_plot) %>%
  arrange(desc(neglog10))

ggplot(enrichment_df_plot_combined, aes(x = neglog10, y = name, fill = as.factor(cluster))) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y") +
  labs(x = "-log10(FDR pvalue)", 
       y = NULL) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave("reports/figures/AnalysisGroup1_4methodoverlap_2clusters_pvalues.pdf", 
       width = 16, height = 8, units = "cm")


