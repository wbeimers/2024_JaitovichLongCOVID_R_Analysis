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

col1 <- viridis_pal(option = "mako")(100)[round(c(0.25, 0.5, 0.75) * 100)]


# plot colors
pie(rep(1, length(col)), col = col , main="") 

color_to_adjust <- col[2]
lighter_shades <- lighten(color_to_adjust, amount = c(0, 0.4))


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
anal <- 7
# comparison (Age, Sex, QoL, BMI)
comp <- "group7_PASCnoPASC"
# formula (1, 2, 3, etc.)
form <- 57
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
  "GO:0030527", "GO:0015629" # transcript
  #"GO:0060589", "GO:0006959" # protein
               )

table(volc_plot$ome[volc_plot$diffexp == "YES"])

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


tp <- ggplot(volc_plot %>%
         filter(metadata_value == "NO") %>%
         filter(ome != omea), 
       aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.05,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value == "NO") %>%
               filter(ome == omea), 
             mapping = aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.5,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value != "NO") %>%
               filter(ome == omea), 
             mapping = aes(effect_size, neglogpvalue,
                           fill = metadata_value,
                           size = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2,
             alpha = 1) +
  scale_fill_manual(values = lighter_shades) +
  scale_size_manual(values = c(0.5, 2), guide = "none") +
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
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
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

pp + lp + tp

ggsave(paste0("reports/figures/Volcano_group_", anal, "_", comp, "_formula", form, "_ALLLLGOterms_col.pdf"), 
       width = 7, height = 2, units = "in")



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

fwrite(fgsea, paste0("data/processed/fgsea_AnalysisGroup7_PASCnoPASC_", omea, "_allGOterms.csv"))


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
              mutate(ome = "t")) %>%
  filter(Cohort %in% c("Healthy", "Acute_fu", "PASC")) %>%
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID",
    Cohort %in% c("PASC") ~ "Long_COVID")) 

## Filter for analysis_group_1 samples, and overlapping biomolecules
## plot:boxplots - significant biomolecules PASCnoPASC ----
signif_feat <- volc_plot %>%
  filter(diffexp == "YES") %>%
  pull(biomolecule_id)

comparison_df_a <- filtered_df_a %>%
  filter(analysis_group_1 == 1) %>%
  filter(biomolecule_id %in% signif_feat)

for (i in unique(comparison_df_a$standardized_name)) {
  
  poi <- i
  
  single_feat_df <- comparison_df_a %>%
    filter(standardized_name == poi)
  
  om <- unique(single_feat_df$ome)
  
  ggplot(single_feat_df, aes(PASCnoPASC, normalized_abundance, fill = PASCnoPASC, color = PASCnoPASC)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_manual(values = col1) +
    scale_color_manual(values = col1) +
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
  
  ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AnalysisGroup7_PreTube_singlebiomolecule_', poi, '_distribution_PASCnoPASC.pdf'), 
         width = 8, height = 6, units = "cm")
}

ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AnalysisGroup1_PreTube_singleprotein_POTEJ_distribution_Cohort.pdf'), 
       width = 6, height = 4, units = "cm")


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








## Dial into specific pathways that are enriched ----
# find a pathway, find proteins, and plot log2fc for both PASCnoPASC and Acute/Healthy

pathway <- "GO:0060589"

pathway_bms <- GO_term_list[[pathway]]

bmol_genes <- biomolecules_metadata %>%
  select(-metadata_id) %>%
  filter(metadata_type %in% c("gene_name", "gene_symbol")) %>%
  group_by(biomolecule_id) %>%
  filter(n() == 1 | metadata_type == "gene_symbol") %>%
  ungroup() %>%
  select(biomolecule_id, metadata_value)

# organize dataframes
fc_barplots <- volc_plot %>%
  filter(biomolecule_id %in% pathway_bms) %>%
  filter(q_value < 0.05) %>%
  left_join(bmol_genes,
            by = "biomolecule_id")

ggplot(fc_barplots, aes(reorder(metadata_value, q_value), effect_size, fill = ome)) + 
  geom_col(position = position_dodge(),
           width = 0.6) + 
  geom_hline(yintercept = 0,
             size = 0.2) +
  scale_fill_manual(values = c(col[1], col[3])) +
  labs(x = NULL,
       y = "Effect Size",
       title = paste(pathway)) +
  scale_y_continuous(expand = c(0,0), limits = c(min(fc_barplots$effect_size) * 1.1, 
                                                 max(fc_barplots$effect_size) * 1.1)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Biomolecule')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.9, 0.9), 
        legend.justification = c("right", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")
  ) 
ggsave(paste0('reports/figures/PASCnoPASC_Biomolecules_All_GO0015629_effectsizes.pdf'), 
       width = 16, height = 6, units = 'cm')



ggplot(fc_barplots, aes(-log10(q_value), effect_size, color = ome)) + 
  geom_point() + 
  geom_text_repel(aes(label = metadata_value.y),
                  size = 2) +
  scale_color_manual(values = c(col[1], col[3])) +
  labs(y = "Effect Size",
       x = "-log10(q value)",
       title = paste(pathway)) +
  scale_y_continuous(expand = c(0,0), limits = c(min(fc_barplots$effect_size) * 1.1, 
                                                 max(fc_barplots$effect_size) * 1.1)) +
  scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Biomolecule')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.9, 0.1), 
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "cm")
  ) 
ggsave(paste0('reports/figures/PASCnoPASC_Biomolecules_All_GO0015629_effectsizes.pdf'), 
       width = 12, height = 6, units = 'cm')






## Lipid Enrichment Analysis ----
# options:
# analysis_group (1, 2, 3, 0)
anal <- 7
# comparison (Age, Sex, QoL, BMI)
comp <- "group7_PASCnoPASC"
# formula (1, 2, 3, etc.)
form <- 57
# ome
omea <- "lipid"


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


# choose a lipid metadata to use to filter ("Unsaturation_Level", "Lipid_Category", "Lipid_Class", 
# "Main_Class", "NumFattyAcylCarbons", "NumFattyAcylUnsaturations", "Sub_Class")
opt <- "Main_Class"
# set up Lipid term list
biomolecules_metadata_lipid <- biomolecules_metadata %>%
  filter(metadata_type == opt) %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id)) 
#%>% filter(metadata_value %in% GO_set$GO_term)

lipid_term_list <- split(biomolecules_metadata_lipid$GO_terms,
                      biomolecules_metadata_lipid$metadata_value)

# Rank of biomolecules by effect size and pvalue
# only do proteins or transcripts one at a time
gene_list <- volc_plot %>%
  filter(ome == omea)

gene_list <- gene_list %>%
  mutate(rank = effect_size * neglogpvalue)

gene_list <- setNames(gene_list$rank,
                      gene_list$biomolecule_id)


# do enrichment
fgsea <- fgsea(pathways = lipid_term_list, 
               stats    = gene_list,
               minSize  = 5,
               maxSize  = 1000) 

fwrite(fgsea, paste0("data/processed/fgsea_AnalysisGroup7_PASCnoPASC_", omea, "_", opt, ".csv"))



# color by lipids stuff
opti <- c("Glycerophosphocholines", "Phosphosphingolipids") # lipid

biomolecules_metadata_lipid <- biomolecules_metadata %>%
  filter(metadata_type == opt) %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id)) %>%
  filter(metadata_value %in% opti)


volc_plot <- volc_plot %>%
  left_join(biomolecules_metadata_lipid, by = c("biomolecule_id" = "GO_terms")) %>% 
  mutate(metadata_value = ifelse(is.na(metadata_value), "NO", metadata_value))


lp <- ggplot(volc_plot %>%
         filter(metadata_value == "NO") %>%
         filter(ome != omea), 
       aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.05,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value == "NO") %>%
               filter(ome == omea), 
             mapping = aes(size = diffexp),
             shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 0.5,
             stroke = 0.2) +
  geom_point(data = volc_plot %>%
               filter(metadata_value != "NO") %>%
               filter(ome == omea), 
             mapping = aes(effect_size, neglogpvalue,
                           fill = metadata_value,
                           size = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2,
             alpha = 1) +
  scale_fill_manual(values = lighter_shades) +
  scale_size_manual(values = c(0.5, 2), guide = "none") +
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
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
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
lp
ggsave(paste0("reports/figures/Volcano_group_", anal, "_", comp, "_formula", form, "_LipidUnsaturationterms_col.pdf"), 
       width = 7, height = 6, units = "cm")

geom_only_lp <- lp +
  theme_void() +  
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none"
  )
geom_only_lp
ggsave(paste0("reports/figures/Volcano_group_", anal, "_", comp, "_formula", form, "_LipidUnsaturationterms_col_geom.png"), 
       width = 7, height = 6, dpi = 300, units = "cm")


## plot:combined fgsea results all omes ----

lipid_gsea <- fread("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Enrichment/AnalysisGp1/PASC_noPASC/GSEA_enrichment_summary_analysisGroup1_noBatch1_PASC_noPASC_v2.csv")
protein_gsea <- fread("data/processed/fgsea_AnalysisGroup7_PASCnoPASC_protein_allGOterms.csv")
transcript_gsea <- fread("data/processed/fgsea_AnalysisGroup7_PASCnoPASC_transcript_allGOterms.csv")

# split between positive NES/negative NES, and filter for <0.05 qvalue
lipid_gsea <- lipid_gsea %>%
  filter(`FDR q-val` < 0.05) %>%
  filter(Comparison == "PASC_noPASC") %>%
  rename(name = Term,
         pval = `NOM p-val`,
         padj = `FDR q-val`,
         log2err = `FWER p-val`,
         leadingEdge = Lead_genes,
         pathway = Name) %>%
  select(-`Tag %`, -`Gene %`, -Comparison) %>%
  mutate(size = 1,
         ome = "l")

protein_gsea <- protein_gsea %>%
  filter(padj < 0.05) %>%
  mutate(ome = "p")

transcript_gsea <- transcript_gsea %>%
  filter(padj < 0.05) %>%
  mutate(ome = "t")



all_gsea <- lipid_gsea %>%
  bind_rows(protein_gsea) %>%
  bind_rows(transcript_gsea)

# find ones that are in both protein and transcript
all_gsea_shared <- all_gsea %>%
  filter(duplicated(pathway)| duplicated(pathway, fromLast = TRUE))


# # exclude small sets wholly encompassed by a larger leading edge
# # 1: Convert each string to a list of values
# leading_lists <- str_split(all_gsea$leadingEdge, "\\|")
# 
# # 2: Convert to a list of sets for easy comparison
# leading_sets <- map(leading_lists, ~ unique(.x))  
# 
# # Step 3: Create a logical vector marking which rows are *not* subsets of any other row
# keep <- sapply(seq_along(leading_sets), function(i) {
#   current <- leading_sets[[i]]
#   others <- leading_sets[-i]
#   is_subset <- sapply(others, function(x) all(current %in% x))
#  !any(is_subset)
#})

# Step 4: Filter the original dataframe
#filtered_all_gsea <- all_gsea[keep, ]


## Keep specific interesting sets to plot
keep <- c("GO:0005615", "GO:0030545", "GO:0005539", "GO:0005198", "GO:0005509", "GO:0004175", "GO:0006955",
              "GO:0030246", "GO:0008233", "GO:0035639", "GO:0035556", "GO:0008092", "GO:0140993", "GO:0007017",
          "GO:0060589", "GO:0006959", "GO:0180051", "GO:0005788", "GO:0005539", "GO:0032993", "GO:0015629")
filtered_all_gsea <- all_gsea %>%
  filter(ome %in% c("l", "t") | (ome == "p" & pathway %in% keep))


plo <- filtered_all_gsea %>%
  mutate(signif = case_when(
    padj <= 0.05 & padj > 0.01 ~ "*",
    padj <= 0.01 & padj > 0.0001 ~ "**",
    padj <= 0.0001 & padj > 0.000001 ~ "***",
    padj <= 0.000001 ~ "****"
  ))

p_plot <- ggplot(plo %>%
                   filter(ome == "p"), aes(x = NES, y = reorder(name, NES, decreasing = T), fill = ome)) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y") +
  scale_fill_manual(values = c(col[1])) +
  geom_text(aes(label = signif, x = -0.1), hjust = 1.2, size = 3) +
  labs(x = "NES", 
       y = NULL) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "none",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 

l_plot <- ggplot(plo %>%
                   filter(ome == "l"), aes(x = NES, y = reorder(name, NES, decreasing = T), fill = ome)) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y") +
  scale_fill_manual(values = c(col[2])) +
  geom_text(aes(label = signif, x = -0.1), hjust = 1.2, size = 3) +
  labs(x = "NES", 
       y = NULL) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "none",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 

t_plot <- ggplot(plo %>%
                   filter(ome == "t"), aes(x = NES, y = reorder(name, NES, decreasing = T), fill = ome)) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y") +
  scale_fill_manual(values = c(col[3])) +
  geom_text(aes(label = signif, x = -0.1), hjust = 1.2, size = 3) +
  labs(x = "NES", 
       y = NULL) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "none",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 

p_plot + l_plot + t_plot

ggsave("reports/figures/PASCnoPASC_Enrichment_allomes_barplot_filtered_separate_7.pdf", 
       width = 7, height = 2, units = "in")





##* STRING Network Analysis ----
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)



## Make list of genenames in pathway
path_name <- GO_terms %>%
  filter(GO_term == pathway) %>%
  pull(name)
genenames <- fc_barplots %>%
  mutate(gene_map = str_extract(metadata_value, "^[^_]+")) %>%
  select(gene_map) %>%
  distinct() %>%
  pull(gene_map)
writeLines(genenames, paste0("data/processed/GOsignificantnames", path_name, "_PASCnoPASC.txt"))





string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

mapped_genes <- string_db$map(fc_barplots, 
                              "metadata_value", 
                              removeUnmappedRows = TRUE)

# Get STRING interaction network for those proteins
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Create a graph object
graph <- graph_from_data_frame(interactions, directed = FALSE)

# Add attributes
V(graph)$gene <- mapped_genes$metadata_value[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$ome <- mapped_genes$ome[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$effect_size <- mapped_genes$effect_size[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$neglogpvalue <- mapped_genes$neglogpvalue[match(V(graph)$name, mapped_genes$STRING_id)]

# Convert to tidygraph format
graph_tbl <- as_tbl_graph(graph)

# Plot with ggraph
set.seed(202)
ggraph(graph_tbl, layout = "fr") +  # "fr" is Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.25, color = "gray60", edge_width = 0.2) +
  geom_node_point(aes(fill = effect_size, size = neglogpvalue ), 
                  shape = 21, 
                  color = "black", 
                  stroke = 0.2) +
  geom_node_text(aes(label = gene), 
                 repel = TRUE, 
                 size = 1.75,
                 lineheight = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Effect Size") +
  scale_size_continuous(name = "-log10(p value)", range = c(1, 4)) +
  theme_void() +
  theme(legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.position = "inside",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        plot.title = element_text(size = 7)
  ) +
  ggtitle(direc)
ggsave(paste0('reports/figures/noLC_LC_STRINGnetwork_', path_name, "_", direc, '.pdf'), 
       width = 16, height = 8, units = 'cm')



##* POTEE, POTEI, POTEJ Interactors ----
pote_int <- read_lines("data/processed/POTEE_POTEI_POTEJ_string_interactions_50interactors.txt")

bmol_genes <- biomolecules_metadata %>%
  select(-metadata_id) %>%
  filter(metadata_type %in% c("gene_name", "gene_symbol")) %>%
  group_by(biomolecule_id) %>%
  filter(n() == 1 | metadata_type == "gene_symbol") %>%
  ungroup() %>%
  select(biomolecule_id, metadata_value)

# organize dataframe
string_plots <- volc_plot %>%
  filter(q_value < 0.05) %>%
  left_join(bmol_genes,
            by = "biomolecule_id") %>%
  filter(metadata_value %in% pote_int)


# string mapping
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

mapped_genes <- string_db$map(string_plots, 
                              "metadata_value", 
                              removeUnmappedRows = TRUE)

# Get STRING interaction network for those proteins
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Create a graph object
graph <- graph_from_data_frame(interactions, directed = FALSE)

# Add attributes
V(graph)$gene <- mapped_genes$metadata_value[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$ome <- mapped_genes$ome[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$effect_size <- mapped_genes$effect_size[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$neglogpvalue <- mapped_genes$neglogpvalue[match(V(graph)$name, mapped_genes$STRING_id)]

# Convert to tidygraph format
graph_tbl <- as_tbl_graph(graph)

# Plot with ggraph
set.seed(202)
ggraph(graph_tbl, layout = "fr") +  # "fr" is Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.25, color = "gray60", edge_width = 0.2) +
  geom_node_point(aes(fill = effect_size, size = neglogpvalue ), 
                  shape = 21, 
                  color = "black", 
                  stroke = 0.2) +
  geom_node_text(aes(label = gene), 
                 repel = TRUE, 
                 size = 1.75,
                 lineheight = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Effect Size") +
  scale_size_continuous(name = "-log10(p value)", range = c(1, 4)) +
  theme_void() +
  theme(legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.position = "inside",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        plot.title = element_text(size = 7)
  ) +
  ggtitle(direc)
ggsave('reports/figures/noLC_LC_STRINGnetwork_POTEEPOTEIPOTEJ_50interactions.pdf', 
       width = 10, height = 6, units = 'cm')




## RNAseq cell cycle ----
s_phase_genes <- c("MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", 
                   "PCNA", "TYMS", "RRM1", "RRM2", "CDK2", "CCNE1", 
                   "CCNE2", "CHEK1")

bmol_genes <- biomolecules_metadata %>%
  select(-metadata_id) %>%
  filter(metadata_type %in% c("gene_name", "gene_symbol")) %>%
  group_by(biomolecule_id) %>%
  filter(n() == 1 | metadata_type == "gene_symbol") %>%
  ungroup() %>%
  select(biomolecule_id, metadata_value)


# transcript volc_plot with gene names
cell_cycle_transcripts <- volc_plot %>%
  filter(ome == "transcript") %>%
  left_join(bmol_genes,
            by = "biomolecule_id") %>%
  filter(metadata_value %in% s_phase_genes)

# volc plot
ggplot(cell_cycle_transcripts, 
             aes(effect_size, neglogpvalue)) + 
  geom_point(shape = 21,
             color = "black",
             fill = "lightgray",
             alpha = 1,
             stroke = 0.2) +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_x_continuous(limits = c(-max(abs(cell_cycle_transcripts$effect_size)), max(abs(cell_cycle_transcripts$effect_size)))) +
  geom_text_repel(data = cell_cycle_transcripts, aes(label = metadata_value), size = 2) +
  xlab(paste("Effect Size", comp)) +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")) 
ggsave(paste0("reports/figures/Volcano_group1_protein_cellcyclemarkers_col.pdf"), 
       width = 4, height = 3, units = "in")




## RNAseq leukocyte populations ----
library(immunedeconv)

bmol_genes <- biomolecules_metadata %>%
  select(-metadata_id) %>%
  filter(metadata_type %in% c("gene_name", "gene_symbol")) %>%
  group_by(biomolecule_id) %>%
  filter(n() == 1 | metadata_type == "gene_symbol") %>%
  ungroup() %>%
  select(biomolecule_id, metadata_value)

rnaseq_matrix <- filtered_df_t %>%
  left_join(bmol_genes, by = "biomolecule_id") %>%
  select(sample_id, metadata_value, Counts) %>%
  filter(!is.na(metadata_value)) %>%
  pivot_wider(names_from = sample_id,
              values_from = Counts) %>%
  column_to_rownames(var = "metadata_value")
  

##* calculate TPM for RNAseq data ----
# get gene lengths
library(biomaRt)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
transcript_lengths <- getBM(attributes = c("hgnc_symbol", "start_position", "end_position"), mart = mart) %>%
  mutate(length_kb = ((end_position - start_position) / 1000)) %>% 
  filter(if_all(everything(), ~ . != "")) %>%
  filter(hgnc_symbol %in% rownames(rnaseq_matrix)) %>%
  group_by(hgnc_symbol) %>%
  slice_max(length_kb, with_ties = FALSE) %>%
  ungroup()

gene_lengths <- setNames(transcript_lengths$length_kb, 
                         transcript_lengths$hgnc_symbol)

# get rnaseq matrix organized
rnaseq_matrix1 <- rnaseq_matrix %>%
  rownames_to_column(var = "ID") %>%
  filter(ID %in% transcript_lengths$hgnc_symbol) %>%
  column_to_rownames(var = "ID")

# do TPM calc
calculate_tpm <- function(counts, gene_lengths) {
  # Compute RPK
  rpk <- sweep(counts, 1, gene_lengths, "/")
  
  # Compute per-sample scaling factor
  scaling_factors <- colSums(rpk)
  
  # Compute TPM
  tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6
  return(tpm)
}

tpm_matrix <- calculate_tpm(counts = rnaseq_matrix1, 
                            gene_lengths = gene_lengths)

##* Plot cell types ----

immune_cell_types <- deconvolute(tpm_matrix, "quantiseq") %>%
  pivot_longer(cols = -cell_type,
               names_to = "sample_id",
               values_to = "proportion") %>%
  mutate(sample_id = as.integer(sample_id)) %>%
  left_join(metadata, by = "sample_id") %>% 
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID",
    Cohort %in% c("PASC", "PASC_fu") ~ "Long_COVID")) 

ggplot(immune_cell_types, aes(proportion, factor(sample_id), fill = cell_type)) + 
  geom_col(position = "stack",
           width = 0.6) + 
  scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = "proportion",
       y = "sample_id") +
  #scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(title = 'Biomolecule')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "right", 
        legend.justification = c("right", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")
  ) +
  facet_wrap(PASCnoPASC ~ ., scales = "free_y")
ggsave(paste0('reports/figures/PASCnoPASC_Transcript_ImmuneCellTypes.pdf'), 
       width = 12, height = 24, units = 'cm')




