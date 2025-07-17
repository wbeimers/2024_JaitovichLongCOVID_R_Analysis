#### Overview ####

# 1. Healthy vs Acute linear model (accounting for Age, Sex, BMI variables)
# 2. See if any differentially expressed biomolecules are retained in the PASC-noPASC comparsion we did (in Analysis Group 1)


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
library(ROTS)
library(patchwork)
library(pheatmap)
library(broom)
library(fgsea)


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


col1 <- viridis_pal(option = "rocket")(100)[round(c(0.25, 0.5, 0.75) * 100)]


# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                               FROM proteomics_measurement')
lipidomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                              FROM lipidomics_measurements')
transcriptomics <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, sample_id, Counts, normalized_counts
                                    FROM rnaseq_measurements")
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_id, rawfile_name, Sample, sample_id, run_type, ome_id, keep, batch
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')

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


# combined lipid and protein (no transcript in Acute samples)
df_pl <- df_p %>%
  mutate(ome = "p") %>%
  bind_rows(df_l %>%
              mutate(ome = "l"))

filtered_df_pl <- filtered_df_p %>%
  mutate(ome = "p") %>%
  bind_rows(filtered_df_l %>%
              mutate(ome = "l")) %>%
  filter(batch != 1) %>% # remove samples from Batch 1
  mutate(PASCnoPASC = case_when(
    Cohort %in% c("Acute") ~ "Acute_COVID",
    Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID",
    Cohort %in% c("PASC") ~ "Long_COVID")) 


## FILTER OUT LIPID BATCH 1 SAMPLES????


## Acute vs Healthy Samples ----
# Get grouping vector and data matrix
# Groups: Acute, Acute_fu, Acute_NC, Healthy, PASC, PASC_fu
# PASC_Cohort: first, second
# PG_change_collection_cutoff: 0, 1
group1 <- "Acute"
group2 <- "Healthy"

diffexp_df <- filtered_df_pl %>%
  filter(Cohort %in% c(group1, group2)) %>%
  select(standardized_name, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = "sample_id", values_from = "normalized_abundance") %>%
  tibble::column_to_rownames(var = "standardized_name") %>%
  select(where(~ all(!is.na(.))))

grouping_vector <- metadata %>%
  filter(sample_id %in% colnames(diffexp_df)) %>%
  arrange(match(sample_id, colnames(diffexp_df))) %>%
  pull(Cohort)


##* ROTS DEA ----
# run ROTS
results <- ROTS(data = diffexp_df, 
                groups = grouping_vector, 
                B = 5000, 
                K = 500, 
                seed = 1234)

summary(results, fdr = 0.05)

# output df
volc_plot <- data.frame(results$logfc) %>%
  tibble::rownames_to_column(var = "standardized_name") %>%
  mutate(logfc = results.logfc) %>%
  dplyr::select(-results.logfc) %>%
  mutate(pvalue = results$pvalue) %>%
  mutate(qvalue = results$FDR)

bmol_ids <- filtered_df_pl %>%
  select(biomolecule_id) %>%
  distinct() %>%
  pull(biomolecule_id)

volc_plot <- volc_plot %>%
  mutate(neglogpvalue = -log10(pvalue)) %>%
  mutate(diffexp = case_when(
    logfc >= 0.263 & qvalue <= 0.05 ~ "UP",
    logfc <= -0.263 & qvalue <= 0.05 ~ "DOWN",
    T ~ "NO"
  )) %>%
  left_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "standardized_name") %>%
  filter(biomolecule_id %in% bmol_ids)

fwrite(volc_plot, paste0("data/processed/ProteinLipid_ROTS_Volcano_", group1, "_", group2, "_nobatch1.csv"))

volc_plot <- fread(paste0("data/processed/AllPlates_Samples_Volcano_", group1, "_", group2, ".csv"))

counts <- volc_plot %>%
  filter(diffexp %in% c("DOWN", "UP", "NO")) %>%
  count(diffexp)

ggplot(volc_plot, aes(logfc, neglogpvalue)) + 
  geom_point(aes(size = diffexp, fill = as.factor(omics_id), alpha = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2) +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_fill_manual(values = c(col[1], col[2])) +
  scale_size_manual(values = c(1.5, 0.5, 1.5), guide = "none") +
  scale_alpha_manual(values = c(0.8, 0.2, 0.8), guide = "none") +
  scale_x_continuous(limits = c(-max(abs(volc_plot$logfc)), max(abs(volc_plot$logfc)))) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab("Log2(FC) Acute/Healthy") +
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
        legend.key.size = unit(0.25, "cm")) +
  #facet_wrap(. ~ factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")),
  #           scales = "free",
  #           ncol = 2) +
  #geom_text(data = volc_plot %>% distinct(method),
  #          aes(x = -3.5, y = 8, label = method),
  #          size = 3, 
  #          color = "black", 
  #          hjust = 0) +
  geom_text(data = counts[counts$diffexp == "DOWN",], 
            aes(x = -3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "UP",], 
            aes(x = 3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "NO",], 
            aes(x = 0, y = Inf, 
                label = paste(n)),
            hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) 
ggsave(paste0("reports/figures/ProteinLipid_ROTS_50pconditionMissingImputed_Volcano_small_", group1, "_", group2, "_omics_id_nobatch1.pdf"), 
       width = 8, height = 6, units = "cm")


## Similar Acute vs PASC features ----
# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

pvalues <- dbGetQuery(con, 'SELECT biomolecule_id, analysis_group, test, comparison, formula, predictor, effect_size, eta_squared, lratio, p_value, q_value
                            FROM pvalues')
formulas <- dbGetQuery(con, 'SELECT *
                             FROM formula_table')
dbDisconnect(con)

# options:
# analysis_group (1, 2, 3, 0)
anal <- 7
# comparison (Age, Sex, QoL, BMI)
comp <- "group7_PASCnoPASC"
# formula (1, 2, 3, etc.)
form <- 57

volc_plot_PASC <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 & effect_size > 0 ~ "UP",
    q_value <= 0.05 & effect_size < 0 ~ "DOWN",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id)


# options:
# analysis_group (1, 2, 3, 0)
anal <- 6
# comparison (Age, Sex, QoL, BMI)
comp <- "Acute_nonAcute"
# formula (1, 2, 3, etc.)
form <- 53

volc_plot_Acute <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(effect_size = -(effect_size)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 & effect_size > 0 ~ "UP",
    q_value <= 0.05 & effect_size < 0 ~ "DOWN",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id) 

# find overlap of UP and DOWN features between both Acute and PASC

UP_df <- volc_plot_PASC %>%
  filter(diffexp == "UP") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "UP"),
             by = "biomolecule_id")

DOWN_df <- volc_plot_PASC %>%
  filter(diffexp == "DOWN") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "DOWN"),
             by = "biomolecule_id")

UP_DOWN_together <- UP_df %>%
  select(biomolecule_id, ome.x, diffexp.x) %>%
  bind_rows(DOWN_df %>%
              select(biomolecule_id, ome.x, diffexp.x))

# make euler plots of up and down overlaps
library(eulerr)
up_overlaps <- list(
  UP_Acute = volc_plot_Acute %>%
    filter(diffexp == "UP") %>%
    pull(biomolecule_id),
  UP_PASC = volc_plot_PASC %>%
    filter(diffexp == "UP") %>%
    pull(biomolecule_id)
)

euler_ov <- euler(up_overlaps)

pdf("reports/figures/EulerPlot_AcutePASC_UP_Overlap.pdf", width = 2, height = 2)
plot(euler_ov, 
     fills = c(col[4], col[5], col[6]), 
     edges = T,
     quantities = T)
dev.off()


down_overlaps <- list(
  DOWN_Acute = volc_plot_Acute %>%
    filter(diffexp == "DOWN") %>%
    pull(biomolecule_id),
  DOWN_PASC = volc_plot_PASC %>%
    filter(diffexp == "DOWN") %>%
    pull(biomolecule_id)
)

euler_ov <- euler(down_overlaps)

pdf("reports/figures/EulerPlot_AcutePASC_DOWN_Overlap.pdf", width = 2, height = 2)
plot(euler_ov, 
     fills = c(col[4], col[5], col[6]), 
     edges = T,
     quantities = T)
dev.off()


##* SingleBiomoleculePlots-AcutePASCOverlap ----

comparison_df_pl <- filtered_df_pl %>%
  filter(Cohort %in% c("Healthy", "Acute", "PASC")) %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  filter(PG_change_collection_cutoff == 0)

for (i in unique(comparison_df_pl$standardized_name)) {
  
  poi <- i
  
  single_feat_df <- comparison_df_pl %>%
    filter(standardized_name == poi) %>%
    mutate(Cohort = factor(Cohort, levels = c("Acute", "PASC", "Healthy")))
  
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
  
  ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AcutePASCcomparison_PreTube_singlebiomolecule_', poi, '_distribution_Cohort.pdf'), 
         width = 8, height = 6, units = "cm")
}


##* GO:terms Proteins ----

# select limited GO_set
#GO_set <- fread("data/metadata/GOtermset_HealthyAcutePASCnoPASC_Proteins_BPonly_5to100_95coverage_20overlap_UP.csv")

# select only one of the GO types
#GO_set <- GO_terms %>%
#  filter(namespace == "biological_process")


## Acute PASC Overlap
# Set up GO term
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                                          FROM biomolecules_metadata')
GO_terms <- dbGetQuery(con, 'SELECT *
                             FROM GO_terms')
dbDisconnect(con)


# Make background list of all proteins in our study
background_list <- filtered_df_p %>%
  select(biomolecule_id) %>%
  distinct() %>%
  pull(biomolecule_id)

## Katie/Salma Function for testing enrichment ##
enrichment <- function(set, reference_sets, background){
  # Make variables
  nset <- length(set)
  nbackground <- length(background)
  
  # Initialize output dataframe
  output <- data.frame(reference = names(reference_sets), 
                       pvalue = rep(1, length(names(reference_sets))), 
                       fdr_pvalue = rep(NA, length(names(reference_sets))),
                       enrichment_ratio = rep(NA, length(names(reference_sets))),
                       stringsAsFactors = FALSE)
  
  # Loop through each reference set
  for (i in 1:nrow(output)){
    # Calculate observed hits (features in 'set' present in the reference set)
    observed_hits <- length(intersect(set, reference_sets[[output[i, "reference"]]]))
    
    # Calculate expected hits by chance (total features in reference set * total features in 'set' / total features in background)
    expected_hits <- length(reference_sets[[output[i, "reference"]]]) * nset / nbackground
    
    # Calculate enrichment ratio
    enrichment_ratio <- observed_hits / expected_hits
    
    # Calculate Background Hits
    hitsBackground <- length(intersect(background, reference_sets[[output[i,1]]]))
    
    # Calculate p-value using hypergeometric distribution
    p_value <- phyper(observed_hits - 1, hitsBackground, 
                      nbackground - hitsBackground, 
                      nset, lower.tail = FALSE)
    
    # Update output dataframe with results
    output[i, "pvalue"] <- p_value
    output[i, "enrichment_ratio"] <- enrichment_ratio
  }
  
  # Adjust p-value for multiple testing using Benjamini-Hochberg method
  output$fdr_pvalue <- p.adjust(output$pvalue, method = "BH")
  
  # Return the output dataframe
  return(output)
} 

## UP GENES
# GO term list
GO_set <- fread("data/metadata/GOtermset_HealthyAcutePASCnoPASC_Proteins_BPonly_3to50_95coverage_10overlap_UP.csv")

biomolecule_metadata_GO <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  separate_rows(metadata_value, sep = ";") %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id), .groups = "drop") #%>% 
  #filter(metadata_value %in% GO_set$GO_term)

GO_term_list <- split(biomolecule_metadata_GO$GO_terms,
                      biomolecule_metadata_GO$metadata_value)

# Gene list
gene_list <- UP_DOWN_together %>%
  filter(diffexp.x == "UP") %>%
  filter(ome.x == "protein") %>%
  pull(biomolecule_id)

# gene_list <- HighOpposite_together %>%
#   filter(ome.x == "protein") %>%
#   pull(biomolecule_id)

# Enrichment
enrichment_df_UP <- enrichment(gene_list, GO_term_list, background_list) %>%
  left_join(GO_terms, by = c("reference" = "GO_term")) %>%
  mutate(rank = enrichment_ratio * -log10(pvalue))

# enrichment_df_opposite <- enrichment(gene_list, GO_term_list, background_list) %>%
#   left_join(GO_terms, by = c("reference" = "GO_term")) %>%
#   mutate(rank = enrichment_ratio * -log10(pvalue))


## DOWN GENES
# GO term list
GO_set <- fread("data/metadata/GOtermset_HealthyAcutePASCnoPASC_Proteins_BPonly_3to50_95coverage_10overlap_DOWN.csv")

biomolecule_metadata_GO <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  separate_rows(metadata_value, sep = ";") %>%
  group_by(metadata_value) %>%
  summarize(GO_terms = unique(biomolecule_id), .groups = "drop") #%>% 
  #filter(metadata_value %in% GO_set$GO_term)

GO_term_list <- split(biomolecule_metadata_GO$GO_terms,
                      biomolecule_metadata_GO$metadata_value)

# Gene list
gene_list <- UP_DOWN_together %>%
  filter(diffexp.x == "DOWN") %>%
  filter(ome.x == "protein") %>%
  pull(biomolecule_id)

# Enrichment
enrichment_df_DOWN <- enrichment(gene_list, GO_term_list, background_list) %>%
  left_join(GO_terms, by = c("reference" = "GO_term")) %>%
  mutate(rank = enrichment_ratio * -log10(pvalue))



## ALL
# Gene list
#gene_list <- UP_DOWN_together %>%
#  filter(ome.x == "protein") %>%
#  pull(biomolecule_id)

# Enrichment
#enrichment_df_ALL <- enrichment(gene_list, GO_term_list, background_list) %>%
#  left_join(GO_terms, by = c("reference" = "GO_term"))







## Plot enrichment barplots
enrichment_df_UP_plot <- enrichment_df_UP %>%
  mutate(neglog10 = -log10(fdr_pvalue)) %>%
  filter(neglog10 >= 1.301) %>%
  mutate(direction = "UP")

enrichment_df_DOWN_plot <- enrichment_df_DOWN %>%
  mutate(neglog10 = -log10(fdr_pvalue)) %>%
  filter(neglog10 >= 1.301) %>%
  mutate(direction = "DOWN")

enrichment_df_plot_combined <- enrichment_df_UP_plot %>%
  bind_rows(enrichment_df_DOWN_plot) %>%
  mutate(signif = case_when(
    fdr_pvalue <= 0.05 & fdr_pvalue > 0.01 ~ "*",
    fdr_pvalue <= 0.01 & fdr_pvalue > 0.001 ~ "**",
    fdr_pvalue <= 0.001 & fdr_pvalue > 0.00001 ~ "***",
    fdr_pvalue <= 0.00001 ~ "****"
  )) %>%
  arrange(desc(neglog10))

ggplot(enrichment_df_plot_combined, aes(x = enrichment_ratio, y = reorder(name, enrichment_ratio), fill = direction)) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y") +
  scale_fill_manual(values = c(col[4], col[5])) +
  geom_text(aes(label = signif, x = 0.1), size = 4) +
  labs(x = "Enrichment Ratio", 
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
ggsave("reports/figures/AC_LC_nC_overlap_UPDOWN_pvalues.pdf", 
       width = 6, height = 8, units = "cm")


##* Volcano Plots - Colored by Overlap ----

volc_plot_Acute_plot <- volc_plot_Acute %>%
  mutate(shared = case_when(
    biomolecule_id %in% UP_df$biomolecule_id ~ "UP",
    biomolecule_id %in% DOWN_df$biomolecule_id ~ "DOWN",
    T ~ "NO"
  )) %>%
  filter(ome %in% c("protein", "lipid"))

counts <- volc_plot_Acute_plot %>%
  filter(shared %in% c("UP", "DOWN")) %>%
  count(shared)

ac_p <- ggplot(volc_plot_Acute_plot, aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = shared, fill = shared, alpha = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2) +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  #scale_fill_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
  scale_fill_manual(values = c(col1[3], "lightgray", col1[1])) +
  scale_size_manual(values = c(1.5, 0.5, 1.5), guide = "none") +
  scale_alpha_manual(values = c(0.8, 0.2, 0.8), guide = "none") +
  scale_x_continuous(limits = c(-max(abs(volc_plot_Acute_plot$effect_size)), max(abs(volc_plot_Acute_plot$effect_size)))) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab("Effect Size") +
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
        legend.key.size = unit(0.25, "cm")) +
  geom_text(data = counts[counts$shared == "UP",], 
            aes(x = 1, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 2, show.legend = FALSE) +
  geom_text(data = counts[counts$shared  == "DOWN",], 
            aes(x = -1, y = Inf, 
                label = paste(n)),
            hjust = 0.5, vjust = 1.5, size = 2, show.legend = FALSE) 
ac_p
ggsave("reports/figures/ProteinLipid_LRM_Acute_noAcute_Volcano_small_sharedPASC_nobatch1.pdf", 
       width = 8, height = 6, units = "cm")


volc_plot_PASC_plot <- volc_plot_PASC %>%
  mutate(shared = case_when(
    biomolecule_id %in% UP_df$biomolecule_id ~ "UP",
    biomolecule_id %in% DOWN_df$biomolecule_id ~ "DOWN",
    T ~ "NO"
  )) %>%
  filter(ome %in% c("protein", "lipid"))

counts <- volc_plot_PASC_plot %>%
  filter(shared %in% c("UP", "DOWN")) %>%
  count(shared)

lc_p <- ggplot(volc_plot_PASC_plot, aes(effect_size, neglogpvalue)) + 
  geom_point(aes(size = shared, fill = shared, alpha = diffexp),
             shape = 21,
             color = "black",
             stroke = 0.2) +
  #geom_vline(xintercept = c(-0.263, 0.263), 
  #          col = "black",
  #           size = 0.2) +
  #geom_hline(yintercept = -log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_fill_manual(values = c(col1[3], "lightgray", col1[2])) +
  scale_size_manual(values = c(1.5, 0.5, 1.5), guide = "none") +
  scale_alpha_manual(values = c(0.8, 0.2, 0.8), guide = "none") +
  scale_x_continuous(limits = c(-max(abs(volc_plot_PASC_plot$effect_size)), max(abs(volc_plot_PASC_plot$effect_size)))) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab("Effect Size") +
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
        legend.key.size = unit(0.25, "cm")) +
  geom_text(data = counts[counts$shared == "UP",], 
            aes(x = 1, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 2, show.legend = FALSE) +
  geom_text(data = counts[counts$shared  == "DOWN",], 
            aes(x = -1, y = Inf, 
                label = paste(n)),
            hjust = 0.5, vjust = 1.5, size = 2, show.legend = FALSE) 
lc_p
ggsave("reports/figures/ProteinLipid_LRM_LC_noLC_Volcano_small_sharedPASC_nobatch1.pdf", 
       width = 8, height = 6, units = "cm")

ac_p + lc_p
ggsave("reports/figures/ProteinLipid_LRM_LC_noLC_Volcano_small_sharedBOTHPLOTS_nobatch1_gorup7.pdf", 
       width = 7, height = 2, units = "in")

##* Heatmap - overlapping biomolecules ----

pal <- c("#E78AC3", 
         '#AA3377', 
         "#B3B3B3", 
         '#229100', 
         '#66CCEE', 
         '#4477AA')

expression_matrix <- filtered_df_pl %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  filter(Cohort %in% c("Healthy", "Acute", "PASC")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  dplyr::select(biomolecule_id, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  tibble::column_to_rownames(var = "biomolecule_id") %>%
  select(where(~ !any(is.na(.))))

sample_annot <- filtered_df_pl %>%
  ungroup() %>%
  select(sample_id, 
         Cohort, 
         #Age, 
         #Sex, 
         #BMI, 
         #SF.36.QOL.Score, 
         #PASC_Cohort, 
         #PG_change_collection_cutoff
  ) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id")

dat <- filtered_df_pl %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  select(ome, biomolecule_id) %>%
  distinct()

row_annot <- data.frame(ome = dat$ome,
                        row.names = dat$biomolecule_id) %>%
  mutate(ome1 = case_when(
    ome == "p" ~ "protein",
    ome == "l" ~ "lipid",
    ome == "t" ~ "transcript"
  )) %>%
  select(-ome) %>%
  rename(Ome = ome1)

k <- 2

pheat <- pheatmap(expression_matrix,
                  #color = viridis(24, direction = 1, option = "plasma"),
                  color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(15)),
                  breaks = c(-4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4),
                  cluster_rows = T,
                  #cutree_rows = k, 
                  #gaps_row = T,
                  cluster_cols = T,
                  treeheight_row = 0,
                  treeheight_col = 10,
                  show_rownames = F,
                  show_colnames = F,
                  border_color = NA,
                  scale = "row",
                  annotation_row = row_annot,
                  annotation_col = sample_annot,
                  annotation_colors = list(Cohort = c(Acute = pal[1], Healthy = pal[4], PASC = pal[5]),
                                           Ome = c(protein = col[1], lipid = col[2], transcript = col[3])),
                  filename = paste0("reports/figures/Heatmap_plots/Heatmap_AcutePASCoverlaptogether_biomolecules.png"),
                  width = 6,
                  height = 7.5)

cluster_assignments <- cutree(pheat$tree_row, k = k)
fwrite(as.data.frame(cluster_assignments) %>%
         tibble::rownames_to_column(var = "biomolecule_id"), paste0("reports/figures/Heatmap_plots/Clusters_AllOmes_AnalysisGroup1_4methodoverlap_2clusters.csv"))



##* PCA ----

# take normalized abundance values for PCA
pca_set <- filtered_df_pl %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  filter(PASCnoPASC %in% c("Acute_COVID", "No_COVID", "Long_COVID")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  dplyr::select(biomolecule_id, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  select(where(~ !any(is.na(.))))

t_pca_set <- as.data.frame(t(pca_set))

# turn row into colnames
colnames(t_pca_set) <- as.character(t_pca_set[1, ])
# Remove the row
t_pca_set <- t_pca_set[-1, ]
#change to numeric
t_pca_set <- data.frame(lapply(t_pca_set, as.numeric), row.names = rownames(t_pca_set))
# re-add sample_id
t_pca_set <- tibble::rownames_to_column(t_pca_set, "sample_id") %>%
  mutate(sample_id = as.integer(sample_id))

pca_score <- prcomp(t_pca_set[,c(2:ncol(t_pca_set))],
                    scale. = T)
summary(pca_score)

# make variables to plot
# scree
explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
variance <- data.frame(proportion = explained_variance,
                       PC = 1:(ncol(pca_set)-1))

# pca scores
scores <- as.data.frame(pca_score$x) %>%
  mutate(sample_id = t_pca_set$sample_id) %>%
  left_join(filtered_df_pl %>% 
              select(sample_id, Cohort, PASCnoPASC, Age, Sex, BMI, unique_patient_id, Collection_date) %>%
              distinct(), 
            by = "sample_id")

# loadings
loadings <- pca_score$rotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = 'biomolecule_id')

# pca plot
p <- ggplot(scores, aes(PC1, PC2, fill = PASCnoPASC)) + 
  geom_point(shape = 21,
             size = 1.5,
             alpha = 0.9,
             color = "black",
             stroke = 0.1) +
  stat_ellipse(aes(color = PASCnoPASC), 
               geom = "path", 
               show.legend = FALSE,
               linewidth = 0.2) +
  #geom_text_repel(aes(label = Sample), size = 2) +
  #scale_fill_manual(values = pal) +
  scale_fill_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
  scale_color_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
  xlab(paste("PC1", round(variance$proportion[1]*100, 2))) +
  ylab(paste("PC2", round(variance$proportion[2]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "inside", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )

p
ggsave(paste0("reports/figures/PCA_plots/PCA_AC_LC_nC_OverlapBiomolecules_proteinlipid_nobatch1_pretube_group7.pdf"), 
       width = 8, height = 5, units = "cm")

# loadings plot
loadings_1 <- loadings %>%
  mutate(biomolecule_id = as.integer(str_remove(biomolecule_id, "^X"))) %>%
  left_join(biomolecules %>%
              select(biomolecule_id, standardized_name),
            by = "biomolecule_id")

p1 <- ggplot(loadings_1, aes(PC1, PC2)) + 
  geom_point(shape = 21,
             size = 1,
             color = "black",
             stroke = 0.1) +
  geom_text_repel(aes(label = standardized_name), size = 2) +
  xlab(paste("PC1", round(variance$proportion[1]*100, 2))) +
  ylab(paste("PC2", round(variance$proportion[2]*100, 2))) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "right", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )

p1
ggsave("reports/figures/PCA_plots/Loadings_AC_LC_nC_OverlapBiomolecules_proteinlipid_nobatch1_pretube.pdf", 
       width = 8, height = 6, units = "cm")



## Dial into specific pathways that are enriched ----
# find a pathway, find proteins, and plot log2fc for both PASCnoPASC and Acute/Healthy

pathway <- "GO:0003676"

pathway_bms <- GO_term_list[[pathway]]

# organize dataframes
volc_plot_Acute_1 <- volc_plot_Acute %>%
  filter(biomolecule_id %in% pathway_bms) %>%
  filter(q_value < 0.05)

volc_plot_PASC_1 <- volc_plot_PASC %>%
  filter(biomolecule_id %in% pathway_bms) %>%
  filter(q_value < 0.05)


fc_barplots <- volc_plot_Acute_1 %>%
  bind_rows(volc_plot_PASC_1) %>%
  filter(ome == "protein")
# filter by shared biomolecules
shared_ids <- fc_barplots %>%
  distinct(comparison, biomolecule_id) %>%
  group_by(biomolecule_id) %>%
  summarise(n_cat = n_distinct(comparison), .groups = "drop") %>%
  filter(n_cat == 2) %>%
  pull(biomolecule_id)
fc_barplots <- fc_barplots %>%
  filter(biomolecule_id %in% shared_ids) %>%
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "Protein_names") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id")


ggplot(fc_barplots, aes(reorder(metadata_value, q_value), effect_size, fill = comparison)) + 
  geom_col(position = position_dodge(),
           width = 0.8) + 
  geom_hline(yintercept = 0,
             linewidth = 0.2) +
  scale_fill_manual(values = c(pal1[5], pal1[3])) +
  labs(x = NULL,
       y = 'Effect Size') +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title = 'Comparison')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.95, 0.95), 
        legend.justification = c("right", "top"),
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")
  ) 
ggsave(paste0('reports/figures/AcutePASC_OppositeBiomolecules_proteinlipid_GO0007219_foldchanges.pdf'), 
       width = 16, height = 6, units = 'cm')




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
  select(metadata_value) %>%
  distinct() %>%
  pull(metadata_value)
writeLines(genenames, paste0("data/processed/GOsignificantnames", path_name, ".txt"))




## CHOOSE DIRECTION "Acute_nonAcute" or "PASC_noPASC"
direc <- "PASC_noPASC"


string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

mapped_genes <- string_db$map(UPDOWN_genes, 
                              "metadata_value", 
                              removeUnmappedRows = TRUE)

# Get STRING interaction network for those proteins
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Create a graph object
graph <- graph_from_data_frame(interactions, directed = FALSE)

# Add attributes
V(graph)$gene <- mapped_genes$metadata_value[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$effect_size <- mapped_genes$effect_size[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$neglogpvalue <- mapped_genes$neglogpvalue[match(V(graph)$name, mapped_genes$STRING_id)]

# Convert to tidygraph format
graph_tbl <- as_tbl_graph(graph)

# Plot with ggraph
set.seed(1997)
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
ggsave(paste0('reports/figures/AC_noLC_LC_STRINGnetwork_', path_name, "_", direc, '.pdf'), 
       width = 5, height = 6, units = 'cm')


##* STRING: ALL SHARED GENES ----
## CHOOSE DIRECTION "Acute_nonAcute" or "group7_PASCnoPASC"
direc <- "Acute_nonAcute"

# find shared
UPDOWN_genes <- UP_DOWN_together %>%
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "gene_name") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id") %>%
  pull(metadata_value)
UPDOWN_genes <- UPDOWN_genes[!is.na(UPDOWN_genes)]


# organize dataframes
volc_plot_Acute_1 <- volc_plot_Acute %>%
  filter(q_value < 0.05)

volc_plot_PASC_1 <- volc_plot_PASC %>%
  filter(q_value < 0.05)

shared_genes <- volc_plot_Acute_1 %>%
  bind_rows(volc_plot_PASC_1) %>%
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "gene_name") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id") %>%
  filter(ome == "protein") %>%
  filter(metadata_value %in% UPDOWN_genes) %>%
  filter(comparison == direc)


string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

mapped_genes <- string_db$map(shared_genes, 
                              "metadata_value", 
                              removeUnmappedRows = TRUE)

# Get STRING interaction network for those proteins
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Create a graph object
graph <- graph_from_data_frame(interactions, directed = FALSE)

# Add attributes
V(graph)$gene <- mapped_genes$metadata_value[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$effect_size <- mapped_genes$effect_size[match(V(graph)$name, mapped_genes$STRING_id)]
V(graph)$neglogpvalue <- mapped_genes$neglogpvalue[match(V(graph)$name, mapped_genes$STRING_id)]

# Convert to tidygraph format
graph_tbl <- as_tbl_graph(graph)

# Plot with ggraph
set.seed(1997)
ggraph(graph_tbl, layout = "fr") +  # "fr" is Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.25, color = "gray60", edge_width = 0.2) +
  geom_node_point(aes(fill = effect_size, size = neglogpvalue ), 
                  shape = 21, 
                  color = "black", 
                  stroke = 0.2) +
  geom_node_text(aes(label = gene), 
                 repel = TRUE, 
                 size = 1.75,
                 lineheight = 0.2,
                 segment.size = 0.2) +
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
ggsave(paste0('reports/figures/AC_noLC_LC_STRINGnetwork_ALLSHAREDGENES_', direc, '_small_group7.pdf'), 
       width = 7, height = 8, units = 'in')







## Opposite Acute vs PASC features ----
# find overlap of UP vs DOWN features between both Acute and PASC

HighPASC_df <- volc_plot_PASC %>%
  filter(diffexp == "UP") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "DOWN"),
             by = "biomolecule_id")

HighAcute_df <- volc_plot_PASC %>%
  filter(diffexp == "DOWN") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "UP"),
             by = "biomolecule_id")

HighOpposite_together <- HighPASC_df %>%
  select(biomolecule_id, ome.x, diffexp.x, diffexp.y) %>%
  bind_rows(HighAcute_df %>%
              select(biomolecule_id, ome.x, diffexp.x, diffexp.y))


## SingleBiomoleculePlots-AC->LC->noLC ----

comparison_df_pl <- filtered_df_pl %>%
  filter(PASCnoPASC %in% c("No_COVID", "Acute_COVID", "Long_COVID")) %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  filter(PG_change_collection_cutoff == 0)

for (i in unique(comparison_df_pl$standardized_name)) {
  
  poi <- i
  
  single_feat_df <- comparison_df_pl %>%
    filter(standardized_name == poi) %>%
    mutate(Cohort = factor(PASCnoPASC, levels = c("No_COVID", "Acute_COVID", "Long_COVID")))
  
  om <- unique(single_feat_df$ome)
  
  ggplot(single_feat_df, aes(PASCnoPASC, normalized_abundance, fill = PASCnoPASC, color = PASCnoPASC)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
    scale_color_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = NULL,
         y = "Log2 Abundance") +
    scale_y_continuous(expand = c(0,0), limits = c(min(single_feat_df$normalized_abundance) / 1.1, 
                                                   max(single_feat_df$normalized_abundance) * 1.1)) +
    theme_classic() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size = 7),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none"
    )
  
  ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AC_LC_noLC_comparison_PreTube_singlebiomolecule_', poi, '_distribution_PASCnoPASC_group7.pdf'), 
         width = 6, height = 4, units = "cm")
}


## Figure 3, 4 boxplots ----
bimol <- c("Q86TL0", "P35626", "Q8WVC0", "Q8N7H5")


for (i in bimol) {
  
  poi <- i
  
  single_feat_df <- comparison_df_pl %>%
    filter(standardized_name == poi) %>%
    mutate(Cohort = factor(PASCnoPASC, levels = c("No_COVID", "Acute_COVID", "Long_COVID")))
  
  om <- unique(single_feat_df$ome)
  
  x_p <- ggplot(single_feat_df, aes(PASCnoPASC, normalized_abundance, fill = PASCnoPASC, color = PASCnoPASC)) + 
    geom_jitter(alpha = 0.5, 
                width = 0.1, 
                size = 0.2) +
    geom_boxplot(width = 0.4, 
                 alpha = 0.25, 
                 outliers = F,
                 size = 0.2) +
    scale_fill_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
    scale_color_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
    ggtitle(paste(poi, "Abundance")) +
    labs(x = NULL,
         y = "Log2 Abundance") +
    scale_y_continuous(expand = c(0,0), limits = c(min(single_feat_df$normalized_abundance) / 1.1, 
                                                   max(single_feat_df$normalized_abundance) * 1.1)) +
    theme_classic() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size = 6),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none"
    )
  
  name <- paste0(i, "_p")       
  assign(name, x_p)
  
}

Q86TL0_p + P35626_p + Q8WVC0_p + Q8N7H5_p + plot_layout(ncol = 4)

ggsave(paste0('reports/figures/SingleProteinPlots/', om, '_AC_LC_noLC_comparison_PreTube_singlebiomolecule_', poi, '_distribution_figure3_v2.pdf'), 
       width = 18, height = 4, units = "cm")


#### Acute -> PASC -> Acute_fu -> Healthy Differences ----
##* ANOVA ----
# start with filtered_df_pl for only proteins and lipids (due to acute), also no batch 1

# means
APAH_means <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id, Cohort) %>%
  summarise(mean_value = mean(normalized_abundance), .groups = "drop")

# anova and post-hoc tukey hsd for each biomolecule/comparison
APAH_comparison <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  do({
    model <- aov(normalized_abundance ~ Cohort, data = .)
    tidy_result <- tidy(model)
    tukey <- TukeyHSD(model)
    tukey_df <- as.data.frame(tukey$Cohort)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$biomolecule_id <- unique(.$biomolecule_id)
    tukey_df
  }) %>%
  ungroup()


anova_results <- filtered_df_pl %>%
  filter(Cohort %in% c("Acute", "PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  group_by(biomolecule_id) %>%
  do(tidy(aov(normalized_abundance ~ Cohort, data = .))) %>%
  filter(term == "Cohort") 



# find the proportion of significant biomolecules for each comparison

APAH_signif <- APAH_comparison %>%
  group_by(comparison) %>%
  summarise(
    n_significant = sum(`p adj` < 0.05, na.rm = TRUE),
    n_total = n(),
    proportion_significant = n_significant / n_total
  ) %>%
  arrange(desc(proportion_significant))


##* Healthy vs Acute_fu ttest----





















#### Rank-based AC, no LC, LC comparison ----
# take the 522 (or however) differentially expressed proteins and still rank them by -log10(p)*effectsize
# use AC vs no LC and LC vs no LC and add ranks together to get total rank

UP_rank <- UP_df %>%
  mutate(rank = ((effect_size.x * neglogpvalue.x) + (effect_size.y * neglogpvalue.y)))
  
DOWN_rank <- DOWN_df %>%
  mutate(rank = ((effect_size.x * neglogpvalue.x) + (effect_size.y * neglogpvalue.y)))

UP_DOWN_rank <- UP_rank %>%
  select(biomolecule_id, ome.x, diffexp.x, rank) %>%
  bind_rows(DOWN_rank %>%
              select(biomolecule_id, ome.x, diffexp.x, rank))

UP_DOWN_list <- setNames(UP_DOWN_rank$rank,
                         UP_DOWN_rank$biomolecule_id)


# do enrichment
fgsea <- fgsea(pathways = GO_term_list, 
               stats    = UP_DOWN_list,
               minSize  = 5,
               maxSize  = 1000) %>%
  left_join(GO_terms %>%
              select(GO_term, name),
            by = c("pathway" = "GO_term"))

fwrite(fgsea, paste0("data/processed/fgsea_AnalysisGroup7_AC_noLC_LC_allGOterms.csv"))


plo <- fgsea %>%
  filter(padj < 0.05) %>%
  mutate(signif = case_when(
    padj <= 0.05 & padj > 0.01 ~ "*",
    padj <= 0.01 & padj > 0.0001 ~ "**",
    padj <= 0.0001 & padj > 0.000001 ~ "***",
    padj <= 0.000001 ~ "****"
  ))

 # exclude small sets wholly encompassed by a larger leading edge
 # 1: Convert each string to a list of values
 leading_lists <- str_split(plo$leadingEdge, "\\|")
 
 # 2: Convert to a list of sets for easy comparison
 leading_sets <- map(leading_lists, ~ unique(.x))  

 # Step 3: Create a logical vector marking which rows are *not* subsets of any other row
 keep <- sapply(seq_along(leading_sets), function(i) {
   current <- leading_sets[[i]]
   others <- leading_sets[-i]
   is_subset <- sapply(others, function(x) all(current %in% x))
  !any(is_subset)
})

 #Step 4: Filter the original dataframe
filtered_plo <- plo[keep, ]


# Now do manual checking to exclude redundant go terms
stay <- c("GO:0060589", "GO:0005798", "GO:0004674", "GO:0005096", "GO:0032271", "GO:0019902", 
          "GO:0003676", "GO:0006397", "GO:0006955", "GO:1990904", "GO:0051093"
          )

filtered_plo <- filtered_plo %>%
  filter(pathway %in% stay)



ggplot(filtered_plo, aes(x = NES, y = reorder(name, NES, decreasing = T))) +
  geom_col(size = 1.5,
           alpha = 0.8,
           orientation = "y",
           fill = "lightgray") +
  #scale_fill_manual(values = c(col[2], col[1], col[3])) +
  geom_text(aes(label = signif, x = -0.1), hjust = 1.2, size = 3) +
  labs(x = "NES", 
       y = NULL) +
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
        legend.position = "none",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm"))
ggsave("reports/figures/AnalysisGroup2_AC_noLC_LC_fgseaenrichment_barplot_filtered_group7.pdf", 
       width = 6, height = 6, units = "cm")




#### figure:overlap volcanos with pca ----

ac_p + lc_p + p

ggsave(paste0("reports/figures/ProteinLipid_LRM_ALL_Volcanos_andPCA_small_sharedPASC_nobatch1_group7.pdf"), 
       width = 7, height = 2, units = "in")






fdsa <- metadata %>%
  filter(Cohort != "Acute_NC",
         Cohort != "PASC_fu")
length(unique(fdsa$unique_patient_id))



rawfiles_l <- rawfiles %>%
  filter(run_type == "Sample") %>%
  filter(!grepl("NC", Sample)) %>%
  filter(keep == "1") %>%
  filter(!grepl("\\.2", Sample))
table(rawfiles_l$ome_id)
lipid_sampl <- rawfiles_l %>%
  filter(ome_id == 2) %>%
  pull(sample_id)


metad_lip <- metadata %>%
  filter(sample_id %in% lipid_sampl)
length(unique(metad_lip$unique_patient_id))


transc_sampl <- unique(transcriptomics$sample_id)
metad_transc <- metadata %>%
  filter(sample_id %in% transc_sampl)
length(unique(metad_transc$unique_patient_id))
table(metad_transc$Cohort)
table(metadata$Cohort)











#### Supplemental Table 1 - Overlapping biomolecules (w/ fold changes) ####

UP_df <- volc_plot_PASC %>%
  filter(diffexp == "UP") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "UP"),
             by = "biomolecule_id")

DOWN_df <- volc_plot_PASC %>%
  filter(diffexp == "DOWN") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "DOWN"),
             by = "biomolecule_id")

SHARED_bmols_df <- UP_df %>%
  select(biomolecule_id, standardized_name.x, comparison.x, effect_size.x, q_value.x,
         comparison.y, effect_size.y, q_value.y,
         ome.x, diffexp.x) %>%
  bind_rows(DOWN_df %>%
              select(biomolecule_id, standardized_name.x, comparison.x, effect_size.x, q_value.x,
                     comparison.y, effect_size.y, q_value.y,
                     ome.x, diffexp.x)) %>%
  rename(standardized_name = standardized_name.x,
         ome = ome.x,
         diffexp = diffexp.x) %>%
  mutate(comparison.x = if_else(comparison.x == "group7_PASCnoPASC", "NoCOVID_LongCOVID", comparison.x)) %>%
  mutate(comparison.y = if_else(comparison.y == "Acute_nonAcute", "NoCOVID_AcuteCOVID", comparison.y))

write_csv(SHARED_bmols_df, "data/temp/TableS1_SharedChangedBiomolecules.csv")



















