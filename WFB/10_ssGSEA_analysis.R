#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(RSQLite)
library(ggrepel)

library(GSVA)
library(mogsa)
library(msigdbr)
library(UniProt.ws)
library(pheatmap)


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

proteomics <- dbGetQuery(con, 'SELECT *
                               FROM proteomics_measurement')
lipidomics <- dbGetQuery(con, 'SELECT *
                              FROM lipidomics_measurements')
transcriptomics <- dbGetQuery(con, "SELECT *
                                    FROM rnaseq_measurements")
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_id, rawfile_name, Sample, sample_id, run_type, ome_id, keep, batch
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')

dbDisconnect(con)



## Make proteomics/lipidomics/transcriptomics dataframes ----



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



## GSVA ssGSEA organization ----
## * find gene sets ----
## ** MsigDB gene set ----
m_df <- msigdbr(species = "Homo sapiens",
                collection = "H")
m_gs <- m_df %>%
  split(x = .$gene_symbol, f = .$gs_name)


## * Proteomics ----
## organize into wide format, can choose subsets of samples, or just do all
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
  filter(biomolecule_id %in% biomolecules_p) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  filter(Cohort %in% c("Acute", "PASC", "Healthy", "Acute_fu"))
length(unique(filtered_df_p$sample_id))

filtered_df_p_ssgsea <- filtered_df_p %>%
  left_join(biomolecules %>%
              select(biomolecule_id, standardized_name),
            by = "biomolecule_id") %>%
  select(normalized_abundance, sample_id, standardized_name) %>%
  pivot_wider(names_from = sample_id,
              values_from = normalized_abundance)

# map uniprot ids to gene symbols
filtered_df_p_ssgsea_symbols <- mapUniProt(
  from = "UniProtKB_AC-ID",
  to = "Gene_Name",
  query = filtered_df_p_ssgsea$standardized_name)

# switch the standardized_name column to the gene symbols
filtered_df_p_ssgsea <- filtered_df_p_ssgsea %>%
  separate_rows(standardized_name, sep = ";") %>%
  left_join(filtered_df_p_ssgsea_symbols %>%
              rename(gene_name = To,
                     standardized_name = From),
            by = "standardized_name") %>%
  select(-standardized_name) %>%
  filter(!is.na(gene_name)) %>%
  relocate(gene_name, .before = everything()) %>%
  group_by(gene_name) %>%
  filter(n() < 2) %>% # get rid of the ones with multiple mappings
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()
  
## map ssGSEA for proteomics ----
gsvaPar <- gsvaParam(filtered_df_p_ssgsea,
                     m_gs)

ss_scores <- gsva(gsvaPar, verbose=FALSE) %>%
  as.data.frame()


## heatmap
sample_annot <- filtered_df_p %>%
  filter(!is.na(Cohort)) %>%
  select(sample_id, Cohort, Sex, Age) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id")

ss_heatmap <- pheatmap(ss_scores,
                       annotation_col = sample_annot)

## heatmap, average by cohort for samples
ss_scores_cohort_average <- ss_scores %>%
  rownames_to_column(var = "gene_set") %>%
  pivot_longer(cols = -gene_set,
               names_to = "sample_id",
               values_to = "ssgsea_abundance") %>%
  mutate(sample_id = as.integer(sample_id)) %>%
  left_join(metadata %>%
              select(sample_id, Cohort),
            by = "sample_id") %>%
  group_by(Cohort, gene_set) %>%
  summarize(avg_ssgsea_zscore = mean(ssgsea_abundance)) %>%
  pivot_wider(names_from = Cohort,
              values_from = avg_ssgsea_zscore) %>%
  column_to_rownames(var = "gene_set")

ss_heatmap_cohort_average <- pheatmap(ss_scores_cohort_average)
