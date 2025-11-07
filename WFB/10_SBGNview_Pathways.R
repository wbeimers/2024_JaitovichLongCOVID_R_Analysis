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



# ## Make proteomics/lipidomics/transcriptomics dataframes ----
# 
# # proteomics
# # Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# # subset rawfiles to only include sample proteomics runs
# rawfiles_p <- rawfiles %>%
#   filter(ome_id == 1) %>%
#   select(-keep) %>%
#   filter(grepl('Sample', run_type))
# 
# df_p <- proteomics %>%
#   inner_join(rawfiles_p, by = 'rawfile_id') %>%
#   inner_join(metadata, by = 'sample_id')
# 
# biomolecules_p <- biomolecules %>%
#   filter(omics_id == 1) %>%
#   filter(keep == "1") %>%
#   pull(biomolecule_id)
# 
# filtered_df_p <- df_p %>%
#   filter(biomolecule_id %in% biomolecules_p)
# 
# 
# # lipidomics
# rawfiles_l <- rawfiles %>%
#   filter(ome_id == 2) %>%
#   filter(keep == "1") %>%
#   filter(run_type == "Sample") %>%
#   select(-keep)
# 
# df_l <- lipidomics %>%
#   inner_join(rawfiles_l, by = 'rawfile_id') %>%
#   inner_join(metadata, by = 'sample_id')
# 
# biomolecules_l <- biomolecules %>%
#   filter(omics_id == 2) %>%
#   filter(keep == "1") %>%
#   pull(biomolecule_id)
# 
# filtered_df_l <- df_l %>%
#   filter(biomolecule_id %in% biomolecules_l)
# 
# 
# # combined lipid and protein (no transcript in Acute samples)
# df_pl <- df_p %>%
#   mutate(ome = "p") %>%
#   bind_rows(df_l %>%
#               mutate(ome = "l"))
# 
# filtered_df_pl <- filtered_df_p %>%
#   mutate(ome = "p") %>%
#   bind_rows(filtered_df_l %>%
#               mutate(ome = "l")) %>%
#   filter(batch != 1) %>% # remove samples from Batch 1
#   mutate(PASCnoPASC = case_when(
#     Cohort %in% c("Acute") ~ "Acute_COVID",
#     Cohort %in% c("Healthy", "Acute_fu") ~ "No_COVID",
#     Cohort %in% c("PASC") ~ "Long_COVID"))




## read in overlapping biomolecules data ----
shared_bmols_df <- read_csv("data/temp/TableS1_SharedChangedBiomolecules.csv")



## Export to Reactome for Pathway Analysis ----
## Choose only significant proteins from pasc vs no pasc comparison and export as uniprot id list
protein_list <- shared_bmols_df %>%
  filter(ome == "protein") %>%
  separate_rows(standardized_name, sep = ";") %>%
  pull(standardized_name)
write_lines(protein_list, "data/processed/ACUTEvsHEALTHYvsPASC_overlap_significant_proteins_forreactome.txt")





## PANTHER protein levels ----
# find a comparison like pasc vs nopasc and take all proteins and go to effect size
pathway_data <- shared_bmols_df %>%
  select(standardized_name, effect_size.y, diffexp, ome) %>%
  filter(ome == "protein") %>%
  separate_rows(standardized_name, sep = ";")


## * make the pathway plot for the pathways ----
library(SBGNview)
data("sbgn.xmls")

data("pathways.info")
pathways <- findPathways(c("rho"))
head(pathways)
pathway <- pathways$pathway.id[16]

# translate uniprot to other id and take effect size of proteins for your pathway
nodes <- pathway_data %>%
  select(standardized_name, effect_size.y) %>%
  column_to_rownames(var = "standardized_name") 


# for protein!!
gene.data.pc <- changeDataId(
  data.input.id = nodes,
  input.type = "UNIPROT",
  output.type = "entrez",
  mol.type = "gene"
)




SBGNview.obj <- SBGNview(
  gene.data = gene.data.pc, 
  gene.id.type = "entrez",
  input.sbgn = pathway,
  output.file = "reactome_rhogtpaseactivatewave_PROTEINS_acutevshealthyvspascoverlap", 
  output.formats = c("png")
) 

print(SBGNview.obj)



