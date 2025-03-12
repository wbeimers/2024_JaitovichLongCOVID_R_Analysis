
## Proteomics Biomolecule Metadata ----
# 1. read in metadata
# 2. remove biomolecule id to update
# 3. duplicate rows for NPS






# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
dbDisconnect(con)

biomolecule_metadata <- fread("data/metadata/protein_uniprot_metadata.csv") %>%
  dplyr::select(-biomolecule_id)

biom_df <- biomolecule_metadata %>%
  left_join(biomolecules, by = "standardized_name") %>%
  filter(keep == "1")


con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
dbWriteTable(con, "proteomics_metadata", biomolecule_metadata, overwrite = T)
dbDisconnect(con)


