#### Editing the metadata tables #### 


library(DBI)
library(RSQLite)


#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### List tables in db #####

dbGetInfo(con)
dbListTables(con)
dbReadTable(con, "biomolecules")


#### Fetch biomolecules table and biomolecules metadata tables, to be combined ####

biomolecules <- dbReadTable(con, "biomolecules")
#lipidomics_metadata <- dbReadTable(con, "lipidomics_metadata")
proteomics_metadata <- dbReadTable(con, "proteomics_metadata")
transcriptomics_metadata <- dbReadTable(con, "rnaseq_metadata")
transcriptomics_metadata2 <- dbReadTable(con, "rnaseq_measurements")

dbDisconnect(con)

lipidomics_metadata <- read.csv("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/db_formatted_tables/biomolecules_metadata_v2.csv")

#### Create correct format for biomolecules_metadata table ####
## Biomolecules metadata should contain following headers: metadata_id, biomolecule_id,
## metadata_type, and metdata_value.   
names(lipidomics_metadata)

lipid_metadata_forDB <- rbind(data.frame(metadata_id = NA, 
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "lipid_category", 
                                         metadata_value = lipidomics_metadata$LipidCategory),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "lipd_class", 
                                         metadata_value = lipidomics_metadata$LipidClass),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "RT", 
                                         metadata_value = lipidomics_metadata$RT),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "m/z", 
                                         metadata_value = lipidomics_metadata$mz),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "adduct", 
                                         metadata_value = lipidomics_metadata$Adduct),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "main_class", 
                                         metadata_value = lipidomics_metadata$MainClass),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "sub_class", 
                                         metadata_value = lipidomics_metadata$SubClass),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "number_of_carbons", 
                                         metadata_value = lipidomics_metadata$NumFattyAcylCarbons),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "number_of_unsaturations", 
                                         metadata_value = lipidomics_metadata$NumFattyAcylUnsaturations),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "average_unsaturations", 
                                         metadata_value = lipidomics_metadata$AverageUnsaturations),
                              data.frame(metadata_id = NA,
                                         biomolecule_id = lipidomics_metadata$biomolecule_id, 
                                         metadata_type = "unsaturation_level", 
                                         metadata_value = lipidomics_metadata$UnsaturationLevel))


proteomics_metadata <- merge(biomolecules, proteomics_metadata, by = 'standardized_name')

names(proteomics_metadata)

protein_metadata_forDB <- rbind(data.frame(metadata_id = NA, 
                                           biomolecule_id = proteomics_metadata$biomolecule_id, 
                                           metadata_type = "gene_symbol", 
                                           metadata_value = proteomics_metadata$gene_name),
                                data.frame(metadata_id = NA,
                                           biomolecule_id = proteomics_metadata$biomolecule_id, 
                                           metadata_type = "GO_terms", 
                                           metadata_value = proteomics_metadata$GO_terms),
                                data.frame(metadata_id = NA,
                                           biomolecule_id = proteomics_metadata$biomolecule_id, 
                                           metadata_type = "entry_name", 
                                           metadata_value = proteomics_metadata$Entry.Name),
                                data.frame(metadata_id = NA,
                                           biomolecule_id = proteomics_metadata$biomolecule_id, 
                                           metadata_type = "protein_names", 
                                           metadata_value = proteomics_metadata$Protein.names))

names(transcriptomics_metadata)
names(transcriptomics_metadata2)

RNA_metadata <- merge(transcriptomics_metadata, transcriptomics_metadata2[,c(1:3,6)], by = 'standardized_name')

head(RNA_metadata)
RNA_metadata<- unique(RNA_metadata)
RNA_metadata<- RNA_metadata[RNA_metadata$biomolecule_id.x == RNA_metadata$biomolecule_id.y,]


RNA_metadata_forDB <- rbind(data.frame(metadata_id = NA, 
                                       biomolecule_id = RNA_metadata$biomolecule_id.x, 
                                       metadata_type = "gene_symbol", 
                                       metadata_value = RNA_metadata$SYMBOL),
                            data.frame(metadata_id = NA,
                                       biomolecule_id = RNA_metadata$biomolecule_id.x, 
                                       metadata_type = "GO_terms", 
                                       metadata_value = RNA_metadata$GO_terms),
                            data.frame(metadata_id = NA,
                                       biomolecule_id = RNA_metadata$biomolecule_id.x, 
                                       metadata_type = "gene_name", 
                                       metadata_value = RNA_metadata$GENENAME),
                            data.frame(metadata_id = NA,
                                       biomolecule_id = RNA_metadata$biomolecule_id.x, 
                                       metadata_type = "entrez_id", 
                                       metadata_value = RNA_metadata$standardized_name))


## combine all metadata 

biomolecules_metadata <- rbind(protein_metadata_forDB, RNA_metadata_forDB, lipid_metadata_forDB)

biomolecules_metadata$metadata_id <- seq(1:nrow(biomolecules_metadata))

##### Add biomolecules_metadata to Db and remove other metadata tables ####


## Establish connection to SQLite db 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

## Add biomolecules_metadata 

dbWriteTable(con, "biomolecules_metadata", biomolecules_metadata, overwrite = T)

dbGetQuery(con, 'PRAGMA table_info(biomolecules_metadata)')

## remove old metadata table 
dbRemoveTable(con, "lipidomics_metadata")
dbRemoveTable(con, "proteomics_metadata")
dbRemoveTable(con, "rnaseq_metadata")
dbRemoveTable(con, "biomolecule_metadata")

## check 

dbListTables(con)

## disconnect
dbDisconnect(con)

## Change made 3/17/2025 at 2:15pm 