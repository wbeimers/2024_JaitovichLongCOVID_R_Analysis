#### make a transcriptomics_measurements table for the sql database ----

library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)
library(data.table)


## Add RNAseq data ----
# read in report
rnaseq <- fread("data/RNA_seq/For Coon Group RPM normalized read counts PASC RNA-Seq.csv",
                data.table = F)

rnaseq_long <- rnaseq %>%
  dplyr::select(-V1) %>%
  pivot_longer(cols = -c("ENTREZID", "SYMBOL", "GENENAME"),
               names_to = "Sample",
               values_to = "Counts") %>%
  mutate(Sample = gsub("^X|\\.BAM$", "", Sample)) %>%
  mutate(Sample = case_when(
    Sample == 'PC2001' ~ '2001.58',
    Sample == 'PC2002' ~ '2002.33',
    Sample == 'PC2003' ~ '2003.16',
    Sample == 'PC2004' ~ '2004.51',
    Sample == 'PC2005' ~ '2005.66',
    Sample == 'PC2006' ~ '2006.95',
    Sample == 'PC2007' ~ '2007.57',
    Sample == 'PC2008' ~ '2008.53',
    Sample == 'PC2009' ~ '2009.56',
    Sample == 'PC2010' ~ '2010.11',
    Sample == 'PC2011' ~ '2011.101',
    Sample == 'PC2012' ~ '2012.7',
    Sample == 'PC2013' ~ '2013.24',
    Sample == 'PC2015' ~ '2015.77',
    Sample == 'PC2017' ~ '2017.9',
    Sample == 'PC2018' ~ '2018.67',
    Sample == 'PC2019' ~ '2019.61',
    T ~ Sample)) %>%
  mutate(biomolecule_id = 9595 + cumsum(!duplicated(ENTREZID))) %>%
  mutate(measurement_id = row_number()) %>%
  rename(standardized_name = ENTREZID)


# Add rnaseq table to the sql database
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbWriteTable(con, "rnaseq_measurements", rnaseq_long, append = T)

dbDisconnect(con)


## Update Omes table ----
rna_ome <- data.frame(omics_id = 3,
                      omics_name = "Transcriptomics")

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbWriteTable(con, "omes", rna_ome, append = T, row.names = F)

dbDisconnect(con)


## Update biomolecules table ----
# make data frame to append to the end

rnaseq_long_bm <- rnaseq_long %>%
  select(biomolecule_id, standardized_name) %>%
  distinct() %>%
  mutate(omics_id = 3) %>%
  mutate(keep = 1)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbWriteTable(con, "biomolecules", rnaseq_long_bm, append = T, row.names = F)

dbDisconnect(con)



## Put new biomolecule IDs in the rnaseq table

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
transcriptomics <- dbGetQuery(con, 'SELECT standardized_name, SYMBOL, GENENAME, Sample, Counts, measurement_id
                                    FROM rnaseq_measurements')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                                 FROM biomolecules')
dbDisconnect(con)

transcriptomics1 <- transcriptomics %>%
  mutate(standardized_name = as.character(standardized_name)) %>%
  left_join(biomolecules %>% dplyr::select(biomolecule_id, standardized_name), by = "standardized_name")

transcriptomics1 <- transcriptomics1 %>%
  dplyr::select(measurement_id, standardized_name, biomolecule_id, everything())


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "rnaseq_measurements", transcriptomics1, append = T, row.names = F)
dbDisconnect(con)


## RNAseq GO Terms ----

transcript_go <- fread("data/metadata/transcripts_go_terms.csv") %>%
  mutate(standardized_name = as.character(standardized_name))

transcript_go <- biomolecules %>%
  filter(omics_id == 3) %>%
  dplyr::select(biomolecule_id, standardized_name) %>%
  right_join(transcript_go, by = "standardized_name")

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
dbWriteTable(con, "rnaseq_metadata", transcript_go, row.names = F, overwrite = T)
dbDisconnect(con)
