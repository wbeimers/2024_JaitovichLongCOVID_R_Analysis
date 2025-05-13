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
  mutate(measurement_id = row_number()) %>%
  rename(standardized_name = ENTREZID) %>%
  mutate(standardized_name = as.character(standardized_name)) %>%
  left_join(biomolecules %>% 
              filter(omics_id == 3) %>%
              select(standardized_name, biomolecule_id),
            by = "standardized_name") %>%
  left_join(rawfiles %>% 
              select(Sample, sample_id) %>%
              distinct(),
            by = "Sample") %>%
  mutate(normalized_counts = log2(Counts)) %>%
  select(-Sample) %>%
  select(measurement_id, standardized_name, biomolecule_id, SYMBOL, GENENAME, sample_id, Counts, normalized_counts)


hist(rnaseq_long %>%
       filter(biomolecule_id == 25000) %>%
       pull(normalized_counts))



con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

# drop old patient metadata table before adding new one
dbExecute(con, "DROP TABLE IF EXISTS rnaseq_measurements")

# Write the new patient_metadata table to the database
dbWriteTable(con, "rnaseq_measurements", rnaseq_long, append = F, overwrite = T)

dbDisconnect(con)







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


## RNAseq counts filtering ----
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
transcriptomics <- dbGetQuery(con, 'SELECT *
                                    FROM rnaseq_measurements')
dbDisconnect(con)


# Also filter for completeness
# Show how many counts<10 there are for each gene
counts_summary <- transcriptomics %>%
  group_by(biomolecule_id) %>%
  summarise(countsunder10 = mean(Counts < 10), .groups = 'drop')

hist(counts_summary$countsunder10)

# Make a list of IDs to keep where there are at least 25% samples counts < 10
ids_to_exclude <- counts_summary %>%
  group_by(biomolecule_id) %>%
  filter(countsunder10 >= 0.75) %>% 
  pull(biomolecule_id)

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
biomolecules <- dbGetQuery(con, 'SELECT *
                                    FROM biomolecules')
dbDisconnect(con)

biomolecules1 <- biomolecules %>%
  mutate(keep = case_when(
    biomolecule_id %in% ids_to_exclude ~ "0",
    T ~ keep
  ))

asdfasdfs <- biomolecules1 %>%
  filter(omics_id == 3)

table(asdfasdfs$keep)


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

# drop old patient metadata table before adding new one
dbExecute(con, "DROP TABLE IF EXISTS biomolecules")

# Write the new patient_metadata table to the database
dbWriteTable(con, "biomolecules", biomolecules1, append = F, overwrite = T)

dbDisconnect(con)




## RNAseq normalized_counts Inf to 0 ----
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")
transcriptomics <- dbGetQuery(con, 'SELECT *
                                    FROM rnaseq_measurements')
dbDisconnect(con)

transcriptomics1 <- transcriptomics %>% 
  mutate(normalized_counts = case_when(
    is.infinite(normalized_counts) ~ 0,
    T ~ normalized_counts
  )) # replace infinite values (0 counts) with 0



con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

# drop old patient metadata table before adding new one
dbExecute(con, "DROP TABLE IF EXISTS rnaseq_measurements")

# Write the new patient_metadata table to the database
dbWriteTable(con, "rnaseq_measurements", transcriptomics1, append = F, overwrite = T)

dbDisconnect(con)




## Drop extra columns ----

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbExecute(con, "ALTER TABLE rnaseq_measurements DROP COLUMN SYMBOL")
dbExecute(con, "ALTER TABLE rnaseq_measurements DROP COLUMN GENENAME")

dbDisconnect(con)

