#### Exploring pvalues with and without Batch 1 lipids ####

library(DBI)
library(RSQLite)
library(pheatmap)

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### pvalues and biomolecule data from db #####

dbGetInfo(con)
dbListTables(con)

#### Fetch biomolecules table and biomolecules_metadata tables ####
pvalues_QOL <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, effect_size, eta_squared, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '1'
           ")

pvalues_QOL_noBatch1 <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, effect_size, eta_squared, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '22'
           ")
pvalues_PASC <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, effect_size, eta_squared, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '14'
           ")

pvalues_PASC_noBatch1 <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, effect_size, eta_squared, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '35'
           ")

biomolecules_metadata <- dbReadTable(con, "biomolecules_metadata")

dbDisconnect(con)

#### Plot pvalue correlations across the different tables #### 

plot(-log(pvalues_PASC$q_value)*pvalues_PASC$effect_size, -log(pvalues_PASC_noBatch1$q_value)*pvalues_PASC_noBatch1$effect_size, col = pvalues_PASC$omics_id, pch = 19)
table(pvalues_PASC$q_value < 0.05, pvalues_PASC$omics_id)
table(pvalues_PASC_noBatch1$q_value < 0.05, pvalues_PASC$omics_id)


plot(-log(pvalues_QOL$q_value)*pvalues_QOL$effect_size, -log(pvalues_QOL_noBatch1$q_value)*pvalues_QOL_noBatch1$effect_size, col = pvalues_QOL$omics_id, pch = 19)
table(pvalues_QOL$q_value < 0.05, pvalues_QOL$omics_id)
table(pvalues_QOL_noBatch1$q_value < 0.05, pvalues_QOL$omics_id)

pvalues_PASC$rank <- -log(pvalues_PASC$q_value)*-sign(pvalues_PASC$effect_size)
pvalues_PASC_noBatch1$rank <- -log(pvalues_PASC_noBatch1$q_value)*-sign(pvalues_PASC_noBatch1$effect_size)
pvalues_QOL$rank <- -log(pvalues_QOL$q_value)*sign(pvalues_QOL$effect_size)
pvalues_QOL_noBatch1$rank <- -log(pvalues_QOL_noBatch1$q_value)*sign(pvalues_QOL_noBatch1$effect_size)

plot(pvalues_PASC_noBatch1$rank, pvalues_QOL_noBatch1$rank)
cor(pvalues_QOL_noBatch1$effect_size, -pvalues_PASC_noBatch1$effect_size)
plot(pvalues_QOL_noBatch1$effect_size, -pvalues_PASC_noBatch1$effect_size)

#### 

sig_biomolecules <- pvalues_PASC_noBatch1$biomolecule_id[pvalues_PASC_noBatch1$q_value < 0.05 & abs(pvalues_PASC_noBatch1$effect_size) >0.5]

sig_biomolecules_pvalues_table <- pvalues_PASC_noBatch1[match(sig_biomolecules, pvalues_PASC_noBatch1$biomolecule_id), ]
sig_biomolecules_pvalues_table[order(sig_biomolecules_pvalues_table$rank), ]
x <- biomolecules_metadata[match(sig_biomolecules, biomolecules_metadata$biomolecule_id), ]

#### Looking through the Compound Discoverer viewer to determine if we can figure out more info on lipid unknowns #### 

## Unknown RT 11.62 mz 466.29814 
### this unknown appears to have ok quant

## Unknown RT 12.72 mz 223.09644 
###This unknownh as poor quant. The feature appears to be bimodal - the RT 12.8 is the better quantified feature 
### Recommend setting keep to 0 

## Unknown RT 12.897 mz 263.08887
### This unkown has good quant. No MS/MS data. has same elution profile as mz 223.09637

## Unknown RT 14.60 mz 532.27958
### Quant appears ok. Similar elution profile as mz 516.30559

## Unknown RT 22.86 mz 1596.19541 
### This appears to be an in-source heterodimer. Fragment ions are 810.70 and 786.59 which are the primary MS1 ions at this time point. 
### Recommend setting keep to 0

## Unknown RT 24.035 mz 834.57681 
### Originally identified as PC 40:6, but the the RT is not consistent with this identification (based on degreaser plots)




