library(DBI)
library(RSQLite)
library(pheatmap)

colors2 <- c("#D66127", "#199D77", "#7670B2", "#FBB469")

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### pvalues and biomolecule data from db #####

dbGetInfo(con)
dbListTables(con)

#### Fetch biomolecules table and biomolecules_metadata tables ####

pvalues_PASC_noBatch1 <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, effect_size, eta_squared, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '35'
           ")

biomolecules_metadata <- dbReadTable(con, "biomolecules_metadata")

dbDisconnect(con)
#### number of significant features ####

sig_by_ome <- table(sig = pvalues_PASC_noBatch1$q_value < 0.05, direction = pvalues_PASC_noBatch1$effect_size >0, pvalues_PASC_noBatch1$omics_id)

str(sig_by_ome)

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_biomolecules_sig_withPASC_horiz.pdf", height = 7, width = 9)

x <- barplot(cbind(c(sig_by_ome[2,1,1], sig_by_ome[1,1,1], sig_by_ome[1,2,1], sig_by_ome[2,2,1]),
              c(sig_by_ome[2,1,2], sig_by_ome[1,1,2], sig_by_ome[1,2,2], sig_by_ome[2,2,2]), 
              c(sig_by_ome[2,1,3], sig_by_ome[1,1,3], sig_by_ome[1,2,3], sig_by_ome[2,2,3])),
        col = c("blue", "lightblue", "pink", "red"), las = 1, xlim = c(-10000, 10000), horiz = T)
text( c(rep(8000, 3), rep(-8000, 3)),c(x,x), c(sig_by_ome[2,2,1:3], sig_by_ome[2,1,1:3]))
dev.off()

#### Plot Volcano plots for each ome #### 

##### Lipids #####
pvalues_PASC_noBatch1$omics_id == 2 
names(pvalues_PASC_noBatch1)


table(biomolecules_metadata$metadata_type)
table(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "Lipid_Class"])
sum(table(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "Lipid_Class"]))

as.numeric(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "NumFattyAcylCarbons"]) %% 2 == 1

## SMs are enriched in PASC
SMs <- biomolecules_metadata$biomolecule_id[which(biomolecules_metadata$metadata_value == "SM" & biomolecules_metadata$metadata_type == "Lipid_Class")]
SMs_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in% SMs

PCs <- biomolecules_metadata$biomolecule_id[which(biomolecules_metadata$metadata_value == "PC" & biomolecules_metadata$metadata_type == "Lipid_Class")]
Odd_chain <- biomolecules_metadata$biomolecule_id[which(biomolecules_metadata$metadata_type == "NumFattyAcylCarbons")] 
Odd_chain <- Odd_chain[which(as.numeric(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "NumFattyAcylCarbons"]) %% 2 == 1)] 

Odd_chain[Odd_chain %in% PCs]

Odd_chain_PCs_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in% Odd_chain[Odd_chain %in% PCs]
Odd_chain_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in% Odd_chain


Cer <- biomolecules_metadata$biomolecule_id[which(biomolecules_metadata$metadata_value == "Cer[NS]" & biomolecules_metadata$metadata_type == "Lipid_Class")]
Cer_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in% Cer


## plotting limits 
max(pvalues_PASC_noBatch1$effect_size) #1.76
min(pvalues_PASC_noBatch1$effect_size) #-1.84
max(-log(pvalues_PASC_noBatch1$q_value)) #45.3

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_lipids_sig_withPASC.pdf", height = 7, width = 9)
plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2, names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2, names(pvalues_PASC_noBatch1) == "q_value"]), 
    pch = 19, col = colors2[1], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

text(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & -log(pvalues_PASC_noBatch1$q_value) >4, names(pvalues_PASC_noBatch1)== "effect_size"],
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & -log(pvalues_PASC_noBatch1$q_value) >4, names(pvalues_PASC_noBatch1)== "q_value"]),
      pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & -log(pvalues_PASC_noBatch1$q_value) >4, names(pvalues_PASC_noBatch1)== "standardized_name"])
points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & SMs_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & SMs_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
     pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")
dev.off()



pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_lipids_sig_withPASC_odd_chain.pdf", height = 7, width = 9)


plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2, names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2, names(pvalues_PASC_noBatch1) == "q_value"]), 
     pch = 19, col = colors2[1], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & Odd_chain_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
       -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & Odd_chain_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
       pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

legend("topright", c("lipids", "odd-chain lipids"), col =c(colors2[1], "black"), pch = 19)

dev.off()



pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_lipids_sig_withPASC_odd_chain_noUnknowns.pdf", height = 7, width = 9)

plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & !grepl("Unknown", pvalues_PASC_noBatch1$standardized_name), names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & !grepl("Unknown", pvalues_PASC_noBatch1$standardized_name), names(pvalues_PASC_noBatch1) == "q_value"]), 
     pch = 19, col = colors2[1], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & Odd_chain_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
       -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & Odd_chain_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
       pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

legend("topright", c("lipids", "odd-chain lipids"), col =c(colors2[1], "black"), pch = 19)
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_lipids_sig_withPASC_SMs_noUnknowns.pdf", height = 7, width = 9)

plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & !grepl("Unknown", pvalues_PASC_noBatch1$standardized_name), names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & !grepl("Unknown", pvalues_PASC_noBatch1$standardized_name), names(pvalues_PASC_noBatch1) == "q_value"]), 
     pch = 19, col = colors2[1], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & SMs_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
       -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 2 & SMs_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
       pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")

legend("topright", c("lipids", "SM lipids"), col =c(colors2[1], "black"), pch = 19)
dev.off()



##### proteins ##### 

#### 
GO <- "GO:0015629"
GO <- "GO:0046983"

Plot_GO  <- function(GO, type = "proteins"){
  #table(biomolecules_metadata$metadata_type)
  #table(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "GO_terms"])
  
  biomolecules_GO <- biomolecules_metadata$biomolecule_id[grep(GO, biomolecules_metadata$metadata_value)]
 
if(type == "proteins"){
  GO_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in% biomolecules_GO 
  
  plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1, names(pvalues_PASC_noBatch1) == "effect_size"], 
       -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1, names(pvalues_PASC_noBatch1) == "q_value"]), 
       pch = 19, col = colors2[2], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")
  
  points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & GO_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
         -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & GO_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
         pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")
  
  legend("topright", c("proteins", paste(GO)), col =c(colors2[2], "black"), pch = 19)  
}    
  
  
  if(type == "transcripts"){
    biomolecules_GO <- biomolecules_metadata$biomolecule_id[grep(GO, biomolecules_metadata$metadata_value)]
    GO_pvalues <- pvalues_PASC_noBatch1$biomolecule_id %in%   biomolecules_GO
    
    plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 3, names(pvalues_PASC_noBatch1) == "effect_size"], 
         -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 3, names(pvalues_PASC_noBatch1) == "q_value"]), 
         pch = 19, col = colors2[3], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")
    
    points(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 3 & GO_pvalues, names(pvalues_PASC_noBatch1) == "effect_size"], 
           -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 3 & GO_pvalues, names(pvalues_PASC_noBatch1) == "q_value"]), 
           pch = 19, col = "black", xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")
    
    legend("topright", c("transcripts", paste(GO)), col =c(colors2[3], "black"), pch = 19)  
  } 
  
}



GeneName_proteins <- biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "gene_name"][match(pvalues_PASC_noBatch1$biomolecule_id, biomolecules_metadata$biomolecule_id[biomolecules_metadata$metadata_type == "gene_name"])]

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC.pdf", height = 7, width = 9)

plot(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1, names(pvalues_PASC_noBatch1) == "effect_size"], 
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1, names(pvalues_PASC_noBatch1) == "q_value"]), 
     pch = 19, col = colors2[2], xlim = c(-2,2), las = 1, ylab = "-log(q-value)", xlab = "effect size", bty = "l")


text(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) > 5 & pvalues_PASC_noBatch1$effect_size> 0.6, names(pvalues_PASC_noBatch1)== "effect_size"],
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) >5 & pvalues_PASC_noBatch1$effect_size> 0.6, names(pvalues_PASC_noBatch1)== "q_value"]),
     GeneName_proteins[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) >5 & pvalues_PASC_noBatch1$effect_size> 0.6])


text(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) > 5 & pvalues_PASC_noBatch1$effect_size <  -0.6, names(pvalues_PASC_noBatch1)== "effect_size"],
     -log(pvalues_PASC_noBatch1[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) >5 & pvalues_PASC_noBatch1$effect_size < -0.6, names(pvalues_PASC_noBatch1)== "q_value"]),
     GeneName_proteins[pvalues_PASC_noBatch1$omics_id == 1 & -log(pvalues_PASC_noBatch1$q_value) >5 & pvalues_PASC_noBatch1$effect_size < -0.6])

dev.off()



pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC_GOlipid_transport.pdf", height = 7, width = 9)

# lipid transport
Plot_GO("GO:0006869", "proteins")
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC_GOSphingolipid_delta4_desaturase_activity.pdf", height = 7, width = 9)

# sphingolipid delta-4 desaturase activity 
Plot_GO("GO:0042284", "proteins")
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC_GOhumoral_immune_response.pdf", height = 7, width = 9)

# humoral immune response
Plot_GO("GO:0006959", "proteins")

dev.off()


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC_GOglycosaminoglycan_binding.pdf", height = 7, width = 9)

# glycosaminoglycan binding
Plot_GO("GO:0005539", "proteins")

dev.off()


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/09_proteins_sig_withPASC_GOglycosaminoglycan_binding.pdf", height = 7, width = 9)
#  
Plot_GO("GO:0030545", "proteins")

dev.off()



#  
Plot_GO("GO:0030218", "transcripts")

length(unique(biomolecules_metadata$metadata_value[biomolecules_metadata$metadata_type == "Protein_names"]))
