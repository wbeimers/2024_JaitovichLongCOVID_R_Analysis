#### Exploring cross-ome correlation matrix ####

library(DBI)
library(RSQLite)
library(pheatmap)

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

#### pvalues and biomolecule data from db #####

dbGetInfo(con)
dbListTables(con)

#### Fetch biomolecules table and biomolecules_metadata tables ####
pvalues <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, lratio, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '1'
           ")

biomolecules_metadata <- dbReadTable(con, "biomolecules_metadata")

dbDisconnect(con)

#### load the cross-ome correlation matrices #### 

lipid_lipid <- read.csv("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/processed_data/BJA/cross_ome_correlations/Kendall-cross-ome-corrs_lipid_lipid.csv")

protein_lipid <- read.csv("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/processed_data/BJA/cross_ome_correlations/Kendall-cross-ome-corrs_protein_lipid.csv")

protein_rna <- read.csv("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/processed_data/BJA/cross_ome_correlations/Kendall-cross-ome-corrs_protein_rnaseq.csv")

#### preliminary plots/exploration of lipid-lipid clusters #### 
str(lipid_lipid)
lipid_lipid[1:10,1:10]

# first column should be the row name
rownames(lipid_lipid)<-lipid_lipid[,1]
lipid_lipid<- as.matrix(lipid_lipid[,-1])

matching_names <- pvalues[match(rownames(lipid_lipid), pvalues$biomolecule_id),]
lipid_class <- biomolecules_metadata[biomolecules_metadata$metadata_type == "Lipid_Class",] 
matching_names <- cbind(matching_names, lipid_class[match(rownames(lipid_lipid),lipid_class$biomolecule_id), ])
matching_names$metadata_value[is.na(matching_names$metadata_value)] <- "unidentified"

sd(lipid_lipid)*3 #~0.49
table(abs(lipid_lipid)>0.49) ##100389

## filter correlation matrix by at least 2 features with 0.49 Tau per row/column 
filter_lipid_lipid <- colSums(abs(lipid_lipid) >  0.49) > 1

table(filter_lipid_lipid) ## true = 1460

pheatmap(lipid_lipid[filter_lipid_lipid,filter_lipid_lipid], show_rownames = F, show_colnames = F)


## annotations for heatmap   
annotation_col <- data.frame(lipid_class = as.numeric(as.factor(matching_names$metadata_value)), sig_QOL = as.numeric(matching_names$q_value < 0.05))
row.names(annotation_col) <- rownames(lipid_lipid)
annotation_col <- annotation_col[filter_lipid_lipid, ]

pheatmap(lipid_lipid[filter_lipid_lipid,filter_lipid_lipid], show_rownames = F, show_colnames = F, annotation_row = annotation_col)
# looks like a lot of different lipids are in the middle cluster (perhaps HDL or LDL particles?), interestingly the features significant with QOL are not in the middle cluster

## looking at some of the sig with QOL and nearest neighbors 
i = matching_names[matching_names$q_value <0.05,]$biomolecule_id[1]

for(i in matching_names[matching_names$q_value <0.05,]$biomolecule_id) {
  near_neighbors <- lipid_lipid[,grep(i,rownames(lipid_lipid))] 
  rownames(lipid_lipid)[which(near_neighbors >0.4)]
  pheatmap(lipid_lipid[which(near_neighbors >0.4),which(near_neighbors >0.4)], labels_row = matching_names$standardized_name[which(near_neighbors >0.4)])

}

#### Protein - lipid correlations #### 

str(protein_lipid)

# first column should be the row name
rownames(protein_lipid)<-protein_lipid[,1]
protein_lipid<- as.matrix(protein_lipid[,-1])

matching_names_lipids <- pvalues[match(rownames(protein_lipid), pvalues$biomolecule_id),]
matching_names_lipids <- cbind(matching_names_lipids, lipid_class[match(rownames(protein_lipid),lipid_class$biomolecule_id), ])
matching_names_lipids$metadata_value[is.na(matching_names_lipids$metadata_value)] <- "unidentified"

matching_names_protein <- pvalues[match(sub("X", "", colnames(protein_lipid)), pvalues$biomolecule_id),]
protein_genename <- biomolecules_metadata[biomolecules_metadata$metadata_type == "gene_name",] 
matching_names_protein <- cbind(matching_names_protein, protein_genename[match(matching_names_protein$biomolecule_id, protein_genename$biomolecule_id), ])

#dim(matching_names_protein)
#dim(protein_lipid)

sd(protein_lipid)*3 #~0.265
table(abs(protein_lipid)>0.27) ##41634

## filter correlation matrix by at least 2 features with 0.49 Tau per row/column 
filter_protein <- colSums(abs(protein_lipid) >  0.3) > 1
filter_lipid <- rowSums(abs(protein_lipid) > 0.3) > 1

table(filter_lipid) ## true = 621
table(filter_protein) ## 2722

pheatmap(protein_lipid[filter_lipid,filter_protein], show_rownames = F, show_colnames = F)

## Annotation rows and columns 

annotation_row <- data.frame(lipid_class = as.numeric(as.factor(matching_names_lipids$metadata_value)), sig_QOL = as.numeric(matching_names_lipids$q_value < 0.05))
row.names(annotation_row) <- rownames(protein_lipid)
#annotation_col <- annotation_col[filter_lipid_lipid, ]

annotation_col <- data.frame(sig_QOL = as.numeric(matching_names_protein$q_value < 0.05) )
row.names(annotation_col) <- paste("X", matching_names_protein$biomolecule_id, sep="")

## plot heatmap
pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/protein_lipid_cor_heatmap_0p4atleast2.pdf", width = 40)

pheatmap(protein_lipid[filter_lipid,filter_protein], labels_col = matching_names_protein$metadata_value[filter_protein], labels_row = matching_names_lipids$standardized_name[filter_lipid], annotation_row = annotation_row, annotation_col = annotation_col)

dev.off()


#### 
hist(-log(pvalues$q_value[pvalues$omics_id == 1]))
hist(-log(pvalues$q_value[pvalues$omics_id == 2]))
hist(-log(pvalues$q_value[pvalues$omics_id == 3]))

table(pvalues$omics_id[pvalues$q_value < 0.05])
table(pvalues$omics_id)
####

#### exploring protein-rnaseq correlations #### 


str(protein_rna)

# first column should be the row name
rownames(protein_rna)<-protein_rna[,1]
protein_rna<- as.matrix(protein_rna[,-1])

matching_names_rna <- pvalues[match(rownames(protein_rna), pvalues$biomolecule_id),]
rna_genesymbol <- biomolecules_metadata[biomolecules_metadata$metadata_type == "gene_symbol",] 
matching_names_rna <- cbind(matching_names_rna, rna_genesymbol[match(rownames(protein_rna),rna_genesymbol$biomolecule_id), ])

# matching_names_protein <- pvalues[match(sub("X", "", colnames(protein_lipid)), pvalues$biomolecule_id),]
# matching_names_protein <- cbind(matching_names_protein, protein_genename[match(matching_names_protein$biomolecule_id, protein_genename$biomolecule_id), ])

sd(protein_rna)*3 #~0.15
table(abs(protein_rna)>0.3) ##23097

## filter correlation matrix by at least 2 features with 0.49 Tau per row/column 
filter_protein <- colSums(abs(protein_rna) >  0.30) > 0
filter_rna <- rowSums(abs(protein_rna) > 0.30) > 0

table(filter_rna) ## true = 21
table(filter_protein) ## 12

pheatmap(protein_rna[filter_rna,filter_protein], show_rownames = F, show_colnames = F)

## Annotation rows and columns 

annotation_row <- data.frame(sig_QOL = as.numeric(matching_names_rna$q_value < 0.05))
row.names(annotation_row) <- rownames(protein_rna)
#annotation_col <- annotation_col[filter_lipid_lipid, ]

annotation_col <- data.frame(sig_QOL = as.numeric(matching_names_protein$q_value < 0.05) )
row.names(annotation_col) <- paste("X", matching_names_protein$biomolecule_id, sep="")

## plot heatmap
pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/protein_rna_cor_heatmap_0p27atleast1.pdf", width = 15, height = 22)

pheatmap(protein_rna[filter_rna,filter_protein], labels_col = matching_names_protein$metadata_value[filter_protein], labels_row = matching_names_rna$metadata_value[filter_rna], annotation_row = annotation_row, annotation_col = annotation_col)

dev.off()


