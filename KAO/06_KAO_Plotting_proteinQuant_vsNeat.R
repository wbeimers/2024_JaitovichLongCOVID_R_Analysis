### Plot comparing Neat prep to Seer prep #### 
library(DBI)
library(RSQLite)
library(VennDiagram)

colors <- c("#72AF82", "#CF5E94", "#C07A9E", "#40ACE0", "#294D81")
colors2 <- c("#D66127", "#199D77", "#7670B2", "#FBB469")

### for ASMS poster June 2025 

neat <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Proteomics/20250224_NeatRePrep/Spectronaut19/20250301_131543_20250228_WFB_JaitovichLongCOVID_NeatRePrep_dDIA_stringent/20250228_WFB_JaitovichLongCOVID_NeatRePrep_dDIA_stringent_Report_WFB_Report (Normal).tsv"

#### read in files #### 
neat_df <- read.delim(neat)
names(neat_df)

neat_df_2 <- neat_df[,names(neat_df)=='PG.Quantity'|names(neat_df) == 'R.FileName' |names(neat_df) == 'PG.GroupLabel']
neat_df_2<-neat_df_2[!duplicated(neat_df_2),]
rm(neat_df)


quantified_PG_neat <- table(neat_df_2$PG.Quantity > 0, neat_df_2$R.FileName)
quantified_PG_neat[1,]

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbListTables(con)

patient_metadata<- dbReadTable(con, "patient_metadata")

protein <- dbGetQuery(con, "SELECT proteomics_measurement.biomolecule_id, proteomics_measurement.standardized_name, proteomics_measurement.rawfile_id, raw_abundance, normalized_abundance, patient_metadata.sample_id, patient_metadata.Cohort, patient_metadata.PG_IDs, patient_metadata.PG_change_collection_cutoff 
           FROM proteomics_measurement
           INNER JOIN rawfiles_all ON rawfiles_all.rawfile_id = proteomics_measurement.rawfile_id
           INNER JOIN patient_metadata ON patient_metadata.sample_id = rawfiles_all.sample_id
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = proteomics_measurement.biomolecule_id
           WHERE rawfiles_all.ome_id = 1
           AND rawfiles_all.keep = '1'
           AND biomolecules.keep = '1'
           ")


dbDisconnect(con)

#### mapping same samples neat and seer ####

neat_sample <- sub("20250225_WFB_JaitovichLongCOVID_NeatRerun_", "",colnames(quantified_PG_neat))
neat_sample <- sub("-",".",neat_sample, fixed=T)
quantified_PG_neat <- as.data.frame(t(quantified_PG_neat), neat_sample)
quantified_PG_neat$V3 <- as.numeric(quantified_PG_neat)

matched_samples <- merge(quantified_PG_neat, patient_metadata, by.x = 0, by.y = "Sample", all.x = T)

#### Plotting neat vs. seer number of protein groups #####

mean((as.numeric(matched_samples$PG_IDs) - as.numeric(matched_samples$V3) )/ as.numeric(matched_samples$V3))

table(matched_samples$PG_change_collection_cutoff ==0)


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/06_SeervsNeat_PG.pdf", height = 5, width = 9)
xdims <- barplot(as.numeric(matched_samples$PG_IDs)[order(matched_samples$PG_change_collection_cutoff)], ylim = c(0, 7500),border = NA,col = colors2[1],pch = 15, las = 1)
barplot(as.numeric(matched_samples$V3)[order(matched_samples$PG_change_collection_cutoff)], add = T, col = colors2[4], pch =19, las = 1,  border = NA )

abline(v= xdims[35,]+0.5, ) ## marking change in collection tube lot 
legend("topright", c("Proteograph XT", "Neat"), col = colors2[c(1,4)], pch = 15,bty = "n")
dev.off()

#### pca of protein data marking before after tube change ####
names(protein)

protein_wide_sample <- reshape(protein[,-c(2:4)], timevar = "biomolecule_id", v.names = "normalized_abundance",
                       idvar = "sample_id", direction = "wide" )


pca <- prcomp(protein_wide_sample[,-c(1:4)], scale = T)

plot(pca)
pca$sdev^2 / sum(pca$sdev^2)

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/06_PCA_proteins_prePostTubeChange.pdf", height = 5, width = 5)
plot(pca$x[,1], pca$x[,2], col = c(colors2[1],"gray20")[protein_wide_sample$PG_change_collection_cutoff+1], bg = colors2[1], pch = 21, bty = "l", 
     xlab = "principal component 1 (30.5%)", ylab = "principal component 2 (8.3%)", main = "Principal component analysis of protein data")
legend("topright", c("all samples", "new lot collection tube"), pt.bg = rep(colors2[1],2), col = c(colors2[1], "gray20"), pch = 21)
dev.off()

#### Venn Diagram of shared protein groups between seer and neat #### 
names(protein)
protein_wide <- reshape(protein[,-c(1,3,5, 7:9)], timevar = "sample_id", v.names = "raw_abundance",
                        idvar = "standardized_name", direction = "wide" )

row.names(protein_wide) <- protein_wide$standardized_name

neat_wide <- reshape(neat_df_2, timevar = "R.FileName", v.names = "PG.Quantity", idvar = "PG.GroupLabel", direction = "wide")

row.names(neat_wide) <- neat_wide$PG.GroupLabel

#count overlap
table(protein_wide$standardized_name %in% neat_wide$PG.GroupLabel)
length(protein_wide$standardized_name)

table(neat_wide$PG.GroupLabel %in% protein_wide$standardized_name)
length(neat_wide$PG.GroupLabel)

# plot Venn Diagram

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/06_proteins_Venn_seervNeat.pdf", height = 4, width = 4)

venn.plot <- draw.pairwise.venn(length(protein_wide$standardized_name), length(neat_wide$PG.GroupLabel), table(neat_wide$PG.GroupLabel %in% protein_wide$standardized_name)[2], fill = colors2[c(1,4)],
                                  category= c("Proteograph XT", "Neat"), main = "Protein overlap")
dev.off()

#### Correlation Neat vs. Seer #### 

## correlation of mean abundance 

seer_mean <- log2(rowMeans(protein_wide[protein_wide$standardized_name %in% neat_wide$PG.GroupLabel, -1], na.rm = T))
neat_mean <- log2(rowMeans(neat_wide[neat_wide$PG.GroupLabel %in% protein_wide$standardized_name, -1], na.rm = T))

matched_proteins <- merge(seer_mean, neat_mean, by = 0)

cor(matched_proteins[,2], matched_proteins[,3], method = "pearson")


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/06_proteins_cor_seervNeat.pdf", height = 4, width = 4)
plot(matched_proteins[,2], matched_proteins[,3], ylab = "Neat average log2(quant)", xlab = "Proteograph XT average log2(quant)", pch =19, col = colors2[1], bty = "l", las = 1)
legend("topleft", "Pearson's r = 0.55")
dev.off()
