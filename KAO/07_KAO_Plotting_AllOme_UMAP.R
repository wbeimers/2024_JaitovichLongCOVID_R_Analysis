### Plot comparing Neat prep to Seer prep #### 
library(DBI)
library(RSQLite)
library(VennDiagram)
install.packages("umap")
library(umap)

colors <- c("#72AF82", "#CF5E94", "#C07A9E", "#40ACE0", "#294D81")
colors2 <- c("#D66127", "#199D77", "#7670B2")

### for ASMS poster June 2025 


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

##### Creating df of lipid raw abundance ###### 
lipids <- dbGetQuery(con, "SELECT lipidomics_measurements.biomolecule_id, lipidomics_measurements.standardized_name,normalized_abundance, patient_metadata.sample_id, patient_metadata.Cohort
           FROM lipidomics_measurements
           INNER JOIN rawfiles_all ON rawfiles_all.rawfile_id = lipidomics_measurements.rawfile_id
           INNER JOIN patient_metadata ON patient_metadata.sample_id = rawfiles_all.sample_id
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles_all.ome_id = 2
           AND patient_metadata.analysis_group_1 = 1
           AND rawfiles_all.keep = '1'
           AND biomolecules.keep = '1'
           ")

protein <- dbGetQuery(con, "SELECT proteomics_measurement.biomolecule_id, proteomics_measurement.standardized_name, normalized_abundance, patient_metadata.sample_id, patient_metadata.Cohort 
           FROM proteomics_measurement
           INNER JOIN rawfiles_all ON rawfiles_all.rawfile_id = proteomics_measurement.rawfile_id
           INNER JOIN patient_metadata ON patient_metadata.sample_id = rawfiles_all.sample_id
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = proteomics_measurement.biomolecule_id
           WHERE rawfiles_all.ome_id = 1
           AND patient_metadata.analysis_group_1 = 1
           AND rawfiles_all.keep = '1'
           AND biomolecules.keep = '1'
           ")


transcripts <- dbGetQuery(con, "SELECT rnaseq_measurements.biomolecule_id, rnaseq_measurements.standardized_name, normalized_counts, patient_metadata.sample_id, patient_metadata.Cohort 
           FROM rnaseq_measurements
           INNER JOIN patient_metadata ON patient_metadata.sample_id = rnaseq_measurements.sample_id
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = rnaseq_measurements.biomolecule_id
           AND patient_metadata.analysis_group_1 = 1
           AND biomolecules.keep = '1'
           ")


dbDisconnect(con)

#### wide formats for merging #### 
names(transcripts)

protein_wide_sample <- reshape(protein[,-c(2)], timevar = "biomolecule_id", v.names = "normalized_abundance",
                               idvar = "sample_id", direction = "wide" )

lipids_wide_sample <- reshape(lipids[,-c(2)], timevar = "biomolecule_id", v.names = "normalized_abundance",
                               idvar = "sample_id", direction = "wide" )

transcripts_wide_sample <- reshape(transcripts[,-c(2)], timevar = "biomolecule_id", v.names = "normalized_counts",
                               idvar = "sample_id", direction = "wide" )

protein_lipid <- merge(protein_wide_sample, lipids_wide_sample[,-2], by = "sample_id")
all_ome <- merge(protein_lipid, transcripts_wide_sample[,-2], by = "sample_id")

all_ome_t <- t(all_ome[,-c(1:2)])

str(all_ome_t)
length(protein_wide_sample)

group <- c(rep(1, length(protein_wide_sample[-c(1:2)])), rep(2, length(lipids_wide_sample[-c(1:2)])), rep(3, length(transcripts_wide_sample[-c(1:2)])))

predictive_features <- c(match('normalized_counts.37456', names(all_ome[-c(1:2)])),
                         match('normalized_counts.21322', names(all_ome[-c(1:2)])),
                         match('normalized_counts.26264', names(all_ome[-c(1:2)])),
                         match('normalized_counts.37825', names(all_ome[-c(1:2)])),
                         match('normalized_abundance.300', names(all_ome[-c(1:2)])),
                         match('normalized_abundance.2226', names(all_ome[-c(1:2)])),
                         match('normalized_abundance.1814', names(all_ome[-c(1:2)])),
                         match('normalized_abundance.40012', names(all_ome[-c(1:2)])),
                         match('normalized_abundance.40088', names(all_ome[-c(1:2)])))


write.csv(all_ome, "D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome.csv")                         

#### UMAP #### 

set.seed(123)
umap_allome <- umap(all_ome_t, metric = "pearson", min_dist = 0.8)

head(umap_allome$layout)

umap.defaults

plot(umap_allome$layout, col = colors2[group])

umap_allome_v2 <- umap(all_ome_t, metric = "pearson2", min_dist = 0.8)

#pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/07_Umap_pearson2_allOme_highlightPredictors.pdf", height = 6, width = 6)
plot(umap_allome_v2$layout, col = colors2[group], pch=19, las = 1, bty = "l", ylab = "UMAP2", xlab = "UMAP1", cex = 0.5)
points(umap_allome_v2$layout[predictive_features,], col = "black", bg = colors2[c(3,3,3,3,1,1,1,2,2)], pch = 21, cex = 2)
text(umap_allome_v2$layout[predictive_features,], c("RNR1", "HBE1", "USP32P1", "XIST", "RGPD3", "ADAMTS2" , "FCN3", "Lipid RT 12.9 m/z 263.09", "Lipid RT 13.3 m/z 623.34"))
dev.off()

#### plotting XIST ### 

plot(all_ome_t[predictive_features[9],]~ as.factor(all_ome$Cohort))

#### Correlations #### 

