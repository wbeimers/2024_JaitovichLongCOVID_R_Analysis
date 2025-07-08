
#### Exploring the lipidomics measurements #### 
# I noticed there was some strange correlation relationships in the lipid-lipid correlation matrix
# and there are very few lipids that reach high significance in the analysis group QOL association

library(DBI)
library(RSQLite)
library(pheatmap)

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")


##### Creating df of lipid raw abundance ###### 
lipids <- dbGetQuery(con, "SELECT lipidomics_measurements.biomolecule_id, lipidomics_measurements.standardized_name, rawfiles_all.batch,raw_abundance, normalized_abundance, patient_metadata.sample_id, patient_metadata.Cohort
           FROM lipidomics_measurements
           INNER JOIN rawfiles_all ON rawfiles_all.rawfile_id = lipidomics_measurements.rawfile_id
           INNER JOIN patient_metadata ON patient_metadata.sample_id = rawfiles_all.sample_id
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles_all.ome_id = 2
           AND patient_metadata.analysis_group_1 = 1
           AND rawfiles_all.keep = '1'
           AND biomolecules.keep = '1'
           ")

pvalues <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, lratio, p_value, q_value, biomolecules.standardized_name, biomolecules.omics_id
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           WHERE biomolecules.keep = '1'  
           AND biomolecules.omics_id = '2'
           AND pvalues.analysis_group = '1'
           AND pvalues.formula = '1'
           ")


dbDisconnect(con)

#### Long to wide #### 
names(lipids) 

#long to wide for normalized abundance only 
lipids_wide <- reshape(lipids[,-c(3:4,7)], timevar = "sample_id", v.names = "normalized_abundance",
               idvar = "biomolecule_id", direction = "wide" )

sample_info <- unique(lipids[,c(3,6,7)])

##### Correlation matrix #### 
# I want to look at correlation matrix with and without batch1 data

lipid_wide_matrix <- lipids_wide[,-c(1:2)]
str(lipid_wide_matrix)
lipid_cor_all <- cor(t(lipid_wide_matrix), method = "kendall")
sample_cor_all <- cor(lipid_wide_matrix, method = "kendall")

lipid_cor_sans_batch1 <- cor(t(lipid_wide_matrix[,sample_info$batch != 1]), method = "kendall")


pheatmap(lipid_cor_all, show_rownames = F, show_colnames = F, k = 50)

pheatmap(lipid_cor_sans_batch1, show_rownames = F, show_colnames = F, k = 50)

hist(lipid_cor_all)
hist(lipid_cor_sans_batch1)

delta<- lipid_cor_sans_batch1-lipid_cor_all
hist(delta)

filter_delta <- rowSums(abs(delta) >0.25) >0
table(filter_delta)

pheatmap(delta[filter_delta, filter_delta], labels_row = lipids_wide[,2])

###### which lipids are significantly different in batch 1 ######

pvalues_batch1 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 1]), unname(x[sample_info$batch != 1]), var.equal = T)$p.value)
pvalues_batch2 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 2]), unname(x[sample_info$batch != 2]), var.equal = T)$p.value)
pvalues_batch3 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 3]), unname(x[sample_info$batch != 3]), var.equal = T)$p.value)
pvalues_batch4 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 4]), unname(x[sample_info$batch != 4]), var.equal = T)$p.value)
pvalues_batch5 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 5]), unname(x[sample_info$batch != 5]), var.equal = T)$p.value)
pvalues_batch6 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 6]), unname(x[sample_info$batch != 6]), var.equal = T)$p.value)
pvalues_batch7 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 7]), unname(x[sample_info$batch != 7]), var.equal = T)$p.value)
pvalues_batch8 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 8]), unname(x[sample_info$batch != 8]), var.equal = T)$p.value)
pvalues_batch9 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 9]), unname(x[sample_info$batch != 9]), var.equal = T)$p.value)
pvalues_batch10 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 10]), unname(x[sample_info$batch != 10]), var.equal = T)$p.value)

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_pvalues_hist.pdf", height = 10)
par(mfrow = c(4,3))
hist(pvalues_batch1)
hist(pvalues_batch2)
hist(pvalues_batch3)
hist(pvalues_batch4)
hist(pvalues_batch5)
hist(pvalues_batch6)
hist(pvalues_batch7)
hist(pvalues_batch8)
hist(pvalues_batch9)
hist(pvalues_batch10)
dev.off()

pvalues_batches <- cbind(pvalues_batch1, pvalues_batch2, pvalues_batch3, pvalues_batch4, pvalues_batch5,
                         pvalues_batch6, pvalues_batch7, pvalues_batch8, pvalues_batch9, pvalues_batch10)

qvalues_batches <- apply(pvalues_batches, 2, p.adjust, method = "BH")

table(rowSums(pvalues_batches < 0.05))

colSums(pvalues_batches < 0.05)

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_qvalueLessThan0p05.pdf")
barplot(colSums(qvalues_batches < 0.05), las = 1, names = seq(1,10), main = "# of lipids with batch-specific q<0.05")
dev.off()


lipids_wide[rowSums(pvalues_batches < 0.05)>6, 2]

par(mar = c(5,20,4,2))
barplot(-log(pvalues_batches[,1][order(pvalues_batches[,1])][1:12]), horiz = T, names = lipids_wide[rev(order(pvalues_batches[,1])),2][1:12], las = 1)
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lowqvalue_batch1.pdf", width = 20, height = 20)

par(mfrow = c(4,4))
for(i in 1:16){
  boxplot(unname(unlist((lipid_wide_matrix[order(qvalues_batches[,1]), ][i,]))) ~ sample_info$batch, main = paste(lipids_wide[order(qvalues_batches[,1]),2][i]),
          las = 1, ylab = "log2(Abundance)", xlab = "Batch")
  
}

dev.off()


## pvalues when batch1 is excluded 
pvalues2_batch2 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 2]), unname(x[sample_info$batch != 2|1]), var.equal = T)$p.value)
pvalues2_batch3 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 3]), unname(x[sample_info$batch != 3|1]), var.equal = T)$p.value)
pvalues2_batch4 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 4]), unname(x[sample_info$batch != 4|1]), var.equal = T)$p.value)
pvalues2_batch5 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 5]), unname(x[sample_info$batch != 5|1]), var.equal = T)$p.value)
pvalues2_batch6 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 6]), unname(x[sample_info$batch != 6|1]), var.equal = T)$p.value)
pvalues2_batch7 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 7]), unname(x[sample_info$batch != 7|1]), var.equal = T)$p.value)
pvalues2_batch8 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 8]), unname(x[sample_info$batch != 8|1]), var.equal = T)$p.value)
pvalues2_batch9 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 9]), unname(x[sample_info$batch != 9|1]), var.equal = T)$p.value)
pvalues2_batch10 <- apply(lipid_wide_matrix, 1, function(x) t.test(unname(x[sample_info$batch == 10]), unname(x[sample_info$batch != 10|1]), var.equal = T)$p.value)


pvalues2_batches <- cbind(pvalues2_batch2, pvalues2_batch3, pvalues2_batch4, pvalues2_batch5,
                         pvalues2_batch6, pvalues2_batch7, pvalues2_batch8, pvalues2_batch9, pvalues2_batch10)

qvalues2_batches <- apply(pvalues2_batches, 2, p.adjust, method = "BH")


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_pvalues_hist_withoutBatch1.pdf", height = 10)
par(mfrow = c(3,3))
hist(pvalues2_batch2)
hist(pvalues2_batch3)
hist(pvalues2_batch4)
hist(pvalues2_batch5)
hist(pvalues2_batch6)
hist(pvalues2_batch7)
hist(pvalues2_batch8)
hist(pvalues2_batch9)
hist(pvalues2_batch10)
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_qvalueLessThan0p05_sansBatch1.pdf", width = 5)
par(mfrow = c(2,1))
barplot(colSums(qvalues_batches < 0.05), las = 1, names = seq(1,10), ylim = c(0, 700), main = "# of lipids with batch-specific q<0.05")
barplot(c(0,colSums(qvalues2_batches < 0.05)), las = 1, names = seq(1,10), ylim = c(0, 700), main = "q<0.05 when Batch 1 is removed")
dev.off()

# what lipids are varying between batches 
table(rowSums(qvalues2_batches < 0.05))
lipids_wide[rowSums(qvalues2_batches < 0.05)>2, 2]

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_varyingByBatch.pdf", width = 10, height = 10)
par(mfrow = c(4,4))
for(i in 1:12){
  boxplot(unname(unlist((lipid_wide_matrix[rowSums(qvalues2_batches < 0.05)>3, ][i,]))) ~ sample_info$batch, main = paste(lipids_wide[rowSums(qvalues2_batches < 0.05)>3,2][i]),
          las = 1, ylab = "log2(Abundance)", xlab = "Batch")

}
dev.off()



pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/lipid_batches_varyingByBatch.pdf", width = 10, height = 10)
par(mfrow = c(4,4))
for(i in 1:16){
  boxplot(unname(unlist((lipid_wide_matrix[rowSums(pvalues_batch1 < 0.05)>5, ][i,]))) ~ sample_info$batch, main = paste(lipids_wide[rowSums(pvalues_batch1 < 0.05)>5,2][i]),
          las = 1, ylab = "log2(Abundance)", xlab = "Batch")
  
}
dev.off()


#### summarized batch pca to find most variying lipids #### 

mean_by_batch <- t(aggregate(t(lipid_wide_matrix), by = list(sample_info$batch), mean))
mean_by_batch <- mean_by_batch[-1,]

pca_batch_sansBatch1 <- prcomp(t(mean_by_batch[,-1]), scale = T)

pca_batch <- prcomp(t(mean_by_batch), scale = T)

plot(pca_batch_sansBatch1$x)
text(pca_batch_sansBatch1$x[,1], pca_batch_sansBatch1$x[,2], c(2:10), pos = 2)

plot(pca_batch$x)
text(pca_batch$x[,1], pca_batch$x[,2], c(1:10), pos = 2)

lipids_wide[rev(order(pca_batch_sansBatch1$rotation[,1])),2][1:10]
lipids_wide[rev(order(pca_batch_sansBatch1$rotation[,2])),2][1:10]

barplot(pca_batch_sansBatch1$rotation[,2][order(pca_batch_sansBatch1$rotation[,2])])
