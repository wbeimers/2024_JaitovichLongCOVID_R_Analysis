library(DBI)
library(RSQLite)
library(pheatmap)

colors <- c("#72AF82", "#CF5E94", "#C07A9E", "#40ACE0", "#294D81")
colors2 <- c("#D66127", "#199D77", "#7670B2")

### for ASMS poster June 2025 

## from 07_KAO_Plotting_allOme_UMAP.R script

df<- read.csv( "D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome.csv")                         

#### correlation #### 

df_t <- t(df[,-c(1:2)])

df_corr<- cor(df[,-c(1:3)], method = "kendall")

write.csv(df_corr, "D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr.csv")

# read in 
df_corr <- read.csv("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr.csv", row.names = 1 )

#### plotting correlations #### 

filter <- rowSums(df_corr > abs(0.5)) > 5 
table(filter)

pheatmap(df_corr[filter,filter], kmeans_k = 50)

hclust_corr <- hclust(dist(df_corr[filter,filter]))

clusters <- cutree(as.hclust(hclust_corr), k = 50)
write.csv(clusters, "D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr_clusters.csv")

predictive_features <- c(match('normalized_counts.37456', names(clusters)),
                         match('normalized_counts.21322', names(clusters)),
                         match('normalized_counts.26264', names(clusters)),
                         match('normalized_counts.37825', names(clusters)),
                         match('normalized_abundance.300', names(clusters)),
                         match('normalized_abundance.2226', names(clusters)),
                         match('normalized_abundance.1814', names(clusters)),
                         match('normalized_abundance.40012', names(clusters)),
                         match('normalized_abundance.40088', names(clusters)))

filter2 <- rowSums(abs(df_corr[filter, filter][,clusters == 47]) > 0.3) > 5 
table(clusters)
pheatmap(df_corr[filter,filter][clusters == 47, filter2 ])

##### long format #### 

# read in 
df_corr <- read.csv("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr.csv", row.names = 1 )
row.names(df_corr) <- sub("normalized_abundance.","", row.names(df_corr))
row.names(df_corr) <- sub("normalized_counts.","", row.names(df_corr))

colnames(df_corr) <- sub("normalized_abundance.","", colnames(df_corr))
colnames(df_corr) <- sub("normalized_counts.","", colnames(df_corr))

clusters<- read.csv("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr_clusters.csv")


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)|lower.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

x<-Sys.time()

df_corr_long <- flattenCorrMatrix(df_corr)

y<-Sys.time()
x-y

write.csv(df_corr_long, "D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/all_ome_kendall_corr_long.csv")

rm(df_corr_long)
