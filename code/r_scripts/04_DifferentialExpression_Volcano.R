####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(perm)
library(ggrepel)


#Colors#
#Make a classic palette
col <- brewer.pal(8, "Set2") 

pal <- c("#66C2A5",
         "#FFD92F",
         "#8DA0CB",
         "#FC8D62",
         "#A6CEE3",
         "#E78AC3",
         "#A6D854",
         "#FDB462",
         "#B3B3B3",
         "#B2DF8A")

#Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),"white", colorRampPalette(col1)(100))

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files#
spec_PG_NPs <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")
spec_PG_NPs_50 <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs_50%.csv")
groups <- read.csv("data/metadata/sample_groups.csv")
file_info <- read.csv("data/metadata/AllPlates_Samples_file_info.csv")
spec_PG_NPs <- spec_PG_NPs[-1]
spec_PG_NPs_50 <- spec_PG_NPs_50[-1]

#change sample colnames
colnames <- colnames(spec_PG_NPs)
colnames <- sapply(strsplit(colnames, "_"), function(x) paste(x[5], collapse = "_"))
colnames[colnames == "NA"] <- "allgenes"
colnames(spec_PG_NPs) <- colnames


####DiffExp####
#split samples into groups to compare
pos_runs_list <- groups[groups$Set == "Acute", ]
pos_runs_list <- pos_runs_list$Sample

neg_runs_list <- groups[groups$Set == "PASC", ]
neg_runs_list <- neg_runs_list$Sample

#Here's where to change between datasets
pos_df <- spec_PG_NPs[, colnames(spec_PG_NPs) %in% pos_runs_list]
pos_df <- cbind(spec_PG_NPs["allgenes"], pos_df)

neg_df <- spec_PG_NPs[, colnames(spec_PG_NPs) %in% neg_runs_list]
neg_df <- cbind(spec_PG_NPs["allgenes"], neg_df)

#filter out NAs
filter_na <- rowSums(is.na(pos_df[,-1])) >= 97
pos_df <- pos_df[!filter_na,]

filter_na <- rowSums(is.na(neg_df[,-1])) >= 184
neg_df <-neg_df[!filter_na,]

all_df <- merge(pos_df, neg_df, by = "allgenes")
write.csv(all_df, file = "data/processed/AllPlates_NoFilter_Samples_PASCvsHealthy_Matrix.csv")

all_fc = NULL
for (i in 1:nrow(all_df)) {
  
  #get gene name
  gene <- all_df[i,1]
  
  #perform statistical test (ttest)
  NP_ttest <- t.test(x = as.numeric(all_df[i,colnames(pos_df)[-1]]),
                     y = as.numeric(all_df[i,colnames(neg_df)[-1]]),
                     paired = F,
                     var.equal = F)
  
  #perform permutation test
  NP_perm <- permTS(x = as.numeric(all_df[i,colnames(pos_df)[-1]]),
                    y = as.numeric(all_df[i,colnames(neg_df)[-1]]),
                    alternative = "two.sided",
                    method = "pclt",
                    control = permControl(nmc = 10000, seed = 500, digits = 5, tsmethod = "central"))
  
  #gives p.value of ttest
  t_pv <- NP_ttest$p.value
  
  #gives p.value of perm test
  p_pv <- NP_perm$p.value
  
  #gives the fold change between the two groups
  #fc <- 2^(mean(as.numeric(all_df[i,na.omit(neg_runs_list)])) - mean(as.numeric(all_df[i,na.omit(pos_runs_list)])))
  fc <- NP_ttest$estimate[2]/NP_ttest$estimate[1]
  
  all_fc <- rbind(all_fc, c(gene, fc, t_pv, p_pv))
  
}


all_fc <- as.data.frame(all_fc)
all_fc$`mean of y` <- as.numeric(all_fc$`mean of y`)
all_fc$V3 <- as.numeric(all_fc$V3)
all_fc$V4 <- as.numeric(all_fc$V4)
all_fc$t_padj <- p.adjust(all_fc$V3, method = "BH")
all_fc$p_padj <- p.adjust(all_fc$V4, method = "BH")
all_fc$log2fc <- log2(all_fc$`mean of y`)
all_fc$log10_t_padj <- -log10(all_fc$t_padj)
all_fc$log10_p_padj <- -log10(all_fc$p_padj)
all_fc$diffexp <- "NO"
all_fc$diffexp[all_fc$log2fc > 1 & all_fc$t_padj < 0.05] <- "UP"
all_fc$diffexp[all_fc$log2fc < -1 & all_fc$t_padj < 0.05] <- "DOWN"
colnames(all_fc) <- c("gene", "fc", "t_pv", "p_pv", "t_padj", "p_padj", "log2(fc)", "-log10(t_padj)", "-log10(p_padj)", "diffexp")


ggplot(all_fc, aes(`log2(fc)`, `-log10(t_padj)`, color = factor(diffexp), size = factor(diffexp))) + 
  geom_point() +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  scale_color_manual(values = c(col[1], col[8], col[5])) +
  scale_size_manual(values = c(3,1,3)) +
  scale_x_continuous(limits = c(-max(abs(all_fc$`log2(fc)`)), max(abs(all_fc$`log2(fc)`))), breaks = seq(-7, 7, by = 1)) +
  #geom_text_repel(data = subset(all_fc, diffexp == "DOWN"), aes(label = gene), size = 4) +
  ggtitle("PASC vs Acute") +
  xlab("Log2 Fold Change (PASC/Acute)") +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
  guides(fill = guide_legend(title = "Differetially Expressed")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) 
ggsave("reports/figures/AllPlates_NoFilter_Samples_PASCvsAcute_Volcano.pdf", width = 24, height = 16, units = "cm")

write.csv(all_fc, file = "data/processed/AllPlates_NoFilter_Samples_PASCvsAcute_DE_list.csv")

