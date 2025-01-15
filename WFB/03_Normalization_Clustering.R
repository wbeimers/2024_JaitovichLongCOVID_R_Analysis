####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(missForest)
library(VIM)
library(viridis)


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
col1 <- colorRampPalette(col)(16)

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files#
spec_PG_NPs <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs.csv")
spec_PG_NPs_50 <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs_50%.csv")
groups <- read.csv("data/metadata/sample_groups.csv")
file_info <- read.csv("data/metadata/AllPlates_Samples_file_info.csv")
spec_PG_NPs <- spec_PG_NPs[-1]
spec_PG_NPs_50 <- spec_PG_NPs_50[-1]


####DataTransformation####
#Look at distribution of protein quant values in the unfiltered and filtered for 50% missingness of PGs
table(!is.na(spec_PG_NPs$X20240612_WFB_Plate_01_100))
table(!is.na(spec_PG_NPs_50$X20240612_WFB_Plate_01_100))

hist(spec_PG_NPs$X20240612_WFB_Plate_01_100,
     breaks = 100)
hist(log2(spec_PG_NPs$X20240612_WFB_Plate_01_100),
     breaks = 100)
hist(spec_PG_NPs_50$X20240612_WFB_Plate_01_100,
     breaks = 100)
hist(log2(spec_PG_NPs_50$X20240612_WFB_Plate_01_100),
     breaks = 100)

#Log2 transform the protein groups with at least 50% presence in the data
spec_PG_NPs_50_tf <- log2(spec_PG_NPs_50[-1])
spec_PG_NPs_50_tf <- data.frame(ID = spec_PG_NPs_50$PG.ProteinGroups, spec_PG_NPs_50_tf, stringsAsFactors = FALSE)

#look at distributions of missing values
matrixplot(spec_PG_NPs_50)
histMiss(c(spec_PG_NPs_50_tf$X20240612_WFB_Plate_01_100, spec_PG_NPs_50_tf$X20240626_WFB_Plate_04_10),
         only.miss = F)

####Imputation####
#write csv files to go impute in a different program
write.csv(spec_PG_NPs_50_tf, file = "data/processed/PG_Matrix_AllPlates_Samples_NPs_50%_tf.csv")

#read in imputed file from perseus
spec_PG_NPs_50_tf_imp <- read.csv("data/processed/PG_Matrix_AllPlates_Samples_NPs_50%_tf_imp.csv")

#check imputation
hist(spec_PG_NPs_50_tf$X20240612_WFB_Plate_01_100,
     breaks = 100)
hist(spec_PG_NPs_50_tf_imp$N..X20240612_WFB_Plate_01_100,
     breaks = 100)


####PCA####
pca_set <- spec_PG_NPs_50_tf_imp

pca_colnames <- colnames(pca_set)
pca_colnames <- sapply(strsplit(pca_colnames, "_"), function(x) paste(x[5], collapse = "_"))
pca_colnames[pca_colnames == "NA"] <- "allgenes"
colnames(pca_set) <- pca_colnames

t_pca_set <- as.data.frame(t(pca_set))

#get them into shape
colnames(t_pca_set) <- as.character(t_pca_set[401, ])
# Remove the row
t_pca_set <- t_pca_set[-401, ]
t_pca_set_n <- rownames(t_pca_set)
t_pca_set_1 <- as.data.frame(lapply(t_pca_set, function(x) as.numeric(as.character(x))))
rownames(t_pca_set_1) <- t_pca_set_n
t_pca_set <- t_pca_set_1
t_pca_set <- tibble::rownames_to_column(t_pca_set, "Samples")

#merge metadata
t_pca_set <- merge(t_pca_set, groups, by.x = "Samples", by.y = "Sample", all.x = T)
t_pca_set <- merge(t_pca_set, file_info, by = "Samples", all.x = T)
t_pca_set$plate <- sapply(strsplit(t_pca_set$Runs, "_"), function(x) paste(x[3:4], collapse = "_"))

pca_score <- prcomp(t_pca_set[,c(2:6088)],
                    scale. = T)
summary(pca_score)

#autoplot
pca_plot <- autoplot(pca_score,
                     data = t_pca_set,
                     color = "plate",
                     fill = "plate",
                     loadings = F,
                     loadings.label = F,
                     scale = 0,
                     frame = F,
                     frame.type = "t",
                     frame.color = "plate") +
  geom_point(aes(fill = plate), shape = 21, size = 4, color = "black") +
  #stat_ellipse(aes(fill = Set), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  #scale_fill_manual(values = col1) +
  #scale_color_manual(values = col1) +
  scale_color_viridis(option = "plasma", discrete = T) +
  scale_fill_viridis(option = "plasma", discrete = T) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) 
pca_plot
dev.off()
ggsave("reports/figures/AllPlates_Samples_50percentImputed_plate_PCA.pdf", width = 32, height = 16, units = "cm")

####scree plot####
plot(pca_score$sdev^2)



