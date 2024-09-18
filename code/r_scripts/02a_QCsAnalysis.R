####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)


#Colors#
#Make a classic palette
col <- brewer.pal(8, "Set2") 

#Make another palette
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
spec_PG_QCs <- read.csv("data/processed/PG_Matrix_AllPlates_QCs.csv")
file_info <- read.csv("data/metadata/AllPlates_QCs_file_info.csv")
spec_PG_QCs <- spec_PG_QCs[-1]
file_info <- file_info[-1]
file_info <- file_info[-15,]
file_info$order <- seq_along(file_info$mtime)

names(spec_PG_QCs) <- sub(pattern = "X", replacement = "", x = names(spec_PG_QCs))
names(spec_PG_QCs) <- sub(pattern = ".x", replacement = "", x = names(spec_PG_QCs))
names(spec_PG_QCs) <- sub(pattern = ".y", replacement = "", x = names(spec_PG_QCs))



####Combine NPA and NPB search QCs together so I don't have to deal with both####
#use this to take either NPA or NPB for all samples for each pg, depending on which set is more complete
#make a filter to choose whichever row has fewer missing values for NPA or NPB. If same, choose NPA.
x <- apply(spec_PG_QCs[,grep(".x", names(spec_PG_QCs))], 1, function(x) table(unlist(unname(x)) > 0)[1])
x[is.na(x)] <-0

y <- apply(spec_PG_QCs[,grep(".y", names(spec_PG_QCs))], 1, function(x) table(unlist(unname(x)) > 0)[1])
y[is.na(y)] <-0

filter_xy <- x>=y
table(is.na(filter_xy))

#separate NPA and NPB into different dataframes
x_df <- spec_PG_QCs[,c(1,grep(".x", names(spec_PG_QCs)))]
y_df <- spec_PG_QCs[,c(1,grep(".y", names(spec_PG_QCs)))]

#apply filter to choose correct rows for either dataframe
x_df_filter <- x_df[filter_xy,]
y_df_filter <- y_df[!filter_xy,]
names(x_df_filter) <- sub(pattern = ".x", replacement = "", x = names(x_df_filter))
names(y_df_filter) <- sub(pattern = ".y", replacement = "", x = names(y_df_filter))
x_df_filter_1 <- x_df_filter[,c(1, order(names(x_df_filter)[-1])+1)]
y_df_filter_1 <- y_df_filter[,c(1, order(names(y_df_filter)[-1])+1)]

spec_PG_QCs <- rbind(x_df_filter, y_df_filter)



####plot:RunOrdervsProteinNumber####
lengths <- colSums(!is.na(spec_PG_QCs))
lengths <- data.frame(experiment = names(lengths),
                      number = lengths)
lengths$plate <- sapply(strsplit(lengths$experiment, "_"), function(x) paste(x[3:4], collapse = "_"))
lengths <- merge(lengths, file_info, by.x = "experiment", by.y = "Runs", all.x = TRUE)
lengths <- lengths[-51,]

ggplot(lengths, aes(order, number)) + 
  geom_point(aes(color = plate), size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  #scale_fill_manual(values = mycolors) +
  scale_color_manual(values = pal) +
  ggtitle("Protein Group Counts") +
  xlab("Run Order") +
  ylab("Protein Group IDs") +
  ylim(4000, 7000) +
  #scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title="Set")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave("reports/figures/AllPlates_QCs_PGvsRunOrder.pdf", width = 32, height = 16, units = "cm")



####CVs Analysis####
#make function to calculate rsd
rsd <- function(featureRow, index = 1:length(featureRow)){sd(featureRow[index], na.rm = T)/mean(featureRow[index], na.rm = T) * 100}

#add column of RSD across all QCs
spec_PG_QCss <- cbind(spec_PG_QCs, rsd_all = apply(spec_PG_QCs[,grepl("QC",colnames(spec_PG_QCs))], 1, rsd))
summQCss <- summary(spec_PG_QCss$rsd_all)

#add columns of RSDs within plates
spec_PG_QCss$rsd_plate_01 <- apply(spec_PG_QCss[, grepl("Plate_01",colnames(spec_PG_QCs))], 1, rsd)

names(spec_PG_QCss)
hist(spec_PG_QCss$rsd_plate_01,
     breaks = 100)
abline(v = summQCss[3], col = "red")

####plot rsd vs quant####
QCs_means <- apply(spec_PG_QCs[-1], 1, mean, na.rm = TRUE)
QCs_medians <- apply(spec_PG_QCs[-1], 1, median, na.rm = TRUE)

QCs_quant_rsd <- cbind(spec_PG_QCs$PG.ProteinGroups, QCs_means, QCs_medians, spec_PG_QCss$rsd_all)
QCs_quant_rsd <- as.data.frame(QCs_quant_rsd)
QCs_quant_rsd$QCs_means <- as.numeric(QCs_quant_rsd$QCs_means)
QCs_quant_rsd$QCs_medians <- as.numeric(QCs_quant_rsd$QCs_medians)
QCs_quant_rsd$V4 <- as.numeric(QCs_quant_rsd$V4)
hist(QCs_quant_rsd$V4,
     breaks = 100)
plot(log2(QCs_quant_rsd$QCs_medians), QCs_quant_rsd$V4)

ggplot(QCs_quant_rsd, aes(log2(QCs_medians), V4)) + 
  geom_point(size = 3, alpha = 0.05) + 
  #geom_smooth( se = FALSE, color = "red") +
  #scale_fill_manual(values = mycolors) +
  #scale_color_manual(values = pal) +
  ggtitle("QCs %RSD vs Median PG Quant") +
  xlab("log2(Median PG Quantity)") +
  ylab("%RSD") +
  #ylim(4000, 7000) +
  #scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title="Set")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave("reports/figures/AllPlates_QCs_PGQuantvsRSD.pdf", width = 32, height = 16, units = "cm")


#boxplot of quant for each in run order
spec_PG_QCs_long <- spec_PG_QCs[-1] %>%
  pivot_longer(cols = everything(), names_to = "Group", values_to = "Value")

spec_PG_QCs_long$plate <- sapply(strsplit(spec_PG_QCs_long$Group, "_"), function(x) paste(x[3:4], collapse = "_"))

ggplot(spec_PG_QCs_long, aes(Group, log2(Value), fill = plate)) + 
  geom_boxplot() + 
  #scale_fill_manual(values = mycolors) +
  scale_fill_manual(values = pal) +
  ggtitle("Post-Normalization PG Quantity Distribution") +
  xlab("Run Order") +
  ylab("log2(PG Quantity)") +
  #ylim(4000, 7000) +
  #scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title="Plate")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave("reports/figures/AllPlates_QCs_PostNormalizationQuantBoxplotvsRunOrder.pdf", width = 32, height = 16, units = "cm")



#boxplot of RSD for each run in run order



####PCA####
#only take complete cases for ease of analysis
spec_PG_QCs_100 <- spec_PG_QCs[complete.cases(spec_PG_QCs),]

#Log2 transform the protein groups with at least 50% presence in the data
spec_PG_QCs_100_tf <- log2(spec_PG_QCs_100[-1])
spec_PG_QCs_100_tf <- data.frame(ID = spec_PG_QCs_100$PG.ProteinGroups, spec_PG_QCs_100_tf, stringsAsFactors = FALSE)

#PCA
pca_set <- spec_PG_QCs_100_tf

pca_colnames <- colnames(pca_set)
pca_colnames <- sapply(strsplit(pca_colnames, "_"), function(x) paste(x[3:6], collapse = "_"))
pca_colnames[pca_colnames == "NA_NA_NA_NA"] <- "allgenes"
colnames(pca_set) <- pca_colnames

t_pca_set <- as.data.frame(t(pca_set))

#get them into shape
colnames(t_pca_set) <- as.character(t_pca_set[1, ])
# Remove the row
t_pca_set <- t_pca_set[-1, ]
t_pca_set <- mutate_all(t_pca_set, as.numeric)
t_pca_set <- tibble::rownames_to_column(t_pca_set, "Samples")

#merge metadata
t_pca_set$plate <- sapply(strsplit(t_pca_set$Samples, "_"), function(x) paste(x[1:2], collapse = "_"))

pca_score <- prcomp(t_pca_set[,c(2:3740)],
                    scale. = T)
summary(pca_score)

#autoplot
pca_plot <- autoplot(pca_score,
                     size = 5,
                     data = t_pca_set,
                     color = "plate",
                     loadings = F,
                     loadings.label = F,
                     scale = 0,
                     frame = F,
                     frame.type = "t",
                     frame.color = "plate") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  #scale_color_gradient(low = "#FFFFD4", high = "#CC4C02") +
  #scale_fill_gradient(low = "#FFFFD4", high = "#CC4C02") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) 
pca_plot
dev.off()
ggsave("reports/figures/AllPlates_QCs_completecases_plate_PCA.pdf", width = 24, height = 16, units = "cm")



