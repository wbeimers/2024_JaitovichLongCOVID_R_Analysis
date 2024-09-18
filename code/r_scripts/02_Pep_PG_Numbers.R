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
spec_PG_NPs <- read.csv("data/processed/PG_Matrix_Plate_01-05_Samples_NPs.csv")
spec_pep_NPs <- read.csv("data/processed/pep_Matrix_Plate_01-05_Samples_NPs.csv")
groups <- read.csv("data/metadata/sample_groups.csv")
file_info <- read.csv("data/metadata/AllPlates_Samples_file_info.csv")
spec_PG_NPs <- spec_PG_NPs[-1]
file_info <- file_info[-1]


####plot:proteinvsmissing####
#Calc percent of samples each protein was detected in
percent_non_na <- rowMeans(!is.na(spec_PG_NPs)) * 100
non_na_counts_ordered <- percent_non_na[order(-percent_non_na)]

sum(percent_non_na == 100) #1215
sum(percent_non_na >= 50) #6095

#find percent between 
p <- rowSums(!is.na(spec_PG_NPs[, -c(1, 2)]))
op <- p[order(-p)]
fiftyp <- op[1:7608]
divided <- 100*fiftyp/40
p.mv <- mean(divided)
print(p.mv)

missing <- data.frame(Order = seq_along(non_na_counts_ordered), Value = non_na_counts_ordered)
missing$color <- "Y"

ggplot(missing, aes(Order, Value)) + 
  #geom_point(size = 3) + 
  geom_area(fill = col[8], alpha = 0.5) + # Shade the area under the line
  geom_line(color = col[3], linewidth = 2) + # Add the line on top
  geom_vline(xintercept = 1215, linetype = "solid", color = "black", linewidth = 1) +
  geom_vline(xintercept = 6095, linetype = "solid", color = col[2], linewidth = 1) + 
  geom_hline(yintercept = 50, linetype = "solid", color = col[2], linewidth = 1) +   
  #scale_fill_manual(values = col[5]) +
  ggtitle("Protein Group Counts") +
  xlab("Protein Groups") +
  ylab("% of Samples Quantified") +
  #ylim(0, 20) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,1)) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 101)) +
  guides(fill = guide_legend(title="Impairment\nStatus")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave("reports/figures/AllPlates_Samples_PGvsMissing.pdf", width = 24, height = 16, units = "cm")



####plot:RunOrdervsProteinNumber####
lengths <- colSums(!is.na(spec_PG_NPs))
barplot(lengths,las=2,border = NA)

lengths <- data.frame(experiment = names(lengths),
                      number = lengths)
lengths$runs <- sapply(strsplit(lengths$experiment, "_"), function(x) paste(x[5], collapse = "_"))
lengths <- merge(lengths, groups, by.x = "runs", by.y = "Sample", all.x = TRUE)
lengths <- merge(lengths, file_info, by.x = "runs", by.y = "Samples", all.x = TRUE)
lengths$runs <- factor(lengths$runs)

ggplot(lengths, aes(order, number, fill = factor(Set))) + 
  geom_col() + 
  #scale_fill_manual(values = mycolors) +
  scale_fill_manual(values = col) +
  ggtitle("Protein Group Counts") +
  xlab("Run Order") +
  ylab("Protein Group IDs") +
  #ylim(0, 20) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,1)) +
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
ggsave("reports/figures/AllPlates_Sample_PGvsRunOrder_Set_Bar.pdf", width = 32, height = 16, units = "cm")

