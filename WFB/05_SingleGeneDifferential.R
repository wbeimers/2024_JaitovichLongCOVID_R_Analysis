####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggsignif)


#Colors#
#Make a classic palette
col <- brewer.pal(10, "Set2") 

#Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),"white", colorRampPalette(col1)(100))

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files#
all_df <- read.csv("data/processed/Plate_01-02_Samples_AcutevsPASC_Matrix.csv")
groups <- read.csv("data/metadata/sample_groups.csv")
file_info <- read.csv("data/metadata/plate_01-02_sample_file_info.csv")
file_info <- file_info[-1]



####Single Gene Distributions####
#select gene
g_one <- dplyr::filter(all_df, all_df$allgenes == "Q9H6X2")

#ggplot
g_one_plot <- data.frame(Samples = unlist(names(g_one[-1])),
                         Values = unlist(g_one[-1]),
                         Group = c(rep("Acute_fu", 17), rep("PASC", 187)))

ggplot(g_one_plot, aes(Group, log2(Values), fill = factor(Group))) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.06, alpha = 0.5) +
  geom_point() +
  geom_signif(comparisons = list(c("Acute_fu", "PASC")), 
              map_signif_level = TRUE,
              tip_length = 0.01,
              size = 0.6,
              textsize = 8,
              vjust = 0.5) +
  scale_fill_manual(values = c(col[2], col[5])) +
  ggtitle(paste("Intensity Distribution of", g_one[1,1])) +
  xlab(NULL) +
  ylab("log2(Intensity) (LFQ)") +
  #ylim(2, 8) +
  #scale_y_continuous(expand = expansion(0,0.1)) +
  guides(fill = guide_legend(title="Group")) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 32, face = 'bold', hjust = 0.5), 
        axis.title = element_text(face = 'bold'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) 
ggsave("reports/figures/AllPlates_Q9H6X2_PASCvsAcuteFU_ProteinIntensityDistribution.pdf", width = 24, height = 16, units = "cm")


table(grepl("P01579", spec_PG_NPs$allgenes))
