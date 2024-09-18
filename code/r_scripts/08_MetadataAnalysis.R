####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)

#Colors#
#Make a classic palette
col <- brewer.pal(10, "Set2") 

#Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),"white", colorRampPalette(col1)(100))

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files
groups <- read.csv("data/metadata/sample_groups.csv")



####plot:groups####
counts <- table(groups$Set)
pie(counts,
    labels = counts,
    col = col,
    main = "Group")
legend("topright",
       legend = unique(counts),
       fill = col,
       bty = "n")
dev.off()




category_counts <- as.data.frame(table(groups$Set))
colnames(category_counts) <- c("Category", "Count")

ggplot(category_counts, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = col) +
  theme_void() +
  ggtitle("PASC vs Acute") +
  guides(fill = guide_legend(title = "Group")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) 
