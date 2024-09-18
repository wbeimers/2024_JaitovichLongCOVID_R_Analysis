####First set proper working directory under Session tab####
library(devtools)
library(tidyverse)
library(eulerr)


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


#Files#
spec_PG_NPA <- read.csv(file = "data/processed/PG_Matrix_AllPlates_QCsSamples_NPA.csv")
spec_PG_NPB <- read.csv(file = "data/processed/PG_Matrix_AllPlates_QCsSamples_NPB.csv")


PG_NPA <- spec_PG_NPA$PG.ProteinGroups
PG_NPB <- spec_PG_NPB$PG.ProteinGroups

sets <- list(NPA = PG_NPA, NPB = PG_NPB)

#euler
fit <- euler(sets)

plot(fit, 
     fills = list(fill = c(pal[5],pal[8]), alpha = 0.5),
     labels = list(cex = 1.5),
     edges = list(col = "black"),
     quantities = T)

