#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(missForest)
library(VIM)
library(viridis)
library(RSQLite)
library(ggrepel)
library(pheatmap)
library(ROTS)
library(data.table)
library(enrichR)
library(enrichplot)
setEnrichrSite("Enrichr")
library(fgsea)
library(org.Hs.eg.db)


# Colors #
# Make a classic palette
col <- brewer.pal(8, "Set2") 

#pal <- c("#66C2A5",
#         "#FFD92F",
#         "#8DA0CB",
#         "#FC8D62",
#         "#A6CEE3",
#         "#E78AC3",
#         "#A6D854",
#         "#FDB462",
#         "#B3B3B3",
#         "#B2DF8A")

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')

# Make a Custom Gradient
col1 <- colorRampPalette(col)(16)

# plot colors
pie(rep(1, length(col)), col = col , main="") 


# Files #
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

metadata <- dbGetQuery(con, "SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples
                           FROM patient_metadata")
dbDisconnect(con)



#### plot:cohort pie chart ####
counts <- table(metadata$Cohort)
pdf("reports/figures/CohortSampleNumbers.pdf",
    width = 5,
    height = 4)
pie(counts,
    labels = counts,
    col = pal,
    main = "Cohorts")
legend("topright",
       legend = unique(counts),
       fill = pal,
       bty = "n")
dev.off()


#### plot:variable changes from PASC to PASC_fu ####
PASC_paired_metadata <- metadata %>%
  filter(!is.na(Paired_samples)) %>%
  filter(Cohort %in% c("PASC", "PASC_fu")) %>% 
  separate(Paired_samples, into = c("Sample_names", "timepoint"), sep = "_") %>%
  mutate(`SF.36.QOL.Score` = as.numeric(`SF.36.QOL.Score`))

ggplot(PASC_paired_metadata, 
       aes(timepoint, `SF.36.QOL.Score`, group = Sample_names)) +
  geom_point(size = 0.2,
             alpha = 0.5) +
  geom_line(linewidth = 0.2,
            alpha = 0.2) +
  ylim(c(0, 900)) +
  labs(x = "Measurement", 
       y = "QOL Score") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom") 
ggsave(paste0('reports/figures/Metadata_QOLchange_PASC_PASCfu.pdf'), 
       width = 4, height = 3, units = "cm")


#### check age+sex+bmi combos ####
metadata1 <- metadata %>%
  mutate(patientid = paste0(Age, "_", Sex, "_", BMI))


