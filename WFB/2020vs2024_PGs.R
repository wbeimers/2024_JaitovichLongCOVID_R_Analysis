#### import old proteomics data from 2020 study to compare depth


library(devtools)
library(tidyverse)
library(DBI)
library(RSQLite)



#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Covid-19 Study DB.sqlite")


#### Pull Data ####

dbListTables(con)


df_proteins <- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Charlson_score, SOFA
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

table(!is.na(unique(df_proteins$biomolecule_id)))



#Do the sampe calculations for current study data (remove PGs containing more than 50% missing values)

#6100 protein groups

#make bar plot

PG_2020vs2024 <- data.frame(PGs = c(517, 6100),
                            Study = c("2020", "2024"))

ggplot(PG_2020vs2024, aes(Study, PGs)) + 
  geom_col(fill = col[7]) + 
  #scale_fill_manual(values = mycolors) +
  #scale_fill_manual(values = col) +
  #ggtitle("Protein Group Counts") +
  xlab("Study Year") +
  ylab("Protein Groups") +
  #ylim(0, 20) +
  scale_y_continuous(expand = c(0,0)) +
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
ggsave("reports/figures/50PGs_2020vs2024_Study.pdf", width = 24, height = 16, units = "cm")




