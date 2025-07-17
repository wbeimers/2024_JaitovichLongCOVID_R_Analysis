#### Distributions of participants ####

### Plots for ASMS poster June 2025

install.packages("beeswarm")


library(DBI)
library(RSQLite)
library(beeswarm)

colors <- c("#72AF82", "#CF5E94", "#C07A9E", "#40ACE0", "#294D81")

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

dbListTables(con)

patient_metadata <- dbReadTable(con, "patient_metadata")

dbDisconnect(con)

#### Plot age by group #### 

# remove Acute_NC 
patient_metadata <- patient_metadata[patient_metadata$Cohort != "Acute_NC",]
patient_metadata$Cohort <- relevel(as.factor(patient_metadata$Cohort), "Healthy") 


pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/age_by_group.pdf", height = 4, width = 9)
beeswarm(patient_metadata$Age ~ patient_metadata$Cohort, las = 1,
         ylab = "Age", xlab ="", bty = "l", col = colors, pch =19)
dev.off()

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/BMI_by_group.pdf", height = 4, width = 9)
beeswarm(patient_metadata$BMI ~ patient_metadata$Cohort, las = 1,
         ylab = "BMI", xlab ="", bty = "l", col = colors, pch =19)
dev.off()

table(patient_metadata$Sex, patient_metadata$Cohort)

pdf("D:/2024_LongCovid/2024_JaitovichLongCOVID_R_Analysis/KAO/plots/Sex_by_group.pdf", height = 9, width = 8)
par(mfrow = c(3,2))
for (i in 1:5){ 
pie(table(patient_metadata$Sex, patient_metadata$Cohort)[,i], col = c(colors[i], "white"), main = levels(patient_metadata$Cohort)[i])
}
dev.off()