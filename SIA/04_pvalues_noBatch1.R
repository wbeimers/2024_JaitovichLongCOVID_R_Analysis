########################### load libraries ################################

library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(purrr) 
library(MASS)
library(lmtest)
library(ggplot2)
library(gridExtra)
library(effectsize)
library(tidyr)
library(lme4)
library(partR2)


########################### Exporting recent version of database ################################

source("H:/Coon Lab/Scripts/R/LongCOVID/XX_Export_Database_Files.R")
export_database_data("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite")

########################### Important Functions ################################

connect_db <- function(db_path) {
  dbConnect(RSQLite::SQLite(), dbname = db_path)
}

query_data <- function(con, table) {
  dbGetQuery(con, paste0("SELECT * FROM ", table))
}



#interaction terms fixed model

perform_lrt <- function(data, formula_null, formula_test, predictor, return = 'both') {
  results <- lapply(unique(data$biomolecule_id), function(biomolecule_id) {
    
    # Subset for biomolecule
    subset_data <- data[data$biomolecule_id == biomolecule_id, ]
    
    # Fit models
    set.seed(1007)
    model_null <- tryCatch(lm(formula_null, data = subset_data), error = function(e) NULL)
    set.seed(1007)
    model_test <- tryCatch(lm(formula_test, data = subset_data), error = function(e) NULL)
    
    # Initialize results
    lrt_lratio <- NA
    lrt_pvalue <- NA
    effect_size <- NA
    eta_sq <- NA
    
    # Compute LRT and effect size if both models valid
    if (!is.null(model_null) && !is.null(model_test)) {
      lrt <- tryCatch(lmtest::lrtest(model_null, model_test), error = function(e) NULL)
      if (!is.null(lrt)) {
        lrt_lratio <- lrt$Chisq[2]
        lrt_pvalue <- lrt$`Pr(>Chisq)`[2]
      }
      
      # Get effect size (coefficient)
      coefs <- tryCatch(summary(model_test)$coefficients, error = function(e) NULL)
      if (!is.null(coefs)) {
        matching_coef <- grep(paste0("^", predictor), rownames(coefs), value = TRUE)
        if (length(matching_coef) > 0) {
          effect_size <- coefs[matching_coef[1], "Estimate"]
        }
      }
      
      # Compute eta squared from ANOVA table (to support interaction terms)
      eta_result <- tryCatch(effectsize::eta_squared(anova(model_test), partial = TRUE), error = function(e) NULL)
      if (!is.null(eta_result) && predictor %in% eta_result$Parameter) {
        eta_sq <- eta_result$Eta2[eta_result$Parameter == predictor]
      }
    }
    
    # Return results for this biomolecule
    data.frame(
      biomolecule_id = biomolecule_id,
      predictor = predictor,
      effect_size = effect_size,
      eta_squared = eta_sq,
      lratio = lrt_lratio,
      p_value = lrt_pvalue
    )
  })
  
  # Combine all results
  return(do.call(rbind, results))
}


#no_r2 not bootstrapping
perform_lrt_mixed <- function(data, formula_null, formula_test, predictor) {
  results <- lapply(unique(data$biomolecule_id), function(biomolecule_id) {
    
    # Subset data
    subset_data <- data[data$biomolecule_id == biomolecule_id, ]
    
    # Fit mixed models with random intercept
    set.seed(1007)
    model_null <- tryCatch(
      lmer(update(formula_null, . ~ . + (1 | unique_patient_id)), data = subset_data, REML = FALSE),
      error = function(e) NULL
    )
    
    set.seed(1007)
    model_test <- tryCatch(
      lmer(update(formula_test, . ~ . + (1 | unique_patient_id)), data = subset_data, REML = FALSE),
      error = function(e) NULL
    )
    
    # Initialize return values
    lrt_lratio <- NA
    lrt_pvalue <- NA
    effect_size <- NA
    
    if (!is.null(model_null) && !is.null(model_test)) {
      # LRT
      lrt <- tryCatch(anova(model_null, model_test), error = function(e) NULL)
      if (!is.null(lrt) && nrow(lrt) == 2) {
        lrt_lratio <- lrt$Chisq[2]
        lrt_pvalue <- lrt$`Pr(>Chisq)`[2]
      }
      
      # Effect size
      coefs <- tryCatch(summary(model_test)$coefficients, error = function(e) NULL)
      if (!is.null(coefs)) {
        matching_coef <- grep(paste0("^", predictor), rownames(coefs), value = TRUE)
        if (length(matching_coef) > 0) {
          effect_size <- coefs[matching_coef[1], "Estimate"]
        }
      }
    }
    
    # Return all in one row
    data.frame(
      biomolecule_id = biomolecule_id,
      predictor = predictor,
      effect_size = effect_size,
      lratio = lrt_lratio,
      p_value = lrt_pvalue
    )
  })
  
  do.call(rbind, results)
}


# Function for paired t-test 
run_paired_ttest <- function(data) {
  data_wide <- data %>%
    pivot_wider(names_from = Cohort, values_from = normalized_abundance)
  
  if (all(c("PASC", "PASC_fu") %in% colnames(data_wide))) {
    t_result <- tryCatch(
      t.test(data_wide$`PASC`, data_wide$`PASC_fu`, paired = TRUE),
      error = function(e) NULL
    )
    
    fold_change <- log2(mean(data_wide$`PASC`, na.rm = TRUE) / mean(data_wide$`PASC_fu`, na.rm = TRUE))
    
    if (!is.null(t_result)) {
      return(data.frame(
        p_value = t_result$p.value,
        fold_change = fold_change
      ))
    }
  }
  
  return(data.frame(p_value = NA, fold_change = NA))
}

############################ Pull data from DB ###################################

db_path <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite"

con <- connect_db(db_path)

dbListTables(con)
dbListFields(con, "biomolecules")

metadata <- query_data(con, "patient_metadata")
biomolecules <- query_data(con, "biomolecules")
biomolecules_metadata <- query_data(con, "biomolecules_metadata")
rawfiles <- query_data(con, "rawfiles_all")
formulas <- query_data(con, "formula_table")
pvalues <- query_data(con, "pvalues")

df_proteins <- query_data(con, "proteomics_measurement")

df_lipids <- query_data(con, "lipidomics_measurements")

df_transcripts <- query_data(con, "rnaseq_measurements")


dbDisconnect(con)

############################ Sorting and Binding all omes data to metadata ###################################

df_proteins <- df_proteins %>%
  inner_join(rawfiles %>% 
               dplyr::select(-c(timestamp, run_type, rawfile_name)) %>%
               rename(keep_rawfiles = keep), by = "rawfile_id") %>%
  inner_join(metadata %>%
               dplyr::select(-Sample), by = "sample_id") %>%
  inner_join(biomolecules %>%
               dplyr::select(-standardized_name) %>%
               rename(keep_biomolecules = keep), by = "biomolecule_id") %>%
  filter(keep_rawfiles == "1", keep_biomolecules == "1")


df_lipids <- df_lipids %>%
  inner_join(rawfiles %>%
               dplyr::select(-timestamp, -run_type, -rawfile_name) %>%
               rename(keep_rawfiles = keep), by = "rawfile_id") %>%
  inner_join(metadata %>%
               dplyr::select(-Sample), by = "sample_id") %>%
  inner_join(biomolecules %>%
               dplyr::select(-standardized_name) %>%
               rename(keep_biomolecules = keep), by = "biomolecule_id") %>%
  filter(keep_rawfiles == "1", keep_biomolecules == "1",
         batch != 1)


df_transcripts <- df_transcripts%>%
  rename(normalized_abundance = normalized_counts) %>%
  inner_join(metadata, by = "sample_id") %>%
  inner_join(biomolecules %>% 
               dplyr::select(-standardized_name) %>%
               rename(keep_biomolecules = keep), by = "biomolecule_id") %>%
  filter(keep_biomolecules == "1")                             


common_columns <- Reduce(intersect, list(colnames(df_proteins), 
                                         colnames(df_lipids), 
                                         colnames(df_transcripts)))

# Select and reorder columns that are matching only
df_proteins <- df_proteins %>% dplyr::select(all_of(common_columns))
df_lipids <- df_lipids %>% dplyr::select(all_of(common_columns))
df_transcripts <- df_transcripts %>% dplyr::select(all_of(common_columns))


#Bind all omes
df <- bind_rows(df_proteins, df_lipids, df_transcripts)


#Add more metadata/grouping columns 
df <- df %>%
  mutate(
    SF.36.QOL.Score = ifelse(is.na(SF.36.QOL.Score) | SF.36.QOL.Score == "NULL", 900, as.numeric(SF.36.QOL.Score)),
    PASC_noPASC = if_else(grepl("PASC", Cohort), "PASC", "noPASC"),
    Acute_Healthy = case_when(
      grepl("^Acute$", Cohort) ~ "Acute",
      grepl("^Healthy$", Cohort) ~ "Healthy",
      TRUE ~ NA_character_
    ),
    Acute_nonAcute = case_when(
      grepl("^Acute$", Cohort) ~ "Acute",
      grepl("^Healthy$", Cohort) | grepl("Acute_fu", Cohort) ~ "nonAcute",
      TRUE ~ NA_character_
    )
  )

#remove duplicated columns if any
colnames(df)
which(duplicated(names(df)))
names(df)[duplicated(names(df))]
df <- df[, !duplicated(names(df))]


#Subset Analysis Groups
df_subset1 <- df[df$analysis_group_1 == 1, ]
df_subset2 <- df[df$analysis_group_2 == 1, ]
df_subset3 <- df[df$analysis_group_3 == 1, ]
df_subset4 <- df[df$analysis_group_4 == 1, ]
df_subset5 <- df[df$analysis_group_5 == 1, ]
df_subset6 <- df[df$analysis_group_6 == 1, ]
df_subset7 <- df[df$analysis_group_7 == 1, ]


######################################  ANALYSIS GROUP 1 ################################

##########################  Calculating pvalues - LRT ################################

##### P-values for QoL ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_QoL1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score"
)

df_pvalues_QoL1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "QoL", formula = "22")

df_pvalues_QoL1 <- df_pvalues_QoL1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoL1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )



df_pvalues_QoL1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_QoL1), by =1 ), df_pvalues_QoL1)


hist(df_pvalues_QoL1$p_value, breaks =100, main = "Histogram of pvalues - QoL")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoL1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_age1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "Age", formula = "23")

df_pvalues_age1 <- df_pvalues_age1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_age1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age1), by =1 ), df_pvalues_age1)

hist(df_pvalues_age1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_sex1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)

df_pvalues_sex1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "Sex", formula = "24")

df_pvalues_sex1 <- df_pvalues_sex1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_sex1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex1), by =1 ), df_pvalues_sex1)

hist(df_pvalues_sex1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")

lrt_results_bmi1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "BMI", formula = "25")

df_pvalues_bmi1 <- df_pvalues_bmi1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi1), by =1 ), df_pvalues_bmi1)

hist(df_pvalues_bmi1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### Pvalues w/ QoL and Age interaction #####

rownames(summary(lm(formula_test, data = test))$coefficients)

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score+ Age + Sex + BMI + Age:Sex")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score+ Age + Sex + BMI + SF.36.QOL.Score:Age + Age:Sex")

lrt_results_QoLAge1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score:Age"
)


df_pvalues_QoLAge1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "QoL/Age interaction", formula = "26")

df_pvalues_QoLAge1 <- df_pvalues_QoLAge1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoLAge1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_QoLAge1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_QoLAge1), by =1 ), df_pvalues_QoLAge1)

hist(df_pvalues_QoLAge1$p_value, breaks = 100, main = 'Histogram of pvalues - QoL/Age interaction')


## Append interaction Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoLAge1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### Pvalues w/ Age and Sex interaction #####

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score+ Age + Sex + BMI + SF.36.QOL.Score:Age")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score+ Age + Sex + BMI + SF.36.QOL.Score:Age + Age:Sex")

lrt_results_AgeSex1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age:Sex"
)

df_pvalues_AgeSex1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "Age/Sex interaction", formula = "27")

df_pvalues_AgeSex1 <- df_pvalues_AgeSex1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_AgeSex1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_AgeSex1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_AgeSex1), by =1 ), df_pvalues_AgeSex1)

hist(df_pvalues_AgeSex1$p_value, breaks = 100, main = 'Histogram of pvalues - Age/Sex interaction')


## Append interaction Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_AgeSex1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



########################## Saving Histograms - Group 1 ##########################

# Create individual plots
plot_qol <- ggplot(df_pvalues_QoL1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "steelblue") + 
  ggtitle("QoL")

plot_age <- ggplot(df_pvalues_age1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex <- ggplot(df_pvalues_sex1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi <- ggplot(df_pvalues_bmi1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")

plot_qol_age <- ggplot(df_pvalues_QoLAge1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#FFDAC1") + 
  ggtitle("QoL/Age Interaction")

plot_age_sex <- ggplot(df_pvalues_AgeSex1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#FFB7C5") + 
  ggtitle("Age/Sex Interaction")

# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group1_v2noBatch1noBatch1_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_qol, plot_age, plot_sex, plot_bmi, plot_qol_age, plot_age_sex, ncol = 2)
dev.off()


######################################  ANALYSIS GROUP 2 ################################

##########################  Calculating pvalues - LRT ################################
##### P-values for QoL ###### 

formula_null <- as.formula("normalized_abundance ~  BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + BMI")


lrt_results_QoL2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score"
)

df_pvalues_QoL2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "QoL", formula = "28")

df_pvalues_QoL2 <- df_pvalues_QoL2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoL2, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

lrt_results_QoL2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(lrt_results_QoL2), by =1 ), lrt_results_QoL2)

hist(df_pvalues_QoL2$p_value, breaks =100, main = "Histogram of pvalues - QoL")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoL2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + BMI")

lrt_results_bmi2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "BMI", formula = "29")

df_pvalues_bmi2 <- df_pvalues_bmi2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi2, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi2), by =1 ), df_pvalues_bmi2)

hist(df_pvalues_bmi2$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 2 ##########################


# Create individual plots
plot_qol <- ggplot(df_pvalues_QoL2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "steelblue") + 
  ggtitle("QoL")

plot_bmi <- ggplot(df_pvalues_bmi2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")

# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group2noBatch1noBatch1_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_qol, plot_bmi, ncol = 2)
dev.off()


######################################  ANALYSIS GROUP 2 - different formulas ################################


##### P-values for QoL ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_QoL2_2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score"
)

df_pvalues_QoL2_2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "QoL", formula = "30")

df_pvalues_QoL2_2 <- df_pvalues_QoL2_2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoL2_2 , by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_QoL2_2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_QoL2_2), by =1 ), df_pvalues_QoL2_2)

hist(df_pvalues_QoL2_2$p_value, breaks =100, main = "Histogram of pvalues - QoL")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoL2_2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_age2_2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age2_2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "Age", formula = "31")

df_pvalues_age2_2 <- df_pvalues_age2_2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age2_2, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age2_2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age2_2), by =1 ), df_pvalues_age2_2)

hist(df_pvalues_age2_2$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age2_2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_sex2_2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex2_2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "Sex", formula = "32")

df_pvalues_sex2_2 <- df_pvalues_sex2_2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex2_2, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex2_2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex2_2), by =1 ), df_pvalues_sex2_2)

hist(df_pvalues_sex2_2$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex2_2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")

lrt_results_bmi2_2 <- perform_lrt(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi2_2 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "BMI", formula = "33")

df_pvalues_bmi2_2 <- df_pvalues_bmi2_2 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi2_2, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi2_2 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi2_2), by =1 ), df_pvalues_bmi2_2)

hist(df_pvalues_bmi2_2$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi2_2, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



########################## Saving Histograms - Group 2 - Age and Sex Added to the Model ##########################

plot_qol <- ggplot(df_pvalues_QoL2_2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "steelblue") + 
  ggtitle("QoL")

plot_age <- ggplot(df_pvalues_age2_2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex <- ggplot(df_pvalues_sex2_2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi <- ggplot(df_pvalues_bmi2_2, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")

# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group2_AgeSexAddedtotheModelnoBatch1noBatch1_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_qol, plot_age, plot_sex, plot_bmi, ncol = 2)
dev.off()




########################## Paired T-test - Group 2 ##########################

#Per patient result (Paired_T_Test)
df_summary <- df_subset2 %>%
  group_by(unique_patient_id, Cohort) %>%
  summarize(mean_abundance = mean(normalized_abundance, na.rm = TRUE), .groups = "drop")

View(df_wide)
# 2. Pivot wider to get paired values side-by-side
df_wide <- df_summary %>%
  pivot_wider(names_from = Cohort, values_from = mean_abundance)

# 3. Run the paired t-test
Paired_t_result <- t.test(df_wide$PASC, df_wide$PASC_fu, paired = TRUE)



#Per biomoelcule_id result (Paired_T_Test)

paired_results <- df_subset2 %>%
  group_by(biomolecule_id, unique_patient_id, Cohort) %>%
  summarize(mean_abundance = mean(normalized_abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cohort, values_from = mean_abundance) %>%
  filter(!is.na(PASC), !is.na(PASC_fu)) %>%
  group_by(biomolecule_id) %>%
  summarize(
    p_value = tryCatch(t.test(PASC, PASC_fu, paired = TRUE)$p.value, error = function(e) NA),
    effect_size = log2(mean(PASC, na.rm = TRUE) / mean(PASC_fu, na.rm = TRUE)),
    .groups = "drop"
  )

# 2. Wrap into consistent format
df_paired_pvalues <- data.frame(
  biomolecule_id = paired_results$biomolecule_id,
  analysis_group = "2",  # or any group number you're working with
  test = "paired_test",
  comparison = "PASC_vs_PASC_fu",
  formula = "34"
) %>%
  left_join(paired_results, by = "biomolecule_id") %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

# 3. Add unique pvalue_id that continues from previous entries in pvalues
df_paired_pvalues <- cbind(
  pvalue_id = seq(nrow(pvalues) + 1, length.out = nrow(df_paired_pvalues), by = 1),
  df_paired_pvalues
)

# 4. (Optional) View histogram
hist(df_paired_pvalues$p_value, breaks = 100, main = "Histogram of Paired t-test p-values")


## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_paired_pvalues, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

########################## PASCvsnoPASC - Group 1 ###############################

##### P-values for PASCnoPASC ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ PASC_noPASC + Age + Sex + BMI")


lrt_results_PASCnoPASC <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "PASC_noPASC"
)

df_pvalues_PASCnoPASC <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "PASC_noPASC", formula = "35")

df_pvalues_PASCnoPASC <- df_pvalues_PASCnoPASC %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_PASCnoPASC, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_PASCnoPASC <- cbind(pvalue_id = row.names(df_pvalues_PASCnoPASC), df_pvalues_PASCnoPASC)

hist(df_pvalues_PASCnoPASC$p_value, breaks =100, main = "Histogram of pvalues - PASC no PASC")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_PASCnoPASC, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ PASC_noPASC + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ PASC_noPASC + Age + Sex + BMI")


lrt_results_age1_1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age1_1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "Age", formula = "36")

df_pvalues_age1_1 <- df_pvalues_age1_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age1_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age1_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age1_1), by =1 ), df_pvalues_age1_1)

hist(df_pvalues_age1_1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age1_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ PASC_noPASC + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ PASC_noPASC + Age + Sex + BMI")


lrt_results_sex1_1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex1_1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "Sex", formula = "37")

df_pvalues_sex1_1 <- df_pvalues_sex1_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex1_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex1_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex1_1), by =1 ), df_pvalues_sex1_1)

hist(df_pvalues_sex1_1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex1_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ PASC_noPASC + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ PASC_noPASC + Age + Sex + BMI")

lrt_results_bmi1_1 <- perform_lrt(
  data = df_subset1,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi1_1 <- data.frame(biomolecule_id = unique(df_subset1$biomolecule_id), analysis_group = "1", test = "LR_test", comparison = "BMI", formula = "38")

df_pvalues_bmi1_1 <- df_pvalues_bmi1_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi1_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi1_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi1_1), by =1 ), df_pvalues_bmi1_1)

hist(df_pvalues_bmi1_1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi1_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 1 PascnoPASC ##########################

# Create individual plots
plot_PAScnoPASC <- ggplot(df_pvalues_PASCnoPASC, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "maroon") + 
  ggtitle("QoL")

plot_age1_1 <- ggplot(df_pvalues_age1_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex1_1 <- ggplot(df_pvalues_sex1_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi1_1 <- ggplot(df_pvalues_bmi1_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")


# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group1_PASCnoPASCnoBatch1_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_PAScnoPASC, plot_age1_1, plot_sex1_1, plot_bmi1_1, ncol = 2)
dev.off()



########################## Analysis Group 2 - Random intercept ###############################


######################################  ANALYSIS GROUP 2 - different formulas - mixed effect model ################################


##### P-values for QoL ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI  + (1 | unique_patient_id)")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI  + (1 | unique_patient_id)")

lrt_results_QoL2_3 <- perform_lrt_mixed(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score"
)

df_pvalues_QoL2_3 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "QoL", formula = "39")

df_pvalues_QoL2_3 <- df_pvalues_QoL2_3 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoL2_3 , by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_QoL2_3 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_QoL2_3), by =1 ), df_pvalues_QoL2_3)

hist(df_pvalues_QoL2_3$p_value, breaks =100, main = "Histogram of pvalues - QoL")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoL2_3, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Sex + BMI + (1 | unique_patient_id)")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI + (1 | unique_patient_id)")


lrt_results_age2_3 <- perform_lrt_mixed(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age2_3 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "Age", formula = "40")

df_pvalues_age2_3 <- df_pvalues_age2_3 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age2_3, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age2_3 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age2_3), by =1 ), df_pvalues_age2_3)

hist(df_pvalues_age2_3$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age2_3, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + BMI + (1 | unique_patient_id)")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI + (1 | unique_patient_id)")


lrt_results_sex2_3 <- perform_lrt_mixed(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex2_3 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "Sex", formula = "41")

df_pvalues_sex2_3 <- df_pvalues_sex2_3 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex2_3, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex2_3 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex2_3), by =1 ), df_pvalues_sex2_3)

hist(df_pvalues_sex2_3$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex2_3, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex  + (1 | unique_patient_id)")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI  + (1 | unique_patient_id)")

lrt_results_bmi2_3 <- perform_lrt_mixed(
  data = df_subset2,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi2_3 <- data.frame(biomolecule_id = unique(df_subset2$biomolecule_id), analysis_group = "2", test = "LR_test", comparison = "BMI", formula = "42")

df_pvalues_bmi2_3 <- df_pvalues_bmi2_3 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi2_3, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi2_3 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi2_3), by =1 ), df_pvalues_bmi2_3)

hist(df_pvalues_bmi2_3$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi2_3, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



########################## Saving Histograms - Group 2 - Age and Sex Added to the Model ##########################

plot_qol2_3 <- ggplot(df_pvalues_QoL2_3, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "steelblue") + 
  ggtitle("QoL")

plot_age2_3 <- ggplot(df_pvalues_age2_3, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex2_3 <- ggplot(df_pvalues_sex2_3, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi2_3 <- ggplot(df_pvalues_bmi2_3, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")

# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group2_AgeSexAddedtotheModel_MixedEffectnoBatch1_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_qol2_3, plot_age2_3, plot_sex2_3, plot_bmi2_3, ncol = 2)
dev.off()



########################################### Analysis Group 4 ##########################################
########################## Acute/Healthy - Group 4 ###############################

##### P-values for Acute/Healthy ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_Healthy + Age + Sex + BMI")


lrt_results_AcuteHealthy <- perform_lrt(
  data = df_subset4,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Acute_Healthy"
)

df_pvalues_AcuteHealthy <- data.frame(biomolecule_id = unique(df_subset4$biomolecule_id), analysis_group = "4", test = "LR_test", comparison = "Acute_Healthy", formula = "43")

df_pvalues_AcuteHealthy <- df_pvalues_AcuteHealthy %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_AcuteHealthy, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_AcuteHealthy <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_AcuteHealthy), by =1 ), df_pvalues_AcuteHealthy)

hist(df_pvalues_AcuteHealthy$p_value, breaks =100, main = "Histogram of pvalues - Acute/Healthy")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_AcuteHealthy, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ Acute_Healthy + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_Healthy + Age + Sex + BMI")

lrt_results_age4_1 <- perform_lrt(
  data = df_subset4,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age4_1 <- data.frame(biomolecule_id = unique(df_subset4$biomolecule_id), analysis_group = "4", test = "LR_test", comparison = "Age", formula = "44")

df_pvalues_age4_1 <- df_pvalues_age4_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age4_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age4_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age4_1), by =1 ), df_pvalues_age4_1)

hist(df_pvalues_age4_1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age4_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ Acute_Healthy + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_Healthy + Age + Sex + BMI")


lrt_results_sex4_1 <- perform_lrt(
  data = df_subset4,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex4_1 <- data.frame(biomolecule_id = unique(df_subset4$biomolecule_id), analysis_group = "4", test = "LR_test", comparison = "Sex", formula = "45")

df_pvalues_sex4_1 <- df_pvalues_sex4_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex4_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex4_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex4_1), by =1 ), df_pvalues_sex4_1)

hist(df_pvalues_sex4_1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex4_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ Acute_Healthy + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ Acute_Healthy + Age + Sex + BMI")

lrt_results_bmi4_1 <- perform_lrt(
  data = df_subset4,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)

df_pvalues_bmi4_1 <- data.frame(biomolecule_id = unique(df_subset4$biomolecule_id), analysis_group = "4", test = "LR_test", comparison = "BMI", formula = "46")

df_pvalues_bmi4_1 <- df_pvalues_bmi4_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi4_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi4_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi4_1), by =1 ), df_pvalues_bmi4_1)

hist(df_pvalues_bmi4_1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi4_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 4 Acute/Healthy ##########################

# Create individual plots
plot_AcuteHealthy <- ggplot(df_pvalues_AcuteHealthy, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "black") + 
  ggtitle("Acute/Healthy")

plot_age4_1 <- ggplot(lrt_results_age4_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex4_1 <- ggplot(df_pvalues_sex4_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi4_1 <- ggplot(df_pvalues_bmi4_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")


# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group4_AcuteHealthy_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_AcuteHealthy, plot_age4_1, plot_sex4_1, plot_bmi4_1, ncol = 2)
dev.off()

########################################### Analysis Group 5 ##########################################
########################## PASC/noPASC - Group 5 ###############################

##### P-values for Acute/Healthy ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_QoL5_1 <- perform_lrt(
  data = df_subset5,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "SF.36.QOL.Score"
)

df_pvalues_QoL5_1 <- data.frame(biomolecule_id = unique(df_subset5$biomolecule_id), analysis_group = "5", test = "LR_test", comparison = "QoL", formula = "47")

df_pvalues_QoL5_1 <- df_pvalues_QoL5_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_QoL5_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_QoL5_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_QoL5_1), by =1 ), df_pvalues_QoL5_1)

hist(df_pvalues_QoL5_1$p_value, breaks =100, main = "Histogram of pvalues - QoL")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_QoL5_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")

lrt_results_age5_1 <- perform_lrt(
  data = df_subset5,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age5_1 <- data.frame(biomolecule_id = unique(df_subset5$biomolecule_id), analysis_group = "5", test = "LR_test", comparison = "Age", formula = "48")

df_pvalues_age5_1 <- df_pvalues_age5_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age5_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age5_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age5_1), by =1 ), df_pvalues_age5_1)

hist(df_pvalues_age5_1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age5_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")


lrt_results_sex5_1 <- perform_lrt(
  data = df_subset5,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex5_1 <- data.frame(biomolecule_id = unique(df_subset5$biomolecule_id), analysis_group = "5", test = "LR_test", comparison = "Sex", formula = "49")

df_pvalues_sex5_1 <- df_pvalues_sex5_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex5_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex5_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex5_1), by =1 ), df_pvalues_sex5_1)

hist(df_pvalues_sex5_1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex5_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ SF.36.QOL.Score + Age + Sex + BMI")

lrt_results_bmi5_1 <- perform_lrt(
  data = df_subset5,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)

df_pvalues_bmi5_1 <- data.frame(biomolecule_id = unique(df_subset5$biomolecule_id), analysis_group = "5", test = "LR_test", comparison = "BMI", formula = "50")

df_pvalues_bmi5_1 <- df_pvalues_bmi5_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi5_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi5_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi5_1), by =1 ), df_pvalues_bmi5_1)

hist(df_pvalues_bmi5_1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi5_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 5 PASC/PASC fu ##########################

# Create individual plots
plot_QoL5_1 <- ggplot(df_pvalues_QoL5_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "steelblue") + 
  ggtitle("QoL")

plot_age5_1 <- ggplot(df_pvalues_age5_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex5_1 <- ggplot(df_pvalues_sex5_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi5_1 <- ggplot(df_pvalues_bmi5_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")


# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group5_QoL_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_QoL5_1, plot_age5_1, plot_sex5_1, plot_bmi5_1, ncol = 2)
dev.off()


########################################### Analysis Group 6 ##########################################
########################## Acute/nonAcute - Group 6 ###############################

##### P-values for Acute/nonAcute ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + Sex + BMI")

lrt_results_AcutenonAcute <- perform_lrt(
  data = df_subset6,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Acute_nonAcute"
)

df_pvalues_AcutenonAcute <- data.frame(biomolecule_id = unique(df_subset6$biomolecule_id), analysis_group = "6", test = "LR_test", comparison = "Acute_nonAcute", formula = "53")

df_pvalues_AcutenonAcute <- df_pvalues_AcutenonAcute %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_AcutenonAcute, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_AcutenonAcute <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_AcutenonAcute), by =1 ), df_pvalues_AcutenonAcute)

hist(df_pvalues_AcutenonAcute$p_value, breaks =100, main = "Histogram of pvalues - Acute/nonAcute")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_AcutenonAcute, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ Acute_nonAcute + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + Sex + BMI")

lrt_results_age6_1 <- perform_lrt(
  data = df_subset6,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age6_1 <- data.frame(biomolecule_id = unique(df_subset6$biomolecule_id), analysis_group = "6", test = "LR_test", comparison = "Age", formula = "54")

df_pvalues_age6_1 <- df_pvalues_age6_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age6_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age6_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age6_1), by =1 ), df_pvalues_age6_1)

hist(df_pvalues_age6_1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age6_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + Sex + BMI")


lrt_results_sex6_1 <- perform_lrt(
  data = df_subset6,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex6_1 <- data.frame(biomolecule_id = unique(df_subset6$biomolecule_id), analysis_group = "6", test = "LR_test", comparison = "Sex", formula = "55")

df_pvalues_sex6_1 <- df_pvalues_sex6_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex6_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex6_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex6_1), by =1 ), df_pvalues_sex6_1)

hist(df_pvalues_sex6_1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex6_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ Acute_nonAcute + Age + Sex + BMI")

lrt_results_bmi6_1 <- perform_lrt(
  data = df_subset6,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)

df_pvalues_bmi6_1 <- data.frame(biomolecule_id = unique(df_subset6$biomolecule_id), analysis_group = "6", test = "LR_test", comparison = "BMI", formula = "56")

df_pvalues_bmi6_1 <- df_pvalues_bmi6_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi4_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi6_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi6_1), by =1 ), df_pvalues_bmi6_1)

hist(df_pvalues_bmi6_1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi6_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 4 Acute/nonAcute ##########################

# Create individual plots
plot_AcutenonAcute <- ggplot(df_pvalues_AcutenonAcute, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "black") + 
  ggtitle("Acute/nonAcute")

plot_age6_1 <- ggplot(lrt_results_age6_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex6_1 <- ggplot(df_pvalues_sex6_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi6_1 <- ggplot(df_pvalues_bmi6_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")


# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group6_AcutenonAcute_noBatch1.pdf", width = 12, height = 10)
grid.arrange(plot_AcutenonAcute, plot_age6_1, plot_sex6_1, plot_bmi6_1, ncol = 2)
dev.off()





########################## PASCvsnoPASC - Group 7 ###############################

##### P-values for PASCnoPASC ###### 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + Sex + BMI")


lrt_results_PASCnoPASC <- perform_lrt(
  data = df_subset7,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "group7_PASCnoPASC"
)

df_pvalues_PASCnoPASC <- data.frame(biomolecule_id = unique(df_subset7$biomolecule_id), analysis_group = "7", test = "LR_test", comparison = "group7_PASCnoPASC", formula = "57")

df_pvalues_PASCnoPASC <- df_pvalues_PASCnoPASC %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_PASCnoPASC, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_PASCnoPASC <- cbind(pvalue_id = row.names(df_pvalues_PASCnoPASC), df_pvalues_PASCnoPASC)

hist(df_pvalues_PASCnoPASC$p_value, breaks =700, main = "Histogram of pvalues - Group 7 PASC no PASC")

#### Establish a connection to the DB
con <- connect_db(db_path)
#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_PASCnoPASC, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


####### P-values w/ age ######

formula_null <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + Sex + BMI")


lrt_results_age7_1 <- perform_lrt(
  data = df_subset7,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Age"
)

df_pvalues_age7_1 <- data.frame(biomolecule_id = unique(df_subset7$biomolecule_id), analysis_group = "7", test = "LR_test", comparison = "Age", formula = "58")

df_pvalues_age7_1 <- df_pvalues_age7_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_age7_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_age7_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age7_1), by =1 ), df_pvalues_age7_1)

hist(df_pvalues_age7_1$p_value, breaks = 100, main = 'Histogram of pvalues - Age')


## Append age Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_age7_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



##### P-values for Sex ###### 

formula_null <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + BMI")
formula_test <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + Sex + BMI")


lrt_results_sex7_1 <- perform_lrt(
  data = df_subset7,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "Sex"
)


df_pvalues_sex7_1 <- data.frame(biomolecule_id = unique(df_subset7$biomolecule_id), analysis_group = "7", test = "LR_test", comparison = "Sex", formula = "59")

df_pvalues_sex7_1 <- df_pvalues_sex7_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_sex7_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_sex7_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_sex7_1), by =1 ), df_pvalues_sex7_1)

hist(df_pvalues_sex7_1$p_value, breaks = 100, main = 'Histogram of pvalues - Sex')


## Append Sex Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_sex7_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 



####### P-values w/ BMI ######

formula_null <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + Sex")
formula_test <- as.formula("normalized_abundance ~ group7_PASCnoPASC + Age + Sex + BMI")

lrt_results_bmi7_1 <- perform_lrt(
  data = df_subset7,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "BMI"
)


df_pvalues_bmi7_1 <- data.frame(biomolecule_id = unique(df_subset7$biomolecule_id), analysis_group = "7", test = "LR_test", comparison = "BMI", formula = "60")

df_pvalues_bmi7_1 <- df_pvalues_bmi7_1 %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_bmi7_1, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )

df_pvalues_bmi7_1 <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_bmi7_1), by =1 ), df_pvalues_bmi7_1)

hist(df_pvalues_bmi7_1$p_value, breaks = 100, main = 'Histogram of pvalues - BMI')

## Append _bmi Pvalues to db 
## Establish a connection to the DB 
con <- connect_db(db_path)

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_bmi7_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 


########################## Saving Histograms - Group 1 PascnoPASC ##########################

# Create individual plots
plot_PAScnoPASC <- ggplot(df_pvalues_PASCnoPASC, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "maroon") + 
  ggtitle("Group7 PASC_noPASC")

plot_age7_1 <- ggplot(df_pvalues_age7_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "grey") + 
  ggtitle("Age")

plot_sex7_1 <- ggplot(df_pvalues_sex7_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#D7BDE2") + 
  ggtitle("Sex")

plot_bmi7_1 <- ggplot(df_pvalues_bmi7_1, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "#C5E1A5") + 
  ggtitle("BMI")


# Combine plots into one PDF
pdf("P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/Data_Analysis/Linear Regression/Combined_Pvalues_Histograms_Group7_PASCnoPASC_noBatch1Lipids.pdf", width = 12, height = 10)
grid.arrange(plot_PAScnoPASC, plot_age7_1, plot_sex7_1, plot_bmi7_1, ncol = 2)
dev.off()


