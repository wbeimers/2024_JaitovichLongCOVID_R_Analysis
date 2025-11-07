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

########################### Important Functions ################################

connect_db <- function(db_path) {
  dbConnect(RSQLite::SQLite(), dbname = db_path)
}

query_data <- function(con, table) {
  dbGetQuery(con, paste0("SELECT * FROM ", table))
}



#LRM - fixed model

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

df_subset <- df[df$Cohort == "PASC" & df$PG_change_collection_cutoff == 0, ]


##########################  Calculating pvalues - LRT ################################
################################ Cohort - PASC ################################

##### P-values for Incubation Period ################# 

formula_null <- as.formula("normalized_abundance ~ Age + Sex + BMI")
formula_test <- as.formula("normalized_abundance ~ delta_time_infection_enrollment + Age + Sex + BMI")


lrt_results_Incubation <- perform_lrt(
  data = df_subset,
  formula_null = formula_null,
  formula_test = formula_test,
  predictor = "delta_time_infection_enrollment"
)

df_pvalues_incubation <- data.frame(biomolecule_id = unique(df_subset$biomolecule_id), analysis_group = "8", test = "LRT_test", comparison = "Incubation_Period", formula = "61")

df_pvalues_incubation <- df_pvalues_incubation %>%
  distinct(biomolecule_id, .keep_all = TRUE) %>%
  inner_join(lrt_results_Incubation, by = "biomolecule_id") %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
  )


df_pvalues_incubation <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_incubation), by =1 ), df_pvalues_incubation)

#pvalues histogram
hist(df_pvalues_incubation$p_value, breaks =100, main = "Histogram of pvalues - Incubaion Period")
hist(df_pvalues_incubation$q_value, breaks =100, main = "Histogram of pvalues - Incubaion Period")

#connect to the DB
con <- connect_db(db_path)
# write table to DB 
dbWriteTable(con, "pvalues", df_pvalues_incubation, append = T)
# check
pvalues <- dbReadTable(con, "pvalues")
# disconnect
dbDisconnect(con) 



#update formulas table

#connect to the DB
con <- connect_db(db_path)
#read formulas
formulas <- query_data(con, "formula_table")
colnames(formulas)

# Create the new row as a small data frame
Updatedformulas <- data.frame(
  X = 61,
  formula = 61,
  test = "LRT_test",
  analysis_group = 8,
  comparison = "Incubation_Period",
  formula_null = "\"normalized_abundance ~ Age + Sex + BMI\"",
  formula_test = "\"normalized_abundance ~ delta_time_infection_enrollment + Age + Sex + BMI\"",
  Notes = "No batch1_Lipids; Analysis group 8 includes PASC only (no PASC fu)"
)

# Append it to the existing formulas table
formulas <- bind_rows(formulas, Updatedformulas)

# write table to DB 
dbWriteTable(con, "formulas", formulas, overwrite = T)
# check
newformulas <- dbReadTable(con, "formulas")
# disconnect
dbDisconnect(con) 
