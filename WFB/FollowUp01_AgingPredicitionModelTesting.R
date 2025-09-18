
################## organize proteomics data for testing aging predictive model ######################


#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(RSQLite)
library(data.table)
library(glmnet)


# Colors #
col <- brewer.pal(8, "Dark2") 

pal1 <- c("#66C2A5",
          "#FFD92F",
          "#8DA0CB",
          "#FC8D62",
          "#A6CEE3",
          "#E78AC3",
          "#A6D854",
          "#FDB462",
          "#B3B3B3",
          "#B2DF8A")

pal <- c('#EE6677', 
         '#AA3377', 
         '#CCBB44', 
         '#228833', 
         '#66CCEE', 
         '#4477AA')


pal <- c("Acute" = "#E78AC3", 
         "Acute_fu" = '#AA3377', 
         "Acute_NC" = "#B3B3B3", 
         "Healthy" = '#229100', 
         "PASC" = '#66CCEE', 
         "PASC_fu" = '#4477AA')


col1 <- viridis_pal(option = "rocket")(100)[round(c(0.25, 0.5, 0.75) * 100)]


# plot colors
pie(rep(1, length(col)), col = col , main="") 



# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT *
                               FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT *
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                             FROM biomolecules_metadata')

dbDisconnect(con)



# proteomics
# Merge rawfiles and proteomics, Combine NPA and NPB measurements by completeness by protein group
# subset rawfiles to only include sample proteomics runs
rawfiles_p <- rawfiles %>%
  filter(ome_id == 1) %>%
  select(-keep) %>%
  filter(grepl('Sample', run_type))

df_p <- proteomics %>%
  inner_join(rawfiles_p, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_p <- biomolecules %>%
  filter(omics_id == 1) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_p <- df_p %>%
  filter(biomolecule_id %in% biomolecules_p) 




## make correct format for python lasso regression check ----
filtered_df_p_python <- filtered_df_p %>%
  filter(PG_change_collection_cutoff == 0) %>% # filter out post-tube change samples for strange behavior
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "gene_name") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id") %>% # add gene name for proper organization
  select(sample_id, Age, Cohort, metadata_value, normalized_abundance) %>%
  group_by(metadata_value) %>%
  mutate(zscore = (normalized_abundance - mean(normalized_abundance)) / sd(normalized_abundance)) %>% # make z score column
  ungroup() %>%
  select(-normalized_abundance) %>%
  pivot_wider(id_cols = c(sample_id, Age, Cohort),
              names_from = metadata_value,
              values_from = zscore)

write_csv(filtered_df_p_python, "data/processed/ProteinsForAgingPredictionModel.csv")




## make correct format for coefficients ----
coef <- read_csv("data/processed/Lehallier2019_ModelCoefficients.csv") %>%
  mutate(Variable = word(Variable, 1, 1, "\\."))
intercept <- coef$Coefficients[1]

write_csv(coef, "data/processed/Lehallier2019_ModelCoefficients_fixed.csv")



## do modeling with coefficients ----
coef <- coef[-1,]

##* match proteins in model to my data ----
model_proteins <- coef$Variable
available_proteins <- intersect(model_proteins, colnames(filtered_df_p_python))

cat("Model uses", length(model_proteins), "proteins, found", length(available_proteins), "in your data\n")

# Subset coeffs and data
coefs <- coef %>% 
  filter(Variable %in% available_proteins) %>%
  group_by(Variable) %>%
  slice_max(order_by = Coefficients, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(ScaledCoefs = Coefficients * 3) #scale coefficients for my data

X <- filtered_df_p_python[, available_proteins, drop = FALSE]


##* apply the model ----
# make named vector of coefficients and make sure X is in the same order
beta <- setNames(coefs$ScaledCoefs, coefs$Variable)
X <- X[, names(beta), drop = FALSE]

# calculate predicted age using matrix multiplication
pred_age <- intercept + as.numeric(as.matrix(X) %*% beta)

# add to data
filtered_df_p_python$Predicted_Age <- pred_age


##* check model performance ----
true_age <- filtered_df_p_python$Age

r2_val <- summary(lm(pred_age ~ true_age))$r.squared
cor_test <- cor.test(true_age, pred_age)

# Print performance
cat("Model Performance:\n",
    sprintf(" R²: %.3f\n Pearson r: %.3f (p=%.2e)\n",
            r2_val, cor_test$estimate, cor_test$p.value))


##* plot pred vs true age ----
ggplot(filtered_df_p_python,
       aes(Age, Predicted_Age, color = Cohort)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  ylim(16,94) +
  xlim(16,94) +
  labs(x = "True Age",
       y = "Predicted Age")


## Make my own LASSO model ----
# Use all 6972 proteins and 330 samples from my data
filtered_df_p_R <- filtered_df_p %>%
  filter(PG_change_collection_cutoff == 0) %>% # filter out post-tube change samples for strange behavior
  left_join(biomolecules_metadata %>%
              select(-metadata_id) %>%
              filter(metadata_type == "gene_name") %>%
              select(biomolecule_id, metadata_value),
            by = "biomolecule_id") %>% # add gene name for proper organization
  select(sample_id, Age, Cohort, metadata_value, normalized_abundance) %>%
  group_by(metadata_value) %>%
  mutate(zscore = (normalized_abundance - mean(normalized_abundance)) / sd(normalized_abundance)) %>% # make z score column
  ungroup() %>%
  select(-normalized_abundance)

allproteins <- unique(filtered_df_p_R$metadata_value)

# protein matrix
filtered_df_p_R_matrix <- filtered_df_p_R %>%
  pivot_wider(names_from = metadata_value,
              values_from = zscore)


##* do glmnet fit ----
cvfit <- cv.glmnet(x = as.matrix(filtered_df_p_R_matrix[,allproteins]),
                   y = filtered_df_p_R_matrix$Age,
                   alpha = 0.5,
                   type.measure = "mse", 
                   nfolds = 20)
print(cvfit)
plot(cvfit)

coefs_glmnet <- as.data.frame(as.matrix(coef.glmnet(cvfit, s = "lambda.min"))) %>%
  rownames_to_column(var = "coef") %>%
  filter(lambda.min != 0)


##* check predicted vs true age with coefficients
# set up named vector of coefficients and intercept
intercept <- coefs_glmnet %>%
  filter(coef == "(Intercept)") %>%
  pull(lambda.min)
coefs_glmnet <- coefs_glmnet[-1,]
beta <- setNames(coefs_glmnet$lambda.min,
                 coefs_glmnet$coef)

# set up data matrix
X <- as.matrix(filtered_df_p_R_matrix[,names(beta)])

# calculate predicted age using matrix multiplication
pred_age <- intercept + as.numeric(X %*% beta)

# add to data
filtered_df_p_R_matrix$Predicted_Age <- pred_age



##* plot pred vs true age ----
ggplot(filtered_df_p_R_matrix,
       aes(Age, Predicted_Age, color = Cohort)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  ylim(16,94) +
  xlim(16,94) +
  labs(x = "True Age",
       y = "Predicted Age")


##* plot top coefficients vs age ----
coc <- "SOST"

single_prot_age_plot <- filtered_df_p_R %>%
  filter(metadata_value == coc)

ggplot(single_prot_age_plot,
       aes(Age, zscore, color = Cohort)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  labs(x = "Age",
       y = "z-score",
       title = paste(coc, "correlation with age"))


## check for agreement with my model and the other ----

##* shared coefficients ----
shared <- intersect(coef$Variable,
          coefs_glmnet$coef)
shared
# 17 shared proteins

mine_filtered <- coefs_glmnet %>%
  filter(coef %in% shared)
other_filtered <- coef %>%
  filter(Variable %in% shared) %>%
  group_by(Variable) %>%
  slice_max(order_by = Coefficients, n = 1, with_ties = FALSE) %>%
  ungroup()

plot(mine_filtered$lambda.min ~ other_filtered$Coefficients)
abline(h = 0, col = "red", lty = 2, lwd = 2)
abline(v = 0, col = "red", lty = 2, lwd = 2)
