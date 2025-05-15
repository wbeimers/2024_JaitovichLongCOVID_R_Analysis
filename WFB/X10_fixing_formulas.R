
## Fix Formula Tables ----
# 1. add 51 and 52 to the formulas table
# 2. replace 43 and 44 with 51 and 52 in the machine_learning_coefficients




# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
formulas <- dbGetQuery(con, 'SELECT *
                           FROM formula_table')
ml <- dbGetQuery(con, 'SELECT *
                           FROM machine_learning_coefficients')
dbDisconnect(con)

# update formulas
upd <- data.frame(formula = c(51,52),
                  test = c("machine_learning_elastic_net", "machine_learning_L1_logistic_regression"),
                  analysis_group = c(1,1),
                  comparison = c("predict_QoL", "predict_PASC_vs_non-PASC"),
                  formula_null = rep(NA_character_, 2),
                  formula_test = c("all biomolecules + Age + Sex + BMI", "all biomolecules + Age + Sex + BMI"),
                  Notes = c("No batch1_Lipids", "No batch1_Lipids"))

formulas1 <- formulas %>%
  bind_rows(upd)

# update ml coefficients
ml1 <- ml %>%
  mutate(formula_id = case_when(
    formula_id == 43 ~ 51,
    formula_id == 44 ~ 52
  ))


# update tables in database
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
dbWriteTable(con, "formula_table", formulas1, overwrite = T)
dbWriteTable(con, "machine_learning_coefficients", ml1, overwrite = T)
dbDisconnect(con)
