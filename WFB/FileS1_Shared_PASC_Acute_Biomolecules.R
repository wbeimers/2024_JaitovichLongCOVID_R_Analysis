library(devtools)
library(tidyverse)
library(RSQLite)








## Similar Acute vs PASC features ----
# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

pvalues <- dbGetQuery(con, 'SELECT biomolecule_id, analysis_group, test, comparison, formula, predictor, effect_size, eta_squared, lratio, p_value, q_value
                            FROM pvalues')
formulas <- dbGetQuery(con, 'SELECT *
                             FROM formula_table')
dbDisconnect(con)

# options:
# analysis_group (1, 2, 3, 0)
anal <- 7
# comparison (Age, Sex, QoL, BMI)
comp <- "group7_PASCnoPASC"
# formula (1, 2, 3, etc.)
form <- 57

volc_plot_PASC <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 & effect_size > 0 ~ "UP",
    q_value <= 0.05 & effect_size < 0 ~ "DOWN",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id)


# options:
# analysis_group (1, 2, 3, 0)
anal <- 6
# comparison (Age, Sex, QoL, BMI)
comp <- "Acute_nonAcute"
# formula (1, 2, 3, etc.)
form <- 53

volc_plot_Acute <- pvalues %>%
  filter(analysis_group == anal) %>%
  filter(comparison == comp) %>%
  filter(formula == form) %>%
  inner_join(biomolecules %>%
               select(biomolecule_id, standardized_name, omics_id),
             by = "biomolecule_id") %>%
  mutate(neglogpvalue = -log10(p_value)) %>%
  mutate(effect_size = -(effect_size)) %>%
  mutate(diffexp = case_when(
    q_value <= 0.05 & effect_size > 0 ~ "UP",
    q_value <= 0.05 & effect_size < 0 ~ "DOWN",
    T ~ "NO"
  )) %>% 
  mutate(ome = case_when(
    omics_id == 1 ~ "protein",
    omics_id == 2 ~ "lipid",
    omics_id == 3 ~ "transcript"
  )) %>%
  select(-omics_id) 

# find overlap of UP and DOWN features between both Acute and PASC

UP_df <- volc_plot_PASC %>%
  filter(diffexp == "UP") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "UP"),
             by = "biomolecule_id")

DOWN_df <- volc_plot_PASC %>%
  filter(diffexp == "DOWN") %>%
  inner_join(volc_plot_Acute %>%
               filter(diffexp == "DOWN"),
             by = "biomolecule_id")

UP_DOWN_figS1 <- UP_df %>%
  mutate(direction = "UP") %>%
  bind_rows(DOWN_df %>%
              mutate(direction = "DOWN"))

write_csv(UP_DOWN_figS1, "data/processed/FileS1_shared_biomolecules.csv")
