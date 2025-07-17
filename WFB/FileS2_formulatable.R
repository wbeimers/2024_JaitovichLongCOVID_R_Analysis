library(devtools)
library(tidyverse)
library(RSQLite)








## Similar Acute vs PASC features ----
# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

formulas <- dbGetQuery(con, 'SELECT *
                            FROM formula_table')

dbDisconnect(con)


write_csv(formulas, "data/metadata/FileS2_FormulaTable.csv")





