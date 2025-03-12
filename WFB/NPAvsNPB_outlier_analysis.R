
## Comparison of NPA or NPB and see if the batch effect holds with just one or both ----


# files
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT standardized_name, rawfile_id, biomolecule_id, raw_abundance, normalized_abundance
                         FROM proteomics_measurement')
biomolecules <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, omics_id, keep 
                           FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_name, Sample, sample_id, ome_id, keep , rawfile_id, run_type
                           FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT Sample, sample_id, Cohort, Age, Sex, BMI, `SF.36.QOL.Score`, PASC_Cohort, Paired_samples
                           FROM patient_metadata')

dbDisconnect(con)



# Make a big dataframe joining only proteomics together
df <- proteomics %>%
  left_join(rawfiles %>%
              dplyr::select(-keep) %>%
              filter(ome_id == 1) %>%
              filter(grepl('Sample', run_type)),
            by = "rawfile_id") %>%
  left_join(metadata %>%
              dplyr::select(-Sample),
            by = "sample_id") %>%
  left_join(biomolecules %>%
              filter(omics_id == 1) %>%
              dplyr::select(biomolecule_id, keep),
            by = "biomolecule_id")
  
filtered_df <- df %>%
  filter(keep == "1")

length(unique(df$standardized_name))
length(unique(filtered_df$standardized_name))



# Make a heatmap for each of these separately (or protein group numbers)
NPA_df <- df %>%
  filter(grepl("NPA", rawfile_name))
NPB_df <- df %>%
  filter(grepl("NPB", rawfile_name))


#### run order vs protein number ----
run_ids_NPA <- NPA_df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  inner_join(df %>% ungroup() %>% select(sample_id, rawfile_id, Cohort), by = 'sample_id') %>%
  distinct() %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())

run_ids_NPB <- NPB_df %>%
  group_by(sample_id) %>%
  summarize(count = sum(!is.na(raw_abundance))) %>%
  inner_join(df %>% ungroup() %>% select(sample_id, rawfile_id, Cohort), by = 'sample_id') %>%
  distinct() %>%
  group_by(sample_id) %>%
  slice_min(rawfile_id, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rawfile_id) %>%
  mutate(order = row_number())

run_ids_combined <-run_ids_NPA %>%
  dplyr::select(sample_id, count) %>%
  left_join(run_ids_NPB,
            by = "sample_id")
  



ggplot(run_ids_combined, aes(count.x, count.y, color = Cohort)) + 
  geom_point() + 
  scale_color_manual(values = pal) +
  xlab('NPA') +
  ylab('NPB') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 8000)) +
  guides(fill = guide_legend(title = 'Study Group')) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(1, 0.05), 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "in")
  ) 
ggsave('reports/figures/AllPlates_Sample_NPB_PGvsRunOrder_Set_Point.pdf', 
       width = 16, height = 6, units = 'cm')




## Heatmap with NAs ----
expression_matrix <- NPA_df %>%
  dplyr::select(standardized_name, raw_abundance, sample_id) %>%
  mutate(raw_abundance = log2(raw_abundance)) %>%
  pivot_wider(names_from = sample_id, values_from = raw_abundance) %>%
  tibble::column_to_rownames(var = "standardized_name")

#Make annotation dataframe for the sample groups
sample_annot <- NPA_df %>%
  ungroup() %>%
  dplyr::select(sample_id, Cohort, Age, Sex, BMI, SF.36.QOL.Score) %>%
  distinct() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(SF.36.QOL.Score = as.numeric(SF.36.QOL.Score))

clust <- expression_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) # column mean

pheatmap(expression_matrix,
         clustering_distance_rows = dist(clust), 
         clustering_distance_cols = dist(t(clust)),
         #color = viridis(24, direction = 1, option = "plasma"),
         color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
         na_col = "white",
         breaks = c(-11, -9, -7, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11),
         cluster_rows = T,
         cluster_cols = T,
         treeheight_row = 0,
         treeheight_col = 10,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         #annotation_row = heatmap_annot,
         annotation_col = sample_annot,
         annotation_colors = list(Cohort = c(Acute = pal[1], Acute_fu = pal[2], Acute_NC = pal[3], Healthy = pal[4], PASC = pal[5], PASC_fu = pal[6])),
         filename = "reports/figures/Unfiltered_NPA_heatmap_Cohort_withNAs.png",
         width = 8,
         height = 8)
