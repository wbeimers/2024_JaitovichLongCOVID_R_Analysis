
library(umap)

pca_set <- filtered_df_pl %>%
  filter(biomolecule_id %in% UP_DOWN_together$biomolecule_id) %>%
  filter(PASCnoPASC %in% c("Acute_COVID", "No_COVID", "Long_COVID")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  dplyr::select(biomolecule_id, normalized_abundance, sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = normalized_abundance) %>%
  select(where(~ !any(is.na(.))))

t_pca_set <- as.data.frame(t(pca_set))

# turn row into colnames
colnames(t_pca_set) <- as.character(t_pca_set[1, ])
# Remove the row
t_pca_set <- t_pca_set[-1, ]
#change to numeric
t_pca_set <- data.frame(lapply(t_pca_set, as.numeric), row.names = rownames(t_pca_set))
# re-add sample_id
t_pca_set <- tibble::rownames_to_column(t_pca_set, "sample_id") %>%
  mutate(sample_id = as.integer(sample_id))



umap_score <- umap(t_pca_set[,c(2:ncol(t_pca_set))])

# pca scores
umap_scores <- as.data.frame(umap_score$layout) %>%
  mutate(sample_id = t_pca_set$sample_id) %>%
  left_join(filtered_df_pl %>% 
              select(sample_id, Cohort, PASCnoPASC, Age, Sex, BMI, unique_patient_id, Collection_date) %>%
              distinct(), 
            by = "sample_id")



# pca plot
ggplot(umap_scores, aes(V1, V2, fill = PASCnoPASC)) + 
  geom_point(shape = 21,
             size = 1.5,
             alpha = 0.9,
             color = "black",
             stroke = 0.1) +
  stat_ellipse(aes(color = PASCnoPASC), 
               geom = "path", 
               show.legend = FALSE,
               linewidth = 0.2) +
  #geom_text_repel(aes(label = Sample), size = 2) +
  #scale_fill_manual(values = pal) +
  scale_fill_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
  scale_color_viridis(discrete = T, option = "rocket", begin = 0.25, end = 0.75) +
  xlab("V1") +
  ylab("V2") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "inside", 
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.25, "cm")
  )

ggsave(paste0("reports/figures/PCA_plots/UMAP_AC_LC_nC_OverlapBiomolecules_proteinlipid_nobatch1_pretube.pdf"), 
       width = 8, height = 5, units = "cm")
