#### libraries/colors/files ####
# libraries #
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggfortify)
library(missForest)
library(VIM)
library(viridis)
library(RSQLite)
library(ggrepel)
library(pheatmap)
library(ROTS)
library(data.table)
library(enrichR)
library(enrichplot)
setEnrichrSite("Enrichr")
library(fgsea)
library(UniProt.ws)
up <- UniProt.ws(taxId = 9606)


# Colors #
# Make a classic palette
col <- brewer.pal(8, "Set2") 

pal <- c("#66C2A5",
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

# Make a Custom Gradient
col1 <- colorRampPalette(col)(16)

# plot colors
pie(rep(1, length(col)), col = col , main="") 


# Files #
# Groups: Acute, Acute_fu, Acute_NC, Healthy, PASC, PASC_fu
# PASC_Cohort: first, second
group1 <- "first"
group2 <- "second"

all_fc <- fread(paste0("data/processed/AllPlates_Samples_Volcano_", group1, "_", group2, ".csv"))



#### Look at ENRICHMENT of proteins ####
#databases
dbs <- listEnrichrDbs()
dbs <- c("GO_Cellular_Component_2023", "GO_Biological_Process_2023", "GO_Molecular_Function_2023")

#switch uniprot_id to gene_id
#taking from volcano plot DE list
up_fc <- filter(all_fc, grepl("UP", diffexp))
down_fc <- filter(all_fc, grepl("DOWN", diffexp))
#do either up list or down list
uniprot_ids <- up_fc$PG.ProteinGroup
gene_ids <- mapIds(org.Hs.eg.db, keys = uniprot_ids, column = "SYMBOL", keytype = "UNIPROT")

#check enrichment of a set of genes
enriched <- enrichr(gene_ids, dbs)
enriched[["GO_Biological_Process_2023"]]

#plot enrichment
plotEnrich(enriched$GO_Cellular_Component_2023, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GOBP")

ggsave("reports/figures/AllPlates_GOCC_GeneEnrichment_PASCvsAcuteFU_UP.pdf", width = 24, height = 16, units = "cm")



####foldchange-based enrichment####
#data("examplePathways")
#data("exampleRanks")

#fgseaRes <- fgsea(pathways = examplePathways, 
#                  stats    = exampleRanks,
#                  minSize  = 15,
#                  maxSize  = 500)

#subset data by fold change
gene_list <- all_fc
gene_list$symbol <- mapIds(org.Hs.eg.db, keys = gene_list$PG.ProteinGroup, column = "SYMBOL", keytype = "UNIPROT")
gene_list <- gene_list %>%
  mutate(neglogpvalue = if_else(is.infinite(neglogpvalue), max(neglogpvalue[is.finite(neglogpvalue)]), neglogpvalue))
gene_list$rank <- gene_list$logfc*gene_list$neglogpvalue



gene_list <- setNames(gene_list$rank, gene_list$symbol)
gene_list <- gene_list[!is.na(names(gene_list))]


#load pathways file
gmt_file <- "data/metadata/c5.go.v2023.2.Hs.symbols.gmt"
pathways <- gmtPathways(gmt_file)

fgsea <- fgsea(pathways = pathways, 
               stats    = gene_list,
               minSize  = 10,
               maxSize  = 500)

#check top pathways
head(fgsea[order(padj), ])

#significant pathways?
sum(fgsea[, padj < 0.01])

#plot one pathway
plotEnrichment(pathways[["GOBP_COMPLEMENT_ACTIVATION"]],
               gene_list) + labs(title="GOBP_COMPLEMENT_ACTIVATION")


#plot top pathways
topPathwaysUp <- fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], gene_list, fgsea, 
              gseaParam=0.5)

#save fgsea
data.table::fwrite(fgsea, file = "data/processed/Plate_01-05_Samples_fgsea_GOpathways.tsv", sep = "\t", sep2 = c("", " ", "")) 






## Enrichment of Unique Proteins in Each Cohort ----
# Make data frames of unique detections from each cohort with unfiltered df to capture all proteins
cohorts <- unique(df$Cohort)

for (i in cohorts) {
  
  x <- df %>%
    group_by(standardized_name) %>%
    filter(
      any(!is.na(raw_abundance[Cohort == i])) &
        all(is.na(raw_abundance[Cohort != i]))
    )
  
  y <- x %>%
    group_by(standardized_name) %>%
    summarise(non_na_count = sum(!is.na(raw_abundance)))
  
  fwrite(y, paste0("data/processed/AllPlates_Samples_", i, "_uniqueIDs.csv"))

  assign(paste0("unique_df_", i), x)
  assign(paste0("unique_df_counts_", i), y)
  
}

# look at enrichment of unique Acute proteins

gene_list <- unique(df$standardized_name)
mapping <- UniProt.ws::select(up, keys = gene_list, keytype = "UniProtKB", columns = "gene_primary")
background <- unlist(strsplit(mapping$Gene.Names..primary., ";"))
write.csv(background, "data/processed/AllPlates_Proteins_All_background.csv")

gene_list <- unique_df_counts_Acute$standardized_name
mapping <- UniProt.ws::select(up, keys = gene_list, keytype = "UniProtKB", columns = "gene_primary")
input <- unlist(strsplit(mapping$Gene.Names..primary., ";"))
write.csv(input, "data/processed/ApplPlatesProteins_PASCvsHealthy_unique_Input.csv")

enriched <- enrichr(input, dbs, background = background)

enrichr <- enriched$GO_Cellular_Component_2023

enrichr <- enrichr %>%
  mutate(`-log10qvalue` = -log10(Adjusted.P.value)) %>%
  mutate(Rank = as.numeric(enrichr$Rank))

ggplot(enrichr, aes(Odds.Ratio, `-log10qvalue`)) +
  geom_point(shape = 16,
             size = 2,
             aes(color = Rank)) +
  geom_text_repel(aes(label = Term), size = 2) +
  scale_color_viridis(option = "plasma") +
  xlab("Odds Ratio") +
  ylab("-log10(q-value)") +
  #xlim(c(-150, 150)) +
  ylim(c(0, 4)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(color = "gray", fill = NA, size = 0.2))
ggsave("reports/figures/AllPlates_Samples_PASCvsHealthyUnique_ENRICHr_GOCC.pdf", 
       width = 18, height = 12, units = "cm")





