####Libraries/Load Files/colors####
#libraries#
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(enrichR)
library(enrichplot)
setEnrichrSite("Enrichr")
library(fgsea)
library(org.Hs.eg.db)


#Colors#
#Make a classic palette
col <- brewer.pal(10, "Set2") 

#Make a Custom Gradient
col1 <- c(rev(colorRampPalette(col)(100)),"white", colorRampPalette(col1)(100))

#plot colors
pie(rep(1, length(col)), col = col , main="") 


#files#
all_fc <- read.csv("data/processed/Plate_01-05_Samples_AcutevsPASC_DE_list.csv")
all_fc <- all_fc[-1]
colnames(all_fc) <- c("gene", "fc", "t_pv", "p_pv", "t_padj", "p_padj", "log2(fc)", "-log10(t_padj)", "-log10(p_padj)", "diffexp")



####Look at ENRICHMENT of proteins####
#databases
dbs <- listEnrichrDbs()
dbs <- c("GO_Cellular_Component_2023", "GO_Biological_Process_2023", "GO_Molecular_Function_2023")

#switch uniprot_id to gene_id
#taking from volcano plot DE list
up_fc <- filter(all_fc, grepl("UP", diffexp))
down_fc <- filter(all_fc, grepl("DOWN", diffexp))
#do either up list or down list
uniprot_ids <- up_fc$gene
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
gene_list$symbol <- mapIds(org.Hs.eg.db, keys = gene_list$gene, column = "SYMBOL", keytype = "UNIPROT")
gene_list$t_rank <- gene_list$`log2(fc)`*(-log10(gene_list$t_pv))
gene_list$p_rank <- gene_list$`log2(fc)`*(-log10(gene_list$p_pv))

gene_list <- setNames(gene_list$t_rank, gene_list$symbol)
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
plotGseaTable(pathways[topPathwaysUp], gene_list, fgsea, 
              gseaParam=0.5)

#save fgsea
data.table::fwrite(fgsea, file = "data/processed/Plate_01-05_Samples_fgsea_GOpathways.tsv", sep = "\t", sep2 = c("", " ", "")) 

