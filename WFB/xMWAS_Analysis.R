
# load package
library(xMWAS)


# files #
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')

proteomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                               FROM proteomics_measurement')
lipidomics <- dbGetQuery(con, 'SELECT biomolecule_id, standardized_name, rawfile_id, raw_abundance, normalized_abundance
                              FROM lipidomics_measurements')
transcriptomics <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name, sample_id, Counts, normalized_counts
                                    FROM rnaseq_measurements")
biomolecules <- dbGetQuery(con, 'SELECT *
                                 FROM biomolecules')
rawfiles <- dbGetQuery(con, 'SELECT rawfile_id, rawfile_name, Sample, sample_id, run_type, ome_id, keep, batch
                             FROM rawfiles_all')
metadata <- dbGetQuery(con, 'SELECT *
                             FROM patient_metadata')

dbDisconnect(con)



## Make proteomics/lipidomics/transcriptomics dataframes ----

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


# lipidomics
rawfiles_l <- rawfiles %>%
  filter(ome_id == 2) %>%
  filter(keep == "1") %>%
  filter(run_type == "Sample") %>%
  select(-keep)

df_l <- lipidomics %>%
  inner_join(rawfiles_l, by = 'rawfile_id') %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_l <- biomolecules %>%
  filter(omics_id == 2) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_l <- df_l %>%
  filter(biomolecule_id %in% biomolecules_l)


# transcriptomics
df_t <- transcriptomics %>%
  filter(!is.na(sample_id)) %>%
  inner_join(metadata, by = 'sample_id')

biomolecules_t <- biomolecules %>%
  filter(omics_id == 3) %>%
  filter(keep == "1") %>%
  pull(biomolecule_id)

filtered_df_t <- df_t %>%
  filter(biomolecule_id %in% biomolecules_t) %>%
  rename(normalized_abundance = normalized_counts)






#### try first for PASC vs no PASC ####
# find overlapping samples
p_ids <- filtered_df_p %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  pull(sample_id) %>%
  unique()
l_ids <- filtered_df_l %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  pull(sample_id) %>%
  unique()
t_ids <- filtered_df_t %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  pull(sample_id) %>%
  unique()

overlap <- as.character(intersect(intersect(p_ids, l_ids), t_ids))



xMWAS_p <- filtered_df_p %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  select(sample_id, normalized_abundance, standardized_name) %>%
  pivot_wider(names_from = sample_id,
              values_from = normalized_abundance) %>%
  select(standardized_name, overlap)

xMWAS_l <- filtered_df_l %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  select(sample_id, normalized_abundance, standardized_name) %>%
  pivot_wider(names_from = sample_id,
              values_from = normalized_abundance) %>%
  select(standardized_name, overlap)

xMWAS_t <- filtered_df_t %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  select(sample_id, normalized_abundance, standardized_name) %>%
  pivot_wider(names_from = sample_id,
              values_from = normalized_abundance) %>%
  select(standardized_name, overlap)

classlabels <- patient_metadata %>%
  filter(Cohort %in% c("PASC", "Acute_fu", "Healthy")) %>%
  filter(PG_change_collection_cutoff == 0) %>%
  select(sample_id, Cohort) %>%
  rename(SampleID = sample_id) %>%
  rename(Class = Cohort) %>%
  mutate(Class = if_else(Class %in% c("Acute_fu", "Healthy"), "noPASC", Class))

write_csv(xMWAS_p, "data/processed/xMWAS_p.csv")
write_csv(xMWAS_l, "data/processed/xMWAS_l.csv")
write_csv(xMWAS_t, "data/processed/xMWAS_t.csv")
write_csv(classlabels, "data/processed/classlabels.csv")

launchShinyApp()

#example dataset that includes metabolome, transcriptome, and cytokine data from the H1N1 mice study (Chandler 2016)
#data(exh1n1)
#data(classlabels) #example classlabels file for case vs control design
#data(classlabels_repeatmeasures) #example classlabels file for repeat measures design

# xMat<-exh1n1$metabolome
# yMat<-exh1n1$transcriptome
# zMat<-exh1n1$cytokine
# wMat<-NA
# classlabels<-exh1n1$classlabels

#Code for reading tab-delimited text files as input data
#currently turned off:
# if(FALSE)
# {
#   fname1<-"/Users/karanuppal/Downloads/OneDrive_1_11-3-2017/gene.txt"
#   fname2<-"/Users/karanuppal/Downloads/OneDrive_1_11-3-2017/clinical.txt"
#   fname3<-"/Users/karanuppal/Downloads/OneDrive_1_11-3-2017/metabolomics.txt"
#   class_fname<-"/Users/karanuppal/Downloads/OneDrive_1_11-3-2017/Classfile.txt"
#   xMat<-read.table(fname1,sep="\t",header=TRUE,row.names=1)
#   yMat<-read.table(fname2,sep="\t",header=TRUE,row.names=1)
#   zMat<-read.table(fname3,sep="\t",header=TRUE,row.names=1)
#   classlabels<-read.table(class_fname,sep="\t",header=TRUE)
#   xMat<-as.data.frame(xMat)
#   yMat<-as.data.frame(yMat)
#   zMat<-as.data.frame(zMat)
#   wMat<-NA
# }
###################

output <- "" #change for your computer

#Please see user manual for description of arguments:
#https://github.com/kuppal2/xMWAS/blob/master/example_manual_tutorial/xMWAS-manual.pdf

#call the run_xmwas() function:
xmwas_res <- run_xmwas(Xome_data = xMWAS_p, Yome_data = xMWAS_l, Zome_data = xMWAS_t, Wome_data = NA, outloc = "xMWAS_output",
                     classlabels = NA, class_fname = NA, xmwasmethod = "pls", plsmode = "regression",
                     max_xvar=10000, #e.g. select top 10000 of the variabels in X dataset based on relative standard deviation; change according to your dataset; you can also use proportion such as round(nrow(xMat)*0.3) to select top 30% of the variables.
                     max_yvar=10000, #select top 10000 of the variabels in Y dataset based on relative standard deviation;  change according to your dataset; you can also use proportion such as round(nrow(yMat)*0.3) to select top 30% of the variables.
                     max_zvar=10000, #select top 10000 variabels in Z dataset based on relative standard deviation;  change according to your dataset; you can also use proportion such as round(nrow(zMat)*0.3) to select top 30% of the variables.
                     max_wvar=10000, #select top 10000 variabels in W dataset based on relative standard deviation;  change according to your dataset; you can also use proportion such as round(nrow(wMat)*0.3) to select top 30% of the variables.
                     rsd.filt.thresh=1,
                     corthresh=0.4, #absolute correlation threshold
                     keepX=1000, #select up to top 1000 variables in the sPLS model; change according to your dataset
                     keepY=1000, #select up to top 1000 variables in the sPLS model; change according to your dataset
                     keepZ=1000, #select up to top 1000 variables in the sPLS model; change according to your dataset
                     keepW=1000, #select up to top 1000 variables in the sPLS model; change according to your dataset
                     pairedanalysis=FALSE, #set to TRUE if repeated measures study design
                     optselect=FALSE, #perform optimal PLS componenet selection; TRUE or FALSE; set to FALSE for exact Pearson correlation calculation using PLS regression
                     rawPthresh=0.05, #p-value threshold for correlation based on Student's t-test
                     numcomps=5, #max number of PLS components to use; set to N-1 (N: number of samples) for exact Pearson correlation calculation using PLS regression
                     net_edge_colors=c("blue","red"),
                     net_node_colors=c("orange", "green","cyan","pink"),
                     Xname="Protein", #change the name of dataset X
                     Yname="Lipid", #change the name of dataset Y
                     Zname="Transcript", #change the name of dataset Z
                     Wname="W", #change the name of dataset W
                     net_node_shape=c("square","circle","triangle","star"),
                     all.missing.thresh=0, #filter based on missing values: set to NA to turn it OFF; otherwise specify a value between: 0 to 1 (e.g. 0.8 to require that at least 80% of the samples have a non-missing value)
                     missing.val=0,
                     seednum=100,label.cex=0.2,vertex.size=6,
                     interactive=FALSE,max_connections=NA,
                     centrality_method="eigenvector", #centrality evaluation method
                     use.X.reference=FALSE,removeRda=TRUE,
                     compare.classes=FALSE, #compare classes: TRUE or FALSE
                     class.comparison.allvar=TRUE,
                     modularity.weighted=TRUE,
                     globalcomparison=TRUE,
                     plot.pairwise=FALSE, #plot results for pairwise comparisons: TRUE or FALSE
                     apply.sparse.class.comparison=FALSE, #perform variable selection in sPLS during class-wise comparison (default: FALSE)
                     layout.type="fr1")

suppressWarnings(try(sink(file=NULL),silent=TRUE))



