
##########  REFERENCE SETS   #########
make_reference_sets <- function(reference_index, id_index, split_str = "; "){
  all_reference <- NULL
  for(i in 1:length(reference_index)){
    all_reference <- append(all_reference, strsplit(reference_index[i], split_str)[[1]])
  }
  
  unique_reference <- as.list(unique(all_reference), stringsAsFactors = F)
  
  reference_sets <- lapply(unique_reference, function(x) id_index[grep(x[1], reference_index, fixed= T)])
  names(reference_sets) <- unique_reference
  
  reference_sets
}


############ Function for testing enrichment ############## 

enrichment <- function(set, reference_sets, background){
  # Initialize output dataframe
  output <- data.frame(reference = names(reference_sets), 
                       pvalue = rep(1, length(names(reference_sets))), 
                       fdr_pvalue = rep(NA, length(names(reference_sets))),
                       enrichment_ratio = rep(NA, length(names(reference_sets))),
                       stringsAsFactors = FALSE)
  
  # Loop through each reference set
  for (i in 1:nrow(output)){
    # Calculate observed hits (features in 'set' present in the reference set)
    observed_hits <- length(intersect(set, reference_sets[[output[i, "reference"]]]))
    
    # Calculate expected hits by chance (total features in reference set * total features in 'set' / total features in background)
    expected_hits <- length(reference_sets[[output[i, "reference"]]]) * length(set) / length(background)
    
    # Calculate enrichment ratio
    enrichment_ratio <- observed_hits / expected_hits
    
    # Calculate p-value using hypergeometric distribution
    p_value <- phyper(observed_hits - 1, length(reference_sets[[output[i, "reference"]]]), 
                      length(background) - length(reference_sets[[output[i, "reference"]]]), 
                      length(set), lower.tail = FALSE)
    
    # Adjust p-value for multiple testing using Benjamini-Hochberg method
    fdr_p_value <- p.adjust(p_value, method = "BH")
    
    # Update output dataframe with results
    output[i, "pvalue"] <- p_value
    output[i, "fdr_pvalue"] <- fdr_p_value
    output[i, "enrichment_ratio"] <- enrichment_ratio
  }
  
  # Return the output dataframe
  return(output)
} 
