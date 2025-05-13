"
    This method answers the question, What is the smallest set of GO terms that will cover, some percent
        of your proteins, while avoiding using the overly broad, high-level GO terms with hundreds to thousands
        of genes in them?

    Method is useful for making a high-level overview of your proteins by separating them into groups that have minimal overlap.
    Also useful for selecting the best set of GO terms to plot after doing GSEA or other enrichment.

    Uses a greedy algorithm that checks for the next GO term that will include the largest number of features not yet included.

    Returns:
    A 4-Tuple of (
        The list of GO terms selected,
        The set of features covered by these GO terms,
        The intersection of features each combo of GO terms,
        The size of the intersection for each combo,
    )

    params:

    features_to_go_terms_dict: A dictionary of unique feature ID to a list (or a semicolon-separated string) of GO terms.

    min/max_features: Only include GO terms between these numbers of features

    target_coverage: float between (0.0, 1.0] indicating % of coverage it will aim to reach

    required_go_terms, go_terms_to_leave_out: force inclusion/exclusion of certain GO terms

    max_overlap: The maximum number of features that are allowed to overlap when a new GO term is added.
    Setting this to 0 means that each GO term will be totally in the features it covers, but you likely cannot reach your target_coverage percentage this way.

    include_bp/cc/mf: Whether to use terms if they fall under one of the 3 high-level namespaces:
        Biological Process, Cellular Component, Molecular Function.
"


## R Version ----
# read tables from database
con <- dbConnect(RSQLite::SQLite(), dbname = 'P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Database/Long Covid Study DB.sqlite')
biomolecules_metadata <- dbGetQuery(con, 'SELECT *
                                          FROM biomolecules_metadata')
GO_terms <- dbGetQuery(con, 'SELECT *
                             FROM GO_terms')
dbDisconnect(con)


# features_to_go_terms
f_t_g_t <- biomolecules_metadata %>%
  filter(metadata_type == "GO_terms") %>%
  select(biomolecule_id, metadata_value) %>%
  filter(biomolecule_id %in% biomolecules_p) # Change how many items to search here
f_t_g_t <- split(f_t_g_t$metadata_value, f_t_g_t$biomolecule_id)
  
# go_term_ontology
g_t_o <- GO_terms %>%
  select(GO_term, namespace) %>%
  deframe()


find_minimum_set_of_go_terms <- function(
    features_to_go_terms,
    min_features,
    max_features,
    target_coverage,
    required_go_terms = NULL,
    go_terms_to_leave_out = NULL,
    max_overlap = 3,
    include_bp = TRUE,
    include_cc = TRUE,
    include_mf = TRUE,
    go_term_ontology = NULL # named vector: GO ID -> namespace ("biological_process", etc.)
) {
  if (is.null(go_term_ontology)) stop("Provide `go_term_ontology` as a named vector: GO ID -> ontology")
  
  # Build GO to features map
  go_to_features <- list()
  for (i in seq_along(features_to_go_terms)) {
    feature <- names(features_to_go_terms)[i]
    go_terms <- features_to_go_terms[[i]]
    if (is.character(go_terms)) {
      go_terms <- unlist(strsplit(go_terms, ";"))
    }
    for (go in go_terms) {
      if (go != "") {
        if (!go %in% names(go_to_features)) {
          go_to_features[[go]] <- character(0)
        }
        go_to_features[[go]] <- unique(c(go_to_features[[go]], feature))
      }
    }
  }
  
  # Filter by size
  filtered_go_terms <- go_to_features[sapply(go_to_features, function(x) {
    length(x) >= min_features && length(x) <= max_features
  })]
  
  # Exclude specific GO terms
  if (!is.null(go_terms_to_leave_out)) {
    filtered_go_terms <- filtered_go_terms[!names(filtered_go_terms) %in% go_terms_to_leave_out]
  }
  
  # Filter by ontology
  ontology_filter <- sapply(names(filtered_go_terms), function(go) {
    ns <- go_term_ontology[go]
    if (is.na(ns)) return(FALSE)
    if (ns == "biological_process") return(include_bp)
    if (ns == "cellular_component") return(include_cc)
    if (ns == "molecular_function") return(include_mf)
    return(FALSE)
  })
  filtered_go_terms <- filtered_go_terms[ontology_filter]
  
  if (length(filtered_go_terms) == 0) stop("No GO terms remaining after filtering.")
  
  covered_features <- character(0)
  required_features <- ceiling(length(features_to_go_terms) * target_coverage)
  selected_go_terms <- character(0)
  
  if (!is.null(required_go_terms)) {
    for (go in required_go_terms) {
      covered_features <- unique(c(covered_features, go_to_features[[go]]))
      selected_go_terms <- c(selected_go_terms, go)
    }
  }
  
  while (length(covered_features) < required_features && length(filtered_go_terms) > 0) {
    # Filter by overlap
    filtered_go_terms <- filtered_go_terms[sapply(filtered_go_terms, function(features) {
      length(intersect(features, covered_features)) <= max_overlap
    })]
    
    if (length(filtered_go_terms) == 0) {
      message("No more GO terms satisfy the max_overlap constraint.")
      break
    }
    
    # Select GO term with most new features
    best_go <- names(which.max(sapply(filtered_go_terms, function(features) {
      length(setdiff(features, covered_features))
    })))
    selected_go_terms <- c(selected_go_terms, best_go)
    covered_features <- unique(c(covered_features, filtered_go_terms[[best_go]]))
    filtered_go_terms[[best_go]] <- NULL
  }
  
  # Calculate intersections
  # intersections <- list()
  # intersection_sizes <- list()
  # for (go1 in selected_go_terms) {
  #   intersections[[go1]] <- list()
  #   intersection_sizes[[go1]] <- list()
  #   for (go2 in selected_go_terms) {
  #     inter <- intersect(go_to_features[[go1]], go_to_features[[go2]])
  #     intersections[[go1]][[go2]] <- inter
  #     intersection_sizes[[go1]][[go2]] <- length(inter)
  #   }
  # }
  
  return(list(
    selected_go_terms = selected_go_terms,
    covered_features = covered_features
    # intersections = intersections,
    # intersection_sizes = intersection_sizes
  ))
}



## run function
min_go_terms <- find_minimum_set_of_go_terms(
  features_to_go_terms = f_t_g_t,
  min_features = 5,
  max_features = 50,
  target_coverage = 0.9,
  max_overlap = 5,
  go_term_ontology = g_t_o,
  include_cc = FALSE
)


GO_set <- GO_terms %>%
  filter(GO_term %in% min_go_terms$selected_go_terms)
  

fwrite(GO_set,
       "data/metadata/GOtermset_PASCnoPASC_Proteins_5to50_90coverage_5overlap.csv")
