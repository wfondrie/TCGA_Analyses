library(tidyverse)
library(SummarizedExperiment)

#' @param study A vector of TCGA studies to retrieve data from.
#' @param GOIs A vector of genes of interest. Leaving as NULL returns expression
#' for all genes

load_expr <- function(study, GOIs = NULL) {
  studyFiles <- paste0("GDCdata/TCGA_", study, ".rda")
  
  dat <- lapply(studyFiles, function(f){ load(f); return(data)})
  
  dat <- do.call(cbind, dat)
  
  # Format Expression Data
  expr <- as_tibble(rownames_to_column(data.frame(assay(dat)), var = "ensembl_gene_id"))
  expr <- full_join(as_tibble(rowData(dat)), expr)
  row.names(expr) <- paste0(expr$external_gene_name, "_", expr$ensembl_gene_id)
  expr <- select(expr, grep("^^(TCGA)", names(expr)))
  expr <- t(expr)
  expr <- as_tibble(rownames_to_column(as.data.frame(expr), var = "barcode"))
  expr$barcode <- gsub("\\.", "-", expr$barcode)
  
  if(!is.null(GOIs)){
    genes <- paste0(GOIs, collapse = "|")
    genes <- paste0("^((", genes, "|)_|barcode)")
    expr <- select(expr, grep(genes, names(expr)))
  }
  
  expr <- full_join(expr, as_tibble(colData(dat)))
  
  return(expr)
}
