# Load libraries ---------------------------------------------------------------
# from CRAN:
library(tidyverse)

# From Bioconductor:
library(TCGAbiolinks)
library(SummarizedExperiment)

# Download Data ----------------------------------------------------------------
# Only try to download data if T
download <- F
dataPath <- "E:/GDCdata" # Data is saved to external HDD

if(download) {
  projects <- getGDCprojects() # list projects
  
  # Function to download project files. Specifically RNA-Seq (HTSeq - FPKM-UQ)
  # Due to frequent errors when downloading, added a while loop and try statement.
  # Tries to download a maximum of 3 times.
  downloadExp <- function(study, gistic) {
    maxAttempts <- 3
    attempt <- 1
    studyNames <- gsub("\\-", "_", study)
    dat <- NULL
    
    while(is.null(dat) && attempt <= maxAttempts) {
      print(paste0("Trying attempt ", attempt, " for ", study, "..."))
      attempt <- attempt + 1
      
      try({
        query <- GDCquery(project = study,
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - FPKM-UQ")
        
        GDCdownload(query, directory = dataPath, method = "api")
        
        dat <- GDCprepare(query,
                          save = T,
                          save.filename = paste0("GDCdata/", studyNames, ".rda"),
                          directory = dataPath,
                          add.gistic2.mut = gistic)
      })
      
    }
    
    return(dat)
  }
  
  # Select only TCGA projects
  tcgaProjects <- projects %>% filter(grepl("^TCGA", project_id)) %>% select(project_id)
  
  # Download files
  expDat <- tcgaProjects %>% 
    group_by(project_id) %>% 
    do(summexp = downloadExp(.$project_id, NULL)) %>%
    mutate(failed = is.null(summexp)) %>%
    filter(!failed)

  
  if(nrow(expDat) < nrow(tcgaProjects)) {
    warning("1 or more projects failed to download. See `failed` variable.") 
    failed <- tcgaProjects %in% filter(!(project_id %in% expDat$project_id))
    }
  
  rm(expDat, expDat2)
  
} 



