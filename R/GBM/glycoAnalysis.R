library(tidyverse)
library(SummarizedExperiment)
library(forcats)
library(biomaRt)

# Load data --------------------------------------------------------------------
load("GDCdata/TCGA_GBM.rda")
gbm <- data
rm(data)

genes <- as.data.frame(rowRanges(gbm))
genes <- genes %>% dplyr::select(ensembl_gene_id, external_gene_name)
expr <- as.data.frame(assay(gbm))
expr$ensembl_gene_id <- row.names(expr)
clin <- as.data.frame(colData(gbm))
clin <- clin %>%
    dplyr::select(barcode, shortLetterCode, days_to_death, days_to_last_follow_up)


# Retrieve Ensembl gn's of glycosylation enzymes -------------------------------
# This makes it possible to rerun w/o internet if needed
if(!file.exists("temp/glycoEnzymes.rda")) {
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Retrieve Glycosidase proteins
gs <- getBM("ensembl_gene_id",
             filters = "go_id", 
             value = "GO:0016798", 
             mart = ensembl)

# Retrieve Glycotransferase proteins
gt <- getBM("ensembl_gene_id", 
             filters = "go_id", 
             value = "GO:0016757", 
             mart = ensembl)

save(gs, gt, file = "temp/glycoEnzymes.rda") 
rm(ensembl)

} else { load("temp/glycoEnzymes.rda") }

allEnz <- c(gs$ensembl_gene_id, gt$ensembl_gene_id)


# Reformat some things ---------------------------------------------------------
enzExpr <- expr %>% 
    filter(ensembl_gene_id %in% allEnz) %>%
    gather(barcode, FPKM, -ensembl_gene_id) %>%
    as_tibble() %>%
    left_join(clin) %>%
    left_join(genes) %>%
    mutate(enzType = ifelse(ensembl_gene_id %in% gs$ensembl_gene_id, "glycosidase", "glycotransferase"))

norm <- enzExpr %>% filter(enzExpr$shortLetterCode == "NT")
tumor <- enzExpr %>% filter(enzExpr$shortLetterCode != "NT")

t_test <- function(gene, fpkm, norm_df, ret) {
    norm_df <- norm_df %>% filter(ensembl_gene_id == gene)
    test <- t.test(log2(fpkm + 1), log2(norm_df$FPKM + 1))
    
    if(ret == "pval") {
         return(test$p.value)
    } else if(ret == "diff") {
         return(test$estimate[2] - test$estimate[1])
    }
}

diffExpr <- tumor %>% 
    group_by(ensembl_gene_id, external_gene_name, shortLetterCode, enzType) %>% 
    summarise(p_val = t_test(ensembl_gene_id[1], FPKM, norm, "pval"),
              diff = t_test(ensembl_gene_id[1], FPKM, norm, "diff")) %>%
    arrange(p_val, desc(diff))

# Some plots -------------------------------------------------------------------
gois <- diffExpr$external_gene_name[1:10]

enzExpr %>% 
    filter(external_gene_name %in% gois & shortLetterCode != "TR") %>%
    ggplot(aes(y = log2(FPKM + 1), x = shortLetterCode)) +
    geom_boxplot() +
    facet_grid( . ~ external_gene_name)


# Survival Curves --------------------------------------------------------------