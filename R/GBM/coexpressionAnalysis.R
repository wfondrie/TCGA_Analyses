library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(WGCNA)
library(fondrie)

key <- read.delim("data/key.txt", header = T)
key <- filter(key, Disease == "GBM")

seqFiles <- list.files("data/RNASeqV2/unc.edu_GBM.IlluminaHiSeq_RNASeqV2.Level_3.1.2.0",
                       pattern = ".rsem.genes.normalized_results",
                       full.names = T)

a <- read.delim(seqFiles[1])

expression <- ldply(seqFiles, function(x){
  df <- read.delim(x,stringsAsFactors = F,header=T)
  fName <- gsub("^.*/","",x)
  keyMatch <- grepl(fName, key$File.Name)
  df$barcode <- key$Barcode[keyMatch]
  df
}, .progress = "text")


expdf <- dcast(expression, barcode ~ gene_id, value.var = "normalized_count")
expmat <- as.matrix(expdf[ ,names(expdf) != "barcode"])
row.names(expmat) <- expdf$barcode

lowCountFilter <- colSums(expmat) < nrow(expmat)

expmat <- expmat[, !lowCountFilter]
expmat <- log2(expmat + 1)

lowVarFilter <- aaply(expmat, 2, var) > 0

expmat <- expmat[, lowVarFilter]

thold <- pickSoftThreshold(expmat, networkType = "signed", verbose = 10)

adjmat <- adjacency(expmat, typ= "signed", power = thold$powerEstimate)

geneTree <- hclust(as.dist(1-adjmat), method = "average")

modules <- cutreeDynamicTree(geneTree, minModuleSize = 200, deepSplit = F)


plotDendroAndColors(geneTree, labels2colors(modules), dendroLabels = F)

moddf <- data.frame(gene_name = row.names(adjmat), module = modules)

lrp1b <- moddf[grep("lrp1b", moddf$gene_name, ignore.case = T),  ]
lrp1b

graph <- adj2graphml(adjmat, nodeAttributes = moddf, threshold = 0.1)

