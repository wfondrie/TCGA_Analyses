library(plyr)
library(ggplot2)


pxFile <- list.files("data/Clinical/Biotab", pattern="drug", full.names=T)
px <- read.delim(pxFile, skip=1, stringsAsFactors = F)
pxAll <- px
px <- px[2:nrow(px),]
px <- px[grepl("^(tem|bev|ava)", px$drug_name, ignore.case = T),]
px$drug_name[grepl("bev|ava", px$drug_name, ignore.case = T)] <- "bevacizumab"
px$drug_name[grepl("^tem", px$drug_name, ignore.case=T)] <- "temozolamide"

patients <- unique(px$bcr_patient_barcode)

keyFile1 <- list.files("data/METADATA/UNC__AgilentG4502A_07_1/", pattern="sdrf", full.names=T)
keyFile2 <- list.files("data/METADATA/UNC__AgilentG4502A_07_2/", pattern="sdrf", full.names=T)
key1 <- read.delim(keyFile1, stringsAsFactors = F)
key2 <- read.delim(keyFile2, stringsAsFactors = F)
key <- merge(key1,key2, all=T)
key <- key[,c("Source.Name","Array.Data.File")]


exFiles1 <- list.files("data/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3", full.names = T)
exFiles2 <- list.files("data/Expression-Genes/UNC__AgilentG4502A_07_2/Level_3", full.names = T)

exFiles <- c(exFiles1, exFiles2)



lrp1b <- ldply(patients, .progress = "text",  function(x){
  keyvalue <- key[grepl(x, key[,1]),]
  keyvalue <- keyvalue$Array.Data.File
  file <- grep(keyvalue, exFiles, value=T)
  if(length(file) > 0){
    exp <- read.delim(file[1], stringsAsFactors = F, header = F)
    vals <- exp[grepl("lrp1b", exp[,1], ignore.case=T),]
  } else {
    vals <- data.frame("V1"=NA, "V2"=NA, "REF"=NA)
  }
  vals$patient <- x
  vals
})

pxKeep <- c("bcr_patient_barcode", "drug_name")

pxSum <- ddply(px, "bcr_patient_barcode", function(x) {
  rows <- x
  if(nrow(x) == 1) {
    fin <- rows[, pxKeep]
  } else {
    rows <- rows[, pxKeep]
    fin <- rows[1,]
    
    if(all(rows$drug_name == "bevacizumab")) { fin$drug_name == "bevacizumab"}
    else if(all(rows$drug_name == "temozolamide")) {fin$drug_name == "temozolamide"}
    else { fin$drug_name == "both" }
  }
  fin
})


names(pxSum)[1] <- "patient"
pData <- merge(pxSum, lrp1b, all = T)
names(pData)[3:4] <- c("gene","log2_ratio")
pData$log2_ratio <- as.numeric(pData$log2_ratio)

##### Patient Clinical Data ##############################################################
pFile <- list.files("data/Clinical/Biotab", patter="patient_gbm.txt", full.names=T)
p <- read.delim(pFile, skip=1, stringsAsFactors = F)
p <- p[p$bcr_patient_barcode %in% pData$patient,]
p$patient <- p$bcr_patient_barcode

keep <- c("patient", "prior_glioma", "gender", "days_to_birth", "race", "vital_status", "days_to_death",
          "karnofsky_performance_score", "age_at_initial_pathologic_diagnosis","histological_type")

p <- p[,keep]
##########################################################################################


bev <- pData$log2_ratio[pData$drug_name == "bevacizumab"]
tem <- pData$log2_ratio[pData$drug_name == "temozolamide"]

test <- t.test(bev,tem)

pBar <- data.frame(x = c(1,1,2,2), y = c(7.3,7.5,7.5,7.3))
p <- ggplot(pData, aes(y=log2_ratio, x = drug_name)) + 
  geom_boxplot() +
  xlab("Drug") +
  ylab("Ratio to standard (log2)")
p <- p + geom_path(data = pBar, aes(x=x, y=y), size = 1) +
  annotate("text", x = 1.5, y = 8, label= paste0("p = ",round(test$p.value, digits = 3))) +
  annotate("text", x = 1, y = -1, label = paste0("n = ", length(bev))) +
  annotate("text", x = 2, y = -1, label = paste0("n = ", length(tem))) 
p

ggsave("temp/drug_test.png", width = 4, height = 4)

