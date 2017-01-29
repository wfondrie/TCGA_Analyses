library(plyr)
library(dplyr)
library(survival)

# Read in clinical data for all patients
clinFile <- list.files("data/Clinical/Biotab", pattern="patient", full.names=T)
clin <- read.delim(clinFile, skip=1, stringsAsFactors = F)
clin <- clin[-1,] # secondary headers are in the first row, so remove them.

# Define what clinical data to keep
keep <- c("bcr_patient_uuid",
          "bcr_patient_barcode",
          "days_to_last_followup",
          "days_to_death")
clin <- clin[,keep]

# Format data for Kaplan Meier survival analysis
clin$time <- as.numeric(clin$days_to_death) # will produce NAs
clin$death <- rep(T, nrow(clin))
clin$death[is.na(clin$time)] <- F
clin$time[is.na(clin$time)] <- clin$days_to_last_followup[is.na(clin$time)]
clin$time <- as.numeric(clin$time)

key <- read.delim("data/key.txt", header = T)
key <- filter(key, Disease == "GBM")

# Read RNA-Seq data
seqFiles <- list.files("data/RNASeqV2/unc.edu_GBM.IlluminaHiSeq_RNASeqV2.Level_3.1.2.0",
                       pattern = ".rsem.genes.normalized_results",
                       full.names = T)

# Reads in the RNA-Seq RSEM for LRP1B and matches the filename to UUID and barcodes
rna <- ldply(seqFiles, function(x) {
  dat <- read.delim(x)
  fName <- gsub("^.*/","",x)
  dat <- dat[grepl("LRP1B", dat[,1]),2]
  keyMatch <- grepl(fName, key$File.Name)
  UUID <- key$UUID[keyMatch]
  barcode <- key$Barcode[keyMatch]
  data.frame(file = fName, 
             bcr_patient_barcode_seq = barcode, 
             LRP1B = dat)
})

rna$bcr_patient_barcode <- substr(rna$bcr_patient_barcode_seq, 1, 12)
rna$log.LRP1B <- log(rna$LRP1B)

# Determining normal and tumor samples by the character 14 of barcodes
# 1 = normal, 0 = tumor
nIdx <- which(substr(rna$bcr_patient_barcode_seq, 14,14) == "1")
tIdx <- which(substr(rna$bcr_patient_barcode_seq, 14,14) == "0")

# Scale RNA RSEM to give Z-score, then calculate Z-test p-value
rnaScale <- function(tumor, normal) {
  normMean <- mean(normal)
  normSD <- sd(normal)
  (tumor - normMean)/normSD
}

rnaTumor <- rna[tIdx,]
rnaTumor$z.score <- rnaScale(rnaTumor$log.LRP1B,rna$log.LRP1B[nIdx])
rnaTumor$p.value <- pnorm(rnaTumor$z.score)

dat <- merge(rnaTumor, clin)
dat$event <- ifelse(dat$p.value <= 0.05, 1, 0)

sdat <- Surv(dat$time, dat$death)
fit <- survfit(sdat~dat$event)
plot(fit, col = c(1:3))

ggsurv(fit)
ggsave("results/kp-05.png")
