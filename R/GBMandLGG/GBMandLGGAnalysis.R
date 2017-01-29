library(tidyverse)
library(SummarizedExperiment)
library(forcats)

load("GDCdata/TCGA_GBM.rda")
gbm <- data
load("GDCdata/TCGA_LGG.rda")
lgg <- data
rm(data)

dat <- cbind(gbm, lgg)
rm(gbm,lgg)

goi <- as_tibble(rowData(dat)[grep("^(LRP1|LRP1B)$", rowData(dat)$external_gene_name), ])
expr <- as.data.frame(t(assay(dat)[row.names(dat) %in% goi$ensembl_gene_id, ]))
names(expr) <- c("LRP1", "LRP1B")
expr$barcode <- row.names(expr)
expr <- as_tibble(expr)
clin <- as_tibble(colData(dat))

btumors <- full_join(expr,clin)
btumors$project_id <- as.factor(btumors$project_id)

naGrade <- btumors %>% filter(is.na(subtype_Grade))

btumors <- btumors %>% 
  select(LRP1, LRP1B, barcode, shortLetterCode, subtype_Grade, 
         days_to_death, days_to_last_follow_up, project_id)

btumors$subtype_Grade <- as.character(btumors$subtype_Grade)
btumors$subtype_Grade[btumors$shortLetterCode == "NT"] <- "NT"

btumors$subtype_Grade <- fct_relevel(btumors$subtype_Grade, c("NT","G2","G3","G4"))

btumors %>% 
  filter(!is.na(subtype_Grade)) %>%
  ggplot(aes(x = log10(LRP1), y = log10(LRP1B), color = subtype_Grade)) +
  scale_color_discrete("Tumor Grade") +
  xlab("LRP1 Expression (Log10(FPKM))") +
  ylab("LRP1B Expression (Log10(FPKM))") +
  geom_point()
ggsave("scatter.png")

btumors %>% 
  gather(Gene, FPKM, LRP1, LRP1B) %>%
  filter(!is.na(subtype_Grade)) %>%
  ggplot(aes(x = Gene, y = log10(FPKM), fill = subtype_Grade)) +
  scale_fill_discrete("Tumor Grade") +
  geom_boxplot()
  
ggsave("box.png")     

clin %>% ggplot(aes(x = days_to_death)) + geom_histogram()
sum(!is.na(clin$days_to_death))
