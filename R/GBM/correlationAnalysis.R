library(plyr)
library(reshape2)
library(Hmisc)
library(qvalue)
library(ggplot2)
library(xlsx)

# These are where the data files are stored
exFiles1 <- list.files("data/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3", full.names = T)
exFiles2 <- list.files("data/Expression-Genes/UNC__AgilentG4502A_07_2/Level_3", full.names = T)

exFiles <- c(exFiles1, exFiles2) # Merges lists into one
#exFiles <- exFiles2

expression <- llply(exFiles, function(x){
  df <- read.delim(x,stringsAsFactors = F,header=F)
  sample <- df[1,2]
  names(df) <- c("gene",sample)
  df <- df[3:nrow(df),]
  df[,2] <- as.numeric(df[,2])
  df
  
}, .progress = "text")

#takes awhile
expdf <- Reduce(function(x, y) merge(x, y, all=TRUE), expression)
remove("expression")

rownames(expdf) <- expdf$gene
exmat <- t(expdf[,2:ncol(expdf)])
remove("expdf")

#m <- cor(exmat) #also takes awhile
m2 <- rcorr(exmat, type = "pearson")
#corrplot(m, method = "circle")

lrp1b <- data.frame(gene = colnames(m2$r), cor =  m2$r[,"LRP1B"], pval = m2$P[,"LRP1B"])

qvals <- qvalue(p=lrp1b$pval)
lrp1b$qval.st <- qvals$qvalues
lrp1b$qval.bh <- p.adjust(p=lrp1b$pval, method = "BH")
lrp1b$lfdr <- qvals$lfdr
save(file = "lrp1b_corrs.rda", lrp1b)
write.csv(lrp1b, "lrp1b_corrs.csv", row.names = F)
write.table(lrp1b, "lrp1b_corrs.txt", row.names = F, quote = F, sep = "\t")

pdata <- data.frame(scale(exmat[,c("LRP1B","WNT3A","EGFR","MGMT","TERT","APC","TP53")]))
pdata$sample <- 1:nrow(pdata)

plot_fun <- function(y,dat){
  #eq <- lm(dat[,y] ~ LRP1B, dat)
  #r <- signif(sqrt(as.numeric(summary(eq)$r.squared)), 3)
  r <- round(cor(dat[,y], dat$LRP1B, method = "pearson"),2)
  rtxt <- paste0("r = ",r)
  
  p <- ggplot(dat, aes(x = LRP1B, y = dat[,y])) + 
    geom_point(color = "grey22", size = 0.5) +
    stat_smooth(method="lm", se = F) +
    #geom_text(aes(x = max(LRP1B)-0.5, y = max(dat[,y])), label = rtxt, size = 4) +
    xlab("LRP1B expression") +
    ylab(paste0(y," expression")) + 
    theme_bw() +
    theme(text = element_text(size=16),
          panel.border = element_rect(size = 1, color = "black"))
  if(r > 0){
    p <- p + geom_label(aes(x = min(LRP1B)+1, y = max(dat[,y])), 
                        label = rtxt, 
                        size = 2.5,
                        fill = "white")
  }else{
    p <- p + geom_label(aes(x = max(LRP1B)-0.9, 
                        y = max(dat[,y])), 
                        label = rtxt, 
                        size = 2.5,
                        fill = "white")  
  }
  p
}

lapply(c("WNT3A","APC","TERT","EGFR","MGMT","TP53"), function(x){
  p <- plot_fun(x,pdata)
  ggsave(paste0("results/",x,"_corr.png"), p, height = 1.5, width = 1.5)
  ggsave(paste0("results/",x,"_corr.pdf"), p, height = 2.5, width = 2.5, useDingbats = F)
  p
})

ggplot(lrp1b, aes(x = cor)) + 
  geom_histogram(binwidth = 0.025) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", size = 1),
        panel.grid = element_blank(),
        text = element_text(size = 20)) +
  xlim(c(-1,1)) +
  xlab("Pearson Correlation Coefficient, r") +
  ylab("Number of Genes") +
  geom_vline(xintercept = c(-0.5,0.5), size = 0.75, color = "blue", linetype = 2)

ggsave("results/cor_histogram.pdf", width = 7.5, height = 3, useDingbats = F)


network <- exmat
