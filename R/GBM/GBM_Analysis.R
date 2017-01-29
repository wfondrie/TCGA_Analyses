source("R/functions/load_expr.R")
library(forcats)

expr <- load_expr("GBM", "LRP1B")

hist(log10(expr$days_to_death))

lmFit <- lm(days_to_death ~ LRP1B_ENSG00000168702, data = expr)

exprSumm <- expr %>% group_by(shortLetterCode) %>% 
  summarize(avg = mean(log10(LRP1B_ENSG00000168702)), st.dev = sd(log10(LRP1B_ENSG00000168702)))

expr$z_LRP1B <- (log10(expr$LRP1B_ENSG00000168702) - exprSumm$avg[exprSumm$shortLetterCode == "NT"]) /
  exprSumm$st.dev[exprSumm$shortLetterCode == "NT"]


lmFit <- lm(days_to_death ~ z_LRP1B, data = expr)
expr %>% ggplot(aes(y = days_to_death, x = z_LRP1B)) +
  geom_point()