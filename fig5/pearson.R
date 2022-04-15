library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(gridExtra)
options(stringsAsFactors = FALSE)


cor_0 <- read.csv('cor_by_severity/Mono_0_spr_cor.5NN.csv')
cor_4 <- read.csv('cor_by_severity/Mono_4-5_spr_cor.5NN.csv')
cor_6 <- read.csv('cor_by_severity/Mono_6-7_spr_cor.5NN.csv')

cor_0 <- cor_0[with(cor_0, order(genes, peaks)),]
cor_4 <- cor_4[with(cor_4, order(genes, peaks)),]
cor_6 <- cor_6[with(cor_6, order(genes, peaks)),]

all(cor_0$genes == cor_4$genes)
all(cor_0$genes == cor_6$genes)
all(cor_0$peaks == cor_4$peaks)
all(cor_0$peaks == cor_6$peaks)

cor_sum <- data.frame('genes' = cor_0$genes,
                      'peaks' = cor_0$peaks,
                      'cor_0' = cor_0$Spearman.cor,
                      'cor_4' = cor_4$Spearman.cor,
                      'cor_6' = cor_6$Spearman.cor)


multiome_data_0 <- read.csv('cor_by_severity/Mono_0_multiome_data.5NN.csv', row.names = 1)
multiome_data_0$current_severity_bin <- '0'
multiome_data_4 <- read.csv('cor_by_severity/Mono_4-5_multiome_data.5NN.csv', row.names = 1)
multiome_data_4$current_severity_bin <- '4-5'
multiome_data_6 <- read.csv('cor_by_severity/Mono_6-7_multiome_data.5NN.csv', row.names = 1)
multiome_data_6$current_severity_bin <- '6-7'




cor_sum$pr_0 <- 0
cor_sum$pr_4 <- 0
cor_sum$pr_6 <- 0

slot_k <- function(gene, peak, multiome_data, ignore_0 = F){
  multiome_data_missna <- multiome_data
  if (ignore_0) {
    multiome_data_missna[multiome_data_missna==0] <- NA
    multiome_data_missna[is.na(multiome_data_missna$current_severity_bin),
                         'current_severity_bin'] <- '0'
  }
  lfit = lm(multiome_data_missna[,peak] ~ multiome_data_missna[,gene])
  lfit.summ <- summary(lfit)
  return(lfit.summ$coefficients[2,1])
}

cor_sum_cp <- cor_sum
cor_sum_cp <- cor_sum_cp[with(cor_sum_cp, order(cor_6)),]

for (i in 1:4979) {
  gene <- cor_sum_cp$genes[i]
  gene <- gsub('-','.',gene)
  gene <- gsub(':','.',gene)
  peak <- cor_sum_cp$peaks[i]
  peak <- gsub('-','.',peak)
  peak <- gsub(':','.',peak)
  
  cor_sum_cp$pr_0[i] <- slot_k(gene, peak, multiome_data_0, ignore_0 = T)
  cor_sum_cp$pr_4[i] <- slot_k(gene, peak, multiome_data_4, ignore_0 = T)
  cor_sum_cp$pr_6[i] <- slot_k(gene, peak, multiome_data_6, ignore_0 = T)
}


write.table(cor_sum_cp, 'cor_by_severity/Mono_cor_sum.5NN.csv',quote=F, sep=',')


