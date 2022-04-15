library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(gridExtra)
options(stringsAsFactors = FALSE)

# NK example plot ----
multiome_data <- read.csv('cor_by_severity/NK_multiome_data.5NN.csv', row.names = 1)
spr_cor <- read.csv('cor_by_severity/NK_spr_cor.5NN.csv', row.names = 1)
spr_cor$gene_ANOVA.p <- 0
spr_cor$peak_ANOVA.p <- 0
rownames(spr_cor) <- c(1:2358)

multiome_data_missna <- multiome_data
multiome_data_missna[multiome_data_missna==0] <- NA
multiome_data_missna[is.na(multiome_data_missna$current_severity_bin),
                     'current_severity_bin'] <- 0


for (i in 1:dim(spr_cor)[1]) {
  gene <- spr_cor$genes[i]
  gene <- gsub('-','.',gene)
  gene <- gsub(':','.',gene)
  peak <- spr_cor$peaks[i]
  peak <- gsub('-','.',peak)
  peak <- gsub(':','.',peak)
  
  res.aov <- aov(as.formula(paste0(gene, ' ~ ', 'current_severity_bin')), 
                 data = multiome_data_missna)
  res.aov <- summary(res.aov)
  spr_cor[i,'gene_ANOVA.p'] <- res.aov[[1]]$`Pr(>F)`[1]
  
  res.aov <- aov(as.formula(paste0(peak, ' ~ ', 'current_severity_bin')), 
                 data = multiome_data_missna)
  res.aov <- summary(res.aov)
  spr_cor[i,'peak_ANOVA.p'] <- res.aov[[1]]$`Pr(>F)`[1]
}

summary(spr_cor$gene_ANOVA.p)
summary(spr_cor$peak_ANOVA.p)

#spr_cor_2 <- spr_cor[(spr_cor$gene_ANOVA.p < 1e-3) & (spr_cor$peak_ANOVA.p < 0.05),]


i <- 3
gene <- spr_cor$genes[i]
gene <- gsub('-','.',gene)
gene <- gsub(':','.',gene)
peak <- spr_cor$peaks[i]
peak <- gsub('-','.',peak)
peak <- gsub(':','.',peak)


lfit = lm(multiome_data_missna[,peak] ~ multiome_data_missna[,gene])
lfit.summ <- summary(lfit)
severity.colors <- c('0'= "#fab8ae", '4-5'="#f25036", '6-7'="#991e0a")

p1 <- ggplot(multiome_data_missna[,c(gene,peak,'current_severity_bin')], 
            aes(x=multiome_data_missna[,gene], y=multiome_data_missna[,peak], 
            col=current_severity_bin)) +
  geom_point(size=2) +
  scale_colour_manual(values = severity.colors) +
  theme_bw() +
  theme(text = element_text(size = 12.5)) +
  geom_abline(intercept = lfit.summ$coefficients[1,1] , slope = lfit.summ$coefficients[2,1], color="gray30") +
  #ggtitle('  ') +
  guides(color=guide_legend(title="Severity")) +
  xlab(gsub('.','-',gene, fixed = TRUE)) +
  ylab(gsub('.','-',peak, fixed = TRUE))


p2 <- ggplot(multiome_data_missna[,c(gene,'current_severity_bin')],
             aes(fill=current_severity_bin, y=multiome_data_missna[,gene], x=current_severity_bin)) + 
  geom_violin(position="dodge", alpha=0.75, draw_quantiles = c(0.5)) + #, outlier.colour="transparent") +
  geom_jitter(height = 0, width = 0.1, size=1.5, alpha=0.6) +
  scale_fill_manual(values = severity.colors) +
  scale_color_manual(values = severity.colors) +
  #scale_fill_viridis(discrete=T, name="") +
  theme_bw()  +
  theme(text = element_text(size = 12.5)) +
  guides(fill=guide_legend(title="Severity")) +
  xlab("WHO Severity Group") +
  ylab(gsub('.','-',gene, fixed = TRUE))


p3 <- ggplot(multiome_data_missna[,c(peak,'current_severity_bin')],
       aes(fill=current_severity_bin, y=multiome_data_missna[,peak], x=current_severity_bin)) + 
  geom_violin(position="dodge", alpha=0.75, draw_quantiles = c(0.5)) + #, outlier.colour="transparent") +
  geom_jitter(height = 0, width = 0.1, size=1.5, alpha=0.6) +
  scale_fill_manual(values = severity.colors) +
  scale_color_manual(values = severity.colors) +
  #scale_fill_viridis(discrete=T, name="") +
  theme_bw()  +
  theme(text = element_text(size = 12.5)) +
  guides(fill=guide_legend(title="Severity")) +
  xlab("WHO Severity Group") +
  ylab(gsub('.','-',peak, fixed = TRUE))


# 14x3.85
grid.arrange(p1, p2, p3, nrow = 1)

# Mono 14 examples plot ----
multiome_data <- read.csv('cor_by_severity/Mono_multiome_data.5NN.csv', row.names = 1)
spr_cor <- read.csv('cor_by_severity/Mono_spr_cor.5NN.csv', row.names = 1)
spr_cor$gene_ANOVA.p <- 0
spr_cor$peak_ANOVA.p <- 0
rownames(spr_cor) <- c(1:4979)

multiome_data_missna <- multiome_data
multiome_data_missna[multiome_data_missna==0] <- NA
multiome_data_missna[is.na(multiome_data_missna$current_severity_bin),
                     'current_severity_bin'] <- 0

for (i in 1:dim(spr_cor)[1]) {
  gene <- spr_cor$genes[i]
  gene <- gsub('-','.',gene)
  gene <- gsub(':','.',gene)
  peak <- spr_cor$peaks[i]
  peak <- gsub('-','.',peak)
  peak <- gsub(':','.',peak)
  
  res.aov <- aov(as.formula(paste0(gene, ' ~ ', 'current_severity_bin')), 
                 data = multiome_data_missna)
  res.aov <- summary(res.aov)
  spr_cor[i,'gene_ANOVA.p'] <- res.aov[[1]]$`Pr(>F)`[1]
  
  res.aov <- aov(as.formula(paste0(peak, ' ~ ', 'current_severity_bin')), 
                 data = multiome_data_missna)
  res.aov <- summary(res.aov)
  spr_cor[i,'peak_ANOVA.p'] <- res.aov[[1]]$`Pr(>F)`[1]
}

summary(spr_cor$gene_ANOVA.p)
summary(spr_cor$peak_ANOVA.p)

#spr_cor_2 <- spr_cor[(spr_cor$gene_ANOVA.p < 1e-3) & (spr_cor$peak_ANOVA.p < 0.05),]

multiome_data_missna <- multiome_data
multiome_data_missna[multiome_data_missna==0] <- NA
multiome_data_missna[is.na(multiome_data_missna$current_severity_bin),
                     'current_severity_bin'] <- 0

i <- 2 #1905 #185 ###278 ##182 ##146 ##4851 ##4976 ##74 ###46 #38 #33 #27 #26 #25
gene <- spr_cor$genes[i]
gene <- gsub('-','.',gene)
gene <- gsub(':','.',gene)
peak <- spr_cor$peaks[i]
peak <- gsub('-','.',peak)
peak <- gsub(':','.',peak)

lfit = lm(multiome_data_missna[,peak] ~ multiome_data_missna[,gene])
lfit.summ <- summary(lfit)
severity.colors <- c('0'= "#fab8ae", '4-5'="#f25036", '6-7'="#991e0a")

p1 <- ggplot(multiome_data_missna[,c(gene,peak,'current_severity_bin')], 
             aes(x=multiome_data_missna[,gene], y=multiome_data_missna[,peak], 
                 col=current_severity_bin)) +
  geom_point(size=2) +
  scale_colour_manual(values = severity.colors) +
  theme_bw() +
  theme(text = element_text(size = 12.5)) +
  geom_abline(intercept = lfit.summ$coefficients[1,1] , slope = lfit.summ$coefficients[2,1], color="gray30") +
  #ggtitle('  ') +
  guides(color=guide_legend(title="Severity")) +
  xlab(gsub('.','-',gene, fixed = TRUE)) +
  ylab(gsub('.','-',peak, fixed = TRUE))


p2 <- ggplot(multiome_data_missna[,c(gene,'current_severity_bin')],
             aes(fill=current_severity_bin, y=multiome_data_missna[,gene], x=current_severity_bin)) + 
  geom_violin(position="dodge", alpha=0.75, draw_quantiles = c(0.5)) + #, outlier.colour="transparent") +
  geom_jitter(height = 0, width = 0.1, size=1.5, alpha=0.6) +
  scale_fill_manual(values = severity.colors) +
  scale_color_manual(values = severity.colors) +
  #scale_fill_viridis(discrete=T, name="") +
  theme_bw()  +
  theme(text = element_text(size = 12.5)) +
  guides(fill=guide_legend(title="Severity")) +
  xlab("WHO Severity Group") +
  ylab(gsub('.','-',gene, fixed = TRUE))


p3 <- ggplot(multiome_data_missna[,c(peak,'current_severity_bin')],
             aes(fill=current_severity_bin, y=multiome_data_missna[,peak], x=current_severity_bin)) + 
  geom_violin(position="dodge", alpha=0.75, draw_quantiles = c(0.5)) + #, outlier.colour="transparent") +
  geom_jitter(height = 0, width = 0.1, size=1.5, alpha=0.6) +
  scale_fill_manual(values = severity.colors) +
  scale_color_manual(values = severity.colors) +
  #scale_fill_viridis(discrete=T, name="") +
  theme_bw()  +
  theme(text = element_text(size = 12.5)) +
  guides(fill=guide_legend(title="Severity")) +
  xlab("WHO Severity Group") +
  ylab(gsub('.','-',peak, fixed = TRUE))

#p1
#p2
#p3
# 14x3.85
grid.arrange(p1, p2, p3, nrow = 1)


# Monocytes correlation by severity ----
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

plot_example_dots <- function(gene, peak, multiome_data, ignore_0 = F){
  multiome_data_missna <- multiome_data
  if (ignore_0) {
    multiome_data_missna[multiome_data_missna==0] <- NA
    multiome_data_missna[is.na(multiome_data_missna$current_severity_bin),
                         'current_severity_bin'] <- '0'
  }
  lfit = lm(multiome_data_missna[,peak] ~ multiome_data_missna[,gene])
  lfit.summ <- summary(lfit)
  ggplot(multiome_data_missna[,c(gene,peak,'current_severity_bin')], 
               aes(x=multiome_data_missna[,gene], y=multiome_data_missna[,peak],
                   col=current_severity_bin)) +
    geom_point(size=2) +
    scale_colour_manual(values = severity.colors) +
    theme_bw() +
    theme(text = element_text(size = 12.5)) +
    geom_abline(intercept = lfit.summ$coefficients[1,1] , slope = lfit.summ$coefficients[2,1], color="gray30") +
    #ggtitle('  ') +
    guides(color=guide_legend(title="Severity")) +
    xlab(gsub('.','-',gene, fixed = TRUE)) +
    ylab(gsub('.','-',peak, fixed = TRUE))
}

pair_screen <- function(i, ignore_0=F){
  gene <- cor_sum$genes[i]
  gene <- gsub('-','.',gene)
  gene <- gsub(':','.',gene)
  peak <- cor_sum$peaks[i]
  peak <- gsub('-','.',peak)
  peak <- gsub(':','.',peak)
  
  
  p1 <- plot_example_dots(gene, peak, multiome_data_0, ignore_0 = ignore_0)
  #p1
  p2 <- plot_example_dots(gene, peak, multiome_data_4, ignore_0 = ignore_0)
  #p2
  p3 <- plot_example_dots(gene, peak, multiome_data_6, ignore_0 = ignore_0)
  #p3
  return(grid.arrange(p1, p2, p3, nrow = 1))
}

severity.colors <- c('0'= "#fab8ae", '4-5'="#f25036", '6-7'="#991e0a")


cor_interest <- cor_sum_cp[(cor_sum_cp$pr_0>0)&(cor_sum_cp$pr_4>0)&(cor_sum_cp$pr_6<0),]

p <- pair_screen(i=1398, ignore_0 = T)



