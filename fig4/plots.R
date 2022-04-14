library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(gridExtra)
options(stringsAsFactors = FALSE)

safe_sd <- function(x) {
  if (length(x) == 1)
    return(0.0)
  return(sd(x))
}

library(randomcoloR)
n <- 6
# group.colors <- distinctColorPalette(n)
# pie(rep(1,n), col=group.colors)
# names(group.colors) <- c('Non', 'True Pair','5NN',
#                                   '10NN', '15NN', '20NN')
# 
# saveRDS(group.colors, 'color.palette.rds')
group.colors <- readRDS('color.palette.rds')


# Smoothing sparsity comparison ----
spar_df <- read.csv('sparsity_sum.csv', stringsAsFactors = F)
spar_df$Smoothing[spar_df$Smoothing == 'No'] <- 'No Smoothing'
spar_df$Smoothing <- factor(spar_df$Smoothing, levels = c('No Smoothing', '5NN', '10NN',
                                                          '15NN', '20NN'))
color_palette <- group.colors[names(group.colors) %in% c('No Smoothing', '5NN', '10NN',
                                                         '15NN', '20NN')]

ggplot(spar_df, aes(x=Feature, y=Sparsity, fill=Smoothing)) + 
  #geom_violin(alpha=0.5) +
  geom_boxplot(position=position_dodge(0.9),
               notch=TRUE,
               width=0.5,
               alpha=0.75,
               outlier.colour="lightgrey",
               outlier.fill="lightgrey",
               outlier.size=0.1) +
  scale_fill_manual(name = "Smoothing", 
                    values = color_palette) +
  xlab('Modality')+
  #facet_wrap(~Dataset) +
  theme_bw()

# Smoothing correlation comparison ----
library(ggpubr)
noSm <- read.csv('correlation/spr.cor.pbmc.2k.NoSm.csv')
summary(noSm$is_pcHiC)
pair <- read.csv('correlation/spr.cor.pbmc.2k.Pair.csv')
cor_2kb <- rbind(noSm, pair)

for (n in c(5,10,15,20)) {
  Sm <- read.csv(paste0('correlation/spr.cor.pbmc.2k.',n,'NN.csv'))
  cor_2kb <- rbind(cor_2kb, Sm)
}

summary(cor_2kb$Method)
cor_2kb$Method <- as.character(cor_2kb$Method)
cor_2kb$Method[cor_2kb$Method == 'SiaNN'] <- 'Non'
cor_2kb$Method[cor_2kb$Method == 'True_pair'] <- 'True Pair'

cor_2kb$Method <- factor(cor_2kb$Method, levels=c('Non', 'True Pair', '5NN', '10NN',
                                                  '15NN', '20NN'))
df_mean <- cor_2kb %>% 
  group_by(Method) %>% 
  summarize(average = mean(Spearman.cor)) %>%
  ungroup()

color_palette <- group.colors

ggplot(cor_2kb, aes(x=Method, y=Spearman.cor, fill=Method)) + 
  #geom_violin(alpha=0.5) +
  geom_boxplot(position=position_dodge(0.5),
               notch=TRUE,
               width=0.7,
               alpha=0.75,
               outlier.colour="lightgrey",
               outlier.fill="lightgrey",
               outlier.size=0.1) +
  geom_point(data = df_mean,
             mapping = aes(x = Method, y = average),
             color="red", alpha=0.75)+
  geom_line(data = df_mean, 
            mapping = aes(x = Method, y = average, group=1))+
  scale_fill_manual(name = "Method", 
                    values = color_palette) +
  #ylim(0, 0.7) +
  xlab(' ')+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "True Pair", paired = T, label.y.npc = 'bottom') +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Non", paired = T) + 
  geom_text(data=df_mean, aes(x=Method, y=average + 0.08, label=round(df_mean$average,3)), col='black', size=3)+
  
  #facet_wrap(~Dataset) +
  theme_bw()

tt <- t.test(noSm$Spearman.cor, Sm$Spearman.cor, paired = T)$p.value

sum(noSm$Spearman.cor < 0)
sum(Sm$Spearman.cor < 0)

mean(abs(noSm$Spearman.cor) - abs(Sm$Spearman.cor)) / mean(abs(noSm$Spearman.cor))


# Example correlation ----
noSm <- read.csv('correlation/spr.cor.pbmc.2k.NoSm.csv')
summary(noSm$is_pcHiC)

sm_5nn <- read.csv('correlation/spr.cor.pbmc.2k.10NN.csv')
noSm <- noSm[with(sm_5nn, order(Spearman.cor, decreasing = T)),]
sm_5nn <- sm_5nn[with(sm_5nn, order(Spearman.cor, decreasing = T)),]

examples <- sm_5nn[(sm_5nn$Spearman.cor * noSm$Spearman.cor) < 0, ]
examples$noSm <- noSm[(sm_5nn$Spearman.cor * noSm$Spearman.cor) < 0, 'Spearman.cor']

data_0 <- read.csv('correlation/pseudo_data.noSm.csv', row.names = 1)
data_5 <- read.csv('correlation/pseudo_data.5NN.csv', row.names = 1)
data_10 <- read.csv('correlation/pseudo_data.10NN.csv', row.names = 1)
data_15 <- read.csv('correlation/pseudo_data.15NN.csv', row.names = 1)
data_20 <- read.csv('correlation/pseudo_data.20NN.csv', row.names = 1)

plot_example_dots <- function(gene, peak, multiome_data, ignore_0 = F, md){
  multiome_data_missna <- multiome_data
  if (ignore_0) {
    multiome_data_missna[multiome_data_missna==0] <- NA
  }
  lfit = lm(multiome_data_missna[,peak] ~ multiome_data_missna[,gene])
  lfit.summ <- summary(lfit)
  multiome_data_missna$md <- md
  ggplot(multiome_data_missna[,c(gene,peak,'md')], 
         aes(x=multiome_data_missna[,gene], y=multiome_data_missna[,peak],
             col = md)) +
    geom_point(size=2, alpha=0.75) +
    scale_colour_manual(values = group.colors) +
    theme_bw() +
    theme(text = element_text(size = 10), legend.position="none") +
    geom_abline(intercept = lfit.summ$coefficients[1,1] , slope = lfit.summ$coefficients[2,1], color="gray30") +
    #ggtitle('  ') +
    guides(color=guide_legend(title="Method")) +
    xlab(gsub('.','-',gene, fixed = TRUE)) +
    ylab(gsub('.','-',peak, fixed = TRUE))
}

i <- 1
gene <- examples$genes[i]
gene <- gsub('-','.',gene)
gene <- gsub(':','.',gene)
peak <- examples$peaks[i]
peak <- gsub('-','.',peak)
peak <- gsub(':','.',peak)



p1 <- plot_example_dots(gene, peak, data_0, ignore_0 = F, md='No Smoothing')
p1

p2 <- plot_example_dots(gene, peak, data_5, ignore_0 = F, md='5NN')
p2

p3 <- plot_example_dots(gene, peak, data_10, ignore_0 = F, md='10NN')
p3

p4 <- plot_example_dots(gene, peak, data_15, ignore_0 = F, md='15NN')
p4

p5 <- plot_example_dots(gene, peak, data_20, ignore_0 = F, md='20NN')
p5

grid.arrange(p1, p2, p3, p4, p5, nrow = 1)


# Correlation increase and pair-wise t-test ----
noSm <- read.csv('correlation/spr.cor.pbmc.2k.NoSm.csv')
noSm <- noSm[noSm$is_pcHiC == 'True',]
pair <- read.csv('correlation/spr.cor.pbmc.2k.Pair.csv')
pair <- pair[pair$is_pcHiC == 'True',]

Sm <- read.csv(paste0('correlation/spr.cor.pbmc.2k.',5,'NN.csv'))
Sm <- Sm[Sm$is_pcHiC == 'True',]

increase_noSm <- Sm
increase_noSm$delta.Spr <- Sm$Spearman.cor - noSm$Spearman.cor

increase_tp <- Sm
increase_tp$delta.Spr <- Sm$Spearman.cor - pair$Spearman.cor

method <- c('5NN','10NN','15NN','20NN')
# test.p.list.no <- c(0,0,0,0)
# names(test.p.list.no) <- method
# test.p.list.no[1] <- t.test(noSm$Spearman.cor, Sm$Spearman.cor, paired = T)$p.value
size.list.no <- c(0,0,0,0)
names(size.list.no) <- method
size.list.no[1] <- mean(Sm$Spearman.cor - noSm$Spearman.cor) / mean(noSm$Spearman.cor)
# 
# test.p.list.tp <- c(0,0,0,0)
# names(test.p.list.tp) <- method
# test.p.list.tp[1] <- t.test(pair$Spearman.cor, Sm$Spearman.cor, paired = T)$p.value
size.list.tp <- c(0,0,0,0)
names(size.list.tp) <- method
size.list.tp[1] <- mean(Sm$Spearman.cor - pair$Spearman.cor) / mean(pair$Spearman.cor)


for (n in 2:4) {
  Sm <- read.csv(paste0('correlation/spr.cor.pbmc.2k.', method[n],'.csv'))
  Sm <- Sm[Sm$is_pcHiC == 'True',]
  Sm$delta.Spr <- Sm$Spearman.cor - noSm$Spearman.cor
  increase_noSm <- rbind(increase_noSm, Sm)
  #test.p.list.no[n] <- t.test(noSm$Spearman.cor, Sm$Spearman.cor, paired = T)$p.value
  size.list.no[n] <- mean(Sm$Spearman.cor - noSm$Spearman.cor) / mean(noSm$Spearman.cor)
  
  
  Sm <- read.csv(paste0('correlation/spr.cor.pbmc.2k.',method[n],'.csv'))
  Sm <- Sm[Sm$is_pcHiC == 'True',]
  Sm$delta.Spr <- Sm$Spearman.cor - pair$Spearman.cor
  increase_tp <- rbind(increase_tp, Sm)
  #test.p.list.tp[n] <- t.test(pair$Spearman.cor, Sm$Spearman.cor, paired = T)$p.value
  size.list.tp[n] <- mean(Sm$Spearman.cor - pair$Spearman.cor) / mean(pair$Spearman.cor)
}
increase_tp$dataset <- 'Comparing with true pair'
increase_noSm$dataset <- 'Comparing with non-smoothing'

increase_df <- rbind(increase_tp, increase_noSm)
increase_df$Method <- factor(increase_df$Method, levels = method)
color_palette <- group.colors[names(group.colors) %in% method]

text_df <- data.frame(size = c(size.list.tp, size.list.no),
                      dataset = rep(c('Comparing with true pair','Comparing with non-smoothing'), each=4),
                      Method = rep(method,2))
text_df$size <- round(text_df$size*100, 2)
text_df$size <- paste0(text_df$size,'%')

ggplot(increase_df, aes(x=Method, y=delta.Spr, fill=Method)) + 
  #geom_violin(alpha=0.5) +
  geom_boxplot(position=position_dodge(0.9),
               notch=TRUE,
               width=0.7,
               alpha=0.75,
               outlier.colour="lightgrey",
               outlier.fill="lightgrey",
               outlier.size=0.1) +
  scale_fill_manual(name = "Method", 
                    values = group.colors) +
  #ylim(0, 0.7) +
  xlab(' ')+
  facet_wrap(~dataset) +
  geom_hline(yintercept = 0.0, linetype='dashed')+
  geom_text(data=text_df, aes(x=Method, y=rep(0.5,8), label=size), col='black', size=3)+
  theme_bw()



# pcHi-C evidence validation ----
distance_order <- c('0-25kb', '25-50kb', '50-75kb',
                    '75-100kb', '100-125kb', '125-150kb')
distance_bin <- function(cor_dat){
  cor_dat$tss_dist <- -cor_dat$tss_dist
  cor_dat$distance <- '0'
  cor_dat$distance[cor_dat$tss_dist<=150000] <- '125-150kb'
  cor_dat$distance[cor_dat$tss_dist<=125000] <- '100-125kb'
  cor_dat$distance[cor_dat$tss_dist<=100000] <- '75-100kb'
  cor_dat$distance[cor_dat$tss_dist<=75000] <- '50-75kb'
  cor_dat$distance[cor_dat$tss_dist<=50000] <- '25-50kb'
  cor_dat$distance[cor_dat$tss_dist<=25000] <- '0-25kb'
  cor_dat$distance <- factor(cor_dat$distance,
                             levels = distance_order)
  print(summary(cor_dat$distance))
  return(cor_dat)
}

cor_pc_0 <- read.csv('correlation/spr.cor.pbmc.150k.NoSm.csv')
cor_pc_0 <- distance_bin(cor_pc_0)

cor_pc_tp <- read.csv('correlation/spr.cor.pbmc.150k.Pair.csv')
cor_pc_tp <- distance_bin(cor_pc_tp)

cor_pc <- rbind(cor_pc_0,cor_pc_tp)

for (i in 3:6) {
  cor_pc_sm <- read.csv(paste0('correlation/spr.cor.pbmc.150k.', methods[i],'.csv'))
  cor_pc_sm <- distance_bin(cor_pc_sm)
  cor_pc <- rbind(cor_pc,cor_pc_sm)
}

cor_pc[is.na(cor_pc)] <- 0
cor_pc$Spearman.cor <- abs(cor_pc$Spearman.cor)

cor_pc$Method[cor_pc$Method == 'SiaNN'] <- 'Non'
cor_pc$Method[cor_pc$Method == 'True_pair'] <- 'True Pair'
cor_pc$Method <- factor(cor_pc$Method, levels = methods)

cor_pc[cor_pc$is_pcHiC=='False', 'is_pcHiC'] <- 'False'
cor_pc[cor_pc$is_pcHiC=='True', 'is_pcHiC'] <- 'True'


cor_pc_sum <- cor_pc %>% 
  group_by(distance, is_pcHiC, Method) %>% 
  summarize(average = mean(Spearman.cor)) %>%
  ungroup()

ggplot(cor_pc, aes(x=is_pcHiC, y=Spearman.cor, fill=Method)) + 
  #geom_violin(alpha=0.5) +
  geom_boxplot(position=position_dodge(0.9),
               notch=TRUE,
               width=0.7,
               alpha=0.75,
               outlier.colour="grey",
               outlier.fill="grey",
               outlier.size=0.75) +
  #guides(fill=guide_legend(title="pcHi-C")) +
  scale_fill_manual(name = "Method", 
                    values = group.colors) +
  xlab('Method')+
  #ylim(0,1)+
  facet_wrap(~distance) +
  theme_bw()



ccr2 <- cor_pc_5[cor_pc_5$genes == 'CCR2',]

# pcHiC evidence mean heatmap ----
methods <- names(group.colors)
cor_pc_0 <- read.csv('correlation/spr.cor.pbmc.150k.NoSm.csv')
cor_pc_0 <- distance_bin(cor_pc_0)

cor_pc_tp <- read.csv('correlation/spr.cor.pbmc.150k.Pair.csv')
cor_pc_tp <- distance_bin(cor_pc_tp)

cor_pc <- rbind(cor_pc_0,cor_pc_tp)

for (i in 3:6) {
  cor_pc_sm <- read.csv(paste0('correlation/spr.cor.pbmc.150k.', methods[i],'.csv'))
  cor_pc_sm <- distance_bin(cor_pc_sm)
  cor_pc <- rbind(cor_pc,cor_pc_sm)
}

cor_pc[is.na(cor_pc)] <- 0
cor_pc$Spearman.cor <- abs(cor_pc$Spearman.cor)

summary(factor(cor_pc$Method))
cor_pc$Method[cor_pc$Method == 'SiaNN'] <- 'Non'
cor_pc$Method[cor_pc$Method == 'True_pair'] <- 'True Pair'
cor_pc$Method <- factor(cor_pc$Method, levels = methods)

cor_pc[cor_pc$is_pcHiC=='False', 'is_pcHiC'] <- 'pcHi-C: False'
cor_pc[cor_pc$is_pcHiC=='True', 'is_pcHiC'] <- 'pcHi-C: True'

cor_pc_sum <- cor_pc %>% 
  group_by(distance,is_pcHiC,Method) %>% 
  summarize(average = mean(Spearman.cor)) %>%
  ungroup()

gp <- ggplot(cor_pc_sum, aes(Method, distance, fill = average))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.35, limit = c(0.19,0.5),
                       space = "Lab", 
                       name="Spearman.cor") +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  ylab("Distance") +
  facet_wrap(~is_pcHiC) +
  coord_fixed()
# Print the heatmap
gp

