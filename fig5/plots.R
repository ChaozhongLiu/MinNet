library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)

safe_sd <- function(x) {
  if (length(x) == 1)
    return(0.0)
  return(sd(x))
}

# JEM covid data UMAP ----
jem_meta <- read.csv('data/JEM_meta.csv',
                     stringsAsFactors = F)

umap_df <- read.csv('SiaNN_UMAP_PBMCmodel.txt',
                    stringsAsFactors = F)
all(jem_meta$X == umap_df$X)
obs_n1 <- sum(!startsWith(jem_meta$X, 'ATAC_'))
tech_label <- rep('ATAC',length(umap_df$index))
tech_label[1:obs_n1] <- 'GEX'
umap_df$Label <- tech_label
umap_df$dataset <- 'Modality'
umap_df$method <- 'PBMC model'

umap_df_tmp <- read.csv('SiaNN_UMAP_BMMCmodel.txt',
                        stringsAsFactors = F)
umap_df_tmp$Label <- tech_label
umap_df_tmp$dataset <- 'Modality'
umap_df_tmp$method <- 'BMMC model'
umap_df <- rbind(umap_df,umap_df_tmp)

labelss <- colnames(jem_meta)[2:4]
labels_name <- c('Sample', 'Cell type', 'Severity')
for (i in 1:3) {
  umap_df_tmp <- read.csv('SiaNN_UMAP_PBMCmodel.txt',
                          stringsAsFactors = F)
  umap_df_tmp$Label <- jem_meta[[labelss[i]]]
  umap_df_tmp$dataset <- labels_name[i]
  umap_df_tmp$method <- 'PBMC model'
  umap_df <- rbind(umap_df,umap_df_tmp)
  
  umap_df_tmp <- read.csv('SiaNN_UMAP_BMMCmodel.txt',
                          stringsAsFactors = F)
  umap_df_tmp$Label <- jem_meta[[labelss[i]]]
  umap_df_tmp$dataset <- labels_name[i]
  umap_df_tmp$method <- 'BMMC model'
  umap_df <- rbind(umap_df,umap_df_tmp)
}

umap_df[umap_df$Label == '0','Label'] <- 'Healthy'
umap_df[umap_df$Label == '4-5','Label'] <- 'Moderate'
umap_df[umap_df$Label == '6-7','Label'] <- 'Fatal'

lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:13] <- lvl_label[c(7:12,14,20,21,23:26)]
lvl_label[14] <- 'GEX'
lvl_label[15] <- 'ATAC'
lvl_label[16] <- 'Healthy'
lvl_label[17] <- 'Moderate'
lvl_label[18] <- 'Fatal'
lvl_label[19:26] <- unique(jem_meta$sample)

umap_df$Label <- factor(umap_df$Label,
                        levels = lvl_label)
umap_df$method <- factor(umap_df$method, levels = c('PBMC model',
                                                    'BMMC model'))

library(randomcoloR)
n <- length(lvl_label)
group.colors <- distinctColorPalette(n)
pie(rep(1,n), col=group.colors)
#group.colors <- c(group.colors, group.colors[1:8])
group.colors.saved <- readRDS('../fig2/color.palette.rds')

group.colors[14] <- group.colors.saved['GEX']
group.colors[15] <- group.colors.saved['ATAC']
group.colors[16] <- "#fab8ae"
group.colors[17] <- "#f25036"
group.colors[18] <- "#991e0a"

names(group.colors) <- lvl_label

gp <- ggplot(umap_df[,2:6], aes(x=X0, y=X1, color=Label)) +
  geom_point(size=0.025,alpha=0.5) + 
  theme_bw()+
  theme(text = element_text(size = 30))+
  scale_color_manual(values=group.colors)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 7.5,alpha=0.95)))+
  facet_grid(method ~ dataset, switch='y')+
  xlab('UMAP_1')+
  ylab('UMAP_2')

gp


# Label trasfer Heatmap ----
tansfer_acc_pbmc_rna <- read.csv('SiaNN_transfer_acc.PBMCmodel.rna.csv')
tansfer_acc_pbmc_rna$Direction <- 'ATAC to RNA'
tansfer_acc_pbmc_atac <- read.csv('SiaNN_transfer_acc.PBMCmodel.atac.csv')
tansfer_acc_pbmc_atac$Direction <- 'RNA to ATAC'
tansfer_acc <- rbind(tansfer_acc_pbmc_rna, tansfer_acc_pbmc_atac)
tansfer_acc$Model <- 'PBMC model'

tansfer_acc_bmmc_rna <- read.csv('SiaNN_transfer_acc.BMMCmodel.rna.csv')
tansfer_acc_bmmc_rna$Direction <- 'ATAC to RNA'
tansfer_acc_bmmc_rna$Model <- 'BMMC model'
tansfer_acc <- rbind(tansfer_acc,tansfer_acc_bmmc_rna)

tansfer_acc_bmmc_atac <- read.csv('SiaNN_transfer_acc.BMMCmodel.atac.csv')
tansfer_acc_bmmc_atac$Direction <- 'RNA to ATAC'
tansfer_acc_bmmc_atac$Model <- 'BMMC model'
tansfer_acc <- rbind(tansfer_acc,tansfer_acc_bmmc_atac)

gp <- ggplot(tansfer_acc, aes(cell_type, Model, fill = SiaNN))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Accuracy") +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  xlab("Cell type") +
  facet_wrap(~Direction) +
  coord_fixed()
# Print the heatmap
gp





