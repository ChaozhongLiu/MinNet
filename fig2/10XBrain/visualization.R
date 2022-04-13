library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)


meta <- read.csv('../data/rna.meta.csv',
                 stringsAsFactors = F, sep='\t')

umap_df <- read.csv('SiaNN_UMAP.txt',
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$seurat_clusters
umap_df$dataset <- 'Clusters'

umap_df_tmp <- read.csv('SiaNN_UMAP.txt',
                        stringsAsFactors = F)
umap_df_tmp$Label <- c(rep('GEX',2423),rep('ATAC',2423))
umap_df_tmp$dataset <- 'Modality'
umap_df <- rbind(umap_df,umap_df_tmp)

umap_df$dataset <- factor(umap_df$dataset, levels = c('Clusters',
                                                      'Modality'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[8] <- 'GEX'
lvl_label[9] <- 'ATAC'

library(randomcoloR)
n <- length(lvl_label)
group.colors <- distinctColorPalette(n)
pie(rep(1,n), col=group.colors)
group.colors.bmmc <- readRDS('../../color.palette.rds')
group.colors[8] <- group.colors.bmmc['GEX']
group.colors[9] <- group.colors.bmmc['ATAC']
names(group.colors) <- lvl_label


umap_df$Label <- factor(umap_df$Label,
                        levels = lvl_label)



gp <- ggplot(umap_df[,2:5], aes(x=X0, y=X1, color=Label)) +
  geom_point(size=0.5, alpha=0.75) + 
  theme_bw()+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=group.colors)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 3, alpha=0.75)))+
  facet_grid(~dataset, switch='y')

gp # 8 x 4












