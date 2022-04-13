setwd('/mnt/hdd/chaozhong/manuscript/fig2/PBMC_data/')
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)

safe_sd <- function(x) {
  if (length(x) == 1)
    return(0.0)
  return(sd(x))
}

#color palette ----
lvl_order <- c("SiaNN",
               "bindSC",
               "Seurat v3",
               "Liger")
color_palette <- c("#7cef61", #SiaNN
                   "#1f77b4", #bindSC
                   "#9566bd", #seurat v3
                   "#fe7f0b") #Liger

# Silhouette Score ----
sil_score <- read.csv('results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score)[1] <- 'Method'
sil_score[sil_score$Method == 'Seurat (UMAP)','Method'] <- 'Seurat v3'
sil_score$Method <- factor(sil_score$Method, levels = lvl_order)


sil_score[sil_score$modality_score<0.95,'modality_score'] <- 0.95
df_summarise <- sil_score %>%
  group_by(Method) %>%
  summarise(
    cell_type_score_sd = safe_sd(cell_type_score),
    cell_type_score = mean(cell_type_score),
    modality_score_sd = safe_sd(modality_score),
    modality_score = mean(modality_score)
  ) %>%
  as.data.frame()

gp <- ggplot(
  data = sil_score, mapping = aes(
    x = cell_type_score,
    y = modality_score,
    color = Method
  )
) +
  geom_point(data = sil_score, alpha = 0.6, size = 4) +
  scale_color_manual(name = "Method", 
                     values = color_palette) + 
  scale_x_continuous(name = "Cell type silhouette score") +
  #facet_wrap(~Dataset) +
  scale_y_continuous(name = "1 - modality silhouette score",limits = c(0.95, 1.0)) +
  theme_bw()

gp


#FOSCTTM score ----
fos_score <- read.csv('results/foscttm_score.csv',stringsAsFactors = F)
colnames(fos_score)[1] <- 'Method'
fos_score[fos_score$Method == 'Seurat','Method'] <- 'Seurat v3'
fos_score[fos_score$Method == 'siaNN','Method'] <- 'SiaNN'
fos_score$Method <- factor(fos_score$Method, levels = lvl_order)

df_summarise <- fos_score %>%
  group_by(Method) %>%
  summarise(
    foscttm_sd = safe_sd(foscttm),
    foscttm = mean(foscttm)
  ) %>%
  as.data.frame()


gp <- ggplot(data = df_summarise, mapping = aes(x = Method, y = foscttm, fill = Method)) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8, position = position_dodge2(width = 0.7)) + 
  geom_errorbar(
    mapping = aes(ymin = foscttm - foscttm_sd, ymax = foscttm + foscttm_sd),
    width = 0.2, position = position_dodge(width = 0.7)
  ) +
  geom_point(
    data = fos_score, size = 0.2,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7),
    show.legend = FALSE
  ) +
  scale_fill_manual(name = "Method", values = color_palette) +
  #facet_wrap(~Dataset) +
  #scale_x_discrete(name = "Dataset") +
  scale_y_continuous(name = "FOSCTTM") +
  theme_bw()

gp


# label transfer accuracy heatmap ----
#10X - 1
transfer_acc <- read.csv('results/label_trasfer_acc.csv',stringsAsFactors = F)
colnames(transfer_acc)[1] <- 'Method'
colnames(transfer_acc)[2] <- 'Total'
transfer_acc[transfer_acc$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc$Method <- factor(transfer_acc$Method, levels = rev(lvl_order))

df_summarise <- transfer_acc %>%
  group_by(Method) %>%
  summarise_at(
    colnames(transfer_acc)[2:21],
    mean
    #Total = mean(Total)
  ) %>%
  as.data.frame()

rownames(df_summarise) <- df_summarise$Method
df_summarise <- df_summarise[,2:21]
df_summarise_1 <- as.data.frame(melt(df_summarise))
df_summarise_1$Method <- rep(rownames(df_summarise),20)
df_summarise_1$Method <- factor(df_summarise_1$Method, levels = rev(lvl_order))

gp <- ggplot(df_summarise_1, aes(variable, Method, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Accuracy") +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  xlab("Cell type") +
  #facet_wrap(~Dataset) +
  coord_fixed()

gp


# UMAP ----
meta <- read.csv('data/rna.meta.csv', sep='\t',
                 stringsAsFactors = F)
file_list <- c('Seurat_umap.csv', 'bindSC_UMAP.txt',
               'SiaNN_UMAP.txt', 'Liger_UMAP.1.txt')

umap_df <- read.csv(paste0('results/raw/',file_list[1]),
                    stringsAsFactors = F)
umap_df$X <- gsub('RNA_','',umap_df$X)
umap_df$X <- gsub('ATAC_','',umap_df$X)

all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$seurat_annotations
umap_df$dataset <- 'Cell type'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:4) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$seurat_annotations
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  colnames(umap_df_tmp) <- colnames(umap_df)
  umap_df <- rbind(umap_df,umap_df_tmp)
}


for (i in 1:4) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- c(rep('GEX',10412),rep('ATAC',10412))
  umap_df_tmp$dataset <- 'Modality'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  colnames(umap_df_tmp) <- colnames(umap_df)
  umap_df <- rbind(umap_df,umap_df_tmp)
}
umap_df$dataset <- factor(umap_df$dataset, levels = c('Cell type',
                                                      'Modality'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:19] <- lvl_label[c(2:11,13:21)]
lvl_label[20] <- 'GEX'
lvl_label[21] <- 'ATAC'

library(randomcoloR)
n <- length(lvl_label)
group.colors <- distinctColorPalette(n)
pie(rep(1,n), col=group.colors)
group.colors.bmmc <- readRDS('../color.palette.rds')
group.colors[20] <- group.colors.bmmc['GEX']
group.colors[21] <- group.colors.bmmc['ATAC']
names(group.colors) <- lvl_label
saveRDS(group.colors, 'color.palette.rds')

umap_df$Label <- factor(umap_df$Label,
                        levels = lvl_label)
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df$method <- factor(umap_df$method, levels = lvl_order)

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.005, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp # 25x11






