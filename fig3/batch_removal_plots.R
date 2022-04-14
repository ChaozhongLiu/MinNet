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
lvl_order <- c("SiaNN","GLUE","bindSC","Seurat v3",
               "Liger")
color_palette <- c(
  "#7cef61", #SiaNN
  "#d72628", #GLUE
  "#1f77b4", #bindSC
  "#9566bd", #seurat v3
  "#fe7f0b") #Liger

group.colors <- readRDS('../fig2/color.palette.rds')

#Silhouette score----
sil_score_1 <- read.csv('../fig2/BMMC_data/results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score_1)[1] <- 'Method'
sil_score_1[sil_score_1$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_1 <- sil_score_1[sil_score_1$Method %in% c('SiaNN','GLUE','Liger','Seurat v3','bindSC'),]
sil_score_1$Method <- factor(sil_score_1$Method, levels = lvl_order)
sil_score_1 <- sil_score_1[,c(1,2,4)]
sil_score_1$Dataset <- '1 - Integration: s1d2 & s3d7'

sil_score_2 <- read.csv('results/silh_score.s3s4.csv',stringsAsFactors = F)
colnames(sil_score_2)[1] <- 'Method'
sil_score_2[sil_score_2$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_2$Method <- factor(sil_score_2$Method, levels = lvl_order)
sil_score_2$Dataset <- '2 - GEX:s3d7 & ATAC:s4d1'

sil_score_3 <- read.csv('results/silh_score.s4s3.csv',stringsAsFactors = F)
colnames(sil_score_3)[1] <- 'Method'
sil_score_3[sil_score_3$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_3$Method <- factor(sil_score_3$Method, levels = lvl_order)
sil_score_3$Dataset <- '3 - GEX:s4d1 & ATAC:s3d7'

sil_score <- rbind(sil_score_1,sil_score_2,sil_score_3)
#sil_score[sil_score$modality_score<0.95,'modality_score'] <- 0.95
df_summarise <- sil_score %>%
  group_by(Dataset, Method) %>%
  summarise(
    cell_type_score_sd = safe_sd(cell_type_score),
    cell_type_score = mean(cell_type_score),
    batch_score_sd = safe_sd(batch_score),
    batch_score = mean(batch_score)
  ) %>%
  as.data.frame()

#df_summarise$Dataset <- factor(df_summarise$Dataset,levels = c('Integration - s1d2 & s3d7',
#                                                               'GEX:s3d7 - ATAC:s4d1',
#                                                               'GEX:s4d1 - ATAC:s3d7'))
gp <- ggplot(
  data = sil_score, mapping = aes(
    x = cell_type_score,
    y = batch_score,
    color = Method
  )
) +
  geom_point(data = sil_score, alpha = 0.6, size = 4) +
  scale_color_manual(name = "Method", 
                     values = color_palette) + 
  scale_x_continuous(name = "Cell type silhouette score") +
  facet_wrap(~Dataset) +
  scale_y_continuous(name = "1 - batch silhouette score",limits = c(0.75, 1.0)) +
  theme_bw()

gp # 


# label transfer accuracy bar plot ----
# transfer_acc_1 <- read.csv('../fig2/BMMC_data/results/label_trasfer_acc.csv',stringsAsFactors = F)[,c('method','total')]
# colnames(transfer_acc_1)[1] <- 'Method'
# transfer_acc_1[transfer_acc_1$Method == 'Seurat','Method'] <- 'Seurat v3'
# transfer_acc_1[transfer_acc_1$Method == 'glue','Method'] <- 'GLUE'
# transfer_acc_1 <- transfer_acc_1[transfer_acc_1$Method %in% c('SiaNN','GLUE','Liger','Seurat v3','bindSC'),]
# transfer_acc_1$Method <- factor(transfer_acc_1$Method, levels = lvl_order)
# transfer_acc_1$Dataset <- 's1d2 & s3d7 integration'

transfer_acc_2 <- read.csv('results/label_trasfer_acc.s3s4.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_2)[1] <- 'Method'
transfer_acc_2[transfer_acc_2$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_2$Method <- factor(transfer_acc_2$Method, levels = lvl_order)
transfer_acc_2$Dataset <- '1 - GEX:s3d7 & ATAC:s4d1'

transfer_acc_3 <- read.csv('results/label_trasfer_acc.s4s3.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_3)[1] <- 'Method'
transfer_acc_3[transfer_acc_3$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_3$Method <- factor(transfer_acc_3$Method, levels = lvl_order)
transfer_acc_3$Dataset <- '2 - GEX:s4d1 & ATAC:s3d7'

transfer_acc <- rbind(transfer_acc_2,transfer_acc_3)

df_summarise <- transfer_acc %>%
  group_by(Dataset,Method) %>%
  summarise(
    total_sd = safe_sd(total),
    total = mean(total)
  ) %>%
  as.data.frame()

gp <- ggplot(data = df_summarise, mapping = aes(x = Method, y = total, fill = Method)) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8, position = position_dodge2(width = 0.7)) + 
  geom_errorbar(
    mapping = aes(ymin = total - total_sd, ymax = total + total_sd),
    width = 0.2, position = position_dodge(width = 0.7)
  ) +
  geom_point(
    data = transfer_acc, size = 0.2,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7),
    show.legend = FALSE
  ) +
  scale_fill_manual(name = "Method", values = color_palette) +
  #scale_x_discrete(name = "Dataset") +
  facet_wrap(~Dataset) +
  scale_y_continuous(name = "Label transfer accuracy") +
  theme_bw()

gp


# UMAP batch integration ----
meta <- read.csv('../fig2/BMMC_data/bmmc.rna.meta.csv',
                 stringsAsFactors = F)
file_list <- list.files(path = "../fig2/BMMC_data/results/raw/", 
                        pattern = '_UMAP.txt')
file_list <- file_list[c(1,2,4,5,6)]

umap_df <- read.csv(paste0('../fig2/BMMC_data/results/raw/',file_list[1]),
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:5) {
  umap_df_tmp <- read.csv(paste0('../fig2/BMMC_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('../fig2/BMMC_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:22] <- lvl_label[c(1:21,24)]
lvl_label[23] <- 's1d2'
lvl_label[24] <- 's3d7'

group.colors.samp <- group.colors[names(group.colors) %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df$method <- factor(umap_df$method, levels = c('SiaNN','GLUE','bindSC',
                                                    'Seurat v3','Liger'))

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.025, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')
  
gp


# UMAP real case test s3s4 ----
meta_2 <- read.csv('bmmc.meta.s3s4.csv',
                   stringsAsFactors = F)
file_list_2 <- list.files(path = "results/raw/", 
                        pattern = '_UMAP.s3s4.txt')

umap_df <- read.csv(paste0('results/raw/',file_list_2[1]),
                    stringsAsFactors = F)
all(rownames(meta_2) == umap_df$X)
umap_df$Label <- meta_2$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]


for (i in 2:5) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list_2[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta_2$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list_2[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list_2[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta_2$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list_2[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:21] <- lvl_label[c(1:20,23)]
lvl_label[22] <- 's3d7'
lvl_label[23] <- 's4d1'
group.colors.samp <- group.colors[names(group.colors) %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df$method <- factor(umap_df$method, levels = c('SiaNN','GLUE','bindSC',
                                                    'Seurat v3','Liger'))

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.025, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=2,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp


# UMAP real case test s4s3 ----
meta_3 <- read.csv('bmmc.meta.s4s3.csv',
                   stringsAsFactors = F)
file_list_3 <- list.files(path = "results/raw/", 
                          pattern = '_UMAP.s4s3.txt')

umap_df <- read.csv(paste0('results/raw/',file_list_3[1]),
                    stringsAsFactors = F)
all(rownames(meta_3) == umap_df$X)
umap_df$Label <- meta_3$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]


for (i in 2:5) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list_3[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta_3$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list_3[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('results/raw/',file_list_3[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta_3$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list_3[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:21] <- lvl_label[c(1:20,23)]
lvl_label[22] <- 's4d1'
lvl_label[23] <- 's3d7'
group.colors.samp <- group.colors[names(group.colors) %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df$method <- factor(umap_df$method, levels = c('SiaNN','GLUE','bindSC',
                                                    'Seurat v3','Liger'))

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.025, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=2,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp



