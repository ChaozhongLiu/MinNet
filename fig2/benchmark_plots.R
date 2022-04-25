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
lvl_order <- c("MinNet","GLUE","bindSC","Seurat v3",
               "Liger", "Liger (UMAP)","Online iNMF",
               "Online iNMF (UMAP)","Liger (UINMF)")
color_palette <- c("#7cef61", #MinNet
                   "#d72628", #GLUE
                   "#1f77b4", #bindSC
                   "#9566bd", #seurat v3
                   "#fe7f0b", #Liger
                   "#feb776", #Liger UMAP
                   "#ac5b04", #Liger online
                   "#785a3b", #Liger online UMAP
                   "#ebe343") #Liger UINMF

meta_1 <- read.csv('BMMC_data/bmmc.rna.meta.csv',
                 stringsAsFactors = F)
meta_1 <- unique(meta_1$cell_type)
meta_2 <- read.csv('BMMC_test/bmmc.rna.meta.csv',
                 stringsAsFactors = F)
meta_2 <- unique(meta_2$cell_type)

meta_3 <- read.csv('Cite_data/cite.rna.meta.csv',
                 stringsAsFactors = F)
meta_3 <- unique(meta_3$cell_type)

meta_4 <- read.csv('Cite_data_2/cite.rna.meta.csv',
                 stringsAsFactors = F)
meta_4 <- unique(meta_4$cell_type)

meta_labels <- unique(c(meta_1,meta_2,meta_3,meta_4))
meta_labels <- c(meta_labels, 'GEX', 'ATAC', 'ADT','s1d2','s3d7','s4d1','s4d8','s4d9')

# library(randomcoloR)
# n <- length(meta_labels) - 8
# group.colors <- distinctColorPalette(n)
# pie(rep(1,n), col=group.colors)
# group.colors <- c(group.colors, group.colors[1:8])
# names(group.colors) <- meta_labels
# 
# saveRDS(group.colors, 'color.palette.rds')
group.colors <- readRDS('color.palette.rds')

#silhouette score----
sil_score_1 <- read.csv('BMMC_data/results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score_1)[1] <- 'Method'
sil_score_1[sil_score_1$Method == 'SiaNN','Method'] <- 'MinNet'
sil_score_1[sil_score_1$Method == 'LigerOnline','Method'] <- 'Online iNMF'
sil_score_1[sil_score_1$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
sil_score_1[sil_score_1$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_1$Method <- factor(sil_score_1$Method, levels = lvl_order)
sil_score_1$Dataset <- '10X Multiome - 1'

sil_score_2 <- read.csv('BMMC_test/results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score_2)[1] <- 'Method'
sil_score_2[sil_score_2$Method == 'SiaNN','Method'] <- 'MinNet'
sil_score_2[sil_score_2$Method == 'LigerOnline','Method'] <- 'Online iNMF'
sil_score_2[sil_score_2$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
sil_score_2[sil_score_2$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_2$Method <- factor(sil_score_2$Method, levels = lvl_order)
sil_score_2$Dataset <- '10X Multiome - 2'

sil_score_3 <- read.csv('Cite_data/results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score_3)[1] <- 'Method'
sil_score_3[sil_score_3$Method == 'SiaNN','Method'] <- 'MinNet'
sil_score_3[sil_score_3$Method == 'Liger','Method'] <- 'Liger (UINMF)'
sil_score_3[sil_score_3$Method == 'LigerINMF','Method'] <- 'Liger'
sil_score_3[sil_score_3$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_3$Method <- factor(sil_score_3$Method, levels = lvl_order)
sil_score_3$Dataset <- 'Cite-seq - 1'

sil_score_4 <- read.csv('Cite_data_2/results/silh_score.csv',stringsAsFactors = F)
colnames(sil_score_4)[1] <- 'Method'
sil_score_4[sil_score_4$Method == 'SiaNN','Method'] <- 'MinNet'
sil_score_4[sil_score_4$Method == 'Liger','Method'] <- 'Liger (UINMF)'
sil_score_4[sil_score_4$Method == 'LigerINMF','Method'] <- 'Liger'
sil_score_4[sil_score_4$Method == 'Seurat','Method'] <- 'Seurat v3'
sil_score_4$Method <- factor(sil_score_4$Method, levels = lvl_order)
sil_score_4$Dataset <- 'Cite-seq - 2'



sil_score <- rbind(sil_score_1,sil_score_2,sil_score_3,sil_score_4)
sil_score[sil_score$modality_score<0.95,'modality_score'] <- 0.95
df_summarise <- sil_score %>%
  group_by(Dataset, Method) %>%
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
  facet_wrap(~Dataset) +
  scale_y_continuous(name = "1 - modality silhouette score",limits = c(0.95, 1.0)) +
  theme_bw()

gp



# ARI score ----
ari_df_1 <- read.csv('BMMC_data/results/ARI_score.csv', stringsAsFactors = F)
ari_df_1 <- ari_df_1[ari_df_1$X.Cluster <= 25,]
ari_df_1 <- ari_df_1[ari_df_1$X.Cluster > 7,]
ari_df_1[ari_df_1$Method == 'SiaNN','Method'] <- 'MinNet'
ari_df_1[ari_df_1$Method == 'Seurat','Method'] <- 'Seurat v3'
ari_df_1$Method <- factor(ari_df_1$Method, levels = lvl_order)

df_summarise_1 <- ari_df_1 %>%
  group_by(Method, X.Cluster) %>%
  summarise(
    ARI_score = mean(ARI)
  ) %>%
  as.data.frame()
df_summarise_1$Dataset <- '10X Multiome - 1'


ari_df_2 <- read.csv('BMMC_test/results/ARI_score.csv', stringsAsFactors = F)
ari_df_2 <- ari_df_2[ari_df_2$X.Cluster <= 25,]
ari_df_2 <- ari_df_2[ari_df_2$X.Cluster > 7,]
ari_df_2[ari_df_2$Method == 'SiaNN','Method'] <- 'MinNet'
ari_df_2[ari_df_2$Method == 'Seurat','Method'] <- 'Seurat v3'
ari_df_2$Method <- factor(ari_df_2$Method, levels = lvl_order)

df_summarise_2 <- ari_df_2 %>%
  group_by(Method, X.Cluster) %>%
  summarise(
    ARI_score = mean(ARI)
  ) %>%
  as.data.frame()
df_summarise_2$Dataset <- '10X Multiome - 2'

ari_df_3 <- read.csv('Cite_data/results/ARI_score.csv', stringsAsFactors = F)
ari_df_3 <- ari_df_3[ari_df_3$X.Cluster <= 25,]
ari_df_3 <- ari_df_3[ari_df_3$X.Cluster > 7,]
ari_df_3[ari_df_3$Method == 'SiaNN','Method'] <- 'MinNet'
ari_df_3[ari_df_3$Method == 'Seurat','Method'] <- 'Seurat v3'
ari_df_3[ari_df_3$Method == 'Liger(UINMF)','Method'] <- 'Liger (UINMF)'
ari_df_3[ari_df_3$Method == 'LigerINMF','Method'] <- 'Liger'
ari_df_3$Method <- factor(ari_df_3$Method, levels = lvl_order)

df_summarise_3 <- ari_df_3 %>%
  group_by(Method, X.Cluster) %>%
  summarise(
    ARI_score = mean(ARI)
  ) %>%
  as.data.frame()
df_summarise_3$Dataset <- 'Cite-seq - 1'

ari_df_4 <- read.csv('BMMC_data/results/ARI_score.csv', stringsAsFactors = F)
ari_df_4 <- ari_df_4[ari_df_4$X.Cluster <= 25,]
ari_df_4 <- ari_df_4[ari_df_4$X.Cluster > 7,]
ari_df_4[ari_df_4$Method == 'SiaNN','Method'] <- 'MinNet'
ari_df_4[ari_df_4$Method == 'Seurat','Method'] <- 'Seurat v3'
ari_df_4[ari_df_4$Method == 'Liger','Method'] <- 'Liger (UINMF)'
ari_df_4[ari_df_4$Method == 'LigerINMF','Method'] <- 'Liger'
ari_df_4$Method <- factor(ari_df_4$Method, levels = lvl_order)

df_summarise_4 <- ari_df_4 %>%
  group_by(Method, X.Cluster) %>%
  summarise(
    ARI_score = mean(ARI)
  ) %>%
  as.data.frame()
df_summarise_4$Dataset <- 'Cite-seq - 2'

ARI_score <- rbind(df_summarise_1,df_summarise_2,df_summarise_3,df_summarise_4)
#ARI_score[ARI_score$modality_score<0.95,'modality_score'] <- 0.95


p<-ggplot(ARI_score, aes(x=X.Cluster, y=ARI_score, group=Method)) +
  geom_line(aes(color=Method)) +
  geom_point(aes(color=Method)) + 
  scale_color_manual(name = "Method", 
                     values = color_palette) + 
  scale_x_continuous(name = "Number of Clusters") +
  facet_wrap(~Dataset) +
  scale_y_continuous(name = "Adjusted Rand Index") +
  theme_bw()
p




#FOSCTTM score ----
fos_score_1 <- read.csv('BMMC_data/results/foscttm_score.csv',stringsAsFactors = F)
colnames(fos_score_1)[1] <- 'Method'
fos_score_1[fos_score_1$Method == 'SiaNN','Method'] <- 'MinNet'
fos_score_1[fos_score_1$Method == 'LigerOnline','Method'] <- 'Online iNMF'
fos_score_1[fos_score_1$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
fos_score_1[fos_score_1$Method == 'Seurat','Method'] <- 'Seurat v3'
fos_score_1[fos_score_1$Method == 'glue','Method'] <- 'GLUE'
fos_score_1$Method <- factor(fos_score_1$Method, levels = lvl_order)
fos_score_1$Dataset <- '10X Multiome'

fos_score_2 <- read.csv('BMMC_test/results/foscttm_score.csv',stringsAsFactors = F)
colnames(fos_score_2)[1] <- 'Method'
fos_score_2[fos_score_2$Method == 'SiaNN','Method'] <- 'MinNet'
fos_score_2[fos_score_2$Method == 'LigerOnline','Method'] <- 'Online iNMF'
fos_score_2[fos_score_2$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
fos_score_2[fos_score_2$Method == 'seurat','Method'] <- 'Seurat v3'
fos_score_2[fos_score_2$Method == 'glue','Method'] <- 'GLUE'
fos_score_2$Method <- factor(fos_score_2$Method, levels = lvl_order)
fos_score_2$Dataset <- '10X Multiome'

fos_score_3 <- read.csv('Cite_data/results/foscttm_score.csv',stringsAsFactors = F)
colnames(fos_score_3)[1] <- 'Method'
fos_score_3[fos_score_3$Method == 'SiaNN','Method'] <- 'MinNet'
fos_score_3[fos_score_3$Method == 'Liger','Method'] <- 'Liger (UINMF)'
fos_score_3[fos_score_3$Method == 'LigerINMF','Method'] <- 'Liger'
fos_score_3[fos_score_3$Method == 'Seurat','Method'] <- 'Seurat v3'
fos_score_3$Method <- factor(fos_score_3$Method, levels = lvl_order)
fos_score_3$Dataset <- 'Cite-seq'

fos_score_4 <- read.csv('Cite_data_2/results/foscttm_score.csv',stringsAsFactors = F)
colnames(fos_score_4)[1] <- 'Method'
fos_score_4[fos_score_4$Method == 'SiaNN','Method'] <- 'MinNet'
fos_score_4[fos_score_4$Method == 'Liger','Method'] <- 'Liger (UINMF)'
fos_score_4[fos_score_4$Method == 'LigerINMF','Method'] <- 'Liger'
fos_score_4[fos_score_4$Method == 'Seurat','Method'] <- 'Seurat v3'
fos_score_4$Method <- factor(fos_score_4$Method, levels = lvl_order)
fos_score_4$Dataset <- 'Cite-seq'

fos_score <- rbind(fos_score_1,fos_score_2,fos_score_3,fos_score_4)

df_summarise <- fos_score %>%
  group_by(Dataset,batch,Method) %>%
  summarise(
    foscttm_sd = safe_sd(foscttm),
    foscttm = mean(foscttm)
  ) %>%
  as.data.frame()


gp <- ggplot(data = df_summarise, mapping = aes(x = batch, y = foscttm, fill = Method)) + 
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
  facet_wrap(~Dataset) +
  #scale_x_discrete(name = "Dataset") +
  scale_y_continuous(name = "FOSCTTM") +
  theme_bw()

gp


# label transfer accuracy bar plot ----
transfer_acc_1 <- read.csv('BMMC_data/results/label_trasfer_acc.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_1)[1] <- 'Method'
transfer_acc_1[transfer_acc_1$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc_1[transfer_acc_1$Method == 'LigerOnline','Method'] <- 'Online iNMF'
transfer_acc_1[transfer_acc_1$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
transfer_acc_1[transfer_acc_1$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_1[transfer_acc_1$Method == 'glue','Method'] <- 'GLUE'
transfer_acc_1$Method <- factor(transfer_acc_1$Method, levels = lvl_order)
transfer_acc_1$Dataset <- '10X Multiome - 1'


transfer_acc_2 <- read.csv('BMMC_test/results/label_trasfer_acc.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_2)[1] <- 'Method'
transfer_acc_2[transfer_acc_2$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc_2[transfer_acc_2$Method == 'LigerOnline','Method'] <- 'Online iNMF'
transfer_acc_2[transfer_acc_2$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
transfer_acc_2[transfer_acc_2$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_2[transfer_acc_2$Method == 'glue','Method'] <- 'GLUE'
transfer_acc_2$Method <- factor(transfer_acc_2$Method, levels = lvl_order)
transfer_acc_2$Dataset <- '10X Multiome - 2'

transfer_acc_3 <- read.csv('Cite_data/results/label_trasfer_acc.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_3)[1] <- 'Method'
transfer_acc_3[transfer_acc_3$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc_3[transfer_acc_3$Method == 'Liger','Method'] <- 'Liger (UINMF)'
transfer_acc_3[transfer_acc_3$Method == 'LigerINMF','Method'] <- 'Liger'
transfer_acc_3[transfer_acc_3$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_3$Method <- factor(transfer_acc_3$Method, levels = lvl_order)
transfer_acc_3$Dataset <- 'Cite-seq - 1'

transfer_acc_4 <- read.csv('Cite_data_2/results/label_trasfer_acc.csv',stringsAsFactors = F)[,c('method','total')]
colnames(transfer_acc_4)[1] <- 'Method'
transfer_acc_4[transfer_acc_4$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc_4[transfer_acc_4$Method == 'Liger','Method'] <- 'Liger (UINMF)'
transfer_acc_4[transfer_acc_4$Method == 'LigerINMF','Method'] <- 'Liger'
transfer_acc_4[transfer_acc_4$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc_4$Method <- factor(transfer_acc_4$Method, levels = lvl_order)
transfer_acc_4$Dataset <- 'Cite-seq - 2'

transfer_acc <- rbind(transfer_acc_1,transfer_acc_2,transfer_acc_3,transfer_acc_4)

df_summarise <- transfer_acc %>%
  group_by(Dataset,Method) %>%
  summarise(
    total_sd = safe_sd(total),
    total = mean(total)
  ) %>%
  as.data.frame()


gp <- ggplot(data = df_summarise, mapping = aes(x = Dataset, y = total, fill = Method)) + 
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
  scale_y_continuous(name = "Label transfer accuracy") +
  theme_bw()

gp


# label transfer accuracy heatmap ----
#10X - 1
transfer_acc <- read.csv('BMMC_data/results/label_trasfer_acc.csv',stringsAsFactors = F)
colnames(transfer_acc)[1] <- 'Method'
colnames(transfer_acc)[2] <- 'Total'
transfer_acc[transfer_acc$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc[transfer_acc$Method == 'LigerOnline','Method'] <- 'Online iNMF'
transfer_acc[transfer_acc$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
transfer_acc[transfer_acc$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc[transfer_acc$Method == 'glue','Method'] <- 'GLUE'
transfer_acc$Method <- factor(transfer_acc$Method, levels = rev(lvl_order))

df_summarise <- transfer_acc %>%
  group_by(Method) %>%
  summarise_at(
    colnames(transfer_acc)[2:24],
    mean
    #Total = mean(Total)
  ) %>%
  as.data.frame()

rownames(df_summarise) <- df_summarise$Method
df_summarise <- df_summarise[,2:24]
df_summarise_1 <- as.data.frame(melt(df_summarise))
df_summarise_1$Method <- rep(rownames(df_summarise),23)
df_summarise_1$Method <- factor(df_summarise_1$Method, levels = rev(lvl_order))
#df_summarise_1$variable <- rep(cell_type_label,each=8)
#df_summarise_1$variable <- factor(df_summarise_1$variable, levels = cell_type_label)
df_summarise_1$Dataset <- '10X Multiome - 1'

#10X - 2
transfer_acc_2 <- read.csv('BMMC_test/results/label_trasfer_acc.csv',
                         stringsAsFactors = F)
colnames(transfer_acc)[1] <- 'Method'
colnames(transfer_acc)[2] <- 'Total'
transfer_acc[transfer_acc$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc[transfer_acc$Method == 'LigerOnline','Method'] <- 'Online iNMF'
transfer_acc[transfer_acc$Method == 'LigerOnline (UMAP)','Method'] <- 'Online iNMF (UMAP)'
transfer_acc[transfer_acc$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc[transfer_acc$Method == 'glue','Method'] <- 'GLUE'
transfer_acc$Method <- factor(transfer_acc$Method, levels = rev(lvl_order))

df_summarise <- transfer_acc %>%
  group_by(Method) %>%
  summarise_at(
    colnames(transfer_acc)[2:22],
    mean
    #Total = mean(Total)
  ) %>%
  as.data.frame()

rownames(df_summarise) <- df_summarise$Method
df_summarise <- df_summarise[,2:22]
df_summarise_2 <- as.data.frame(melt(df_summarise))
df_summarise_2$Method <- rep(rownames(df_summarise),21)
df_summarise_2$Method <- factor(df_summarise_2$Method, levels = rev(lvl_order))
#cts <- as.character(unique(df_summarise_2$variable))
#df_summarise_2$variable <- rep(cell_type_label,each=8)
#df_summarise_2$variable <- factor(df_summarise_2$variable, levels = cell_type_label)
df_summarise_2$Dataset <- '10X Multiome - 2'

df_summarise <- rbind(df_summarise_1,df_summarise_2)
cts <- levels(df_summarise$variable)
cts <- gsub('\\.\\.','+ ',cts)
cts <- gsub('\\.',' ',cts)
cts <- gsub('MK E prog','MK/E prog',cts)
cts <- gsub('ID2 hi myeloid prog','ID2 hi-myeloid prog',cts)
cts <- gsub('G M prog','G/M prog',cts)
cts_old <- levels(df_summarise$variable)
df_summarise$variable <- as.character(df_summarise$variable)
for (i in 1:23) {
  df_summarise[df_summarise$variable == cts_old[i], 'variable'] <- cts[i]
}
df_summarise$variable <- factor(df_summarise$variable, levels=cts)
# Create a ggheatmap
gp <- ggplot(df_summarise, aes(variable, Method, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Accuracy") +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  xlab("Cell type") +
  facet_wrap(~Dataset) +
  coord_fixed()
# Print the heatmap
gp


#Cite - 1
transfer_acc <- read.csv('Cite_data/results/label_trasfer_acc.csv',stringsAsFactors = F)
colnames(transfer_acc)[1] <- 'Method'
colnames(transfer_acc)[2] <- 'Total'
transfer_acc[transfer_acc$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc[transfer_acc$Method == 'Liger','Method'] <- 'Liger (UINMF)'
transfer_acc[transfer_acc$Method == 'LigerINMF','Method'] <- 'Liger'
transfer_acc[transfer_acc$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc$Method <- factor(transfer_acc$Method, levels = rev(lvl_order))

df_summarise <- transfer_acc %>%
  group_by(Method) %>%
  summarise_at(
    colnames(transfer_acc)[2:27],
    mean
    #Total = mean(Total)
  ) %>%
  as.data.frame()

rownames(df_summarise) <- df_summarise$Method
df_summarise <- df_summarise[,2:27]
df_summarise_3 <- as.data.frame(melt(df_summarise))
df_summarise_3$Method <- rep(rownames(df_summarise),26)
df_summarise_3$Method <- factor(df_summarise_3$Method, levels = rev(lvl_order))
df_summarise_3$Dataset <- 'Cite-seq - 1'


#Cite - 2
transfer_acc <- read.csv('Cite_data_2/results/label_trasfer_acc.csv',stringsAsFactors = F)
colnames(transfer_acc)[1] <- 'Method'
colnames(transfer_acc)[2] <- 'Total'
transfer_acc[transfer_acc$Method == 'SiaNN','Method'] <- 'MinNet'
transfer_acc[transfer_acc$Method == 'Liger','Method'] <- 'Liger (UINMF)'
transfer_acc[transfer_acc$Method == 'LigerINMF','Method'] <- 'Liger'
transfer_acc[transfer_acc$Method == 'Seurat','Method'] <- 'Seurat v3'
transfer_acc$Method <- factor(transfer_acc$Method, levels = rev(lvl_order))

df_summarise <- transfer_acc %>%
  group_by(Method) %>%
  summarise_at(
    colnames(transfer_acc)[2:25],
    mean
    #Total = mean(Total)
  ) %>%
  as.data.frame()

rownames(df_summarise) <- df_summarise$Method
df_summarise <- df_summarise[,2:25]
df_summarise_4 <- as.data.frame(melt(df_summarise))
df_summarise_4$Method <- rep(rownames(df_summarise),24)
df_summarise_4$Method <- factor(df_summarise_4$Method, levels = rev(lvl_order))
df_summarise_4$Dataset <- 'Cite-seq - 2'


df_summarise <- rbind(df_summarise_3,df_summarise_4)
cts <- levels(df_summarise$variable)
cts <- gsub('\\.\\.','+ ',cts)
cts <- gsub('\\.',' ',cts)
cts <- gsub('MK E prog','MK/E prog',cts)
cts <- gsub('ID2 hi myeloid prog','ID2 hi-myeloid prog',cts)
cts <- gsub('G M prog','G/M prog',cts)
cts_old <- levels(df_summarise$variable)
df_summarise$variable <- as.character(df_summarise$variable)
for (i in 1:26) {
  df_summarise[df_summarise$variable == cts_old[i], 'variable'] <- cts[i]
}
df_summarise$variable <- factor(df_summarise$variable, levels=cts)

# Create a ggheatmap
gp2 <- ggplot(df_summarise, aes(variable, Method, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Accuracy") +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  xlab("Cell type") +
  facet_wrap(~Dataset) +
  coord_fixed()
# Print the heatmap
gp2



# UMAP visualization BMMC_val -----
lvl_order <- c("MinNet","GLUE","bindSC","Seurat v3",
               "Liger","Online iNMF")
color_palette <- c(
  "#7cef61", #MinNet
  "#d72628", #GLUE
  "#1f77b4", #bindSC
  "#9566bd", #seurat v3
  "#fe7f0b", #Liger
  "#ac5b04") #Liger online

meta <- read.csv('BMMC_data/bmmc.rna.meta.csv',
                 stringsAsFactors = F)
file_list <- list.files(path = "BMMC_data/results/raw/", 
                        pattern = '_UMAP.txt')

umap_df <- read.csv(paste0('BMMC_data/results/raw/',file_list[1]),
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- c(rep('GEX',12163),rep('ATAC',12163))
  umap_df_tmp$dataset <- 'Modality'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}
umap_df$dataset <- factor(umap_df$dataset, levels = c('Cell type',
                                                      'Modality',
                                                      'Batch'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:22] <- lvl_label[c(2:10,12:23,26)]
lvl_label[23] <- 'GEX'
lvl_label[24] <- 'ATAC'
lvl_label[25] <- 's1d2'
lvl_label[26] <- 's3d7'
group.colors.samp <- group.colors[meta_labels %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df[umap_df$method == 'SiaNN','method'] <- 'MinNet'
umap_df[umap_df$method == 'LigerOnline','method'] <- 'Online iNMF'
umap_df$method <- factor(umap_df$method, levels = lvl_order)

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.005, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp # 25x11



# UMAP visualization BMMC_test -----
lvl_order <- c("MinNet","GLUE","bindSC","Seurat v3",
               "Liger","Online iNMF")
color_palette <- c(
  "#7cef61", #SiaNN
  "#d72628", #GLUE
  "#1f77b4", #bindSC
  "#9566bd", #seurat v3
  "#fe7f0b", #Liger
  "#ac5b04") #Liger online

meta <- read.csv('BMMC_test/bmmc.rna.meta.csv',
                 stringsAsFactors = F)
file_list <- list.files(path = "BMMC_test/results/raw/", 
                        pattern = '_UMAP.txt')

umap_df <- read.csv(paste0('BMMC_test/results/raw/',file_list[1]),
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_test/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_test/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:6) {
  umap_df_tmp <- read.csv(paste0('BMMC_test/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- c(rep('GEX',20009),rep('ATAC',20009))
  umap_df_tmp$dataset <- 'Modality'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}
umap_df$dataset <- factor(umap_df$dataset, levels = c('Cell type',
                                                      'Modality',
                                                      'Batch'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:20] <- lvl_label[c(2:10,12:21,25)]
lvl_label[21] <- 'GEX'
lvl_label[22] <- 'ATAC'
lvl_label[23] <- 's4d1'
lvl_label[24] <- 's4d8'
lvl_label[25] <- 's4d9'

group.colors.samp <- group.colors[meta_labels %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))

umap_df[umap_df$method == 'SiaNN','method'] <- 'MinNet'
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df[umap_df$method == 'LigerOnline','method'] <- 'Online iNMF'
umap_df$method <- factor(umap_df$method, levels = lvl_order)

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.005, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp



# UMAP visualization Cite - 1 -----
lvl_order <- c("MinNet","bindSC","Seurat v3",
               "Liger","Liger (UINMF)")
color_palette <- c(
  "#7cef61", #SiaNN
  "#1f77b4", #bindSC
  "#9566bd", #seurat v3
  "#fe7f0b", #Liger
  "#ebe343") #Liger UINMF

meta <- read.csv('Cite_data/cite.rna.meta.csv',
                 stringsAsFactors = F)
file_list <- list.files(path = "Cite_data/results/raw/", 
                        pattern = '_UMAP.txt')

umap_df <- read.csv(paste0('Cite_data/results/raw/',file_list[1]),
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- c(rep('GEX',14826),rep('ADT',14826))
  umap_df_tmp$dataset <- 'Modality'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}
umap_df$dataset <- factor(umap_df$dataset, levels = c('Cell type',
                                                      'Modality',
                                                      'Batch'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:25] <- lvl_label[c(2:10,12:24,27:29)]
lvl_label[26] <- 'GEX'
lvl_label[27] <- 'ADT'
lvl_label[28] <- 's1d2'
lvl_label[29] <- 's3d7'

group.colors.samp <- group.colors[meta_labels %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))

umap_df[umap_df$method == 'SiaNN','method'] <- 'MinNet'
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df[umap_df$method == 'Liger','method'] <- 'Liger (UINMF)'
umap_df[umap_df$method == 'LigerINMF','method'] <- 'Liger'

umap_df$method <- factor(umap_df$method, levels = lvl_order)

gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.005, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 5,alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp #21x11



# UMAP visualization Cite - 2 -----
lvl_order <- c("MinNet","bindSC","Seurat v3",
               "Liger","Liger (UINMF)")
color_palette <- c(
  "#7cef61", #SiaNN
  "#1f77b4", #bindSC
  "#9566bd", #seurat v3
  "#fe7f0b", #Liger
  "#ebe343") #Liger UINMF

meta <- read.csv('Cite_data_2/cite.rna.meta.csv',
                 stringsAsFactors = F)
file_list <- list.files(path = "Cite_data_2/results/raw/", 
                        pattern = '_UMAP.txt')

umap_df <- read.csv(paste0('Cite_data_2/results/raw/',file_list[1]),
                    stringsAsFactors = F)
all(rownames(meta) == umap_df$X)
umap_df$Label <- meta$batch
umap_df$dataset <- 'Batch'
umap_df$method <- strsplit(file_list[1],'_')[[1]][1]

for (i in 2:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data_2/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$batch
  umap_df_tmp$dataset <- 'Batch'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data_2/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- meta$cell_type
  umap_df_tmp$dataset <- 'Cell type'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}

for (i in 1:5) {
  umap_df_tmp <- read.csv(paste0('Cite_data_2/results/raw/',file_list[i]),
                          stringsAsFactors = F)
  umap_df_tmp$Label <- c(rep('GEX',15066),rep('ADT',15066))
  umap_df_tmp$dataset <- 'Modality'
  umap_df_tmp$method <- strsplit(file_list[i],'_')[[1]][1]
  umap_df <- rbind(umap_df,umap_df_tmp)
}
umap_df$dataset <- factor(umap_df$dataset, levels = c('Cell type',
                                                      'Modality',
                                                      'Batch'))
lvl_label <- sort(unique(umap_df$Label))
lvl_label[1:23] <- lvl_label[c(2:9,11:23,27,28)]
lvl_label[24] <- 'GEX'
lvl_label[25] <- 'ADT'
lvl_label[26] <- 's4d1'
lvl_label[27] <- 's4d8'
lvl_label[28] <- 's4d9'
group.colors.samp <- group.colors[meta_labels %in% lvl_label]
umap_df$Label <- factor(umap_df$Label,
                        levels = names(group.colors.samp))

umap_df[umap_df$method == 'SiaNN','method'] <- 'MinNet'
umap_df[umap_df$method == 'Seurat','method'] <- 'Seurat v3'
umap_df[umap_df$method == 'Liger','method'] <- 'Liger (UINMF)'
umap_df[umap_df$method == 'LigerINMF','method'] <- 'Liger'

umap_df$method <- factor(umap_df$method, levels = lvl_order)


gp <- ggplot(umap_df[,2:6], aes(x=UMAP_1, y=UMAP_2, color=Label)) +
  geom_point(size=0.005, alpha=0.35) + 
  theme_bw()+
  theme(text = element_text(size = 20))+
  scale_color_manual(values=group.colors.samp)+
  guides(color = guide_legend(ncol=1,override.aes = list(size = 5, alpha=0.95)))+
  facet_grid(dataset ~ method, switch='y')

gp
