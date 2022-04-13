setwd('/mnt/hdd/chaozhong/manuscript/fig2/Cite_data_2/seurat/')
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(rhdf5)
library(pheatmap)
library(cowplot)
library(RColorBrewer)

rm(list=ls())

set.seed(14)


# RNA-seq processing ----
count_mtx <- readRDS('cite.rna.counts.rds')
cite.rna <- CreateSeuratObject(counts = count_mtx, 
                               project = "CITE",
                               assay = 'RNA')
metadata.rna <- readRDS('cite.rna.meta.rds')
cite.rna <- AddMetaData(cite.rna, metadata.rna)

cite.rna <- NormalizeData(cite.rna, normalization.method = "LogNormalize", scale.factor = 10000)
cite.rna <- FindVariableFeatures(cite.rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cite.rna)
cite.rna <- ScaleData(cite.rna, features = all.genes)

cite.rna <- RunPCA(cite.rna, features = VariableFeatures(object = cite.rna))
#bmmc.rna <- JackStraw(bmmc.rna, num.replicate = 100)
#bmmc.rna <- ScoreJackStraw(bmmc.rna, dims = 1:20)
#JackStrawPlot(bmmc.rna, dims = 1:15)
#ElbowPlot(bmmc.rna) #20 is better

cite.rna <- FindNeighbors(cite.rna, dims = 1:20)
#bmmc.rna <- FindClusters(bmmc.rna, resolution = 0.1)
cite.rna <- RunUMAP(cite.rna, dims = 1:20)
#rna <- RunTSNE(rna, reduction = 'pca', dims=1:15)
DimPlot(cite.rna, reduction = "umap", group.by='cell_type')

saveRDS(cite.rna,'cite.rna.rds')


# ADT processing ----
count_mtx <- readRDS('cite.adt.counts.rds')
cite.adt <- CreateSeuratObject(counts = count_mtx, 
                               project = "CITE",
                               assay = 'ADT')
metadata.adt <- readRDS('cite.adt.meta.rds')
cite.adt <- AddMetaData(cite.adt, metadata.adt)

cite.adt <- NormalizeData(cite.adt, normalization.method = "LogNormalize", scale.factor = 10000)
cite.adt <- FindVariableFeatures(cite.adt, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cite.adt)
cite.adt <- ScaleData(cite.adt, features = all.genes)

cite.adt <- RunPCA(cite.adt, features = VariableFeatures(object = cite.adt))

cite.adt <- FindNeighbors(cite.adt, dims = 1:20)
#bmmc.rna <- FindClusters(bmmc.rna, resolution = 0.1)
cite.adt <- RunUMAP(cite.adt, dims = 1:20)
#rna <- RunTSNE(rna, reduction = 'pca', dims=1:15)
DimPlot(cite.adt, reduction = "umap", group.by='cell_type')

saveRDS(cite.adt,'cite.adt.rds')



# integration ----
cite.rna <- readRDS('cite.rna.rds')
cite.adt <- readRDS('cite.adt.rds')

gene.use <- intersect(rownames(cite.rna),rownames(cite.adt))
cite.rna.homo <- CreateSeuratObject(counts=cite.rna@assays$RNA@counts[gene.use,])
cite.adt.homo <- CreateSeuratObject(counts=cite.adt@assays$ADT@counts[gene.use,])
metadata.rna <- readRDS('cite.rna.meta.rds')
cite.rna.homo <- AddMetaData(cite.rna.homo, metadata.rna)
metadata.adt <- readRDS('cite.adt.meta.rds')
cite.adt.homo <- AddMetaData(cite.adt.homo, metadata.adt)

cite.rna.homo <- ScaleData(cite.rna.homo, features = gene.use, do.scale = TRUE)
cite.rna.homo <- RunPCA(cite.rna.homo, features = gene.use, verbose = FALSE)

transfer.anchors <- FindTransferAnchors(reference = cite.adt.homo, 
                                        query = cite.rna.homo, 
                                        features = gene.use,  
                                        dims = seq(1,15,1), 
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")

# label transfer ----
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = cite.adt.homo$cell_type,
                                     weight.reduction = cite.rna.homo[["pca"]])

cite.rna.homo <- AddMetaData(cite.rna.homo, metadata = celltype.predictions)
pred_df <- cite.rna.homo@meta.data[,c('cell_type','predicted.id')]
colnames(pred_df) <- c('cell_type','prediction')
write.table(pred_df, '../results/raw/Seurat_prediction.txt',sep=',',quote=F)

predictions <- table(cite.rna.homo$cell_type, cite.rna.homo$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions <- predictions[with(predictions,order(Var1)),]
write.csv(predictions, 'seurat_prediction.csv',quote = F)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", 
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1

mean(cite.rna.homo$cell_type == cite.rna.homo$predicted.id)
predictions <- read.csv('seurat_prediction.csv',stringsAsFactors = F)
acc <- predictions[predictions$Var1 == predictions$Var2,]
acc[with(acc,order(Var1)),]
write.table(acc, '../results/raw/Seurat_acc.csv', quote = F, sep=',',row.names=F)



# co-embedding ----
refdata <- GetAssayData(cite.adt.homo, assay = "RNA", slot = "data")[gene.use, ]

imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata, 
                           weight.reduction = cite.rna.homo[["pca"]], 
                           dims=seq(1,15,1))
tmp <- cite.rna.homo
cite.rna.homo[["RNA"]] <- imputation
coembed <- merge(x =cite.adt.homo, y = cite.rna.homo)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = gene.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = gene.use, verbose = FALSE)

pca_coembed <- coembed@reductions$pca@cell.embeddings
dim(pca_coembed)
all(gsub('_1','_2',rownames(pca_coembed)[1:15066]) == rownames(pca_coembed)[15067:30132])
write.csv(pca_coembed, '../results/raw/Seurat_pca.csv',quote = F)

saveRDS(coembed,'coembed.rds')
saveRDS(cite.rna.homo,'cite.rna.homo.rds')

coembed <- RunUMAP(coembed, dims = 1:30)

coembed$tech <- c(rep(0,15066),rep(1,15066))
DimPlot(coembed, group.by = "cell_type")
DimPlot(coembed, group.by = "batch")
DimPlot(coembed, group.by = "tech")

coembed <- readRDS('coembed.rds')
coembed <- RunUMAP(coembed, dims = 1:20)
umap_coembed <- coembed@reductions$umap@cell.embeddings
dim(umap_coembed)
all(gsub('_1','_2',rownames(umap_coembed)[1:15066]) == rownames(umap_coembed)[15067:30132])
write.csv(umap_coembed, '../results/raw/Seurat_UMAP.txt',quote = F)


coembed <- readRDS('coembed.rds')
pca_coembed <- coembed@reductions$pca@cell.embeddings[,1:20]
dim(pca_coembed)
all(gsub('_1','_2',rownames(pca_coembed)[1:15066]) == rownames(pca_coembed)[15067:30132])
write.csv(pca_coembed, '../results/raw/Seurat_embed.pca.txt',quote = F)





