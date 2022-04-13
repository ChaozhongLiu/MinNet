#setwd('/mnt/hdd/chaozhong/manuscript/fig2/BMMC_test/seurat/')
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

# rds is got from humble:~/thesis/BMMC_test/
# ATAC processing ----
count_mtx <- readRDS('bmmc.atac.counts.rds')
bmmc.atac <- CreateSeuratObject(counts = count_mtx, 
                                project = "BMMC",
                                assay = 'peaks')
metadata.atac <- readRDS('bmmc.atac.meta.rds')
bmmc.atac <- AddMetaData(bmmc.atac, metadata.atac)

peaks <- bmmc.atac@assays$peaks@counts
peaks <- peaks[startsWith(rownames(peaks), 'chr'),]
rownames(peaks) <- gsub('chr','',rownames(peaks))

peaks_tmp <- peaks[,1:10000]
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks_tmp,
                                            annotation.file = "../../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                            seq.levels = c(1:22,"X","Y"), 
                                            upstream = 2000, verbose = TRUE)

peaks_tmp <- peaks[,10001:20009]
activity.matrix_tmp <- CreateGeneActivityMatrix(peak.matrix = peaks_tmp,
                                            annotation.file = "../../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                            seq.levels = c(1:22,"X","Y"), 
                                            upstream = 2000, verbose = TRUE)
activity.matrix <- cbind(activity.matrix,  activity.matrix_tmp)


#creat gene activity assay
bmmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
DefaultAssay(bmmc.atac) <- "ACTIVITY"

bmmc.atac <- FindVariableFeatures(bmmc.atac)
bmmc.atac <- NormalizeData(bmmc.atac)
bmmc.atac <- ScaleData(bmmc.atac)

DefaultAssay(bmmc.atac) <- "peaks"
VariableFeatures(bmmc.atac) <- names(which(Matrix::rowSums(bmmc.atac) > 100))
bmmc.atac <- RunLSI(bmmc.atac, n = 50, scale.max = NULL)
bmmc.atac <- RunUMAP(bmmc.atac, reduction = "lsi", dims = 2:50)
DimPlot(bmmc.atac, group.by = "cell_type", label = TRUE, repel = TRUE)

saveRDS(bmmc.atac,'bmmc.atac.rds')


# RNA-seq processing ----
count_mtx <- readRDS('bmmc.rna.counts.rds')
bmmc.rna <- CreateSeuratObject(counts = count_mtx, 
                               project = "BMMC",
                               assay = 'RNA')
metadata.rna <- readRDS('bmmc.rna.meta.rds')
bmmc.rna <- AddMetaData(bmmc.rna, metadata.rna)

bmmc.rna <- NormalizeData(bmmc.rna, normalization.method = "LogNormalize", scale.factor = 10000)
bmmc.rna <- FindVariableFeatures(bmmc.rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bmmc.rna)
bmmc.rna <- ScaleData(bmmc.rna, features = all.genes)

bmmc.rna <- RunPCA(bmmc.rna, features = VariableFeatures(object = bmmc.rna))
#bmmc.rna <- JackStraw(bmmc.rna, num.replicate = 100)
#bmmc.rna <- ScoreJackStraw(bmmc.rna, dims = 1:20)
#JackStrawPlot(bmmc.rna, dims = 1:15)
#ElbowPlot(bmmc.rna) #20 is better

bmmc.rna <- FindNeighbors(bmmc.rna, dims = 1:20)
#bmmc.rna <- FindClusters(bmmc.rna, resolution = 0.1)
bmmc.rna <- RunUMAP(bmmc.rna, dims = 1:20)
#rna <- RunTSNE(rna, reduction = 'pca', dims=1:15)
DimPlot(bmmc.rna, reduction = "umap", group.by='cell_type')

saveRDS(bmmc.rna,'bmmc.rna.rds')


# integration and label transfer ----
bmmc.rna <- readRDS('bmmc.rna.rds')
bmmc.atac <- readRDS('bmmc.atac.rds')

transfer.anchors <- FindTransferAnchors(reference = bmmc.rna, 
                                        query = bmmc.atac, 
                                        features = VariableFeatures(object = bmmc.rna), 
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = bmmc.rna$cell_type, 
                                     weight.reduction = bmmc.atac[["lsi"]])
bmmc.atac <- AddMetaData(bmmc.atac, metadata = celltype.predictions)
pred_df <- bmmc.atac@meta.data[,c('cell_type','predicted.id')]
colnames(pred_df) <- c('cell_type','prediction')
write.table(pred_df,'../results/raw/Seurat_prediction.txt',quote=F,sep=',')


predictions <- table(bmmc.atac$cell_type, bmmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions <- predictions[with(predictions,order(Var1)),]
write.csv(predictions, 'seurat_prediction.csv',quote = F)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", 
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1

mean(bmmc.atac$cell_type == bmmc.atac$predicted.id)
predictions <- read.csv('seurat_prediction.csv',stringsAsFactors = F)
acc <- predictions[predictions$Var1 == predictions$Var2,]
acc[with(acc,order(Var1)),]
write.table(acc, '../results/raw/Seurat_acc.csv', quote = F, sep=',',row.names=F)



# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(bmmc.rna)
refdata <- GetAssayData(bmmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = bmmc.atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
bmmc.atac[["RNA"]] <- imputation
coembed <- merge(x = bmmc.rna, y = bmmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "cell_type"))

umap_coembed <- coembed@reductions$umap@cell.embeddings
dim(umap_coembed)
all(gsub('_1','_2',rownames(umap_coembed)[1:20009]) == rownames(umap_coembed)[20010:40018])
write.csv(umap_coembed, '../results/raw/Seurat_UMAP.csv',quote = F)
saveRDS(coembed,'coembed.rds')

metad <- bmmc.rna@meta.data
write.table(metad, '../bmmc.rna.meta.csv',sep=',',quote=F)

coembed <- readRDS('coembed.rds')
DimPlot(coembed, group.by = c("batch", "cell_type"))

coembed <- readRDS('coembed.rds')
pca_coembed <- coembed@reductions$pca@cell.embeddings[,1:20]
dim(pca_coembed)
all(gsub('_1','_2',rownames(pca_coembed)[1:20009]) == rownames(pca_coembed)[20010:40018])
write.csv(pca_coembed, '../results/raw/Seurat_embed.pca.txt',quote = F)


