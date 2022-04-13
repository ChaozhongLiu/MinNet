library(SeuratData)
# install the dataset and load requirements
#InstallData("cbmc")
library(Seurat)
library(dplyr)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(rhdf5)
#install.packages('rhdf5')
rm(list=ls())


# Integration ----
pbmc.atac <- readRDS('../data/pbmc.atac.rds')
pbmc.rna <- readRDS('../data/pbmc.rna.rds')

pbmc.atac <- RenameCells(pbmc.atac, add.cell.id = 'ATAC')
pbmc.rna <- RenameCells(pbmc.rna, add.cell.id = 'RNA')

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                                     weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

anchor_set <- data.frame(transfer.anchors@anchors)
atac_true <- pbmc.rna$seurat_annotations[c(1:10412) %in% unique(anchor_set[,2])]
atac_pred <- celltype.predictions$predicted.id[c(1:10412) %in% unique(anchor_set[,2])]
mean(atac_true == atac_pred)

pbmc.atac <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                          query = pbmc.atac,
                          weight.reduction = pbmc.atac[["lsi"]], dims = 2:30,
                          store.weights=T)
weight_matrix <- pbmc.atac@tools$TransferData$weights.matrix


pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2

pred_df <- pbmc.atac@meta.data[,c('seurat_annotations','predicted.id')]
colnames(pred_df) <- c('cell_type','prediction')
write.table(pred_df,'../results/raw/Seurat_prediction.txt',quote=F,sep=',')

predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions <- predictions[with(predictions,order(Var1)),]
write.csv(predictions, '../results/raw/seurat_prediction.csv',quote = F)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", 
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) + 
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct", 
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct", 
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2



predictions <- read.csv('seurat_prediction.csv')
acc <- predictions[predictions$Var1 == predictions$Var2,]
acc[with(acc,order(Var1)),]

#co-embedding space ----
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]], 
                           dims = 2:30)
pbmc.atac[["RNA"]] <- imputation

coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

umap_coembed <- coembed@reductions$umap@cell.embeddings
dim(umap_coembed)
all(gsub('RNA','ATAC',rownames(umap_coembed)[1:10412]) == rownames(umap_coembed)[10413:20824])
write.csv(umap_coembed, '../results/raw/seurat_umap.csv',quote = F)
saveRDS(coembed,'coembed.rds')





