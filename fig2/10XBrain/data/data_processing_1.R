library(Seurat)
library(ggplot2)
library(patchwork)
library(rhdf5)

rm(list=ls())

set.seed(14)


# H5AD processing ----
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5('human_brain_3k_filtered_feature_bc_matrix.h5')

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create RNA Seurat object ----
rna <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = "RNA")
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
rna <- subset(rna, 
              subset = nCount_RNA < 25000 &
                nCount_RNA > 1000 &
                percent.mt < 10)

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
#rna <- JackStraw(rna, num.replicate = 100)
#rna <- ScoreJackStraw(rna, dims = 1:20)
#JackStrawPlot(rna, dims = 1:15)
#ElbowPlot(rna) #15 is better
rna <- FindNeighbors(rna, dims = 1:15)
rna <- FindClusters(rna, resolution = 0.1)
rna <- RunUMAP(rna, dims = 1:15)
#rna <- RunTSNE(rna, reduction = 'pca', dims=1:15)
DimPlot(rna, reduction = "umap",pt.size = 0.5)
#DimPlot(rna, reduction = "umap", group.by = 'orig.ident',pt.size = 1.5)


# Create ATAC Seurat object ----
atac <- CreateSeuratObject(counts = atac_counts, assay = "ATAC", project = "ATAC")
atac <- subset(atac, 
               subset = nCount_ATAC < 1e5 &
                 nCount_ATAC > 5e3)

dim(atac)

atac$cell.use <- colnames(atac) %in% colnames(rna)
summary(atac$cell.use)
rna$cell.use <- colnames(rna) %in% colnames(atac)
summary(rna$cell.use)


# Kepp cells consistent ----
rna <- subset(rna, subset = cell.use)
atac <- subset(atac, subset = cell.use)
all(colnames(rna) == colnames(atac))
saveRDS(rna,'rna.rds')
saveRDS(atac,'atac.rds')

# Create activity matrix ----
atac$tech <- "atac"

peaks <- atac@assays$ATAC@data
peaks <- peaks[!startsWith(rownames(peaks), 'chrUn'),]
rownames(peaks) <- gsub('^chr', '', rownames(peaks))
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "../../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)

atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
DefaultAssay(atac) <- "ACTIVITY"
atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)
rm(peaks)
rm(activity.matrix)


DefaultAssay(atac) <- "ATAC"
VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 50))
atac <- RunLSI(atac, n = 50, scale.max = NULL)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)
DimPlot(atac, reduction = "umap") #,group.by = 'predicted.id')

saveRDS(atac,'atac.rds')

# save for SiaNN ----
rna <- readRDS('rna.rds')
rna_mtx <- Matrix::t(rna[['RNA']]@counts)
rna_mtx <- as.matrix(rna_mtx)
h5createFile("rna.brain.h5")
h5write(rna_mtx, 'rna.brain.h5', 'RNA')
rna_meta <- rna@meta.data
write.table(rna_meta, 'rna.meta.csv',quote = F, sep='\t')
rna_gene <- rownames(rna)
write.table(rna_gene, 'rna_gene_name.txt',quote=F,row.names=F,col.names=F)

atac <- readRDS('atac.rds')
atac_mtx <- Matrix::t(atac[['ACTIVITY']]@counts)
atac_mtx <- as.matrix(atac_mtx)
h5createFile("atac.brain.h5")
h5write(atac_mtx, 'atac.brain.h5', 'ATAC')
atac_meta <- atac@meta.data
write.table(atac_meta, 'atac.meta.csv',quote = F, sep='\t')
DefaultAssay(atac) <- 'ACTIVITY'

atac_gene <- rownames(atac)
write.table(atac_gene, 'atac_gene_name.txt',quote=F,row.names=F,col.names=F)




















