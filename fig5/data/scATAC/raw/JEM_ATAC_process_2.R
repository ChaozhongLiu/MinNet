setwd('/mnt/hdd/chaozhong/manuscript/fig5/JEM_covid_data/scATAC/raw')
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(rhdf5)
library(pheatmap)
library(cowplot)
library(RColorBrewer)

rm(list=ls())

set.seed(1)
# ATAC processing merge ----
count_mtx <- readRDS('peaks_mtx/merge.peak.rds')
atac <- CreateSeuratObject(counts = count_mtx, 
                           project = "Covid19",
                           assay = 'peaks')
peaks <- atac@assays$peaks@counts
peaks <- peaks[startsWith(rownames(peaks), 'chr'),]
rownames(peaks) <- gsub('chr','',rownames(peaks))
dim(peaks)

peak_tmp <- peaks[,1:5000]
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peak_tmp,
                                            annotation.file = "../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                            seq.levels = c(1:22,"X","Y"), 
                                            upstream = 2000, verbose = TRUE)

for (i in 1:9) {
  peak_tmp <- peaks[,(5000*i+1):(5000*(i+1))]
  activity_tmp <- CreateGeneActivityMatrix(peak.matrix = peak_tmp,
                                              annotation.file = "../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                              seq.levels = c(1:22,"X","Y"), 
                                              upstream = 2000, verbose = TRUE)
  activity.matrix <- cbind(activity.matrix,  activity_tmp)
}

peak_tmp <- peaks[,50001:51043]
activity_tmp <- CreateGeneActivityMatrix(peak.matrix = peak_tmp,
                                         annotation.file = "../../../biolib/Homo_sapiens.GRCh38.104.gtf", 
                                         seq.levels = c(1:22,"X","Y"), 
                                         upstream = 2000, verbose = TRUE)
activity.matrix <- cbind(activity.matrix,  activity_tmp)

subindex <- sample.int(51043, 25522)
activity.matrix <- activity.matrix[,subindex]
peaks <- peaks[,subindex]
all(colnames(peaks) == colnames(activity.matrix))

VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 1000))
length(VariableFeatures(atac))
peaks_use <- VariableFeatures(atac)
dim(peaks)
rownames(peaks) <- paste0('chr',rownames(peaks))
peaks <- peaks[peaks_use,]
dim(peaks)

h5createFile(paste0("activity_mtx/",'merge',".h5"))
h5write(as.matrix(activity.matrix)[,1:5000], paste0("activity_mtx/",'merge',".h5"), 'ACTIVITY-1')
h5write(as.matrix(activity.matrix)[,5001:10000], paste0("activity_mtx/",'merge',".h5"), 'ACTIVITY-2')
h5write(as.matrix(activity.matrix)[,10001:20000], paste0("activity_mtx/",'merge',".h5"), 'ACTIVITY-3')
h5write(as.matrix(activity.matrix)[,20001:25522], paste0("activity_mtx/",'merge',".h5"), 'ACTIVITY-4')

h5write(as.matrix(peaks[,1:2500]), paste0("activity_mtx/",'merge',".h5"), 'peaks-1')
h5write(as.matrix(peaks[,2501:5000]), paste0("activity_mtx/",'merge',".h5"), 'peaks-2')
h5write(as.matrix(peaks[,5001:7500]), paste0("activity_mtx/",'merge',".h5"), 'peaks-3')
h5write(as.matrix(peaks[,7501:10000]), paste0("activity_mtx/",'merge',".h5"), 'peaks-4')
h5write(as.matrix(peaks[,10001:12500]), paste0("activity_mtx/",'merge',".h5"), 'peaks-5')
h5write(as.matrix(peaks[,12501:15000]), paste0("activity_mtx/",'merge',".h5"), 'peaks-6')
h5write(as.matrix(peaks[,15001:17500]), paste0("activity_mtx/",'merge',".h5"), 'peaks-7')
h5write(as.matrix(peaks[,17501:20000]), paste0("activity_mtx/",'merge',".h5"), 'peaks-8')
h5write(as.matrix(peaks[,20001:22500]), paste0("activity_mtx/",'merge',".h5"), 'peaks-9')
h5write(as.matrix(peaks[,22501:25522]), paste0("activity_mtx/",'merge',".h5"), 'peaks-10')

peaks_name <- rownames(peaks)
genes_name <- rownames(activity.matrix)
write.table(peaks_name, 'peaks_index.txt', quote=F, row.names = F, col.names = F)
write.table(genes_name, 'genes_index.txt', quote=F, row.names = F, col.names = F)

cell_index <- colnames(peaks)
write.table(cell_index, 'cells_index.txt', quote=F, row.names = F, col.names = F)



summary(startsWith(rownames(activity.matrix), 'ENSG0'))

