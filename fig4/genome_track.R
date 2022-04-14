library(SeuratData)
#install.packages('Seurat')
library(Seurat)
library(dplyr)
#BiocManager::install("Rsamtools")
#install.packages("Signac")
library(Signac)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(rhdf5)
#install.packages('rhdf5')
rm(list=ls())


atac <- readRDS('../fig2/PBMC_data/seurat/pbmc.atac.rds')
rna <- readRDS('../fig2/PBMC_data/seurat/pbmc.rna.rds')
Idents(atac) <- 'seurat_annotations'
Idents(rna) <- 'seurat_annotations'
all(colnames(rna) == colnames(atac))

atac[['RNA']] <- rna[['RNA']]

summary(atac$seurat_annotations)

#LEF1
atac <- subset(atac, idents=c('CD8 Naive','CD8 TEM_1','CD8 TEM_2',
                              'CD4 TCM', 'CD4 TEM', 'Treg', 'CD4 Naive'))

summary(atac$seurat_annotations)
atac$seurat_annotations_2 <- factor(atac$seurat_annotations,
                                    levels = c('CD8 Naive','CD4 Naive',
                                               'Treg', 'CD8 TEM_1',
                                               'CD4 TCM', 'CD4 TEM',
                                               'CD8 TEM_2'))
Idents(atac) <- 'seurat_annotations_2'
link2plot <- StringToGRanges(regions = c('chr4-108167286-108172179',
                                         'chr4-108167286-108315389'))
values(link2plot) <- DataFrame(score = c(0.8831908, 0.6128973))

Links(atac[['ATAC']]) <- link2plot

CoveragePlot(
  object = atac,
  region = "chr4-108120000-108315900",
  features = "LEF1",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = TRUE
)
# zoom in
CoveragePlot(
  object = atac,
  region = "chr4-108200000-108315900",
  features = "LEF1",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = TRUE
)

#CCR2
atac$cell.type <- 'none'
atac$cell.type[atac$seurat_annotations %in% c('CD8 Naive','CD8 TEM_1',
                                              'CD8 TEM_2')] <- 'CD8 T'
atac$cell.type[atac$seurat_annotations %in% c('CD4 TCM','CD4 TEM',
                                              'CD4 Naive')] <- 'CD4 T'
atac$cell.type[atac$seurat_annotations %in% c('Naive B','Intermediate B',
                                              'Memory B')] <- 'B cell'
atac$cell.type[atac$seurat_annotations %in% c('Naive B','Intermediate B',
                                              'Memory B')] <- 'B cell'
atac$cell.type[atac$seurat_annotations %in% c('gdT')] <- 'gdT'
atac$cell.type[atac$seurat_annotations %in% c('MAIT')] <- 'MAIT'
atac$cell.type[atac$seurat_annotations %in% c('pDC')] <- 'pDC'
atac$cell.type[atac$seurat_annotations %in% c('cDC')] <- 'cDC'
atac$cell.type[atac$seurat_annotations %in% c('CD14 Mono')] <- 'CD14 Mono'
atac$cell.type[atac$seurat_annotations %in% c('CD16 Mono')] <- 'CD16 Mono'
atac$cell.type <- factor(atac$cell.type)
Idents(atac) <- 'cell.type'

#Idents(atac) <- 'seurat_annotations'

atac <- subset(atac, idents=c('CD8 T','CD4 T','B cell','CD14 Mono',
                              'CD16 Mono', 'cDC', 'gdT', 'MAIT', 'pDC'))

link2plot <- StringToGRanges(regions = c('chr3-46208488-46354308',
                                         'chr3-46318335-46354308',
                                         'chr3-46299654-46354308',
                                         'chr3-46213035-46354308',
                                         'chr3-46312980-46354308',
                                         'chr3-46228635-46354308'))

values(link2plot) <- DataFrame(score = c(0.647090439, 0.634225218,0.612477135,
                                         0.611680250, 0.607241642, 0.605348218))

Links(atac[['ATAC']]) <- link2plot

CoveragePlot(
  object = atac,
  region = "chr3-46205000-46364308",
  features = "CCR2",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = TRUE
)

summary(atac$seurat_annotations)




