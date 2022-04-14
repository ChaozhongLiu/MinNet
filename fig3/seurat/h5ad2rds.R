library(anndata)
#library(stringr)
library(reticulate)
library(Seurat)
library(scater)
use_python("/usr/lib/python3.6")
#remotes::install_github("mojaveazure/seurat-disk")
library(dplyr)



# RNA-seq processing ----
batchs <- c('s1d2','s3d7','s4d1')
for (i in 1:3) {
  rna_anndat <- anndata::read_h5ad(paste0('../data/Multiome_GEX.',batchs[i],'.h5ad'))
  count_mtx <- rna_anndat$layers['counts']
  count_mtx <- Matrix::t(count_mtx)
  saveRDS(count_mtx, paste0('../data/multiome.',batchs[i],'.rna.counts.rds'))
  metadata.rna <- rna_anndat$obs
  saveRDS(metadata.rna, paste0('../multiome.',batchs[i],'.meta.rds'))
  
  atac_anndat <- anndata::read_h5ad(paste0('../data/Multiome_ATAC.peak.',batchs[i],'.h5ad'))
  count_mtx <- atac_anndat$layers['counts']
  count_mtx <- Matrix::t(count_mtx)
  saveRDS(count_mtx, paste0('../data/multiome.',batchs[i],'.atac.counts.rds'))
}

