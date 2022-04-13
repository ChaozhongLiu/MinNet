library(anndata)
library(reticulate)
library(dplyr)
set.seed(42)
use_python("/usr/lib/python3.6")


# RNA-seq processing ----
rna_anndat <- anndata::read_h5ad('../data/Multiome_GEX.test.h5ad')
count_mtx <- rna_anndat$layers['counts']
count_mtx <- Matrix::t(count_mtx)
saveRDS(count_mtx, 'bmmc.rna.counts.rds')
metadata.rna <- rna_anndat$obs
saveRDS(metadata.rna, 'bmmc.rna.meta.rds')
# meta data needed by other algorithms
write.table(metadata.rna, '../bmmc.rna.meta.csv', quote=F, sep=',')


# ATAC-seq processing ----
atac_anndat <- anndata::read_h5ad('../data/Multiome_ATAC.test.h5ad')
count_mtx <- atac_anndat$layers['counts']
count_mtx <- Matrix::t(count_mtx)
saveRDS(count_mtx, 'bmmc.atac.counts.rds')
metadata.atac <- atac_anndat$obs
saveRDS(metadata.atac, 'bmmc.atac.meta.rds')


