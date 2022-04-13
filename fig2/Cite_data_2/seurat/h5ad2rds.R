library(anndata)
#library(stringr)
library(reticulate)
library(scater)
#use_python("/usr/lib/python3.6")
#remotes::install_github("mojaveazure/seurat-disk")
library(dplyr)
set.seed(42)


# RNA-seq processing ----
rna_anndat <- anndata::read_h5ad('../data/Cite_GEX.test.h5ad')
count_mtx <- rna_anndat$layers['counts']
count_mtx <- Matrix::t(count_mtx)
saveRDS(count_mtx, 'cite.rna.counts.rds')
metadata.rna <- rna_anndat$obs
saveRDS(metadata.rna, 'cite.rna.meta.rds')


# ADT processing ----
adt_anndat <- anndata::read_h5ad('../data/Cite_ADT.test.h5ad')
count_mtx <- adt_anndat$layers['counts']
count_mtx <- Matrix::t(count_mtx)
saveRDS(count_mtx, 'cite.adt.counts.rds')
metadata.rna <- adt_anndat$obs
saveRDS(metadata.rna, 'cite.adt.meta.rds')

