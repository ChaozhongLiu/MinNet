#setwd('/mnt/hdd/chaozhong/manuscript/fig2/BMMC_test/bindSC/')
#library(devtools)
#install_github('KChen-lab/bindSC')
#library(argparse)
library(bindSC)
library(Seurat)
library(Signac)


#---------------------------------
cat("[1/4] Reading data...\n")
rna <- readRDS('../seurat/bmmc.rna.rds')
atac <- readRDS('../seurat/bmmc.atac.rds')

rna.so <- CreateSeuratObject(counts = rna@assays$RNA@counts, 
                             meta.data = rna@meta.data)
atac.so <- CreateSeuratObject(counts = atac@assays$peaks@counts, 
                              meta.data = atac@meta.data)
atac2rna.so <- CreateSeuratObject(counts = atac@assays$ACTIVITY@counts, 
                                  meta.data = atac@meta.data)

rna <- FindVariableFeatures(rna,
                            selection.method = "vst",
                            nfeatures = 2500)

hvg <- VariableFeatures(rna) #rownames(rna$var)[rna$var$highly_variable]

hvg <- intersect(hvg, rownames(atac2rna.so))

X.clst <- rna@meta.data$cell_type
Y.clst <- atac@meta.data$cell_type
X.batch <- rna@meta.data$batch
Y.batch <- atac@meta.data$batch

n_cells <- ncol(rna) + ncol(atac)

rm(rna, atac)
gc()  # Reduce memory usage


cat("[2/4] Data preprocessing...\n")
rna.so <- NormalizeData(rna.so)
atac2rna.so <- NormalizeData(atac2rna.so)

X <- GetAssayData(rna.so)[hvg, ]
Z0 <- GetAssayData(atac2rna.so)[hvg, ]
out <- dimReduce(dt1 = X, dt2 = Z0, K = 30)
X <- t(out$dt1)
Z0 <- t(out$dt2)

#atac.so <- RunTFIDF(atac.so)
atac.so <- FindTopFeatures(atac.so, min.cutoff = "q0")
atac.so <- RunLSI(atac.so, n = 50)

Y <- t(Embeddings(atac.so, reduction = "lsi"))

rm(rna.so, atac.so, atac2rna.so)
gc()  # Reduce memory usage



cat("[3/4] Running bindSC...\n")
res <- BiCCA(
  X = X, Z0 = Z0, Y = Y, X.clst = X.clst, Y.clst = Y.clst,
  alpha = 0.5, lambda = 0.5, K = 15,
  temp.path = './', save = FALSE
)



cat("[4/4] Saving results...\n")
write.table(
  res$u, 'bindsc_embed.rna.txt', #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
)
write.table(
  res$r, 'bindsc_embed.atac.txt', #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
)



#UMAP plot ----
umap_plt <- umap(rbind(res$u, res$r))

umap_plt  <- data.frame("UMAP1"=umap_plt$layout[,1],
                        "UMAP2"=umap_plt$layout[,2],
                        "celltype" = c(X.clst,Y.clst),
                        "data" = c(rep("scRNA-seq",ncol(X)),
                                   rep("scATAC-seq",ncol(Y))),
                        "batch" = c(X.batch,Y.batch))

xlim <- c(min(umap_plt$UMAP1), max(umap_plt$UMAP1))
ylim <- c(min(umap_plt$UMAP2), max(umap_plt$UMAP2))

UMAP_plot(meta = umap_plt,  alpha=1,
          color = "celltype",
          xlim = xlim, ylim = ylim,
          mylabel = paletteDiscrete(umap_plt$celltype) )

UMAP_plot(meta = umap_plt,  alpha=1,
          color = "data",
          xlim = xlim, ylim = ylim,
          mylabel = paletteDiscrete(umap_plt$data) )

UMAP_plot(meta = umap_plt,  alpha=1,
          color = "batch",
          xlim = xlim, ylim = ylim,
          mylabel = paletteDiscrete(umap_plt$batch) )







