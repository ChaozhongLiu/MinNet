setwd('/mnt/hdd/chaozhong/manuscript/fig2/Cite_data_2/bindSC/')
#library(devtools)
#install_github('KChen-lab/bindSC')
#library(argparse)
library(bindSC)
library(Seurat)
library(Signac)

#https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/CITE-seq/CITE_seq.html
#---------------------------------
cat("[1/4] Reading data...\n")
rna <- readRDS('../seurat/cite.rna.rds')
adt <- readRDS('../seurat/cite.adt.rds')
#rna.homo <- readRDS('../seurat/cite.rna.homo.rds')
gene.use <- intersect(rownames(rna),rownames(adt))

rna.so <- CreateSeuratObject(counts = rna@assays$RNA@counts, 
                             meta.data = rna@meta.data)
adt.so <- CreateSeuratObject(counts = adt@assays$ADT@counts[gene.use,], 
                             meta.data = adt@meta.data)
rna2adt.so <- CreateSeuratObject(counts=rna.so@assays$RNA@counts[gene.use,])

rna <- FindVariableFeatures(rna,
                            selection.method = "vst",
                            nfeatures = 3000)

hvg <- VariableFeatures(rna) #rownames(rna$var)[rna$var$highly_variable]

X.clst <- adt@meta.data$cell_type
Y.clst <- rna@meta.data$cell_type
X.batch <- adt@meta.data$batch
Y.batch <- rna@meta.data$batch

n_cells <- ncol(rna) + ncol(adt)

rm(rna, adt)
gc()  # Reduce memory usage


cat("[2/4] Data preprocessing...\n")
rna.so <- NormalizeData(rna.so)
adt.so <- NormalizeData(adt.so)
rna2adt.so <- NormalizeData(rna2adt.so)
all(rownames(rna2adt.so) == rownames(adt.so))

X <- GetAssayData(adt.so)
Z0 <- GetAssayData(rna2adt.so)
Y <- GetAssayData(rna.so)[hvg,]

rm(rna.so, adt.so, rna2adt.so)
gc()  # Reduce memory usage



cat("[3/4] Running bindSC...\n")
res <- BiCCA( X = X ,
              Y = Y, 
              Z0 = Z0, 
              X.clst = X.clst,
              Y.clst = Y.clst,
              alpha = 0.1, 
              lambda = 0.7,
              K = 15,
              temp.path  = "./",
              num.iteration = 50,
              tolerance = 0.01,
              save = FALSE,
              parameter.optimize = FALSE,
              block.size = 0)


cat("[4/4] Saving results...\n")
write.table(
  res$u, 'bindsc_embed.adt.txt', #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
)
write.table(
  res$r, 'bindsc_embed.rna.txt', #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
)



#UMAP plot ----
umap_plt <- umap(rbind(res$u, res$r))

umap_plt  <- data.frame("UMAP1"=umap_plt$layout[,1],
                        "UMAP2"=umap_plt$layout[,2],
                        "celltype" = c(X.clst,Y.clst),
                        "data" = c(rep("ADT",ncol(X)),
                                   rep("GEX",ncol(Y))),
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

bindsc_umap <- umap_plt[,c('UMAP1','UMAP2')]
write.table(
  bindsc_umap, 'bindsc_UMAP.txt', #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
)










