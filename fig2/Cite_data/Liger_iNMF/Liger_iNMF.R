#setwd('/mnt/hdd/chaozhong/manuscript/fig2/Cite_data/Liger_iNMF/')
#devtools::install_github('welch-lab/liger')
library(rliger)
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#rand_seed <- as.integer(args[1])
rand_seed <- 1
set.seed(rand_seed)


#---------------------------------
cat("[1/4] Reading data...\n")
rna <- readRDS('../seurat/data/cite.rna.counts.rds')
adt <- readRDS('../seurat/data/cite.adt.counts.rds')

colnames(rna) <- paste(colnames(rna), "RNA", sep = ".")
#rownames(atac$obs) <- paste(rownames(atac$obs), "ATAC", sep = ".")  # Avoid collision
colnames(adt) <- paste(colnames(adt), "ADT", sep = ".")

int.liger <- createLiger(list(
  adt = adt,
  rna = rna
))

rna_cells <- colnames(rna)
adt_cells <- colnames(adt)
all(rna_cells != gsub('/.ADT','.RNA',adt_cells))
n_cells <- ncol(rna) + ncol(adt)

rm(rna, adt)
gc()  # Reduce memory usage



#---------------------------------
cat("[2/4] Data preprocessing...\n")
int.liger <- normalize(int.liger)
int.liger <- selectGenes(int.liger, datasets.use = 2)
int.liger <- scaleNotCenter(int.liger)



#---------------------------------
cat("[3/4] Running LIGER...\n")
int.liger <- optimizeALS(int.liger, k=20, rand.seed = rand_seed)
int.liger <- quantile_norm(int.liger, rand.seed = rand_seed)
int.liger <- runUMAP(int.liger)



#---------------------------------
cat("[4/4] Saving results...\n")
combined_latent <- int.liger@H.norm
missing_cells <- setdiff(
  union(rna_cells, adt_cells),
  rownames(combined_latent)
)  # Because of cell filtering in scaleNotCenter

combined_latent <- rbind(combined_latent, matrix(
  nrow = length(missing_cells),
  ncol = ncol(combined_latent),
  dimnames = list(missing_cells, colnames(combined_latent))
))  # Fill with NA
rna_latent <- combined_latent[rna_cells, ]
rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
adt_latent <- combined_latent[adt_cells, ]
rownames(adt_latent) <- gsub("\\.ADT$", "", rownames(adt_latent))
write.table(
  rna_latent, paste0('LigerINMF_embed.rna.',rand_seed,'.txt'), #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)
write.table(
  adt_latent, paste0('LigerINMF_embed.adt.',rand_seed,'.txt'), #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)


# UMAP space
combined_latent <- int.liger@tsne.coords
missing_cells <- setdiff(
  union(rna_cells, adt_cells),
  rownames(combined_latent)
)  # Because of cell filtering in scaleNotCenter

combined_latent <- rbind(combined_latent, matrix(
  nrow = length(missing_cells),
  ncol = ncol(combined_latent),
  dimnames = list(missing_cells, colnames(combined_latent))
))  # Fill with NA
rna_latent <- combined_latent[rna_cells, ]
rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
adt_latent <- combined_latent[adt_cells, ]
rownames(adt_latent) <- gsub("\\.ADT$", "", rownames(adt_latent))
write.table(
  rna_latent, paste0('LigerINMF_UMAPembed.rna.',rand_seed,'.txt'), #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)
write.table(
  adt_latent, paste0('LigerINMF_UMAPembed.adt.',rand_seed,'.txt'), #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)












