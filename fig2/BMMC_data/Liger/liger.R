#setwd('/mnt/hdd/chaozhong/manuscript/fig2/BMMC_data/Liger/')
#devtools::install_github('welch-lab/liger')
library(rliger)
args <- commandArgs(trailingOnly = TRUE)
#print(args)
rand_seed <- as.integer(args[1])
set.seed(rand_seed)

#https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/walkthrough_rna_atac.html
#---------------------------------
cat("[1/4] Reading data...\n")
rna <- readRDS('../seurat/bmmc.rna.rds')
atac <- readRDS('../seurat/bmmc.atac.rds')
rna <- rna@assays$RNA@counts
atac <- atac@assays$ACTIVITY@counts

#rownames(rna$obs) <- paste(rownames(rna$obs), "RNA", sep = ".")  # Avoid collision
colnames(rna) <- paste(colnames(rna), "RNA", sep = ".")
#rownames(atac$obs) <- paste(rownames(atac$obs), "ATAC", sep = ".")  # Avoid collision
colnames(atac) <- paste(colnames(atac), "ATAC", sep = ".")

int.liger <- createLiger(list(
  atac = atac,
  rna = rna
))

rna_cells <- colnames(rna)
atac_cells <- colnames(atac)
all(rna_cells != gsub('/.ATAC','.RNA',atac_cells))
n_cells <- ncol(rna) + ncol(atac)

rm(rna, atac)
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
  union(rna_cells, atac_cells),
  rownames(combined_latent)
)  # Because of cell filtering in scaleNotCenter

combined_latent <- rbind(combined_latent, matrix(
  nrow = length(missing_cells),
  ncol = ncol(combined_latent),
  dimnames = list(missing_cells, colnames(combined_latent))
))  # Fill with NA
rna_latent <- combined_latent[rna_cells, ]
rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
atac_latent <- combined_latent[atac_cells, ]
rownames(atac_latent) <- gsub("\\.ATAC$", "", rownames(atac_latent))
write.table(
  rna_latent, paste0('Liger_embed.rna.',rand_seed,'.txt'), #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)
write.table(
  atac_latent, paste0('Liger_embed.atac.',rand_seed,'.txt'), #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)


# UMAP space
combined_latent <- int.liger@tsne.coords
missing_cells <- setdiff(
  union(rna_cells, atac_cells),
  rownames(combined_latent)
)  # Because of cell filtering in scaleNotCenter

combined_latent <- rbind(combined_latent, matrix(
  nrow = length(missing_cells),
  ncol = ncol(combined_latent),
  dimnames = list(missing_cells, colnames(combined_latent))
))  # Fill with NA
rna_latent <- combined_latent[rna_cells, ]
rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
atac_latent <- combined_latent[atac_cells, ]
rownames(atac_latent) <- gsub("\\.ATAC$", "", rownames(atac_latent))
write.table(
  rna_latent, paste0('Liger_UMAPembed.rna.',rand_seed,'.txt'), #args$output_rna,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)
write.table(
  atac_latent, paste0('Liger_UMAPembed.atac.',rand_seed,'.txt'), #args$output_atac,
  sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
)







