devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#args <- commandArgs(trailingOnly = TRUE)
#sample_id <- as.integer(args[1])
library(ArchR)
library(pheatmap)
library(Signac)

#ArchR::installExtraPackages()
set.seed(1)
addArchRThreads(threads = 8)
addArchRGenome("hg38")


archr.proj <- loadArchRProject("data/scATAC/raw/merge/")
summary(factor(archr.proj$Sample))
severity <- rep('none',length(archr.proj$Sample))
severity[archr.proj$Sample %in% c('HIP002','HIP015','HIP023',
                                  'HIP043','HIP044','HIP045')] <- '0'
severity[archr.proj$Sample %in% c('555_1','556','558')] <- '4-5'
severity[archr.proj$Sample %in% c('555_2','559')] <- '6-7'
archr.proj$severity <- severity

atac_meta <- read.csv('data/JEM_ATAC_meta.consistent.csv')
summary(archr.proj$cellNames %in% gsub('^ATAC_','',atac_meta$X))
idxSample <- BiocGenerics::which(archr.proj$cellNames %in% gsub('^ATAC_','',atac_meta$X))
cellsSample <- archr.proj$cellNames[idxSample]
archr.proj.use <- archr.proj[cellsSample, ]
archr.proj.use <- archr.proj[gsub('^ATAC_','',atac_meta$X), ]
summary(archr.proj.use$cellNames == gsub('^ATAC_','',atac_meta$X))

# NK
archr.proj.use$CellType <- atac_meta$cell_type
idxSample <- BiocGenerics::which(archr.proj.use$CellType %in% c('NK'))
cellsSample <- archr.proj.use$cellNames[idxSample]
archr.proj.use.ct <- archr.proj.use[cellsSample, ]

feature2plot <- StringToGRanges(regions = c('chr6-32940846-32941346'))

p <- plotBrowserTrack(
  ArchRProj = archr.proj.use.ct, 
  groupBy = "severity", 
  geneSymbol = c('HLA-DPB1'),
  features = c(feature2plot),
  upstream = 150000,
  downstream = 20000
)
  
grid::grid.newpage()
grid::grid.draw(p$`HLA-DPB1`)


## Monocytes
idxSample <- BiocGenerics::which(archr.proj.use$CellType %in% c('CD14 Mono','CD16 Mono'))
cellsSample <- archr.proj.use$cellNames[idxSample]
archr.proj.use.ct <- archr.proj.use[cellsSample, ]

feature2plot <- StringToGRanges(regions = c('chr6-31541012-31541512'))

p <- plotBrowserTrack(
  ArchRProj = archr.proj.use.ct, 
  groupBy = "severity", 
  geneSymbol = c('LST1'),
  features = c(feature2plot),
  upstream = 50000,
  downstream = 20000
)

grid::grid.newpage()
#5x10
grid::grid.draw(p$`LST1`)


