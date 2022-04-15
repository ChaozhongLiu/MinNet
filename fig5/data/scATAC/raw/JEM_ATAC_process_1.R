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


#################################################
# Fragments file downloaded from GSE174072
#################################################

# create arrow file ----
samples <- c('555_1','555_2','556','558','559',
             'HIP002','HIP015','HIP023','HIP043','HIP044','HIP045')
for (sample_id in samples) {
  meta <- read.csv('../sample_list.txt')
  file_name <- meta[sample_id, 'filename']
  sample_name <- meta[sample_id, 'orig.ident']
  TSS_cut <- as.numeric(meta[sample_id, 'TSS'])
  frag_cut <- as.numeric(meta[sample_id, 'fragments'])
  
  
  ArrowFiles <- createArrowFiles(
    inputFiles = paste0('fragments/',file_name),
    sampleNames = sample_name,
    minTSS = TSS_cut, #Dont set this too high because you can always increase later
    minFrags = 10^frag_cut, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE)
  
  #predicting doublets
  doubScores <- addDoubletScores(
    input = paste0(sample_name,'.arrow'),
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
  )
}



#create archR project ----
archr.proj <- ArchRProject(
  ArrowFiles = c(paste0(samples,'.arrow')), 
  outputDirectory = 'merge',
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
paste0("Memory Size = ", round(object.size(archr.proj) / 10^6, 3), " MB")
getAvailableMatrices(archr.proj)

#filter doublets
archr.proj <- filterDoublets(archr.proj, cutEnrich = 2)

#CLustering
#LSI dimensional reduction
archr.proj <- addIterativeLSI(
  ArchRProj = archr.proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30
)

#Harmony remove batch effect
archr.proj <- addHarmony(
  ArchRProj = archr.proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#find cluster
archr.proj <- addClusters(
  input = archr.proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustersHarmony",
  resolution = 0.8
)


#pseudo-bulk replicates ----
archr.proj <- addGroupCoverages(ArchRProj = archr.proj, groupBy = "ClustersHarmony")

#Peak calling ----
pathToMacs2 <- findMacs2()
summary(factor(archr.proj$ClustersHarmony))
archr.proj <- addReproduciblePeakSet(
  ArchRProj = archr.proj, 
  groupBy = "ClustersHarmony", 
  pathToMacs2 = pathToMacs2
)

archr.proj <- addPeakMatrix(archr.proj)
getAvailableMatrices(archr.proj)

saveArchRProject(ArchRProj = archr.proj, outputDirectory = 'merge', load = FALSE)

peak_mtx <- getMatrixFromProject(archr.proj, useMatrix='PeakMatrix')
peakset <- paste0(seqnames(peak_mtx@rowRanges), ':',
                  start(peak_mtx@rowRanges), '-',
                  end(peak_mtx@rowRanges))
peak_mtx <- peak_mtx@assays@data$PeakMatrix
#print(paste0(samples[i],' has ',dim(peak_mtx)[2],' cells'))
rownames(peak_mtx) <- peakset
saveRDS(peak_mtx, paste0('peaks_mtx/','merge.peak.rds'))
