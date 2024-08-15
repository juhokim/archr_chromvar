library(ArchR)

addArchRThreads(threads = 32)

addArchRGenome("mm10")
fragpath = "/home/exouser/projects/sc_atac_seq/data/scATAC_melanoma/atac_fragments.tsv.gz"

ArrowFiles <- createArrowFiles(
  inputFiles = fragpath,
  sampleNames = 'sc_atac_mel',
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "scATAC_mel",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p <- ggAlignPlots(p1, p2, type = "h")

# Peak calling
#pathToMacs2 <- findMacs2()
pathToMacs2 <- "/home/exouser/software/anaconda3/bin/macs2"
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters", force=TRUE)
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters", pathToMacs2 = pathToMacs2)


# chromVAR 
proj <- addPeakMatrix(proj)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif", force = TRUE)

# deviation matrix
deviation_matrix <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
deviation_matrix_data <- assay(deviation_matrix)

