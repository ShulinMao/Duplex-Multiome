isDoublet <- function(X, ident_rate = 0.05){
  # https://www.bioconductor.org/packages/release/bioc/html/scds.html
  require(scds)
  # Getting metadata from the Seurat Object
  nCells <- ncol(X)
  barcodeID <- colnames(X)
  # Creating a SingleCellExperiment Object
  X <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = X@assays$RNA@counts))
  # Scoring Cells
  X <- cxds(X)
  X <- X$cxds_score
  # Defining the return vector
  isDoublet <- rep(FALSE, nCells)
  names(isDoublet) <- barcodeID
  
  isDoublet[X > quantile(X, probs = (1 - ident_rate))] <- TRUE
  #isDoublet[names(boxplot.stats(X)$out)] <- TRUE
  # Return
  return(isDoublet)
}

scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  #require(Matrix)
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)
  
  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mtThreshold) & (librarySize < 2 * mean(librarySize)))
  } else {
    selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 3 * mean(librarySize)))
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}


