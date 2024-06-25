#' @title Preprocess data
#' @description
#' Preprocess the gene expression, location coordinate and cell type proportion matrices,
#' including normalization of count data (use NormalizaData in Seurat package),
#' scaling of coodinates and quality control at both spot and gene levels.
#' @param object STANCE object
#' @param spot.threshold (default 10) filter out the spots whose total expression count 
#' lower than the threshold.
#' @param gene.threshold (default 0.05) filter out low-expressed genes
#'whose zero-expression rate lower than the threshold.
#' @return return STANCE object.
#' 
#' @import Seurat
#' 
#' @export

data_preprocess <- function(object, spot.threshold = 10, gene.threshold = 0.05){
  counts <- object@gene_expression
  prop <- object@cell_type_compositions
  pos <- object@location
  
  ## Normalize count data
  cat('Normalizing the count data...\n')
  counts.normalized <- Seurat::NormalizeData(counts)
  
  ## Scale locations
  pos.use <- pos
  pos.use[,1] <- (pos[,1] - min(pos)) / (max(pos) - min(pos))
  pos.use[,2] <- (pos[,2] - min(pos)) / (max(pos) - min(pos))
  cat('The location matrix has been scaled.\n')
  
  ## Spot QC
  spots.use <- which(colSums(counts) >= spot.threshold)
  counts.use <- counts.normalized[, spots.use]
  pos.use <- pos.use[spots.use, ]
  prop.use <- prop[spots.use,]
  numSpots.removed <- nrow(pos) - nrow(pos.use)
  cat(paste(numSpots.removed, 'spots with gene expression less than', spot.threshold, 'have been removed. \n'))
  
  ## Gene QC
  # Filter out mt genes
  if (length(grep("mt-", rownames(counts.use), ignore.case = T)) > 0) {
    mt_gene_list <- grep("mt-", rownames(counts.use), ignore.case = T)
    counts.use <- counts.use[-mt_gene_list,]
    cat(paste(length(mt_gene_list), 'mt genes have been removed. \n'))
  }
  # Filter out low expressed genes
  ExpR <- apply(counts.use, MARGIN = 1, FUN = function(data.vector){
    non.zero <- sum(data.vector != 0)
    return(non.zero / length(data.vector))
  })
  gene.use <- intersect(names(ExpR[ExpR >= gene.threshold]), rownames(counts.use))
  counts.use <- counts.use[match(gene.use, rownames(counts.use)),]
  numGenes.removed <- nrow(counts) - nrow(counts.use)
  cat(paste(numGenes.removed, 'genes with gene zero expression rate less than', gene.threshold, 'have been removed. \n'))
  
  object@gene_expression <- as.matrix(counts.use)
  object@location <- pos.use
  object@cell_type_compositions <- prop.use
  object@original_location <- pos
  object@original_gene_expression <- counts
  object@original_cell_type_compositions <- prop
  return(object)
}