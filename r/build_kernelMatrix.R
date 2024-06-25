#' @title Build the Gaussian kernel matrix 
#' @description
#' Build the Gaussian kernel matrix to model spatial correlation between spots.
#' The Gaussian kernel bandwidth will be selected automatically:
#' for datasets with 5000 spots or fewer, 
#' non-parametric Sheather-Jones' bandwidths are computed for each gene, 
#' and the median of these values is used as the common bandwidth; 
#' for datasets with more than 5000 spots, Silverman's bandwidths are computed for each gene, 
#' and the median of these gene-specific bandwidths is used as the common bandwidth.
#' @param object the ... object.
#' @return return ... object with the Gaussian kernel matrix.
#' @import stats
#' @import KRLS
#' @export

build_kernelMatrix <- function(object){
  if(is.null(object@original_gene_expression) | is.null(object@original_location) | is.null(object@original_cell_type_compositions)){
    stop('Please run \'data_preprocess()\'  before building the kernel matrix.')
  }
  counts <- object@gene_expression
  pos <- object@location
  n <- nrow(counts)
  
  if (n > 5000) {
    ## Bandwidth selection by Silverman's
    bw_vector <- apply(counts, MARGIN = 1, stats::bw.nrd0)
    bw <- median(na.omit(bw_vector))
    object@bandwidth <- bw
  } else {
    ## Bandwidth selection by SJ's
    bw_vector <- apply(counts, MARGIN = 1, stats::bw.SJ)
    bw <- median(na.omit(bw_vector))
    object@bandwidth <- bw
  }
  
  ## Gaussian kernel
  KK <- KRLS::gausskernel(X = pos, sigma = bw)
  object@kernel <- KK
  return(object)
}