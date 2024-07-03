#' @title Build the Gaussian kernel matrix and compute covariance matrices
#' @description
#' Build the Gaussian kernel matrix to model spatial correlation between spots, and
#' then contruct covariance matrices (Sigma_k) for all the variance components.
#' The Gaussian kernel bandwidth will be selected automatically:
#' for datasets with 5000 spots or fewer,
#' non-parametric Sheather-Jones' bandwidths are computed for each gene,
#' and the median of these values is used as the common bandwidth;
#' for datasets with more than 5000 spots, Silverman's bandwidths are computed for each gene,
#' and the median of these gene-specific bandwidths is used as the common bandwidth.
#' @param object the STANCE object.
#' @return return STANCE object with the Gaussian kernel matrix.
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

  # Number of cell types
  K <- dim(object@cell_type_compositions)[2]

  # List for Sigma_k & Sigma matrices
  Sigma_k.list <- lapply(1:K, function(k, obj = object){
    REDesign <- diag(obj@cell_type_compositions[,k])
    Sigma_k <- REDesign %*% object@kernel %*% REDesign
    return(Sigma_k)
  })
  object@Sigma_k_matrices <- Sigma_k.list

  return(object)
}
