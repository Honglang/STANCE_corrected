#' @title Construct the model matrices
#' @description
#' Construct the important matrices, including design matrices, covariance matrices 
#' and Sigma_k matrices for all the random effects, before fitting STANCE model.
#' @param object the STANCE object.
#' @return return STANCE object.
#' 
#' @export

get_modelMatrices <- function(object){
  if(is.null(object@kernel)){
    stop('Please run \'build_kernelMatrix()\' first.')
  }
  n <- dim(object@gene_expression)[2]
  K <- dim(object@cell_type_compositions)[2]
  
  ## List for the design matrices of all (K+1) random effects
  # For k-th spatial variance components, the design matrix should be
  # the diagonal matrix with k-th column of cell type composition matrix
  # as diagonal elements
  REDesign.list <- lapply(1:K, function(k, obj = object){
    REDesign <- diag(obj@cell_type_compositions[,k])
    return(REDesign)
  })
  # For the error term, the design matrix should be an identity matrix
  #REDesign.list[[K + 1]] <- I_n
  object@Random_effect_design_matrices <- REDesign.list
  
  ## List for covariance matrices of all (K+1) random effects
  # For each of K spatial variance components, the covariance matrix should be
  # the kernel matrix
  Cov.list <- rep(list(object@kernel), K)
  # For the error term, the covariance matrix should be an identity matrix
  #Cov.list[[K + 1]] <- I_n
  object@Covariance_matrices <- Cov.list
  
  # List for Sigma_k & Sigma matrices
  Sigma_k.list <- lapply(1:K, function(k, obj = object){
    Sigma_k <- obj@Random_effect_design_matrices[[k]] %*% 
      Cov.list[[k]] %*% 
      obj@Random_effect_design_matrices[[k]]
    return(Sigma_k)
  })
  object@Sigma_k_matrices <- Sigma_k.list
  
  return(object)
}