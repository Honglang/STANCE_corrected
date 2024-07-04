#' @title Estimate the variance components for STANCE model
#' @description
#' Estimate all the variance components of the full STANCE model for each gene of interest.
#' @param object STANCE object
#' @param Genes_to_est (default NULL) a vector of character strings
#' claiming all the genes of interest to estimate.
#' If NULL, estimation of variance components will perform across all the genes.
#' @param ncores (default 1) an integer value of number of CPU cores to use
#' when running variance component estimation.
#' @return return STANCE object.
#'
#' @import gaston
#'
#' @export

varcomp_est <- function(object, Genes_to_est = NULL, ncores = 1){
  n <- dim(object@gene_expression)[2]
  K <- dim(object@cell_type_compositions)[2]
  numGenes <- dim(object@gene_expression)[1]

  X <- cbind(matrix(1,nrow = n), object@covariates)
  counts <- object@gene_expression

  if(!is.null(Genes_to_est)){
    genes_to_est.list <- intersect(row.names(counts), Genes_to_est)
    if(identical(genes_to_est.list, character(0))){
      stop('None of the Genes_to_est matches the gene list!')
    } else {
      counts.use <- counts[match(genes_to_est.list, row.names(counts)),]
    }
  } else {
    counts.use <- counts
  }

  Est_result <- mclapply(1:nrow(counts.use), function(i) {
    y <- counts.use[i, ]
    ## Full model
    model.full <- gaston::lmm.aireml(Y = y, X = X,
                                     K = object@Sigma_k_matrices,
                                     verbose = FALSE)
    output <- c(model.full$tau, model.full$sigma2)
    names(output) <- c(object@cell_types, "Error_residuals")
    return(output)
  }, mc.cores = ncores)

  Est_result <- do.call(rbind, Est_result)
  rownames(Est_result) <- rownames(counts.use)

  object@VarComp_estimates <- as.data.frame(Est_result)
  return(object)
}
