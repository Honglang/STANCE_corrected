#' @title Perform STANCE individual tests for detecting ctSVGs
#' @description
#' Perform STANCE individual tests to detect ctSVGs for each cell type of interest given.
#' @param object STANCE object
#' @param Cell_types_to_test (default NULL) a vector of character strings
#' claiming all the cell types of interest to test.
#' If NULL, tests will perform across all the cell types.
#' @param Genes_to_test (default NULL) a vector of character strings
#' claiming all the genes of interest to test.
#' If NULL, tests will perform across all the genes.
#' @param correction (default FALSE) if TRUE, perform a bias correction 
#' for the approximation of the scaled chi-square distribution.
#' @param pv.adjust (default "BY") p-value adjustment method, a character string.
#' @return return STANCE object.
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import MM4LMM
#' 
#' @export


Rcpp::cppFunction(depends = "RcppArmadillo",code = '
arma::mat invert(const arma::mat& X) {
  return arma::inv(X);
}
')

runTest2 <- function(object, 
                     Cell_types_to_test = NULL, 
                     Genes_to_test = NULL,
                     pv.adjust = "BY",
                     correction = FALSE){
  if(is.null(object@Covariance_matrices) | is.null(object@Random_effect_design_matrices) | is.null(object@Sigma_k_matrices)){
    stop('Please run \'get_modelMatrices()\' before doing tests.')
  }
  if(is.null(Cell_types_to_test)) {
    Cell_types_to_test <- object@cell_types
  }
  CT_to_test_index <- match(Cell_types_to_test, object@cell_types)
  #if(length(CT_to_test_index) == 1){
  #  if(is.na(CT_to_test_index)) {
  #    stop("No one in Cell_types_to_test matches the cell types in the object! \n")
  #  }
  #} else {
    if(all(is.na(CT_to_test_index))) {
      stop("No one in Cell_types_to_test matches the cell types in the object! \n")
    } else if(any(is.na(CT_to_test_index))) {
      CT_nomatch <- Cell_types_to_test[is.na(CT_to_test_index)]
      cat("The following cell types in Cell_types_to_test does not match any cell types in the object \n")
      cat(CT_nomatch, "\n")
      cat("The rest will be tested...\n")
      CT_to_test_index <- na.omit(CT_to_test_index)
    }
  #}
  
  if(is.null(Genes_to_test)) {
    cat("All genes will be tested...\n")
    Genes_to_test <- row.names(object@gene_expression)
  }
  Genes_to_test_index <- match(Genes_to_test, row.names(object@gene_expression))
  #if(length(Genes_to_test_index) == 1){
  #  if(is.na(Genes_to_test_index)) {
  #    stop("No one in Genes_to_test matches the genes in the object! \n")
  #  }
  #} else {
    if(all(is.na(Genes_to_test_index))) {
      stop("No one in Genes_to_test matches the genes in the object! \n")
    } else if(any(is.na(Genes_to_test_index))) {
      Genes_nomatch <- Genes_to_test[is.na(Genes_to_test_index)]
      cat("The following genes in Genes_to_test does not match any genes in the object \n")
      cat(Genes_nomatch, "\n")
      cat("The rest will be tested...\n")
      Genes_to_test_index <- na.omit(Genes_to_test_index)
    }
  #}
  
  
  n <- dim(object@gene_expression)[2]
  K <- dim(object@cell_type_compositions)[2]
  # Identity matrix
  In <- diag(rep(1,n))
  # one vector
  jn <- matrix(rep(1,n), nrow = n)
  # Design matrix
  X <- cbind(jn, object@covariates)
  Xt <- t(X)
  # Lists for covariance matrices & random effect design matrices
  Var.list <- object@Covariance_matrices
  Var.list[[K + 1]] <- In
  Z.list <- object@Random_effect_design_matrices
  Z.list[[K + 1]] <- In
  
  Test2_result <- lapply(1:length(CT_to_test_index), function(i){
    iCT <- CT_to_test_index[i]
    #if(length(Genes_to_test_index) == 1)
    # Test3_result_single <- sapply(object@gene_expression[Genes_to_test_index,], function(y))
    #if(length(Genes_to_test_index) != 1)
    Sigma_l <- object@Sigma_k_matrices[[iCT]]
    Sigma_k.list_l <- object@Sigma_k_matrices[-iCT]
    Var.list_l <- Var.list[-iCT]
    Z.list_l <- Z.list[-iCT]
    
    if(length(Genes_to_test_index) == 1){
      y <- object@gene_expression[Genes_to_test_index,]
      model.l <- MM4LMM::MMEst(Y = y,
                               Cofactor = object@covariates, 
                               VarList = Var.list_l,
                               ZList = Z.list_l,
                               Verbose = FALSE)
      VarComp.est <- ifelse(is.na(model.l$NullModel$Sigma2),
                            yes = 0,
                            no = model.l$NullModel$Sigma2)
      # V_3 matrix
      V3 <- matrix(rep(0,n*n), ncol = n)
      for (k in 1:(K-1)){
        V3 <- V3 + VarComp.est[k] * Sigma_k.list_l[[k]]
      }
      V3 <- V3 + VarComp.est[K] * In
      diag(V3) <- ifelse(diag(V3) == 0,
                         yes = 1e-10,
                         no = diag(V3))
      V3inv <- invert(V3)
      V3invX <- V3inv %*% X
      XtV3invX <- crossprod(X, V3invX)
      diag(XtV3invX) <- ifelse(diag(XtV3invX) == 0,
                               yes = 1e-10,
                               no = diag(XtV3invX))
      XtV3invX_inv <- invert(XtV3invX)
      P3 <- V3invX %*% tcrossprod(XtV3invX_inv, V3invX)
      P3 <- V3inv - P3
      
      P3Sigma_l <- P3 %*% Sigma_l
      P3Sigma_lP3 <- P3Sigma_l %*% P3
      
      e <- sum(diag(P3Sigma_l)) / 2
      
      if(correction) {
        # List for P0 %*% Sigma_k, k = 1,...,(l-1),(l+1),...K
        P3Sigma_k.list <- lapply(1:(K-1), function(k){
          P3Sigma_k <- P3 %*% Sigma_k.list_l[[k]]
        })
        
        # I_{l,l}
        I_ll <- sum(diag(P3Sigma_l %*% P3Sigma_l)) / 2
        
        # Vector I_{-l, l}
        I_nll <- sum(diag(P3Sigma_lP3)) / 2
        for(k in 1:(K-1)){
          I_nll <- c(I_nll, sum(diag(P3Sigma_l %*% P3Sigma_k.list[[k]])) / 2)
        }
        I_nll <- matrix(I_nll, nrow = K)
        
        # Matrix I_{-l, -l}
        I_nlnl <- matrix(rep(0,K*K), ncol = K)
        I_nlnl[1,1] <- sum(diag(P3 %*% P3))/2
        for (i in 2:K) {
          I_nlnl[i,1] <- sum(diag(P3Sigma_k.list[[i - 1]] %*% P3))/2
          for (j in 2:K) {
            I_nlnl[i,j] <- sum(diag(P3Sigma_k.list[[i - 1]] %*% P3Sigma_k.list[[j -1]]))/2
          }
        }
        I_nlnl[1,2:K] <- t(I_nlnl[2:K,1])
        
        # Matrix I^{-1}_{-l, -l}
        I_nlnl_inv <- invert(I_nlnl)
        
        v <- I_ll - crossprod(I_nll, I_nlnl_inv) %*% I_nll
      } else {
        v <- sum(diag(P3Sigma_l %*% P3Sigma_l)) / 2
      }
      
      # Test statistic
      # U <- (t(y) %*% P3Sigma_lP3 %*% y) / 2
      U <- (crossprod(y, P3Sigma_lP3) %*% y) / 2
      U <- as.numeric(U)
      
      # Scaler
      a <- v / (2 * e)
      # Degree of freedom
      g <- (2 * e^2) / v
      # p-value = P(U/a > chi2_g)
      p_value <- pchisq(U/a, df = g, lower.tail = F)
      
      Test2_result_single <- c(U, p_value)
      names(Test2_result_single) <- c("statistic", "p_value")
    } else {
      Test2_result_single <- apply(object@gene_expression[Genes_to_test_index,], MARGIN = 1,
                                   function(y){
                                     model.l <- MM4LMM::MMEst(Y = y,
                                                              Cofactor = object@covariates, 
                                                              VarList = Var.list_l,
                                                              ZList = Z.list_l,
                                                              Verbose = FALSE)
                                     VarComp.est <- ifelse(is.na(model.l$NullModel$Sigma2),
                                                           yes = 0,
                                                           no = model.l$NullModel$Sigma2)
                                     # V_3 matrix
                                     V3 <- matrix(rep(0,n*n), ncol = n)
                                     for (k in 1:(K-1)){
                                       V3 <- V3 + VarComp.est[k] * Sigma_k.list_l[[k]]
                                     }
                                     V3 <- V3 + VarComp.est[K] * In
                                     diag(V3) <- ifelse(diag(V3) == 0,
                                                        yes = 1e-10,
                                                        no = diag(V3))
                                     V3inv <- invert(V3)
                                     V3invX <- V3inv %*% X
                                     XtV3invX <- crossprod(X, V3invX)
                                     diag(XtV3invX) <- ifelse(diag(XtV3invX) == 0,
                                                              yes = 1e-10,
                                                              no = diag(XtV3invX))
                                     XtV3invX_inv <- invert(XtV3invX)
                                     P3 <- V3invX %*% tcrossprod(XtV3invX_inv, V3invX)
                                     P3 <- V3inv - P3
                                     
                                     P3Sigma_l <- P3 %*% Sigma_l
                                     P3Sigma_lP3 <- P3Sigma_l %*% P3
                                     
                                     e <- sum(diag(P3Sigma_l)) / 2
                                     if(correction) {
                                       # List for P0 %*% Sigma_k, k = 1,...,(l-1),(l+1),...K
                                       P3Sigma_k.list <- lapply(1:(K-1), function(k){
                                         P3Sigma_k <- P3 %*% Sigma_k.list_l[[k]]
                                       })
                                       
                                       # I_{l,l}
                                       I_ll <- sum(diag(P3Sigma_l %*% P3Sigma_l)) / 2
                                       
                                       # Vector I_{-l, l}
                                       I_nll <- sum(diag(P3Sigma_lP3)) / 2
                                       for(k in 1:(K-1)){
                                         I_nll <- c(I_nll, sum(diag(P3Sigma_l %*% P3Sigma_k.list[[k]])) / 2)
                                       }
                                       I_nll <- matrix(I_nll, nrow = K)
                                       
                                       # Matrix I_{-l, -l}
                                       I_nlnl <- matrix(rep(0,K*K), ncol = K)
                                       I_nlnl[1,1] <- sum(diag(P3 %*% P3))/2
                                       for (i in 2:K) {
                                         I_nlnl[i,1] <- sum(diag(P3Sigma_k.list[[i - 1]] %*% P3))/2
                                         for (j in 2:K) {
                                           I_nlnl[i,j] <- sum(diag(P3Sigma_k.list[[i - 1]] %*% P3Sigma_k.list[[j -1]]))/2
                                         }
                                       }
                                       I_nlnl[1,2:K] <- t(I_nlnl[2:K,1])
                                       
                                       # Matrix I^{-1}_{-l, -l}
                                       I_nlnl_inv <- invert(I_nlnl)
                                       
                                       v <- I_ll - crossprod(I_nll, I_nlnl_inv) %*% I_nll
                                     } else {
                                       v <- sum(diag(P3Sigma_l %*% P3Sigma_l)) / 2
                                     }
                                     # Test statistic
                                     # U <- (t(y) %*% P3Sigma_lP3 %*% y) / 2
                                     U <- (crossprod(y, P3Sigma_lP3) %*% y) / 2
                                     U <- as.numeric(U)
                                     
                                     # Scaler
                                     a <- v / (2 * e)
                                     # Degree of freedom
                                     g <- (2 * e^2) / v
                                     # p-value = P(U/a > chi2_g)
                                     p_value <- pchisq(U/a, df = g, lower.tail = F)
                                     
                                     output <- c(U, p_value)
                                     names(output) <- c("statistic", "p_value")
                                     return(output)
                                   })
      Test2_result_single <- t(Test2_result_single)
      Test2_result_single <- as.data.frame(Test2_result_single)
    }
    return(Test2_result_single)
  })
  # Rename the result list
  names(Test2_result) <- object@cell_types[CT_to_test_index]
  
  object@Test_2 <- Test2_result
  return(object)
}

