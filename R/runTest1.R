#' @title Perform STANCE overall tests for detecting SVGs and ctSVGs
#' @param object the STANCE object.
#' @param correction (default FALSE) if TRUE, perform a bias correction
#' for the approximation of the scaled chi-square distribution.
#' @param pv.adjust (default "BY") p-value adjustment method, a character string.
#' @return the STANCE object with overall test results.
#' @export

runTest1 <- function(object, correction = FALSE, pv.adjust = "BY"){
  if(is.null(object@kernel) | is.null(object@Sigma_k_matrices)){
    stop('Please run \'build_kernelMatrix()\' before doing tests.')
  }
  n <- dim(object@gene_expression)[2]
  K <- dim(object@cell_type_compositions)[2]
  numGenes <- dim(object@gene_expression)[1]

  # one vector
  jn <- matrix(rep(1,n), nrow = n)

  # Design matrix
  X <- cbind(jn, object@covariates)
  Xt <- t(X)
  XtX <- crossprod(X,X)
  XtX_inv <- solve(XtX)
  X_XtX_inv_Xt <- X %*% XtX_inv %*% Xt
  # P_tilde matrix
  P0tilde <- diag(n) - X_XtX_inv_Xt

  # Sigma matrix
  Sigma <- matrix(rep(0, n*n), nrow = n)
  for (k in 1:K) {
    Sigma <- Sigma + object@Sigma_k_matrices[[k]]
  }

  # PSigmaP matrix
  P0tildeSigma <- P0tilde %*% Sigma
  P0tildeSigmaP0tilde <- P0tildeSigma %*% P0tilde

  cat('Start running STANCE Stage 1 test...\n')

  # Test 1 statistic & the estimates of beta and sigma2
  Test1_result <- apply(object@gene_expression, MARGIN = 1,
                        function(y, Xm = X){
                          model.null <- lm(y ~ X - 1)
                          beta <- model.null$coefficients
                          sigma2_epsilon <- mean(model.null$residuals^2)
                          Utilde <- crossprod(y, P0tildeSigmaP0tilde) %*% y
                          U <- as.numeric(Utilde) / (2 * sigma2_epsilon^2)
                          output <- c(U, beta, sigma2_epsilon)
                          names(output) <- c("U", paste0("beta", 0:(length(beta)-1)), "sigma2_epsilon")
                          return(output)
                        })

  # Common part of E(U)
  e.share <- sum(diag(P0tildeSigma))
  if(correction == F){
    # Common part of Var(U)
    v.share <- sum(diag(P0tildeSigma %*% P0tildeSigma))
  } else{
    # List for P0 %*% Sigma_k
    P0tildeSigma_k.list <- lapply(1:K, function(k){
      P0tildeSigma_k <- P0tilde %*% object@Sigma_k_matrices[[k]]
    })

    # Information matrix for tau-tau
    # (i,j)-th element = tr(P0 %*% Sigma_i %*% P0 %*% Sigma_j)
    Itt <- matrix(rep(0,K*K), ncol = K)
    for (i in 1:K) {
      for (j in 1:K) {
        Itt[i,j] <- sum(diag(P0tildeSigma_k.list[[i]] %*% P0tildeSigma_k.list[[j]]))
      }
    }

    # Information for tau-sigma
    # k-th element = tr(P0 %*% Sigma_k)
    Its <- c()
    for (k in 1:K) {
      Its[k] <- sum(diag(P0tildeSigma_k.list[[k]]))
    }
    Its <- matrix(Its, nrow = K)

    # Information for sigma-sigma = tr(P_o)
    Iss <- sum(diag(P0tilde))
    Iss <- as.numeric(Iss)

    # Efficient information for tau-tau
    # = Itt - Its %*% Iss %*% t(Its)
    IEtt <- Itt - (tcrossprod(Its, Its) / Iss)

    # Common part of Var(U) = sum of element of IEtt
    v.share <- sum(IEtt)
  }

  # Compute p-value
  p_values <- apply(Test1_result, MARGIN = 2,
                    function(input, e.tilde = e.share, v.tilde = v.share){
                      test_statistic <- input["U"]
                      sigma2_epsilon <- input["sigma2_epsilon"]
                      # E(U)
                      e <- e.tilde / (2 * sigma2_epsilon)
                      # Var(U)
                      v <- v.tilde / (2 * sigma2_epsilon^2)
                      # Scaler
                      a <- v / (2 * e)
                      # Degree of freedom
                      g <- (2 * e^2) / v
                      # p-value = P(U/a > chi2_g)
                      p_value <- pchisq(test_statistic/a,
                                        df = g, lower.tail = F)
                      return(p_value)
                    })

  # FDR control
  p_values.adj <- p.adjust(p_values, method = pv.adjust,
                           n = length(p_values))
  output <- data.frame(statistic = Test1_result["U",],
                       p_value = p_values,
                       p_value_adj = p_values.adj)
  row.names(output) <- row.names(object@gene_expression)
  object@Test_1 <- output
  return(object)
}
