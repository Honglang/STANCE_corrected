#' @useDynLib STANCE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Invert a matrix using RcppArmadillo
#'
#' @param mat A numeric matrix to be inverted
#' @return The inverted matrix
#' @export
#' @examples
#' invert(matrix(c(1, 2, 3, 4), nrow = 2))
invert <- function(mat) {
  stopifnot(is.matrix(mat))
  .Call('_STANCE_invert', mat)
}
