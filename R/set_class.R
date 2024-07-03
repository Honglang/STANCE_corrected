#' @title An S4 class to store key information.
#'
#' @slot gene_expression The gene expression matrix with genes in rows
#' and spots in columns.
#' @slot location The location coordinate matrix of spots.
#' @slot cell_type_compositions The cell type composition matrix
#' to store cell type proportion at each spot.
#' Rows represent spots and columns represent cell types.
#' @slot covariates The covariate design matrix modeling gene expression.
#' @slot cell_types Cell types included in this experiment.
#' @slot kernel The kernel matrix to capture spatial correlation.
#' The default is a Gaussian kernel matrix.
#' @slot bandwidth The bandwidth parameter of Gaussian kernel.
#' @slot Sigma_k_matrices A list for Sigma_k matrices.
#' @slot Test_1 The results of STANCE Test 1 (overall test for SVGs and ctSVGs),
#' containing the test statistic values, p values and possibly adjusted p values
#' for each gene.
#' @slot Test_2 The results of STANCE Test 2 (test for ctSVGs of each cell type of interest),
#' containing the test statistic values, p values and possibly adjusted p values
#' for each gene and each cell type of interest.
#' @slot VarComp_estimates The ReML estimates of all the variance components for each gene.
#' @slot cell_type_top_genes The top significant genes with highest proportion of variance
#' explained by the spatial effect corresponding to the cell type of interest.
#' @slot original_gene_expression Used only after the normalization of expression count data
#' and the quality control steps,
#' so as to store the original gene expression count matrix.
#' @slot original_location Used only after the scaling procedure of location coordinates
#' and the quality control steps,
#' so as to store the original spot location matrix.
#' @slot original_cell_type_compositions Used only the quality control steps,
#' so as to store the original spot location matrix.
#' @import methods
#' @export

setClass("STANCE", slots = list(
  gene_expression = "matrix",
  location = "matrix",
  cell_type_compositions = "matrix",
  covariates = "ANY",
  cell_types = "character",
  kernel = "ANY",
  bandwidth = "numeric",
  Sigma_k_matrices = "ANY",
  Test_1 = "ANY",
  Test_2 = "ANY",
  VarComp_estimates = "ANY",
  cell_type_top_genes = "ANY",
  original_gene_expression = "matrix",
  original_location =  "matrix",
  original_cell_type_compositions =  "matrix"
))
