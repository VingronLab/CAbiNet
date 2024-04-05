
#' @import APL
#' @import RcppEigen
#' @import methods
#' @importFrom Rcpp evalCpp
#' @importMethodsFrom Matrix colSums
#' @importFrom Matrix sparseMatrix
#' @importFrom leiden leiden
#' @importFrom dplyr desc arrange
#' @importFrom ggplot2 ggplot aes geom_point theme_bw
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData rowData
#' @importClassesFrom biclust Biclust
#' @importFrom BiocParallel SerialParam
NULL


# Load python packages
.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
}
