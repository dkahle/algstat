#' algstat : Algebraic statistics in R
#' 
#' algstat is a package for algebraic statistics in R.  Current 
#' applications include exact inference in log-linear models for 
#' contingency table data, analysis of ranked and partially ranked 
#' data, and general purpose tools for multivariate polynomials, 
#' building on the mpoly package.  To aid in the process, algstat 
#' leverages ports to Macaulay2 (through m2r), Bertini, LattE and
#' 4ti2 (through latter).
#' 
#' @import stringr mpoly reshape2 Rcpp lpSolve parallel memoise 
#'   ggplot2 latter m2r 
#' @importFrom stats deriv dmultinom loglin runif sd model.frame
#' @importFrom utils combn download.file
#' @importFrom plyr ddply
#'   
#' @useDynLib algstat
#' @docType package
#' @name algstat
#' @aliases algstat package-algstat
NULL