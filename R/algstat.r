#' algstat : Algebraic statistics in R
#'
#' algstat is a package for algebraic statistics in R.  Current applications
#' include exact inference in log-linear models for contingency table data,
#' analysis of ranked and partially ranked data, and general purpose tools for
#' multivariate polynomials, building on the mpoly package.  To aid in the
#' process, algstat leverages ports to Macaulay2 (through m2r), Bertini (through
#' bertini), LattE and 4ti2 (through latte).
#'
#' @import Rcpp mpoly latte m2r bertini
#' @importFrom stringr str_detect str_c str_dup str_replace str_replace_all
#'   str_split str_sub str_sub<-
#' @importFrom lpSolve lp
#' @importFrom ggplot2 ggplot scale_x_continuous scale_y_continuous
#'   scale_fill_gradient scale_fill_gradient2 qplot theme_bw coord_equal theme
#'   element_blank
#' @importFrom parallel mclapply makeCluster stopCluster parLapply detectCores
#' @importFrom reshape2 melt
#' @importFrom stats deriv dmultinom loglin runif sd
#' @importFrom utils combn download.file
#' @useDynLib algstat
#' @docType package
#' @name algstat
#' @aliases algstat package-algstat
NULL