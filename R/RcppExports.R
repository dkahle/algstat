# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

computeCRsCpp <- function(x, exp, lambda) {
    .Call('_algstat_computeCRsCpp', PACKAGE = 'algstat', x, exp, lambda)
}

computeG2sCpp <- function(x, exp) {
    .Call('_algstat_computeG2sCpp', PACKAGE = 'algstat', x, exp)
}

computeNMsCpp <- function(x, exp) {
    .Call('_algstat_computeNMsCpp', PACKAGE = 'algstat', x, exp)
}

computeUProbsCpp <- function(x) {
    .Call('_algstat_computeUProbsCpp', PACKAGE = 'algstat', x)
}

computeX2sCpp <- function(x, exp) {
    .Call('_algstat_computeX2sCpp', PACKAGE = 'algstat', x, exp)
}

metropolis_hypergeometric_cpp <- function(init, moves, iter, burn, thin) {
    .Call('_algstat_metropolis_hypergeometric_cpp', PACKAGE = 'algstat', init, moves, iter, burn, thin)
}

metropolis_uniform_cpp <- function(init, moves, iter, burn, thin) {
    .Call('_algstat_metropolis_uniform_cpp', PACKAGE = 'algstat', init, moves, iter, burn, thin)
}

rfiberOne <- function(A, b) {
    .Call('_algstat_rfiberOne', PACKAGE = 'algstat', A, b)
}

walk <- function(current, moves, iter, thin) {
    .Call('_algstat_walk', PACKAGE = 'algstat', current, moves, iter, thin)
}

