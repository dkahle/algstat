#' Compute the Lawrence lifting of a configuration matrix
#' 
#' Compute the Lawrence lifting of a configuration matrix
#' 
#' @param A the configuration matrix of the model
#' @return an integer matrix
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso Aoki, Hara, and Takemura (2012), p. 61
#' @export 
#' @examples
#' 
#' 
#' A <- hmat(c(2,2), 1:2)
#' lawrence(A)
#' 
#' 
lawrence <- function(A){
  m <- nrow(A); n <- ncol(A)
  mat <- matrix(0L, nrow = 2*m+n, ncol = 2*n)
  mat[1:m,1:n] <- A
  mat[(m+1):(2*m),(n+1):(2*n)] <- A
  mat[(2*m+1):(2*m+n),1:n] <- intDiag(n)
  mat[(2*m+1):(2*m+n),(n+1):(2*n)] <- intDiag(n)
  mat
}



intDiag <- function(int){
  mat <- as.integer(diag(int))
  dim(mat) <- c(int, int)
  mat
}