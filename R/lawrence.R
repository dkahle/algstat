#' Compute the Lawrence lifting of a configuration matrix
#'
#' Compute the Lawrence lifting of a configuration matrix
#'
#' @param A the configuration matrix of the model
#' @param k kth Lawrence lifting
#' @param extent short or tall version
#' @return an integer matrix
#' @author David Kahle \email{david@@kahle.io}
#' @seealso Aoki, Hara, and Takemura (2012), pp. 61-62
#' @export
#' @examples
#'
#'
#' (A <- hmat(c(2,2), 1:2))
#' plot_matrix(A) # in package latte
#'
#' lawrence(A)
#' plot_matrix(lawrence(A))
#'
#' lawrence(A, 2)
#' plot_matrix(lawrence(A, 2))
#'
#' lawrence(A, extent = "tall")
#' plot_matrix(lawrence(A, extent = "tall"))
#'
#' plot_matrix(A)
#' plot_matrix(lawrence(A))
#' plot_matrix(lawrence(A, 2))
#' plot_matrix(lawrence(A, 3))
#' plot_matrix(lawrence(A, 20))
#'
#' plot_matrix(lawrence(A))
#' plot_matrix(lawrence(A, extent = "tall"))
#'
#'
#' \dontrun{ # requires 4ti2
#'
#' J <- 10 # number of levels
#' A <- rbind(1L, 1:J)
#' markov(A)
#'
#' lawrence(A)
#' plot_matrix(lawrence(A))
#' markov(lawrence(A))
#' # markov(lawrence(A, 2)) # takes unknown long time when J=10
#' zbasis(lawrence(A))
#'
#' }
#'
#' 
lawrence <- function(A, k = 1, extent = c("short", "tall")){
  
  extent <- match.arg(extent)
  if(k > 1 && extent == "tall") stop("if extent = \"tall\", k must be 1.")
  m <- nrow(A); n <- ncol(A)  
  
  if(k == 1 && extent == "short"){
    mat <- matrix(0L, nrow = m+n, ncol = 2*n)
    mat[1:m,1:n] <- A
    mat[(m+1):(m+n),1:n]         <- integer_identity_mat(n)
    mat[(m+1):(m+n),(n+1):(2*n)] <- integer_identity_mat(n)  
  } else if(k == 1 && extent == "tall"){
    mat <- matrix(0L, nrow = 2*m+n, ncol = 2*n)
    mat[1:m,1:n] <- A
    mat[(m+1):(2*m),(n+1):(2*n)] <- A
    mat[(2*m+1):(2*m+n),1:n]         <- integer_identity_mat(n)
    mat[(2*m+1):(2*m+n),(n+1):(2*n)] <- integer_identity_mat(n)  
  }
  
  
  if(k > 1 && extent == "short"){
    mat <- matrix(0L, nrow = k*m + n, ncol = (k+1)*n)
    for(i in 0:(k-1)) mat[i*m + 1:m, i*n + 1:n] <- A
    for(i in 0:k) mat[k*m + 1:n, i*n + 1:n] <- integer_identity_mat(n)
  } 
  
  mat
}



integer_identity_mat <- function(int){
  structure(
    as.integer(diag(int)),
    dim = c(int, int)
  )
}