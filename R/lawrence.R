#' Compute the Lawrence lifting of a configuration matrix
#' 
#' Compute the Lawrence lifting of a configuration matrix
#' 
#' @param A the configuration matrix of the model
#' @param extent short or tall version
#' @return an integer matrix
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso Aoki, Hara, and Takemura (2012), pp. 61-62
#' @export 
#' @examples
#' 
#' 
#' A <- hmat(c(2,2), 1:2)
#' lawrence(A)
#' lawrence(A, extent = "tall")
#' 
#' \dontrun{ # requires 4ti2 installed
#' markov(A)
#' markov(lawrence(A))
#' markov(lawrence(A, extent = "tall"))
#' }
#' 
#' 
lawrence <- function(A, extent = c("short", "tall")){
  
  extent <- match.arg(extent)
  m <- nrow(A); n <- ncol(A)  
  
  if(extent == "short"){
    mat <- matrix(0L, nrow = m+n, ncol = 2*n)
    mat[1:m,1:n] <- A
    mat[(m+1):(m+n),1:n]         <- integerIdentityMatrix(n)
    mat[(m+1):(m+n),(n+1):(2*n)] <- integerIdentityMatrix(n)  
  } else if(extent == "tall"){
    mat <- matrix(0L, nrow = 2*m+n, ncol = 2*n)
    mat[1:m,1:n] <- A
    mat[(m+1):(2*m),(n+1):(2*n)] <- A
    mat[(2*m+1):(2*m+n),1:n]         <- integerIdentityMatrix(n)
    mat[(2*m+1):(2*m+n),(n+1):(2*n)] <- integerIdentityMatrix(n)  
  }
  
  mat
}



integerIdentityMatrix <- function(int){
  mat <- as.integer(diag(int))
  dim(mat) <- c(int, int)
  mat
}