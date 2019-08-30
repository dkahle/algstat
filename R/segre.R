#' Segre Product
#' 
#' Compute the Segre product of an arbitrary number of matrices
#'
#' @param ... A listing of matrices
#'
#' @return A matrix that is the Segre product of the specified matrices.
#' @export
#'
#' @examples
#'
#' A <- B <- C <- matrix(c(1,1,1,2,1,3,1,4,1,5), nrow = 2, ncol = 5)
#' 
#' # two matrices
#' segre(A, B)
#' 
#' # more 
#' segre(A, B, C)
segre <- function(...) {
  Reduce(function(x,y) {
    dim_x <- dim(x)
    dim_y <- dim(y)
    cols <- dim_x[2] * dim_y[2]
    rows <- dim_x[1] + dim_y[1]
    mat <- matrix(0L, nrow = rows, ncol = cols)
    k <- 1
    for (i in 1:dim_x[2]) {
      for (j in 1:dim_y[2]) {
        mat[,k] <- c(x[,i],y[,j])
        k <- k+1
      }
    }
    mat[!duplicated(mat),]
  }, list(...))
}