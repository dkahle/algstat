#' Count the elements of a fiber Ax = b
#' 
#' Count the elements of a fiber Ax = b
#' 
#' \code{countFiber} uses LattE's count function (via algstat's 
#' \code{\link{count}} function) to count the fiber  In many cases,
#' the number of such tables is enormous.  In these cases, instead
#' of giving back an integer \code{countFiber} provides a character
#' string with the integer in it; see examples.
#' 
#' @param A the A matrix of Ax = b
#' @param b the b vector of Ax = b
#' @param dir directory to place the files in, without an ending /
#' @param opts options for count
#' @param quiet show latte output
#' @return an integer
#' @seealso \code{\link{count}}
#' @export
#' @name countFiber
#' @examples
#' 
#' \dontrun{ 
#' 
#' 
#' data(politics)
#' (A <- hmat(c(2,2), list(1, 2)))
#' countTables(politics, A)
#' b <- A %*% tab2vec(politics)
#' countFiber(A, b)
#' 
#' 
#' 
#' 
#' }
#' 
countFiber <- function(A, b, dir = tempdir(), opts = "", quiet = TRUE){
  
  ## basic quantities
  m <- nrow(A); n <- ncol(A)
  
  ## turn into latte format
  bNegA <- unname(cbind(b, -A))
  
  ## add linearity and nonnegative attributes
  attr(bNegA, "linearity") <- 1:m
  attr(bNegA, "nonnegative") <- 1:n
  
  ## make code
  code <- write.latte(bNegA)
  
  ## count
  count(spec = code, dir = dir, opts = opts, quiet = quiet)  
}












#' @param ... ...
#' @export
#' @rdname countFiber
memCountFiber <- memoise::memoise(countFiber)









