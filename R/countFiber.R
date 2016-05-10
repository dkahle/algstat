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
#' @param quiet show latte output
#' @param cache use count (default) or fcount
#' @param ... additional arguments to \code{\link{count}}
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
#' # this counts the number of ways 10 non-negative integers can sum to 100
#' A <- ones(1, 10)
#' countFiber(A, 100, cache = FALSE)
#' system.time(countFiber(A, 100, cache = FALSE))
#' system.time(countFiber(A, 100))
#' system.time(countFiber(A, 100))
#' 
#' 
#' }
#' 
countFiber <- function(A, b, dir = tempdir(), quiet = TRUE, cache = TRUE, ...){
  
  ## basic quantities
  m <- nrow(A); n <- ncol(A)
  
  ## turn into latte format
  bNegA <- unname(cbind(b, -A))
  
  ## add linearity and nonnegative attributes
  attr(bNegA, "linearity") <- 1:m
  attr(bNegA, "nonnegative") <- 1:n
  
  ## make code
  code <- format_latte(bNegA)

  ## count
  f <- if(cache) latter::count else latter::fcount
  f(spec = code, dir = dir, quiet = quiet, ...)  
}











