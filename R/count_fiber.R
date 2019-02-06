#' Count the elements of a fiber Ax = b
#'
#' Count the elements of a fiber Ax = b
#'
#' \code{count_fiber} uses LattE's count function (via algstat's
#' \code{\link{latte_count}} function) to count the fiber  In many cases, the
#' number of such tables is enormous.  In these cases, instead of giving back an
#' integer \code{count_fiber} provides a character string with the integer in
#' it; see examples.
#'
#' @param A the A matrix of Ax = b
#' @param b the b vector of Ax = b
#' @param dir directory to place the files in, without an ending /
#' @param quiet show latte output
#' @param cache use count (default) or fcount
#' @param ... additional arguments to \code{\link{latte_count}}
#' @return an integer
#' @seealso \code{\link{latte_count}}
#' @name count-fiber
#' @examples
#'
#' \dontrun{ requires LattE
#'
#' data(politics); politics
#' 
#' (A <- hmat(c(2,2), list(1, 2)))
#' count_tables(politics, A)
#' b <- A %*% tab2vec(politics)
#' count_fiber(A, b)
#'
#'
#' # this counts the number of ways 10 non-negative integers can sum to 100
#' A <- ones(1, 10)
#' count_fiber(A, 100, cache = FALSE)
#' system.time(count_fiber(A, 100, cache = FALSE))
#' system.time(count_fiber(A, 100))
#' system.time(count_fiber(A, 100))
#'
#'
#' }
#' 
#' 


#' @export
#' @rdname count-fiber
count_fiber <- function(A, b, dir = tempdir(), quiet = TRUE, cache = TRUE, ...){
  
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
  f <- if(cache) latte::latte_count else latte::latte_fcount
  f(spec = code, dir = dir, quiet = quiet, ...)  
}






#' @export
#' @rdname count-fiber
countFiber <- function(...) {
  .Deprecated("count_fiber")
  count_fiber(...)
}





