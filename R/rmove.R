#' Generate random moves on a fiber
#' 
#' Generate random moves on a fiber
#' 
#' @param n the number of observations
#' @param A the configuration matrix of the model defining the fiber
#' @param b integer vector; the vector of sufficient statistics
#' @param x integer vector; a vector in the fiber to be used as an 
#'   initialization point for the "walk" method
#' @param method "sis", "walk", or "hybrid"
#' @param dist "uniform" or "hypergeometric"
#' @param parallel parallelize the sampling?
#' @param format format of the returned moves, "mat", "vec", or "tab"
#' @param dim the dimensions of the table if "tab" is used, see 
#'   \code{\link{vec2tab}}
#' @param ... arguments to pass to \code{\link{rfiber}}, such as \code{parallel 
#'   = TRUE}
#' @return a named numeric vector
#' @author David Kahle \email{david.kahle@@gmail.com}, Ruriko Yoshida 
#'   \email{ruriko.yoshida@@uky.edu}
#' @seealso \code{\link{rfiber}}
#' @export rmove
#' @examples
#' 
#' 
#' data(politics)
#' politics
#' x <- tab2vec(politics)
#' x
#' 
#' A <- hmat(c(2,2), 1:2) # independence model on 2x2 table
#' b <- A %*% x # vector of sufficient statistics
#' rmove(10, A = A, b = b)
#' rmove(10, A = A, x = x, method = "walk")
#' 
#' 
#' 
#' 
#' moves <- rmove(10000, A = A, x = x, method = "walk")
#' unique(t(moves))
#' 
#' 
#' ## simple example
#' A <- hmat(c(2,2), 1:2)
#' x <- c(1, 3, 3, 1)
#' b <- rep.int(4, 4)
#' 
#' rfiber(10, A = A, b = b)
#' (tabs <- rfiber(10, A = A, b = b, format = "vec"))
#' unique(tabs)
#' lapply(unique(tabs), vec2tab, dim = c(2, 2))
#' 
#' 
#' 
#' \dontrun{ # save R CMD check time and requires outside software
#' 
#' 
#' library(microbenchmark)
#' microbenchmark(
#'   rmove(100, A = A, b = b),
#'   rmove(100, A = A, x = x, method = "walk")
#' )
#' 
#' 
#' 
#' markov(A)
#' markov(A, "vec")
#' markov(A, "tab", c(2,2))
#' 
#' rmove(5, A, b)
#' rmove(5, A, b, format = "vec")
#' rmove(5, A, b, format = "tab", dim = c(2,2))
#' 
#' 
#' 
#' 
#' ## politics example
#' data(politics)
#' politics
#' tab2vec(politics)
#' b <- A %*% tab2vec(politics)
#' tabs <- rfiber(1000, A = A, b = b, format = "vec")
#' lapply(unique(tabs), vec2tab, dim = c(2, 2))
#' 
#' rmove(5, A, b)
#' rmove(5, A, b, format = "tab", dim = c(2,2))
#' 
#' 
#' 
#' ## parallelizing
#' system.time(rmove(1e4, A, b)) # ~35s
#' system.time(rmove(1e4, A, b, parallel = TRUE)) # ~7.5s, 8 cores
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ## drugs example
#' data(drugs)
#' A <- hmat(c(2,2,2), 1:3)
#' b <- A %*% tab2vec(drugs)
#' markovMoves <- markov(A)
#' rmove(5, A, b)
#' rmove(5, A, b, format = "tab", dim = c(2,2,2))
#' 
#' 
#' 
#' }
#' 
rmove <- function(n, A, b, x,
  method = c("sis", "walk", "hybrid"), 
  dist = c("uniform", "hypergeometric"), 
  parallel = FALSE,
  format = c("mat" ,"vec", "tab"), dim,
  ...
){
  
  # arg check
  method <- match.arg(method)
  dist   <- match.arg(dist)
  format <- match.arg(format)
  
  # generate 2*n tables
  if(missing(x)){
    tabs <- rfiber(2*n, A = A, b = b, method = method, dist = dist, parallel = parallel, format = "mat", ...)
  } else if(missing(b)){
    tabs <- rfiber(2*n, A = A, x = x, method = method, dist = dist, parallel = parallel, format = "mat", ...)
  }
  
  
  # subtract them and split the results
  moves <- tabs[,1:n] - tabs[,(n+1):(2*n)]
  
  # rename and out
  if(format == "mat"){
    return(moves)
  } else if(format == "vec"){
    moves <- split(moves, rep(1:n, each = ncol(A)))
    names(moves) <- NULL
    return(moves)
  } else if(format == "tab"){
    moves <- split(moves, rep(1:n, each = ncol(A)))
    names(moves) <- NULL
    moves <- lapply(moves, vec2tab, dim = dim)
    return(moves)
  }

}





