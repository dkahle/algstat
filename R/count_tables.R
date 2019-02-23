#' Count Similarly Margined Contingency Tables
#'
#' Count the number of contingency tables with the same marginals as a given
#' table.
#'
#' \code{count_tables} uses LattE's count function (via algstat's
#' \code{\link{latte_count}} function) to count the tables.  In many cases, the
#' number of such tables is enormous.  In these cases, instead of giving back an
#' integer \code{count_tables} provides a character string with the integer in
#' it; see examples.
#'
#' @param table the table of interest
#' @param A the configuration/transpose design matrix
#' @param dir directory to place the files in, without an ending /
#' @param quiet show latte output
#' @param cache use count (default) or fcount
#' @param ... arguments to pass to \code{\link{latte_count}}
#' @return an integer
#' @seealso \code{\link{latte_count}}, \code{\link{count_fiber}}
#' @name count-tables
#' @examples
#' 
#' 
#' \dontrun{ requires LattE
#'
#'
#' data(politics); politics
#' count_tables(politics)
#' (A <- hmat(c(2,2), list(1, 2)))
#' count_tables(politics, A)
#'
#'
#'
#' data(handy); handy
#' count_tables(handy)
#'
#'
#'
#' data(HairEyeColor); HairEyeColor
#' eyeHairColor <- margin.table(HairEyeColor, 2:1)
#' count_tables(eyeHairColor)
#'
#' system.time(count_tables(eyeHairColor)) # it was computed above
#' system.time(count_tables(eyeHairColor)) # it was computed above
#' system.time(count_tables(eyeHairColor, cache = FALSE))
#' system.time(count_tables(eyeHairColor, cache = FALSE))
#'
#'
#' library(gmp)
#' as.bigz(count_tables(eyeHairColor))
#'
#'
#'
#' # notice that even tables with small cells can have huge fibers
#' data(drugs); drugs
#' count_tables(drugs)
#'
#'
#'
#'
#' # 0-1 tables can be very hard and produce very large fibers
#' # the 4x4 table below has 154 elements in its independence fiber
#' # the 5x5 has 16830, and the compute times are on the order of
#' # 1 and 10 seconds, respectively.
#' set.seed(1)
#' n <- 5
#' tab <- matrix(sample(0:1, n^2, replace = TRUE), nrow = n)
#' dimnames(tab) <- list(X = paste0("x", 1:n), Y = paste0("y", 1:n))
#' tab
#' count_tables(tab)
#'
#' count_tables(eyeHairColor, quiet = FALSE)
#'
#'
#'
#' }
#' 


#' @export
#' @rdname count-tables
count_tables <- function(table, 
    A = hmat(dim(table), as.list(1:length(dim(table)))), 
    dir = tempdir(), quiet = TRUE, cache = TRUE, ...
){
  
  ## make column names
  cellVars <- paste0("t", 1L:ncol(A))
  
  
  ## make the sums
  margConds <- unname(apply(A, 1, function(v){
    nonzero_ndcs <- unname(which(v > 0))
    one_ndcs     <- unname(which(v[nonzero_ndcs] == 1))
    terms <- paste(v[nonzero_ndcs], cellVars[nonzero_ndcs])
    terms[one_ndcs] <- str_sub(terms[one_ndcs], 3)
    paste(terms, collapse = " + ")
  }))

  
  
  ## compute the marginals
  marginals <- as.integer(A %*% tab2vec(table))
  
  
  ## make the equalities and inqualities
  margConds <- str_c(margConds, " == ", marginals)
  nonnegConds <- str_c(cellVars, " >= 0")
  
  
  ## count
  f <- if(cache) latte::latte_count else latte::latte_fcount
  f(c(margConds, nonnegConds), dir, quiet, ...)  
}




countTables <- function(...) {
  .Deprecated("count_tables")
  count_tables(...)
}








