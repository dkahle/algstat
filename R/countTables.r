#' Count Similarly Margined Contingency Tables
#' 
#' Count the number of contingency tables with the same marginals as
#' a given table.
#' 
#' \code{countTables} uses LattE's count function (via algstat's 
#' \code{\link{count}} function) to count the tables.  In many
#' cases, the number of such tables is enormous.  In these cases,
#' instead of giving back an integer \code{countTables} provides a
#' character string with the integer in it; see examples.
#' 
#' @param table the table of interest
#' @param A the configuration/transpose design matrix
#' @param dir directory to place the files in, without an ending /
#' @param quiet show latte output
#' @param cache use count (default) or fcount
#' @param ... arguments to pass to \code{\link{count}}
#' @return an integer
#' @seealso \code{\link{count}}, \code{\link{countFiber}}
#' @export
#' @name countTables
#' @examples
#' \dontrun{ 
#' 
#' 
#' data(politics)
#' countTables(politics)
#' (A <- hmat(c(2,2), list(1, 2)))
#' countTables(politics, A)
#' 
#' 
#' 
#' data(handy)
#' countTables(handy)
#' 
#' 
#' 
#' data(HairEyeColor)
#' eyeHairColor <- margin.table(HairEyeColor, 2:1)
#' countTables(eyeHairColor)
#' 
#' system.time(countTables(eyeHairColor)) # it was computed above
#' system.time(countTables(eyeHairColor)) # it was computed above
#' system.time(countTables(eyeHairColor, cache = FALSE))
#' system.time(countTables(eyeHairColor, cache = FALSE))
#' 
#' 
#' library(gmp)
#' as.bigz(countTables(eyeHairColor))
#' 
#' 
#' 
#' # notice that even tables with small cells can have 
#' # huge fibers
#' data(drugs)
#' countTables(drugs)
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
#' countTables(tab)
#' 
#' countTables(eyeHairColor, quiet = FALSE)
#' 
#' 
#' 
#' }
#' 
countTables <- function(table, 
    A = hmat(dim(table), as.list(1:length(dim(table)))), 
    dir = tempdir(), quiet = TRUE, cache = TRUE, ...
){
  
  ## make column names
  cellNames <- paste0("t", colnames(A))
  
  
  ## make the sums
  margConds <- unname(apply(A, 1, function(v){
    nonzero_ndcs <- unname(which(v > 0))
    one_ndcs     <- unname(which(v[nonzero_ndcs] == 1))
    terms <- paste(v[nonzero_ndcs], cellNames[nonzero_ndcs])
    terms[one_ndcs] <- str_sub(terms[one_ndcs], 3)
    paste(terms, collapse = " + ")
  }))

  
  
  ## compute the marginals
  marginals <- as.integer(A %*% tab2vec(table))
  
  
  ## make the equalities and inqualities
  margConds <- paste0(margConds, " == ", marginals)
  nonnegConds <- paste0(cellNames, " >= 0")
  
  
  ## count
  f <- if(cache) latter::count else latter::fcount
  f(c(margConds, nonnegConds), dir, quiet, ...)  
}










