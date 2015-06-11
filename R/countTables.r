#' Count Similarly Margined Contingency Tables
#' 
#' Count the number of contingency tables with the same marginals as a given
#' table.
#' 
#' \code{countTables} uses LattE's count function (via algstat's
#' \code{\link{count}} function) to count the tables.  In many cases, the number
#' of such tables is enormous.  In these cases, instead of giving back an
#' integer \code{countTables} provides a character string with the integer in
#' it; see examples.
#' 
#' @param table the table of interest
#' @param A the configuration/transpose design matrix
#' @param dir directory to place the files in, without an ending /
#' @param opts options for count
#' @param quiet show latte output
#' @return an integer
#' @seealso \code{\link{count}}
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
#' system.time(countTables(eyeHairColor))
#' system.time(countTables(eyeHairColor))
#' system.time(memCountTables(eyeHairColor))
#' system.time(memCountTables(eyeHairColor))
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
#' countTables(eyeHairColor, quiet = FALSE)
#' 
#' 
#' 
#' }
#' 
countTables <- function(table, 
    A = hmat(dim(table), as.list(1:length(dim(table)))), 
    dir = tempdir(), opts = "", quiet = TRUE
){
  
  ## make column names
  cellNames <- paste0("t", colnames(A))
  
  
  ## make the sums
  margConds <- unname(apply(A, 1, function(v){
    paste(paste(v, cellNames), collapse = " + ")
  }))
  
  
  ## compute the marginals
  marginals <- as.integer(A %*% tab2vec(table))
  
  
  ## make the equalities and inqualities
  margConds <- paste0(margConds, " == ", marginals)
  nonnegConds <- paste0(cellNames, " >= 0")
  
  
  ## count
  count(c(margConds, nonnegConds), dir, opts, quiet)  
}












#' @param ... ...
#' @export
#' @rdname countTables
memCountTables <- memoise::memoise(countTables)









