#' Construct a Hierarchical Model Matrix
#'
#' Determine the A matrix associated with a hierarchical model on a contingency
#' table.  In algebraic statistics, the A matrix of a log-linear model is the
#' transpose of the design matrix of the (cell-means parameterized) ANOVA
#' corresponding to the model.
#'
#' @param varlvls a vector containing the number of levels of each variable
#' @param facets the facets generating the hierarchical model, a list of vectors
#'   of variable indices
#' @seealso \code{\link{genmodel}}
#' @return a named matrix
#' @export hmat
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures
#'   on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @examples
#'
#' # 2x2 independence example
#' # following convention, the first index indicates rows
#' varlvls <- c(2,2)
#' facets <- list(1,2)
#' ( A <- hmat(varlvls, facets) )
#'
#' # alternatively:
#' hmat(c(2, 2), 1:2)
#'
#'
#' # LAS example 1.2.11, p.16
#' varlvls <- c(2,2,2,2)
#' facets <- list(c(1,2), c(1,4), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#'
#' cat(format_latte(A))
#' 
hmat <- function(varlvls, facets){

  # set basic variables
  p <- length(varlvls) # number of variables in table (p-way)  
  numFacets <- length(facets)
  
  
  # make cells
  varsNlvls <- lapply(as.list(varlvls), function(x) number2Glyph(1:x))
  cellsDf   <- expand.grid(rev(varsNlvls))[,p:1]  
  
  
  # make colnames
  colNames <- apply(cellsDf, 1, paste, collapse = "")
  nCells   <- length(colNames)
  
  
  # make rownames
  configsPerFacet <- vapply(facets, function(facet) prod(varlvls[facet]), numeric(1))
  totalRows       <- sum(configsPerFacet)
  
  facetConfigs <- lapply(facets, function(facet){
  	tmpDf <- expand.grid(rev.default(varsNlvls[facet]))[,length(facet):1, drop = FALSE]  
  	names(tmpDf) <- facet
  	tmpDf
  })
  
  rowNames <- unlist(lapply(facets, function(facet){
    facetConfigs <- expand.grid(rev(varsNlvls[facet]))[,length(facet):1, drop = FALSE]  
    facetConfigsMat <- matrix('+', nrow = nrow(facetConfigs), ncol = p)
    facetConfigsMat[,facet] <- as.matrix(facetConfigs)
    facetConfigs <- apply(facetConfigsMat, 1, paste, collapse = "")
    facetConfigs
  }))
  
  
  # make A
  A <- matrix(0L, nrow = totalRows, ncol = nCells, dimnames = list(rowNames, colNames))

  
  # put in 0's and 1's
  for(k in 1:totalRows){
    ndcsToMatch   <- which(strsplit(rowNames[k],"")[[1]] != "+")
    configToMatch <- gsub('\\+', "", rowNames[k])
    A[k,] <- vapply(strsplit(colNames, ""), function(l){
      paste(l[ndcsToMatch], collapse = "") == configToMatch
    }, logical(1)) + 0
  }
  
  
  # convert to ints
  intA <- as.integer(A)
  intA <- matrix(intA, nrow = nrow(A), ncol = ncol(A), dimnames = list(rowNames, colNames))

  
  # return
  intA
}





number2Glyph <- function(n) c(0:9, letters, LETTERS)[n+1]

glyph2Number <- function(g){
  x <- 0:62
  names(x) <- c(0:9, letters, LETTERS)
  unname(x[g])
}







