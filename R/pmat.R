#' Construct a Model Matrix for Poisson Regression
#' 
#' Determine the A matrix associated with a hierarchical model on a
#' contingency table for Poisson Regression.  
#' 
#' @param levels a list containing the number of levels of each
#'   variable
#' @param facets the facets generating the hierarchical model, a
#'   list of vectors of variable indices

#' @return The configuration matrix for the given levels and facets
#' 
#' @examples
#' # Single Covariate
#' levels <- 1:5
#' facets <- list(1)
#' pmat(levels, facets)
#' 
#' # multiple covariates, each has levels 1, ..., 5
#' levels <- list(1:5, 1:5)
#' facets <- list(1, 2, c(1,2))
#' pmat(levels, facets)
#' @export pmat

pmat <- function(levels, facets){

  #Small function to make single covariate configuration matrix
  func <- function(x) unname(rbind(rep(1, length(x)), x))
  
  if(!is.list(levels)){
    numCovariates <- 1
    fullMat <- func(levels)
  }else{
  numCovariates <- length(levels)
  #Make single covariate configuration matrix for each covariate
  matList <- lapply(levels, func)
  #Full heirarchicial config matrix with all interactions included
  fullMat <- do.call(kprod, matList)
  }
  expCov <- 1:numCovariates

  #Checking heirarchical sturcture of facets
  if(any(sapply(facets, length) > 1)){
    longListElts <- facets[which(sapply(facets, length) > 1)]
    
    uniqueVals <- unique(unlist(longListElts))
    
    heirarc <- as.list(c(uniqueVals, longListElts))
    
    facets <- union(heirarc, facets)
  }
  
  #All possible combinations of covariates (powerset like) to be compared to facets
   if(length(expCov) == 1) {
     facetList <- list(expCov)
   }else{
     facetList <- list(integer(0))
   for(i in seq_along(expCov)){
     facetList <- c(facetList, lapply(facetList, function(x) c(x,expCov[i])))
   }
     facetList <- facetList[-1]
   }
  #return the configuration matrix which includes only the elements need for the heirarchical model
  fullMat[c(TRUE, facetList %in% type.convert(facets)),]
}


