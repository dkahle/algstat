#' Construct a Model Matrix for Poisson Regression
#' 
#' Determine the A matrix associated with a hierarchical model on a
#' contingency table for Poisson Regression.  
#' 
#' @param levels a vector containing the number of levels of each
#'   variable
#' @param facets the facets generating the hierarchical model, a
#'   list of vectors of variable indices

#' @return a matrix
#' @export pmat

pmat <- function(levels, facets){

  levels <- levels[levels != 1]
  num_covariates <- length(levels)
  mat_list <- list()
  
  for(i in 1:num_covariates){
    mat_list[[i]] <- matrix(c(rep(1,levels[i]), 1:levels[i]), nrow = 2, byrow = T)
  }
  
  return(do.call(kprod, mat_list))

}


