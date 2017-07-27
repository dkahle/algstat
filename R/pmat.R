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
  
  #Basic Variables
  num_covariates <- length(levels)
  cov_index <- 1:num_covariates
  num_elements <- prod(levels)
  num_rows <- length(facets)
  
  #Base level combos 
  varsNlvls <- lapply(as.list(levels), function(x) number2Glyph(1:x))
  baseLvls  <- expand.grid(rev(varsNlvls))[,num_covariates:1]
  
  #Make initial matrix
  A <- matrix(0L, nrow = num_rows, ncol = num_elements)
  
  for(i in 1:num_rows){
    if(length(facets[[i]]) == 1){

      A[i, ] <- baseLvls[,facets[[i]]]
    }else{
      inter_terms <- cov_index %in% facets[[i]]
      sub_mat <- as.matrix(baseLvls[ ,inter_terms])
      class(sub_mat) <- "numeric"
      A[i,] <- apply(sub_mat, 1, prod)
    }
  }
  return(rbind(rep(1, num_elements), A))
}

number2Glyph <- function(n) c(0:9, letters, LETTERS)[n+1]
