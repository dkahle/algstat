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

  #########Old Setup#########
  #Setup the levels
  #levels <- levels[levels != 1]
  #num_covariates <- length(levels)
  #mat_list <- list()
  #exp_cov <-1:num_covariates
 #############################
  #Small function to make single covariate configuration matrix
  func <- function(x){
    rbind(rep(1, length(x)), x)
  }
  
  ######### New Setup ########
  if(is.vector(levels)){
    num_covariates <- 1
    full_mat <- func(levels)
  }else{
  num_covariates <- ncol(levels)
  #Make single covariate configuration matrix for each covariate
  mat_list <- alply(levels, 2, func)
  #Full heirarchicial config matrix with all interactions included
  full_mat <- do.call(kprod, mat_list)
  }
  exp_cov <- 1:num_covariates
  
  #Checking heirarchical sturcture of facets
  if(any(sapply(facets, length) > 1)){
    long_list_elts <- facets[which(sapply(facets, length)>1)]
    
    unique_vals <- unique(unlist(long_list_elts))
    
    heirarc <- as.list(c(unique_vals, long_list_elts))
    
    facets <- union(heirarc, facets)
  }
  
  ########Old#######
  #List of config matrices, one for each covariate
  #for(i in exp_cov){
  #  mat_list[[i]] <- matrix(c(rep(1,levels[i]), 1:levels[i]), nrow = 2, byrow = TRUE)
  #}
  ################
  
  
  #All possible combinations of covariates (powerset like) to be compared to facets
   if(length(exp_cov) == 1) {
     facet_list <- list(exp_cov)
   }else{
     facet_list <- list(integer(0))
   for(i in seq_along(exp_cov)){
     facet_list <- c(facet_list, lapply(facet_list, function(x) c(x,exp_cov[i])))
   }
     facet_list <- facet_list[-1]
   }
  #return the configuration matrix which includes only the elements need for the heirarchical model
  return(full_mat[c(TRUE, facet_list %in% facets),])
}


