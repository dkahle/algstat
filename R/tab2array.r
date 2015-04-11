#' Table to array conversion
#'
#' Convert a table into an array.
#' 
#' @param tab a table
#' @return an array
#' @export tab2array
#' @seealso \code{\link{tab2vec}}, \code{\link{array2tab}}
#' @examples
#' 
#' data(handy)
#' handy
#' tab2array(handy)
#' 
#' data(Titanic)
#' Titanic
#' tab2array(Titanic)
#' 
#'
tab2array <- function(tab){
  stopifnot(is.table(tab))
  a <- tab
  attributes(a)[c("class","dimnames")] <- NULL
  a
}

















#' Array to table conversion
#'
#' Convert an array into a table.
#' 
#' @param array an array
#' @return an tabl
#' @export array2tab
#' @seealso \code{\link{tab2array}}
#' @examples
#' 
#' mat <- matrix(rpois(25,2), 5, 5) 
#' array2tab(mat)
#' 
#' data(handy)
#' array2tab(handy)
#' 
#'
array2tab <- function(array){
  
  if(is.table(array)){
    message("this is already a table!")
    return(array)
  }
  
  stopifnot(is.array(array))
  
  ## initialize list
  dimnamesList <- as.list(1:length(dim(array)))
  
  ## give names as X1, X2, and so on
  names(dimnamesList) <- paste0("X", 1:length(dimnamesList))
  
  ## populate list with X1-1, X1-2, and so on
  for(k in 1:length(dimnamesList)){
    dimnamesList[[k]] <- paste0("X", k, "-", 1:dim(array)[k])
  }
  
  ## set dimnames
  dimnames(array) <- dimnamesList
  
  ## return
  array
}


