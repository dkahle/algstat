#' Array to vector conversion
#' 
#' Convert an array into a vector.
#' 
#' This function converts an array (or a multi-way contingency table) into a
#' vector, using a consistent ordering of the cells.  The ordering of the cells
#' is lexicographical and cannot be specified by the user.
#' 
#' @param tab an array of counts
#' @return a named integer vector.  the names correspond to the cell indices in
#'   the table.
#' @export tab2vec
#' @seealso \code{\link{vec2tab}}
#' @examples
#' 
#' a <- array(1:24, c(2,3,4))
#' tab2vec(a)
#' 
#' data(Titanic)
#' tab2vec(Titanic)
#' Titanic[1,1,1,1]
#' Titanic[1,1,1,2]
#' 
#' 
tab2vec <- function(tab){
  
  # if is a vector, return
  if(is.null(dim(tab))){
    tab <- as.vector(tab)
    names(tab) <- 1:length(tab)
    return(tab)
  }
  
  # if it's a vector already, short-circuit
  if(length(dim(tab)) == 1){
    tab <- as.vector(tab)
    names(tab) <- 1:length(tab)
    return(tab)
  }
  
  # otherwise, rotate and class
  u <- aperm(tab, length(dim(tab)):1)
  if(class(tab[1]) == "numeric") u <- as.vector(u)
  if(class(tab[1]) == "integer") u <- as.integer(u)  
  
  # create cell indices
  tmpdf <- expand.grid( 
    rev.default(lapply(dim(tab), function(x) 1:x)) 
  )[,length(dim(tab)):1]
  
  # assign them as names to u
  names(u) <- apply(tmpdf, 1, paste, collapse = ',')
  
  # return
  u
}







#' Vector to array conversion
#'
#' Convert a vector into an array given a set of dimensions; it therefore simply wraps \code{aperm} and \code{array}.
#'
#' This function converts an array (or a multi-way contingency table) into a vector, using a consistent ordering of the cells.  The ordering of the cells is lexicographical and cannot be specified by the user.
#' 
#' @param vec a vector
#' @param dim the desired array dimensions, oftentimes a vector of the number of levels of each variable in order
#' @return an array
#' @export vec2tab
#' @seealso \code{\link{tab2vec}}, \code{\link{aperm}}, \code{\link{array}}
#' @examples
#' 
#' data(Titanic)
#' Titanic
#' tab2vec(Titanic)
#' vec2tab(tab2vec(Titanic), dim(Titanic))
#' vec2tab(tab2vec(Titanic), dim(Titanic)) == Titanic
#' all(vec2tab(tab2vec(Titanic), dim(Titanic)) == Titanic)
#' 
#'
vec2tab <- function(vec, dim){
  aperm(
    array(vec, rev(dim)),
    length(dim):1
  )
}


















