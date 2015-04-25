#' Fill a table with a number
#' 
#' Fill a table with a number
#' 
#' @param tab a contingency table
#' @param fill the number to fill the contingency table with
#' @return a named numeric vector
#' @export tabFill
#' @examples
#' 
#' Titanic
#' tabFill(Titanic)
#' tabFill(Titanic, 1)
#' tabFill(Titanic, 0)
#' 
#' nCells <- prod(dim(Titanic))
#' tabFill(Titanic, rpois(nCells, 5))
#' 
#' 
tabFill <- function(tab, fill = 1L){
  tab[] <- fill
  tab
}