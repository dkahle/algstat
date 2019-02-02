#' Lp Norm
#' 
#' Compute the Lp norm of a vector.
#' 
#' @param x x
#' @param p p
#' @return ...
#' @export lpnorm
#' @examples
#' 
#' lpnorm(1:10)
#' lpnorm(matrix(1:25, 5, 5))
#' lpnorm(split(1:25, rep(1:5, each = 5)))
#' 
#' lpnorm(1:10, 1)
#' lpnorm(matrix(1:25, 5, 5), 1)
#' lpnorm(split(1:25, rep(1:5, each = 5)), 1)
#' 
#' lpnorm(rnorm(10), 0)
#' lpnorm(matrix(rnorm(25), 5, 5), 0)
#' lpnorm(split(rnorm(25), rep(1:5, each = 5)), 0)
#' 
#' lpnorm(-5:5, Inf)
#' lpnorm(matrix(-25:-1, 5, 5), Inf)
#' lpnorm(split(-25:-1, rep(1:5, each = 5)), Inf)
#' 
lpnorm <- function(x, p = 2){
  if(is.vector(x) && !is.list(x)){
  	if(p == Inf) return(max(abs(x)))
    if(p >= 1) return( sum(abs(x)^p)^(1/p) )
    if(0 <= p && p < 1) return( sum(abs(x)^p) )
  }
  if(is.matrix(x)) return(apply(x, 1, lpnorm, p))
  if(is.list(x)) return(sapply(x, lpnorm, p))
  NA  
}








#' Random Spectral Data
#'
#' Generate spectral data for testing purposes.
#' 
#' @param nVoters number of voters voting
#' @param nObjects number of objects up for selection
#' @param kSelected number of objects selected by each voter
#' @return ...
#' @export rvotes
#' @examples
#' rvotes(100, 10, 3)
#'
rvotes <- function(nVoters, nObjects, kSelected){
  t(replicate(nVoters, sort(sample(nObjects, kSelected))))
}







#' Create an upper triangular matrix
#'
#' Create an upper triangular matrix.
#' 
#' @param x a vector
#' @return ...
#' @seealso \code{\link{lower}}
#' @export upper
#' @examples
#' upper(1:3)
#' lower(1:3)
#'
#' upper(1:6)
#' lower(1:6)
#' 
#' upper(rnorm(6))
#'
upper <- function(x){
  l <- length(x)
  p <- round( (1+sqrt(1+8*l))/2 )	
  ndcs <- unlist(sapply(as.list(0:(p-2)), function(k){
    k*p + (k+2):p
  }))  
  
  if(all(is.integer(x))){
    m <- as.integer( matrix(rep(0, p^2), p, p) )
    m[ndcs] <- as.integer(x)
  } else {
    m <- matrix(rep(0, p^2), p, p)
    m[ndcs] <- x
  }
  
  dim(m) <- c(p,p)
  t(m)
}




#' Create a lower triangular matrix
#'
#' Create a lower triangular matrix.
#' 
#' @param x a vector
#' @return ...
#' @export lower
#' @seealso \code{\link{upper}}
#' @examples
#' upper(1:3)
#' lower(1:3)
#'
#' upper(1:6)
#' lower(1:6)
#' 
#' upper(rnorm(6))
#'
lower <- function(x) t(upper(x))










#' Vector Projection
#'
#' Project a vector onto the column space of a matrix or the orthogonal complement of the column
#' space of a matrix; the null space of A transpose.
#' 
#' @param A a matrix
#' @param x a vector
#' @return ...
#' @name project-onto
#' @seealso \code{\link{qr.fitted}}
#' @examples
#' 
#' A <- diag(5)[,1:2]
#' x <- 1:5
#' project_onto(A, x)
#' project_onto_perp(A, x)
#' 


#' @rdname project-onto
#' @export
project_onto <- function(A, x) qr.fitted(qr(A), x) 


#' @rdname project-onto
#' @export
project_onto_perp <- function(A, x) as.vector(diag(length(x))%*%x - qr.fitted(qr(A), x) )













#' Multinomial Coefficient
#'
#' Compute the multinomial coefficient.
#'
#' This function computes the multinomial coefficient by computing the factorial
#' of each number on a log scale, differencing log(n!) - sum(log(x!)), and then
#' exponentiating.  It then checks to see if this is an integer; if it's not, it
#' issues a warning.
#'
#' @param n an integer
#' @param x a vector of integers
#' @return ...
#' @export mchoose
#' @examples
#'
#' mchoose(6, c(2,2,1,1))
#'
#'
#'
#' 
mchoose <- function(n, x){
  lboth <- lfactorial(c(n,x))
  out <- exp(lboth[1] - sum(lboth[-1]))
  if(out != round(out)) warning("log-retransformed multinomial coefficient != rounded value.")
  as.integer(round(out))
}






















timeStamp <- function(){
  timeStamp <- as.character(Sys.time())
  timeStamp <- chartr("-", "_", timeStamp)
  timeStamp <- chartr(" ", "_", timeStamp)
  timeStamp <- chartr(":", "_", timeStamp)
  timeStamp
}






is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol




file.path2 <- function(...){
  dots <- list(...)
  if(.Platform$OS.type == "unix"){
    sep <- "/"
  } else {
    sep <- "\\"
  }
  paste0(dots, collapse = sep)
}



is.formula <- function(x) class(x) == "formula"





is.mac <- function() grepl("darwin", R.version$platform)
is.win <- function() .Platform$OS.type == 'windows'
is.unix <- function() .Platform$OS.type == "unix"







sample_int_between <- function(l, u){
  if(l == u) return(l)
  sample(l:u, 1)
}





capitalize <- function(s){
  if(length(s) > 1) return(vapply(s, capitalize, character(1)))  
  str_c(toupper(str_sub(s, 1, 1)), str_sub(s, 2))
}












`%notin%` <- function(elem, set){
  if(length(elem) > 1) return(vapply(elem, `%notin%`, logical(1), set = set))
  !(elem %in% set)
}
















#' Test whether an mpoly object is linear.
#'
#' Test whether an mpoly object is linear.
#'
#' @param x an mpoly or mpolyList object
#' @return a logical vector
#' @examples
#' 
#' 
#' is.linear(mp("0"))
#' is.linear(mp("x + 1"))
#' is.linear(mp("x + y"))
#' is.linear(mp(c("0", "x + y")))
#' 
#' is.linear(mp("x + x y"))
#' is.linear(mp(c("x + x y", "x")))
#' 
#' 
is.linear <- function(x){
  
  stopifnot(is.mpoly(x) || is.mpolyList(x))
  
  if(is.mpolyList(x)) return(sapply(x, is.linear))
  
  all(
    vapply(x, function(term){
      if(all(length(term) <= 2)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }, logical(1))
  )  
}













#' Fill a table with a number
#' 
#' Fill a table with a number
#' 
#' @param tab a contingency table
#' @param fill the number to fill the contingency table with
#' @param ... ...
#' @return a named numeric vector
#' @name tab-fill
#' @examples
#' 
#' Titanic
#' tab_fill(Titanic)
#' tab_fill(Titanic, 1L)
#' tab_fill(Titanic, 0L)
#' 
#' nCells <- prod(dim(Titanic))
#' tab_fill(Titanic, rpois(nCells, 5))
#' 




#' @export
#' @rdname tab-fill
tab_fill <- function(tab, fill = 1L){
  tab[] <- fill
  tab
}


#' @export
#' @rdname tab-fill
tabFill <- function (...) {
  .Deprecated("tab_fill")
  tab_fill(...)
}








#' @importFrom latte tab2vec
#' @export
latte::tab2vec




#' @importFrom latte vec2tab
#' @export
latte::vec2tab










