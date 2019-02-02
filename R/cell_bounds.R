#' Compute the lower and upper bounds
#'
#' Compute the lower and upper bounds of a cell in a vectorized contingency
#' table.
#'
#' @param A the configuration matrix of the model defining the fiber
#' @param b solution vector in Ax = b (vector of sufficient statistics)
#' @param strategy "lp" for linear programming and "ip" for integer programming.
#'   (both use \code{\link{lp}})
#' @param start reduce the size of the matrix A, i.e. use A[,start:ncol(A)]
#'   instead of A
#' @param messaging TRUE for messages
#' @param ... ...
#' @return a named numeric vector
#' @author Ruriko Yoshida \email{ruriko.yoshida@@uky.edu}, David Kahle
#'   \email{david.kahle@@gmail.com}
#' @name cell-bounds
#' @examples
#'
#'
#' A <- hmat(c(2,2), 1:2)
#' b <- rep.int(4, 4)
#' cell_bounds(A, b, messaging = TRUE)
#'
#'
#' 




#' @export
#' @rdname cell-bounds
cell_bounds <- function(A, b, strategy = c("lp", "ip"), start = 1, messaging = FALSE){ 
  
  ## arg check
  stopifnot(all(is.wholenumber(A)))
  stopifnot(all(is.wholenumber(b)))
  
  strategy <- match.arg(strategy)
  A <- A[,start:ncol(A), drop = FALSE]
  
  ## set dimensions of A
  m <- nrow(A); n <- ncol(A)
  
  ## setup objective function, see ?lpSolve::lp
  f.obj <- c(1L, rep(0L, n-1))
  
  ## setup directions
  f.dir <- c(rep("=", m), rep(">=", n))
  
  ## setup constraints
  f.con <- rbind(A, diag(n)) # diag(n) not int
  
  ## setup right hand side
  f.rhs <- c(b, rep(0L, n))

  if(messaging){
    message("(min/max)imize x[1]")
    message("s.t.")

    matLines <- apply(unname(f.con), 1, paste, collapse = " ")
   
    xs <- paste0("x[", format(1:ncol(f.con)), "]")
    xLines <- rep(str_dup(" ", nchar(xs[1])), length(matLines))
    xLines[1:length(xs)] <- xs

    message(
      paste(
        str_c(
          matLines, "   ", 
          xLines,   "   ", 
          format(f.dir, justify = "right"), "   ",
          format(f.rhs, justify = "right")
        ), 
        collapse = "\n"
      )
    )
    message("")
  }
  
  ## run lp
  lowerOut <- lp("min", f.obj, f.con, f.dir, f.rhs, all.int = (strategy == "ip"))    
  upperOut <- lp("max", f.obj, f.con, f.dir, f.rhs, all.int = (strategy == "ip"))    
  
  ## make list
  out <- vector("integer", 2)
  names(out) <- c("lower", "upper")
  
  ## populate list
  if(lowerOut$status == 2){
    out["lower"] <- NA
  } else {
    out["lower"] <- as.integer(ceiling(lowerOut$solution[1]))
  }
  
  if(upperOut$status == 2){
    out["upper"] <- NA
  } else {
    out["upper"] <- as.integer(floor(upperOut$solution[1]))
  }
  
  out
}








#' @export
#' @rdname cell-bounds
cellBounds <- function(...) {
  .Deprecated("cell_bounds")
  cell_bounds(...)
}





