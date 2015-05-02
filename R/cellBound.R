#' Compute the lower or upper bound
#' 
#' Compute the lower or upper bound of a cell in a vectorized contingency table.
#' 
#' @param A the configuration matrix of the model defining the fiber
#' @param b solution vector in Ax = b (vector of sufficient statistics)
#' @param side bound side to compute, "lower" or "upper"
#' @param strategy "lp" for linear programming and "ip" for integer programming.  
#' (both use \code{\link{lp}}.)
#' @return an integer
#' @author Ruriko Yoshida \email{ruriko.yoshida@@uky.edu}, David Kahle \email{david.kahle@@gmail.com}
#' @export 
#' @examples
#' 
#' 
#' A <- hmat(c(2,2), 1:2)
#' b <- rep.int(4, 4)
#' cellBound(A, b, "lower") # 0
#' cellBound(A, b, "upper") # 4
#' 
#' 
#' 
#' 
cellBound <- function(A, b, side = c("lower", "upper"), strategy = c("lp", "ip")){ 
  
  ## arg check
  side   <- match.arg(side)
  strategy <- match.arg(strategy)
  
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
  
  ## run lp
  out <- lp(
    ifelse(side == "lower", "min", "max"), 
    f.obj, f.con, f.dir, f.rhs, 
    all.int = (strategy == "ip")
  )
  
  ## clean output and return
  if(out$status == 2) return(NA)
  
  if(side == "lower"){
    o <- as.integer(ceiling(out$solution[1]))
  } else if(side == "upper"){
    o <- as.integer(floor(out$solution[1]))
  }
  
  o
}


