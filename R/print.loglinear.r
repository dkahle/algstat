#' Pretty printing of loglinear's output
#'
#' Pretty printing of loglinear's output.
#' 
#' @param x an object of class loglinear
#' @param digits digits to round to
#' @param ... additional parameters
#' @usage \method{print}{loglinear}(x, digits = 4, ...)
#' @return Invisible string of the printed object.
#' @seealso [loglinear()]
#' @export
#' @examples
#' 
#' # see ?loglinear
#' 
#' 
print.loglinear <- function(x, digits = 4, ...){
	
  # argument checking and basic variable setting
  stopifnot(class(x) == "loglinear")
  
  # print call
  cat("Call:\n")  
  print(x$call)
  cat("\n")
  
  # print method
  cat("Fitting method:\n")  
  if(x$method == "mcmc") cat("Metropolis-Hastings using Markov moves\n")
  if(x$method == "ipf") cat("Iterative proportional fitting (with stats::loglin)\n")
  cat("\n")
  
  # print method
  cat("MCMC details:\n")  
  cat(paste0("N = ", x$iter, " samples (after thinning), burn in = ", x$burn, ", thinning = ", x$thin))
  cat("\n\n")  
  
  # print statistics
  statMat <- cbind(
    "stat" = round(x$statistic, digits), 
    "P(>= stat)" = round(x$p.value, digits),
    "se" = round(x$p.value.std.err, digits), 
    "mid-P-Value" = round(x$mid.p.value, digits)
  )
  statMat[1,1] <- ""
  statMat <- cbind(
    format(c("P(table)", "Pearson X^2", "Likelihood G^2", 
      "Freeman-Tukey", "Cressie-Read", "Neyman X^2"), justify = "right"),
    statMat
  )
  statMat <- apply(statMat, 2, format, justify = "left")
  statMat <- rbind(
    c("", "stat", "p.value", "se", "mid.p.value"),
    statMat
  )
  statMat <- apply(statMat, 2, format, justify = "right") 
  apply(statMat, 1, function(x){cat(x); cat("\n")})
  invisible()  
}



