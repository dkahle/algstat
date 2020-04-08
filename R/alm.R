#' Algebraic least squares regression
#'
#' Algebraic least squares regression
#'
#' @param form A formula specifying the basis monomials, passed to
#'   [stats::model.matrix]
#' @param data The dataset on which to apply the method, passed to
#'   [stats::model.matrix]
#' @param ... Additional parameters, presently discarded
#' @return A numeric vector of model parameter estimates
#' @name alm
#' @examples
#'
#' n <- 101
#' df <- data.frame(x = seq(0, 1, length.out = n))
#' df$y <- 3 - 2*df$x + rnorm(n, 0, sd = .15)
#' 
#' with(df, plot(x, y))
#' 
#' (mod <- alm(~ x + y, data = df))
#' (y_coef <- coef(mod)[["y"]])
#' coef(mod) / y_coef
#' 
#' lm(y ~ x, data = df)
#' 
#' 
#' # example 2
#' df <- rvnorm(100, mp("x^2 + y^2 - 1"), sd = .1)
#' plot(df, asp = 1)
#' df <- as.data.frame(df)
#' mod <- alm(~ poly(x, 2, raw = TRUE) + poly(y, 2, raw = TRUE), data = df)
#' mod
#' 
#' 
  
  
  
  
  
  
  
  
  
#' @export
#' @rdname alm
alm <- function(form, data, ...){
    
  mat <- model.matrix(form, data)
  e <- eigen(crossprod(mat), symmetric = TRUE)
  
  ests <- structure(
    e$vectors[, which.min(e$values)],
    .Names = colnames(mat)
  )
  
  structure(
    .Data = list(
      coefficients = ests,
      eigen = e
    ),
    class = "alm"
  )
   
}




#' @export
#' @rdname alm
print.alm <- function(x, ...) {
  print(x$coefficients, ...)
}






