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
#'
#' ## fitting a linear form
#' ########################################
#' 
#' n <- 101
#' df <- data.frame(x = seq(0, 1, length.out = n))
#' df$y <- 3 - 2*df$x + rnorm(n, 0, sd = .15)
#' 
#' plot(df)
#' 
#' (mod <- alm(~ x + y, data = df))
#' (y_coef <- coef(mod)[["y"]])
#' coef(mod) / y_coef
#' 
#' lm(y ~ x, data = df)
#' 
#' 
#' ## fitting a quadratic form
#' ########################################
#' 
#' df <- rvnorm(100, mp("x^2 + (2 y)^2 - 1"), sd = .05, output = "tibble")
#' ggplot(df, aes(x, y)) + geom_point() + coord_equal()
#' p <- function(x, deg) {
#'   out <- poly(x, deg, raw = TRUE, simple = TRUE)
#'   names <- paste0(sprintf("%s^", deparse(substitute(x))), 1:deg)
#'   names <- gsub("\\^1", "", names)
#'   colnames(out) <- names
#'   out
#' }
#' mod <- alm(~ p(x,2) + p(y,2) + I(x*y), data = df)
#' mod
#' str(mod)
#' 
#' poly <- paste(
#'   coef(mod), 
#'   c("1", "x", "x^2", "y", "y^2", "x y"), 
#'   collapse = " + "
#' )
#' poly <- mp(poly)
#' ggvariety(poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   geom_point(aes(x, y), data = df, alpha = .1, inherit.aes = FALSE) +
#'   coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))
#'   
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






