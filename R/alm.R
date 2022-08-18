#' Algebraic least squares regression
#'
#' Algebraic least squares regression
#'
#' @param formula A formula specifying the basis monomials, passed to
#'   [stats::model.matrix]
#' @param data The dataset on which to apply the method, passed to
#'   [stats::model.matrix]
#' @param tol Tolerance in determining whether the design matrix is full rank
#'   (imposed on eigenvalues)
#' @param lambda A ridge regression penalty
#' @param ... Additional parameters, presently discarded
#' @param x parameter to pass to [print()]
#' @return A numeric vector of model parameter estimates
#' @name alm
#' @examples
#'
#'
#' ## fitting a linear formula
#' ########################################
#'
#' n <- 101
#' df <- data.frame(x = seq(0, 1, length.out = n))
#' df$y <- 3 - 2*df$x + rnorm(n, 0, sd = .15)
#'
#' plot(df)
#'
#' (mod <- alm(~ x + y, data = df))
#' sum(coef(mod)^2) # euclidean norm = 1
#'
#' X <- model.matrix(~ x + y, data = df)
#' XtX <- crossprod(X) # = t(X) %*% X
#' (e <- eigen(XtX) )
#' # eigen(crossprod(X) + .01*diag(3)) # ridge penalty
#' coef(mod) /  coef(mod)[["y"]]
#'
#'
#' (best_fit_poly <- as.mpoly(mod) )
#' 
#'
#' # compare to lm
#' lm(y ~ x, data = df)
#' 
#' (ests <- structure(
#'   qr.Q(qr(XtX))[, which.min(e$values)], 
#'   .Names = colnames(X)
#' ))
#' ests / ests[["y"]]
#' 
#' X <- cbind(1, df$x)
#' y <- df$y
#' solve(crossprod(X), t(X) %*% y)
#' 
#'
#'
#' ## fitting a quadratic formula
#' ########################################
#'
#' df <- rvnorm(100, mp("x^2 + (2 y)^2 - 1"), sd = .05, output = "tibble")
#'
#'
#' (p <- mp("1 + x + y + x^2 + x y + y^2"))
#' str(p)
#' lapply(monomials(p), print, silent = TRUE)
#'
#' (p_monos <- monomials(p))
#' (mono_names <- sapply(monomials(p), print, silent = TRUE))
#' fs <- as.function(p_monos, varorder = c("x", "y"))
#' head(t(apply(df, 1, fs)))
#' fit <- alm(~ t(apply(df, 1, fs)) - 1, data = df)
#' betas <- coef(fit)
#' names(betas) <- mono_names
#'
#' library("ggplot2")
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
#' (poly <- mp(poly))
#' ggvariety(poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   geom_point(aes(x, y), data = df, alpha = .1, inherit.aes = FALSE) +
#'   coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))
#'
#'
#'
#' ## ridge penalty
#' ########################################
#'
#' # n = 2
#' (df <- data.frame(x = 0:1, y = 0:1))
#' lm(y ~ x, data = df)
#' lm(y ~ x + I(x^2), data = df)
#' (X <- model.matrix(y ~ x + I(x^2), data = df))
#' (XtX <- crossprod(X))
#' eigen(XtX) # 2 nonzero eigenvalues
#' 
#' (mod <- alm(~ x + y, data = df))
#' as.mpoly(mod)
#' alm(~ x + y, data = df)$eigen
#' alm(~ x + y, data = df, lambda = .01)
#' alm(~ x + I(x^2) + y, data = df)
#' alm(~ x + I(x^2) + y, data = df, lambda = .01)
#' 
  
  
  
  
  
  
  
  
  
#' @export
#' @rdname alm
alm <- function(formula, data, tol = 1e-10, lambda, ...){
    
  
  X <- model.matrix(formula, data)
  XtX <- crossprod(X)
  if (!missing(lambda)) XtX <- XtX + lambda*diag(nrow(XtX))
  e <- eigen(XtX, symmetric = TRUE)
  
  if (any(e$values <= tol)) {
    warning(
      "  Very small eigenvalues detected;", 
      "\n  the design matrix is numerically singular.",
      call. = FALSE
    )
  }
  
  ests <- structure(
    e$vectors[, which.min(e$values)],
    .Names = colnames(X)
  )
  
  structure(
    .Data = list(
      coefficients = ests,
      eigen = e,
      X = X,
      XtX = XtX
    ),
    class = "alm"
  )
   
}




#' @export
#' @rdname alm
print.alm <- function(x, ...) {
  print(x$coefficients, ...)
}






