#' Projection onto a variety
#'
#' A R-based implementation of the gradient descent homotopies
#'
#'
#' @param x0 Atomic vector, the point to be projected.
#' @param poly An mpoly object, typically created with [mp()].
#' @param dt The t-mesh size for the homotopy.
#' @param varorder A character vector specifying the variable order to pass to
#'   [mpoly::as.function.mpoly()].
#' @param n_correct The number of Newton correction iterations to use.
#' @param al A numeric vector of length 2; the patch to do projective
#'   calculations over.
#' @return A numeric vector the same length as \code{x0}.
#' @references Griffin, Z. and J. Hauenstein (2015). Real solutions to systems
#'   of polynomial equations and parameter continuation. \emph{Advances in
#'   Geometry} 15(2), pp.173--187.
#' @author David Kahle
#' @examples
#'
#'
#' ## basic usage
#' ########################################
#'
#' x0 <- c(1,1)
#' p <- mp("x^2 + y^2 - 1")
#' (x0_proj <- project_onto_variety(x0, p))
#' 
#' as.function(p)(x0_proj)
#' sqrt(2)/2
#' 
#' cbind(t(x0), t(x0_proj)) %>% 
#'   as.data.frame() %>% as_tibble() %>% 
#'   purrr::set_names(c("x", "y", "x_proj", "y_proj")) -> df
#' 
#' ggvariety(p) + coord_equal() +
#'   ggplot2::geom_segment(
#'     aes(x, y, xend = x_proj, yend = y_proj), 
#'     data = df, inherit.aes = FALSE
#'   )
#' 
#' 
#' 
#' ## options
#' ########################################
#' 
#' x0 <- c(1,1)
#' p <- mp("x^2 + y^2 - 1")
#' project_onto_variety(x0, p, message = TRUE)
#' project_onto_variety(x0, p, dt = .25, message = TRUE)
#' 
#' 
#' 
#' ## more complex example
#' ########################################
#' 
#' p <- mp("(x^2 + y^2)^2 - 2 (x^2 - y^2)")
#' ggvariety(p, c(-2, 2), n = 201) + coord_equal()
#' 
#' x0 <- c(.05, .10)
#' (x0_proj <- project_onto_variety(x0, p))
#' 
#' cbind(t(x0), t(x0_proj)) %>% 
#'   as.data.frame() %>% as_tibble() %>% 
#'   purrr::set_names(c("x", "y", "x_proj", "y_proj")) -> df
#' 
#' ggvariety(p, c(-2, 2)) + coord_equal() +
#'   ggplot2::geom_segment(
#'     aes(x, y, xend = x_proj, yend = y_proj), 
#'     data = df, inherit.aes = FALSE
#'   )
#' 
#' 
#' 
#' 
#' ## projecting a dataset - grid
#' ########################################
#' 
#' library("ggplot2")
#' library("dplyr")
#' 
#' (p <- lissajous(3, 3, 0, 0))
#' ggvariety(p, n = 251) + coord_equal()
#' 
#' set.seed(1)
#' (s <- seq(-1, 1, .25))
#' n <- length(s)
#' grid <- expand.grid(x = s, y = s)
#' grid$x <- jitter(grid$x)
#' grid$y <- jitter(grid$y)
#' 
#' ggplot(grid, aes(x, y)) + geom_point() + coord_equal()
#' 
#' grid %>% as.matrix() %>% 
#'   apply(1, function(x0) project_onto_variety(x0, p)) %>% t() %>% 
#'   as.data.frame() %>% as_tibble() %>% 
#'   purrr::set_names(c("x_proj", "y_proj")) ->
#'   grid_proj
#'   
#' ggvariety(p, n = 251) + coord_equal() +
#'   geom_segment(
#'     aes(x, y, xend = x_proj, yend = y_proj), 
#'     data = bind_cols(grid, grid_proj), inherit.aes = FALSE
#'   ) +
#'   geom_point(aes(x, y), data = grid, inherit.aes = FALSE)
#'   
#'   
#' # here's what happens when you use a naive implementation
#' f <- as.function(p, varorder = c("x", "y"))
#' naive_project_to_variety <- function(x0) {
#'   optim(x0, function(.) f(.)^2, method = "BFGS")$par
#' }
#' 
#' grid %>% 
#'   select(x, y) %>% as.matrix() %>% 
#'   apply(1, naive_project_to_variety) %>% t() %>% 
#'   as.data.frame() %>% as_tibble() %>% 
#'   purrr::set_names(c("x_proj", "y_proj")) %>% as_tibble() ->
#'   grid_proj2
#'   
#' df <- bind_rows(
#'   bind_cols(grid, grid_proj) %>% mutate(method = "homotopy"),
#'   bind_cols(grid, grid_proj2) %>% mutate(method = "naive")
#' )
#'   
#' ggvariety(p, n = 251) +
#'   geom_segment(
#'     aes(x, y, xend = x_proj, yend = y_proj), 
#'     data = df, inherit.aes = FALSE
#'   ) +
#'   geom_point(aes(x, y), data = grid, inherit.aes = FALSE) +
#'   coord_equal() +
#'   facet_grid(. ~ method)
#' 
#' 
#' ## projecting a dataset - rvnorm
#' ########################################
#' 
#' library("ggplot2")
#' library("dplyr")
#' 
#' \dontrun{ requires stan
#' 
#' (p <- lissajous(3, 3, 0, 0))
#' ggvariety(p, n = 251) + coord_equal()
#' 
#' set.seed(1)
#' (samps <- rvnorm(1e4, p, sd = .025, output = "tibble"))
#' 
#' ggplot(samps, aes(x, y)) + 
#'   geom_point(aes(color = chain)) + 
#'   coord_equal() + 
#'   facet_wrap(~ chain)
#'   
#' ggplot(samps, aes(x, y)) + 
#'   geom_bin2d(binwidth = .03*c(1,1)) + 
#'   coord_equal()
#'   
#' # cut down on draws for time
#' subsamps <- samps %>% sample_n(500)
#' ggplot(subsamps, aes(x, y)) + geom_point() + coord_equal()
#' 
#' subsamps %>% 
#'   select(x, y) %>% 
#'   as.matrix() %>% 
#'   apply(1, function(x0) project_onto_variety(x0, p)) %>% t() %>% 
#'   as.data.frame() %>% as_tibble() %>% 
#'   purrr::set_names(c("x_proj", "y_proj")) %>% 
#'   bind_cols(subsamps, .) ->
#'   subsamps
#'   
#' ggvariety(p, n = 251) + coord_equal() +
#'   geom_segment(
#'     aes(x, y, xend = x_proj, yend = y_proj), 
#'     data = subsamps, inherit.aes = FALSE
#'   ) +
#'   geom_point(
#'     aes(x, y, color = factor(chain)), 
#'     data = subsamps, inherit.aes = FALSE
#'   )
#' 
#' }
#'



#' @export
project_onto_variety <- function(
  x0, poly, dt = .05, varorder = vars(poly), 
  n_correct = 2, al = rnorm(length(x0)), 
  message = FALSE, tol = sqrt(.Machine$double.eps)
) {
  
  n_vars <- length(varorder)
  
  ts <- seq(1, 0, -dt)
  vn <- c(x0, al[1], 0)
  
  # polynomial as a function
  g <- poly
  gfunc <- as.function(g, varorder = varorder, silent = TRUE)
  
  # jacobian
  dg <- deriv(g, var = varorder)
  dgfunc <- as.function(dg, varorder = varorder, silent = TRUE)
  
  # hessian
  ddg <- lapply(dg, deriv, var = varorder)
  ddgfunc_list <- lapply(ddg, as.function, varorder = varorder, silent = TRUE)
  ddgfunc <- function(x) sapply(ddgfunc_list, function(f) f(x))
  
  
  Ha <- function(v, t) {
    # v = (x, y, la0, la1)
    x <- v[1:2]; la <- v[3:4]
    c(
      gfunc(x) - t*gfunc(x0),
      as.numeric(cbind((x-x0), dgfunc(x)) %*% la),
      la[1] + la[2]*al[2] - al[1]
    )
  }
  # Ha(vn, .99)
  
  JHa <- function(v, t) {
    x <- v[1:n_vars]
    la0 <- v[n_vars+1] 
    la1 <- v[n_vars+2] 
    al1 <- al[2]
    rbind(
      c(dgfunc(x), 0, 0),
      cbind(
        la0*diag(n_vars) + la1*ddgfunc(x),
        x - x0,
        dgfunc(x)
      ),
      c(rep(0, n_vars), 1, al1)
    )
  }
  # JHa(vn, .99)
  
  Ht <- function(v, t) c(-gfunc(x0), 0, 0, 0)
  
  
  for (i in 2:length(ts)) {
    
    # predict
    vn <- vn + as.numeric(solve(JHa(vn, t = ts[i-1])) %*% Ht(vn, t = ts[i-1])) * dt
    if (message) message(paste(round(vn, 5), collapse = " "))
    
    # correct
    for (. in 1:n_correct) {
      vnp1 <- as.numeric( vn - solve(JHa(vn, t = ts[i])) %*% Ha(vn, t = ts[i]) )
      vn <- vnp1
      if (message) message("  ", paste(round(vn, 5), collapse = " "))
    }
    
  }
  
  resid <- gfunc(vn[1:n_vars])
  if (abs(resid) >= tol) warning(sprintf("Tolerance not met (residual = %f).", resid))
  
  vn[1:n_vars]
}







