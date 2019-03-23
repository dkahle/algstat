#' Plot a variety
#'
#' Uses geom_contour() and ggplot() to plot an mpoly object representing a
#' variety in 2D space.
#'
#' @param mp an mpoly object
#' @param xlim vector representing x bounds
#' @param ylim vector representing y bounds
#' @param n number of mesh points in each dimension
#' @param nx number of mesh points in the abcissa (x)
#' @param ny number of mesh points in the ordinate (y)
#' @param ... additional parameters
#' @name ggvariety
#' @return A ggplot object containing variety plot
#' @author Phillip Hossu, Ryan Hebdon, Chong Sun, Grant Innerst, David Kahle
#' @examples
#'
#' ## basic usage
#' ##################################################
#'
#' ggvariety("y - x")
#' ggvariety("y - x^2")
#' ggvariety("x^2 + y^2 - 1")
#' ggvariety(c("x^2 + y^2 - 1", "y - x"))
#' 
#' 
#' 
#' ## setting limits
#' ##################################################
#' 
#' ggvariety("y - x^2")
#' ggvariety("y - x^2", xlim = c(-2,2), ylim = c(-2,2))
#' 
#' 
#' ## ggplot2 styling
#' ##################################################
#' 
#' library("ggplot2")
#' ggvariety("x^2 + y^2 - 1") + coord_equal()
#' ggvariety("x^2 + y^2 - 1") + coord_equal() + theme_bw()
#' ggvariety("x^2 + y^2 - 1") + coord_equal() + theme_classic()
#' ggvariety("x^2 + y^2 - 1") + coord_equal() + theme_void()
#' 
#' ggvariety(c("x^2 + y^2 - 1", "(x^2 + y^2)^3 - 4 x^2 y^2")) +
#'   coord_equal() + theme_void() +
#'   scale_color_manual(values = c("red", "blue"), guide = FALSE)
#'
#'
#'
#' ## possible issues
#' ##################################################
#' 
#' # at a low level, ggvariety() uses grDevices::contourLines()
#' # to numerically detect zero crossings. this is an imperfect process,
#' # so you may see gaps where none exist. as a general strategy, upping
#' # the number of sampled points on the grid is recommended.
#' # the below are commented to cut check time; they run
#' 
#' # ggvariety("y^2 - x^3 - x^2") + coord_equal()
#' # ggvariety("y^2 - x^3 - x^2", n = 201) + coord_equal()
#' # ggvariety(mp(c("x^2 + y^2 - 1", "y - x")))
#' 



#' @rdname ggvariety
#' @export
ggvariety <- function(
  mp, 
  xlim = c(-1, 1), 
  ylim = c(-1, 1), 
  n = 101, 
  nx = n, 
  ny = n, 
  ...
) {
  
  # define/wipe variables
  x <- NULL; rm(x)
  y <- NULL; rm(y)
  z <- NULL; rm(z)
  . <- NULL; rm(.)
  `p(x,y)` <- NULL; rm(`p(x,y)`)
  poly <- NULL; rm(poly)
  
  # check if mp argument was mpoly obj
  if (!is.mpoly(mp)) mp <- mpoly::mp(mp)
  p <- if (inherits(mp, "mpolyList")) length(mp) else 1L
  f <- as.function(mp, varorder = c("x", "y"), vector = FALSE, silent = TRUE)

  if (p == 1L) {
    
    tibble(
      "x" = seq(xlim[1], xlim[2], length.out = nx), 
      "y" = seq(ylim[1], ylim[2], length.out = ny)
    ) %>% 
      cross_df() %>% 
      mutate("z" = f(x, y)) %>% 
      ggplot(aes(x, y, "z" = z)) + 
        geom_contour(breaks = 0)
    
  } else { # several varieties
    
    # make xy grid
    list(
      "x" = seq(xlim[1], xlim[2], length.out = nx), 
      "y" = seq(ylim[1], ylim[2], length.out = ny)
    ) %>% 
      cross() ->
      xy_combos
    
    # evaluate polys at xy grid
    xy_combos %>% 
      map(~ f(.x$x, .x$y)) %>% 
      do.call(rbind, .) %>% 
      as.data.frame() %>% as_tibble() %>% 
      set_names(str_c("p", 1:p)) ->
      poly_vals_tb
    
    # combine grid and eval'd tibbles into one wide tibble, reshape
    xy_poly_tb <- xy_combos %>% 
      bind_rows() %>% 
      bind_cols(poly_vals_tb) %>% 
      gather("poly", "p(x,y)", -x, -y)
    
    # make plot
    xy_poly_tb %>% 
      ggplot(aes(x, y, z = `p(x,y)`, color = poly)) +
        geom_contour(breaks = 0)    
  }
  
}


# look at mesh, adaptively refine to add more points to 
# df near places where you think a crossing occurs



