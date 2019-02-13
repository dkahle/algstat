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
#' ggvariety(mp("y - x"))
#' ggvariety(mp("y - x^2"))
#' ggvariety(mp("x^2 + y^2 - 1"))
#' 
#' 
#' ## setting limits
#' ##################################################
#' 
#' ggvariety(mp("y - x^2"))
#' ggvariety(mp("y - x^2"), xlim = c(-2,2), ylim = c(-2,2))
#' 
#' 
#' ## ggplot2 styling
#' ##################################################
#' 
#' library("ggplot2")
#' ggvariety(mp("x^2 + y^2 - 1")) + coord_equal()
#' ggvariety(mp("x^2 + y^2 - 1")) + coord_equal() + theme_bw()
#' ggvariety(mp("x^2 + y^2 - 1")) + coord_equal() + theme_void()
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
#' ggvariety("y^2 - x^3 - x^2") + coord_equal()
#' ggvariety("y^2 - x^3 - x^2", n = 201) + coord_equal()
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
  
  # check if mp argument was mpoly obj
  if (!is.mpoly(mp)) mp <- mpoly::mp(mp)
  f <- as.function(mp, varorder = c("x", "y"), vector = FALSE, silent = TRUE)
  
  # make tibble / ggplot
  tibble(
    "x" = seq(xlim[1], xlim[2], length.out = nx), 
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  ) %>% 
    cross_df() %>% 
    mutate("z" = f(x, y)) %>% 
    ggplot(aes(x, y, "z" = z)) + 
      geom_contour(breaks = 0)
}


# look at mesh, adaptively refine to add more points to 
# df near places where you think a crossing occurs



