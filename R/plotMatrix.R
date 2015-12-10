#' Plot a matrix
#' 
#' plotMatrix is a R variant of Matlab's spy function.
#' 
#' @param A a matrix
#' @return a ggplot object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @export 
#' @examples
#' 
#' # the no-three-way interaction configuration
#' A <- hmat(c(3,3,3), subsets(1:3, 2))
#' plotMatrix(A)
#' plotMatrix(markov(A))
#' 
plotMatrix <- function(A){
  low <- if(any(A < 0)){
    low <- "blue"; high <- "red"
    fillScale <- scale_fill_gradient2(
      low = "blue", mid = "grey80", 
      high = "red", midpoint = 0, guide = FALSE, space = "Lab"
    )
  } else {
    fillScale <- scale_fill_gradient(
      low = "white", high = "black", guide = FALSE
    )
  }
  df <- expand.grid(x = 1:ncol(A), y = 1:nrow(A))
  B <- t(A)
  df$A <- as.integer(B[,ncol(B):1])
  qplot(x, y, data = df, fill = A, geom = "tile") +
    fillScale +
    theme_bw() + coord_equal() + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(
      line = element_blank(), text = element_blank(),
      plot.margin = grid::unit(c(.5, .5, -.5, -.5), "lines")
    )
}