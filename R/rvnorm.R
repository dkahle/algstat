#' The Variety Normal Distribution
#'
#' Un-normalized density function and random generation for the variety normal
#' distribution with mean equal to \code{poly} and "standard deviation" equal to
#' \code{sd}. Please see details for caveats.
#' 
#' If the variety you are interested in is connected, this strategy should work 
#' well out of the box.  If it isn't, you'll likely need to rely on running 
#' multiple chains, and it is very likely, if not probable, that the sampling 
#' will be biased to one or more of those components and down-sample others.
#' Question: what is the relative likelihood of each component, or an equal unit 
#' of length, on different components? How does this generalize to more 
#' varieties of varying dimensions?
#'
#' @param n the number of draws desired from each chain
#' @param poly an mpoly object
#' @param sd the "standard deviation"
#' @param chains the number of chains to run for theh generation
#' @param warmup the warmup
#' @param cores number of cores to use
#' @param keep_warmup discard warmup?
#' @param output \code{"simple"}, \code{"more"}, \code{"stanfit"}
#' @param normalized should the polynomial be gradient-normalized?
#' @param inject_direct directly specify printed polynomial to string inject
#'   into the stan code. requires you specify vars, numerator, and denominator
#' @param verbose message user?
#' @param vars a character vector of the indeterminates in the distribution
#' @param numerator,denominator a character(1) containing the printed numerator
#'   of the variety normal distribution
#' @param w a box window (-w,w) of the same dimension as the number of variables
#' @param ... additional parameters
#' @name rvnorm
#' @return A matrix whose rows are the individual draws from the distribution
#' @author David Kahle
#' @examples
#' 
#' \dontrun{ contains runs rstan
#'
#' library("ggplot2"); theme_set(theme_minimal())
#' library("dplyr")
#'
#' ## 0d variety in 1d
#' ########################################
#' 
#' p <- mp("x")
#' samps <- rvnorm(2000, p, sd = 1)
#' head(samps, breaks = "scott")
#' hist(samps, breaks = seq(-5, 5, .1))
#' ggplot(as.data.frame(samps), aes(x)) + geom_histogram(bins = 100)
#' qqnorm(samps); curve(1*x, col = "red", add = TRUE)
#' shapiro.test(sample(samps, 5000L))
#' 
#' p <- mp("(x-3) x (x+3)")
#' (samps <- rvnorm(2000, p, sd = 1, output = "tibble"))
#' ggplot(samps, aes(x)) + geom_histogram(bins = 100)
#' ggplot(samps, aes(x = iter, y = x)) + geom_line(aes(color = factor(chain)))
#' # poor mixing
#' 
#' 
#' ## 1d variety in 2d
#' ########################################
#' 
#' # circle
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(2000, p, sd = .01)
#' head(samps)
#' plot(samps, asp = 1)
#' ggplot(as.data.frame(samps), aes(x, y)) + geom_point() + coord_equal()
#'   
#' # this illustrates that the sampling is clearly efficiently moving around the
#' # support of the distribution
#' (samps <- rvnorm(2000, p, sd = .05, output = "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point() + coord_equal()
#' ggplot(samps, aes(x, y, color = iter)) +
#'   geom_point(size = .5) + geom_path(alpha = .2) + 
#'   coord_equal() + facet_wrap(~ factor(chain))
#'   
#'   
#'   
#'   
#'   
#' # two circles
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' (samps <- rvnorm(2000, p, .01, "tibble", keep_warmup = TRUE))
#' ggplot(samps, aes(x, y)) + geom_point() + coord_equal()
#' ggplot(samps, aes(x, y, color = factor(chain))) + geom_point() + coord_equal()
#' 
#' # it's only through the multiple chains that we see both components
#' # if we initialized one point one each component, we could essentially
#' # guarantee we've seen them all, as opposed to this, which is more probabilistic
#' ggplot(samps, aes(x, y, color = iter)) + 
#'   geom_point(size = .5) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#' 
#' (samps <- rvnorm(2000, p, .2, "tibble", keep_warmup = TRUE))
#' ggplot(samps, aes(x, y, color = iter)) + 
#'   geom_point(size = .5) + geom_path(alpha = .2) +  
#'   coord_equal() + facet_wrap(~ factor(chain))
#' 
#' 
#' 
#' # alpha curve
#' p <- mp("y^2 - (x^3 + x^2)")
#' (samps <- rvnorm(2000, p, .025, "tibble", w = 1.15))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' 
#' 
#' 
#' # heart
#' p <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
#' (samps <- rvnorm(2000, p, .01, "tibble", keep_warmup = TRUE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'   
#' ggplot(samps, aes(x, y, color = iter)) + 
#'   geom_point(size = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#' 
#' 
#' 
#' 
#' # iva 4e pp10-11, positions of "hand" in two-rod system.
#' # first rod is of length 2, second of length 1
#' # this provides random configurations of the system
#' # (remember: the conditions are not met exactly, so both bars are only 
#' # approximately the lengths they should be.)
#' p <- mp("(x^2 + y^2 - 4)^2 + ((x - z)^2 + (y - w)^2 - 1)^2")
#' (samps <- rvnorm(2000, p, .01, "tibble"))
#' ggplot(samps, aes(z, w)) + geom_point(size = .5) + coord_equal()
#'   
#' samps %>% 
#'   sample_n(25) %>% mutate(config = factor(1:25)) %>% 
#'   ggplot() +
#'     geom_hline(yintercept = 0, color = "gray50") +
#'     geom_vline(xintercept = 0, color = "gray50") +
#'     geom_point(aes(0, 0)) +
#'     geom_segment(aes(x = 0, y = 0, xend = x, yend = y)) +
#'     geom_point(aes(x, y)) +
#'     geom_segment(aes(x = x, y = y, xend = z, yend = w)) +
#'     geom_point(aes(z, w)) +
#'     coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
#'     facet_wrap(~ config)
#'   
#'   
#'   
#' # encircled clover - normalized vs unnormalized
#' p <- mp("(x^2 + y^2 - 1) ((x^2 + y^2)^3 - 4 x^2 y^2)")
#' samps_normalized <- rvnorm(2000, p, .01, "tibble")
#' ggplot(samps_normalized, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' 
#' ggplot(samps_normalized, aes(x, y, color = iter)) + 
#'   geom_point(size = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#'  
#' samps_unnormalized <- rvnorm(2000, p, .01, "tibble", normalized = FALSE)
#' ggplot(samps_unnormalized, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' 
#' samps <- bind_rows(
#'   mutate(samps_normalized, strategy = "normalized"),
#'   mutate(samps_unnormalized, strategy = "unnormalized")
#' )
#' ggplot(samps, aes(x, y)) + 
#'   geom_point(size = .5) + coord_equal() +
#'   facet_wrap(~ strategy)
#' 
#' 
#' 
#' # jon's curve
#' p <- mp(
#'   "0.6575807899*x^6 - 0.8710143351*x^5 + 1.97274237*x^4*y^2 + 
#'   0.5151773876*x^4*y - 2.19037755*x^4 - 1.74202867*x^3*y^2 - 
#'   1.101761115*x^3*y + 1.697857343*x^3 + 1.97274237*x^2*y^4 + 
#'   1.030354775*x^2*y^3 - 4.066060954*x^2*y^2 - 0.2963687857*x^2*y + 
#'   2.251360555*x^2 - 0.8710143351*x*y^4 - 1.101761115*x*y^3 + 
#'   0.1583742363*x*y^2 + 1.654420569*x*y + 0.04312151611*x + 
#'   0.6575807899*y^6 + 0.5151773876*y^5 - 1.875683404*y^4 - 
#'   1.179380775*y^3 + 1.342115975*y^2 + 0.3931953974*y - 0.1446830451", 
#'   stars_only = TRUE
#' )
#' (samps <- rvnorm(2000, p, .01, "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' 
#' ggplot(samps, aes(x, y, color = iter)) + 
#'   geom_point(size = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#' 
#' 
#' 
#' # elliptic curve with a = -1 and b = 1
#' # https://en.wikipedia.org/wiki/Elliptic_curve
#' p <- mp("y^2 - (x^3 - x + 1)")
#' (samps <- rvnorm(2000, p, .01, "tibble", w = 3))
#' ggplot(samps, aes(x, y)) + geom_point() + coord_equal()
#'   
#' ggplot(samps, aes(x, y, color = iter)) +
#'   geom_point() + geom_path(alpha = .2) + coord_equal() +
#'   facet_wrap(~ factor(chain))
#' 
#' 
#' 
#' 
#' ## surface examples
#' ########################################
#' 
#' library("plotly")
#' 
#' 
#' ## 1d in 3d
#' #########################
#' 
#' # twisted cubic
#' p <- mp("(y - x^2)^2 + (z - x^3)^2")
#' (samps <- rvnorm(2000, p, .01, "tibble", keep_warmup = TRUE, w = 5))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 2, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' 
#' ## 2d in 3d
#' #########################
#' 
#' # sphere
#' p <- mp("x^2 + y^2 + z^2 - 1")
#' (samps <- rvnorm(2000, p, .01, "tibble", keep_warmup = TRUE))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' # torus
#' p <- mp("(x^2 + y^2 + z^2 + 2^2 - 1^2)^2 - 4 2^2 (x^2 + y^2)")
#' (samps <- rvnorm(2000, p, .01, "tibble"))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' # double torus
#' (p <- mp("-0.01 + 4 x^2 - 20 x^3 + 41 x^4 - 44 x^5 + 26 x^6 - 8 x^7 + x^8 - 
#'          4 x y^2 + 10 x^2 y^2 - 8 x^3 y^2 + 2 x^4 y^2 + y^4 + z^2"))
#' (samps <- rvnorm(2000, p, .01, "tibble", keep_warmup = TRUE))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' 
#' # heart
#' p <- mp("(x^2 + 2.25 y^2 + z^2 - 1)^3 - x^2 z^3 - .1125 x^2 z^3")
#' samps <- rvnorm(2000, p, .01, "tibble")
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' 
#' 
#' 
#' # geisha
#' p <- mp("x^2 y z + x^2 z^2 - y^3 z - y^3")
#' (samps <- rvnorm(2000, p, .01, "tibble", w = 5))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' 
#' # whitney
#' p <- mp("x^2 - y^2 z")
#' (samps <- rvnorm(2000, p, .05, "tibble", w = 5))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "lines+markers",
#'   line = list(width = 1.5), marker = list(size = 3),
#'   split = ~factor(chain), opacity = .2
#' )
#' 
#' 
#' 
#' 
#' # columpius
#' p <- mp("x^3 y + x z^3 + y^3 z + z^3 + 7 z^2 + 5 z")
#' (samps <- rvnorm(2000, p, .05, "tibble", w = 5, warmup = 200))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#' 
#' 
#' 
#' 
#' # something else
#' p <- mp("x^2 - y^2 z^2 + z^3")
#' (samps <- rvnorm(4000, p, .01, "tibble", w = 5, warmpup = 100))
#' 
#' plot_ly(
#'   samps, x = ~x, y = ~y, z = ~z, 
#'   type = "scatter3d", mode = "markers",
#'   marker = list(size = 1, color = "black")
#' )
#' 
#'
#' }
#' 













#' @rdname rvnorm
#' @export
rvnorm <- function(
  n, 
  poly, 
  sd, 
  output = "simple", 
  normalized = TRUE,
  chains = parallel::detectCores(), 
  warmup = 1000, 
  cores = parallel::detectCores(), 
  keep_warmup = FALSE, 
  inject_direct = FALSE, 
  verbose = FALSE,
  w, 
  vars, 
  numerator, 
  denominator, 
  ...
) {
  
  # cran guard
  chain <- NULL; rm(chain)
  num <- NULL; rm(num)
  denom <- NULL; rm(denom)
  ng <- NULL; rm(ng)
  lp__ <- NULL; rm(lp__)
  iter <- NULL; rm(iter)
  
  # check args
  if (inject_direct) {
    
    if(missing(vars) || missing(numerator) || missing(denominator)) 
      stop("if inject_direct = TRUE, you must specify vars, numerator, and denominator")
    
  } else {
    
    if (!is.mpoly(poly)) poly <- mp(poly)
    poly <- mpoly:::reorder.mpoly(poly, varorder = sort(mpoly::vars(poly)))
    vars <- mpoly::vars(poly)

    numerator <- print(poly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L) %>% str_replace_all("[*]{2}", "^")
    if (normalized) {
      denominator <- if (length(vars) > 1) Reduce(`+`, gradient(poly)^2) else gradient(poly)^2
      denominator <- print(denominator, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L) %>% str_replace_all("[*]{2}", "^")
      denominator <- glue::glue("sqrt({denominator})")
    } else {
      denominator <- "1"
    }
    
  }
  
  # make stan code
  if (missing(w)) {
    parms <- glue::glue("real {vars};") %>% str_c(collapse = "\n  ") 
  } else {
    parms <- glue::glue("real<lower=-{w},upper={w}> {vars};") %>% str_c(collapse = "\n  ") 
  }

  stan_code <- glue::glue("
parameters {
  {{parms}}
} 

transformed parameters {
  real num = {{numerator}};
  real denom = {{denominator}};
  real ng = num / denom;
} 

model {
  target += normal_lpdf(ng | 0.00, {{sd}});
}
    ",  .open = "{{", .close = "}}"
  )
  if(verbose) cat(stan_code, "\n")
 
  # compile and run stan code
  fit <- rstan::stan(
    model_code = stan_code, 
    data = list("zero" = 0), 
    chains = chains, 
    iter = n, 
    warmup = warmup, 
    control = list("adapt_delta" = .999, "max_treedepth" = 20L),
    cores = cores
  )

  # return
  if (output == "stanfit") {
    
    return(fit)
    
  } else if(output == "tibble") {
    
    fit %>% 
      rstan::extract(permuted = FALSE, inc_warmup = keep_warmup) %>% 
      purrr::array_branch(2L) %>% 
      purrr::imap(
        ~ tibble::as_data_frame(.x) %>% mutate(chain = .y)
      ) %>% 
      bind_rows() %>% 
      mutate(
        iter = rep((as.integer(!keep_warmup)*warmup+1):n, chains),
        chain = as.integer(str_sub(chain, 7L))
      )
    
  } else {
    
    fit %>% 
      rstan::extract(permuted = FALSE, inc_warmup = keep_warmup) %>% 
      purrr::array_branch(2L) %>% 
      purrr::imap(
        ~ tibble::as_data_frame(.x) %>% mutate(chain = .y)
      ) %>% 
      bind_rows() %>% 
      mutate(
        iter = rep((as.integer(!keep_warmup)*warmup+1):n, chains),
        chain = as.integer(str_sub(chain, 7L))
      ) %>% 
      dplyr::select(-num, -denom, -ng, -lp__, -chain, -iter) %>% 
      as.matrix()
    
  }
  
}

