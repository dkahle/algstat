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
#' @param n The number of draws desired from each chain.
#' @param poly An mpoly object.
#' @param sd The "standard deviation" component of the normal kernel.
#' @param output \code{"simple"}, \code{"more"}, \code{"stanfit"}.
#' @param chains The number of chains to run for the random number generation,
#'   see [stan()].
#' @param cores The number of CPU cores to distribute the chains across, see
#'   [stan()].
#' @param warmup Number of warmup iterations in [stan()].
#' @param keep_warmup If \code{TRUE}, the MCMC warmup steps are included in the
#'   output.
#' @param normalized If \code{TRUE}, the polynomial is gradient-normalized. This
#'   is highly recommended.
#' @param inject_direct Directly specify printed polynomial to string inject
#'   into the stan code. Requires you specify \code{vars}, \code{numerator}, and
#'   \code{denominator}.
#' @param verbose \code{TRUE} or \code{FALSE}; determines level of messaging.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param numerator,denominator A character(1) containing the printed numerator
#'   of the variety normal distribution.
#' @param w A box window (-w,w) of the same dimension as the number of
#'   variables.
#' @param refresh The \code{refresh} argument of [stan()], which governs how
#'   much information is provided to the user while sampling.
#' @param ... Additional parameters to pass to [stan()].
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class [stanfit-class].
#' @author David Kahle
#' @examples
#'
#' \dontrun{ contains runs rstan
#'
#' ## basic usage
#' ########################################
#'
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(2000, p, sd = .1)
#' head(samps)
#'
#' samps <- rvnorm(2000, p, sd = .1, output = "tibble")
#' ggplot(samps, aes(x, y)) + geom_point() + coord_equal()
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
  chains = 4L, 
  cores = getOption("mc.cores", 1L), 
  warmup = floor(n/2), 
  keep_warmup = FALSE, 
  inject_direct = FALSE, 
  verbose = FALSE,
  w, 
  vars, 
  numerator, 
  denominator, 
  refresh,
  ...
) {
  
  # cran guard
  chain <- NULL; rm(chain)
  num <- NULL; rm(num)
  denom <- NULL; rm(denom)
  ng <- NULL; rm(ng)
  lp__ <- NULL; rm(lp__)
  iter <- NULL; rm(iter)
  if (missing(refresh)) if (verbose) refresh <- max(iter/10L, 1) else refresh <- 0L
  if(!missing(refresh)) stopifnot(is.numeric(refresh) && length(refresh) == 1L)
  
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
  if (verbose) cat(stan_code, "\n")

  # compile and run stan code
  fit <- rstan::stan(
    model_code = stan_code, 
    data = list("zero" = 0), 
    chains = chains, 
    iter = n, 
    warmup = warmup, 
    control = list("adapt_delta" = .999, "max_treedepth" = 20L),
    cores = cores,
    refresh = refresh,
    ...
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
    
  } else if (output == "simple") { 
    
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

