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
#' @param n The number of draws desired from each chain after warmup.
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
#' @param thin [stan()] \code{thin} parameter.
#' @param normalized If \code{TRUE}, the polynomial is gradient-normalized. This
#'   is highly recommended.
#' @param inject_direct Directly specify printed polynomial to string inject
#'   into the stan code. Requires you specify \code{vars}, \code{numerator}, and
#'   \code{denominator}.
#' @param verbose \code{TRUE} or \code{FALSE}; determines level of messaging.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param numerator,denominator A character(1) containing the printed numerator
#'   of the variety normal distribution.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. A If a single number, a box window (-w,w) is applied to all
#'   variables.
#' @param refresh The \code{refresh} argument of [stan()], which governs how
#'   much information is provided to the user while sampling.
#' @param code_only If \code{TRUE}, will only formulate and return Stan code.
#' @param ... Additional parameters to pass to [stan()].
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class [stanfit-class].
#' @author David Kahle
#' @examples
#'
#' \dontrun{ runs rstan
#'
#' library("ggplot2")
#' rstan::rstan_options(auto_write = TRUE) # helps avoid recompiles
#'
#' ## basic usage
#' ########################################
#'
#' # single polynomial
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(2000, p, sd = .1)
#' head(samps)
#' str(samps) # 2000 * (4 chains)
#' plot(samps, asp = 1)
#'
#' # returning a data frame
#' (samps <- rvnorm(2000, p, sd = .1, output = "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y, color = `g[1]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#'
#' # more than one polynomial, # vars >= # eqns
#' p <- mp(c("x^2 + y^2 + z^2 - 1", "z"))
#' samps <- rvnorm(500, p, sd = .1, output = "tibble")
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[2]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, z, color = `g[2]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#'
#' # overdetermined system, vars < # eqns (not yet supported)
#' p <- mp(c("x", "y", "x + y"))
#' samps <- rvnorm(500, p, sd = .1, output = "tibble")
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[3]`)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ## using refresh to get more info
#' ########################################
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 100)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
#' rvnorm(2000, p, sd = .1, "tibble", refresh = -1)
#'
#'
#'
#' ## many chains in parallel
#' ########################################
#'
#' options(mc.cores = parallel::detectCores())
#' p <- mp("x^2 + (4 y)^2 - 1")
#' samps <- rvnorm(250, p, sd = .01, "tibble", verbose = TRUE, chains = 8)
#' ggplot(samps, aes(x, y)) + geom_point() + coord_equal()
#'
#'
#'
#' ## windowing for unbounded varieties
#' ########################################
#' # windowing is needed for unbounded varieties
#' # in the following, look at the parameters block
#'
#' p <- mp("x y - 1") # unbounded variety
#'
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
#'
#' rvnorm(1e3, p, sd = .01, "tibble", w = 1.15, code_only = TRUE)
#'
#' window <- list("x" = c(-1.5, 1.25), "y" = c(-2, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", w = window, code_only = TRUE)
#'
#' window <- list("x" = c(-1.5, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", w = window, code_only = TRUE)
#'
#'
#'
#' ## the importance of normalizing
#' ########################################
#' # one of the effects of the normalizing is to stabilize variances, making
#' # them roughly equivalent globally over the variety.
#'
#' # lemniscate of bernoulli
#' p <- mp("(x^2 + y^2)^2 - 2 (x^2 - y^2)")
#'
#' # normalized, good
#' (samps <- rvnorm(2000, p, .025, "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'
#' # unnormalized, bad
#' (samps <- rvnorm(2000, p, .025, "tibble", normalized = FALSE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'
#'
#'
#' ## semi-algebraic sets
#' ########################################
#' # inside the semialgebraic set x^2 + y^2 <= 1
#' # this is the same as x^2 + y^2 - 1 <= 0, so that
#' # x^2 + y^2 - 1 + s^2 == 0 for some slack variable s
#' # this is the projection of the sphere into the xy-plane.
#'
#' p <- mp("1 - (x^2 + y^2) - s^2")
#' samps <- rvnorm(1e4, p, sd = .01, "tibble", chains = 8)
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d() + coord_equal()
#'
#' ggplot(sample_n(samps, 2e3), aes(x, y, color = s)) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' # alternative representation
#' # x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) == -1
#' # so that s^2 (x^2 + y^2 - 1) + 1 == 0
#' # while this strategy works in theory, it doesn't work
#' # so well in practice, since s^2 is unbounded.
#' # it's gradient is also more complicated.
#'
#'
#' ## keeping the warmup / the importance of multiple chains
#' ########################################
#'
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' ggvariety(p, xlim = c(-3,3)) + coord_equal()
#'
#' samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, keep_warmup = TRUE, w = 5)
#' ggplot(samps, aes(x, y, color = iter)) +
#'   geom_point(size = 1, alpha = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#'
#'
#' ## ideal-variety correspondence considerations
#' ########################################
#'
#' p <- mp("x^2 + y^2 - 1")
#'
#' samps_1 <- rvnorm(250, p^1, sd = .1, output = "tibble", chains = 8)
#' samps_2 <- rvnorm(250, p^2, sd = .1, output = "tibble", chains = 8)
#' samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
#' samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
#' samps <- bind_rows(mget(apropos("samps_")))
#' samps$power <- rep(seq_along(apropos("samps_")), each = 2000)
#'
#' ggplot(samps, aes(x, y, color = num < 0)) +
#'   geom_point(size = .5) +
#'   coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
#'   facet_wrap(~ power)
#'
#'
#' ## neat examples
#' ########################################
#' # an implicit Lissajous region, view in separate window large
#'
#' # x = cos(m t + p)
#' # y = sin(n t + q)
#' (p <- lissajous(3, 2,  -pi/2, 0))
#' (p <- lissajous(4, 3,  -pi/2, 0))
#' (p <- lissajous(5, 4,  -pi/2, 0))
#' (p <- lissajous(3, 3,  0, 0))
#' (p <- lissajous(5, 5,  0, 0))
#' (p <- lissajous(7, 7,  0, 0))
#' ggvariety(p, n = 201) + coord_equal()
#'
#' p <- plug(p, "x", mp(".5 x"))
#' p <- plug(p, "y", mp(".5 y"))
#'
#' # algebraic set
#' samps <- rvnorm(5e3, p, sd = .01, "tibble", chains = 8, refresh = 100)
#' ggplot(samps, aes(x, y, color = factor(chain))) +
#'   geom_point(size = .5) + coord_equal()
#'
#' # semi-algebraic set
#' samps_normd <- rvnorm(5e3, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   normalized = TRUE, refresh = 100
#' )
#' samps_unormd <- rvnorm(5e3, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   normalized = FALSE, refresh = 100
#' )
#'
#' bind_rows(
#'   samps_normd %>% mutate(normd = TRUE),
#'   samps_unormd %>% mutate(normd = FALSE)
#' ) %>%
#'   ggplot(aes(x, y)) +
#'     geom_point() +
#'     facet_wrap(~ normd) +
#'     coord_equal()
#'
#' bind_rows(
#'   samps_normd %>% mutate(normd = TRUE),
#'   samps_unormd %>% mutate(normd = FALSE)
#' ) %>%
#'   ggplot(aes(x, y)) +
#'     geom_bin2d(bins = 100) +
#'     facet_wrap(~ normd) +
#'     coord_equal()
#'
#'
#' samps_normd %>%
#'   ggplot(aes(x, y)) +
#'     geom_bin2d(bins = 200) +
#'     facet_wrap(~ factor(chain)) +
#'     coord_equal()
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
  chains = 4L, 
  warmup = floor(n/2), 
  keep_warmup = FALSE, 
  thin = 1L,
  inject_direct = FALSE, 
  verbose = FALSE,
  cores = getOption("mc.cores", 1L),
  normalized = TRUE,
  w, 
  vars, 
  numerator, 
  denominator, 
  refresh,
  code_only = FALSE,
  ...
) {
  
  # cran guard
  chain <- NULL; rm(chain)
  num <- NULL; rm(num)
  denom <- NULL; rm(denom)
  ng <- NULL; rm(ng)
  lp__ <- NULL; rm(lp__)
  iter <- NULL; rm(iter)
  . <- NULL; rm(.)
  
  if (!is.mpolyList(poly)) poly <- structure(list(poly), class = "mpolyList")
  
  vars <- mpoly::vars(poly)
  
  n_eqs  <- length(poly)
  n_vars <- length(vars)
  # if (n_eqs > n_vars) stop("Overdetermined systems not yet supported.")
  
  if (missing(refresh)) if (verbose) refresh <- max(n/10L, 1) else refresh <- 0L
  if (!missing(refresh)) stopifnot(is.numeric(refresh) && length(refresh) == 1L)
  
  
  
  ## create stan code
  ########################################
  
  needs_to_compile_new_model <- TRUE
  
  if (needs_to_compile_new_model) {
    
    printed_polys <- mpoly:::print.mpolyList(poly, silent = TRUE, stars = TRUE, plus_pad = 0, times_pad = 0) %>% str_replace_all("\\*\\*", "^")
    printed_polys <- str_c(printed_polys, collapse = ", ")
    
    d <- get("deriv.mpoly", asNamespace("mpoly"))
    p <- get("print.mpoly", asNamespace("mpoly"))
      
    printed_jac <- array("", dim = c(n_eqs, n_vars))
    for (i in 1:n_eqs) {
      for (j in 1:n_vars) {
        if (normalized) {
          printed_jac[i,j] <- p(d(poly[[i]], vars[j]), silent = TRUE, stars = TRUE, plus_pad = 0L, times_pad = 0L)
        } else {
          printed_jac[i,j] <- if (i == j) "1" else "0"
        }
      }
    }
    printed_jac <- printed_jac %>% 
      apply(1L, str_c, collapse = ", ") %>% 
      str_c("      [", ., "]", collapse = ", \n") %>% 
      str_c("[\n", ., "\n    ]") %>% 
      str_replace_all("\\*\\*", "^")
    
    
    # set variables
    if (missing(w)) {
      parms <- glue::glue("real {vars};") 
    } else {
      if (is.numeric(w) && length(w) == 1L) {
        parms <- glue::glue("real<lower=-{w},upper={w}> {vars};") 
      } else if (is.list(w)) {
        stopifnot(all(names(w) %in% vars))
        parms <- vector("character", n_vars)
        for (var_ndx in seq_along(vars)) {
          parms[var_ndx] <- if (vars[var_ndx] %in% names(w)) {
            var_ndx_in_w <- which(names(w) == vars[var_ndx])
            glue::glue("real<lower={w[[var_ndx_in_w]][1]},upper={w[[var_ndx_in_w]][2]}> {vars[var_ndx]};") 
          } else {
            glue::glue("real {vars[var_ndx]};") 
          }
        }
      } else {
        stop("bound parameter misspecified, see ?rvnorm.", call. = FALSE)
      }
    }
    parms <- parms %>% str_c(collapse = "\n    ") 

    # variance expression
    if (is.numeric(sd) && is.vector(sd) && length(sd) == 1L) {
      Si_exp <- if (n_vars >= n_eqs) {
        glue::glue("{sd^2} * tcrossprod(J)") # = "{sd^2} * J * J'"
      } else {
        glue::glue("{sd^2} * (tcrossprod(J) + diag_matrix(rep_vector(1, {n_eqs})))")
      }
    } else if (is.numeric(sd) && is.vector(sd)) {
      stop("This sd not yet supported.")
    } else if (is.numeric(sd) && is.matrix(sd)) {
      stop("This sd not yet supported.")
      Si_exp <- glue::glue("quad_form_sym(sd, J)") # = J' * sd * J
    }
    
    stan_code <- glue::glue("
      parameters {
        {{parms}}
      } 
      
      transformed parameters {
        vector[{{n_eqs}}] g = [{{printed_polys}}]';
        matrix[{{n_eqs}},{{n_vars}}] J = {{printed_jac}};
        matrix[{{n_eqs}},{{n_eqs}}] Si = {{Si_exp}};
      } 
      
      model {
        target += multi_normal_lpdf(g | rep_vector(0, {{n_eqs}}), Si);
      }
    ",  .open = "{{", .close = "}}"
    )
    if (verbose) cat(stan_code, "\n")  
    if (code_only) return(stan_code)

    # compile code
    if (!verbose) message("Compiling model... ", appendLF = FALSE)
    model <- rstan::stan_model("model_code" = stan_code, ...)
    if (!verbose) message("done.")
    
  }
  
  
  ## sample from the distribution
  ######################################## 
  
  fit <- rstan::sampling(
    "object" = model,
    # "data" = data, # will be list of coefficients?
    "chains" = chains,
    "iter" = n + warmup,
    "warmup" = warmup,
    "thin" = thin,
    "refresh" = refresh,
    "control" = list("adapt_delta" = .999, "max_treedepth" = 20L),
    "cores" = cores,
    ...
  )
  
  
  ## return, parse output as desired 
  ########################################
  
  if (output == "stanfit") return(fit)
    
  if (output == "tibble") {
    
    samps <- fit %>%
      rstan::extract(permuted = FALSE, inc_warmup = keep_warmup) %>% 
      purrr::array_branch(2L) %>% 
      purrr::imap(
        ~ tibble::as_tibble(.x) %>% mutate(chain = .y)
      ) %>% 
      bind_rows() %>% 
      mutate(
        iter = rep((as.integer(!keep_warmup)*warmup+1):(n+warmup), chains),
        chain = as.integer(str_sub(chain, 7L))
      )
    
    return(samps)
    
  } 
  
  if (output == "simple") { 

    samps <- fit %>% 
      rstan::extract(permuted = FALSE, inc_warmup = keep_warmup) %>% 
      purrr::array_branch(2L) %>% 
      purrr::imap(
        ~ tibble::as_tibble(.x) %>% mutate(chain = .y)
      ) %>% 
      bind_rows() %>% 
      mutate(
        iter = rep((as.integer(!keep_warmup)*warmup+1):(n+warmup), chains),
        chain = as.integer(str_sub(chain, 7L))
      ) %>% 
      # dplyr::select(-num, -denom, -ng, -lp__, -chain, -iter) %>% 
      .[1:n_vars] %>% 
      as.matrix()
    
    return(samps)
    
  }
  
}

