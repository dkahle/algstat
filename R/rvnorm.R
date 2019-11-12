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
#' library("tidyverse")
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
#' 
#' ggplot(samps, aes(x, y, color = g)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'   
#' ggplot(samps, aes(x, y)) + 
#'   stat_density2d(
#'     aes(fill = stat(density)), 
#'     geom = "raster", contour = FALSE
#'    ) + 
#'   coord_equal()
#'
#'
#'
#' # more than one polynomial, # vars > # eqns, underdetermined system
#' p <- mp(c("x^2 + y^2 + z^2 - 1", "z"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[2]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, z, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#'
#'
#' # more than one polynomial, # vars < # eqns, overdetermined system
#' p <- mp(c("3 x", "3 y", "2 x + 2 y", "3 (x^2 + y)", "3 (x^2 - y)"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#' 
#' samps %>% 
#'   select(x, y, starts_with("g")) %>% 
#'   pivot_longer(starts_with("g"), "equation", "value") %>% 
#'   ggplot(aes(x, y, color = value)) + geom_point() +
#'     scale_color_gradient2(mid = "gray80") + coord_equal() +
#'     facet_wrap(~ equation)
#'     
#'
#' ## using refresh to get more info
#' ########################################
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
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
#' (samps <- rvnorm(1e4, p, sd = .01, "tibble", verbose = TRUE, chains = 8))
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .01*c(1,1)) + coord_equal()
#' # decrease sd to get more uniform sampling
#'
#'
#' ## windowing for unbounded varieties
#' ########################################
#' # windowing is needed for unbounded varieties
#' # in the following, look at the parameters block
#'
#' p <- mp("x y - 1") # unbounded variety, 1 poly
#' p <- mp(c("x y - 1", "y - x")) # 2 polys
#'
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = 1.15)
#'
#' window <- list("x" = c(-1.5, 1.25), "y" = c(-2, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = window)
#'
#' window <- list("x" = c(-1.5, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = window)
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
#' (samps <- rvnorm(2000, p, .05, "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' # unnormalized, bad
#' (samps <- rvnorm(2000, p, .05, "tibble", normalized = FALSE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
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
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' ggplot(sample_n(samps, 2e3), aes(x, y, color = s)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' # alternative representation
#' # x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) + 1 == 0
#' # note that it's gradient is more complicated.
#' p <- mp("s^2 (x^2 + y^2 - 1) + 1")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, w = 2, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#'
#' ## keeping the warmup / the importance of multiple chains
#' ########################################
#'
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' ggvariety(p, xlim = c(-3,3)) + coord_equal()
#'
#' # notice the migration of chains initialized away from the distribution
#' # (it helps to make the graphic large on your screen)
#' samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, keep_warmup = TRUE)
#' ggplot(samps, aes(x, y, color = iter)) +
#'   geom_point(size = 1, alpha = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#'   
#' samps <- rvnorm(2500, p, sd = .05, "tibble", chains = 8, keep_warmup = TRUE)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + 
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
#' # samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
#' # samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
#' samps <- bind_rows(mget(apropos("samps_[1-4]")))
#' samps$power <- rep(seq_along(apropos("samps_[1-4]")), each = 2000)
#'
#' ggplot(samps, aes(x, y, color = g < 0)) +
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
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .02*c(1,1)) + coord_equal()
#'   
#' ggplot(samps, aes(x, y, color = factor(chain))) +
#'   geom_point(size = .5) + coord_equal() +
#'   facet_wrap(~ factor(chain))
#'
#' # semi-algebraic set
#' samps_normd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   normalized = TRUE, refresh = 100
#' )
#' samps_unormd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   normalized = FALSE, refresh = 100
#' )
#'
#' bind_rows(
#'   samps_normd  %>% mutate(normd = TRUE),
#'   samps_unormd %>% mutate(normd = FALSE)
#' ) %>%
#'   ggplot(aes(x, y)) +
#'     geom_bin2d(binwidth = .05*c(1,1)) + 
#'     facet_grid(normd ~ chain) +
#'     coord_equal()
#'
#'
#' ggplot(samps_normd, aes(x, y)) + 
#'   geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
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
  g <- NULL; rm(g)
  ndg <- NULL; rm(ndg)
  ng <- NULL; rm(ng)
  lp__ <- NULL; rm(lp__)
  iter <- NULL; rm(iter)
  . <- NULL; rm(.)
  
  
  if (is.mpoly(poly)) {
    n_eqs <- 1L
  } else if (is.character(poly) || is.mpolyList(poly)) {
    n_eqs <- length(poly)
  } else {
    stop("`poly` should be either a character vector, mpoly, or mpolyList.", call. = FALSE)
  }

  if (missing(refresh)) if (verbose) refresh <- max(ceiling(n/10), 1L) else refresh <- 0L
  if (!missing(refresh)) stopifnot(is.numeric(refresh), length(refresh) == 1L)  
  
  
  mpoly_to_stan <- function (mpoly) {
    p <- get("print.mpoly", asNamespace("mpoly"))
    p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L) %>% 
      str_replace_all("[*]{2}", "^")
  }
  
  mpolyList_to_stan <- function (mpolyList) {
    p <- get("print.mpolyList", asNamespace("mpoly"))
    p(mpolyList, silent = TRUE, stars = TRUE, plus_pad = 0, times_pad = 0) %>% 
      str_replace_all("\\*\\*", "^") %>% 
      str_c(collapse = ", ")
  }
  
  d <- get("deriv.mpoly", asNamespace("mpoly"))
  

  # create stan code --------------------------------------------------------

  needs_to_compile_model <- TRUE
  
  if (needs_to_compile_model) {
    
    
    if (n_eqs == 1L) {
    # single polynomial provided
      
      
      
      if (!is.mpoly(poly)) poly <- mp(poly)
      if (missing(vars)) vars <- mpoly::vars(poly)
      n_vars <- length(vars)
      reorder.mpoly <- get("reorder.mpoly", asNamespace("mpoly"))
      poly <- reorder.mpoly(poly, varorder = sort(vars))
      
      g_string <- mpoly_to_stan(poly)
      
      if (normalized) {
        if (n_vars > 1) {
          grad <- deriv(poly, var = mpoly::vars(poly))
          ndg_sq <- Reduce(`+`, grad^2) 
        } else {
          ndg_sq <- gradient(poly)^2
        }
        ndg_sq_string <- mpoly_to_stan(ndg_sq)
        ndg_string <- glue::glue("sqrt({ndg_sq_string})")
      } else {
        ndg_string <- "1"
      }
        
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
      parms <- parms %>% str_c(collapse = "\n  ") 
      
      stan_code <- glue::glue("
data {
  real<lower=0> si;
}
      
parameters {
  {{parms}}
} 
      
transformed parameters {
  real g = {{g_string}};
  real ndg = {{ndg_string}};
} 
      
model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
      ",  .open = "{{", .close = "}}"
      )
      if (verbose) cat(stan_code, "\n")  
      if (code_only) return(stan_code)

      # compile code
      if (!verbose) message("Compiling model... ", appendLF = FALSE)
      model <- rstan::stan_model("model_code" = stan_code, ...)
      if (!verbose) message("done.")
      
      
    
    } else {
    # multiple polynomials provided
      
      vars <- mpoly::vars(poly)
      n_vars <- length(vars)
      # if (n_eqs > n_vars) stop("Overdetermined systems not yet supported.")
      
    
      printed_polys <- mpolyList_to_stan(poly)
        
      printed_jac <- array("", dim = c(n_eqs, n_vars))
      for (i in 1:n_eqs) {
        for (j in 1:n_vars) {
          if (normalized) {
            printed_jac[i,j] <- mpoly_to_stan(d(poly[[i]], vars[j]))
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
      parms <- parms %>% str_c(collapse = "\n  ") 
  
      # normalization matrix
      if (is.numeric(sd) && is.vector(sd) && length(sd) == 1L) {
        
        gbar_string <- if (n_vars == n_eqs) {
          "J \\ g"
        } else if (n_vars > n_eqs) {
          "J' * ((J*J') \\ g)"
        } else {
          "(J'*J) \\ (J'*g)"
        }
        
      } else stop("This sd not yet supported.") 
       
      
      # write stan code
      stan_code <- glue::glue("
data {
  real<lower=0> si;
}
      
parameters {
  {{parms}}
} 
        
transformed parameters {
  vector[{{n_eqs}}] g = [{{printed_polys}}]';
  matrix[{{n_eqs}},{{n_vars}}] J = {{printed_jac}};
} 
        
model {
  target += normal_lpdf(0.00 | {{gbar_string}}, si);
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
  }
  
  

  # run stan sampler --------------------------------------------------------
  
  fit <- rstan::sampling(
    "object" = model,
    "data" = list("si" = sd),
    "chains" = chains,
    "iter" = n + warmup,
    "warmup" = warmup,
    "thin" = thin,
    "refresh" = refresh,
    "control" = list("adapt_delta" = .999, "max_treedepth" = 20L),
    "cores" = cores,
    ...
  )
  
  
  

  # parse output and return -------------------------------------------------

  if (output == "stanfit") return(fit)
  
  convert_mat_to_tibble_add_chain <- function(mat, chain_no) {
    tibble::as_tibble(mat) %>% mutate(chain = chain_no)
  }
    
  if (output == "tibble") {
    
    samps <- fit %>%
      rstan::extract(permuted = FALSE, inc_warmup = keep_warmup) %>% 
      purrr::array_branch(2L) %>% 
      purrr::imap(convert_mat_to_tibble_add_chain) %>% 
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
      purrr::imap(convert_mat_to_tibble_add_chain) %>% 
      bind_rows() %>% 
      mutate(
        iter = rep((as.integer(!keep_warmup)*warmup+1):(n+warmup), chains),
        chain = as.integer(str_sub(chain, 7L))
      ) %>% 
      dplyr::select(-g, -ndg, -lp__, -chain, -iter) %>%
      # .[1:n_vars] %>% 
      as.matrix()
    
    return(samps)
    
  }
  
}






# 
# 
# 
# stanfit_to_tibble <- function(stanfit, inc_warmup) {
#   
# }
# 
# 








