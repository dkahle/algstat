#' Sample from the fiber of a contingency table
#'
#' Sample from the fiber of a contingency table
#'
#' @param n the number of observations
#' @param A the configuration matrix of the model defining the fiber
#' @param b integer vector; the vector of sufficient statistics
#' @param x integer vector; a vector in the fiber to be used as an
#'   initialization point for the "walk" method
#' @param method "sis", "walk", or "hybrid"
#' @param dist "uniform" or "hypergeometric"
#' @param parallel parallelize the sampling?
#' @param includeRejections should rejections be returned? (note: rejection
#'   tables aren't completed in the SIS procedure.)
#' @param format format of the returned moves, "mat", "vec", or "tab"
#' @param dim the dimensions of the table if "tab" is used, see [vec2tab()]
#' @param thin thinning argument (if using walk method)
#' @param ... ...
#' @return a named numeric vector
#' @author David Kahle \email{david@@kahle.io}, Ruriko Yoshida
#'   \email{ruriko.yoshida@@uky.edu}
#' @export rfiber
#' @examples
#'
#'
#' A <- hmat(c(2,2), 1:2)
#' b <- rep.int(4, 4)
#'
#' rfiber(10, A = A, b = b)
#' rfiber(10, A = A, b = b, format = "vec")
#' rfiber(10, A = A, b = b, format = "tab", dim = c(2, 2))
#'
#' set.seed(1)
#' (tab <- rfiber(1, A = A, b = b, format = "tab", dim = c(2, 2))[[1]])
#' x <- tab2vec(tab)
#'
#' \dontrun{ # requires 4ti2
#'
#' rfiber(100, A = A, x = x, method = "walk")
#'
#' library(microbenchmark)
#' microbenchmark(
#'   rfiber(100, A = A, b = b),
#'   rfiber(100, A = A, x = x, method = "walk")
#' )
#'
#' # the distribution of the samples is seemingly near-uniform
#' tabs <- rfiber(1e4, A, b, format = "vec", parallel = TRUE)
#' table(sapply(tabs, paste, collapse = " "))
#'
#'
#'
#' tabs <- rfiber(1000, A = A, b = b, format = "vec")
#' unique(tabs)
#' lapply(unique(tabs), vec2tab, dim = c(2, 2))
#'
#' table(sapply(tabs, paste, collapse = " " ))
#'
#'
#'
#'
#' data(politics); politics
#' tab2vec(politics)
#' b <- A %*% tab2vec(politics)
#' tabs <- rfiber(1000, A = A, b = b)
#' unique(t(tabs))
#' table(apply(tabs, 2, paste, collapse = " " )) # roughly uniform
#'
#'
#'
#'
#'
#' # poisson regression example
#' J <- 5
#' A <- rbind(1L, 1:J)
#' b <- c(5, 15)
#'
#'
#'
#'
#' system.time(rfiber(1e4, A = A, b = b))
#' system.time(rfiber(1e4, A = A, b = b, parallel = TRUE))
#'
#' }
#' 
rfiber <- function(n, ..., A, b, x, 
  method = c("sis", "walk", "hybrid"),
  dist = c("uniform", "hypergeometric"),
  parallel = FALSE, includeRejections = FALSE, thin = 1,
  format = c("mat", "vec", "tab"), dim
){
   
  ## arg checking
  stopifnot(is.wholenumber(n) && n > 0)
  stopifnot(all(is.wholenumber(A)))
  if(!missing(b)){
    stopifnot(all(is.wholenumber(b)))
    stopifnot(nrow(A) == length(b))
  }
  if(!missing(x)){
    stopifnot(all(is.wholenumber(x)))
    stopifnot(ncol(A) == length(x))
  }
  method <- match.arg(method)
  dist   <- match.arg(dist)
  format <- match.arg(format)
  if(includeRejections && format == "mat"){
    stop("format mat is not available when including rejections, change to \"vec\" or \"tab\".", call. = FALSE)
  }
  
  
  ## 
  if(method == "sis"){
    if(!missing(x) && missing(b)) b <- A %*% tab2vec(x)
    samps <- rfiber_sis(n, A = A, b = b, parallel = parallel, includeRejections = includeRejections, ...)
    if(format == "vec"){
      return(samps)
    } else if(format == "mat"){
      return(simplify2array(samps))
    } else if(format == "tab"){
      return(lapply(samps, vec2tab, dim = dim))
    }    
  } else if(method == "walk"){
    samps <- rfiber_walk(n = n, A = A, x = x, thin = thin, ...)
    if(format == "mat"){
      return(samps)
    } else if(format == "vec"){
      return(FALSE)
    } else if(format == "tab"){
      return(FALSE)
    }
  } else if(method == "hybrid"){
    return(rfiber_hybrid(n, A, b, x))
  }
  
}






rfiber_sis <- function(n, A, b, parallel, cluster, includeRejections = FALSE, ...){
    
  ## incorporate accepted-only samples
  if(includeRejections){
    sampleOne <- function() rfiberOne(A, b)$x
  } else {
    sampleOne <- function() {
      rejectedSample <- TRUE
      while(rejectedSample){
        out <- rfiberOne(A, b)
        rejectedSample <- out$reject
      }
      out$x
    }
  }
  
  ## if not run in parallel
  if(missing(parallel) || !parallel) return(replicate(n, sampleOne(), FALSE)  )
  
  ## compute in parallel  
  if(is.mac() || is.unix()){  ## parallel mac/unix
    
    # max out clusters if mc.cores not set (override default of 2)
    if(is.null(getOption("mc.cores"))) options(mc.cores = detectCores())
    
    # compute and return
    return( mclapply(1:n, function(discardThisArgument) sampleOne()) )    
    
  } else if(is.win()){        ## parallel windows
    
    # make cluster if none
    if(missing(cluster)){    
      message("parallel: distributing load across ", detectCores(), " clusters.")
      algstatCluster <- makeCluster(getOption("mc.cores", detectCores()))
      on.exit(stopCluster(algstatCluster), add = TRUE)
    } else {
      algstatCluster <- cluster
    }
    
    # compute and return
    return( parLapply(algstatCluster, 1:n, function(discardThisArgument) sampleOne()) )    
    
  } else {
    
    warning("parallel is not implemented for this system.")
    return(replicate(n, rfiberOne(A, b)$x, FALSE, ...)  )
    
  }  
  
}







rfiber_walk <- function(n, A, x, moves = zbasis(A), thin = 0, ...){
  walk(x, cbind(moves, -moves), n, thin)
}

rfiber_hybrid <- function(n, A, b, x, ...){
  stop("method not yet implemented.")
}





