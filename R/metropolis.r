#' The Metropolis Algorithm
#'
#' Given a starting table (as a vector) and a collection of moves, run the
#' Metropolis-Hastings algorithm starting with the starting table.
#'
#' See Algorithm 1.1.13 in LAS, the reference below.
#'
#' @param init the initial step
#' @param moves the moves to be used (the negatives will be added); they are
#'   arranged as the columns of a matrix.
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param dist steady-state distribution; "hypergeometric" (default) or
#'   "uniform"
#' @param engine \code{"C++"} or \code{"R"}? (C++ is significantly faster)
#' @name metropolis
#' @return a list
#' @export metropolis
#' @author David Kahle
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures
#'   on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @examples
#'
#' ## basic use
#' ############################################################
#'
#' # move up and down integer points on the line y = 100 - x
#' # sampling from the hypergeometric distribution
#' # note: negative moves are added internally
#' init <- c(10L, 90L)
#' moves <- matrix(c(1,-1), ncol = 1)
#'
#' # it helps running each of these lines several times to get a feel for things
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1)
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1, dist = "uniform")
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1, engine = "R")
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1, engine = "R", dist = "uniform")
#'
#' # a bigger simulation
#' iter <- 1e3
#' out <- metropolis(init, moves, iter = iter, burn = 0, thin = 1)
#' hist(out$steps[1,], breaks = 100)
#'
#' # view convergence through trace plot
#' plot(1:iter, out$steps[1,1:iter])
#'
#' # sampling from the hypergeometric distribution
#' out <- metropolis(init, moves, iter = iter, dist = "uniform")
#' hist(out$steps[1,], breaks = 100)
#' plot(1:iter, out$steps[1,1:iter])
#'
#' # look at autocorrelation
#' acf(out$steps[1,])
#'
#'
#' ## thinning to reduce autocorrelation with the thin argument
#' ############################################################
#'
#' out <- metropolis(init, moves, iter = iter, dist = "uniform", thin = 2500)
#' acf(out$steps[1,])
#' hist(out$steps[1,], breaks = 100)
#'
#'
#'
#' ## burn in with the burn argument
#' ############################################################
#'
#' set.seed(1L)
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1, engine = "R")
#'
#' set.seed(1L)
#' metropolis(init, moves, iter = 10, burn = 0, thin = 1, engine = "R")
#'
#' set.seed(1L)
#' metropolis(init, moves, iter =  5, burn = 5, thin = 1, engine = "R")
#'
#' # this can't be shown as easily with the C++ engine because of seed issues
#'
#'
#' 
metropolis <- function(
  init, 
  moves, 
  iter = 1e3L, 
  burn = 0L, 
  thin = 1L,
  dist = c("hypergeometric","uniform"), 
  engine = c("C++","R")
){

  ## preliminary checking
  ##################################################
  dist <- match.arg(dist)
  if (is.character(engine) && length(character) == 1L && engine == "Cpp") engine <- "C++"
  engine <- match.arg(engine)
  
  if (thin == 0L) {
    message("thin = 1 corresponds to no thinning, resetting thin = 1.")
    thin <- 1L
  }
  

  ## in R
  ##################################################
  if(engine == "R"){

    n_moves <- ncol(moves)
    state   <- matrix(nrow = nrow(moves), ncol = iter)
  
    ## burn in
    #########################
  
    current <- unname(init)
    burn_unifs <- runif(burn)
    main_unifs <- runif(iter*thin)  
    
    message("Running chain (R)... ", appendLF = FALSE)
    
    if (burn > 0) {
      
      for(k in 1L:burn){
  
        move      <- sample(c(-1,1), 1) * moves[,sample(n_moves,1)]
        prop_state <- current + move
      
        if(any(prop_state < 0)){
          prob <- 0
        } else {
          if(dist == "hypergeometric"){
            prob <- exp( sum(lfactorial(current)) - sum(lfactorial(prop_state)) )
          } else { # dist == "uniform"
            prob <- 1
          }
        }
      
        if(burn_unifs[k] < prob) current <- prop_state # else current
      
      }
      
      state[,1L] <- current
      
    } else { # no burn-in
      
      state[,1L] <- current
      
    }
  
    ## main sampler
    #########################
  
    totalRuns <- 0L
    prob_total <- 0
    
    for(k in 2:iter){
    	
    	for(j in 1:thin){
  
        move       <- sample(c(-1L,1L), 1L) * moves[,sample(n_moves,1L)]
        prop_state <- current + move
      
        if(any(prop_state < 0)){
          prob <- 0
        } else {
          if(dist == "hypergeometric"){
            prob <- exp( sum(lfactorial(current)) - sum(lfactorial(prop_state)) )
          } else { # dist == "uniform"
            prob <- 1
          }
        }
        prob_total <- prob_total + min(1, prob)
  
        if(main_unifs[k*(thin-1)+j] < prob) current <- prop_state # else current
        
        totalRuns <- totalRuns + 1L     
      }
  
      state[,k] <- current    
    }
    message("done.")  
    
    ## format output
    out <- list(
      steps = state, 
      moves = moves, 
      accept_prob = prob_total / totalRuns
    )
  

  }
  
  ## in Cpp
  ##################################################
  if (engine == "C++") {
    
    current   <- unname(init)  
    allMoves  <- cbind(moves, -moves)  
    sampler   <- if (dist == "hypergeometric") {
      metropolis_hypergeometric_cpp
    } else {
      metropolis_uniform_cpp
    }
    
    message("Running chain (C++)... ", appendLF = FALSE)  
    if (burn > 0) current <- sampler(current, allMoves, burn, 1)$steps[,burn]
    out       <- sampler(current, allMoves, iter, thin)
    out$moves <- moves
    message("done.")

  }


  ## return output
  ##################################################  

  out[c("steps", "moves", "accept_prob")]
}








#' @rdname metropolis
#' @export
rawMetropolis <- function(init, moves, iter = 1E3, dist = "hypergeometric"){
  metropolis(init, moves, iter, burn = 0, thin = 1, dist = dist) 
}

