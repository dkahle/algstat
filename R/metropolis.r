#' The Metropolis Algorithm
#' 
#' Given a starting table (as a vector) and a collection of moves, 
#' run the Metropolis-Hastings algorithm starting with the starting 
#' table.
#' 
#' See Algorithm 1.1.13 in LAS, the reference below.
#' 
#' @param init the initial step
#' @param moves the moves to be used (the negatives will be added); 
#'   they are arranged as the columns of a matrix.
#' @param suffStats the sufficient statistics of the model. Only used when SIS = TRUE. Defaulted to 0.
#' @param config the configuration matrix that encodes the model. Only used when SIS = TRUE. Defaulted to matrix(0).
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param dist steady-state distribution; "hypergeometric" (default)
#'   or "uniform"
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @param hitAndRun Whether or not to use the discrete hit and run algorithm in
#'   the metropolis algorithm. Defaulted to FALSE
#' @param SIS If TRUE, with a small probability the move will be chosen randomly from the uniform distribution 
#'  on the fiber using Sequential Importance "Like" Sampling methods. Defaulted to FALSE
#' @param nonUniform If TRUE, moves will be chosen adaptively using a move weighting system that uses information 
#' from previous steps. Defaulted to FALSE
#' @param adaptive Option when hitAndRun = TRUE. If adaptive = TRUE, hit and run will choose a proposal distribution adaptively. 
#' Defaulted to FALSE
#' 
#' @name metropolis
#' @return a list
#' @export metropolis
#' @author David Kahle
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). 
#'   \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser 
#'   Verlag AG.
#' @examples
#' 
#' \dontrun{
#' 
#' library(ggplot2); theme_set(theme_bw())
#' 
#' # move up and down integer points on the line y = 100 - x
#' # sampling from the hypergeometric distribution
#' init <- c(10,90)
#' moves <- matrix(c(1,-1), ncol = 1)
#' out <- metropolis(init, moves)
#' qplot(out$steps[1,])
#' 
#' # view convergence through trace plot
#' qplot(1:1000, out$steps[1,])
#' 
#' # sampling from the uniform distribution
#' out <- metropolis(init, moves, dist = "uniform")
#' qplot(out$steps[1,])
#' 
#' # view convergence through trace plot
#' qplot(1:1000, out$steps[1,])
#' 
#' # look at autocorrelation
#' acf(out$steps[1,])
#' # thin
#' out <- metropolis(init, moves, dist = "uniform", thin = 2500)
#' acf(out$steps[1,])
#' qplot(out$steps[1,])
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' data(handy)
#' 
#' exp   <- loglin(handy, as.list(1:2), fit = TRUE)$fit
#' e <- unname(tab2vec(exp))
#' h <- t(t(unname(tab2vec(handy))))
#' chisq <- algstat:::computeX2sCpp(h, e)
#' 
#' out <- loglinear(~ Gender + Handedness, data = handy)
#' chisqs <- algstat:::computeX2sCpp(out$steps, e)
#' 
#' mean(chisqs >= chisq)
#' fisher.test(handy)$p.value
#' 
#' 
#' 
#' 
#' 
#' A <- hmat(c(2,2), as.list(1:2))
#' moves <- markov(A)
#' outC <- metropolis(tab2vec(handy), moves, suffStats = tab2vec(handy) %*% A, config = A, 1e4, engine = "Cpp")
#' str(outC)
#' outR <- metropolis(tab2vec(handy), moves, suffStats = tab2vec(handy) %*% A, config = A, 1e4, engine = "R", thin = 20)
#' str(outR)
#' 
#' # showSteps(out$steps)
#' 
#' 
#' library(microbenchmark)
#' microbenchmark(
#'   metropolis(tab2vec(handy), moves, suffStats = tab2vec(handy) %*% A, config = A,engine = "Cpp"),
#'   metropolis(tab2vec(handy), moves, suffStats = tab2vec(handy) %*% A, config = A,engine = "R")
#' )
#' 
#' # cpp ~ 20-25x faster
#' 
#' 
#' 
#' 
#' # examples using the extra options inside metropolis 
#' 
#' 
#' data("HairEyeColor")
#' tbl <- tab2vec(apply(HairEyeColor, c(1, 2), sum))
#' A <- hmat(c(4,4),1:2)
#' moves <- markov(A)
#' suffStats <- A %*% tbl
#' 
#' # base metropolis algorithm 
#' base <- metropolis(tbl, moves, suffStats, A)
#' 
#' # hit and run option
#' har <- metropolis(tbl, moves, suffStats, A, hitAndRun = TRUE)
#' 
#' # check convergence through trace plots
#' baseStats <- algstat:::computeUProbsCpp(base$steps)
#' harStats <- algstat:::computeUProbsCpp(har$steps)
#' 
#' data <- data.frame(baseStats = baseStats, harStats = harStats, steps = 1:1000)
#' 
#' ggplot(data = data) + geom_line(aes(steps, baseStats)) + 
#' geom_line(aes(steps, harStats), color = "red") + 
#' labs(x = "Steps", y = "UNLL value", title = "Base Algorithm vs. Algorithm with Hit and Run option in red")
#' 
#' 
#' showSteps <- function(steps){
#'   apply(steps, 2, function(x){
#'     x <- format(x)
#'     tab <- vec2tab(x, dim(handy))
#'     message(
#'       paste(
#'         apply(tab, 1, paste, collapse = " "),
#'         collapse = " "
#'       )
#'     )
#'     message("
#' ", appendLF = F)
#'   })
#'   invisible()
#' }
#' # showSteps(out$steps)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' }
#' 
#' 
metropolis <- function(init, moves, suffStats = 0, config = matrix(0), iter = 1E3, burn = 0, thin = 1,
                       dist = c("hypergeometric","uniform"), engine = c("Cpp","R"), 
                       hitAndRun = FALSE, SIS = FALSE, nonUniform = FALSE, adaptive = FALSE
){
  
  ## preliminary checking
  ##################################################
  dist <- match.arg(dist)
  engine <- match.arg(engine)
  if(thin == 0){
    message("thin = 1 corresponds to no thinning, resetting thin = 1.")
    thin <- 1
  }
  
  
  ## in R
  ##################################################
  if(engine == "R"){
    
    nMoves <- ncol(moves)
    state  <- matrix(nrow = nrow(moves), ncol = iter)
    ## run burn-in
    
    current <- unname(init)
    
    message("Running chain (R)... ", appendLF = FALSE)
   
    #Setting up non-uniform move sampling framework
    if(nonUniform == TRUE){
      moveDist <- rep(1,nMoves)
      counter <- nMoves
    }
    
    unifs <- runif(burn)
    if(nonUniform == TRUE){
       moveProb <- runif(burn)
    }
    if(burn > 0) {
      for(k in 1:burn){
        #Hit and Run option
        if(hitAndRun)
        {
          move <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
          workMove <- move[move != 0]
          workCurrent <- current[move != 0]
          workMoves <- -1 * workCurrent / workMove
          lowerBound <- if(any(workMoves < 0)){max(subset(workMoves,subset = workMoves < 0))}else{1}
          upperBound <- if(any(workMoves > 0)){min(subset(workMoves,subset = workMoves > 0))}else{-1}
          
          if(any(workMoves == 0)){
          workPropStatelow <- current + lowerBound * move
          workPropStateup <-  current + upperBound * move
          if(any(workPropStatelow < 0)){
            lowerBound <- 1
          }
          if(any(workPropStateup < 0)){
            upperBound <- -1
          }
          }
          
          multiple  <- sample(lowerBound:upperBound,1)
          if(multiple  == 0){
            multiple  <- 1
          }
          propState <- current + multiple  * move
        }
        
        #Non-uniform move sampling option
        if(nonUniform == TRUE)
        {
          
          for(l in 1:nMoves){
            if(moveProb[k] <= sum(moveDist[1:l])/counter){
              move <- moves[,l]
              whichMove <- l
            }
          }
          move <- sample(c(-1,1), 1) * move
          propState <- current + move
        }
        
        else{
          move      <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
          propState <- current + move
        }
       
        
         if(any(propState < 0)){
          prob <- 0
        } else {
          if(dist == "hypergeometric"){
            prob <- exp( sum(lfactorial(current)) - sum(lfactorial(propState)) )
          } else { # dist == "uniform"
            prob <- 1
          }
        }
        
        if(nonUniform == TRUE){
          if(unifs[k] < prob){
            current  <- propState
            moveDist[whichMove] <- moveDist[whichMove] + 1
            counter <- counter + 1
          }
        }else{
          if(unifs[k] < prob) current <- propState # else current
        }
        
      }
      state[,1] <- current
    }
    
    ## run main sampler
    
    totalRuns <- 0
    probTotal <- 0
    unifs <- runif(iter*thin)
    
    for(k in 2:iter){
      
      for(j in 1:thin){
        
        if(hitAndRun)
        {
          move <- moves[,sample(nMoves,1)]
          workMove <- move[move != 0]
          workCurrent <- current[move != 0]
          workMoves <- (-1 * workCurrent) / workMove
          lowerBound <- if(any(workMoves < 0)){max(subset(workMoves,subset = workMoves < 0))}else{1}
          upperBound <- if(any(workMoves > 0)){min(subset(workMoves,subset = workMoves > 0))}else{-1} 
          
          if(any(workMoves == 0)){
            workPropStatelow <- current + lowerBound * move
            workPropStateup <-  current + upperBound * move
            if(any(workPropStatelow < 0)){
              lowerBound <- 1
            }
            if(any(workPropStateup < 0)){
              upperBound <- -1
            }
          }
         
          multiple  <- sample(lowerBound:upperBound,1)
          if(multiple  == 0){
            multiple  <- 1
          }
          propState <- current + multiple  * move
        }
        if(nonUniform)
        {
          moveProb <- runif(1)
          for(l in 1:nMoves){
            if(moveProb <= sum(moveDist[1:l])/counter){
              move <- moves[,l]
              whichMove <- l
              break()
            }
          }
          move <- sample(c(-1,1), 1) * move
          propState <- current + move
        }
        else{
          move      <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
          propState <- current + move
        }
        
        if(SIS){
          if(runif(1) <= .05){
            propState <- sis_table(config, suffStats)
          }
        }
        
        if(any(propState < 0)){
          prob <- 0
        } else {
          if(dist == "hypergeometric"){
            prob <- exp( sum(lfactorial(current)) - sum(lfactorial(propState)) )
          } else { # dist == "uniform"
            prob <- 1
          }
        }
        probTotal <- probTotal + min(1, prob)
        
        if(nonUniform){
          if(unifs[k*(thin-1)+j] < prob){
            current  <- propState
            moveDist[whichMove] <- moveDist[whichMove] + 1
            counter <- counter + 1
          }
        }else{
          if(unifs[k*(thin-1)+j] < prob) current <- propState # else current
        }
        
        totalRuns <- totalRuns + 1        
      }
      
      state[,k] <- current    
    }
    message("done.")  
    ## format output
    out <- list(
      steps = state, 
      moves = moves, 
      accept_prob = probTotal / totalRuns
    )
    
    
    
    
  }
  
  ## in Cpp
  ##################################################
  if(engine == "Cpp"){
    
    current   <- unname(init)  
    allMoves  <- cbind(moves, -moves)  
    sampler   <- if(dist == "hypergeometric") {
      metropolis_hypergeometric_cpp
    } else {
      metropolis_uniform_cpp
    }
    message("Running chain (C++)... ", appendLF = FALSE)  
    if (burn > 0) current <- sampler(current, allMoves, suffStats, config, burn, 1, hitAndRun, SIS, nonUniform, adaptive)$steps[,burn]
    out       <- sampler(current, allMoves, suffStats, config, iter, thin, hitAndRun, SIS, nonUniform, adaptive)
    out$moves <- moves
    message("done.")
    
  }
  
  ## return output
  ##################################################  

  out[c("steps", "moves", "accept_prob")]
  
  # This is a change to the function.
}








#' @rdname metropolis
#' @export
rawMetropolis <- function(init, moves, iter = 1E3, dist = "hypergeometric", hitAndRun = FALSE, SIS = FALSE, nonUniform = FALSE, adaptive = FALSE){
  metropolis(init, moves, iter, burn = 0, thin = 1, dist = dist, hitAndRun, SIS, nonUniform, adaptive) 
}

