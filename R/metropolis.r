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
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param dist steady-state distribution; "hypergeometric" (default)
#'   or "uniform"
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @param hit_and_run Whether or not to use the discrete hit and run algorithm in
#'   the metropolis algorithm
#' @param adaptive Option inside hit_and_run option. If TRUE, hit and run will choose a proposal distribution adaptively. Defaulted to FALSE.
#' @param SIS If TRUE, with a small probability the move will be chosen randomly from the uniform distribution 
#'  on the fiber using Sequential Importance "Like" Sampling methods. Defaulted to FALSE
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
#' # sampling from the hypergeometric distribution
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
#' outC <- metropolis(tab2vec(handy), moves, 1e4, engine = "Cpp")
#' str(outC)
#' outR <- metropolis(tab2vec(handy), moves, 1e4, engine = "R", thin = 20)
#' str(outR)
#' 
#' # showSteps(out$steps)
#' 
#' 
#' library(microbenchmark)
#' microbenchmark(
#'   metropolis(tab2vec(handy), moves, engine = "Cpp"),
#'   metropolis(tab2vec(handy), moves, engine = "R")
#' )
#' 
#' # cpp ~ 20-25x faster
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
metropolis <- function(init, moves, suff_stats, config, iter = 1E3, burn = 0, thin = 1,
                       dist = c("hypergeometric","uniform"), engine = c("Cpp","R"), 
                       hit_and_run = FALSE, SIS = FALSE, non_uniform = FALSE, adaptive = FALSE
){
  
  ## preliminary checking
  ##################################################
  dist <- match.arg(dist)
  engine <- match.arg(engine)
  if(thin == 0){
    message("thin = 1 corresponds to no thinning, resetting thin = 0.")
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
    if(non_uniform == TRUE){
      move_dist <- rep(1,nMoves)
      counter <- nMoves
    }
    
    unifs <- runif(burn)
    if(non_uniform == TRUE){
       move_prob <- runif(burn)
    }
    if(burn > 0) {
      for(k in 1:burn){
        #Hit and Run option
        if(hit_and_run)
        {
          move <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
          w_move <- move[move != 0]
          w_current <- current[move != 0]
          w_moves <- -1 * w_current / w_move
          lower_bound <- if(any(w_moves < 0)){max(subset(w_moves,subset = w_moves < 0))}else{1}
          upper_bound <- if(any(w_moves > 0)){min(subset(w_moves,subset = w_moves > 0))}else{-1}
          
          if(any(w_moves == 0)){
          w_propStatelow <- current + lower_bound * move
          w_propStateup <-  current + upper_bound * move
          if(any(w_propStatelow < 0)){
            lower_bound <- 1
          }
          if(any(w_propStateup < 0)){
            upper_bound <- -1
          }
          }
          
          c_s <- sample(lower_bound:upper_bound,1)
          if(c_s == 0){
            c_s <- 1
          }
          propState <- current + c_s * move
        }
        
        #Non-uniform move sampling option
        if(non_uniform == TRUE)
        {
          
          for(l in 1:nMoves){
            if(move_prob[k] <= sum(move_dist[1:l])/counter){
              move <- moves[,l]
              which_move <- l
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
        
        if(non_uniform == TRUE){
          if(unifs[k] < prob){
            current  <- propState
            move_dist[which_move] <- move_dist[which_move] + 1
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
        
        if(hit_and_run)
        {
          move <- moves[,sample(nMoves,1)]
          w_move <- move[move != 0]
          w_current <- current[move != 0]
          w_moves <- (-1 * w_current) / w_move
          lower_bound <- if(any(w_moves < 0)){max(subset(w_moves,subset = w_moves < 0))}else{1}
          upper_bound <- if(any(w_moves > 0)){min(subset(w_moves,subset = w_moves > 0))}else{-1} 
        
          #New part 
         # #Option 1 Enumerate tables and 
         # line <- lower_bound:upper_bound
         # #Enumerate tables on the line
         # tables <- matrix(0L, nrow =length(init) , ncol = length(line))
         # for(i in 1:length(line)){
         #   tables[,i] <- current + line[i]*move
         # }
         # probs <- apply(tables, 2, function(x) 1/(sum(lfactorial(x))))
         # prob_dist <- probs / sum(probs)
         # unif <- runif(1)
         # dummy <- prob_dist[1]
         # for(i in 1:(length(prob_dist) -1)){
         #   if(unif < dummy){
         #     propState <- tables[,i]
         #     break()
         #   }
         #   dummy <- dummy + prob_dist[i+1]
         # }
          #Option 2
         # line <- lower_bound:upper_bound
         # w_current <- current
         # unifs2 <- runif(2*length(line))
         # for(i in 1:(2*length(line))){
         #   w_propState <- w_current + sample(c(-1,1), 1)*move
         #   if(any(w_propState < 0)){
         #     prob <- 0
         #   } else {
         #     if(dist == "hypergeometric"){
         #       prob <- exp( sum(lfactorial(w_current)) - sum(lfactorial(w_propState)) )
         #     } else { # dist == "uniform"
         #       prob <- 1
         #     }
         #   }
         #   if(unifs2[i] < prob) w_current <- w_propState # else w_current
         # }
         # propState <- w_current
          
          
          if(any(w_moves == 0)){
            w_propStatelow <- current + lower_bound * move
            w_propStateup <-  current + upper_bound * move
            if(any(w_propStatelow < 0)){
              lower_bound <- 1
            }
            if(any(w_propStateup < 0)){
              upper_bound <- -1
            }
          }
         
          c_s <- sample(lower_bound:upper_bound,1)
          if(c_s == 0){
            c_s <- 1
          }
          propState <- current + c_s * move
        }
        if(non_uniform == TRUE)
        {
          move_prob <- runif(1)
          for(l in 1:nMoves){
            if(move_prob <= sum(move_dist[1:l])/counter){
              move <- moves[,l]
              which_move <- l
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
        
        if(SIS == TRUE){
          if(runif(1) <= .05){
            propState <- sis_table(config, suff_stats)
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
        
        if(non_uniform == TRUE){
          if(unifs[k*(thin-1)+j] < prob){
            current  <- propState
            move_dist[which_move] <- move_dist[which_move] + 1
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
      acceptProb = probTotal / totalRuns
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
    if (burn > 0) current <- sampler(current, allMoves, suff_stats, config, burn, 1, hit_and_run, SIS, non_uniform, adaptive)$steps[,burn]
    out       <- sampler(current, allMoves, suff_stats, config, iter, thin, hit_and_run, SIS, non_uniform, adaptive)
    out$moves <- moves
    message("done.")
    
  }
  
  
  ## return output
  ##################################################  
  
  out[c("steps", "moves", "acceptProb")]
}








#' @rdname metropolis
#' @export
rawMetropolis <- function(init, moves, iter = 1E3, dist = "hypergeometric", hit_and_run = FALSE, SIS = FALSE, non_uniform = FALSE, adaptive = FALSE){
  metropolis(init, moves, iter, burn = 0, thin = 1, dist = dist, hit_and_run) 
}

