#' Fitting generalized linear models with algebraic methods
#' 
#'
#' @param model model specification, either in terms of a configuration matrix or a symbolic 
#'   description of the model to be fitted
#' @param data data, as a data frame of raw data with ordinal discrete covariates
#' @param family a description of the error distirbution and link function used in the model
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @param moves the markov moves for the mcmc (as columns of a
#'   matrix).
#' @param ... ...
#' @return a list containing named elements \itemize{ \item
#'   \code{steps}: an integer matrix whose columns represent
#'   individual samples from the mcmc. \item \code{moves}: the moves
#'   used for the proposal distribution in the mcmc, computed with
#'   4ti2 (note that only the positive moves are given). \item
#'   \code{acceptProb}: the average acceptance probability of the
#'   moves, including the thinned moves. \item \code{call}: the call. 
#'   \item \code{obs}: the summarized data. 
#'    \item \code{A}: the sufficient
#'   statistics computing matrix.
#'   \item \code{sufficientStatistics}: The sufficient statistics of the model.
#'   \item \code{p.value}: the exact p-values of individual tests,
#'   accurate to Monte-Carlo error.  these are computed as the
#'   proportion of samples with statistics equal to or larger than
#'   the oberved statistic. \item \code{mid.p.value}: the mid
#'   p.values, see Agresti pp.20--21. \item \code{sampsStats}: 
#'   the statistics computed for each mcmc
#'   sample. \item \code{cells}: the number of cells in the table. }
#' @examples 
#' 
#'  library(ggplot2);theme_set(theme_bw())
#'  
#'  # generating data and running a poisson regression model
#'    # pick beta 0 and beta 1
#'     b0 <- 1; b1 <- 0.3
#'
#'    # generate data
#'     n <- 5000
#'     x <- sample(1:5, n, replace = T)
#'     y <- rpois(n, lambda = exp(b0 + b1*x))
#'     df <- data.frame(
#'      x = x, 
#'      y = y
#'     )
#'     
#'     # function output
#'     out <- aglm(y ~ x, data = df, family = poisson(), thin = 2000)
#'     
#'     # check convergence through trace plot
#'     qplot(1:10000, out$sampsStats$PRs, geom = "line")
#'   
#'    # compare aglm and glm predictions with the truth
#'     
#'     # model fitting with glm
#'     mod <- glm(y ~ x, data = df, family = poisson())
#'     
#'     # truth
#'     exp(b0 + b1*(1:5))
#'      
#'     # glm predictions
#'     predict(mod, data.frame(x = 1:5), type = "response")
#'     
#'     # aglm predictions
#'     rowMeans(out$steps) / ddply(df, "x", nrow)$V1
#'     
#'     
#'     
#'  # generating data and running a logistic regression model
#'  
#'   # helper functions
#'   link    <- function(p) log(p/(1-p))
#'   invlink <- function(x) 1 / (1+exp(-x))
#'   
#'   
#'   # create a fake data set
#'   
#'   # one covariate
#'   b0 <- 0.5; b1 <- 0.2
#'   
#'   n <- 100
#'   x <- sample(1:5, n, replace = T)
#'   y <- rbinom(n = n, size = 1, prob = invlink(b0 + b1*x))
#'   df <- data.frame(
#'     x = x, 
#'     y = y
#'   )
#'   
#'   # aglm 
#'   out <- aglm(y ~ x, data = df, family = binomial())
#'  
#'   # check convergence through trace plot
#'   qplot(1:10000, out$sampsStats$PRs, geom = "line")
#'   
#'   # using glm
#'   mod <- glm(y ~ x, data = df, family = binomial())
#'   
#'   # truth 
#'   invlink(b0 + b1*out$obs[,1])
#'   
#'   # glm predictions
#'   predict(mod, data.frame(x = c(1:5)), type = "response")
#'   
#'   # aglm predictions 
#'   rowMeans(out$steps) / ddply(df, "x", nrow)$V1
#'   
#' @export 


aglm <- function(model, data, family = poisson(),
                     iter = 1E4, burn = 10000, 
                     thin = 100, engine = c("Cpp","R"), 
                      moves, 
                     ...)
{
  ## set/check args
  ##################################################
  
  engine  <- match.arg(engine)
  argList <- as.list(match.call(expand.dots = TRUE))[-1]
  
  if("formula" %in% names(argList)){
    .Deprecated(msg = 
                  'the formula argument is deprecated, please use "model" instead.'
    )
  }
  
  ## check family and link argument
  ############################################################

  method <- family$family
  link <- family$link
  if (!(method == "poisson" | method == "binomial")) {
    stop("only poisson and logistic regression are implemented")
  }
  
  if (method == "poisson" & !link == "log") {
    stop("for poisson regression, only the log link function is implemented")
  }
      
  if (method == "binomial" & !link == "logit") {
    stop("for logistic regression, only the logit link function is implemented")
  }
  
  ## parse model specification (formula for vector of r_k's)
  ##################################################
  
  modelGivenByMatrix <- ifelse(is.matrix(model), TRUE, FALSE)
  
  if(modelGivenByMatrix){
    A <- model
    data   <- suppressMessages(teshape(data, "tab"))
    init <- tab2vec(data)
    nCells <- length(init)
  } else {
    # if it's a formula, convert to list
    if(is.formula(model)){ 
      ## reshape data
      data <- model.frame(model, data)
      
      # name data
      vars  <- names(data)
      
      ## parse formula
      fString    <- as.character(model)
      response   <- fString[2]
      predString <- fString[3]
      
      
      ## make list of facets
      model <- strsplit(predString, " \\+ ")[[1]]
      model <- strsplit(model, " \\* ")
      
      ## format the data 
      names(data)[names(data) == response] <- "response"
      
      if (method == "binomial") {
        data <- rbind(ddply(data, unique(unlist(model)), summarise, sum = sum(response)), 
                      ddply(data, unique(unlist(model)), summarise, sum = length(response) - sum(response))
        )
      } else {
        data <- ddply(data, unique(unlist(model)), summarise, sum = sum(response))
      }
      
      ## any 0 levels
      if (any(data[,-ncol(data)] == 0)) {
        stop("Cannot have a covariate with a level 0")
      }
      
      if(length(model) == 1){
        
        init  <- data$sum
        nCells <- length(init)
        p     <- 1
        
      } else {
        # if model specifiaction, then make table
        p <- ncol(data) - 1
        init <- data$sum
        nCells <- length(init)
      }
    } 
    
    
    
    # make facets (list of index vecs); if model specified with variable
    # names, convert them to indices
    if(all(unlist(model) %in% vars)){ # variable names      
      varname2index <- 1:p
      names(varname2index) <- vars[vars != response]      
      facets <- lapply(model, function(varsInFacet) varname2index[varsInFacet])
    } else if(all(unlist(model) %in% 1:length(vars))){ # by indices
      facets <- lapply(model, as.integer) # to fix the ~ 1 + 2 case, parsed as chars
    } else {
      stop("Invalid model specification, see ?aglm")
    }
    ## levels (assuming all levels are numeric i.e. (1,2,3,...  not Green, Blue, Red, etc.)
    if(ncol(data) <= 2){ 
      levls <- unique(data[,-ncol(data)])
    } else {
      levls <- lapply(data[,-ncol(data)], unique)
    }
    # make configuration (model) matrix
    A <- pmat(levls, facets)
  }
  
  # check to see if all level configurations are there (need work here)
  lvlsInData <- as.list(as.data.frame(t(expand.grid(levls)))) %in% as.list(as.data.frame(t(data[,-ncol(data)])))
  
  # subset A by levels that are present
  A <- A[,lvlsInData]
  
  # if family = "binomial" compute the lawernce lifting of A
  if (method == "binomial") {
    A <- lawrence(A)
  }
  # find the sufficient statistics
  suffStats <- unname(A %*% init)
  
  ## construct A matrix and compute moves
  ##################################################  
  
  if(missing(moves) && !is.null(getOption("4ti2_path"))){
    
    message("Computing Markov moves (4ti2)... ", appendLF = FALSE)  	
    moves <- markov(A)
    message("done.", appendLF = TRUE)      
    
  } else if(missing(moves) && is.null(getOption("4ti2_path"))){
    
    warning(
      "No moves were provided and 4ti2 is not found.\n",
      "  The resulting chain is likely not connected and strongly autocorrelated.\n",
      "  See ?pois_reg.  Consider using rmove to generate SIS moves in advance.",
      immediate. = TRUE
    )
    message("Computing 1000 SIS moves... ", appendLF = FALSE)    
    moves <- rmove(n = 1000, A = A, b = A %*% tab2vec(data), ...)
    message("done.", appendLF = TRUE)      
    
  } else if(is.character(moves)){
    
    movesMat <- NULL
    stopifnot(all(moves %in% c("lattice", "markov", "groebner", "grobner", "graver", "sis")))
    if("lattice"  %in% moves)  movesMat <- cbind(movesMat,   zbasis(A))
    if("markov"   %in% moves)  movesMat <- cbind(movesMat,   markov(A))
    if("groebner" %in% moves)  movesMat <- cbind(movesMat, groebner(A))
    if("grobner"  %in% moves)  movesMat <- cbind(movesMat, groebner(A))
    if("graver"   %in% moves)  stop("graver not yet implemented.")
    moves <- movesMat
    
  }
  
  stopifnot(is.array(moves))
  
  
  
  ## run metropolis-hastings
  ##################################################  
  init <- unname(init) # init
  out <- metropolis(init, moves, suffStats = suffStats, config = unname(A), iter = iter, burn = burn, thin = thin, 
                    engine = engine, ...)  

  
  u <- t(t(data$sum))
  PR <- computeUProbsCpp(matrix(u))  # unnormd prob; numers LAS 1.1.10
  
  
  ## compute MCMC chi squares
  ##################################################  
  PRs <- computeUProbsCpp(out$steps) # unnormd probs; numers LAS 1.1.10
  
  ## add A matrix, p.value and return
  ##################################################  
  out$call <- match.call()   
  out$obs <- data  
  out$A <- A
  out$sufficientStatistics <- suffStats
  
  out$p.value <- c(
    PR = mean(PRs <= PR)
  )
  
  out$p.value.std.err <- c(
    PR = sqrt(mean(PRs <= PR)*(1-mean(PRs <= PR))/iter)
  )  
  
  out$mid.p.value <- c(
    PR = mean(PRs < PR) + mean(PRs == PR)/2
  )  
  
  
  out$iter       <- iter
  out$burn       <- burn
  out$thin       <- thin
  out$statistic  <- c(PR = PR)
  out$sampsStats <- list(PRs = PRs)
  out$cells      <- nCells
  out$method     <- method
  
  class(out) <- "aglm"
  out
  
  
  
  
  
  
  
}