#' Fitting generalized linear models with algebraic methods
#' 
#'[aglm()] fits a generalized linear model to a dataset (typically a
#' contingency table) and performs an exact conditional test on the fitted model.  
#' The exact test, which is a goodness-of-fit test, is performed via 
#' Monte-Carlo sampling from the conditional distribution of the table 
#' given the sufficient statistics of the model.  In short, inference 
#' is drawn by comparing the statistic of the observed table (the unnormalized log-likelihood)
#' to those of the samples.  The proportion of sampled tables with equal to or
#' more extreme values than the observed table is the resulting p-value.
#' 
#' @param model model specification, either in terms of a configuration matrix or a symbolic 
#'   description of the model to be fitted
#' @param data data, as a data frame of raw data with ordinal discrete covariates
#' @param family a description of the error distirbution and link function used in the model
#' @param control a list of arguments that control the MCMC algorithm
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
#'   the oberved statistic. 
#'   \item \code{mid.p.value}: the mid
#'   p.values, see Agresti pp.20--21. 
#'   \item \code{sampsStats}: 
#'   the statistics computed for each mcmc
#'   sample. 
#'   \item \code{cells}: the number of cells in the table. }
#' @examples 
#' 
#'  # simple poisson regression example from Diaconis et. al. 1998
#'  
#'  data <- data.frame(x = c(0,1,2,3,4), y = c(44,25,21,19,11))
#' 
#'  # perform the fitting
#'  fit_aglm <- aglm(y ~ x, data = data, family = poisson())
#' 
#'  # aglm gof p-value
#'  fit_aglm$p.value
#'  
#'  # compare with glm   
#'  fit_glm  <-  glm(y ~ x, data = data, family = poisson())
#'  
#'  # glm gof p-value
#'  1 - pchisq(fit_glm$deviance, fit_glm$df.residual)
#'  
#'  
#'  # logistic regression example with data from Haberman 1974
#'  
#'  # data in elongated format
#'  x <- c(rep(0,6),rep(1,2),rep(2,4),rep(3,9),rep(4,10),rep(5,20),
#'  rep(6,34),rep(7,42),rep(8,124),rep(9,58),rep(10,77),rep(11,95),
#'  rep(12,360))
#'  
#'  y <- c(rep(0,2),rep(1,10),rep(0,3),rep(1,6),rep(0,5),rep(1,5),
#'  rep(0,7),rep(1,13),rep(0,9),rep(1,25),rep(0,15),rep(1,27),rep(0,49),
#'  rep(1,75),rep(0,29),rep(1,29),rep(0,45),rep(1,32),rep(0,59),
#'  rep(1,36),rep(0,245),rep(1,115))
#'  
#'  data_aglm <- data.frame(x=x,y=y)
#'  data_glm <- data.frame(x = 0:12)
#'  data_glm$y <- cbind(c(4, 2, 4, 6, 5, 13, 25, 27, 75, 29, 32, 36, 115),
#'                      c(2, 0, 0, 3, 5, 7, 9, 15, 49, 29, 45, 59, 245))
#'
#'  fit_aglm <- aglm(y ~ x, data = data_aglm, family = binomial())  
#'  fit_aglm$p.value  
#'  
#'  # compare to glm
#'  
#'  # special formatting of data needed for correct analysis
#'  data_glm <- data.frame(x = 0:12)
#'  data_glm$y <- cbind(c(4, 2, 4, 6, 5, 13, 25, 27, 75, 29, 32, 36, 115),
#'                      c(2, 0, 0, 3, 5, 7, 9, 15, 49, 29, 45, 59, 245))
#'                      
#'  fit_glm <-   glm(y ~ x, data = data_glm, family = binomial())
#'  1 - pchisq(fit_glm$deviance, fit_glm$df.residual)
#'  
#'  
#'  # multiple poisson regression 
#'  
#'  # create a fake dataset 
#'  b0 <- 1; b1 <- 0.5; b2 <- 0.5
#'  
#'  covariates <- expand.grid(1:3,1:4)
#'  
#'  data <- data.frame(
#'    x1 = covariates[,1],
#'    x2 = covariates[,2],
#'    y = rpois(12,exp(b0 + b1*covariates[,1] + b2*covariates[,2]))
#'  )
#'  
#'  # compare to glm
#'  
#'  fit_aglm <- aglm(y ~ x1 + x2, data = data, family = poisson())
#'  fit_glm  <-  glm(y ~ x1 + x2, data = data, family = poisson())
#'  
#'  fit_aglm$p.value
#'  with(fit_glm, 1 - pchisq(deviance, df.residual))
#'  
#' @export 


aglm <- function(model, data, family = poisson(),
                     control = list(...), moves, 
                     ...)
{

  ## set/check args
  ##################################################
  control <- do.call("aglm.control", control)
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
      
      ## extract full model 
      pred_string <- attr(terms(data), "term.labels")
      
      ## name data
      vars  <- names(data)
      response   <- vars[1]
      
      ## make list of facets
      model <- strsplit(pred_string, "\\:")
      
      ## format the data 
      names(data)[names(data) == response] <- "response"

      if (method == "binomial") {
        data <- group_by(data, .dots = unique(unlist(model)))
        success <- summarise(data, sum = sum(response))
        failure <- summarise(data, sum = length(response) - sum(response))
        data <- bind_rows(success, failure)
        
      } else {
        data <- group_by(data, .dots = unique(unlist(model)))
        data <- summarise(data, sum = sum(response))
      }
      
      ## any 0 levels
    #  if (any(data[,-ncol(data)] == 0)) {
    #    stop("Cannot have a covariate with a level 0")
    #  }
      
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
    facets <- lapply(facets, unname)
    ## levels (assuming all levels are numeric i.e. (1,2,3,...  not Green, Blue, Red, etc.)
    if(ncol(data) <= 2){ 
      levels <- unique(data[,-ncol(data)])
    } else {
      levels <- lapply(data[,-ncol(data)], unique)
      levels <- lapply(levels, sort)
    }
    # make configuration (model) matrix
    A <- pmat(levels, facets)
  }
  
 # check to see if all level configurations are there (need work here)
 # lvlsInData <- as.list(as.data.frame(t(expand.grid(levels)))) %in% as.list(as.data.frame(t(data[, -ncol(data)])))
  
  # subset A by levels that are present
 # A <- A[,lvlsInData]
  
  # if family = "binomial" compute the lawernce lifting of A
  if (method == "binomial") {
    A <- lawrence(A)
  }
  # find the sufficient statistics
  suffStats <- unname(A %*% init)
  
  ## construct A matrix and compute moves
  ##################################################  
  
  if(missing(moves) && has_4ti2()){
    
    message("Computing Markov moves (4ti2)... ", appendLF = FALSE)  	
    moves <- markov(A, p = "arb")
    message("done.", appendLF = TRUE)      
    
  } else if(missing(moves) && has_4ti2()){
    
    warning(
      "No moves were provided and has_4ti2() = FALSE.\n",
      "  SIS moves will be used; estimates will likely be biased.\n",
      "  Consider using rmove() to generate SIS moves in advance.",
      immediate. = TRUE
    )
    message("Computing 1000 SIS moves... ", appendLF = FALSE)    
    moves <- rmove(n = 1000, A = A, b = A %*% tab2vec(data), ...)
    message("done.", appendLF = TRUE)      
    
  } else if(is.character(moves)){
    
    movesMat <- NULL
    stopifnot(all(moves %in% c("lattice", "markov", "groebner", "grobner", "graver", "sis")))
    if("lattice"  %in% moves)  movesMat <- cbind(movesMat,   zbasis(A, p = "arb"))
    if("markov"   %in% moves)  movesMat <- cbind(movesMat,   markov(A, p = "arb"))
    if("groebner" %in% moves)  movesMat <- cbind(movesMat, groebner(A, p = "arb"))
    if("grobner"  %in% moves)  movesMat <- cbind(movesMat, groebner(A, p = "arb"))
    if("graver"   %in% moves)  stop("graver not yet implemented.")
    moves <- movesMat
    
  }
  
  stopifnot(is.array(moves))
  
  
  
  ## run metropolis-hastings
  ##################################################  
  init <- unname(init) # init
  out <-
    metropolis(
      init,
      moves,
      iter = control$iter,
      burn = control$burn,
      thin = control$thin,
      engine = control$engine,
      hit_and_run = control$hit_and_run,
      adaptive = control$adaptive
    )
  
  
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
    PR = sqrt(mean(PRs <= PR)*(1-mean(PRs <= PR))/control$iter)
  )  
  
  out$mid.p.value <- c(
    PR = mean(PRs < PR) + mean(PRs == PR)/2
  )  
  
  
  out$iter       <- control$iter
  out$burn       <- control$burn
  out$thin       <- control$thin
  out$statistic  <- c(PR = PR)
  out$sampsStats <- list(PRs = PRs)
  out$cells      <- nCells
  out$method     <- method
  
  class(out) <- "aglm"
  out
}












aglm.control <- function(iter = 10000,
                         burn = 10000,
                         thin = 100, 
                         engine = c("C++", "R"), 
                         hit_and_run = FALSE,
                         adaptive = FALSE
) {
  engine <- match.arg(engine)
list(iter = iter, burn = burn, thin = thin, 
     engine = engine, hit_and_run = hit_and_run, 
     adaptive = adaptive)
}