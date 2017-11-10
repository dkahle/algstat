#' Fit a Logistic Regression model with algebraic methods
#' 
#' 
#' 
#' @param model hierarchical poisson model specification
#' @param data data, as a data frame with raw data with discrete covariates
#' @param init the initialization of the chain. by default, this is
#'   the observed table
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @param method should the expected value (exp) be fit using
#'   iterative proportional fitting (via loglin) or the MCMC as the
#'   average of the steps?
#' @param moves the markov moves for the mcmc (as columns of a
#'   matrix).
#' @param ... ...
#' @return a list containing named elements \itemize{ \item
#'   \code{steps}: an integer matrix whose columns represent
#'   individual samples from the mcmc. \item \code{moves}: the moves
#'   used for the proposal distribution in the mcmc, computed with
#'   4ti2 (note that only the positive moves are given). \item
#'   \code{acceptProb}: the average acceptance probability of the
#'   moves, including the thinned moves. \item \code{param}: the
#'   fitted parameters of the log linear model. \item \code{df}:
#'   parameters per term in the model \item \code{quality}: model
#'   selection statistics AIC, AICc, and BIC. \item 
#'   \code{residuals}: the (unstandardized) pearson residuals (O -
#'   E) / sqrt(E) \item \code{call}: the call. \item \code{obs}: the
#'   contingency table given. \item \code{exp}: the fit contingency
#'   table as an integer array. \item \code{A}: the sufficient
#'   statistics computing matrix (from Tmaker). \item 
#'   \code{p.value}: the exact p-values of individual tests,
#'   accurate to Monte-Carlo error.  these are computed as the
#'   proportion of samples with statistics equal to or larger than
#'   the oberved statistic. \item \code{mid.p.value}: the mid
#'   p.values, see Agresti pp.20--21. \item \code{statistic}: the
#'   pearson's chi-squared (X2), likelihood ratio (G2), 
#'   Freeman-Tukey (FT), Cressie-Read (CR), and Neyman modified
#'   chi-squared (NM) statistics computed for the table given. \item
#'   \code{sampsStats}: the statistics computed for each mcmc
#'   sample. \item \code{cells}: the number of cells in the table.
#'   \item \code{method}: the method used to estimate the table. }
#'   @importFrom stats model.frame
#'  @export log_reg

log_reg <- function(model, data,
                    iter = 1E4, burn = 1000, 
                    thin = 10, engine = c("Cpp","R"), 
                    method = c("ipf", "mcmc"), moves, 
                    hit_and_run = FALSE,
                    SIS = FALSE,
                    non_uniform = FALSE,
                    adaptive = FALSE,
                    ...)
{
  
  ## set/check args
  ##################################################
  
  engine  <- match.arg(engine)
  method  <- match.arg(method)  
  argList <- as.list(match.call(expand.dots = TRUE))[-1]
  
  if("formula" %in% names(argList)){
    .Deprecated(msg = 
                  'the formula argument is deprecated, please use "model" instead.'
    )
  }
  
  
  ## reshape data
  ##################################################
  
  # data   <- suppressMessages(teshape(data, "freq"))
  # p      <- length(dim(data))
  # nCells <- length(data)
  
  ## if a pure array is given, give names for later
  # if(is.array(data) && is.null(dimnames(data))) data <- array2tab(data)
  
  ## other basic objects
  #varsNlevels <- dimnames(data)  
  
  
  
  
  
  ## check for sampling zeros
  ##################################################
  #if(any(data == 0L)) message(
  #  "Care ought be taken with tables with sampling zeros to ensure the MLE exists."
  #)
  
  
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
      data <- rbind(ddply(data, unique(unlist(model)), summarise, sum = sum(response)), 
                    ddply(data, unique(unlist(model)), summarise, sum = length(response) - sum(response))
                    )
      
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
      stop("Invalid model specification, see ?log_reg")
    }
    ## levels (assuming all levels are numeric i.e. (1,2,3,...  not Green, Blue, Red, etc.)
    
    if(ncol(data) <= 2){ 
      lvls <- unique(data[,-ncol(data)])
    } else {
      lvls <- lapply(data[,-ncol(data)], unique)
    }
    # make configuration (model) matrix
    A <- pmat(lvls, facets)
  }
  
  # check to see if all level configurations are there (need work here)
  lvlsInData <- as.list(as.data.frame(t(expand.grid(lvls)))) %in% as.list(as.data.frame(t(data[,-ncol(data)])))
  
  # subset A by levels that are present
  A <- A[,lvlsInData]
  
  # compute the Lawrence lifting of A
  A <- lawrence(A)
  
  # find the sufficient statistics
  suff_stats <- unname(A %*% init)
  
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
  out <- metropolis(init, moves, suff_stats = suff_stats, config = unname(A), iter = iter, burn = burn, thin = thin, 
                    engine = engine, hit_and_run = hit_and_run, SIS = SIS, non_uniform = non_uniform, adaptive = adaptive)  
  
  
  
  ## compute data chi square
  ##################################################  
  # if(modelGivenByMatrix && method == "ipf"){
  #   message(
  #     "Iterative proportional fitting is not yet implemented\n", 
  #     "  for models specified by configuration matrices.\n",
  #     "  Changing to method = \"mcmc\"..."
  #   )
  #   method <- "mcmc"
  # }
  # if(method == "ipf"){
  #   exp <- loglin(data, facets, fit = TRUE, print = FALSE)$fit
  # } else if(method == "mcmc"){
  #   exp <- vec2tab(rowMeans(out$steps), dim(data))
  #   dimnames(exp) <- dimnames(data)
  # }
  # e <- unname(tab2vec(exp))
  
  u <- t(t(data$sum))
  PR <- computeUProbsCpp(matrix(u))  # unnormd prob; numers LAS 1.1.10
  # X2 <- computeX2sCpp(u, e)  
  # G2 <- computeG2sCpp(u, e)    
  # FT <- computeCRsCpp(u, e, -.5)      
  # CR <- computeCRsCpp(u, e, 2/3)
  # NM <- computeNMsCpp(u, e)  
  
  
  ## compute MCMC chi squares
  ##################################################  
  PRs <- computeUProbsCpp(out$steps) # unnormd probs; numers LAS 1.1.10
  # X2s <- computeX2sCpp(out$steps, e)  
  # G2s <- computeG2sCpp(out$steps, e) 
  # FTs <- computeCRsCpp(out$steps, e, -.5)
  # CRs <- computeCRsCpp(out$steps, e, 2/3)  
  # NMs <- computeNMsCpp(out$steps, e)   
  
  # 
  #   ## compute parameters
  #   ##################################################      
  #   if(!modelGivenByMatrix){
  #   # in principle, there should be one parameter for every cell.
  #   # there are prod(dim(data)) cells.
  #   # a good reference is BFH, p. 35 (and to a lesser extent 43)
  #   # the prod(dim(data)[terms[[j]]] - 1) line below is like
  #   # (I - 1) (J - 1) (K - 1)
  #   # CDA p.79 also helpful
  #   dimSatModel <- nCells - 1
  #   degFreedom <- rep.int(0, 2^p) # there are 2^p possible subsets of vars, and
  #                                 # therefore there are 2^p possible terms
  #                                 
  #   # possibleTerms are more "types of terms" as opposed to individual terms
  #   # for example, an entry c(1,3) would refer to all combinations of levels
  #   # of variables 1 and 3; ie (# var 1 levels - 1) * (# var 3 levels - 1)
  #   # individual terms (parameters)
  #   possibleTerms <- subsets(p, include_null = TRUE)
  #   names(possibleTerms) <- sapply(possibleTerms, paste, collapse = " ")
  #   names(possibleTerms)[which(names(possibleTerms) == "")] <- "(Intercept)"    
  #   nVarLvls <- dim(data)
  #   # paramsPerTerm <- lapply(possibleTerms, function(x){
  #   #   if(length(x) == 0) return(1L)
  #   #   prod(nVarLvls[x] - 1)
  #   # })
  #   
  #   
  #   # similarly, there are the terms in the model
  #   termsInModel <- unique(unlist(lapply(
  #     lapply(facets, as.character), # to avoid subsets(2)
  #     subsets, include_null = TRUE), 
  #     recursive = FALSE
  #   ))
  #   termsInModel <- lapply(termsInModel, as.integer)
  #   names(termsInModel) <- sapply(termsInModel, paste, collapse = " ")  
  #   names(termsInModel)[which(names(termsInModel) == "")] <- "(Intercept)"
  #   paramsPerTermInModel <- lapply(termsInModel, function(x){
  #     if(length(x) == 0) return(1L) 
  #     prod(nVarLvls[x] - 1)
  #   })
  #   names(paramsPerTermInModel) <- unname(sapply(termsInModel, function(x){
  #     if(length(x) == 0) return("(Intercept)")
  #     paste(names(dimnames(data))[x], collapse = ".")
  #   }))
  #   nParamsInModel <- sum(unlist(paramsPerTermInModel))
  #   dimModel <- nParamsInModel - 1 # the - 1 accounts for the overall mean
  #   overallAsymptoticDegFreedom <- (dimSatModel - dimModel)
  #   
  # 
  #   # compute the parameters  
  #   log_fit <- exp
  #   log_fit[exp > 0] <- log(exp[exp > 0])  
  #   param <- as.list(rep(NA, length(termsInModel)))
  #   names(param) <- names(paramsPerTermInModel) 
  #   for(k in seq_along(param)){
  #     if(length(termsInModel[[k]]) == 0){
  #       param[[k]] <- mean(log_fit)
  #       log_fit <- log_fit - param[[k]]
  #     } else {
  #       param[[k]] <- apply(log_fit, termsInModel[[k]], mean)
  #       log_fit <- sweep(log_fit, termsInModel[[k]], param[[k]])
  #     }
  #   }
  #   # for every step, fit mle
  #   # then decompose mle
  #   # problem : they all have the same marginals, so the same
  #   # mles!
  #   # idea 1 : sample from the multinomial with the same sample
  #   # size (so different marginals), estimate, then decompose
  #   # idea 2 : bootstrap sample from the table, estimate, decompose
  #   # i think i like idea 2 better.
  #   
  # 
  #   # reorder the param estimates in the order of subsets
  #   # so you have the intercept, then all first order terms, and so on
  #   goodOrder <- sapply(
  #     c("(Intercept)", subsets(names(dimnames(data)))),
  #     paste, collapse = "."
  #   )
  #   param <- param[goodOrder[goodOrder %in% names(param)]]
  #   out$param <- param
  #   
  #   }
  #   
  ## compute residuals and model selection, agresti p.81, 216, 324
  ##################################################  
  # out$residuals <- exp
  # out$residuals[exp > 0] <- 
  #   (data[exp > 0] - exp[exp > 0]) / sqrt(exp[exp > 0])
  # 
  # if(!modelGivenByMatrix){
  #   k <- nParamsInModel  # = number of params 
  #   n <- sum(data)       # = sample size
  #   L <- dmultinom(u, sum(u), e, TRUE) # maximized log-likelihood
  #   BIC  <- log(n)*k - 2*L
  #   AIC  <-      2*k - 2*L
  #   AICc <- AIC + 2*k*(k+1)/(n-k-1)
  #   out$df <- paramsPerTermInModel
  #   out$quality <- c(AIC = AIC, AICc = AICc, BIC = BIC)
  # }
  
  ## add A matrix, p.value and return
  ##################################################  
  out$call <- match.call()   
  out$obs <- data  
  # out$exp <- exp
  # out$A <- A
  
  out$p.value <- c(
    PR = mean(PRs <= PR)
  )
  
  # out$p.value <- c(
  #   PR = mean(PRs <= PR),   
  #   X2 = mean(X2s >= X2), 
  #   G2 = mean(G2s >= G2),   
  #   FT = mean(FTs >= FT),
  #   CR = mean(CRs >= CR),
  #   NM = mean(NMs >= NM)
  # )
  
  out$p.value.std.err <- c(
    PR = sqrt(mean(PRs <= PR)*(1-mean(PRs <= PR))/iter)
  )  
  
  # out$p.value.std.err <- c(
  #   PR = sqrt(mean(PRs <= PR)*(1-mean(PRs <= PR))/iter), 
  #   X2 = sqrt(mean(X2s >= X2)*(1-mean(X2s >= X2))/iter), 
  #   G2 = sqrt(mean(G2s >= G2)*(1-mean(G2s >= G2))/iter),   
  #   FT = sqrt(mean(FTs >= FT)*(1-mean(FTs >= FT))/iter),
  #   CR = sqrt(mean(CRs >= CR)*(1-mean(CRs >= CR))/iter), 
  #   NM = sqrt(mean(NMs >= NM)*(1-mean(NMs >= NM))/iter)     
  # )  
  
  out$mid.p.value <- c(
    PR = mean(PRs < PR) + mean(PRs == PR)/2
  )  
  
  # out$mid.p.value <- c(
  #   PR = mean(PRs < PR) + mean(PRs == PR)/2,
  #   X2 = mean(X2s > X2) + mean(X2s == X2)/2, 
  #   G2 = mean(G2s > G2) + mean(G2s == G2)/2,
  #   FT = mean(FTs > FT) + mean(FTs == FT)/2,
  #   CR = mean(CRs > CR) + mean(CRs == CR)/2,
  #   NM = mean(NMs > NM) + mean(NMs == NM)/2    
  # )  
  
  out$iter       <- iter
  out$burn       <- burn
  out$thin       <- thin
  out$statistic  <- c(PR = PR)
  # out$statistic  <- c(PR = PR, X2 = X2, G2 = G2, FT = FT, CR = CR, NM = NM)
  out$sampsStats <- list(PRs = PRs)
  # out$sampsStats <- list(PRs = PRs, X2s = X2s, G2s = G2s, FTs = FTs, CRs = CRs, NMs = NMs)
  out$cells      <- nCells
  out$method     <- method
  
  class(out) <- "logistic"
  out
}