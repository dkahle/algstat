#' Fit a hierarchical Poisson regression models with algebraic
#' methods
#' 
#' Fit a hierarchical Poisson regression models with algebraic
#' methods
#' 
#' @param model hierarchical log-linear model specification
#' @param data data, typically as a table but can be in different 
#'   formats.  see \code{\link{teshape}}
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
#' @export
#' @author David Kahle
#' @seealso \code{\link{loglin}}, \code{\link{loglm}}, 
#'   \code{\link{metropolis}}
#' @references Diaconis, P. and B. Sturmfels (1998). Algebraic 
#'   Algorithms for Sampling from Conditional Distributions. 
#'   \emph{The Annals of Statistics} 26(1), pp.363-397.
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). 
#'   \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser 
#'   Verlag AG.
#' @references Aoki, S., H. Hara, and A. Takemura (2012). 
#'   \emph{Markov Bases in Algebraic Statistics}, Springer.
#' @references Agresti, A. (2002). \emph{Categorical Data Analysis},
#'   Basel: John Wiley & Sons, 2ed.
#' @references Agresti, A. (1992). A Survey of Exact Inference for 
#'   Contingency Tables \emph{Statistical Science} 7(1), pp.131-153.
#' @references Read, T. and Cressie, N. (1998). 
#'   \emph{Goodness-of-Fit Statistics for Discrete Multivariate 
#'   Data}, Springer-Verlag.
#' @examples
#' 
#' \dontrun{
#' 
#' 
#' ## handedness introductory example
#' ############################################################
#' 
#' data(handy)
#' 
#' (out <- loglinear(~ Gender + Handedness, data = handy))
#' 
#' # you can also specify the same model using variable indices...   
#' (out <- loglinear(~ 1 + 2, data = handy))
#'   
#' # ... or as a list of facets given by indices
#' (out <- loglinear(list(1, 2), data = handy))
#' 
#' # ... or as a list of facets given by name
#' (out <- loglinear(list("Gender", "Handedness"), data = handy))
#' 
#' # ... and even via a pre-computed configuration matrix
#' # this method does come with somewhat reduced output
#' A <- hmat(c(2, 2), 1:2)
#' (out <- loglinear(A, data = handy))
#' 
#' 
#' 
#' # loglinear performs the same tasks as loglin and loglm,
#' # but loglinear gives the exact test p values and more goodness-of-fit statistics
#' stats::loglin(handy, list(1, 2))
#' MASS::loglm(~ Gender + Handedness, data = handy)
#' # loglm is just a wrapper of loglin  
#'   
#' # we can check loglinear's output with
#' fisher.test(handy)$p.value
#' out$p.value
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
#' }
#' 
#' 
apoisson <- function(model, data, 
                      init = tab2vec(data), 
                      iter = 1E4, burn = 1000, 
                      thin = 10, engine = c("Cpp","R"), 
                      method = c("ipf", "mcmc"), moves, 
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
  
  data   <- suppressMessages(teshape(data, "tab"))
  p      <- length(dim(data))
  nCells <- length(data)
  
  ## if a pure array is given, give names for later
  if(is.array(data) && is.null(dimnames(data))) data <- array2tab(data)
  
  ## other basic objects
  varsNlevels <- dimnames(data)  
  vars        <- names(varsNlevels)
  
  
  
  
  ## check for sampling zeros
  ##################################################
  if(any(data == 0L)) message(
    "Care ought be taken with tables with sampling zeros to ensure the MLE exists."
  )


  ## parse model specification (formula for vector of r_k's)
  ##################################################
  
  # assume that the model given is the vector
  # of variable levels


  ## construct A matrix and compute moves
  ##################################################  
  
  if(missing(moves) && !is.null(getOption("4ti2_path"))){
   
    randomly_generate_moves <- TRUE
    
  } else if(missing(moves) && is.null(getOption("4ti2_path"))){
    
    randomly_generate_moves <- TRUE
    
  } else if(is.character(moves)){

    stop("not yet implemented.")    
    
  }
  


  ## run metropolis-hastings
  ##################################################  
  init <- unname(init) # init
  # out <- metropolis(init, moves, iter = iter, burn = burn, thin = thin, engine = engine)  
  out <- metropolis(init, moves, iter = iter, burn = burn, thin = thin, engine = engine)  


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
  u <- t(t(unname(tab2vec(data))))
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
  out$exp <- exp
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

  class(out) <- "poisson"
  out
}




