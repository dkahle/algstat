sis_table <- function(configMat, suffStatistics){
  # need to check if configMat is a mat, suffStatistics is vector/mat, 
  # length of suffStatistics is same as number of rows of configMat
  
  # matricies and vectors to work with!
  workA <- configMat
  workSuff <- suffStatistics
  tblElts <- ncol(configMat)
  numConstraints <- nrow(configMat)
  tbl <- vector(mode = "numeric", length = tblElts)
  
  for(i in 1:tblElts){
    constr <- unname(rbind(cbind(rep(1,numConstraints), workSuff, workA),
                           cbind(rep(0,tblElts), rep(0,tblElts), diag(-1,tblElts))))
    objfun <- vector(mode = "numeric", length = tblElts)
    objfun[i] <- -1
    minLp <- lpcdd(constr, objfun)
    maxLp <- lpcdd(constr, objfun, minimize = FALSE)
    
    if(minLp[1] == "Optimal" && maxLp[1] == "Optimal"){
      minimum <- as.numeric(unname(minLp[4]))
      maximum <- as.numeric(unname(maxLp[4]))
      tbl[i] <- if(isTRUE(all.equal(minimum, maximum))){minimum}
      else{sample(minimum:maximum, 1)}
    } else { tbl[i] <- 0 }
    # update constraints and sufficient statistics
    index <- which(workA[,i] == 1)
    workA[index,i] <- 0
    workSuff[index] <- workSuff[index] - tbl[i]
  }
  return(tbl)
}