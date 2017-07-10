sis_table <- function(config_mat, suff_statistics){
  #Need to check if config_mat is a mat, suff_statistics is vector/mat, 
  #length of suff_statistics is same as number of rows of config_mat
  
  #Matricies and vectors to work with!
  work_A <- config_mat
  work_suff <- suff_statistics
  tbl_elts <- ncol(config_mat)
  num_const <- nrow(config_mat)
  tbl <- vector(mode = "numeric", length = tbl_elts)
  
  for(i in 1:tbl_elts){
    constr <- unname(rbind(cbind(rep(1,num_const), work_suff, work_A),
                           cbind(rep(0,tbl_elts), rep(0,tbl_elts), diag(-1,tbl_elts))))
    objfun <- vector(mode = "numeric", length = tbl_elts)
    objfun[i] <- -1
    min_lp <- lpcdd(constr, objfun)
    max_lp <- lpcdd(constr, objfun, minimize = FALSE)
    
    if(min_lp[1] == "Optimal" && max_lp[1] == "Optimal"){
      print(minimum <- as.numeric(unname(min_lp[4])))
      print(maximum <- as.numeric(unname(max_lp[4])))
      tbl[i] <- if(isTRUE(all.equal(minimum, maximum))){minimum}
      else{sample(minimum:maximum, 1)}
    } else { tbl[i] <- 0 }
    #Update constraints and sufficient statistics
    index <- which(work_A[,i] == 1)
    work_A[index,i] <- 0
    work_suff[index] <- work_suff[index] - tbl[i]
  }
  return(tbl)
}