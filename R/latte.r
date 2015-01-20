#' Write a latte file to disk
#'
#' \code{write.latte} writes a latte-formatted file a file.
#' 
#' @param mat a matrix; for example the output of \code{\link{hmat}}
#' @param file the name of the file
#' @return an invisible form of the saved output.
#' @export write.latte
#' @seealso \code{\link{write.latte}}
#' @examples
#' \dontrun{
#'
#' mat <- matrix(sample(9), 3, 3)
#' mat
#' write.latte(mat, "foo.hrep")
#' file.show("foo.hrep")
#' read.latte("foo.hrep")
#' read.latte("foo.hrep", "Ab")
#' 
#' attr(mat, "linearity") <- c(1, 3)
#' attr(mat, "nonnegative") <- 2
#' mat
#' write.latte(mat, "foo.hrep")
#' file.show("foo.hrep")
#' read.latte("foo.hrep")
#' 
#' file.remove("foo.hrep")
#'
#' }
#' 
write.latte <- function(mat, file){  
  
  ## construct file in latte format
  ## e.g. "3 3\n5 6 8\n7 3 4\n2 9 1"
  ## and cat("3 3\n5 6 8\n7 3 4\n2 9 1")
  out <- paste(nrow(mat), ncol(mat))
  out <- paste0(out, "\n")
  out <- paste0(out, 
    paste(
      apply(unname(mat), 1, paste, collapse = " "), 
      collapse = "\n"
    )
  )
  
  
  ## add linearity lines, if present
  if("linearity" %in% names(attributes(mat))){
    linLines <- attr(mat, "linearity")
    linLineToAdd <- paste(
      "linearity", 
      length(linLines), 
      paste(linLines, collapse = " ")
    )
    out <- paste(out, linLineToAdd, sep = "\n")
  }
  
  
  ## add nonnegative lines, if present
  if("nonnegative" %in% names(attributes(mat))){
    nnegLines <- attr(mat, "nonnegative")
    nnegLineToAdd <- paste(
      "nonnegative", 
      length(nnegLines), 
      paste(nnegLines, collapse = " ")
    )
    out <- paste(out, nnegLineToAdd, sep = "\n")
  }
  
  
  ## save it to disk
  writeLines(out, con = file)
  
  
  ## return invisible output
  invisible(out)
}






















#' Read a latte file from disk
#'
#' \code{read.latte} reads a latte-formatted file a file.
#' 
#' @param file the name of the file
#' @param format the format of the read file, as a latte-style matrix [b -A] 
#'   or as matrix A and a solution b
#' @return an invisible form of the saved output.
#' @export read.latte
#' @seealso \code{\link{read.latte}}
#' @examples
#' \dontrun{
#'
#' mat <- matrix(sample(9), 3, 3)
#' mat
#' write.latte(mat, "foo.hrep")
#' file.show("foo.hrep")
#' read.latte("foo.hrep")
#' read.latte("foo.hrep", "Ab")
#' 
#' attr(mat, "linearity") <- c(1, 3)
#' attr(mat, "nonnegative") <- 2
#' mat
#' write.latte(mat, "foo.hrep")
#' file.show("foo.hrep")
#' read.latte("foo.hrep")
#' 
#' file.remove("foo.hrep")
#'
#' }
#' 
read.latte <- function(file, format = c("mat", "Ab")){  
  
  ## check args
  format <- match.arg(format)
  
  
  ## read in file
  ## e.g. [1] "3 3"   "7 2 1" "6 3 4" "9 8 5"
  contents <- readLines(file)
  
  
  ## eliminate dimensions
  ## e.g. [1] "7 2 1" "6 3 4" "9 8 5"
  dim <- as.integer(str_split(str_trim(contents[1]), " ")[[1]])
  
  
  ## split and parse, result is list
  matRows <- lapply(
    str_split(str_trim(contents[2:(1+dim[1])]), " "), 
    function(x) as.integer(x[nchar(x) > 0])
  )
  
  
  ## put into matrix
  mat <- matrix(unlist(matRows), nrow = dim[1], byrow = TRUE)
  
  
  ## check for linearity or nonnegative lines
  if(length(contents) > 1 + dim[1]){
    
    # isolate the added lines
    addedLines <- contents[-(1:(1 + dim[1]))]
    
    # look for linearity, assume only one such line
    linQ <- str_detect(addedLines, "linearity")
    if(any(linQ)){
      linLine  <- addedLines[linQ]
      linStuff <- str_replace(linLine, "linearity ", "")
      linNdcs  <- as.integer(str_split(linStuff, " ")[[1]][-1])
      attr(mat, "linearity") <- linNdcs
    }
    
    # look for nonnegative, assume only one such line
    nnegQ <- str_detect(addedLines, "nonnegative")
    if(any(nnegQ)){
      nnegLine  <- addedLines[nnegQ]
      nnegStuff <- str_replace(nnegLine, "nonnegative ", "")
      nnegNdcs  <- as.integer(str_split(nnegStuff, " ")[[1]][-1])
      attr(mat, "nonnegative") <- nnegNdcs
    }
    
  }
  
  
  ## format
  if(format == "mat"){
    return(mat)
  } else if(format == "Ab"){
    return(list(A = -mat[,-1,drop=FALSE], b = mat[,1]))
  }
}





