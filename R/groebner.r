#' Compute a Grobner basis with 4ti2
#'
#' A Grobner basis of a matrix A is computed with the groebner function of 4ti2, obtained with the LattE-integrale bundle.
#' 
#' @param mat a matrix; for example the output of \code{\link{hmat}}
#' @param format how the moves should be returned (if "mat", moves are columns)
#' @param dim the dimension to be used in vec2tab if format = "tab" is used, oftentimes a vector of the number of levels of each variable in order
#' @param all if TRUE, all moves (+ and -) are given.  if FALSE, only the + moves are given
#' @param dir directory to place the files in, without an ending /
#' @param opts options for groebner
#' @param quiet show 4ti2 output
#' @return a matrix containing the Grobner basis as its columns (for easy addition to tables)
#' @seealso \code{\link{markov}}
#' @export groebner
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @examples
#' 
#' \dontrun{
#' # these examples require having extra software installed
#'
#'
#'
#' 
#' # 2x2 independence example
#' # following convention, the first index indicates rows
#' varlvls <- c(2,2)
#' facets <- list(1,2)
#' ( A <- hmat(varlvls, facets) )
#' groebner(A)
#' groebner(A, "vec")
#' groebner(A, "tab", varlvls)
#' groebner(A, "tab", varlvls, all = TRUE)
#' groebner(A, quiet = FALSE)
#'
#'
#' }
#' 
groebner <- function(mat, format = c("mat", "vec", "tab"), dim = NULL,
  all = FALSE, dir = tempdir(), opts = "-parb", quiet = TRUE
){
  
  ## check for 4ti2
  if(is.null(getOption("markovPath"))){
    stop(
      "algstat doesn't know where groebner is (or any other 4ti2 programs),\n",
      "  and so can't compute a groebner basis.  see ?setMarkovPath", call. = FALSE
    )
  }
  
  
  ## check args
  format <- match.arg(format)
  if(format == "tab" && missing(dim)){
    stop('if format = "tab" is specified, dim must be also.', call. = FALSE) 
  }
  
  
  ## make dir to put 4ti2 files in (within the tempdir) timestamped
  dir2 <- file.path2(dir, timeStamp())
  suppressWarnings(dir.create(dir2))
  
  
  ## make 4ti2 file
  if(!missing(mat)) write.latte(mat, file.path2(dir2, "groebnerCode.mat"))
  
  
  ## switch to temporary directory
  oldWd <- getwd()
  setwd(dir2)
  on.exit(setwd(oldWd), add = TRUE)
  
    
  ## run 4ti2
  if(is.mac() || is.unix()){
      
    system2(
      file.path2(getOption("markovPath"), "groebner"),
      paste(opts, file.path2(dir2, "groebnerCode.mat")),
      stdout = "groebnerOut", stderr = FALSE
    )
      
  } else if(is.win()){ 
      
    matFile <- file.path2(dir2, "groebnerCode.mat")
    matFile <- chartr("\\", "/", matFile)
    matFile <- str_c("/cygdrive/c", str_sub(matFile, 3)) 
      
    system2(
      "cmd.exe",
      paste(
        "/c env.exe", 
        file.path(getOption("markovPath"), "groebner"), 
        opts, matFile
      ), stdout = "groebnerOut", stderr = FALSE
    )
      
  }
    
  
  ## print if requested
  if(!quiet) cat(readLines("groebnerOut"), sep = "\n")

  
  ## figure out what files to keep them, and make 4ti2 object
  basis <- t(read.latte(paste0("groebnerCode.mat", ".gro")))
  
  
  ## fix case of no basis
  basisDim <- dim(basis)
  noBasisFlag <- FALSE
  if(any(basisDim == 0)){
    noBasisFlag <- TRUE
    warning("groebner basis empty, returning 0's.", call. = FALSE)
    basisDim[basisDim == 0] <- 1L
    basis <- rep(0L, prod(basisDim))
    dim(basis) <- basisDim    
  }
  
  
  ## format
  if(all && !noBasisFlag) basis <- cbind(basis, -basis)
  
  
  # out
  if(format == "mat"){
    return(basis)
  } else {        
    lbasis <- as.list(rep(NA, ncol(basis)))
    for(k in 1:ncol(basis)) lbasis[[k]] <- basis[,k]
    if(format == "vec") return(lbasis)     
    if(format == "tab") return(lapply(lbasis, vec2tab, dim = dim))    
  }
}




