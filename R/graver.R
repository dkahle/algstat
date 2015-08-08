#' Compute a Graver basis with 4ti2
#' 
#' A Graver basis of a matrix A is computed with the graver function of 4ti2,
#' obtained with the LattE-integrale bundle.
#' 
#' @param mat a matrix; for example the output of \code{\link{hmat}}
#' @param format how the moves should be returned (if "mat", moves are columns)
#' @param dim the dimension to be used in vec2tab if format = "tab" is used,
#'   oftentimes a vector of the number of levels of each variable in order
#' @param all if TRUE, all moves (+ and -) are given.  if FALSE, only the +
#'   moves are given
#' @param dir directory to place the files in, without an ending /
#' @param opts options for graver
#' @param quiet show 4ti2 output
#' @param dbName the name of the model in the Markov bases database,
#'   http://markov-bases.de, see examples
#' @return a matrix containing the Graver basis as its columns (for easy
#'   addition to tables)
#' @export
#' @name graver
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures
#'   on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
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
#' graver(A)
#' graver(A, "vec")
#' graver(A, "tab", varlvls)
#' graver(A, "tab", varlvls, all = TRUE)
#' 
#' 
#' 
#' 
#' # 3x3 independence example
#' # following convention, the first index indicates rows
#' varlvls <- c(3,3)
#' facets <- list(1,2)
#' ( A <- hmat(varlvls, facets) )
#' markov(A)
#' graver(A)
#' graver(A, "vec")
#' graver(A, "tab", varlvls)
#' graver(A, "tab", varlvls, TRUE)
#' 
#' 
#' 
#' 
#' 
#' 
#' # markov basis vs graver basis
#' varlvls <- c(3,3,3)
#' facets <- list(c(1,2), c(1,3), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#' str(zbasis(A)) #   8 elements
#' str(markov(A)) #  81 elements
#' str(graver(A)) # 795 elements
#' 
#' 
#' 
#' 
#' 
#' 
#' # LAS example 1.2.12, p.16  (no 3-way interaction)
#' varlvls <- c(2,2,2,2)
#' facets <- list(c(1,2), c(1,4), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#' graver(A)
#' graver(A, "tab", varlvls) # hard to understand
#' # tableau(graver(A), varlvls) # error?
#' 
#' 
#' 
#' 
#' 
#' # memoising graver - faster recomputing times
#' # (memorizing computed graver bases)
#' A <- hmat(c(3,3,2), 1:3)
#' system.time(graver(A))
#' system.time(graver(A))
#' system.time(memGraver(A))
#' system.time(memGraver(A))
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
#' markov(diag(1, 10))
#' graver(diag(1, 10))
#' graver(diag(1, 10), "vec")
#' graver(diag(1, 10), "vec", all = TRUE)
#' 
#' }
#' 
graver <- function(mat, format = c("mat", "vec", "tab"), dim = NULL,
  all = FALSE, dir = tempdir(), opts = "-p=gmp", quiet = TRUE,
  dbName
){
  
  ## check for 4ti2
  if(is.null(getOption("markovPath"))){
    stop(
      "algstat doesn't know where graver is (or any other 4ti2 programs),\n",
      "  and so can't compute a Graver basis.  see ?setMarkovPath", call. = FALSE
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
  if(!missing(mat)) write.latte(mat, file.path2(dir2, "PROJECT.mat"))

  
  ## switch to temporary directory
  oldWd <- getwd()
  setwd(dir2)
  on.exit(setwd(oldWd), add = TRUE)
  
  
  ## create/retrieve graver basis
  if(missing(dbName)){
    
    ## run 4ti2 if needed
    if(is.mac() || is.unix()){
      
      system2(
        file.path2(getOption("markovPath"), "graver"),
        paste(opts, file.path2(dir2, "PROJECT")),
        stdout = "graverOut", stderr = FALSE
      )
      
    } else if(is.win()){ 
      
      matFile <- file.path2(dir2, "PROJECT")
      matFile <- chartr("\\", "/", matFile)
      matFile <- str_c("/cygdrive/c", str_sub(matFile, 3)) 
      
      system2(
        "cmd.exe",
        paste(
          "/c env.exe", 
          file.path(getOption("markovPath"), "graver"), 
          opts, matFile
        ), stdout = "graverOut", stderr = FALSE
      )
      
    }
    
    
    if(!quiet) cat(readLines("graverOut"), sep = "\n")
    
  } else { # if the model name is specified
    
    download.file(
      paste0("http://markov-bases.de/data/", dbName, "/", dbName, ".gra"),
      destfile = "PROJECT.gra" # already in tempdir
    )
    
  }
  
  
  ## figure out what files to keep them, and make 4ti2 object
  basis <- t(read.latte("PROJECT.gra"))
  
  
  ## fix case of no basis
  basisDim <- dim(basis)
  noBasisFlag <- FALSE
  if(any(basisDim == 0)){
    noBasisFlag <- TRUE
    warning("graver basis empty, returning 0's.", call. = FALSE)
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






#' @param ... ...
#' @export
#' @rdname graver
memGraver <- memoise::memoise(graver)






