#' Compute a Markov basis with 4ti2
#'
#' A Markov basis of a matrix A is computed with the markov function of 4ti2, obtained with the LattE-integrale bundle.
#' 
#' @param mat a matrix; for example the output of \code{\link{hmat}}
#' @param format how the moves should be returned (if "mat", moves are columns)
#' @param dim the dimension to be used in vec2tab if format = "tab" is used, oftentimes a vector of the number of levels of each variable in order
#' @param all if TRUE, all moves (+ and -) are given.  if FALSE, only the + moves are given
#' @param dir directory to place the files in, without an ending /
#' @param opts options for markov
#' @param quiet show 4ti2 output
#' @param dbName the name of the model in the markov bases database, http://markov-bases.de, see examples
#' @return a matrix containing the Markov basis as its columns (for easy addition to tables)
#' @export 
#' @name markov
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
#' markov(A)
#' markov(A, "vec")
#' markov(A, "tab", varlvls)
#' markov(A, "tab", varlvls, all = TRUE)
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
#' markov(A, "vec")
#' markov(A, "tab", varlvls)
#' markov(A, "tab", varlvls, TRUE)
#' 
#' 
#' 
#' 
#' # LAS example 1.2.1, p.12 (2x3 independence)
#' varlvls <- c(2,3)
#' facets <- list(1, 2)
#' ( A <- hmat(varlvls, facets) )
#' markov(A, "tab", varlvls)
#' # Prop 1.2.2 says that there should be 
#' 2*choose(2, 2)*choose(3,2) # = 6
#' # moves.
#' markov(A, "tab", varlvls, TRUE)
#' 
#'
#' 
#' 
#' 
#' # LAS example 1.2.12, p.17  (no 3-way interaction)
#' varlvls <- c(2,2,2)
#' facets <- list(c(1,2), c(1,3), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#' markov(A)
#' tableau(markov(A), dim = varlvls)
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
#' markov(A)
#' markov(A, "tab", varlvls) # hard to understand
#' tableau(markov(A), varlvls)
#' 
#' 
#'
#' 
#' 
#' # memoising markov - faster recomputing times
#' # (memorizing computed markov bases)
#' A <- hmat(c(4,4,4), 1:3)
#' system.time(markov(A))
#' system.time(markov(A))
#' system.time(memMarkov(A))
#' system.time(memMarkov(A))
#' 
#' 
#' 
#'
#' 
#' # using the markov bases database, must be connected to internet
#' # A <- markov(dbName = "ind3-3")
#' B <- markov(hmat(c(3,3), list(1,2)))
#' # all(A == B)
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
#' markov(diag(1, 10), "vec")
#' markov(diag(1, 10), "vec", all = TRUE)
#'
#' }
#' 
markov <- function(mat, format = c("mat", "vec", "tab"), dim = NULL,
  all = FALSE, dir = tempdir(), opts = "-parb", quiet = TRUE,
  dbName
){
  
  ## check for 4ti2
  if(is.null(getOption("markovPath"))){
    stop(
      "algstat doesn't know where markov is (or any other 4ti2 programs),\n",
      "  and so can't compute a markov basis.  see ?setMarkovPath", call. = FALSE
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
  if(!missing(mat)) write.latte(mat, file.path2(dir2, "markovCode.mat"))


  ## switch to temporary directory
  oldWd <- getwd()
  setwd(dir2)
  on.exit(setwd(oldWd), add = TRUE)
  
  
  ## create/retrieve markov basis
  if(missing(dbName)){
  	
    ## run 4ti2 if needed
    if(is.mac() || is.unix()){
      
      system2(
        file.path2(getOption("markovPath"), "markov"),
        paste(opts, file.path2(dir2, "markovCode.mat")),
        stdout = "markovOut", stderr = FALSE
      )
      
    } else if(is.win()){ 
      
      matFile <- file.path2(dir2, "markovCode.mat")
      matFile <- chartr("\\", "/", matFile)
      matFile <- str_c("/cygdrive/c", str_sub(matFile, 3)) 
      
      system2(
        "cmd.exe",
        paste(
          "/c env.exe", 
          file.path(getOption("markovPath"), "markov"), 
          opts, matFile
        ), stdout = "markovOut", stderr = FALSE
      )
      
    }

    
    if(!quiet) cat(readLines("markovOut"), sep = "\n")
    
  } else { # if the model name is specified
  	
    download.file(
      paste0("http://markov-bases.de/data/", dbName, "/", dbName, ".mar"),
      destfile = "markovCode.mat.mar" # already in tempdir
    )
    
  }
  
  
  ## figure out what files to keep them, and make 4ti2 object
  basis <- t(read.latte(paste0("markovCode.mat", ".mar")))
  
  
  ## fix case of no basis
  basisDim <- dim(basis)
  noBasisFlag <- FALSE
  if(any(basisDim == 0)){
    noBasisFlag <- TRUE
    warning("markov basis empty, returning 0's.", call. = FALSE)
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







#' @export
#' @rdname markov
memMarkov <- memoise::memoise(markov)






