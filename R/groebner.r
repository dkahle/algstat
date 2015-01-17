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
#' \dontrun{
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
#' groebner(A, "tab", varlvls, TRUE)
#'
#'
#'
#' 
#' # 3x3 independence example
#' # following convention, the first index indicates rows
#' varlvls <- c(3,3)
#' facets <- list(1,2)
#' ( A <- hmat(varlvls, facets) )
#' groebner(A)
#' groebner(A, "vec")
#' groebner(A, "tab", varlvls)
#' groebner(A, "tab", varlvls, TRUE)
#' 
#' 
#' 
#' 
#' # LAS example 1.2.1, p.12 (2x3 independence)
#' varlvls <- c(2,3)
#' facets <- list(1, 2)
#' ( A <- hmat(varlvls, facets) )
#' groebner(A, "tab", varlvls)
#' # Prop 1.2.2 says that there should be 
#' 2*choose(2, 2)*choose(3,2) # = 6
#' # moves.
#' groebner(A, "tab", varlvls, TRUE)
#' 
#'
#' 
#' 
#' 
#' # LAS example 1.2.12, p.17  (no 3-way interaction)
#' varlvls <- c(2,2,2)
#' facets <- list(c(1,2), c(1,3), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#' groebner(A)
#' tableau(groebner(A), dim = varlvls)
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
#' groebner(A)
#' groebner(A, "tab", varlvls) # hard to understand
#' tableau(groebner(A), varlvls)
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
#'
#' 
#'
#'
#' 
#' 
#'
#' 
#' groebner(diag(1, 10))
#'
#' }
#' 
groebner <- function(mat, format = c("mat", "vec", "tab"), dim = NULL,
                   all = FALSE, dir = tempdir(), opts = "-parb", quiet = TRUE
){
  
  format <- match.arg(format)
  
  ## redirect with the special case of when the identity is given
  if((nrow(mat) == ncol(mat)) && all(mat == diag(1, nrow(mat)))){
    warning("the identity matrix was supplied to groebner, returning 0's.")
    if(format == "mat") return(matrix(0, nrow = nrow(mat), ncol = 1))
    if(format == "vec") return(list(matrix(0, nrow = nrow(mat), ncol = 1)))
    if(format == "tab") return(vec2tab(matrix(0, nrow = nrow(mat), ncol = 1), dim))    
  }
  
  ## make dir to put 4ti2 files in (within the tempdir) timestamped
  dir2 <- file.path2(dir, timeStamp())
  suppressWarnings(dir.create(dir2))
  
  
  ## define a function to write the code to a file
  formatAndWriteFile <- function(mat, codeFile = "groebnerCode.mat"){
    
    # line numbers up, e.g. 1 1 0 0
    out <- paste(nrow(mat), ncol(mat))
    out <- paste0(out, "\n")
    out <- paste0(out, 
                  paste(apply(unname(mat), 1, paste, collapse = " "), 
                        collapse = "\n")
    )    
    
    # write code file
    writeLines(out, con = file.path2(dir2, codeFile))
    invisible(out)
  }	
  
  
  ## make 4ti2 file
  if(!missing(mat)) formatAndWriteFile(mat)
  
  
  ## switch to temporary directory
  oldWd <- getwd()
  setwd(dir2)
  
  
  ## compute basis with 4ti2
  if(.Platform$OS.type == "unix"){
    
    system2(
      file.path2(getOption("markovPath"), "groebner"),
      paste(opts, file.path2(dir2, "groebnerCode.mat")),
      stdout = "groebnerOut", stderr = FALSE
    )
    
  } else { # windows 
    
    matFile <- file.path2(dir2, "groebnerCode.mat")
    matFile <- chartr("\\", "/", matFile)
    matFile <- str_c("/cygdrive/c", str_sub(matFile, 3)) 
    
    system2(
      "cmd.exe",
      paste(
        "/c env.exe", 
        file.path(getOption("groebnerPath"), "groebner"), 
        opts, matFile
      ), stdout = "groebnerOut", stderr = FALSE
    )
    
  }
  
  
  if(!quiet) cat(readLines("groebnerOut"), sep = "\n")
  
  
  ## figure out what files to keep them, and make 4ti2 object
  basis <- readLines(paste0("groebnerCode.mat", ".gro"))
  basis <- basis[-1]
  basis <- lapply(
    str_split(str_trim(basis), " "), 
    function(x) as.integer(x[nchar(x) > 0])
  )
  if(all) basis <- c(basis, lapply(basis, function(x) -x))
  
  
  ## migrate back to original working directory
  setwd(oldWd)
  
  
  
  # out
  if(format == "mat"){
    basis <- matrix(unlist(basis), ncol = length(basis[[1]]), byrow = TRUE) 
    return(t(basis))
  } else if(format == "vec") {
    return(basis)
  } else { # format == "tab"
    return(lapply(basis, vec2tab, dim = dim))
  }
}




