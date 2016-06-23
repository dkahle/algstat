#' Set paths to external functions
#' 
#' These functions set the path to external programs either by (1) 
#' passing them a character string or (2) using 
#' \code{\link{file.choose}}.
#' 
#' @param path a character string, the path to bertini (for
#'   Bertini), markov (say, for 4ti2), and count (say, for LattE)
#' @return an invisible character string, the path found
#' @name setPaths
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @examples
#' 
#' \dontrun{ # the below code requires suggested external software
#' 
#' 
#' ## for Bertini
#' getOption("bertini")
#' setBertiniPath() 
#' 
#' ## for LattE
#' getOption("latte")
#' setLattePath() 
#' 
#' ## for 4ti2 (typically the same as LattE)
#' getOption("4ti2")
#' set4ti2Path() 
#' 
#' 
#' 
#' ## each of these functions can be used statically as well
#' (`4ti2Path` <- getOption("4ti2"))
#' set4ti2Path("/path/to/4ti2/directory")
#' getOption("4ti2")
#' set4ti2Path(`4ti2Path`) # undoes example
#' 
#' 
#' 
#' }
#' 
NULL













#' @rdname setPaths
#' @export 
setBertiniPath <- function(path){

  if(missing(path) && interactive()){
  	
    bertiniPath <- dirname(file.choose())
    if(is.win() && str_detect(bertiniPath,"C:/")){
      bertiniPath <- str_replace(bertiniPath, "C:/", "/cygdrive/c/")
    }    
    options(bertini = bertiniPath)
    return(invisible(bertiniPath))
    
  } else if(!missing(path)){
  	
    options(bertini = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}











#' @rdname setPaths
#' @aliases setLattEPath
#' @export 
setLattePath <- function(path){

  if(missing(path) && interactive()){
  	
    lattePath <- dirname(file.choose())
    if(is.win() && str_detect(lattePath,"C:/")){
      lattePath <- str_replace(dirname(lattePath), "C:/", "/cygdrive/c/")
    }    
    options(latte = lattePath)
    return(invisible(lattePath))
    
  } else if(!missing(path)){
  	
    options(latte = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}


#' @export 
setLattEPath <- setLattePath











#' @rdname setPaths
#' @export 
set4ti2Path <- function(path){

  if(missing(path) && interactive()){

    `4ti2_path` <- dirname(file.choose())
    if(is.win() && str_detect(`4ti2_path`,"C:/")){
      `4ti2_path` <- str_replace(`4ti2_path`, "C:/", "/cygdrive/c/")
    }
    options(`4ti2` = `4ti2_path`)
    return(invisible(`4ti2_path`))

  } else if(!missing(path)){

    options(`4ti2` = path)
    return(invisible(path))

  } else {
    stop(
      "If the session is not interactive, a path must be specified.",
      call. = FALSE
    )
  }
  
}
