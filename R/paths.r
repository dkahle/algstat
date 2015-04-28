#' Set paths to external functions
#' 
#' These functions set the path to external programs either by (1) passing them a character string or (2) using \code{\link{file.choose}}. 
#' 
#' @param path a character string, the path to m2 (for Macaulay2), bertini (for Bertini), markov (say, for 4ti2), and count (say, for LattE)
#' @return an invisible character string, the path found
#' @name setPaths
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @examples
#' 
#' \dontrun{ # the below code requires suggested external software
#' 
#' ## for Macaulay2
#' getOption("m2Path")
#' setM2Path() 
#' 
#' ## for Bertini
#' getOption("bertiniPath")
#' setBertiniPath() 
#' 
#' ## for LattE
#' getOption("lattePath")
#' setLattePath() 
#' 
#' ## for 4ti2 (typically the same as LattE)
#' getOption("markovPath")
#' setMarkovPath() 
#' 
#' 
#' 
#' ## each of these functions can be used statically as well
#' (markovPath <- getOption("markovPath"))
#' setMarkovPath("/path/to/4ti2/directory")
#' getOption("markovPath")
#' setMarkovPath(markovPath) # undoes example
#' 
#' 
#' 
#' }
#' 
NULL








#' @rdname setPaths
#' @export 
setM2Path <- function(path){

  if(missing(path) && interactive()){
  	
    m2Path <- dirname(file.choose())
    if(is.win() && str_detect(m2Path,"C:/")){
      m2Path <- str_replace(m2Path, "C:/", "/cygdrive/c/")
    }
    options(m2Path = m2Path)
    return(invisible(m2Path))
    
  } else if(!missing(path)){
  	
    options(m2Path = path)
    return(invisible(path))
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}







#' @rdname setPaths
#' @export 
setBertiniPath <- function(path){

  if(missing(path) && interactive()){
  	
    bertiniPath <- dirname(file.choose())
    if(is.win() && str_detect(bertiniPath,"C:/")){
      bertiniPath <- str_replace(bertiniPath, "C:/", "/cygdrive/c/")
    }    
    options(bertiniPath = bertiniPath)
    return(invisible(bertiniPath))
    
  } else if(!missing(path)){
  	
    options(bertiniPath = path)
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
    options(lattePath = lattePath)
    return(invisible(lattePath))
    
  } else if(!missing(path)){
  	
    options(lattePath = path)
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
setMarkovPath <- function(path){

  if(missing(path) && interactive()){
  	
    markovPath <- dirname(file.choose())
    if(is.win() && str_detect(markovPath,"C:/")){
      markovPath <- str_replace(markovPath, "C:/", "/cygdrive/c/")
    }        
    options(markovPath = markovPath)
    return(invisible(markovPath))
    
  } else if(!missing(path)){
  	
    options(markovPath = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}
