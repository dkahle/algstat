.onAttach <- function(...) {
  
  packageStartupMessage('  Please cite algstat! See citation("algstat") for details.')  
  
  if(is.mac()){ ## find the path on a mac	

    mac_search_and_set("bertini", "bertini", "bertini")
    
  } else if(is.win()){ ## find the path on a pc  	
  	
  	if(!any(str_detect(tolower(list.files("C:\\")), "cygwin"))){
  	  psm("  Cygwin is required to run most of algstat on a Windows platform.")
  	  psm("  It needs to be in your C:\\ drive, but wasn't found.")  	  
  	  return(invisible())
  	}
    
    if(!whereis_is_accessible()){ # check for whereis, return if not found
      psm(
        "  The whereis function was not found, so algstat can't find the required exe's.\n",
        "  Try setting the path with setBertiniPath()."
      )
      return()
    }
    

  	win_search_and_set("bertini")
    
    
  } #else {  ## find the path on a unix-type machine
    # note that this is done after os-x with an else
  	
    
    
  #}
  
  ## check for programs
  startup_check_for_program("bertini")    
  
}













.onUnload <- function (libpath) {
  library.dynam.unload("algstat", libpath)
}













# each program (suite of executables) is identified by four names
# optionName : the program name (e.g. latte)
# longName   : the program pretty printed (e.g. LattE)
# execName   : the defining executable of that package (e.g. count); 
#              this is what algstat looks for when looking for the program
# setFun     : the helper function to set that path (e.g. setLatteFun)


longName <- function(optionName){
  switch(optionName,
    bertini = "Bertini", 
    "4ti2" = "4ti2", 
    latte = "LattE"
  )
}

execName <- function(optionName){
  switch(optionName,
    bertini = "bertini", 
    "4ti2" = "markov", 
    latte = "count"
  )
}

setFun <- function(optionName){
  switch(optionName,
    bertini = "setBertiniPath()", 
    "4ti2" = "set4ti2Path()", 
    latte = "setLattePath()"
  )
}










psm  <- packageStartupMessage
psms <- function(fmt, ...) packageStartupMessage(sprintf(fmt, ...))



startup_check_for_program <- function(optionName){
  
  longName <- longName(optionName)
  setFun <- setFun(optionName)
  
  if(is.null(getOption(optionName))){
    psms("  %s not found. Set the location with %s", longName, setFun)
    return(invisible(FALSE))
  }
    
  if(length(list.files(getOption(optionName))) == 0){
    psms("  %s appears to be installed, but it's not where it was expected.", longName)
    psms("  Suggestion : run %s", setFun)
    return(invisible(FALSE))    
  }	
    
  invisible(TRUE)
  
}





program_not_found_stop <- function(optionName){
  
  longName <- longName(optionName)
  setFun <- setFun(optionName)
  
  if(is.null(getOption(optionName))){
    stop(sprintf("  %s not found. Set the location with %s", longName, setFun))
    return(invisible(FALSE))
  }
  
}








setOption <- function(optionName, value){
  eval(parse(text = sprintf('options("%s" = "%s")', optionName, value)))
}












whereis_is_accessible <- function() unname(Sys.which("whereis")) != ""

win_find <- function(s){ 
  wexe <- unname(Sys.which("whereis"))
  x <- system(paste(wexe, s), TRUE)
  str_sub(x, nchar(s)+2)
}

win_search_and_set <- function(optionName){
  
  # search
  x <- win_find(execName(optionName))
  if(str_detect(x, "/")) setOption(optionName, dirname(x))

}














# this seems too slow to load every time, so the below wraps it
# it reduces the search space where this function is then used
mac_find <- function(exec, where){
  
  # query the system and clean attributes
  query <- sprintf("find %s -name %s", where, exec)
  finding <- suppressWarnings(system(query, intern = TRUE, ignore.stderr = TRUE))
  attributes(finding) <- NULL
  
  # get the bin first
  path <- finding[str_detect(finding, paste0("bin/", exec))][1]
  
  # bertini isn't in a bin directory
  if(is.na(path)) path <- finding[1]
  
  # return
  path
}


mac_search_and_set <- function(exec, baseName, optionName){
  
  # look in applications and home
  apps <- list.files("/Applications")
  home <- list.files("~")

  # check for main dir name
  ndxApps  <- which(str_detect(tolower(apps), baseName))
  ndxHome  <- which(str_detect(tolower(home), baseName))   
  
  if(length(ndxApps) == 0 && length(ndxHome) == 0){
    return(NULL)
  } else if(length(ndxApps) != 0){
    path <- mac_find(exec, paste0("/Applications/", apps[ndxApps]))
  } else if(length(ndxHome) != 0){
    path <- mac_find(exec, paste0("~/", apps[ndxApps]))
  }  
  
  setOption(optionName, dirname(path))
  
}














