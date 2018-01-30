.onAttach <- function(...) {
  
  packageStartupMessage('Please cite algstat! See citation("algstat") for details.')  
  
}


.onUnload <- function (libpath) {
  library.dynam.unload("algstat", libpath)
}

