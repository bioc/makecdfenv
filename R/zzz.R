##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib <- function(libname, pkgname, where) {
  ## require(affy, quietly=TRUE) || stop("Cannot load without package \"affy\"")  
  library.dynam("makecdfenv", pkgname, libname)
}

