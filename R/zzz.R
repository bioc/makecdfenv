##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib <- function(libname, pkgname, where) {
  require("utils")
  library.dynam("makecdfenv", pkgname, libname)
}

