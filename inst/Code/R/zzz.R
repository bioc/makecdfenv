.onLoad <- function(libname, pkgname){
      where = asNamespace(pkgname)
      data(list = pkgname, package = pkgname, envir = where)
      packageStartupMessage(AnnotationDbi:::annoStartupMessages(pkgname))
}
