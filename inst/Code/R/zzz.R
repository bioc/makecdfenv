.First.lib <- function(libname, pkgname){
      path = .path.package(pkgname)
      where <- as.environment(match(paste("package:", pkgname, sep = ""),search()))
      data(list="@PKGNAME@", package=pkgname, envir = where)
}
