.First.lib <- function(libname, pkgname){
      path = .path.package(pkgname)
      where <- as.environment(match(paste("package:", pkgname, sep = ""),search()))
      load(file.path(path, "data", "@PKGNAME@.rda"), envir = where)
}
