getInfoInFile <- function(filename, type, unit, property, compress=NULL) {
  if (is.null(compress) && (type == "CDF"))
    compress <- getOption("BioC")$affy$compress.cdf
  if (is.null(compress) && (type == "CEL"))
    compress <- getOption("BioC")$affy$compress.cel
  .Call("getInfo", as.character(filename), as.character(type),
  as.character(unit), as.character(property), as.integer(compress),
  PACKAGE="makecdfenv")
}
