read.cdffile <- function(file, compress=getOption("BioC")$affy$compress.cdf) {

  ff <- new("Cdf")

  ## ---------
  ## the extra operation on the string are done to match what is done
  ## is 'whatcdf'
  tmp <- getInfoInAffyFile(file, "CDF", unit="Chip", property="Name", compress=compress)

  
  tmp <- substr(tmp, 1, nchar(tmp)-2)
  ##we will use cleancdfname later
  ##tmp <- gsub("_","",tmp)
  ##tmp <- tolower(tmp)
  ff@cdfName <- tmp
  ##----------
  tmp <- .Call("readCDFfile", as.character(file),
               as.integer(3), as.integer(compress))
  tmp[tmp == ""] <- NA
  mydim <- dim(tmp)

  tmp <- factor(tmp)
  
  ff@name <- array(as.integer(tmp), mydim)

  ff@name.levels <- levels(tmp)
  rm(tmp)
  gc()
  tmp <- .Call("readCDFfile", as.character(file),
               as.integer(7), as.integer(compress))
  tmp[tmp == ""] <- NA
  mydim <- dim(tmp)
  
  tmp <- factor(tmp)
  
  ff@pbase <- array(as.integer(tmp), mydim)

  ff@pbase.levels <- levels(tmp)
  rm(tmp)
  gc()
  tmp <- .Call("readCDFfile", as.character(file),
               as.integer(8), as.integer(compress))
  tmp[tmp == ""] <- NA
  mydim <- dim(tmp)
  tmp <- factor(tmp)

  ff@tbase <- array(as.integer(tmp), mydim)

  ff@tbase.levels <- levels(tmp)
  rm(tmp)
  gc()
  
  tmp <- .Call("readCDFfile", as.character(file),
               as.integer(9), as.integer(compress))
  tmp[tmp == ""] <- NA
  mydim <- dim(tmp)

  ff@atom <- array(as.integer(tmp), mydim)

  gc() 

  return(ff)
}
