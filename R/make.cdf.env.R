
make.cdf.env <- function(filename,cdf.path=getwd(),env.name="tmpenv",
                         verbose=TRUE){
  ##read in the cdf file into a CDF object
  if(verbose) cat("Reading CDF file.\n")
  cdf <- read.cdffile(paste(path.expand(cdf.path),filename,sep="/"))

  if(verbose) cat("Creating environment\n")
  ##extract number or rows and number of columns
  sizex <- dim(cdf@name)[1]
  sizey <- dim(cdf@name)[2]
  
  ##where are the pm and mm
  pm.or.mm <- as.vector(pmormm(cdf))
  pmindex <- pm.or.mm
  mmindex <- !pm.or.mm
  rm(pm.or.mm)
  
  ## get ready to go
  ## start the new environment
  ## position in vector will correspond to the position on the array
  ## given ncol and nrow we can figure out the 2D (x,y) position
  ## now we assign 
  assign(env.name,new.env(hash=TRUE))
  genenames <- name.levels(cdf) ##so that we only use name.level method once
  n <- length(genenames) ##number of genes
  probenames <-  as.vector(cdf@name)
  position <- 1:length(probenames)
  numbers <- as.vector(cdf@atom)##number relates to position in genome
  tmp<- split(c(position,  pmindex,   mmindex,   numbers),
              c(probenames,probenames,probenames,probenames))
  rm(pmindex,mmindex,position,numbers,probenames)
  names(tmp) <- genenames[as.numeric(names(tmp))]
  ## this list will be converted to an environment using multiassign
  if(verbose)
    cat("Wait for about",round(length(tmp)/100),"dots")
  tmp <- lapply(tmp,function(x){
    if(verbose)
      if(runif(1)<0.01) cat(".")
    n <- length(x)/4
    position <- x[1:n] 
    pmindex <- which(x[(n+1):(2*n)]==1)
    mmindex <- which(x[(2*n+1):(3*n)]==1)
    numbers <- x[(3*n+1):(4*n)]

    ##get pm and mm and put in order according to numbers
    pm <- position[pmindex][order(numbers[pmindex])]
    mm <- position[mmindex][order(numbers[mmindex])]
    if(length(mm)==0)
      return(cbind(pm=pm,mm=rep(NA,length(pm))))
    else
      return(cbind(pm=pm,mm=mm))
  })
  if(verbose) cat("\n")
  multiassign(names(tmp),tmp,get(env.name))

  syms = list(
    XY2I  = paste("xy2i = function(x,y) {y*", sizex, "+x+1}", sep=""),
    I2XY  = paste("i2xy = function(i) {cbind(i%%", sizex, "-1, i%/%", sizex, ")}", sep=""),
    SIZEX = paste(sizex),
    SIZEY = paste(sizey),
    SIZEI = paste(sizex*sizey))
  
  return(list(env=get(env.name), syms=syms))
}

## package maker
make.cdf.package<- function(filename,
                            packagename = NULL,
                            cdf.path = getwd(),
                            package.path = getwd(),
                            author = "The Bioconductor Project",
                            maintainer = "The Bioconductor Project <bioconductor@stat.math.ethz.ch>",
                            version = library(help=makecdfenv)$info[[2]][[2]][2],
                            force   = FALSE,
                            verbose = TRUE) {
  
  ## If no packagename given change CDF filename to packagename.
  ## cleancdfname() is defined in package 'affy'.
  if(is.null(packagename))
    packagename <- cleancdfname(sub("\.cdf$", "", filename, ignore.case=TRUE))

  cdf <- make.cdf.env(filename,cdf.path=cdf.path)
  assign(packagename, cdf$env)

  home  <- .path.package("makecdfenv")
  symbols = append(cdf$syms,
    list(PKGNAME     = packagename,
         VERSION     = version,
         AUTHOR      = author,
         TODAY       = date(),
         CDFFILENAME = filename,
         MAINTAINER  = maintainer))

  ## see in the subdirectory 'Code' of this package for all the
  ## files that go into the package!
  res = createPackage(packagename,
                destinationDir = path.expand(package.path),
                originDir  = file.path(home, "Code"),
                symbolValues = symbols,
                force = force, quiet = !verbose)

  for (fn in c("TITLE", "INDEX", "data/00Index"))
    copySubstitute(file.path(file.path(home, "Code", fn)),
                   file.path(res$pkgdir, fn),
                   symbols)
  
  ## save an XDR file with the environment
  save(list = packagename, file = file.path(res$pkgdir, "data", paste(packagename,".rda",sep="")))
  
  return(TRUE)
}








