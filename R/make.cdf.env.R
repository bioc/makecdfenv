##------------------------------------------------------------
## (C) Rafael Irizarry, Wolfgang Huber 2003
##
## Modifications by B. M. Bolstad for binary CDF files Feb 2005
## Note that the binary parser ignores the cdffile object
## altogether and creates much of the structure in the c code
##
##------------------------------------------------------------
make.cdf.env <- function(filename,
                         cdf.path = getwd(),
                         compress = FALSE,
                         return.env.only = TRUE,
                         verbose = TRUE) {
  stopifnot(is.logical(verbose), is.logical(return.env.only), is.logical(compress),
            is.character(cdf.path), is.character(filename))
  stopifnot(length(filename)==1, length(cdf.path)==1, length(return.env.only)==1,
            length(verbose)==1, length(compress)==1)


  isCDFXDA <- function(filename){
    return (.Call("CheckCDFXDA",filename,PACKAGE="affyio"))
  }



  if (!isCDFXDA(file.path(path.expand(cdf.path),filename))){
    # a "classic" CDF file (ie text)

    ## read in the cdf file into a CDF object
    if(verbose)
      cat("Reading CDF file.\n")
    cdf <- read.cdffile(file.path(path.expand(cdf.path),filename), compress=compress)
    
    if(verbose)
      cat("Creating CDF environment\n")
  
    ## extract number or rows and number of columns
    sizex <- dim(cdf@name)[1]
    sizey <- dim(cdf@name)[2]
  
    ## where are the pm and mm
    pm.or.mm <- as.vector(pmormm(cdf))
    pmindex <- pm.or.mm
    mmindex <- !pm.or.mm
    rm(pm.or.mm)
    
    ## get ready to go
    ## start the new environment
    ## position in vector will correspond to the position on the array
    ## given ncol and nrow we can figure out the 2D (x,y) position
    ## now we assign 
    
    env = new.env(hash=TRUE, parent=emptyenv())
    
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
    
    if(verbose)
      cat("\n")
    
    multiassign(names(tmp),tmp, env)
    
    syms = list(
      I2XY  = paste("y*", sizex, "+x+1", sep=""),
      XY2I  = paste("r=cbind((i-1)%%",sizex,",(i-1)%/%",sizex,"); colnames(r)=c('x','y'); return(r)", sep=""),
      SIZEX = paste(sizex),
      SIZEY = paste(sizey),
      SIZEI = paste(sizex*sizey))
    
    ## return:
    if (return.env.only) {
      rv = env
    } else {
      rv = list(env=env, syms=syms)
    }
  } else {
    ## Binary CDF file("ReadCDFFile",
    tmp <- .Call("ReadCDFFile",file.path(path.expand(cdf.path),filename),PACKAGE="affyio")
    sizex <- tmp[[1]][1]
    sizey <- tmp[[1]][2]
    tmp <- tmp[[2]][order(names(tmp[[2]]))]

    env = new.env(hash=TRUE, parent=emptyenv())
    multiassign(names(tmp),tmp, env)
    syms = list(
      I2XY  = paste("y*", sizex, "+x+1", sep=""),
      XY2I  = paste("r=cbind((i-1)%%",sizex,",(i-1)%/%",sizex,"); colnames(r)=c('x','y'); return(r)", sep=""),
      SIZEX = paste(sizex),
      SIZEY = paste(sizey),
      SIZEI = paste(sizex*sizey))
    
    ## return:
    if (return.env.only) {
      rv = env
    } else {
      rv = list(env=env, syms=syms)
    }
    
  }



    
  return(rv)
}

##------------------------------------------------------------
## package maker
##------------------------------------------------------------
make.cdf.package<- function(filename,
                            packagename = NULL,
                            cdf.path = getwd(),
                            package.path = getwd(),
                            compress = FALSE,
                            author = "The Bioconductor Project",
                            maintainer = "Biocore Data Team <biocannotation@lists.fhcrc.org>",
                            version = packageDescription("makecdfenv", field="Version"),
                            species = NULL,
                            unlink  = FALSE,
                            verbose = TRUE) {
  
  ## If no packagename given change CDF filename to packagename.
  ## cleancdfname() is defined in package 'affy'.
  if(is.null(packagename))
    packagename <- cleancdfname(sub("\\.cdf$", "", filename, ignore.case=TRUE))

  if(is.null(species))
    stop("A species name must be specified, using the correct format:\n Example - Homo_sapiens\n")
  pkgroot <- sub("cdf$","", packagename)

  cdf <- make.cdf.env(filename, cdf.path=cdf.path, compress=compress, return.env.only=FALSE)
  assign(packagename, cdf$env)

  ## Add small env containing nrow and ncol
  dimenv <- new.env(hash = TRUE, parent = emptyenv())
  multiassign(c("NROW", "NCOL"), as.numeric(c(cdf$syms$SIZEX, cdf$syms$SIZEY)), dimenv)
  dimenvname <- paste(sub("cdf$", "", packagename), "dim", sep = "")
  assign(dimenvname, dimenv)

  home  <- .path.package("makecdfenv")
  symbols = append(cdf$syms,
    list(PKGNAME     = packagename,
         VERSION     = version,
         AUTHOR      = author,
         TODAY       = date(),
         CDFFILENAME = filename,
         MAINTAINER  = maintainer,
         SPECIES     = species,
         PKGROOT     = pkgroot,
         DIMENVNAME  = dimenvname))

  ## see in the subdirectory 'Code' of this package for all the
  ## files that go into the package!
  res = createPackage(packagename,
                destinationDir = path.expand(package.path),
                originDir  = file.path(home, "Code"),
                symbolValues = symbols,
                unlink = unlink, quiet = !verbose)

  ## save an XDR file with the environment
  save(list = c(packagename, dimenvname),
       file = file.path(res$pkgdir, "data", paste(packagename,".rda",sep="")))

  pkgdir <- file.path(package.path, packagename)

  if(verbose)
      cat("\n\n\nREADME PLEASE:\n",
          "A source package has now been produced in\n", pkgdir, ".\n",
          "Before using this package it must be installed via 'R CMD INSTALL'\n",
          "at a terminal prompt (or DOS command shell).\nIf you are using ",
          "Windows, you will need to get set up to install packages.\n",
          "See the 'R Installation and Administration' manual, specifically\n",
          "Section 6 'Add-on Packages' as well as 'Appendix E: The Windows Toolset'\n",
          "for more information.\n\n\n", sep="")

}
