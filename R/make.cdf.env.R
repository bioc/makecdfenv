make.cdf.env <- function(filename,cdf.path=getwd(),env.name="tmpenv",
                         verbose=TRUE){
  ##read in the cdf file into a CDF object
  if(verbose) cat("Reading CDF file.\n")
  cdf <- read.cdffile(paste(path.expand(cdf.path),filename,sep="/"))

  if(verbose) cat("Creating environment\n")
  ##extract number or rows and number of columns
  nrow <- dim(cdf@name)[1]
  ncol <- dim(cdf@name)[2]
  
  ##where are the pm and mm
  pm.or.mm <- as.vector(pmormm(cdf))
  pmindex <- pm.or.mm
  mmindex <- !pm.or.mm
  rm(pm.or.mm)
  
  ##get ready to go
  ##start the new environment
  ##position in vector will correspont to the position on the array
  ##given ncol and nrow we can figure out the 2D (x,y) position
  ##now we assing 
  assign(env.name,new.env(hash=TRUE))
  genenames <- name.levels(cdf) ##so that we only use name.level method once
  n <- length(genenames) ##number of genes
  probenames <-  as.vector(cdf@name)
  position <- 1:length(probenames)
  numbers <- as.vector(cdf@atom)##number relates to position in genome
  tmp<- split(c(position,pmindex,mmindex,numbers),
              c(probenames,probenames,probenames,probenames))
  rm(pmindex,mmindex,position,numbers,probenames)
  names(tmp) <- genenames[as.numeric(names(tmp))]
  ## this list will be converted to an environment using multiassign
  if(verbose) cat("Wait for about",round(length(tmp)/100),"dots")
  tmp <- lapply(tmp,function(x){
    if(verbose) if(sample(1:100)==1) cat(".")
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
  return(get(env.name))
}

##very rough package maker. 
make.cdf.package<- function(filename,packagename=NULL,
                            cdf.path=getwd(),package.path=getwd(),
                            author="The Bioconductor Project",
                            maintainer="The Bioconductor Project <bioconductor@stat.math.ethz.ch>",
                            version=library(help=makecdfenv)$info[[2]][[2]][2],
                            verbose=TRUE){
  ##if no packagename given change CDF filename to packagename
  if(is.null(packagename)){
    packagename <- gsub("_","",tolower(filename))
    packagename <- gsub("\\.","",packagename)
    packagename <- gsub("-","",packagename)
    packagename <- gsub("\ ","",packagename)
  }

  fullpath<- paste(path.expand(package.path),packagename,sep="/")
  if(file.exists(fullpath))
    stop(paste(fullpath," exists. Please remove first.\n",sep=""))
  
  dir.create(fullpath)
  dir.create(paste(fullpath,"data",sep="/"))
  dir.create(paste(fullpath,"man",sep="/"))
  dir.create(paste(fullpath,"R",sep="/"))

  ##TITLE FILE
  write(packagename,file=paste(fullpath,"TITLE",sep="/"))

  ##DESCRIPTION FILE
  write(paste("Package: ",packagename,"\n",
              "Title: ",packagename,"\n",
              "Version: ",version,"\n",
              "Created: ",date(),"\n",
              "Author: ",author,"\n",
              "Description: A package containing an evironmet represeting the ",filename," file.\n",
              "Maintainer: ",maintainer,"\n",
              "License: LGPL\n",sep=""),
        file=paste(fullpath,"DESCRIPTION",sep="/"))

  ##INDEX file
  write(paste(packagename,"environment containing the location probe set membership mapping",sep="\t"),file=paste(fullpath,"INDEX",sep="/"))
  
  ##data/00Index file
  write(paste(packagename,"environment containing the location probe set membership mapping",sep="\t"),file=paste(fullpath,"data","00Index",sep="/"))
  
  ##zzz.R
  write(paste(".First.lib <- function(libname, pkgname){\n",
              "\tpath = .path.package(pkgname)\n",
              "\twhere <- as.environment(match(paste(\"package:\", pkgname, sep = \"\"),search()))\n",
              paste("\tload(file.path(path, \"data\", \"",
                    paste(packagename,".rda",sep=""),
                    "\") ,envir = where)\n",sep=""),
              "}\n",sep=""),file=paste(fullpath,"R","zzz.R",sep="/"))
  
  ##Rd file
  write(paste("\\name{",packagename,"}\n",
              "\\alias{",packagename,"}\n",
              "\\non_function{}\n",
              "\\title{",packagename,"}\n",
              "\\description{environment describind CDF file}\n",
              "\\keyword{datasets}\n",sep=""),
        file=paste(fullpath,"man",paste(packagename,".Rd",sep=""),sep="/"))
  

  ###the environmet
  tmp <- make.cdf.env(filename,cdf.path=cdf.path)

  assign(packagename,tmp)
  save(list=packagename,file=paste(package.path,"/",packagename,"/data/",packagename,".rda",sep=""))
  TRUE
}








