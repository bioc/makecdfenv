## Laurent 2002 tneruaL
## Class definitions for 'Cdf'
setClass("Cdf", representation(cdfName="character",
                                 name="matrix",
                                 name.levels="character",
                                 pbase="matrix",
                                 pbase.levels="character",
                                 tbase="matrix",
                                 tbase.levels="character",
                                 atom="matrix"))
  ## name.levels method
  if( !isGeneric("name.levels") ) {
    setGeneric("name.levels", function(object)
               standardGeneric("name.levels"))
  } else
    cat("name.levels is already generic, could be a problem.\n")

   setMethod("name.levels","Cdf",function(object) object@name.levels)

  if( !isGeneric("name.levels<-") )
      setGeneric("name.levels<-", function(object, value)
                 standardGeneric("name.levels<-"))
  setReplaceMethod("name.levels", "Cdf", function(object, value){
    object@name.levels <- value
    object
    })

  ## pbase method
  if( !isGeneric("pbase") ) {
    setGeneric("pbase", function(object)
               standardGeneric("pbase"))
  } else
    cat("pbase is already generic, could be a problem.\n")

   setMethod("pbase","Cdf",function(object) object@pbase)

  if( !isGeneric("pbase<-") )
      setGeneric("pbase<-", function(object, value)
                 standardGeneric("pbase<-"))
  setReplaceMethod("pbase", "Cdf", function(object, value){
    object@pbase <- value
    object
    })

 ##pbase.levels method
  if( !isGeneric("pbase.levels") ) {
    setGeneric("pbase.levels", function(object)
               standardGeneric("pbase.levels"))
  } else
    cat("pbase.levels is already generic, could be a problem.\n")

   setMethod("pbase.levels","Cdf",function(object) object@pbase.levels)

  if( !isGeneric("pbase.levels<-") )
      setGeneric("pbase.levels<-", function(object, value)
                 standardGeneric("pbase.levels<-"))
  setReplaceMethod("pbase.levels", "Cdf", function(object, value){
    object@pbase.levels <- value
    object
    })

  ## tbase method
  if( !isGeneric("tbase") ) {
    setGeneric("tbase", function(object)
               standardGeneric("tbase"))
  } else
    cat("tbase is already generic, could be a problem.\n")

   setMethod("tbase","Cdf",function(object) object@tbase)

  if( !isGeneric("tbase<-") )
      setGeneric("tbase<-", function(object, value)
                 standardGeneric("tbase<-"))
  setReplaceMethod("tbase", "Cdf", function(object, value){
    object@tbase <- value
    object
    })

  ## tbase.levels method
  if( !isGeneric("tbase.levels") ) {
    setGeneric("tbase.levels", function(object)
               standardGeneric("tbase.levels"))
  } else
    cat("tbase.levels is already generic, could be a problem.\n")

   setMethod("tbase.levels","Cdf",function(object) object@tbase.levels)

  if( !isGeneric("tbase.levels<-") )
      setGeneric("tbase.levels<-", function(object, value)
                 standardGeneric("tbase.levels<-"))
  setReplaceMethod("tbase.levels", "Cdf", function(object, value){
    object@tbase.levels <- value
    object
    })

  ## printing method for 'Cdf'
  setMethod("show", "Cdf",
            function(object) {
              ##print("Cdf object\n")
              cat("Cdf object:\n")
              cat("cdfName=", object@cdfName)
              cat("dim(name)=",dim(object@name)[1],"x",
                  dim(object@name)[2]," features\n",sep="")
              cat("nb. affyid=", length(object@name.levels), "\n", sep="")
            })



