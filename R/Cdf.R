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
    setGeneric("name.levels", function(object)
               standardGeneric("name.levels"))

   setMethod("name.levels","Cdf",function(object) object@name.levels)

  setGeneric("name.levels<-", function(object, value)
                 standardGeneric("name.levels<-"))
  setReplaceMethod("name.levels", "Cdf", function(object, value){
    object@name.levels <- value
    object
    })

  ## pbase method
    setGeneric("pbase", function(object)
               standardGeneric("pbase"))

   setMethod("pbase","Cdf",function(object) object@pbase)

  setGeneric("pbase<-", function(object, value)
                 standardGeneric("pbase<-"))
  setReplaceMethod("pbase", "Cdf", function(object, value){
    object@pbase <- value
    object
    })

 ##pbase.levels method
    setGeneric("pbase.levels", function(object)
               standardGeneric("pbase.levels"))

   setMethod("pbase.levels","Cdf",function(object) object@pbase.levels)

  setGeneric("pbase.levels<-", function(object, value)
                 standardGeneric("pbase.levels<-"))
  setReplaceMethod("pbase.levels", "Cdf", function(object, value){
    object@pbase.levels <- value
    object
    })

  ## tbase method
    setGeneric("tbase", function(object)
               standardGeneric("tbase"))

    setMethod("tbase","Cdf",function(object) object@tbase)

  setGeneric("tbase<-", function(object, value)
                 standardGeneric("tbase<-"))
  setReplaceMethod("tbase", "Cdf", function(object, value){
    object@tbase <- value
    object
    })

  ## tbase.levels method
    setGeneric("tbase.levels", function(object)
               standardGeneric("tbase.levels"))

   setMethod("tbase.levels","Cdf",function(object) object@tbase.levels)

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



