\name{read.cdffile}
\alias{read.cdffile}
\title{Read a CDF file}
\description{Read the data contained in a CDF file}
\usage{read.cdffile(file, compress=FALSE)}
\arguments{
  \item{file}{the name of the CDF file}
  \item{compress}{whether the file is compressed or not}
}
\value{
  Returns a \code{\link{Cdf-class}} object.
}
\warnings{
  A \code{Cdf} object is not a cdf environment, which will be needed for
  the computation of expression values from the probe intensities in a
  \code{AffyBatch}. 
}
\details{This function is intended for use by
  \code{\link{make.cdf.env}}, which in turn is called by
  \code{\link{make.cdf.package}}. User may not have much
  benefit from calling this function directly.
  
  In order to save memory, the name corresponding to each value in the
  CEL is a factor. As in R \code{factor} objects
  cannot be also of type \code{matrix}, the names corresponding to the indices
  were stored in a vector of type \code{character} called
  \code{name.levels}.
  The same thing was done with the pbase and cbase information (more for
  consistency than by conviction that some memory
  could be saved here).
}
\author{Laurent Gautier (laurent@cbs.dtu.dk)}
\seealso{\code{\link{Cdf-class}}}
\keyword{file}
\examples{
fn <- file.path(.path.package("makecdfenv"), "data", "Hu6800.CDF.gz")
mycdf <- read.cdffile(fn, compress=TRUE)
mycdf
}  
