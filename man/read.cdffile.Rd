\name{read.cdffile}
\alias{read.cdffile}
\title{Read a CDF file}
\description{
  Read the data contained in a CDF file.
}
\usage{
read.cdffile(file, compress=TRUE)
}
\arguments{
  \item{file}{the name of the CDF file}
  \item{compress}{whether the file is compressed or not...}
}
\value{
  Return a \code{\link[affy:Cdf-class]{Cdf}} object.
}
\warnings{
  A \code{Cdf} object is not a cdf environment, which will be needed for
  the computation of expression values from the probe intensities in a
  \code{AffyBatch}. \code{Cdf} is for you if you are doing probe level work.
}
\details{
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
\seealso{
  \code{\link[affy:Cdf-class]{Cdf}}
}

\keyword{file}

