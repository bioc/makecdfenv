\name{make.cdf.env}
\alias{make.cdf.env}
\title{CDF Environment Maker}
\description{Reads an Affymetrix chip description file (CDF) and creates an
environment used as a hash table for the probe set mapping to location.}
\usage{
make.cdf.env(filename,
             cdf.path = getwd(),
             compress = FALSE,
             return.env.only = TRUE,
             verbose = TRUE)
}
\arguments{
  \item{filename}{Character. Filename of the CDF file -
    \bold{without} the path prefix!}
  \item{cdf.path}{Character. Path to the CDF file.}
  \item{compress}{Logical. If \code{TRUE}, CDF file is compressed.}
  \item{return.env.only}{Logical. If \code{TRUE} (the default),
    then the function returns an environment. Otherwise, a list with
    two elements, the first being the environment, and the second being
    a data structure that contains additional information needed for
    the package builder (see details, and vignette).}
  \item{verbose}{Logical. If \code{TRUE}, messages are shown.}
}

\value{Depending on the argument \code{return.env.only},
  either the environment, or a list with two elements,
  \code{env} and \code{syms}. See details.}

\details{Normally, this function should not be called directly. The
  preferred way to handle CDF information is to use
  \code{\link{make.cdf.package}} to build a package, and to
  install it into R. The CDF information can then either be invoked
  automatically by the package \code{affy}, or can be loaded manually by
  calling, for example, \code{library(hgu133a)}.
  
  Some R installations (typically on Windows) do not offer all the tools
  that are necessary for package building. In such situations, this
  function may be called directly. Please see the vignette for details
  - type:

  \code{openVignette("makecdfenv")}

  \bold{Return values:}  
  \code{env} is an environment, used as a hash table.
  For every probe set name we have a matrix with 2 columns. The first
  column contains the PM locations and the second column the MM
  locations. For PM only chips the MM column will have NAs.
  
  \code{syms} is a list that contains chip-specific (i.e., CDF-file
  specific) information that can be used in the construction of the help
  files for the CDF package.
}
\seealso{\code{\link{make.cdf.package}}}
\author{Rafael A. Irizarry, Wolfgang Huber}
\example{
  env <- make.cdf.env("Hu6800.CDF.gz",
           cdf.path=file.path(.path.package("makecdfenv"), "data"),
           compress=TRUE)
  length(ls(env))
  get("U53347_at", env)	   
}
\keyword{manip}
