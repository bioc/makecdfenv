\name{make.cdf.package}
\alias{make.cdf.package}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{CDF Environment Package Maker}
\description{
This function reads an Affymetrix chip description file (CDF) and
creates an R package that when loaded has the CDF environment available
for use.}
\usage{
make.cdf.package(filename, packagename=NULL,
                 cdf.path = getwd(), package.path = getwd(),
                 author = "The Bioconductor Project",
                 maintainer = "The Bioconductor Project <bioconductor@stat.math.ethz.ch>",
                 version = library(help = makecdfenv)$info[[2]][[2]][2],
                 verbose = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{filename}{a character containing the filename of the CDF file. We
    assume the fullpath is *not* given.}
  \item{packagename}{a character containg the name wanted for the package}
  \item{cdf.path}{The path needed to find the CDF file.}
 \item{package.path}{The path where the package will be created.}
  \item{author}{character denoting what to put in the author field.}
  \item{maintainer}{character denoting what to put in the maintainer field.}
  \item{version}{character denoting what to put in the version
    field. Should be a of the form \code{x.x.x}.}
  \item{verbose}{logical. If \code{TRUE} messages are shown.}
}
\details{
  By default the package name will be the name of the CDF file made
  lower case and with special characters removed (i.e. only
  alpha-numeric).
  For this reason the full path should not be given in
  the \code{filename} argument, instead the user should specify it in
  with the \code{cdf.path}
  argument.
  
  In general one would want to use the name given in by
  \code{cleancdfname(abatch@cdfName)} with \code{abatch} an
  \code{\link[affy]{AffyBatch}} object obtained, for example, using
  \code{\link[affy]{ReadAffy}}. This is the package name that the
  \code{affy} package looks for by default.

  If the user has a CEL file, called \code{filename}, the recommended
  package name for the environment is
  \code{cleancdfname(whatcdf(filename))}. This usually coincides with
  the default.
}
\value{
  If it works it will return \code{TRUE}. A package will be created as well.
}
\author{Rafael A. Irizarry}
\seealso{\code{\link{make.cdf.env}}}

\keyword{manip}% at least one, from doc/KEYWORDS




