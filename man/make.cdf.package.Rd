\name{make.cdf.package}
\alias{make.cdf.package}
\title{CDF Environment Package Maker}
\description{
This function reads an Affymetrix chip description file (CDF) and
creates an R package that when loaded has the CDF environment available
for use.}
\usage{
make.cdf.package(filename,
  packagename  = NULL,
  cdf.path     = getwd(),   
  package.path = getwd(),
  compress     = FALSE,
  author       = "The Bioconductor Project",
  maintainer   = "The Bioconductor Project <bioconductor@stat.math.ethz.ch>",
  version      = packageDescription("makecdfenv", field="Version"),
  unlink       = FALSE,
  verbose      = TRUE)
}
\arguments{
  \item{filename}{Character. Filename of the CDF file -
    \bold{without} the path prefix!}
  \item{packagename}{Character. Name wanted for the package.}
  \item{cdf.path}{Character. Path to the CDF file.}
  \item{package.path}{Character. Path where the package will be created.}
  \item{compress}{Logical. If \code{TRUE}, CDF file is compressed.}
  \item{author}{Character. What to put in the \code{author} field of the package.}
  \item{maintainer}{Character. What to put in the \code{maintainer} field of the package.}
  \item{version}{Character. What to put in the version
   field. Should be a of the form \code{x.x.x}.}
  \item{unlink}{Logical. If \code{TRUE}, and a package directory exists
   already in \code{package.path}, that is overwritten.}
  \item{verbose}{Logical. If \code{TRUE} messages are shown.}
}
\details{The function is called for its side effect, creating a package.
By default the package name will be the name of the CDF file made lower
case and with special characters removed (i.e. only
alpha-numeric). 
  
In general one would want to use the name given in by
\code{cleancdfname(abatch@cdfName)} with \code{abatch} an
\code{\link[affy]{AffyBatch}} object obtained, for example, using
\code{\link[affy]{ReadAffy}}. This is the package name that the
\code{affy} package looks for by default.

If the user has a CEL file, called \code{filename}, the recommended
package name for the environment is
\code{cleancdfname(whatcdf(filename))}. This usually coincides with
the default.

Please see the vignette for more details.
}
\value{
If success, the function returns the name of the created package.
}
\seealso{\code{\link{make.cdf.env}}}
\author{Rafael A. Irizarry, Wolfgang Huber}
\examples{
  pkgpath <- tempdir()
  make.cdf.package("Hu6800.CDF.gz",
           cdf.path=system.file("extdata", package="makecdfenv"),
           compress=TRUE,
	   package.path = pkgpath)
  dir(pkgpath)	   
}

\keyword{manip}




