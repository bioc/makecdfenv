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
  author       = "The Bioconductor Project",
  maintainer   = "The Bioconductor Project <bioconductor@stat.math.ethz.ch>",
  version      = library(help=makecdfenv)$info[[2]][[2]][2],
  force        = FALSE, verbose = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{filename}{A character containing the filename of the CDF file -
    \bold{without} the path prefix!}
  \item{packagename}{A character containing the name wanted for the package.}
  \item{cdf.path}{A character with the path to the CDF file.}
  \item{package.path}{A character with the path where the package will be created.}
  \item{author}{A character denoting what to put in the author field.}
  \item{maintainer}{A character denoting what to put in the maintainer field.}
  \item{version}{A character denoting what to put in the version
   field. Should be a of the form \code{x.x.x}.}
  \item{force}{Logical. If \code{TRUE}, and a package directory exists
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
\example{
  ## Presently, the example is commented out since including
  ## the file HG_U95Av2.CDF would add another 30 MB to the package.
  ## Compression would put this down to 5 MB but is not supported on
  ## all platforms. If you have an example CDF file that is small
  ## enough to be included, please let me know - wh
   
  ## make.cdf.package("HG_U95Av2.CDF")
}

\keyword{manip}




