\name{make.cdf.env}
\alias{make.cdf.env}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{CDF Environment Maker}
\description{
Reads in an Affymetrix chip description file (CDF) and creates an
environment containing a hash table with a location to probe set
membership mapping.}
\usage{
make.cdf.env(filename, cdf.path = getwd(), env.name = "tmpenv",
             verbose = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{filename}{a character containing the filename of the CDF file. We
    assume the fullpath is *not* given.}
  \item{cdf.path}{The path needed to find the CDF file.}
  \item{env.name}{A character containing an environment name to use
    within the function. User should  not need to bother with this argument.}
  \item{verbose}{logical. If \code{TRUE} messages are shown.}
}
\value{An environment containing a hash table. For every probe set name
  we have a matrix with 2 columns. The first column contains the PM
  locations and the second column the MM locations. For PM only chips
  the MM column will have NAs.}
\author{Rafael A. Irizarry}
\seealso{\code{\link{make.cdf.package}}}
\keyword{manip}% at least one, from doc/KEYWORDS
