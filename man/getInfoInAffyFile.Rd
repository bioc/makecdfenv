\name{getInfoInFile}
\alias{getInfoInFile}
\title{get information from Affymetrix data files}
\description{
  Get specific information stored in Affymetrix data files (CEL or CDF).
}
\usage{
getInfoInFile(filename, type, unit, property, compress = NULL)
}
\arguments{
  \item{filename}{a file name }
  \item{type}{\code{"CEL"} or \code{"CDF"} are the only known types.}
  \item{unit}{the 'unit' to find the information}
  \item{property}{the 'property' of interest}
  \item{compress}{a boolean}
}
\details{
  The data files seem to have a structure of 'units'. A 'unit' start
  with something between square brackets (ex: "[UNIT1]"). The parameter
  \code{unit} let one specify in which unit the information of interest
  is located. The 'properties' are constituted of an identifier, the
  sign '=' and the value for the property (ex: "numCells=12000").

  If you are planning to use this function, I assume you know the what
  is in CEL and CDF files in details and that you know what are doing.
}
\value{
  a character with the value of the 'property'.
}
\note{
  The code goes through the file until the rigth 'unit' is found. Then
  it looks for the right 'property' further down (eventually going
  through the next units if the property is not found).
}
\seealso{\code{read.celfile}, \code{read.cdffile}, \code{whatcdf}}
\keyword{manip}
