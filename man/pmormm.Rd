\name{pmormm}
\alias{pmormm}
\title{PM or MM}
\description{
  Determing  if intensities on an array from the cdffile are PM, MM, or neither.
}
\usage{
  pmormm(cdf)
}
\arguments{
  \item{cdf}{A \code{\link[affy:Cdf-class]{Cdf}} object}
}
\details{
  The intensities on a chip fall in three categories: perfect match
  (PM), mismatch (MM) or unknown (No information about them is contained
  in the CDF file, yet they were found informative (details to come)).
  PM are coded as TRUE, MM as FALSE and unknown as NA.
}
\value{A matrix (see section 'details').}
\author{L. Gautier <laurent@cbs.dtu.dk>}
\keyword{manip}
