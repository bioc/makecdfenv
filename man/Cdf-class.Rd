\name{Cdf-class}
\docType{class}
\alias{Cdf-class}
\alias{atom}
\alias{atom<-}
\alias{name.levels}
\alias{name.levels<-}
\alias{pbase}
\alias{pbase<-}
\alias{pbase.levels}
\alias{pbase.levels<-}
\alias{tbase}
\alias{tbase<-}
\alias{tbase.levels}
\alias{tbase.levels<-}
\alias{atom,Cdf-method}
\alias{atom<-,Cdf-method}
\alias{name.levels,Cdf-method}
\alias{name.levels<-,Cdf-method}
\alias{pbase,Cdf-method}
\alias{pbase<-,Cdf-method}
\alias{pbase.levels,Cdf-method}
\alias{pbase.levels<-,Cdf-method}
\alias{tbase,Cdf-method}
\alias{tbase<-,Cdf-method}
\alias{tbase.levels,Cdf-method}
\alias{tbase.levels<-,Cdf-method}
\alias{show,Cdf-method}


\title{Class Cdf}
\description{A class describing the content a Chip Description File.}
\section{Creating Objects}{
  \code{  new('Cdf',}\cr
  \code{    cdfName      = ...., # Object of class character}\cr
  \code{    name         = ...., # Object of class matrix}\cr
  \code{    name.levels  = ...., # Object of class character}\cr
  \code{    pbase        = ...., # Object of class matrix}\cr
  \code{    pbase.levels = ...., # Object of class character}\cr
  \code{    tbase        = ...., # Object of class matrix}\cr
  \code{    tbase.levels = ...., # Object of class character}\cr
  \code{    atom         = ...., # Object of class matrix}\cr
  \code{  )}}
\section{Slots}{
  \describe{
    \item{\code{cdfName}:}{the CDF name tag. Used to link with the right
      CEL files}
    \item{\code{name}:}{Object of class "matrix" of 'factors' for the gene names corresponding to the probes.}
    \item{\code{name.levels}:}{Object of class "character" containing
      the levels corresponding to \code{name}.} 
    \item{\code{pbase}:}{Object of class "matrix" of pbase levels.}
    \item{\code{pbase.levels}:}{Object of class "character"  containing
      the levels corresponding to \code{pbase}. }
    \item{\code{tbase}:}{Object of class "matrix" of tbase levels.}
    \item{\code{tbase.levels}:}{Object of class "character" containing
      the levels corresponding to \code{tbase}.  }
    \item{\code{atom}:}{Object of class "matrix" of atom or probe numbers.}
  }
}

\details{
  The class \code{Cdf} was designed to store the content of Affymetrix's
  Chip Definition Files (CDF). This early class is no longer widely used
  in the package affy. Environments (called \code{cdfenvs}) are
  preferred to have an efficient mapping between probe set identifiers
  and indexes (needed to access/subset particular probe intensities).
  Unless one needs to access every information contained in a
  CDF file, the \code{cdfenvs} will be preferred.

  The following comments are only relevant to someone with interest in
  what is in a CDF file.
 The name associated
 to each probe is not unique, it corresponds to the gene name. It is
 very common to have a name repeated 40 times (20 perfect matches and 20
 mismatches). Probes can be grouped by pairs: a perfect
 match (PM) probe has its mismatch (MM) counterpart. The two probes in a
 pair differ by one base (usually located in the middle of the
 sequence). The information relative to this particular base are stored
 in \code{pbase} and \code{tbase}. At a given position (x,y), having
 pbase[x,y] == tbase[x,y] means having a MM while having pbase[x,y] !=
 tbase[x,y] means having a perfect mismatch if and only if one of the
 bases is A while the other is TRUE, or one is G and the other is C. The
 function 
 \code{\link{pmormm}} returns \code{TRUE} for PMs and \code{FALSE} for
 MMs.

 To know more about the file structure of a CDF file, one has to refer
 to the parsing code (mostly in C).
}
\seealso{
\code{\link{read.cdffile}},
\code{\link{make.cdf.env}},
\code{\link{make.cdf.package}}
}
\author{L. Gautier <laurent@cbs.dtu.dk>}

\section{Methods}{
  \describe{
    \item{atom}{(Cdf): An accesor function for the
      \code{ato,} slot. }
    \item{atom<-}{(Cdf):  A replacement function for
      \code{atom} slot. }
    \item{name.levels}{(Cdf):  An accesor function for the
      \code{name.levels} slot. }
    \item{name.levels<-}{(Cdf): A replacement function for
      \code{name.levels} slot. }
    \item{pbase}{(Cdf):  An accesor function for the
      \code{pbase} slot. }

    \item{pbase<-}{(Cdf): A replacement function for
      \code{pbase} slot. }
    \item{pbase.levels}{(Cdf):  An accesor function for the
      \code{pbase.levels} slot. }
    \item{pbase.levels<-}{(Cdf): A replacement function for
      \code{pbase.levels} slot. }
    \item{show}{(Cdf):  renders information about the \code{Cdf}
      object in a concise way on stdout. }
    \item{tbase}{(Cdf):  An accesor function for the
      \code{tbase} slot. }
    \item{tbase<-}{(Cdf): A replacement function for
      \code{tbase} slot. }
    \item{tbase.levels}{(Cdf):  An accesor function for the
      \code{tbase.levels} slot. }
    \item{tbase.levels<-}{(Cdf): A replacement function for
      \code{tbase.levels} slot. }
  }
}
\keyword{methods}
