% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_data.R
\name{as.sir_vectorised}
\alias{as.sir_vectorised}
\title{Convert MIC or Disk Diffusion to SIR, vectorised over antimicrobials}
\usage{
as.sir_vectorised(mic, mo, ab, accept_ecoff = FALSE, ...)
}
\arguments{
\item{mic}{vector of MIC values}

\item{mo}{vector of microorganism names}

\item{ab}{vector of antibiotic names}

\item{accept_ecoff}{if TRUE, ECOFFs will be used when no clinical breakpoints are available}

\item{...}{additional arguments that are passed to AMR::as.sir}
}
\value{
S3 sir values
}
\description{
The AMR::as.sir function is not vectorised over antimicrobials. This function
provides vectorisation over antimicrobials. Due to the overhead of running
AMR::as.sir, this function tries to be efficient by only running AMR::as.sir
as little as necessary.
}
\examples{
mic <- c("<0.25", "8", "64", ">64")
mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI")
ab <- c("AMK", "AMK", "AMK", "AMK")
as.sir_vectorised(mic, mo, ab)
# using different microorganisms and antibiotics
mic <- c("<0.25", "8", "64", ">64")
mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_PROTS_MRBL", "B_PROTS_MRBL")
ab <- c("AMK", "AMK", "CIP", "CIP")
as.sir_vectorised(mic, mo, ab)
}
