% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{scanNoise}
\alias{scanNoise}
\title{Get intensity threshold for a given signal-to-noise ratio}
\usage{
scanNoise(scan, sn, rho)
}
\arguments{
\item{scan}{A matrix with mz and intensity values}

\item{sn}{The required minimum signal-to-noise to be considered a real
signal}

\item{rho}{The modifier to use for the second lowest ion special case}
}
\value{
A numeric with the intensity of the first ion that gets accepted as
a true signal
}
\description{
This function takes a scan and finds the first non-noise ion in it based on
the dynamic noise level algorithm described by Xu and Freitas (2009). It
returns the intensity of that ion.
}
\references{
Xu, H., & Freitas, M. A. (2009). A dynamic noise level algorithm
for spectral screening of peptide MS/MS spectra. BMC Bioinformatics, 11,
436–436. doi:10.1186/1471-2105-11-436
}

