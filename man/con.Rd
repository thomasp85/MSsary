% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MsData.R, R/MsList.R, R/generics.R
\docType{methods}
\name{con,MsData-method}
\alias{con}
\alias{con,MsData-method}
\alias{con,MsList-method}
\title{Get the MsConnections object}
\usage{
\S4method{con}{MsData}(object)

\S4method{con}{MsList}(object, i, type = "MsConnection")

con(object, ...)
}
\arguments{
\item{object}{The object containing the connection}
}
\value{
An MsConnections object
}
\description{
This function returns the MsConnections object. Only for internal use.
}
\section{Methods (by class)}{
\itemize{
\item \code{MsData}: Get connection from MsData object

\item \code{MsList}: Get connections from an MsList
}}
\keyword{internal}

