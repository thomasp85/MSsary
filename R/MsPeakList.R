#' @include aaa.R
#' @include generics.R
#' @include MsList.R
#' 
NULL

#' Class to handle detected peaks
#' 
#' This class is a container for one or several peaks from MsData object(s).
#' 
#' @slot connections A list of MsConnections
#' 
#' @slot info A data.frame with information on the contained peaks
#' 
#' @slot data A matrix containing the raw data pertaining to the peaks
#' 
#' @slot mapping A matrix with the mapping of peaks to the rows in the 
#' @@data matrix
#' 
#' @family MSsary-classees
#' 
setClass(
    'MsPeakList',
    contains = 'MsList'
)
