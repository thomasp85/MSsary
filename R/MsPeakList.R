################################################################################
# TODO: scans, ions and chroms method

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

### METHODS

#' @describeIn MsPeakList Short summary of object
#' 
#' @param object An MsPeakList object
#' 
setMethod(
    'show', 'MsPeakList',
    function(object) {
        cat('An MsPeakList object with', length(object), 'peaks\n')
    }
)
#' @describeIn MsPeakList Get the names of the scans
#' 
#' @param x An MsPeakList object
#' 
setMethod(
    'names', 'MsPeakList',
    function(x) {
        n <- callNextMethod()
        paste0(n, x@info$peakID)
    }
)

#' @describeIn MsPeakList Subset an MsList object
#' 
setMethod(
    '[', c('MsPeakList', 'numeric', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        callNextMethod()
    }
)
#' @describeIn MsPeakList Subset an MsList object
#' 
setMethod(
    '[', c('MsPeakList', 'logical', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        x[which(i)]
    }
)

#' @describeIn meltMS melt an MsPeakList
#' 
setMethod(
    'meltMS', 'MsPeakList',
    function(object) {
        if(any(isEmpty(object))) {
            warning('Empty scans removed')
        }
        scanIndex <- getElementIndex(object@mapping)
        meltData <- data.frame(peak=msInfo(object)$peakID[scanIndex], object@data)
        meltData
    }
)

#' @describeIn MsPeakList Create a plot
#' 
#' @import ggplot2
#' 
setMethod(
    'msPlot', 'MsPeakList',
    function(object, type='2d', collapse=length(object) > 12,  ...) {
        if(type == '2d') {
            data <- meltMS(object)
            p <- ggplot(data=data) + theme_bw()
            p <- p + geom_line(aes(x=retentionTime, y=intensity, group=peak))
            if(!collapse) {
                p <- p + facet_wrap(~peak, scales='free')
            }
            p <- p + scale_x_continuous('Retention time (sec)') + scale_y_continuous('Intensity')
            p
        } else if(type=='3d') {
            if(!require(rgl)) {
                stop('rgl package needed')
            }
            if(!rgl.useNULL()) open3d()
            clear3d()
            skip <- par3d(skipRedraw = TRUE)
            on.exit(par3d(skip))
            data <- msData(object)
            info <- msInfo(object)
            data <- lapply(1:length(data), function(i) {rbind(cbind(data[[i]], rep(info$mzMean[i], nrow(data[[i]]))), c(NA,NA,NA))})
            data <- do.call(rbind, data)
            lines3d(x=data[,1], y=data[,3], z=data[,2])
            decorate3d(xlab='Retention time (sec)', ylab='m/z', zlab='Intensity', aspect=TRUE, box=FALSE, cex=3)
        }
    }
)


#' @describeIn MsPeakList Get information about the scan
#' 
#' @return \code{msInfo}: A data.frame with one row containing information
#' about the scan
#' 
setMethod(
    'msInfo', 'MsPeakList',
    function(object) {
        callNextMethod()
    }
)

#' @describeIn MsPeakList Get the scan data
#' 
#' @return \code{msData}: A list of matrices with m/z values in the first 
#' column and intensity in the second
#' 
setMethod(
    'msData', 'MsPeakList',
    function(object) {
        callNextMethod()
    }
)

