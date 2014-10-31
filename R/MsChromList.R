################################################################################
# TODO: ions function for MsChromList
#

#' @include aaa.R
#' @include generics.R
#' @include MsList.R
#' 
NULL

#' Class to handle extracted ion chromatograms
#' 
#' This class is a container for one or several XIC from MsData object(s).
#' 
#' @slot connections A list of MsConnections
#' 
#' @slot info A data.frame with information on the contained XIC
#' 
#' @slot data A matrix containing the raw data pertaining to the XIC
#' 
#' @slot mapping A matrix with the mapping of XIC to the rows in the 
#' @@data matrix
#' 
#' @family MSsary-classees
#' 
setClass(
    'MsChromList',
    contains = 'MsList'
)

### METHODS

#' @describeIn MsChromList Short summary of object
#' 
#' @param object An MsChromList object
#' 
setMethod(
    'show', 'MsChromList',
    function(object) {
        cat('An MsChromList object with', length(object), 'chromatograms\n')
    }
)

#' @describeIn MsChromList Subset an MsList object
#' 
#' @param x An MsChromList object
#' 
setMethod(
    '[', c('MsChromList', 'numeric', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        callNextMethod()
    }
)
#' @describeIn MsChromList Subset an MsList object
#' 
setMethod(
    '[', c('MsChromList', 'logical', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        x[which(i)]
    }
)

#' @describeIn meltMS melt an MsChromList
#' 
#' @importFrom reshape2 melt
#' 
setMethod(
    'meltMS', 'MsChromList',
    function(object) {
        if(any(isEmpty(object))) {
            warning('Empty chromatograms removed')
        }
        chromIndex <- getElementIndex(object@mapping)
        meltData <- data.frame(chrom=msInfo(object)$name[chromIndex], object@data)
        meltData <- melt(meltData, id.vars=1:3, variable.name='type', value.name='intensity')
        meltData
    }
)

#' @describeIn MsChromList Create a plot
#' 
#' @param type Which type of plots should be plotted
#' 
#' @param collapse Should 'BPC' and 'TIC' be overlayed or on separate plots, 
#' TRUE means overlayed.
#' 
#' @import ggplot2
#' 
setMethod(
    'msPlot', 'MsChromList',
    function(object, type=c('TIC', 'BPC'), collapse=FALSE, ...) {
        data <- meltMS(object)
        p <- ggplot(data[data$type %in% type,]) + theme_bw()
        if(collapse) {
            p <- p + geom_line(aes(x=retentionTime, y=intensity, colour=chrom, linetype=type))
            p <- p + scale_linetype('Type')
        } else {
            p <- p + geom_line(aes(x=retentionTime, y=intensity, colour=chrom))
            p <- p + facet_grid(type ~ ., scales='free_y')    
        }
        p <- p + scale_y_continuous('Intensity') + scale_x_continuous('Retention time (sec)')
        p <- p + scale_colour_brewer('Chromatogram', palette = 'Set1')
        p <- p + theme(legend.position='bottom', legend.direction='vertical')
        p
    }
)

#' Extract the scans that make up a chromatogram
#' 
setMethod(
    'scans', 'MsChromList',
    function(object, ...) {
        chroms <- msData(object)
        lapply(1:length(chroms), function(i) {
            con <- con(object, i)
            acqNum <- chroms[[i]][, 'acquisitionNum']
            scans <- con$getScans(acqNum)
            info <- dbGetQuery(con$sary(), paste0('SELECT * FROM header WHERE acquisitionNum IN (', paste(acqNum, collapse=', '), ')'))
            mapping <- getListMapping(scans, 1)
            mapping <- cbind(mapping, matrix(isCentroided(scans), dimnames = list(NULL, 'mode')))
            scans <- do.call(rbind, scans)
            colnames(scans) <- c('mz', 'intensity')
            new('MsScanList', connections=list(con), info=info, data=scans, mapping=mapping)
        })
    }
)

#' @describeIn MsChromList Get information about the chromatograms
#' 
#' @return \code{msInfo}: A data.frame with one row containing information
#' about the chromatograms
#' 
setMethod(
    'msInfo', 'MsChromList',
    function(object) {
        callNextMethod()
    }
)

#' @describeIn MsChromList Get the chromatogram data
#' 
#' @return \code{msData}: A list of matrices with acquisitionNum, retentionTime,
#' totIonCurrent and basePeakIntensity for each scan in the chromatogram
#' 
setMethod(
    'msData', 'MsChromList',
    function(object) {
        callNextMethod()
    }
)
