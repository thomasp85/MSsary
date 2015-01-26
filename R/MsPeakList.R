################################################################################
# TODO: 

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
        meltData <- data.frame(peak=uNames(object)[scanIndex], object@data, stringsAsFactors=FALSE)
        meltData
    }
)

annotatePeak <- function(object, precursor) {
    if(!precursor) return(NULL)
    res <- list()
    info <- msInfo(object)
    for(i in 1:length(object)) {
        data <- as.data.frame(object[[i]])
        data$peak <- rownames(info)[i]
        
        data$precursor <- FALSE
        
        acqNum <- getAcqNum(con(object, i), seqNum=BETWEEN(info$scanStart[i], info$scanEnd[i]))
        acqNum <- getAcqNum(con(object, i), precursorScanNum=IN(acqNum[[1]]), precursorMZ=BETWEEN(info$mzMin[i], info$mzMax[i]))
        if(length(acqNum[[1]]) != 0) {
            if(precursor) {
                prec <- con(object, i)$getHeader(acqNum[[1]])
                prec$precursorRT <- con(object, i)$getHeader(prec$precursorScanNum)$retentionTime
                ionMatch <- outer(data$intensity, prec$precursorIntensity, '==') * outer(data$retentionTime, prec$precursorRT, '==')
                data$precursor <- apply(ionMatch, 1, any)
            }
        }
        
        data <- data[data$precursor, , drop=FALSE]
        if(!precursor) data$precursor <- NULL
        
        res[[i]] <- data
    }
    do.call(rbind, res)
}

#' @describeIn MsPeakList Create a plot
#' 
#' @import ggplot2
#' 
setMethod(
    'msPlot', 'MsPeakList',
    function(object, type='2d', collapse=length(object) > 12, simple=FALSE, precursor=!simple, mzscatter=!simple, context=!simple, ...) {
        if(type == '2d') {
            data <- meltMS(object)
            ann <- annotatePeak(object, precursor)
            p <- ggplot(data=data) + theme_bw()
            if(collapse) {
                p <- p + geom_line(aes(x=retentionTime, y=intensity, group=peak, colour=peak))
            } else {
                p <- p + geom_line(aes(x=retentionTime, y=intensity, group=peak))
            }
            p <- p + scale_x_continuous('Retention time (sec)') + scale_y_continuous('Intensity')
            if(precursor & sum(ann$precursor) != 0) {
                p <- p + geom_point(aes(x=retentionTime, y=intensity, colour=precursor), data=ann[ann$precursor,], size=I(3))
                p <- p + scale_colour_manual('', values='forestgreen', labels='Precursor')
            }
            if(length(object) == 1 && any(mzscatter, context)) {
                if(mzscatter) {
                    info <- msInfo(object)
                    scanRange <- c(info$scanStart, info$scanEnd)
                    scanRange <- diff(scanRange)*0.1*c(-1, 1) + scanRange
                    mzRange <- c(info$mzMin, info$mzMax)
                    mzRange <- diff(mzRange)*0.1*c(-1, 1) + mzRange
                    ionData <- ions(con(object, 1, 'MsData'), msLevel=info$msLevel, seqNum=BETWEEN(scanRange[1], scanRange[2]), mz=BETWEEN(mzRange[1], mzRange[2]))
                    ionData <- data.frame(ionData[[1]])
                    ionData$peak <- 'Ion scatter'
                    ionMatch <- outer(ionData$intensity, data$intensity, '==') * outer(ionData$retentionTime, data$retentionTime, '==')
                    ionData$selected <- apply(ionMatch, 1, any)
                    p <- p + geom_point(aes(x=retentionTime, y=mz), data=ionData, colour='grey')
                    p <- p + geom_point(aes(x=retentionTime, y=mz), data=ionData[ionData$selected, ])
                    p <- p + facet_grid(peak~., scales='free_y')
                    p <- ggplotGrob(p)
                    p$heights[p$layout$t[grep('panel', p$layout$name)]] <- lapply(c(3,1), grid::unit, 'null')
                    p$layout$b[p$layout$name=='ylab'] <- 3
                    p <- gtable_add_grob(p, textGrob('m/z', rot=90), 5, 2)
                }
                if(context) {
                    cChrom <- chroms(con(object, 1, 'MsData'), msLevel=1)
                    cData <- data.frame(cChrom[[1]])
                    cData$name <- 'BPC'
                    cPlot <- ggplot(cData, aes(x=retentionTime, y=BPC)) + theme_bw()
                    cPlot <- cPlot + geom_line(colour='grey')
                    cPlot <- cPlot + geom_line(aes(x=retentionTime, y=intensity), data=data)
                    cPlot <- cPlot + facet_grid(name ~ .) + scale_y_continuous('') + scale_x_continuous('') + theme(axis.title.x=element_blank())
                    cPlot <- ggplotGrob(cPlot)
                    if(!inherits(p, 'gtable')) {
                        p <- p + facet_grid(peak~., scales='free_y')
                        p <- ggplotGrob(p)
                    }
                    p <- rbindGtable(cPlot, p, size='max')
                    p$heights[p$layout$t[grep('panel', p$layout$name)[1:2]]] <- lapply(c(1,3), grid::unit, 'null')
                    p$layout$t[p$layout$name=='ylab'] <- 3
                    p <- p[-c(5,7),]
                }
                plot(p)
                invisible(p)
            } else {
                if(!collapse) {
                    p <- p + facet_wrap(~peak, scales='free')
                }
                p
            }
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
            data <- lapply(1:length(data), function(i) {rbind(data[[i]], c(NA,NA,NA))})
            data <- do.call(rbind, data)
            lines3d(x=data[,1], y=data[,3], z=data[,2])
            decorate3d(xlab='Retention time (sec)', ylab='m/z', zlab='Intensity', aspect=TRUE, box=FALSE, cex=3)
        }
    }
)

#' Get scans from peak
#' 
setMethod(
    'scans', 'MsPeakList',
    function(object) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            args$object <- con(object, i, 'MsData')
            args$seqNum <- BETWEEN(info$scanStart[i], info$scanEnd[i])
            args$msLevel <- info$msLevel[i]
            res[[i]] <- do.call(scans, args)
        }
        res
    }
)

#' Get chroms with optional expansion
#' 
setMethod(
    'chroms', 'MsPeakList',
    function(object, mzExpand=0, rtExpand=0) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            args$object <- con(object, i, 'MsData')
            args$retentionTime <- BETWEEN(min(object[[i]][, 'retentionTime'])-rtExpand, max(object[[i]][, 'retentionTime'])+rtExpand)
            args$mz <- BETWEEN(info$mzMin[i]-mzExpand, info$mzMax[i]+mzExpand)
            args$msLevel <- info$msLevel[i]
            res[[i]] <- do.call(chroms, args)
        }
        do.call(c, res)
    }
)
#' Get ions with optional expansion
#' 
setMethod(
    'ions', 'MsPeakList',
    function(object, mzExpand=0, rtExpand=0, ...) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            args$object <- con(object, i, 'MsData')
            args$retentionTime <- BETWEEN(min(object[[i]][, 'retentionTime'])-rtExpand, max(object[[i]][, 'retentionTime'])+rtExpand)
            args$mz <- BETWEEN(info$mzMin[i]-mzExpand, info$mzMax[i]+mzExpand)
            args$msLevel <- info$msLevel[i]
            res[[i]] <- do.call(ions, args)
        }
        do.call(c, res)
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

