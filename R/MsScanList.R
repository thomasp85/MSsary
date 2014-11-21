################################################################################
# TODO: Ordering of plots doesn't follow scan order
#

#' @include aaa.R
#' @include generics.R
#' @include MsList.R
#' 
NULL

#' Class to handle specific scans
#' 
#' This class is a container for one or several scans from MsData object(s).
#' 
#' @slot connections A list of MsConnections
#' 
#' @slot info A data.frame with information on the contained scans
#' 
#' @slot data A matrix containing the raw data pertaining to the scans
#' 
#' @slot mapping A matrix with the mapping of scans to the rows in the 
#' @@data matrix
#' 
#' @family MSsary-classees
#' 
setClass(
    'MsScanList',
    contains = 'MsList'
)
setMethod(
    'initialize', 'MsScanList',
    function(.Object, connections, info, data, mapping) {
        .Object@connections <- connections
        .Object@info <- info
        .Object@data <- data
        .Object@mapping <- mapping
        if(!'mode' %in% colnames(mapping) && nrow(data) != 0) {
            mode <- isCentroided(msData(.Object))
            .Object@mapping <- cbind(.Object@mapping, matrix(mode, dimnames=list(NULL, 'mode')))
        }
        .Object
    }
)

### METHODS

#' @describeIn MsScanList Short summary of object
#' 
#' @param object An MsScanList object
#' 
setMethod(
    'show', 'MsScanList',
    function(object) {
        cat('An MsScanList object with', length(object), 'scans\n')
    }
)
#' @describeIn MsScanList Get the names of the scans
#' 
#' @param x An MsScanList object
#' 
setMethod(
    'names', 'MsScanList',
    function(x) {
        n <- callNextMethod()
        paste0(n, x@info$acquisitionNum)
    }
)

#' @describeIn MsScanList Subset an MsList object
#' 
setMethod(
    '[', c('MsScanList', 'numeric', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        callNextMethod()
    }
)
#' @describeIn MsScanList Subset an MsList object
#' 
setMethod(
    '[', c('MsScanList', 'logical', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        x[which(i)]
    }
)

#' @describeIn MsScanList Get the mode of the scan
#' 
setMethod(
    'scanMode', 'MsScanList',
    function(object) {
        ifelse(as.logical(object@mapping[,'mode']), 'centroid', 'profile')
    }
)

#' @describeIn meltMS melt an MsScanList
#' 
setMethod(
    'meltMS', 'MsScanList',
    function(object) {
        if(any(isEmpty(object))) {
            warning('Empty scans removed')
        }
        scanIndex <- getElementIndex(object@mapping)
        meltData <- data.frame(sample=msInfo(object)$acquisitionNum[scanIndex], mode=scanMode(object)[scanIndex], object@data)
        meltData
    }
)

#' @describeIn MsScanList Create a plot
#' 
#' @import ggplot2
#' 
setMethod(
    'msPlot', 'MsScanList',
    function(object, simple = FALSE,  ...) {
        data <- meltMS(object)
        p <- ggplot(data) + theme_bw()
        if(any(scanMode(object) == 'centroid')) {
            p <- p + geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity), data=subset(data, mode=='centroid'))
        }
        if(any(scanMode(object) == 'profile')) {
            p <- p + geom_line(aes(x=mz, y=intensity), data=subset(data, mode=='profile'))
        }
        if(length(object) == 1 && !simple) {
            scanNum <- msInfo(object)$acquisitionNum
            if(msInfo(object)$msLevel == 1) {
                children <- dbGetQuery(con(object, 1)$sary(), paste0('SELECT precursorMZ FROM header WHERE precursorScanNum == ', scanNum))
                if(nrow(children) != 0) {
                    fragmentScans <- annotateChildren(data, children, scanMode(object))
                    p <- p + geom_point(aes(x=mz, y=intensity, colour=I('black')), data=fragmentScans) + scale_colour_discrete('', breaks='black', labels='Precursor ions')
                }
            }
        }
        p <- p + facet_grid(sample ~ ., scales='free_y')
        p <- p + scale_y_continuous('Intensity') + scale_x_continuous('m/z')
        p
    }
)

#' Get chroms with optional expansion
#' 
setMethod(
    'chroms', 'MsScanList',
    function(object, mz, expand=0) {
        info <- msInfo(object)
        raw <- isRaw(object)
        res <- list()
        for(i in 1:length(object)) {
            args <- list(raw=raw[i])
            args$object <- con(object, i, 'MsData')
            args$retentionTime <- BETWEEN(info$retentionTime[i]-expand, info$retentionTime[i]+expand)
            if(!missing(mz)) {
                if(!rangeFilter(mz)) stop('mz filter must define an interval')
                args$mz <- mz[i]
            }
            args$msLevel <- info$msLevel[i]
            res[[i]] <- do.call(chroms, args)
        }
        do.call(c, res)
    }
)
#' Get ions with optional expansion
#' 
setMethod(
    'ions', 'MsScanList',
    function(object, mz, expand=0, ...) {
        info <- msInfo(object)
        raw <- isRaw(object)
        res <- list()
        for(i in 1:length(object)) {
            args <- list(raw=raw[i])
            args$object <- con(object, i, 'MsData')
            args$retentionTime <- BETWEEN(info$retentionTime[i]-expand, info$retentionTime[i]+expand)
            if(!missing(mz)) {
                if(!rangeFilter(mz)) stop('mz filter must define an interval')
                args$mz <- mz[i]
            }
            args$msLevel <- info$msLevel[i]
            res[[i]] <- do.call(ions, args)
        }
        do.call(c, res)
    }
)

#' Get peaks intersecting with scan
#' 
setMethod(
    'peaks', 'MsScanList',
    function(object, parent=FALSE, ...) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            if(parent) {
                pSeqNum <- dbGetQuery(con(object, i, 'sary'), paste0('SELECT seqNum FROM header WHERE acquisitionNum == ', info$precursorScanNum[i]))$seqNum
                if(length(pSeqNum) == 1) {
                    res[[i]] <- peaks(con(object, i, 'MsData'), seqNum=OVERLAPS(pSeqNum), mz=OVERLAPS(info$precursorMZ[i]), msLevel=info$msLevel-1)
                }
            } else {
                res[[i]] <- peaks(con(object, i, 'MsData'), seqNum=OVERLAPS(info$seqNum[i]), msLevel=info$msLevel[i])                
            }
        }
        res
    }
)

#' @describeIn MsScanList Get the next scan
#' 
#' @return \code{nextScan}: An MsScanList object with the next scan
#' 
setMethod(
    'nextScan', 'MsScanList',
    function(object, sameLevel=TRUE) {
        scInfo <- msInfo(object)
        queries <- mapply(function(x, l) {
            query <- paste0('SELECT * FROM header WHERE seqNum == (SELECT min(seqNum) FROM header WHERE seqNum > ', x)
            if(sameLevel) {
                query <- paste(query, 'AND msLevel ==', l)
            }
            paste0(query, ')')
        }, scInfo$seqNum, scInfo$msLevel, SIMPLIFY = TRUE)
        getScanList(object, queries)
    }
)

#' @describeIn MsScanList Get previous scan
#' 
#' @return \code{previousScan}: An MsScanList object with the previous scan
#' 
setMethod(
    'previousScan', 'MsScanList',
    function(object, sameLevel=TRUE) {
        scInfo <- msInfo(object)
        queries <- mapply(function(x, l) {
            query <- paste0('SELECT * FROM header WHERE seqNum == (SELECT max(seqNum) FROM header WHERE seqNum < ', x)
            if(sameLevel) {
                query <- paste(query, 'AND msLevel ==', l)
            }
            paste0(query, ')')
        }, scInfo$seqNum, scInfo$msLevel, SIMPLIFY = TRUE)
        getScanList(object, queries)
    }
)

#' @describeIn MsScanList Get the scan containing the parent ion for this scan
#' 
#' @return \code{parent}: An MsScanList object with the parent scan
#' 
setMethod(
    'parent', 'MsScanList',
    function(object) {
        scInfo <- msInfo(object)
        queries <- mapply(function(x) {
            query <- paste0('SELECT * FROM header WHERE acquisitionNum == ', x)
        }, scInfo$precursorScanNum, SIMPLIFY = TRUE)
        getScanList(object, queries)
    }
)

#' @describeIn MsScanList Get fragmentation scans from this scan
#' 
#' @return \code{children}: A list of MsScanList objects with the children scans
#' 
setMethod(
    'children', 'MsScanList',
    function(object) {
        lapply(1:length(object), function(i) {
            singleScan <- object[i]
            query <- paste0(
                'SELECT * FROM header WHERE precursorScanNum == ', 
                msInfo(singleScan)$acquisitionNum
            )
            dropEmpty(getScanList(singleScan, query))
        })
    }
)

#' @describeIn MsScanList Get all scans with the same parent
#' 
#' @return \code{siblings}: A list of MsScanList objects with the sibling scans
#' 
setMethod(
    'siblings', 'MsScanList',
    function(object) {
        lapply(1:length(object), function(i) {
            singleScan <- object[i]
            query <- paste0(
                'SELECT * FROM header WHERE precursorScanNum == (SELECT acquisitionNum FROM header WHERE acquisitionNum == ', 
                msInfo(singleScan)$precursorScanNum,
                ')'
            )
            dropEmpty(getScanList(singleScan, query))
        })
    }
)

#' @describeIn MsScanList Get information about the scan
#' 
#' @return \code{msInfo}: A data.frame with one row containing information
#' about the scan
#' 
setMethod(
    'msInfo', 'MsScanList',
    function(object) {
        callNextMethod()
    }
)

#' @describeIn MsScanList Get the scan data
#' 
#' @return \code{msData}: A list of matrices with m/z values in the first 
#' column and intensity in the second
#' 
setMethod(
    'msData', 'MsScanList',
    function(object) {
        callNextMethod()
    }
)

