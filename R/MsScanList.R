################################################################################
# TODO: Ordering of plots doesn't follow scan order
#       writeMGF for MsScanList
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
    contains = 'MsList',
    validity = function(object) {
        if(!('mode' %in% colnames(object@mapping)) & length(object) != 0) {
            return('MsScanList must contain the mode of the spectrum')
        }
        return(TRUE)
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
    function(object, extraData) {
        info <- msInfo(object)
        if(any(isEmpty(object))) {
            warning('Empty scans removed')
        }
        scanIndex <- getElementIndex(object@mapping)
        meltData <- data.frame(sample=uNames(object)[scanIndex], mode=scanMode(object)[scanIndex], object@data)
        if(!missing(extraData)) {
            meltData <- cbind(meltData, info[scanIndex, extraData])
        }
        meltData
    }
)
annotateIons <- function(object, precursor=FALSE, parent=FALSE, peaks=FALSE, masses=0, tolerance=20) {
    if(!any(precursor, parent, peaks, masses!=0)) return(NULL)
    info <- msInfo(object)
    mode <- scanMode(object)
    
    res <- list()
    
    for(i in 1:length(object)) {
        scan <- as.data.frame(object[[i]])
        scan$sample <- rownames(info)[i]
        
        if(mode[i] == 'profile') {
            scan <- scan[which(diff(sign(diff(scan$intensity)))==-2)+1, ]
        }
        scan$peak <- FALSE
        scan$type <- NA_character_
        scan$mass <- FALSE
        
        if(peaks) {
            p <- peaks(object[i])[[1]]
            if(length(p) != 0) {
                pInfo <- msInfo(p)
                lx <- outer(scan$mz, pInfo$mzMin, ">=" )
                hx <- outer(scan$mz, pInfo$mzMax, "<=" )
                scan$peak <- apply(lx*hx, 1, sum) > 0
            }
        }
        
        if(parent && info$precursorMZ[i] != 0) {
            parentInd <- which.min(abs(scan$mz - info$precursorMZ[i]))
            if(abs(scan$mz[parentInd] - info$precursorMZ[i]) < info$precursorMZ[i]*tolerance/1e6) {
                scan$type[parentInd] <- 'Parent'
            }
        }
        
        if(precursor) {
            acqNum <- getAcqNum(con(object, i), precursorScanNum=info$acquisitionNum[i], raw=isRaw(object)[i])
            if(length(acqNum) != 0) {
                childHead <- con(object, i)$getHeader(acqNum[[1]], raw=isRaw(object)[i])
                diffs <- abs(outer(scan$mz, childHead$precursorMZ, `-`))
                childInd <- apply(diffs, 2, which.min)
                scan$type[childInd] <- 'Precursor'
            }
        }
        
        if(masses != 0) {
            nLab <- min(nrow(scan), masses)
            labelInds <- order(scan$intensity, decreasing = T)[1:nLab]
            scan$mass[labelInds] <- TRUE
        }
        
        scan <- scan[scan$peak | !is.na(scan$type) | scan$mass, ]
        if(!peaks) scan$peak <- NULL
        if(!parent && !precursor) scan$type <- NULL
        if(masses == 0) scan$mass <- NULL
        
        res[[i]] <- scan
    }
    do.call(rbind, res)
}

#' @describeIn MsScanList Create a plot
#' 
#' @import ggplot2
#' 
setMethod(
    'msPlot', 'MsScanList',
    function(object, simple = FALSE, parent=!simple, precursor=!simple, peaks=!simple, masses=ifelse(simple, 0, 10), context=!simple, tolerance=20, ...) {
        data <- meltMS(object)
        ann <- annotateIons(object, precursor, parent, peaks, masses, tolerance)
        p <- ggplot(data) + theme_bw()
        if(any(scanMode(object) == 'centroid')) {
            p <- p + geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity), data=subset(data, mode=='centroid'), color=ifelse(peaks, 'grey', 'black'))
        }
        if(any(scanMode(object) == 'profile')) {
            p <- p + geom_line(aes(x=mz, y=intensity), data=subset(data, mode=='profile'), color=ifelse(peaks, 'grey', 'black'))
        }
        if(peaks && any(ann$peak)) {
            p <- p + geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity, color=peak), data=ann[ann$peak,, drop=FALSE])
            p <- p + scale_colour_manual('', values='steelblue', labels='Peaks')
        }
        if(parent | precursor) {
            p <- p + geom_point(aes(x=mz, y=intensity, shape=type), data=ann[!is.na(ann$type), ], color='forestgreen', size=I(3))
            p <- p + scale_shape('')
        }
        if(masses != 0) {
            p <- p + geom_text(aes(x=mz, y=intensity, label=format(mz)), data=ann[ann$mass, , drop=FALSE], size=2.5, vjust=-0.6, hjust=0.5)
        }
        p <- p + facet_grid(sample ~ ., scales='free_y')
        p <- p + scale_y_continuous('Intensity') + scale_x_continuous('m/z')
        if(length(object) == 1 & context) {
            suppressMessages({
                oPlot <- ggplotGrob(p + scale_y_continuous(''))
            })
            if(msInfo(object)$msLevel == 1) {
                cChrom <- chroms(con(object, 1, 'MsData'), msLevel=1)
                cData <- data.frame(cChrom[[1]])
                cData$name <- 'TIC'
                cPlot <- ggplot(cData, aes(x=retentionTime, y=TIC)) + theme_bw()
                cPlot <- cPlot + geom_line(colour='grey')
                prec <- cData[which.min(abs(cData$retentionTime-msInfo(object)$retentionTime)), , drop=FALSE]
                cPlot <- cPlot + geom_point(data=prec)
                cPlot <- cPlot + facet_grid(name ~ .) + scale_y_continuous('') + scale_x_continuous('Retention time (sec)')
            } else {
                cScan <- parent(object)
                cData <- meltMS(cScan)
                if(scanMode(cScan) == 'profile') {
                    cData <- cData[which(diff(sign(diff(cData$intensity)))==-2)+1, ]
                }
                cData$name <- 'Precursor scan'
                prec <- cData[which.min(abs(cData$mz-msInfo(object)$precursorMZ)), , drop=FALSE]
                cPlot <- ggplot(cData) + theme_bw()
                cPlot <- cPlot + geom_vline(aes(xintercept=mz), data=prec, linetype='dashed', colour='darkgrey')
                cPlot <- cPlot + geom_segment(aes(x=mz, xend=mz, y=intensity, yend=0), colour='grey')
                cPlot <- cPlot + geom_segment(aes(x=mz, xend=mz, y=intensity, yend=0), data=prec)
                cPlot <- cPlot + facet_grid(name ~ .) + scale_y_continuous('') + scale_x_continuous('')
            }
            cPlot <- ggplotGrob(cPlot)
            p <- rbindGtable(cPlot, oPlot, size='max')
            p <- gtable_add_grob(p, textGrob('Intensity', rot=90), 3, 2, 9)
            p$heights[p$layout$t[grep('panel', p$layout$name)]] <- lapply(c(1,3), grid::unit, 'null')
            p <- p[-c(5,7),]
            plot(p)
            invisible(p)
        } else {
            p
        }
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

