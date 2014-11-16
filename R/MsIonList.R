################################################################################
# TODO: 
#

#' @include aaa.R
#' @include generics.R
#' @include MsList.R
#' 
NULL

#' Class to handle raw ion reads
#' 
#' This class is a container for ion data from a slice of an MsData object.
#' Unlike other MsList subclasses this class can only contain one slice.
#' 
#' @slot connections A list of length 1 with the MsConnections object
#' 
#' @slot info A data.frame with information on the slice
#' 
#' @slot data A matrix containing the raw data pertaining to the ions in the 
#' slice
#' 
#' @slot mapping A matrix with the mapping
#' 
#' @family MSsary-classees
#' 
setClass(
    'MsIonList',
    contains = 'MsList'
)

### METHODS

#' @describeIn MsIonList Short summary of object
#' 
#' @param object An MsIonList object
#' 
setMethod(
    'show', 'MsIonList',
    function(object) {
        cat('An MsIonList object with', length(object), 'ion sets\n')
    }
)
#' @describeIn MsIonList Get the names of the scans
#' 
#' @param x An MsIonList object
#' 
setMethod(
    'names', 'MsIonList',
    function(x) {
        n <- callNextMethod()
        n <- paste0(n, x@info$msLevel, ', ', floor(x@info$minRT), '-', ceiling(x@info$maxRT), ', ', floor(x@info$minMZ), '-', ceiling(x@info$maxMZ))
        n
    }
)

#' @describeIn MsIonList Create a plot
#' 
#' @param type Which type of plots should be plotted
#' 
#' @import ggplot2 gtable grid
#' @importFrom RColorBrewer brewer.pal
#' 
setMethod(
    'msPlot', 'MsIonList',
    function(object, type='2d', simple=FALSE, precursors=!simple, ...) {
        if(length(object) != 1) stop('Only plotting of single ionsets supported atm.')
        data <- data.frame(object@data)
        cScale <- brewer.pal(9, 'YlOrRd')
        if(type == '2d') {
            data <- binIons(data, ...)
            data$row <- 2
            data$col <- 1
            p <- ggplot(data, aes(x=retentionTime, y=mz, fill=intensity)) + theme_bw()
            p <- p + geom_raster()
            p <- p + scale_y_continuous('m/z') + scale_x_continuous('Retention time (sec)')
            p <- p + scale_fill_gradientn('Intensity', colours=cScale, na.value=rgb(0,0,0,0))
            if(precursors) {
                info <- msInfo(object)
                precursorData <- dbGetQuery(
                    con(object, 1)$sary(), 
                    paste0('SELECT prec.mz AS mz, orig.retentionTime AS retentionTime FROM 
                               (SELECT precursorScanNum, precursorMZ AS mz FROM currentHeader WHERE precursorScanNum IN (SELECT acquisitionNum FROM currentHeader WHERE retentionTime BETWEEN ', info$minRT, ' AND ', info$maxRT, ') AND precursorMZ BETWEEN ', info$minMZ, ' AND ', info$maxMZ, ')
                           AS prec 
                           JOIN currentHeader AS orig ON orig.acquisitionNum = prec.precursorScanNum')
                )
                if(nrow(precursorData) != 0) {
                    precursorData$row <- 2
                    precursorData$col <- 1
                    p <- p + geom_point(aes(colour=I('black'), fill=NULL), data=precursorData) + scale_colour_discrete('', breaks='black', labels='Precursor ions')
                }
            }
            if(!simple) {
                scanData <- data %>% 
                    group_by(mz) %>% 
                    summarise(TIC=sum(intensity, na.rm=T), col=2, row=2)
                chromData <- data.frame(object@data) %>% 
                    group_by(retentionTime) %>% 
                    summarise(TIC=sum(intensity, na.rm=T), col=1, row=1)
                p <- p + geom_line(aes(y=TIC, fill=NULL), data=chromData)
                p <- p + geom_segment(aes(yend=mz, x=TIC, xend=0, fill=NULL), data=scanData)
                p <- p + facet_grid(row~col, scales='free')
                
                p <- p + theme(strip.background=element_blank(), strip.text=element_blank())
                
                gt <- ggplot_gtable(ggplot_build(p))
                emptyGrobInd <- which(gt$layout$t == 4 & gt$layout$l == 6)
                gt$grobs[[emptyGrobInd]] <- grob()
                panels <- gt$layout$t[grep("panel", gt$layout$name)]
                gt$heights[panels] <- lapply(c(1,3), unit, "null")
                panels <- gt$layout$l[grep("panel", gt$layout$name)]
                gt$widths[panels] <- lapply(c(3, 3, 1, 1), unit, "null")
                gt$layout$r[gt$layout$name=='xlab'] <- 4
                gt$layout$t[gt$layout$name=='ylab'] <- 6
                gt <- gtable_add_grob(gt, textGrob('Intensity'), 8, 6)
                gt <- gtable_add_grob(gt, textGrob('Intensity', rot=90), 4, 2)
                
                plot(gt)
                invisible(gt)
            } else {
                p
            }
        } else if(type == '3d') {
            if(!require(rgl)) {
                stop('rgl package needed')
            }
            if(!rgl.useNULL()) open3d()
            clear3d()
            skip <- par3d(skipRedraw = TRUE)
            on.exit(par3d(skip))
            nIons <- nrow(data)
            maxInt <- max(data$intensity)
            cScale <- colorRampPalette(cScale)(256)
            endColours <- cScale[ceiling(256*data$intensity/maxInt)]
            colours <- rep(cScale[1], 2*nIons)
            colours[seq(1, by=2, length.out=nIons)] <- endColours
            segments3d(
                x = rep(data$retentionTime, rep(2, nIons)),
                y = rep(data$mz, rep(2, nIons)),
                z = rbind(rep(0, nIons), data$intensity),
                color = rbind(rep(cScale[1], nIons), endColours)
            )
            decorate3d(xlab='Retention time (sec)', ylab='m/z', zlab='Intensity', aspect=TRUE, box=FALSE, cex=3)
            bg3d(col='#abd9e9')   # #e0f3f8
        }
    }
)

#' Extract the scans that make up an ion list
#' 
setMethod(
    'scans', 'MsIonList',
    function(object, ...) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            iCon <- con(object, i, 'MsData')
            res[[i]] <- scans(iCon, msLevel=info$msLevel[i], retentionTime=BETWEEN(info$minRT[i], info$maxRT[i]))
        }
        res
    }
)

#' Extract chromatograms from an ion list
#' 
setMethod(
    'chroms', 'MsIonList',
    function(object, ...) {
        data <- lapply(msData(object), function(x) {
            chrom <- data.frame(x) %>%
                group_by(retentionTime) %>%
                select(intensity) %>%
                summarise(TIC=sum(intensity), BPC=max(intensity)) %>%
                arrange(retentionTime)
            as.matrix(chrom)
        })
        for(i in unique(object@mapping[, 'conIndex'])) {
            elIndex <- object@mapping[, 'conIndex'] == i
            acqNum <- dbGetPreparedQuery(
                con(object, i)$sary(), 
                'SELECT acquisitionNum, retentionTime 
                FROM header 
                WHERE retentionTime >= $minRT 
                    AND retentionTime <= $maxRT 
                    AND msLevel == $msLevel',
                bind.data=msInfo(object)[elIndex, , drop=FALSE]
            )
            data[elIndex] <- lapply(data[elIndex], function(x) {
                as.matrix(merge(acqNum, x))
            })
        }
        info <- data.frame(
            msLevel=msInfo(object)$msLevel,
            nScan=sapply(data, nrow), 
            maxTIC=sapply(data, function(x) max(x[,'TIC'])), 
            maxBPC=sapply(data, function(x) max(x[,'BPC'])), 
            minRT=msInfo(object)$minRT, 
            maxRT=msInfo(object)$maxRT,
            minMZ=msInfo(object)$minMZ,
            maxMZ=msInfo(object)$maxMZ
        )
        mapping <- getListMapping(data, object@mapping[, 'conIndex'])
        data <- do.call(rbind, data)
        new('MsChromList', connections=object@connections, info=info, data=data, mapping=mapping)
    }
)

#' Extract peaks from an MsIonList
#' 
setMethod(
    'peaks', 'MsIonList',
    function(object, overlap=FALSE) {
        info <- msInfo(object)
        res <- list()
        for(i in 1:length(object)) {
            iCon <- con(object, i, 'MsData')
            retRange <- ifelse(overlap, OVERLAPS(info$minRT[i], info$maxRT[i]), BETWEEN(info$minRT[i], info$maxRT[i]))
            mzRange <- ifelse(overlap, OVERLAPS(info$minMZ[i], info$maxMZ[i]), BETWEEN(info$minMZ[i], info$maxMZ[i]))
            res[[i]] <- peaks(iCon, msLevel=info$msLevel[i], retentionTime=retRange, mz=mzRange)
        }
        res
    }
)

#' @describeIn MsIonList Get information about the ion slice
#' 
#' @return \code{msInfo}: A data.frame with one row containing information
#' about the ion slice
#' 
setMethod(
    'msInfo', 'MsIonList',
    function(object) {
        callNextMethod()
    }
)

#' @describeIn MsIonList Get the ion data
#' 
#' @return \code{msData}: A list of matrices with acquisitionNum, retentionTime,
#' totIonCurrent and basePeakIntensity for each scan in the chromatogram
#' 
setMethod(
    'msData', 'MsIonList',
    function(object) {
        callNextMethod()
    }
)

