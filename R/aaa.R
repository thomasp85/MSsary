#' Asses the mode of a spectrum
#' 
#' This function checks whether a spectrum is present in profile or centroid
#' mode using a simple heuristic. It checks the m/z differences between 
#' consecutive ions and if the 0.25 quantile is below 0.025 it is taken as a 
#' profile spectrum
#' 
#' @param scans A list of matrices containing scan data
#' 
#' @return A logical vector with same length as the input containing TRUE if the
#' corresponding spectrum is centroided
#' 
#' @noRd
#' 
isCentroided <- function(scans) {
    if(class(scans) == 'matrix') scans <- list(scans)
    
    sapply(scans, function(x) {
        if(length(x) == 0) return(NA)
        if(nrow(x) < 10) return(TRUE)
        
        quantile(diff(x[,1]), 0.25) > 0.025
    })
}

#' Create a matrix of start and end indexes for an rbinded scan list
#' 
#' This function calculates the corresponding start and end indexes for each
#' scan in a list when the list is rbinded. Additionally it attaches the 
#' connection index to the output.
#' 
#' @param data A list of matrices with scan info
#' 
#' @param conIndex The index for the connection for each scan
#' 
#' @return A matrix with number of rows corresponding to the length of data, 
#' with columns 'start', 'end' and 'conIndex'
#' 
#' @noRd
#' 
getListMapping <- function(data, conIndex) {
    dataLengths <- sapply(data, nrow)
    mapping <- createIntervals(dataLengths, list(conIndex=conIndex))
    return(mapping)
}

#' Subset an MsList mapping matrix
#' 
#' This function recalculates the start and end intervals after subsetting the
#' data part of an MsList
#' 
#' @param mapping The mapping matrix from an MsList subclass
#' 
#' @param i The indexes of the elements that are going to be kept
#' 
#' @return A matrix with the same columns as mapping and number of rows 
#' corresponding to the length of i
#' 
#' @noRd
#' 
subsetMapping <- function(mapping, i) {
    extraData <- mapping[i, !(colnames(mapping) %in% c('start', 'end')), drop=FALSE]
    elementLengths <- getElementLengths(mapping)
    createIntervals(elementLengths[i], extraData)
}

#' Calculate indexing interval based on lengths
#' 
#' This function calculates the start and end indexes based on a vector of
#' lengths and optinally attaches the data given in extraData.
#' 
#' @param lengths A vector of integers
#' 
#' @param extraData Either a matrix with number of row corresponding to the 
#' length of lengths, or a list with named elements of type vector. Values in 
#' the vectors will be recycled if necessry.
#' 
#' @return A matrix withnumber of rows corresponding to the length of lengths,
#' containing the columns 'start' and 'end' as well as whatever was passed on
#' through extraData
#' 
#' @noRd
#' 
createIntervals <- function(lengths, extraData) {
    mapping <- matrix(NA, ncol=2, nrow=length(lengths))
    currentIndex <- 0
    for(i in 1:length(lengths)) {
        if(lengths[i] == 0) {
            start <- NA
            end <- NA
        } else {
            start <- currentIndex+1
            end <- currentIndex+lengths[i]
            currentIndex <- end
        }
        mapping[i, ] <- c(start, end)
    }
    colnames(mapping) <- c('start', 'end')
    if(!missing(extraData)) {
        if(class(extraData) == 'list') {
            extraData <- mapply(function(name, data) {
                data <- rep(data, length.out = nrow(mapping))
                matrix(data, dimnames=list(NULL, name))
            }, names(extraData), extraData, SIMPLIFY=FALSE)
            extraData <- do.call(cbind, extraData)
        }
        if(class(extraData) != 'matrix') {
            stop('Can only handle extraData of class \'matrix\' or \'list\'')
        }
        mapping <- cbind(mapping, extraData)
    }
    mapping
}

#' Convert a mapping matrix into a vector of element legnths
#' 
#' This function calculates the legnth of each element based on a matrix with
#' 'start' and 'end' columns
#' 
#' @param mapping A mapping matrix from a MsList
#' 
#' @return A vector of integers corresponding to the length of each element in
#' the mapping
#' 
#' @noRd
#' 
getElementLengths <- function(mapping) {
    elementLengths <- apply(mapping, 1, function(x) {
        x['end'] - x['start'] + 1
    })
    elementLengths[is.na(elementLengths)] <- 0
    elementLengths
}

#' Convert a mapping matrix into a vector of indexes
#' 
#' This function creates a vector with the same length as the data slot in an 
#' MsList subclass based on its mapping matrix. The vector contains for each
#' index the corresponding element index
#' 
#' @param mapping A mapping matrix as stored in an MsList subclass
#' 
#' @return A vector with length corresponding to the maximum 'end' value in the 
#' mapping matrix (i.e. the number of rows of the data slot in the object)
#' 
#' @noRd
#' 
getElementIndex <- function(mapping) {
    elementLengths <- getElementLengths(mapping)
    rep(1:nrow(mapping), elementLengths)
}

#' Extract scans based on a query returning acquisitionNum'
#' 
#' This function iterates over the scans in a MsScanList and extract new scan(s)
#' from their connections based on the corresponding query. The query have to 
#' return a data.frame with an acquisitionNum column (though the number of rows
#' can be zero)
#' 
#' @param object An MsScanList object
#' 
#' @param query A vector of strings with valid SQL queries with the same length
#' as object. The query should in general extract from the header table and 
#' include the acquisitionNum column
#' 
#' @return A new MsScanList object
#' 
#' @noRd
#' 
getScanList <- function(object, query) {
    scInfo <- msInfo(object)
    info <- list()
    data <- list()
    for(i in 1:nrow(scInfo)) {
        info[[i]] <- dbGetQuery(con(object, i)$sary(), query[i])
        if(nrow(info[[i]]) == 0) {
            info[[i]][1,] <- NA
            data[[i]] <- list(matrix(ncol=2, nrow=0))
        } else {
            data[[i]] <- con(object, i)$getScans(info[[i]]$acquisitionNum)
        }
    }
    info <- do.call(rbind, info)
    data <- unlist(data, recursive=FALSE)
    mapping <- getListMapping(data, object@mapping[, 'conIndex'])
    data <- do.call(rbind, data)
    colnames(data) <- c('mz', 'intensity')
    new('MsScanList', connections=object@connections, info=info, data=data, mapping=mapping)
}

#' Get the acquisition numbers based on a filter
#' 
#' This function returns the acquisition numbers from an MsData object that 
#' passes a set of filters based on the metadata available.
#' 
#' @param con An MsConnections object
#' 
#' @param seqNum A matrix of intervals
#' 
#' @param retentionTime A matrix of intervals
#' 
#' @param TIC A matrix of intervals
#' 
#' @param BPC A matrix of intervals
#' 
#' @param nPeaks A matrix of intervals
#' 
#' @param msLevels An integer vector
#' 
#' @param isParent Logical
#' 
#' @param isChild Logical
#' 
#' @return A vector of acquisition numbers
#' 
#' @noRd
#' 
getAcqNum <- function(con, ..., raw=FALSE) {
    table <- ifelse(raw, 'header', 'currentHeader')
    filter <- 'none'
    query <- lapply(list(...), toFilter)
    if(length(query) != 0) {
        filter <- 'some'
        query <- mapply(fillSQL, filter=query, name=names(query), SIMPLIFY=FALSE)
        query <- do.call(combineSQL, query)
    }
    switch(
        filter,
        none = list(as.integer(dbGetQuery(con$sary(), paste0('SELECT acquisitionNum FROM ', table))$acquisitionNum)),
        some = lapply(query, function(x) {
            as.integer(dbGetQuery(con$sary(), paste0('SELECT acquisitionNum FROM ', table, ' ', x))$acquisitionNum)
        })
    )
}

#' Extract a continuous range of acquisitionNums from the same MS level
#' 
#' This function works as getAcqNum, but are constrained to only extract a 
#' single continuous interval from the same MS level such as needed when 
#' extracting peaks etc.
#' 
getContAcqNum <- function(con, seqNum, retentionTime, msLevel) {
    msLevel <- toFilter(msLevel)
    if(msLevel$type != 'EQUALS') {
        stop('Only one msLevel can be selected')
    }
    if(!missing(seqNum) && !missing(retentionTime)) {
        stop('Either use seqNum or retentionTime - not both')
    }
    arguments <- list(con=con, msLevel=msLevel)
    if(!missing(seqNum)) {
        if(!rangeFilter(seqNum)) {
            stop('seqNum must specifiy an interval: Use either \'BETWEEN\', \'ABOVE\' or \'BELOW\'')
        }
        arguments$seqNum <- seqNum
    }
    if(!missing(retentionTime)) {
        if(!rangeFilter(retentionTime)) {
            stop('retentionTime must specifiy an interval: Use either \'BETWEEN\', \'ABOVE\' or \'BELOW\'')
        }
        arguments$retentionTime <- retentionTime
    }
    do.call(getAcqNum, arguments)
}

#' Get peakIDs based on a filter
#' 
#' Similar to getAcqNum but used to extract peakID instead
#' 
getPeakIds <- function(con, seqNum, retentionTime, mz, ...) {
    filters <- 'none'
    locQuery <- list()
    if(!missing(retentionTime)) {
        nFilter <- rtToSeqFilter(retentionTime, con)
        locQuery <- c(locQuery, list(fillSQL(nFilter, 'scanStart', 'scanEnd')))
    }
    if(!missing(seqNum)) {
        locQuery <- c(locQuery, list(fillSQL(seqNum, 'scanStart', 'scanEnd')))
    }
    if(!missing(mz)) {
        locQuery <- c(locQuery, list(fillSQL(mz, 'mzMin', 'mzMax')))
    }
    infoQuery <- lapply(list(...), toFilter)
    
    if(length(locQuery) != 0) {
        filters <- 'location'
        locQueryStrings <- do.call(combineSQL, locQuery)
    }
    if(length(infoQuery) != 0) {
        filters <- 'info'
        infoQuery <- query <- mapply(fillSQL, filter=infoQuery, name=names(infoQuery), SIMPLIFY=FALSE)
        infoQueryStrings <- do.call(combineSQL, infoQuery)
    }
    if(length(locQuery) != 0 && length(infoQuery) != 0) {
        filters <- 'both'
        if(length(locQueryStrings) < length(infoQueryStrings)) {
            locQueryStrings <- do.call(combineSQL, c(locQuery, nSQL=length(infoQueryStrings)))
        } else {
            infoQueryStrings <- do.call(combineSQL, c(infoQuery, nSQL=length(locQueryStrings)))
        }
    }
    
    switch(
        filters,
        none = list(as.integer(dbGetQuery(con$sary(), 'SELECT peakID FROM peakLoc')$peakID)),
        location = lapply(locQueryStrings, function(x) {
            as.integer(dbGetQuery(con$sary(), paste0('SELECT peakID FROM peakLoc ', x))$peakID)
        }),
        info = lapply(infoQueryStrings, function(x) {
            as.integer(dbGetQuery(con$sary(), paste0('SELECT peakID FROM peakInfo ', x))$peakID)
        }),
        both = {
            locID <- lapply(locQueryStrings, function(x) {
                as.integer(dbGetQuery(con$sary(), paste0('SELECT peakID FROM peakLoc ', x))$peakID)
            })
            infoID <- lapply(infoQueryStrings, function(x) {
                as.integer(dbGetQuery(con$sary(), paste0('SELECT peakID FROM peakInfo ', x))$peakID)
            })
            mapply(intersect, locID, infoID, SIMPLIFY=FALSE)
        }
    )
}

#' Create a set of intervals based on a two column matrix
#' 
#' This helper function builds up an SQL subexpression that selects 'name' in
#' the intervals defined by the first and second columns in the interval matrix.
#' Column one defines the lower bounds and column 2 the upper bound (bound 
#' included). Each row in interval defines an interval and all intervals are 
#' OR'ed together so the union of intervals are obtained. The returned string is
#' surrounded by parenthesis so it is ready to be combined with other strings in
#' a WHERE clause
#' 
#' @param name The name of the column for which the interval should be created
#' 
#' @param interval a matrix with 2 columns. Column 1 sets the lower bound, 
#' column 2 the upper. NA in either means an open bound.
#' 
#' @return A string
#' 
#' @noRd
#' 
constructInterval <- function(name, interval) {
    query <- c()
    for(i in 1:nrow(interval)) {
        int <- interval[i,]
        if(is.na(int[1])) {
            if(is.na(int[2])) {
                # Skip empty interval
            } else {
                query[i] <- paste0('(', name, '<=', int[2], ')')
            }
        } else {
            if(is.na(int[2])) {
                query[i] <- paste0('(', name, '>=', int[1], ')')
            } else {
                query[i] <- paste0('(', name, '>=', int[1], ' AND ', name, '<=', int[2], ')')
            }
        }
    }
    query <- query[!is.na(query)]
    paste0('(', paste(query, collapse=' OR '), ')')
}

#' Create profile data from ions
#' 
#' This function takes a matrix of ion triples and bins the data to create well
#' formated data
#' 
#' @import dplyr
#' @importFrom reshape2 acast melt
#' 
binIons <- function(ions, fun=max, mzBins=500, rtBins=mzBins) {
    rtRange <- range(ions[, 'retentionTime'])
    mzRange <- range(ions[, 'mz'])
#    rtBins <- ifelse(missing(rtBins), length(unique(ions[, 'retentionTime']))-1, rtBins)
    binWidthRT <- diff(rtRange)/rtBins
    binWidthMZ <- diff(mzRange)/mzBins
    rtBreaks <- seq(rtRange[1]-binWidthRT/2, rtRange[2]+binWidthRT/2, binWidthRT)
    mzBreaks <- seq(mzRange[1], mzRange[2], binWidthMZ)
    
    rtBin <- cut(ions[, 'retentionTime'], rtBreaks, include.lowest=TRUE)
    mzBin <- cut(ions[, 'mz'], mzBreaks, include.lowest=TRUE)
    
    ans <- data.frame(ions, rtBin=rtBin, mzBin=mzBin) %>%
        group_by(rtBin, mzBin) %>%
        select(intensity) %>%
        summarise(intensity=fun(intensity))
    
#    allComb <- expand.grid(rtBin=levels(rtBin), mzBin=levels(mzBin))
#    ans <- merge(allComb, as.data.frame(ans), all.x=TRUE)
    
    ans <- acast(ans, mzBin ~ rtBin, fill=NA, value.var='intensity', drop=FALSE)
    gaps <- findGaps(ans)
    if(!is.null(gaps)) {
        for(i in 1:length(gaps)) {
            gap <- gaps[[i]]
            bounds <- range(gap) + c(-1, 1)
            ans[, gap] <- fillGaps(ans[, bounds[1]], ans[, bounds[2]], length(gap))
        }
    }
    ans <- melt(ans, varnames=c('mzBin', 'rtBin'), value.name='intensity')
    
    data.frame(
        retentionTime = (rtBreaks[1:rtBins]+binWidthRT/2)[match(ans$rtBin, levels(rtBin))],
        mz = (mzBreaks[1:mzBins]+binWidthMZ/2)[match(ans$mzBin, levels(mzBin))],
        intensity = ans$intensity
    )
}

#' Interpolate gaps between RT bins
#' 
fillGaps <- function(from, to, width) {
    naRows <- which(is.na(from) & is.na(to))
    from[is.na(from)] <- 0
    to[is.na(to)] <- 0
    ans <- t(mapply(seq, from, to, MoreArgs=list(length.out=width+2)))
    ans <- ans[, -c(1, width+2), drop=FALSE]
    ans[naRows, ] <- NA
    ans
}
#' Find gaps in RT bins
#' 
findGaps <- function(data) {
    naCols <- apply(data, 2, function(x) all(is.na(x)))
    if(!any(naCols)) return(NULL)
    breaks <- c(0, which(diff(which(naCols)) != 1), length(which(naCols)))
    sapply(seq(length(breaks) - 1), function(i) which(naCols)[(breaks[i] + 1):breaks[i+1]])
}

#' Create a data.frame with precursor scans
#' 
annotateChildren <- function(data, children, mode) {
    if(mode == 'profile') {
        data <- data[which(diff(sign(diff(data$intensity)))==-2)+1, ]
    }
    diffs <- abs(outer(data$mz, children$precursorMZ, `-`))
    cInd <- apply(diffs, 2, which.min)
    data[cInd,]
}

#' Extract mzR file location from sary file
#' 
getMzrPath <- function(saryPath) {
    db <- dbConnect(dbDriver('SQLite'), saryPath)
    rawPath <- dbGetQuery(db, 'SELECT location FROM mzR')$location
    dbDisconnect(db)
    rawPath
}
#' Return a call in which all of the arguments which were supplied
#' or have presets are specified by their full names and supplied
#' or default values.
#'  
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition.
#'  See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be 
#'  included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @param doEval logical, defaults to TRUE. Should function arguments be 
#'  evaluated in the returned call or not?
#'
#' @return An object of class call. 
#' @author fabians
#' @seealso \code{\link[base]{match.call}}
#' @references \href{http://stackoverflow.com/questions/3478923/displaying-the-actual-parameter-list-of-the-function-during-execution}{Stack Overflow answer}
#' 
expand.call <- function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE, doEval=TRUE) {
    safeDeparse <- function(expr){
        #rm line breaks, whitespace             
        ret <- paste(deparse(expr), collapse="")
        return(gsub("[[:space:]][[:space:]]+", " ", ret))
    }
    
    call <- .Internal(match.call(definition, call, expand.dots))
    
    #supplied args:
    ans <- as.list(call)
    if(doEval) {
        for(i in 2:length(ans)) {
            ans[[i]] <- eval(ans[[i]])
        }
    }
    
    #possible args:
    frmls <- formals(safeDeparse(ans[[1]]))
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
}

#' Deparse call into a single text string
#' 
callToString <- function(call) {
    paste(sub('^\\s+', '', deparse(call)), collapse='')
}