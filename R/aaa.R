#' Get the current version of a package
#' 
#' This function reads the version from a package from the description file
#' 
#' @param package The package to get the version for
#' 
#' @return A string with the version
#' 
#' @noRd
#' 
getVersion <- function(package = 'MSsary') {
    desc <- readLines(system.file('DESCRIPTION', package=package))
    vers <- desc[grep('Version:', desc, ignore.case = TRUE)]
    vers <- sub('Version: ', '', vers, ignore.case = TRUE)
    return(vers)
}
.MSsary_version <- getVersion()

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
getAcqNum <- function(con, seqNum, retentionTime, TIC, BPC, nPeaks, msLevels, isParent, isChild) {
    query <- c()
    if(!missing(seqNum)) {
        if(class(seqNum) != 'matrix') {
            seqNum <- matrix(seqNum, ncol=2, byrow=TRUE)
        }
        query <- c(query, constructInterval('seqNum', seqNum))
    }
    if(!missing(retentionTime)) {
        if(class(retentionTime) != 'matrix') {
            retentionTime <- matrix(retentionTime, ncol=2, byrow=TRUE)
        }
        query <- c(query, constructInterval('retentionTime', retentionTime))
    }
    if(!missing(TIC)) {
        if(class(TIC) != 'matrix') {
            TIC <- matrix(TIC, ncol=2, byrow=TRUE)
        }
        query <- c(query, constructInterval('TIC', TIC))
    }
    if(!missing(BPC)) {
        if(class(BPC) != 'matrix') {
            BPC <- matrix(BPC, ncol=2, byrow=TRUE)
        }
        query <- c(query, constructInterval('BPC', BPC))
    }
    if(!missing(nPeaks)) {
        if(class(nPeaks) != 'matrix') {
            nPeaks <- matrix(nPeaks, ncol=2, byrow=TRUE)
        }
        query <- c(query, constructInterval('nPeaks', nPeaks))
    }
    if(!missing(msLevels)) {
        query <- c(query, paste0('(msLevel IN (', paste(msLevels, collapse=', '),'))'))
    }
    if(!missing(isParent)) {
        subQuery <- '(SELECT precursorScanNum FROM header)'
        if(isParent) {
            query <- c(query, paste0('(acquisitionNum IN ', subQuery, ')'))
        } else {
            query <- c(query, paste0('(acquisitionNum NOT IN ', subQuery, ')'))
        }
    }
    if(!missing(isChild)) {
        if(isChild) {
            query <- c(query, paste0('(precursorScanNum != 0)'))
        } else {
            query <- c(query, paste0('(precursorScanNum == 0)'))
        }
    }
    query <- paste(query, collapse=' AND ')
    if(query == '') {
        query <- 'SELECT acquisitionNum FROM header'
    } else {
        query <- paste0('SELECT acquisitionNum FROM header WHERE (', query, ')')
    }
    dbGetQuery(con$sary(), query)$acquisitionNum
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

createChromNames <- function(object, info) {
    name <- file_path_sans_ext(basename(con(object)$saryFile))
    name <- paste0(name, ': ', floor(info$minRT), '-', ceiling(info$maxRT))
    if(!is.null(info$minMZ)) {
        name <- paste0(name, ', ', floor(info$minMZ), '-', ceiling(info$maxMZ))
    }
    name
}