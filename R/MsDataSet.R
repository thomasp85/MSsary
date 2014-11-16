################################################################################
# TODO: members table should be queriable during extraction
#

#' @include MsData.R
#' @include MsSet.R
#' @include generics.R
#' @include aaa.R
NULL

#' Store and operate on sets of samples
#' 
#' This object handles multiple MsData sets and facilitates batch analysis on
#' these. 
#' 
setClass(
    'MsDataSet',
    slots = list(
        setSary='MsSet',
        connections='list'
    ),
    validity = function(object) {
        if(!all(sapply(object@connections, class) == 'MsData')) {
            return('All connections must be MsData objects')
        }
        return(TRUE)
    }
)


## METHODS
setMethod(
    'length', 'MsDataSet',
    function(x) {
        length(x@connections)
    }
)
setMethod(
    'show', 'MsDataSet',
    function(object) {
        cat('An MsDataSet with', length(object), 'members')
    }
)
setMethod(
    '[[', c('MsDataSet', 'numeric', 'missing'),
    function(x, i, j) {
        x@connections[[i]]
    }
)
setMethod(
    '[[', c('MsDataSet', 'character', 'missing'),
    function(x, i, j) {
        i <- which(names(x@connections) == i)
        x[[i]]
    }
)
setMethod(
    'scans', 'MsDataSet',
    function(object, ...) {
        scans <- lapply(object@connections, scans, ...)
        for(i in 1:length(scans)) {
            names(scans[[i]]@connections) <- names(scans)[i]
        }
        names(scans) <- NULL
        do.call(c, scans)
    }
)
setMethod(
    'chroms', 'MsDataSet',
    function(object, ...) {
        chroms <- lapply(object@connections, chroms, ...)
        for(i in 1:length(chroms)) {
            names(chroms[[i]]@connections) <- names(chroms)[i]
        }
        names(chroms) <- NULL
        do.call(c, chroms)
    }
)
setMethod(
    'peaks', 'MsDataSet',
    function(object, ...) {
        peaks <- lapply(object@connections, peaks, ...)
        for(i in 1:length(peaks)) {
            names(peaks[[i]]@connections) <- names(peaks)[i]
        }
        names(peaks) <- NULL
        do.call(c, peaks)
    }
)
setMethod(
    'ions', 'MsDataSet',
    function(object, ...) {
        ions <- lapply(object@connections, ions, ...)
        for(i in 1:length(ions)) {
            names(ions[[i]]@connections) <- names(ions)[i]
        }
        names(ions) <- NULL
        do.call(c, ions)
    }
)

setMethod(
    'detectPeaks', 'MsDataSet',
    function(object, method, ...) {
        for(i in 1:length(object)) {
            detectPeaks(object[[i]], method, ...)
        }
        invisible(NULL)
    }
)


## CONSTRUCTORS

#' Create an MsDataSet from a list of raw files
#' 
#' @export
#' 
createMsDataSet <- function(rawFiles, setName, memberNames, ref, design, force=FALSE) {
    saryFiles <- paste0(file_path_sans_ext(basename(rawFiles), compression = T), '_', setName)
    connections <- list()
    for(i in 1:length(rawFiles)) {
        connections[[i]] <- createMsData(rawFiles[i], saryFiles[i], force=force)
    }
    saryFiles <- file.path(dirname(rawFiles), paste0(saryFiles, '.sary'))
    
    setFile <- file.path(dirname(rawFiles[1]), paste0(setName, '.saryset'))
    if(force) {
        if(unlink(setFile) == 1) {
            stop('Cannot remove ', setFile)
        }
    }
    if(file.exists(setFile)) {
        stop('Saryset file already exists at: ', setFile)
    }
    file.create(setFile)
    
    saryset <- MsSet(setFile)
    
    if(missing(memberNames)) {
        memberNames <- basename(file_path_sans_ext(saryFiles, compression = T))
    }
    names(connections) <- memberNames
    if(missing(ref)) {
        ref <- memberNames[floor(length(memberNames)/2)]
        message(ref, ' is set as reference')
    }
    if(is.integer(ref) || is.numeric(ref)) {
        ref <- memberNames[ref]
    }
    memberInfo <- data.frame(saryLocation=saryFiles, memberName=memberNames, referenceSample=as.integer(memberNames==ref))
    memberTableDef <- memberTableFormat
    if(!missing(design)) {
        if(class(design) != 'data.frame') stop('Design must be supplied as a data.frame')
        if(nrow(design) != length(saryFiles)) stop('Dimensions of design must match number of samples')
        testCon <- dbConnect(dbDriver('SQLite'), ':memory:')
        designDef <- sapply(mtcars, function(x) {dbDataType(testCon, x)})
        dbDisconnect(testCon)
        memberTableDef <- c(memberTableDef, paste(names(designDef), designDef))
        memberInfo <- cbind(memberInfo, design)
    }
    
    saryset$addTable('members', memberTableDef)
    saryset$addData('members', memberInfo)
    
    saryset$addTable('groupLoc', groupRtree, rtree=TRUE)
    saryset$addTable('groupInfo', groupTableFormat)
    saryset$addTable('groups', groupMemberFormat)
    saryset$addTable('rtCor', rtCorTableFormat)
    
    new('MsDataSet', setSary=saryset, connections=connections)
}
#' Load an already created MsDataSet from a saryset file
#' 
#' @export
#' 
loadMsDataSet <- function(setFile, ...) {
    saryset <- MsSet(setFile)
    
    info <- saryset$getInfo()
    
    connections <- list()
    for(i in 1:nrow(info)) {
        connections[[i]] <- loadMsData(info$saryLocation[i], ...)
    }
    names(connections) <- info$memberName
    
    new('MsDataSet', setSary=saryset, connections=connections)
}