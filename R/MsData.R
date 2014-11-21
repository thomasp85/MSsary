################################################################################
# TODO: Register identification object with specific methods
#           A unique ID defined for the id class is merged with either peaks or
#           links - this id is carried around and included in exports and can be
#           used to query the id object at a later time.
#       Create empty MsList when criteria is not met.
#
#       Copy method for MsData
#


#' @include MsConnections.R
#' @include tableFormats.R
#' @include generics.R
#' @include peakMethods.R
#' @include aaa.R
#' 
NULL

#' Main class for interacting with MS data
#' 
#' This class is the main entry point for interacting with MS data in the MSsary
#' package. It is designed to seamlessly integrate data from raw data files and
#' additional results from various algorithms. Furthermore it provides a mean to
#' commit persistent changes to the raw data (e.g. filtering or recalibration),
#' that will be available across sessions in a non-destructive manner.
#' 
#' @details
#' Behind the scene an SQLite-based file is created that handles all deviations
#' and additions to the raw data. The file is named after the raw data file but
#' has a .sary file extension to illustrate its relationship to this package.
#' Because it is just an SQLite database file the data generated with MSsary can
#' easily be accessed from other programs.
#' 
#' @section Constructors:
#' Objects of class MsData can be created either from scratch with 
#' \code{\link{createMsData}} or from an already created sary-file 
#' \code{\link{loadMsData}}.
#' 
#' @slot connections An MsConnections reference class object
#' 
#' @family MSsary-classes
#' 
#' @seealso \code{\link{createMsData}} \code{\link{loadMsData}}
#' 
setClass(
    'MsData',
    slots=list(
        connections='MsConnections'
    )
)

### METHODS

#' @describeIn MsData Short summary of content
#' 
#' @param object An MsData object
#' 
setMethod(
    'show', 'MsData',
    function(object) {
        cat('An MsData object\n')
        cat('\n')
        cat('         Raw file:', basename(con(object)$rawFile), '\n')
        cat('        Sary file:', basename(con(object)$saryFile), '\n')
        cat('\n')
        cat('  Number of scans:', length(object), '\n')
        nPeaks <- con(object)$nPeaks()
        if(nPeaks != 0) {
            cat('  Number of peaks:', nPeaks, '\n')
        }
        cat('\n')
        cat('Last edit history:\n')
        cat('\n')
        print(tail(editHistory(object)[,1:3]))
    }
)

#' @describeIn MsData Number of scans in the raw data
#' 
#' @param x An MsData object
#' 
setMethod(
    'length', 'MsData',
    function(x) {
        con(x)$length()
    }
)

setMethod(
    '==', c('MsData', 'MsData'),
    function(e1, e2) {
        con(e1) == con(e2)
    }
)
setMethod(
    '!=', c('MsData', 'MsData'),
    function(e1, e2) {
        !(e1 == e2)
    }
)

#' @describeIn con Get connection from MsData object
#' 
setMethod(
    'con', 'MsData',
    function(object) {
        return(object@connections)
    }
)

#' @describeIn MsData Get the edit history of the object
#' 
setMethod(
    'editHistory', 'MsData',
    function(object) {
        dbGetQuery(con(object)$sary(), 'SELECT * FROM history')
    }
)

#' Extract scans from an MsData object
#' 
#' This method returns one or several scans based on a range of selection 
#' criteria. See the arguments section for a full list.
#' 
setMethod(
    'scans', 'MsData',
    function(object, ..., raw=FALSE) {
        acqNum <- unlist(getAcqNum(con(object), ..., raw=raw))
        if(length(acqNum) == 0) {
            stop('No scans match criteria')
        }
        scInfo <- con(object)$getHeader(acqNum, raw=raw)
        scData <- con(object)$getScans(acqNum, raw=raw)
        mapping <- getListMapping(scData, 1, raw=as.integer(raw))
        mapping <- cbind(mapping, matrix(isCentroided(scData), dimnames = list(NULL, 'mode')))
        scData <- do.call(rbind, scData)
        colnames(scData) <- c('mz', 'intensity')
        new('MsScanList', connections=list(object), info=scInfo, data=scData, mapping=mapping)
    }
)

#' Extract chromatograms from an MsData object
#' 
#' This function is used to get chromatographic data from an MsData object. If
#' nothing besides the object is passed the full chromatogram is extracted, but
#' this can be altered by passing in scan, mz and retention time constraints.
#' The returned MsScanList object contains both TIC and BPC so there is no need 
#' to specify this during creation
#' 
setMethod(
    'chroms', 'MsData',
    function(object, ..., mz, raw=FALSE) {
        args <- list(...)
        acqNum <- getContAcqNum(con(object), ..., raw=raw)
        msLevels <- toFilter(args$msLevel)$data
        msLevels <- rep(msLevels, length.out=length(acqNum))
        if(missing(mz)) {
            scanNums <- sort(unique(do.call('c', acqNum)))
            data <- con(object)$getHeader(scanNums, raw=raw)[, c('acquisitionNum', 'retentionTime', 'totIonCurrent', 'basePeakIntensity')]
            data <- lapply(acqNum, function(x) {data[data$acquisitionNum %in% x,]})
        } else {
            if(!rangeFilter(mz)) {
                stop('mz must specifiy an interval: Use either \'BETWEEN\', \'ABOVE\' or \'BELOW\'')
            }
            mzFunc <- toFunction(mz)
            if(length(acqNum) > length(mzFunc)) {
                mzFunc <- mzFunc[rep(1:length(mzFunc), length.out=length(acqNum))]
            } else if(length(acqNum) < length(mzFunc)) {
                acqNum <- rep(acqNum, length.out=length(mzFunc))
            }
            data <- con(object)$extractIC(acqNum, mzFunc, raw=raw)
        }
        info <- lapply(1:length(data), function(i) {
            data.frame(
                msLevel=msLevels[i],
                nScan=nrow(data[[i]]), 
                maxTIC=max(data[[i]]$totIonCurrent), 
                maxBPC=max(data[[i]]$basePeakIntensity), 
                minRT=min(data[[i]]$retentionTime), 
                maxRT=max(data[[i]]$retentionTime)
            )
        })
        info <- do.call(rbind, info)
        if(!missing(mz)) {
            mzRange <- toRange(mz)
            mzRange <- mzRange[rep(1:nrow(mzRange), length.out=nrow(info)), , drop=FALSE]
            colnames(mzRange) <- c('minMZ', 'maxMZ')
            info <- cbind(info, mzRange)
        }
        mapping <- getListMapping(data, 1, raw=as.integer(raw))
        data <- as.matrix(do.call(rbind, data))
        cNames <- c('acquisitionNum', 'retentionTime', 'TIC', 'BPC')
        colnames(data) <- cNames
        new('MsChromList', connections=list(object), info=info, data=data, mapping=mapping)
    }
)

#' Extract ions from an MsData object
#' 
setMethod(
    'ions', 'MsData',
    function(object, ..., mz, raw=FALSE) {
        args <- list(...)
        acqNum <- getContAcqNum(con(object), ..., raw=raw)
        msLevels <- toFilter(args$msLevel)$data
        msLevels <- rep(msLevels, length.out=length(acqNum))
        
        if(missing(mz)) {
            mzFunc <- NULL
        } else {
            if(!rangeFilter(mz)) {
                stop('mz must specifiy an interval: Use either \'BETWEEN\', \'ABOVE\' or \'BELOW\'')
            }
            mzFunc <- toFunction(mz)
            if(length(acqNum) > length(mzFunc)) {
                mzFunc <- mzFunc[rep(1:length(mzFunc), length.out=length(acqNum))]
            } else if(length(acqNum) < length(mzFunc)) {
                acqNum <- rep(acqNum, length.out=length(mzFunc))
            }
        }
        data <- con(object)$extractIons(acqNum, mzFunc, raw=raw)
        info <- mapply(function(x, msLevel) {
            data.frame(
                msLevel=msLevel,
                minRT=min(x[, 'retentionTime']), 
                maxRT=max(x[, 'retentionTime']), 
                minMZ=min(x[, 'mz']),
                maxMZ=max(x[, 'mz']),
                minINT=min(x[, 'intensity']),
                maxINT=max(x[, 'intensity'])
            )
        }, x=data, msLevel=msLevels, SIMPLIFY=FALSE)
        info <- do.call(rbind, info)
        mapping <- getListMapping(data, 1, raw=as.integer(raw))
        data <- do.call(rbind, data)
        new('MsIonList', connections=list(object), info=info, data=data, mapping=mapping)
    }
)

#' Extrack peaks from an MsData object
#' 
setMethod(
    'peaks', 'MsData',
    function(object, ...) {
        ids <- unlist(getPeakIds(con(object), ...))
        info <- con(object)$getPeaks(ids)
        data <- lapply(info$peak, matrix, ncol=2, byrow=TRUE)
        info$peak <- NULL
        mapping <- getListMapping(data, 1)
        data <- do.call(rbind, data)
        colnames(data) <- c('retentionTime', 'intensity')
        new('MsPeakList', connections=list(object), info=info, data=data, mapping=mapping)
    }
)

#' Peak detection in MsData object
#' 
setMethod(
    'detectPeaks', 'MsData',
    function(object, method, ...) {
        arguments <- list(...)
        dataSubset <- names(arguments) %in% names(formals(getContAcqNum))
        acqNum <- unlist(do.call(getContAcqNum, c(con=con(object), arguments[dataSubset])))
        
        arguments <- arguments[!dataSubset]
        arguments$name <- method
        arguments$scans <- con(object)$getScans(acqNum)
        arguments$info <- con(object)$getHeader(acqNum)
        
        peaks <- do.call(peakMethods$useMethod, arguments)
        funcCall <- expand.call(call=sys.call(-1), expand.dots = TRUE)
        prevID <- dbGetQuery(con(object)$sary(), 'SELECT IFNULL(MAX(rowid), 0) AS max FROM peakInfo')$max
        con(object)$addData('peakInfo', peaks[, c('msLevel', 'length', 'mzMean', 'maxHeight', 'area', 'peak')])
        ids <- dbGetQuery(con(object)$sary(), paste0('SELECT rowid AS id FROM peakInfo WHERE rowid > ', prevID))$id
        con(object)$addData('peakLoc', data.frame(peakID=ids, peaks[, c('scanStart', 'scanEnd', 'mzMin', 'mzMax')]))
        con(object)$addData(
            'history', 
            data.frame(
                time = as.character(Sys.time()),
                operation = 'Peaks detected',
                MSsary_version = as.character(packageVersion('MSsary')),
                call = callToString(funcCall),
                augPackage = peakMethods$getPackage(method),
                augPackVersion = peakMethods$getVersion(method),
                stringsAsFactors=FALSE
            )
        )
        invisible(TRUE)
    }
)

### CONSTRUCTORS

#' Create an MsData object from a raw MS data file
#' 
#' This function takes an mzR-compliant data file and create the initial, 
#' minimal sary database file and set up the link to both the database and raw
#' file.
#' 
#' @param rawFile The path to an mzR-compliant MS data file
#' 
#' @return An MsData object
#' 
#' @seealso \code{\link{loadMsData}} \code{\linkS4class{MsData}}
#' 
#' @examples
#' file <- 'test'
#' msdata <- createMsData(file)
#' 
#' @importFrom tools file_path_sans_ext file_path_as_absolute
#' @importFrom mzR openMSfile close header
#' 
#' @export
#' 
createMsData <- function(rawFile, saryName, force=FALSE) {
    if(missing(saryName)) {
        dbFile <- paste0(file_path_sans_ext(rawFile,compression = T), '.sary')
    } else {
        dbFile <- file.path(dirname(rawFile), paste0(saryName, '.sary'))
    }
    if(force) {
        if(unlink(dbFile) == 1) {
            stop('Cannot remove ', dbFile)
        }
    }
    if(file.exists(dbFile)) {
        stop('Sary file already exists at: ', dbFile)
    }
    file.create(dbFile)
    
    connection <- MsConnections(rawFile, dbFile)
    
    connection$addTable('header', headerTableFormat)
    
    raw <- openMSfile(rawFile)
    headerInfo <- header(raw)
    close(raw)
    
    connection$addData('header', headerInfo)
    
    connection$addTable('history', historyTableFormat)
    funcCall <- match.call()
    connection$addData(
        'history', 
        data.frame(
            time = as.character(Sys.time()),
            operation = 'Sary created',
            MSsary_version = as.character(packageVersion('MSsary')),
            call = callToString(funcCall),
            stringsAsFactors=FALSE
        )
    )
    connection$addTable('mzR', 'location TEXT NOT NULL')
    connection$addData('mzR', data.frame(location=file_path_as_absolute(rawFile)))
    
    connection$addTable('scans', scanTableFormat)
    
    connection$addTable('peakInfo', peakTableFormat)
    connection$addTable('peakLoc', peakRtree, rtree=TRUE)
    
    dbGetQuery(
        connection$sary(), 
        'CREATE VIEW currentRawHeader AS SELECT * FROM header WHERE acquisitionNum NOT IN (SELECT scanNum FROM scans)'
    )
    dbGetQuery(
        connection$sary(), 
        'CREATE VIEW currentModHeader AS SELECT h.seqNum, h.acquisitionNum, h.msLevel, h.polarity, s.peaksCount, s.totIonCurrent, h.retentionTime, s.basePeakMZ, s.basePeakIntensity, h.collisionEnergy, h.ionisationEnergy, s.lowMZ, s.highMZ, h.precursorScanNum, h.precursorMZ, h.precursorCharge, h.precursorIntensity, h.mergedScan, h.mergedResultScanNum, h.mergedResultStartScanNum, h.mergedResultEndScanNum FROM header AS h JOIN scans AS s ON h.acquisitionNum=s.scanNum AND s.remove == 0'
    )
    dbGetQuery(
        connection$sary(),
        'CREATE VIEW currentHeader AS SELECT * FROM currentRawHeader UNION SELECT * FROM currentModHeader ORDER BY seqNum'
    )
    
    new('MsData', connections=connection)
}

#' Load in an already created sary file as an MsData object
#' 
#' This function creates a new MsData object from an already created sary file,
#' together with the accompanying raw data file. The function automatically 
#' checks the sary file and see if it complies with the format. By default it 
#' also does a light check on whether the data in the raw file corresponds to 
#' the data in the sary file.
#' 
#' @param rawFile The path to the raw datafile
#' 
#' @param saryFile The path to the sary database file
#' 
#' @param testIntegrity Should the raw and sary file be compared to each other
#' 
#' @param testSize The number of randomly selected scans to use for comparing
#' raw and sary if \code{testIntegrity=TRUE}
#' 
#' @return An MsData object
#' 
#' @seealso \code{\link{createMsData}} \code{\linkS4class{MsData}}
#' 
#' @export
#' 
loadMsData <- function(saryFile, rawFile, testIntegrity=TRUE, testSize=10) {
    if(missing(rawFile)) {
        rawFile <- getMzrPath(saryFile)
        if(!file.exists(rawFile)) {
            rawFile <- file.path(dirname(saryFile), basename(rawFile))
            if(!file.exists(rawFile)) {
                stop('Cannot infer raw data location')
            }
        }
    }
    connection <- MsConnections(rawFile, saryFile)
    
    if(!connection$verify()) {
        stop(saryFile, ' does not look like a saryFile')
    }
    if(testIntegrity) {
        if(!connection$validate(testSize)) {
            stop(saryFile, ' and ', rawFile, ' does not match')
        }
    }
    new('MsData', connections=connection)
}