################################################################################
# TODO: Add analysis history to sary file - include timestamp and MSsary version
#
#       Implement R*tree for storing peaks
#
#       Register identification object with specific methods
#


#' @include MsConnections.R
#' @include tableFormats.R
#' @include generics.R
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
    function(object, id, ...) {
        if(length(list(...)) != 0) {
            acqNum <- getAcqNum(con(object), ...)
            if(!missing(id)) {
                acqNum <- id[id %in% acqNum]
            }
        } else {
            acqNum <- id
        }
        if(length(acqNum) == 0) {
            stop('No scans match criteria')
        }
        scInfo <- dbGetQuery(con(object)$sary(), paste0('SELECT * FROM header WHERE acquisitionNum IN (', paste(acqNum, collapse=', '), ')'))
        scData <- con(object)$getScans(scInfo$acquisitionNum)
        if(class(scData) == 'matrix') scData <- list(scData)
        mapping <- getListMapping(scData, 1)
        mapping <- cbind(mapping, matrix(isCentroided(scData), dimnames = list(NULL, 'mode')))
        scData <- do.call(rbind, scData)
        colnames(scData) <- c('mz', 'intensity')
        new('MsScanList', connections=list(con(object)), info=scInfo, data=scData, mapping=mapping)
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
    function(object, seqNum, retentionTime, mzRange, msLevel = 1) {
        #browser()
        if(missing(seqNum)) {
            if(missing(retentionTime)) {
                acqNum <- list(getAcqNum(con(object), msLevels = msLevel))
            } else {
                if(class(retentionTime) != 'matrix') {
                    retentionTime <- matrix(retentionTime, ncol=2, byrow=TRUE)
                }
                acqNum <- lapply(1:nrow(retentionTime), function(i) {
                    getAcqNum(con(object), retentionTime=retentionTime[i,], msLevels=msLevel)
                })
            }
        } else {
            if(!missing(retentionTime)) {
                warning('retention time window ignored')
            }
            if(class(seqNum) != 'matrix') {
                seqNum <- matrix(seqNum, ncol=2, byrow=TRUE)
            }
            acqNum <- lapply(1:nrow(seqNum), function(i) {
                getAcqNum(con(object), seqNum=seqNum[i,], msLevels=msLevel)
            })
        }
        if(missing(mzRange)) {
            scanNums <- sort(unique(do.call('c', acqNum)))
            data <- dbGetQuery(con(object)$sary(), paste0('SELECT acquisitionNum, retentionTime, totIonCurrent, basePeakIntensity FROM header WHERE acquisitionNum IN (', paste(scanNums, collapse=', '), ')'))
            #browser()
            data <- lapply(acqNum, function(x) {data[data$acquisitionNum %in% x,]})
        } else {
            if(class(mzRange) != 'matrix') {
                mzRange <- matrix(mzRange, ncol=2, byrow=TRUE)
            }
            if(length(acqNum) > nrow(mzRange)) {
                mzRange <- mzRange[rep(1:nrow(mzRange), length.out=length(acqNum)), ]
            } else if(length(acqNum) < nrow(mzRange)) {
                acqNum <- rep(acqNum, length.out=nrow(mzRange))
            }
            data <- con(object)$extractIC(acqNum, mzRange)
        }
        info <- lapply(1:length(data), function(i) {
            data.frame(
                name=NA, 
                msLevel=msLevel,
                nScan=nrow(data[[i]]), 
                maxTIC=max(data[[i]]$totIonCurrent), 
                maxBPC=max(data[[i]]$basePeakIntensity), 
                minRT=min(data[[i]]$retentionTime), 
                maxRT=max(data[[i]]$retentionTime)
            )
        })
        info <- do.call(rbind, info)
        if(!missing(mzRange)) {
            colnames(mzRange) <- c('minMZ', 'maxMZ')
            info <- cbind(info, mzRange)
        }
        info$name <- createChromNames(object, info)
        mapping <- getListMapping(data, 1)
        data <- as.matrix(do.call(rbind, data))
        cNames <- c('acquisitionNum', 'retentionTime', 'TIC', 'BPC')
        colnames(data) <- cNames
        new('MsChromList', connections=list(con(object)), info=info, data=data, mapping=mapping)
    }
)

#' Extract ions from an MsData object
#' 
setMethod(
    'ions', 'MsData',
    function(object, seqNum, retentionTime, mzRange, msLevel = 1, SIMPLIFY=TRUE) {
        if(missing(seqNum)) {
            if(missing(retentionTime)) {
                acqNum <- list(getAcqNum(con(object), msLevels = msLevel))
            } else {
                if(class(retentionTime) != 'matrix') {
                    retentionTime <- matrix(retentionTime, ncol=2, byrow=TRUE)
                }
                acqNum <- lapply(1:nrow(retentionTime), function(i) {
                    getAcqNum(con(object), retentionTime=retentionTime[i,], msLevels=msLevel)
                })
            }
        } else {
            if(!missing(retentionTime)) {
                warning('retention time window ignored')
            }
            if(class(seqNum) != 'matrix') {
                seqNum <- matrix(seqNum, ncol=2, byrow=TRUE)
            }
            acqNum <- lapply(1:nrow(seqNum), function(i) {
                getAcqNum(con(object), seqNum=seqNum[i,], msLevels=msLevel)
            })
        }
        if(missing(mzRange)) mzRange <- NULL
        ions <- con(object)$extractIons(acqNum, mzRange)
        ions <- lapply(ions, function(x) {
            new(
                'MsIonList',
                connections=list(con(object)),
                info=data.frame(
                    msLevel=msLevel,
                    minRT=min(x[, 'retentionTime']), 
                    maxRT=max(x[, 'retentionTime']), 
                    minMZ=min(x[, 'mz']),
                    maxMZ=max(x[, 'mz']),
                    minINT=min(x[, 'intensity']),
                    maxINT=max(x[, 'intensity'])
                ),
                data=x,
                mapping=matrix(c(1, nrow(x), 1), nrow=1, dimnames = list(NULL, c('start', 'end', 'conIndex')))
            )
        })
        if(SIMPLIFY && length(ions) == 1) ions <- ions[[1]]
        ions
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
#' @importFrom tools file_path_sans_ext
#' @importFrom mzR openMSfile close header
#' 
#' @export
#' 
createMsData <- function(rawFile) {
    dbFile <- paste0(file_path_sans_ext(rawFile,compression = T), '.sary')
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
    funcCall <- sys.call(0)
    connection$addData(
        'history', 
        data.frame(
            time = as.character(Sys.time()),
            operation = 'Sary created',
            MSsary_version = .MSsary_version,
            call = deparse(funcCall)
        )
    )
    
    connection$addTable('scans', scanTableFormat)
    
    connection$addTable('peaks', peakTableFormat)
    
    dbGetQuery(
        connection$sary(), 
        'CREATE VIEW currentRawHeader AS SELECT * FROM header WHERE acquisitionNum NOT IN (SELECT scanNum FROM scans WHERE retentionTime IS NULL)'
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
loadMsData <- function(rawFile, saryFile, testIntegrity=TRUE, testSize=10) {
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