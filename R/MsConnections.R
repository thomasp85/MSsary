#' A class to handle MS I/O
#' 
#' This class takes care of the communication with the sary database and the raw
#' datafile. Its main responsibility is to ensure a stable connection across
#' sessions without causing the session to crash (This would normally happen if
#' an mzRramp object was accessed in a different session than the one it is 
#' created in). It relies on knowledge of the location of the raw and database 
#' files and will throw an error if these files are moved.
#' 
#' TODO: Schema for updating file location vs. creating a new object (which is 
#' cheap). Decide whether this class also should handle writing to the sary db.
#' 
#' @field rawFile The path to the location of the MS raw data file
#' 
#' @field ramp Object of class mzRramp. A connection to the MS raw data file
#' 
#' @field saryFile The path to the location of the sary database file
#' 
#' @field saryDB An SQLiteConnection object. A connection to the sary database
#' 
#' @seealso \code{\linkS4class{MsData}}
#' 
#' @importFrom DBI dbConnect dbDriver dbGetQuery dbDisconnect dbListTables dbExistsTable
#' @importFrom RSQLite dbGetPreparedQuery dbIsValid
#' @importFrom mzR isInitialized openMSfile header
#' 
MsConnections <- setRefClass(
    'MsConnections',
    fields=list(
        rawFile='character',
        ramp='mzRramp',
        saryFile='character',
        saryDB='SQLiteConnection'
    ),
    methods=list(
        initialize = function(raw, sary, ...) {
            if(class(raw) != 'character' || class(sary) != 'character') {
                stop('Paths must be given as strings')
            }
            if(!file.exists(raw)) {
                stop('No file at specified location:', raw)
            }
            if(!file.exists(sary)) {
                stop('No file at specified location:', sary)
            }
            tryCatch({
                ramp <<- openMSfile(raw)
            },
            error = function(e) {
                stop('Cannot create connection to', raw, '- failed with error:\n', e)
            })
            tryCatch({
                db <- dbConnect(dbDriver('SQLite'), sary)
                dbGetQuery(db, 'pragma schema_version')
                saryDB <<- db
            },
            error = function(e) {
                stop('Cannot create connection to', sary, '- failed with error:\n', e)
            })
            rawFile <<- raw
            saryFile <<- sary
            callSuper(...)
        },
        finalize = function() {
            close(ramp)
            dbDisconnect(saryDB)
        },
        verify = function() {
            'Verify that the default tables are present and their structure conforms'
            
            presentTables <- dbListTables(sary())
            tablesConform <- all(c('header', 'scans', 'peakInfo') %in% presentTables)
            
            headerCols <- dbGetQuery(sary(), 'PRAGMA table_info(header)')$name
            scansCols <- dbGetQuery(sary(), 'PRAGMA table_info(scans)')$name
            peaksCols <- dbGetQuery(sary(), 'PRAGMA table_info(peakInfo)')$name
            
            trueHeaderCols <- sapply(strsplit(headerTableFormat, ' '), function(x) x[1])
            trueScansCols <- sapply(strsplit(scanTableFormat, ' '), function(x) x[1])
            truePeaksCols <- sapply(strsplit(peakTableFormat, ' '), function(x) x[1])
            
            headerConform <- all(trueHeaderCols %in% headerCols)
            scansConform <- all(trueScansCols %in% scansCols)
            peaksConform <- all(truePeaksCols %in% peaksCols)
            
            if(all(c(tablesConform, headerConform, scansConform, peaksConform))) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        },
        validate = function(subset=10) {
            'Checks that the header data matches between raw and sary'
            
            saryLength <- dbGetQuery(sary(), 'SELECT Count(*) FROM header')
            rawLength <- base::length(mzR())
            if(saryLength != rawLength) return(FALSE)
            
            randIndex <- sample(1:rawLength, size = subset, replace = FALSE)
            randRaw <- header(mzR(), randIndex)
            randRaw <- randRaw[order(randRaw$seqNum),]
            rownames(randRaw) <- NULL
            
            randSary <- dbGetQuery(
                sary(),
                paste0(
                    'SELECT * FROM header WHERE seqNum IN (',
                    paste(randIndex, collapse=', '),
                    ')'
                )
            )
            return(isTRUE(all.equal(randRaw, randSary)))
        },
        mzR = function() {
            'Get the mzRramp connection in a safe manner'
            
            if(!file.exists(rawFile)) {
                stop('Raw MS data file no longer present at:', rawFile)
            }
            tryCatch({
                isInitialized(ramp)
            },
            error = function(e) {
                ramp <<- openMSfile(rawFile)
            })
            
            return(ramp)
        },
        sary = function() {
            'Get the sary database connection in a safe manner'
            
            if(!file.exists(saryFile)) {
                stop('Database file no longer present at:', saryFile)
            }
            if(!dbIsValid(saryDB)) {
                saryDB <<- dbConnect(dbDriver('SQLite'), saryFile)
            }
            
            return(saryDB)
        },
        addTable = function(name, definition, rtree=FALSE) {
            'Creates the table \'name\' with the columns defined in \'definition\'.
            The definition can be a character vector that will be concatenated with \', \'.'
            
            if(dbExistsTable(sary(), name)) stop('Table already exists')
            if(rtree) {
                newTable <- paste0(
                   'CREATE VIRTUAL TABLE ',
                   name,
                   ' USING rtree(',
                   paste(definition, collapse=', '),
                   ')'
                )
            } else {
                newTable <- paste0(
                    'CREATE TABLE ', 
                    name, 
                    ' (',
                    paste(definition, collapse=', '),
                    ')'
                )
            }
            
            dbGetQuery(sary(), newTable)
        },
        addData = function(name, data) {
            'Adds data to an already created table. \'name\' specifies the name of the
            table, \'data\' is a data.frame containing the data to be added.'
            
            if(!dbExistsTable(sary(), name)) stop('Table don\'t exists')
            
            tableDef <- dbGetQuery(
                sary(), 
                paste0('PRAGMA table_info(', name,')')
            )
            nameMatch <- names(data) %in% tableDef$name
            if(any(!nameMatch)) {
                warning(sum(!nameMatch), ' Columns ignored')
            }
            
            dbGetPreparedQuery(
                sary(), 
                paste0(
                    'INSERT INTO ', 
                    name, 
                    ' (',
                    paste(names(data)[nameMatch], collapse=', '),
                    ') VALUES (',
                    paste(paste0('$', names(data)[nameMatch]), collapse=', '),
                    ')'
                ),
                bind.data=data
            )
        },
        updateData = function(name, data, key) {
            'Updates column in already existing data'
            
            if(!dbExistsTable(sary(), name)) stop('Table don\'t exists')
            
            tableDef <- dbGetQuery(
                sary(), 
                paste0('PRAGMA table_info(', name,')')
            )
            nameMatch <- names(data) %in% tableDef$name
            if(any(!nameMatch)) {
                warning(sum(!nameMatch), ' Columns ignored')
            }
            if(key %in% names(data)[!nameMatch]) {
                stop('key doesn\'t exist')
            }
            keyMatch <- names(data) %in% key
            nameMatch <- nameMatch & !keyMatch
            
            dbGetPreparedQuery(
                sary(), 
                paste0(
                    'UPDATE ', 
                    name, 
                    ' SET ',
                    paste(paste0(names(data)[nameMatch], ' = $', names(data)[nameMatch]), collapse=', '),
                    ' WHERE ',
                    paste(paste0(names(data)[keyMatch], ' = $', names(data)[keyMatch]), collapse=', ')
                ),
                bind.data=data
            )
        },
        setData = function(name, data, key) {
            'Insert or update data depending on the existence of key'
            
            keys <- dbGetQuery(
                sary(),
                paste0(
                    'SELECT ',
                    key,
                    ' FROM ',
                    name
                )
            )
            existing <- keys %in% data[, key]
            if(sum(existing) > 0) {
                updateData(name, data[existing,], key)
            }
            if(sum(!existing) > 0) {
                addData(name, data[!existing,])
            }
        },
        getHeader = function(ids, raw=FALSE) {
            header <- dbGetQuery(
                sary(), 
                paste0(
                    'SELECT * FROM ', 
                    ifelse(raw, 'header', 'currentHeader'), 
                    ' WHERE acquisitionNum IN (', 
                    paste(ids, collapse=', '), 
                    ')'
                )
            )
            header[match(header$acquisitionNum, ids),]
        },
        getScans = function(ids, raw=FALSE) {
            'Transparently extract scans, getting modified scans if they exists 
            or read from the raw file if they don\'t'
            
            if(raw) {
                seqNum <- dbGetQuery(sary(), paste0('SELECT seqNum, acquisitionNum FROM header WHERE acquisitionNum IN (', paste(ids, collapse=', '), ')'))
                seqNum <- seqNum[match(ids, seqNum$acquisitionNum), ]
                p <- mzR::peaks(mzR(), seqNum$seqNum)
            } else {
                p <- list()
                scans <- dbGetQuery(sary(), paste0('SELECT * FROM scans WHERE scanNum IN (', paste(ids, collapse=', '), ')'))
                if(sum(scans$remove == 0) > 0) {
                    p1 <- lapply(scans$scan[scans$remove == 0], function(s) {
                        matrix(unserialize(s), ncol=2)
                    })
                    p[match(scans$scanNum[scans$remove==0], ids)] <- p1
                }
                if(!all(ids %in% scans$scanNum)) {
                    nIds <- ids[!(ids %in% scans$scanNum)]
                    p2 <- getScans(nIds, raw=TRUE)
                    p[match(nIds, ids)] <- p2
                }
            }
            
            if(class(p) == 'matrix') {
                p <- list(p)
            }
            p
        },
        setScans = function(ids, scans) {
            'Set scans to new values'
            ans <- data.frame(scanNum=ids, peaksCount=NA, totIonCurrent=NA, basePeakMZ=NA, basePeakIntensity=NA, lowMZ=NA, highMZ=NA, remove=0, scan=NA)
            for(i in 1:base::length(ids)) {
                cScan <- scans[[i]]
                nPeaks <- nrow(cScan)
                ans$peaksCount[i] <- nPeaks
                ans$totIonCurrent[i] <- sum(cScan[, 2])
                bpInd <- which.max(cScan[, 2])
                ans$basePeakMZ[i] <- cScan[bpInd, 1]
                ans$basePeakIntensity[i] <- cScan[bpInd, 2]
                ans$lowMZ[i] <- cScan[1, 1]
                ans$highMZ[i] <- cScan[nPeaks, 1]
                ans$scan[i] <- list(serialize(as.vector(cScan), NULL))
            }
            addData('scans', ans)
        },
        removeScans = function(ids, value=TRUE) {
            'Remove scans from the set'
            
            setData('scans', data.frame(scanNum=ids, remove=as.integer(value)), 'scanNum')
        },
        resetScans = function(ids) {
            'Reset the state of scans back to its original value'
            if(missing(ids)) {
                dbGetQuery(sary(), 'DELETE FROM scans')
            } else {
                dbGetQuery(sary(), paste0('DELETE FROM scans WHERE scanNum IN (', paste(ids, collapse=', '), ')'))
            }
        },
        getPeaks = function(ids) {
            'Extract peaks from database'
            
            peaks <- dbGetQuery(
                sary(),
                paste0(
                    'SELECT * FROM (SELECT * FROM peakLoc WHERE peakID IN (',
                    paste(ids, collapse=', '),
                    ')) JOIN peakInfo USING (peakID)'
                )
            )
            peaks$peak <- lapply(peaks$peak, unserialize)
            peaks
        },
        extractIC = function(acqNum, mzwin, raw=FALSE) {
            scanNum <- sort(unique(do.call('c', acqNum)))
            scans <- getScans(scanNum, raw=raw)
            scanInfo <- getHeader(scanNum, raw=raw)[, c('acquisitionNum', 'retentionTime')]
            XIC <- list()
            for(i in 1:base::length(acqNum)) {
                scIndex <- match(acqNum[[i]], scanNum)
                scData <- lapply(scans[scIndex], function(x) {
                    xic <- x[mzwin[[i]](x[,1]), 2]
                    c(sum(xic), max(xic))
                })
                scData <- do.call(rbind, scData)
                scData <- cbind(scanInfo[scIndex, ], scData)
                names(scData)[3:4] <- c('totIonCurrent', 'basePeakIntensity')
                XIC[[i]] <- scData
            }
            XIC
        },
        extractIons = function(acqNum, mzwin=NULL, raw=raw) {
            scanNum <- sort(unique(do.call('c', acqNum)))
            scans <- getScans(scanNum, raw=raw)
            rt <- getHeader(scanNum, raw=raw)$retentionTime
            ions <- list()
            for(i in 1:base::length(acqNum)) {
                scIndex <- match(acqNum[[i]], scanNum)
                rtCurrent <- rep(rt[scIndex], sapply(scans[scIndex], nrow))
                ionData <- cbind(do.call(rbind, scans[scIndex]), rtCurrent)
                colnames(ionData) <- c('mz', 'intensity', 'retentionTime')
                if(!is.null(mzwin)){
                    ionData <- ionData[mzwin[[i]](ionData[, 'mz']), ]
                }
                ionData <- ionData[ionData[, 'intensity'] != 0,]
                ions[[i]] <- ionData
            }
            ions
        },
        show = function() {
            cat('An MsConnections object with connection to the following files')
            cat('\n\n')
            cat('MS raw data:  ', rawFile, '\n')
            cat('Sary database:', saryFile, '\n')
        },
        length = function() {
            as.numeric(dbGetQuery(sary(), 'SELECT Count(acquisitionNum) FROM currentHeader'))
        },
        nPeaks = function() {
            as.numeric(dbGetQuery(sary(), 'SELECT Count(OID) FROM peakInfo'))
        }
    )
)

setMethod(
    '==', c('MsConnections', 'MsConnections'),
    function(e1, e2) {
        (e1$rawFile == e2$rawFile) && (e1$saryFile == e2$saryFile)
    }
)
setMethod(
    '!=', c('MsConnections', 'MsConnections'),
    function(e1, e2) {
        !(e1 == e2)
    }
)