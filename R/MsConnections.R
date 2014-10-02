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
#' @importFrom DBI dbConnect dbDriver dbGetQuery dbCommit dbDisconnect dbListTables dbExistsTable
#' @importFrom RSQLite dbBeginTransaction dbGetPreparedQuery isIdCurrent
#' @importFrom mzR isInitialized openMSfile header
#' 
MsConnections <- setRefClass(
    'MsConnections',
    fields=list(
        rawFile='character',
        ramp='mzRramp',
        saryFile='character',
        saryDB='SQLiteConnection'),
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
            
            presentTables <- dbListTables(.self$sary())
            tablesConform <- all(c('header', 'scans', 'peaks') %in% presentTables)
            
            headerCols <- dbGetQuery(.self$sary(), 'PRAGMA table_info(header)')$name
            scansCols <- dbGetQuery(.self$sary(), 'PRAGMA table_info(scans)')$name
            peaksCols <- dbGetQuery(.self$sary(), 'PRAGMA table_info(peaks)')$name
            
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
            
            saryLength <- dbGetQuery(.self$sary(), 'SELECT Count(*) FROM header')
            rawLength <- length(.self$ramp())
            if(saryLength != rawLength) return(FALSE)
            
            randIndex <- sample(1:rawLength, size = subset, replace = FALSE)
            randRaw <- header(.self$ramp(), randIndex)
            randSary <- dbGetQuery(
                .self$sary(),
                paste0(
                    'SELECT * FROM header WHERE seqNum IN (',
                    paste(randIndex, collapse=', '),
                    ')'
                )
            )
            return(isTRUE(all.equal(randRaw, randSary)))
        },
        raw = function() {
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
            if(!isIdCurrent(saryDB)) {
                saryDB <<- dbConnect(dbDriver('SQLite'), saryFile)
            }
            
            return(saryDB)
        },
        addTable = function(name, definition) {
            'Creates the table \'name\' with the columns defined in \'definition\'.
            The definition can be a character vector that will be concatenated with \', \'.'
            
            if(dbExistsTable(.self$sary(), name)) stop('Table already exists')
            
            newTable <- paste0(
                'CREATE Table ', 
                name, 
                '(',
                paste(definition, collapse=', '),
                ')'
            )
            dbGetQuery(.self$sary(), newTable)
        },
        addData = function(name, data) {
            'Adds data to an already created table. \'name\' specifies the name of the
            table, \'data\' is a data.frame containing the data to be added.'
            
            if(!dbExistsTable(.self$sary(), name)) stop('Table don\'t exists')
            
            tableDef <- dbGetQuery(
                .self$sary(), 
                paste0('PRAGMA table_info(', name,')')
            )
            nameMatch <- names(data) %in% tableDef$name
            if(any(!nameMatch)) warning(sum(!nameMatch), ' Columns ignored')
            
            dbBeginTransaction(.self$sary())
            dbGetPreparedQuery(
                .self$sary(), 
                paste0(
                    'INSERT INTO ', 
                    name, 
                    ' VALUES (',
                    paste(paste0('$', names(data)[nameMatch]), collapse=', '),
                    ')'
                ),
                bind.data=data
            )
            dbCommit(.self$sary())
        },
        show = function() {
            cat('An MsConnections object with connection to the following files')
            cat('\n\n')
            cat('MS raw data:  ', rawFile, '\n')
            cat('Sary database:', saryFile, '\n')
        }
        )
    )