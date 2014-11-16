#' Class to handle MsSet connection
#' 
#' This class mainly handles continuous sqlite connection and I/O
#' 
MsSet <- setRefClass(
    'MsSet',
    fields = list(
        sarysetFile = 'character',
        sarysetDB = 'SQLiteConnection'
    ),
    methods = list(
        initialize = function(saryset, ...) {
            if(class(saryset) != 'character') {
                stop('Path must be given as strings')
            }
            if(!file.exists(saryset)) {
                stop('No file at specified location:', saryset)
            }
            tryCatch({
                db <- dbConnect(dbDriver('SQLite'), saryset)
                dbGetQuery(db, 'pragma schema_version')
                sarysetDB <<- db
            },
            error = function(e) {
                stop('Cannot create connection to', saryset, '- failed with error:\n', e)
            })
            sarysetFile <<- saryset
            callSuper(...)
        },
        finalize = function() {
            dbDisconnect(sarysetDB)
        },
        saryset = function() {
            'Get the saryset database connection in a safe manner'
            
            if(!file.exists(sarysetFile)) {
                stop('Database file no longer present at:', sarysetFile)
            }
            if(!dbIsValid(sarysetDB)) {
                sarysetDB <<- dbConnect(dbDriver('SQLite'), sarysetFile)
            }
            
            return(sarysetDB)
        },
        addTable = function(name, definition, rtree=FALSE) {
            'Creates the table \'name\' with the columns defined in \'definition\'.
            The definition can be a character vector that will be concatenated with \', \'.'
            
            if(dbExistsTable(saryset(), name)) stop('Table already exists')
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
            
            dbGetQuery(saryset(), newTable)
        },
        addData = function(name, data) {
            'Adds data to an already created table. \'name\' specifies the name of the
            table, \'data\' is a data.frame containing the data to be added.'
            
            if(!dbExistsTable(saryset(), name)) stop('Table don\'t exists')
            
            tableDef <- dbGetQuery(
                saryset(), 
                paste0('PRAGMA table_info(', name,')')
            )
            nameMatch <- names(data) %in% tableDef$name
            if(any(!nameMatch)) {
                warning(sum(!nameMatch), ' Columns ignored')
            }
            
            dbGetPreparedQuery(
                saryset(), 
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
            
            if(!dbExistsTable(saryset(), name)) stop('Table don\'t exists')
            
            tableDef <- dbGetQuery(
                saryset(), 
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
                saryset(), 
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
                saryset(),
                paste0(
                    'SELECT ',
                    key,
                    ' FROM ',
                    name
                )
            )
            existing <- keys %in% data[, key]
            updateData(name, data[existing,], key)
            addData(name, data[!existing,])
        },
        getInfo = function() {
            dbGetQuery(saryset(), 'SELECT * FROM members')
        }
    )
)