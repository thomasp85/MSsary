################################################################################
# TODO: Add unique, split, subset
#

#' @include generics.R
#' @include aaa.R
#' @include MsData.R
#' 
NULL

#' Virtual class for MSsary list like objects
#' 
#' This class defines the storage and general methods for extracted MS data
#' within the MSsary framework. The class is virtual and must be sublassed for
#' use with e.g. scans (\code{\linkS4class{MsScanList}}). Data is stored 
#' efficiently in a single matrix, and the class itself takes care of indexing 
#' individual elements. Furthermore it retains a link to the MsConnections from
#' where the data was extracted.
#' 
#' @slot connections A list of MsConnections objects
#' 
#' @slot info A data.frame with information pertaining to each element 
#' (metadata)
#' 
#' @slot data A matrix holding the actual data
#' 
#' @slot mapping A matrix giving start and end indexes in @@data for each row in
#' @@info. Additionally it also contains the index for the MsConnections object
#' in the @@connections list
#' 
setClass(
    'MsList',
    slots = list(
        connections = 'list',
        info = 'data.frame',
        data = 'matrix',
        mapping = 'matrix'
    ),
    contains = 'VIRTUAL',
    validity = function(object) {
        if(nrow(object@info) != nrow(object@mapping)) {
            return('The mapping and info must contain one row per element')
        }
        if(!all(sapply(object@connections, class) == 'MsData')) {
            return('All connections must be MsConnections objects')
        }
        if(!all(c('start', 'end', 'conIndex', 'memberIndex', 'raw') %in% colnames(object@mapping))) {
            return('Mapping must contain the columns \'start\', \'end\' and \'conIndex\' ')
        }
        return(TRUE)
    }
)

### METHODS

#' @describeIn con Get connections from an MsList
#' 
setMethod(
    'con', 'MsList',
    function(object, i, type='MsConnection') {
        theCon <- object@connections[[object@mapping[i, 'conIndex']]]
        switch(
            type,
            
            MsConnection=con(theCon),
            MsData=theCon,
            sary=con(theCon)$sary(),
            mzR=con(theCon)$mzR()
        )
    }
)

#' @describeIn MsList Number of elements in the object
#' 
#' @param x An MsList object
#' 
setMethod(
    'length', 'MsList',
    function(x) {
        nrow(x@info)
    }
)

#' @describeIn MsList Get the names of the elements
#' 
setMethod(
    'names', 'MsList',
    function(x) {
        if(length(x@connections) == 1) {
            conNames <- ''
        } else {
            conNames <- names(x@connections)
            saryNames <- sapply(x@connections, function(x) {basename(file_path_sans_ext(con(x)$saryFile))})
            if(is.null(conNames)) {
                conNames <- saryNames
            } else {
                conNames[conNames==''] <- saryNames[conNames=='']
            }
            conNames <- paste0(conNames, ': ')
        }
        conNames[x@mapping[,'conIndex']]
    }
)

#' @describeIn MsList Subset an MsList object
#' 
setMethod(
    '[', c('MsList', 'numeric', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        i <- i[i != 0]
        if(length(i) == 0) {
            connections <- list()
            info <- x@info[0,]
            data <- x@data[0,]
            mapping <- x@mapping[0,]
        } else {
            if(all(i <= 0)) {
                i <- which(!1:length(x) %in% abs(i))
            } else if(any(i < 0)) {
                stop('only 0\'s may be mixed with negative subscripts')
            }
            info <- x@info[i, ]
            dataIndex <- mapply(function(start, end) {
                if(is.na(start)) return(integer())
                seq(start, end)
            }, x@mapping[i, 'start'], x@mapping[i, 'end'], SIMPLIFY=FALSE)
            dataIndex <- do.call(c, dataIndex)
            data <- x@data[dataIndex,]
            conIndex <- unique(x@mapping[i, 'conIndex'])
            connections <- x@connections[conIndex]
            mapping <- subsetMapping(x@mapping, i)
            mapping[, 'conIndex'] <- match(x@mapping[, 'conIndex'], conIndex)[i]
        }
        new(class(x), connections=connections, info=info, data=data, mapping=mapping)
    }
)
#' @describeIn MsList Subset an MsList object
#' 
setMethod(
    '[', c('MsList', 'logical', 'missing', 'missing'),
    function(x, i, j, ..., drop=TRUE) {
        x[which(i)]
    }
)
#' @describeIn MsList Extract raw data
#' 
setMethod(
    '[[', c('MsList', 'numeric', 'missing'),
    function(x, i, j) {
        elementInd <- x@mapping[i, c('start', 'end')]
        x@data[elementInd[1]:elementInd[2],]
    }
)
#' @describeIn MsList Combine multiple lists
#' 
setMethod(
    'c', 'MsList',
    function(x, y, ..., recursive = FALSE) {
        if(missing(y)) return(x)
        
        if(class(x) != class(y)) stop('Only MsList with the same subclass can be combined')
        if(length(y) == 0) {
            newList <- x
        } else if(length(x) == 0) {
            newList <- y
        } else {
            connections <- x@connections
            mapping <- y@mapping
            for(i in 1:length(y@connections)) {
                conElements <- mapping[, 'conIndex'] == i
                exist <- sapply(connections, function(x) {
                    x == y@connections[[i]]
                })
                if(any(exist)) {
                    mapping[conElements, 'conIndex'] <- which(exist)
                } else {
                    connections <- append(connections, y@connections[i])
                    mapping[conElements, 'conIndex'] <- length(connections)
                }
            }
            mapping[, c('start', 'end')] <- mapping[, c('start', 'end')] + nrow(x@data)
            mapping <- rbind(x@mapping, mapping)
            data <- rbind(x@data, y@data)
            rownames(data) <- NULL
            info <- rbind(x@info, y@info)
            newList <- new(class(x), connections=connections, info=info, data=data, mapping=mapping)
        }
        c(newList, ...)
    }
)

#' @describeIn MsList Sort an MsList by a property
#' 
#' @export
#' 
setMethod(
    'sort', 'MsList',
    function(x, decreasing = FALSE, by, ...) {
        if(!(by %in% names(x@info))) {
            stop(by, ' is not a property of the elements')
        }
        newOrder <- order(x@info[, by], decreasing = decreasing)
        x[newOrder]
    }
)

#' @describeIn MsList See if scans are empty
#' 
setMethod(
    'isEmpty', 'MsList',
    function(object) {
        as.logical(is.na(object@mapping[, 'start']))
    }
)
#' @describeIn MsList Removes empty elements
#' 
setMethod(
    'dropEmpty', 'MsList',
    function(object) {
        object[!isEmpty(object)]
    }
)

#' @describeIn MsList Get the info for each element as a data.frame
#' 
setMethod(
    'msInfo', 'MsList',
    function(object) {
        info <- object@info
        rownames(info) <- names(object)
        info
    }
)

#' @describeIn MsList Get the data as a list of matrices
#' 
setMethod(
    'msData', 'MsList',
    function(object) {
        if(length(object) == 0) return(list(NULL))
        index <- mapply(seq, object@mapping[, 'start'], object@mapping[, 'end'], SIMPLIFY=FALSE)
        lapply(index, function(x) object@data[x,])
    }
)