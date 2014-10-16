################################################################################
# TODO: Add unique, split and c methods
#       Add `[[` that returns a matrix

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
        if(!all(sapply(object@connections, class) == 'MsConnections')) {
            return('All connections must be MsConnections objects')
        }
        if(!all(c('start', 'end', 'conIndex') %in% colnames(object@mapping))) {
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
    function(object, i) {
        object@connections[[object@mapping[i, 'conIndex']]]
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
            conIndex <- x@mapping[i, 'conIndex']
            connections <- x@connections[conIndex]
            mapping <- subsetMapping(x@mapping, i)
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
        object@info
    }
)

#' @describeIn MsList Get the data as a list of matrices
#' 
setMethod(
    'msData', 'MsList',
    function(object) {
        if(length(object) == 0) return(list(NULL))
        index <- getElementIndex(object@mapping)
        res <- list()
        res[unique(index)] <- lapply(split(as.data.frame(object@data), index), as.matrix)
        res
    }
)