#' @export
#' 
ABOVE <- function(x) {
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'ABOVE'
    ans$data <- x
    ans$sql <- c()
    for(i in 1:length(x)) {
        if(is.na(x[i])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0("PLCHLDR > ", x[i])
        }
    }
    ans
}
#' @export
#' 
BELOW <- function(x) {
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'BELOW'
    ans$data <- x
    ans$sql <- c()
    for(i in 1:length(x)) {
        if(is.na(x[i])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0("PLCHLDR2 < ", x[i])
        }
    }
    ans
}
#' @export
#' 
EQUALS <- function(x) {
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'EQUALS'
    ans$data <- x
    ans$sql <- c()
    for(i in 1:length(x)) {
        if(is.na(x[i])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0("PLCHLDR == ", x[i])
        }
    }
    ans
}
#' @export
#' 
DIFFERS <- function(x) {
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'DIFFERS'
    ans$data <- x
    ans$sql <- c()
    for(i in 1:length(x)) {
        if(is.na(x[i])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0("PLCHLDR != ", x[i])
        }
    }
    ans
}
#' @export
#' 
IN <- function(...) {
    args <- list(...)
    if(all(sapply(args, length) == 1)) args <- list(unlist(args))
    
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'IN'
    ans$data <- args
    ans$sql <- c()
    for(i in 1:length(args)) {
        if(is.na(args[[i]][1])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0('PLCHLDR IN (', paste(args[[i]], collapse=', '), ')')
        }
    }
    ans
}
#' @export
#' 
NOTIN <- function(...) {
    args <- list(...)
    if(all(sapply(args, length) == 1)) args <- list(unlist(args))
    
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'NOTIN'
    ans$data <- args
    ans$sql <- c()
    for(i in 1:length(args)) {
        if(is.na(args[[i]][1])) {
            ans$sql[i] <- ""
        } else {
            ans$sql[i] <- paste0('PLCHLDR NOT IN (', paste(args[[i]], collapse=', '), ')')
        }
    }
    ans
}
#' @export
#' 
BETWEEN <- function(x, y) {
    if(length(x) < length(y)) {
        x <- x[rep(1:length(x), length.out=length(y))]
    }
    if(length(y) < length(x)) {
        y <- y[rep(1:length(y), length.out=length(x))]
    }
    
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'BETWEEN'
    ans$data <- list(x=x, y=y)
    ans$sql <- c()
    for(i in 1:length(x)) {
        nas <- is.na(c(x[i], y[i]))
        if(any(nas)) {
            if(!nas[1]) {
                ans$sql[i] <- paste0('PLCHLDR >= ', x[i])
            } else if(!nas[2]) {
                ans$sql[i] <- paste0('PLCHLDR2 <= ', y[i])
            } else {
                ans$sql[i] <- ""
            }
        } else {
            ans$sql[i] <- paste0("PLCHLDR >= ", x[i], ' AND PLCHLDR2 <= ', y[i])
        }
    }
    ans
}
#' @export
#' 
OVERLAPS <- function(x, y=x) {
    if(length(x) < length(y)) {
        x <- x[rep(1:length(x), length.out=length(y))]
    }
    if(length(y) < length(x)) {
        y <- y[rep(1:length(y), length.out=length(x))]
    }
    
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'OVERLAPS'
    ans$data <- list(x=x, y=y)
    ans$sql <- c()
    for(i in 1:length(x)) {
        nas <- is.na(c(x[i], y[i]))
        if(any(nas)) {
            if(!nas[1]) {
                ans$sql[i] <- paste0('PLCHLDR2 >= ', x[i])
            } else if(!nas[2]) {
                ans$sql[i] <- paste0('PLCHLDR <= ', y[i])
            } else {
                ans$sql[i] <- ""
            }
        } else {
            ans$sql[i] <- paste0("PLCHLDR2 >= ", x[i], ' AND PLCHLDR <= ', y[i])
        }
    }
    ans
}
#' @export
#' 
OUTSIDE <- function(x, y) {
    if(length(x) < length(y)) {
        x <- x[rep(1:length(x), length.out=length(y))]
    }
    if(length(y) < length(x)) {
        y <- y[rep(1:length(y), length.out=length(x))]
    }
    
    ans <- list()
    class(ans) <- 'MSsaryFilter'
    ans$type <- 'OUTSIDE'
    ans$data <- list(x=x, y=y)
    ans$sql <- c()
    for(i in 1:length(x)) {
        nas <- is.na(c(x[i], y[i]))
        if(any(nas)) {
            if(!nas[1]) {
                ans$sql[i] <- paste0('PLCHLDR <= ', x[i])
            } else if(!nas[2]) {
                ans$sql[i] <- paste0('PLCHLDR2 >= ', y[i])
            } else {
                ans$sql[i] <- ""
            }
        } else {
            ans$sql[i] <- paste0("PLCHLDR <= ", x[i], ' OR PLCHLDR2 >= ', y[i])
        }
    }
    ans
}

#' @export
#' 
BOTTOM <- function(x, what) {
    res <- list()
    class(res) <- 'MSsaryLimit'
    res$type <- 'BOTTOM'
    res$sql <- paste0(' ORDER BY ', what, ' LIMIT ', x)
    res$data <- list(x=x, what=what)
    res
}

#' @export
#' 
TOP <- function(x, what) {
    res <- list()
    class(res) <- 'MSsaryLimit'
    res$type <- 'TOP'
    res$sql <- paste0(' ORDER BY ', what, ' DESC LIMIT ', x)
    res$data <- list(x=x, what=what)
    res
}

mutateFilter <- function(filter, newdata) {
    do.call(filter$type, as.list(newdata))
}

fillSQL <- function(filter, name, name2=name) {
    filter$sql <- gsub('PLCHLDR2', name2, filter$sql)
    filter$sql <- gsub('PLCHLDR', name, filter$sql)
    filter
}
rangeFilter <- function(filter) {
    filter$type %in% c('BETWEEN', 'ABOVE', 'BELOW')
}
toRange <- function(filter) {
    if(!rangeFilter(filter)) {
        stop('Filter must construct an interval')
    }
    switch(
        filter$type,
        
        BETWEEN = cbind(filter$data$x, filter$data$y),
        ABOVE = cbind(filter$data, NA),
        BELOW = cbind(NA, filter$data)
    )
}
length.MSsaryFilter <- function(x) {
    length(x$sql)
}
`[.MSsaryFilter` <- function(x, i, j, drop) {
    i <- i%%length(x)
    if(i == 0) i <- length(x)
    args <- x$data
    if(class(args) == 'list') {
        args[] <- lapply(args, function(x) {
            x[i]
        })
    } else {
        args <- args[i]
    }
    mutateFilter(x, args)
}
print.MSsaryFilter <- function(object) {
    cat(object$type, 'filter:\n')
    cat(paste(object$sql, collapse='\n'))
}
show.MSsaryFilter <- function(object) {
    print(object)
}

toFilter <- function(x) {
    if(class(x) != 'MSsaryFilter' && class(x) != 'MSsaryLimit') {
        x <- EQUALS(x)
    }
    x
}
toFunction <- function(filter) {
    switch(
        filter$type,
        
        ABOVE = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else value > x
        }),
        BELOW = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else value < x
        }),
        EQUALS = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else value == x
        }),
        DIFFER = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else value != x
        }),
        IN = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else value %in% x
        }),
        NOTIN = lapply(filter$data, function(xin) {
            x <- xin
            function(value) if(is.na(x)) rep(TRUE, length(value)) else !(value %in% x)
        }),
        BETWEEN = mapply(function(xin, yin) {
            x <- xin
            y <- yin
            function(value, value2=value) {
                if(is.na(x) && is.na(y)) {
                    return(rep(TRUE, length(value)))
                }
                if(is.na(x)) {
                    return(value2 <= y)
                }
                if(is.na(y)) {
                    return(value >= x)
                }
                value2 <= y & value >= x
            }
        }, xin=filter$data$x, yin=filter$data$y),
        OVERLAPS = mapply(function(xin, yin) {
            x <- xin
            y <- yin
            function(value, value2=value) {
                if(is.na(x) && is.na(y)) {
                    return(rep(TRUE, length(value)))
                }
                if(is.na(x)) {
                    return(value <= y)
                }
                if(is.na(y)) {
                    return(value2 >= x)
                }
                value <= y & value2 >= x
            }
        }, xin=filter$data$x, yin=filter$data$y),
        OUTSIDE = mapply(function(xin, yin) {
            x <- xin
            y <- yin
            function(value, value2=value) {
                if(is.na(x) && is.na(y)) {
                    return(rep(TRUE, length(value)))
                }
                if(is.na(x)) {
                    return(value > y)
                }
                if(is.na(y)) {
                    return(value2 < x)
                }
                value > y | value2 < x
            }
        }, xin=filter$data$x, yin=filter$data$y)
    )
}

combineSQL <- function(..., limit, nSQL) {
    filters <- list(...)
    if(missing(nSQL)) {
        nSQL <- 0
        if(length(filters) != 0) {
            nSQL <- max(sapply(filters, function(x) {length(x$sql)}))
        }
        if(!missing(limit)) {
            nSQL <- max(nSQL, length(limit$sql))
        }
    }
    if(missing(limit)) {
        limit <- rep('', nSQL)
    } else {
        limit <- paste0(' ', rep(limit$sql, length.out=nSQL))
    }
    res <- sapply(1:nSQL, function(i) {
        query <- c()
        for(f in filters) {
            nF <- length(f$sql)
            fInd <- rep(1:nF, nSQL)
            sql <- f$sql[fInd[i]]
            if(sql != '') {
                query <- c(query, paste0('(', sql, ')'))
            }
        }
        query <- paste(query, collapse=' AND ')
        if(query != '') {
            query <- paste('WHERE', query)
        }
        paste0(query, limit[i])
    })
    res
}

rtToSeqFilter <- function(f, con) {
    f <- toFilter(f)
    rts <- f$data
    if(class(rts) == 'list') {
        rts[] <- lapply(rts, function(x) {
            seqs <- dbGetPreparedQuery(con$sary(), 'SELECT seqNum FROM currentHeader ORDER BY ABS($rt - retentionTime) LIMIT 1', bind.data=data.frame(rt=x))
            x[!is.na(x)] <- seqs[!is.na(x)]
            x
        })
    } else {
        seqs <- dbGetPreparedQuery(con$sary(), 'SELECT seqNum FROM currentHeader ORDER BY ABS($rt - retentionTime) LIMIT 1', bind.data=data.frame(rt=rts))
        rts[!is.na(rts)] <- seqs[!is.na(rts)]
    }
    mutateFilter(f, rts)
}