#' @include methodStore.R
#' 
NULL

#' Object to handle all peak detection methods
#' 
#' This peakMethods object is the sole instance of a peakMethodStore reference 
#' class. Its responsibility is to index all registered peak detection 
#' algorithms and perform input and output check upon calling the functions. 
#' This object should not be interacted with directly but through detectPeaks 
#' methods for the different object from MSsary.
#' 
#' @field functions A named list with all registered peak detection functions
#' @field requirements A named list with requirements for all reqistered 
#' functions
#' @field package A named list with the package from where the function was 
#' registered
#' @field version A named list with the version of the package stated in the 
#' package field
#' 
prepareMethodStore <- setRefClass(
    'prepareMethodStore',
    contains='methodStore',
    methods=list(
        testInput=function(name, scans, info) {
            "Test the provided inputs to see if they are sound and corresponds
            to the requirements set when the function was reqistered"
            
            callSuper(name, scans, info)
            
            if(!is.null(requirements[name]$mode)) {
                randInd <- sample(1:length(scans), 10)
                mode <- ifelse(all(isCentroided(scans[randInd])), 'centroid', 'profile')
                if(mode != requirements[name]$mode) {
                    stop('Scans are of mode ', mode, ', should be ', requirements[name]$mode)
                }
            }
        },
        useMethod=function(name, scans, info, ...) {
            "Runs the function given by \'name\' with the supplied arguments.
            Input and output are automatically tested to ensure conformance."
            
            testInput(name, scans, info)
            res <- getMethod(name)(scans, info, ...)
            
            if(class(res) != 'list') {
                stop('Output must be in the form of a list')
            }
            if(length(res) != length(scans)) {
                stop('Output length must match input name')
            }
            sapply(res, function(x) {
                if(class(x) != 'matrix') stop('Scans must be returned in matrix form')
                if(ncol(x) != 2) stop('Scans must contain two columns')
                if(any(diff(x[,1]) < 0)) stop('Ions must be sorted by m/z')
            })
            res
        },
        registerMethod=function(name, fun, req=list()) {
            "Register a new function with the given \'name\' and optionally sets
             the requirements"
            
            if(!all(c('scans', 'info') %in% names(formals(fun)))) {
                stop('Function arguments must include \'scans\' and \'info\'')
            }
            
            callSuper(name, fun, req)
        }
    )
)
prepareMethods <- prepareMethodStore()



localMax <- function(scans, info) {
    lapply(scans, function(x) {
        x[which(diff(sign(diff(x$intensity)))==-2)+1, ]
    })
}
ionFilter <- function(scans, info, filters) {
    lapply(scans, function(x) {
        for(i in names(filters)) {
            x <- switch(
                i,
                
                threshold = x[x[, 2] > filters[[i]], ],
                top = x[sort(order(x[, 2], decreasing = TRUE)[1:filters[[i]]]), ],
                mz = x[x[, 1] >= filters[[i]][1] & x[, 1] <= filters[[i]][2], ],
                rmMZ = x[x[, 1] <= filters[[i]][1] | x[, 1] >= filters[[i]][2], ],
                DNL = x[x[, 2] >= scanNoise(x, filters[[i]][1], filters[[i]][2]), ]
            )
        }
        x
    })
}