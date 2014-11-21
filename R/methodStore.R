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
setRefClass(
    'methodStore',
    fields=list(
        functions='list',
        requirements='list',
        package='list',
        packVersion='list'
    ),
    methods=list(
        initialize=function() {
            "Initialize the object to be empty"
            
            functions <<- list()
            requirements <<- list()
        },
        hasMethod=function(name) {
            "Check for registered functions of the given name"
            
            any(names(functions) == name)
        },
        getMethod=function(name) {
            "Get the function of the given name"
            
            functions[[name]]
        },
        methodArgs=function(name) {
            names(formals(getMethod(name)))
        },
        getPackage=function(name) {
            "Get the package from where a function originates"
            
            package[[name]]
        },
        getVersion=function(name) {
            "Get the version of the package from where a function originates"
            
            packVersion[[name]]
        },
        testInput=function(name, scans, info) {
            "Test the provided inputs to see if they are sound and corresponds
            to the requirements set when the function was reqistered"
            
            if(!hasMethod(name)) {
                stop('No method of name \"', name, '\" available')
            }
            if(!is.null(requirements[name]$extraInfo)) {
                extraInfo <- requirements[name]$extraInfo %in% names(info)
                if(!all(extraInfo)) {
                    stop('Missing scan information: ', paste(requirements[name]$extraInfo[!extraInfo], collapse=', '))
                }
            }
        },
        useMethod=function(name, scans, info, ...) {
            "Runs the function given by \'name\' with the supplied arguments.
            Input and output are automatically tested to ensure conformance."
            
            stop('useMethod must be overridden for methodStore subclasses')
        },
        registerMethod=function(name, fun, req=list()) {
            "Register a new function with the given \'name\' and optionally sets
            the requirements"
            
            if(any(names(functions) == name)) {
                stop('Name already exists')
            }
            if(class(requirements) != 'list') {
                stop('Requirements must be in the form of a list')
            }
            
            functionName <- sub('\\(\\)$', '', deparse(match.call(call=sys.call(-1))['fun']))
            packageName <- sub('.*:', '', getAnywhere(functionName)$where[1])
            
            functions[[name]] <<- fun
            requirements[[name]] <<- req
            package[[name]] <<- packageName
            packVersion[[name]] <<- as.character(packageVersion(packageName))
        }
    )
)