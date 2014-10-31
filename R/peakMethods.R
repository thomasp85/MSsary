################################################################################
# TODO: Automatically convert minScan and maxScan to correct index on return

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
peakMethodStore <- setRefClass(
    'peakMethodStore',
    fields=list(
        functions='list',
        requirements='list',
        package='list',
        version='list'
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
        getPackage=function(name) {
            "Get the package from where a function originates"
            
            package[[name]]
        },
        getVersion=function(name) {
            "Get the version of the package from where a function originates"
            
            version[[name]]
        },
        testInput=function(name, scans, info) {
            "Test the provided inputs to see if they are sound and corresponds
             to the requirements set when the function was reqistered"
            
            if(!is.null(requirements[name]$mode)) {
                randInd <- sample(1:length(scans), 10)
                mode <- ifelse(all(isCentroided(scans[randInd])), 'centroid', 'profile')
                if(mode != requirements[name]$mode) {
                    stop('Scans are of mode ', mode, ', should be ', requirements[name]$mode)
                }
            }
            if(length(unique(info$msLevel)) != 1) {
                stop('Scans come from different MS levels')
            }
            if(!is.null(requirements[name]$msLevel) && info$msLevel[1] != requirements[name]$msLevel) {
                stop('MS level is ', info$msLevel[1], ', should be ', requirements[name]$msLevel)
            }
            if(sort(info$seqNum) != info$seqNum) {
                stop('Scans doesn\'t seem to be consecutive')
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
            
            resultNames <- c('mzMean', 'mzMin', 'mzMax', 'scanStart', 'scanEnd', 'length', 'area', 'maxHeight', 'peak')
            if(!hasMethod(name)) {
                stop('No method of name \"', name, '\" available')
            }
            testInput(name, scans, info)
            res <- getMethod(name)(scans, info, ...)
            if(class(res) != 'data.frame') {
                stop('Output should be of class data.frame')
            }
            if(all(resultNames %in% names(res))) {
                stop('Output is missing columns: ', resultNames[!resultNames %in% names(res)])
            }
            res
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
            
            functionName <- deparse(substitute(fun))
            packageName <- sub('.*:', '', getAnywhere(functionName)$where[1])
            packageVersion <- getVersion(packageName)
            
            functions[[name]] <<- fun
            requirements[[name]] <<- req
            package[[name]] <<- packageName
            version[[name]] <<- packageVersion
        }
    )
)
peakMethods <- peakMethodStore()

#' Runs the massifquant algorithm
#' 
#' Massifquant is a Kalman filter (KF) based feature detection for XC-MS data in
#' centroid mode.
#' 
#' @param scans A list of matrices with two columns. First column holds the mz
#' values for a single scan, second column holds the intensities.
#' @param info A data.frame with information for each scan. Should minimally 
#' hold a retentionTime, seqNum and msLevel column.
#' @param ppm The minimum estimated parts per million mass resolution a feature 
#' must possess.
#' @param minIntensity The minimum threshold for the maximum intensity of a 
#' feature that must be met.
#' @param minScans The minimum feature length in number of scans that a peak
#' must have.
#' @param consecMissedLim Integer: Suggested values:(1,2,3). While a feature is 
#' in the proces of being detected by a Kalman Filter, the Kalman Filter may not 
#' find a predicted centroid in every scan. After 1 or more consecutive failed 
#' predictions, this setting informs Massifquant when to stop a Kalman Filter 
#' from following a candidate feature.
#' @param criticalValue Numeric: Suggested values: (0.1-3.0). This setting helps 
#' determine the the Kalman Filter prediciton margin of error. A real centroid 
#' belonging to a bonafide feature must fall within the KF prediction margin of 
#' error. Much like in the construction of a confidence interval, criticalVal 
#' loosely translates to be a multiplier of the standard error of the prediction 
#' reported by the Kalman Filter. If the features in the XC-MS sample have a 
#' small mass deviance in ppm error, a smaller critical value might be better 
#' and vice versa.
#' @param combine Integer: set to 1 if apply t-test union on segmentation; set 
#' to 0 if no t-test to be applied on chromatographically continous features 
#' sharing same m/z range. Explanation: With very few data points, sometimes a 
#' Kalman Filter stops tracking a feature prematurely. Another Kalman Filter is 
#' instantiated and begins following the rest of the signal. Because tracking is 
#' done backwards to forwards, this algorithmic defect leaves a real feature 
#' divided into two segments or more. With this option turned on, the program 
#' identifies segmented features and combines them (merges them) into one with a 
#' two sample t-test. The potential danger of this option is that some truly 
#' distinct features may be merged.
#' @param checkBack Integer: set to 1 if turned on; set to 0 if turned off. The 
#' convergence of a Kalman Filter to a feature's precise m/z mapping is very 
#' fast, but sometimes it incorporates erroneous centroids as part of a feature 
#' (especially early on). The "scanBack" option is an attempt to remove the 
#' occasional outlier that lies beyond the converged bounds of the Kalman 
#' Filter. The option does not directly affect identification of a feature 
#' because it is a postprocessing measure; it has not shown to be a extremely 
#' useful thus far and the default is set to being turned off.
#' 
#' @return A data.frame with a row for each detected peak. The data.frame has 
#' the following columns: mz, mzmin, mzmax, scmin, scmax, length, intensity, 
#' maxint and peak. The peak column holds the sequential intensity values for
#' each peak.
#' 
#' @details This algorithm's performance has been tested rigorously on high 
#' resolution LC/{OrbiTrap, TOF}-MS data in centroid mode. Simultaneous kalman 
#' filters identify features and calculate their area under the curve. The 
#' default parameters are set to operate on a complex LC-MS Orbitrap sample. 
#' Users will find it useful to do some simple exploratory data analysis to find 
#' out where to set a minimum intensity, and identify how many scans an average 
#' feature spans. The "consecMissedLimit" parameter has yielded good performance 
#' on Orbitrap data when set to (2) and on TOF data it was found best to be at 
#' (1). This may change as the algorithm has yet to be tested on many samples. 
#' The "criticalValue" parameter is perhaps most dificult to dial in 
#' appropriately and visual inspection of peak identification is the best 
#' suggested tool for quick optimization. The "ppm" and "checkBack" parameters 
#' have shown less influence than the other parameters and exist to give users 
#' flexibility and better accuracy.
#' 
#' @author Christopher Conley and Thomas Lin Pedersen
#' 
#' @references
#' Conley, C. J., Smith, R., Torgrip, R. J. O., Taylor, R. M., Tautenhahn, R., & 
#' Prince, J. T. (2014). Massifquant: open-source Kalman filter-based XC-MS 
#' isotope trace feature detection. Bioinformatics. 
#' doi:10.1093/bioinformatics/btu359
#' 
massifquant <- function(scans, info, ppm, minIntensity, minScans, consecMissedLim, criticalVal, combine, checkBack) {
    res <- mqCpp(
        scans=scans,
        scantime=info$retentionTime,
        minIntensity=minIntensity,
        minCentroids=minScans,
        consecMissedLim=consecMissedLim,
        ppm=ppm,
        criticalVal=criticalVal,
        segs=combine,
        scanBack=checkBack
    )
    
    peakInd <- names(res) == 'peak'
    resDf <- as.data.frame(res[, -peakInd])
    peaks <- lapply(res$peak, function(x) {serialize(x, NULL)})
    resDf$peak <- peaks
    resDf
}