#' @include MsConnections.R
#' @include tableFormats.R
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
    
    connection$addTable('scans', scanTableFormat)
    
    connection$addTable('peaks', peakTableFormat)
    
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
    connection <- MsConnections(rawFile, dbFile)
    
    if(!verifySary(connection)) {
        stop(saryFile, ' does not look like a saryFile')
    }
    if(testIntegrity) {
        if(!connection$validate(testSize)) {
            stop(saryFile, ' and ', rawFile, 'does not match')
        }
    }
    new('MsData', connections=connection)
}