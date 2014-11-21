.onLoad <- function(libname, pkgname) {
    # Spectrum processing registration
    prepareMethods$registerMethod('localMax', localMax, list(mode='profile'))
    prepareMethods$registerMethod('ionFilter', ionFilter, list(mode='centroid'))
    # Peak detection registration
    peakMethods$registerMethod('massifquant', massifquant, list(mode='centroid'))
}

# importFrom utils findMatches
# 
.DollarNames.MsConnections <- function(x, pattern) {
    findMatches(pattern, getRefClass(class(x))$methods())
}
.DollarNames.peakMethodStore <- function(x, pattern) {
    findMatches(pattern, getRefClass(class(x))$methods())
}
.DollarNames.MsSet <- function(x, pattern) {
    findMatches(pattern, getRefClass(class(x))$methods())
}