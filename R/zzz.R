.onLoad <- function(libname, pkgname) {
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