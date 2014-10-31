.onLoad <- function(libname, pkgname) {
    # Peak detection registration
    peakMethods$registerMethod('massifquant', massifquant, list(mode='centroid'))
}