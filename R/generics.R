#' Get the MsConnections object
#' 
#' This function returns the MsConnections object. Only for internal use.
#' 
#' @param object The object containing the connection
#' 
#' @return An MsConnections object
#' 
#' @keywords internal
#' 
#' @export
#' 
setGeneric(
    'con', def = function(object, ...) {standardGeneric('con')}
)

#' Get a list of all edits on an object
#' 
#' This function returns the edit history as stored in the history table in the
#' sary database.
#' 
#' @param object An object containing an MsConnections object
#' 
#' @return A data.frame with edits from first to last
#' 
#' @export
#' 
setGeneric(
    'editHistory', def = function(object) {standardGeneric('editHistory')}
)

#' Create ggplot2 ready data from object
#' 
#' This is an internal function that converts objects to long-format data.frames
#' ready for use in ggplot2 functions.
#' 
#' @param object An object to be converted to a data.frame
#' 
#' @return A data.frame with the relevant information
#' 
#' @keywords internal
#' 
setGeneric(
    'meltMS', def = function(object) {standardGeneric('meltMS')}
)

#' Generic plot function for MSsary objects
#' 
#' This function provides the general plotting interface in MSsary. All MSsary 
#' objects will have this methods, often with several different outputs 
#' depending on other parameters
#' 
#' @param object An object of an MSsary class
#' 
#' @param ... additional parameters passed along
#' 
#' @return The data.frame used to construct the ggplot
#' 
#' @export
#' 
setGeneric(
    'msPlot', def = function(object, ...) {standardGeneric('msPlot')}
)

#' Get unique names of elements
#' 
#' This function ensures unique element names even in the case of duplicated
#' elements
#' 
#' @param object A list-like object from MSsary
#' 
#' @return A character vector of the same length as object
#' 
setGeneric(
    'uNames', def = function(object) {standardGeneric('uNames')}
)

#' Test whether the elements of an object are empty
#' 
#' This function applies to list-like objects such as MsScanList and tests
#' whether the individual elements of the object are empty (i.e. doesn't contain
#' any ions). It returns a logical vector of the same length as the object. If
#' \code{TRUE} the corresponding element is empty
#' 
#' @param object A list-like object from MSsary
#' 
#' @return A logical vector of the same length as object
#' 
#' @export
#' 
setGeneric(
    'isEmpty', def = function(object) {standardGeneric('isEmpty')}
)
#' Subsets an object to remove empty elements
#' 
#' This function removes empty elements from a list-like MSsary object.
#' 
#' @param object A list-like object from MSsary
#' 
#' @return An object of the same class as the one passed to the method
#' 
#' @export
#' 
setGeneric(
    'dropEmpty', def = function(object) {standardGeneric('dropEmpty')}
)

#' Test whether the elements of an object are based on raw value
#' 
#' This function applies to list-like objects such as MsScanList and tests
#' whether the individual elements of the object are based on raw unaltered
#' values or modified ones. Note that this function just test for strictly raw-
#' ness in the sense that it tests whether \code{raw=TRUE} was set during 
#' construction. If the raw data has not been modified it is possible to have
#' raw and unraw MsLists that are numerically identical.
#' 
#' @param object A list-like object from MSsary
#' 
#' @return A logical vector of the same length as object
#' 
#' @export
#' 
setGeneric(
    'isRaw', def = function(object) {standardGeneric('isRaw')}
)

#' Extract scans from an object
#' 
#' This function allows the user to extract scans from a range of objects into
#' an MsScanList object
#' 
#' @param object The object to extract scans from
#' 
#' @param ... Aditional parameters depending on the class of @@object
#' 
#' @return An MsScanList object
#' 
#' @export
#' 
setGeneric(
    'scans', def = function(object, ...) {standardGeneric('scans')}
)

#' Extract chromatograms from an object
#' 
#' This function allows the user to extract chromatograms from a range of 
#' objects into an MsChromList object
#' 
#' @param object The object to extract chromatograms from
#' 
#' @param ... Aditional parameters depending on the class of @@object
#' 
#' @return An MsChromList object
#' 
#' @export
#' 
setGeneric(
    'chroms', def = function(object, ...) {standardGeneric('chroms')}
)

#' Extract ions from an object
#' 
#' This function allows the user to extract ions from a range of 
#' objects into an MsIonList object
#' 
#' @param object The object to extract ions from
#' 
#' @param ... Aditional parameters depending on the class of @@object
#' 
#' @return An MsIonList object
#' 
#' @export
#' 
setGeneric(
    'ions', def = function(object, ...) {standardGeneric('ions')}
)

#' Extract peaks from an object
#' 
#' This function allows the user to extract peaks from a range of 
#' objects into an MsPeakList object
#' 
#' @param object The object to extract ions from
#' 
#' @param ... Aditional parameters depending on the class of @@object
#' 
#' @return An MsPeakList object
#' 
#' @export
#' 
setGeneric(
    'peaks', def = function(object, ...) {standardGeneric('peaks')}
)

#' Get and set the mode of scans in an MsScanList
#' 
#' This function returns or sets the mode of the scans in an MsScanList, either
#' 'centroid' or 'profile'. The mode is internally calculated using a heuristic
#' that looks at the m/z difference between consecutive ions.
#' 
#' @param object An MsScanList object
#' 
#' @param value A vector with 'centroid' or 'profile' values
#' 
#' @return An MsScanList
#' 
#' @export
#' 
setGeneric(
    'scanMode', def = function(object) {standardGeneric('scanMode')}
)
#' @rdname scanMode
#' 
#' @export
#' 
setGeneric(
    'scanMode<-', def = function(object, value) {standardGeneric('scanMode<-')}
)

#' Methods for navigating scans
#' 
#' These methods facilitates creating new MsScanList objects from an already
#' defined one. They utilise the inherent relationship between scans in LC-MS/MS
#' data.
#' 
#' @param object An MsScanList object
#' 
#' @param sameLevel Should only scans of the same level as the starting scan be 
#' considered?
#' 
#' @return Either an MsScanList object or a list of these
#' 
#' @rdname MsScanList-navigators
#' @name MsScanList-navigators
#' 
#' @seealso \code{\linkS4class{MsScanList}}
#' 
NULL

#' @rdname MsScanList-navigators
#' 
#' @details
#' \code{nextScan}: Get the next consecutive scan
#' 
#' @export
#' 
setGeneric(
    'nextScan', def = function(object, sameLevel=TRUE) {standardGeneric('nextScan')}
)
#' @rdname MsScanList-navigators
#' 
#' @details
#' \code{previousScan}: Get the previous consecutive scan
#' 
#' @export
#' 
setGeneric(
    'previousScan', def = function(object, sameLevel=TRUE) {standardGeneric('previousScan')}
)
#' @rdname MsScanList-navigators
#' 
#' @details
#' \code{parent}: Get the parent scan of a fragmentation scan
#' 
#' @export
#' 
setGeneric(
    'parent', def = function(object) {standardGeneric('parent')}
)
#' @rdname MsScanList-navigators
#' 
#' @details
#' \code{children}: Get fragmentation scans from a parent scan
#' 
#' @export
#' 
setGeneric(
    'children', def = function(object) {standardGeneric('children')}
)
#' @rdname MsScanList-navigators
#' 
#' @details
#' \code{siblings}: Get all fragmentation scans from a fragmentation scans 
#' parent
#' 
#' @export
#' 
setGeneric(
    'siblings', def = function(object) {standardGeneric('siblings')}
)

#' Methods for accessing data in MsList objects
#' 
#' These methods makes it possible to extract the data stored in MsList 
#' objects.
#' 
#' @param object An MsList object
#' 
#' @rdname MsList-accessors
#' @name MsList-accessors
#' 
#' @seealso \code{\linkS4class{MsList}}
#' 
NULL

#' @rdname MsList-accessors
#' 
#' @details
#' \code{msInfo}: Get information about the elements
#' 
#' @return
#' \code{msInfo}: A data.frame
#' 
#' @export
#' 
setGeneric(
    'msInfo', def = function(object) {standardGeneric('msInfo')}
)
#' @rdname MsList-accessors
#' 
#' @details
#' \code{msData}: Get the data stored for each element
#' 
#' @return
#' \code{msInfo}: A list of matrices with the data for each element
#' 
#' @export
#' 
setGeneric(
    'msData', def = function(object) {standardGeneric('msData')}
)

#' Detect peaks in chromatographic MS data
#' 
#' @export
#' 
setGeneric(
    'detectPeaks', def = function(object, ...) {standardGeneric('detectPeaks')}
)
#' Preprocess scans
#' 
#' @export
#' 
setGeneric(
    'prepareScans', def = function(object, ...) {standardGeneric('prepareScans')}
)