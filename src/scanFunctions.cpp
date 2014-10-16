#include <Rcpp.h>
using namespace Rcpp;

//' XIC from a list of scans
//' 
//' This function extracts total and base ions from each scan in a list based on
//' a min and max m/z value given
//' 
//' @param scans A list of matrices with scan data as returned by the msData 
//' method for MsScanList
//' 
//' @param mzmin A numeric giving the minimum mz value to include
//' 
//' @param mzmax A numeric giving the maximum mz value to include
//' 
//' @return A matrix with a length equal to the length of scans, and with two 
//' columns. The first column hold the summed intensities, the second the 
//' maximum intensity for each scan
//' 
//' @noRd
//' 
//[[Rcpp::export]]
NumericMatrix getXIC(List scans, double mzmin, double mzmax) {
    int nScans = scans.size();
    NumericMatrix res(nScans, 2);
    
    for(int i=0; i<nScans; i++) {
        NumericMatrix currentScan = scans[i];
        NumericMatrix::Column mz = currentScan( _, 0);
        NumericVector intensities = currentScan( _, 1);
        NumericVector XICints = intensities[mz >= mzmin & mz <= mzmax];
        res(i, 0) = sum(XICints);
        res(i, 1) = max(XICints);
    }
    return(res);
}