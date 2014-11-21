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

//' Get intensity threshold for a given signal-to-noise ratio
//' 
//' This function takes a scan and finds the first non-noise ion in it based on
//' the dynamic noise level algorithm described by Xu and Freitas (2009). It 
//' returns the intensity of that ion.
//' 
//' @param scan A matrix with mz and intensity values
//' @param sn The required minimum signal-to-noise to be considered a real 
//' signal
//' @param rho The modifier to use for the second lowest ion special case
//' 
//' @return A numeric with the intensity of the first ion that gets accepted as
//' a true signal
//' 
//' @references Xu, H., & Freitas, M. A. (2009). A dynamic noise level algorithm 
//' for spectral screening of peptide MS/MS spectra. BMC Bioinformatics, 11, 
//' 436–436. doi:10.1186/1471-2105-11-436
//' 
//[[Rcpp::export]]
NumericVector scanNoise(NumericMatrix scan, double sn, double rho) {
    NumericVector intensity = scan(_,1);
    std::sort(intensity.begin(), intensity.end());
    double iHat;
    double SNR;
    NumericVector res(1);
    
    int i = 0;
    double xSum = 0;
    double ySum = 0;
    double xxSum = 0;
    double xySum = 0;
    double slope;
    double intercept;
    
    for(NumericVector::iterator it = intensity.begin(); it != intensity.end(); ++it) {
        if(i == 1) {
            iHat = *(it-1) * (1 + rho);
        } else if(i > 1) {
            slope = ((i+1) * xySum - xSum * ySum) / ((i+1) * xxSum - xSum * xSum);
            intercept = (ySum - slope * xSum) / i;
            iHat = slope * (i+1) + intercept;
        }
        if(i != 0) {
            SNR = (*it) / iHat;
            if(SNR > sn) {
                res[0] = *it;
                break;
            }
        }
        
        xSum += i;
        ySum += *it;
        xxSum += i*i;
        xySum += i*(*it);
        i++;
    }
    return res;
}