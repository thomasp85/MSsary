#include "OpOverload.h"
#include "DataStore.h"

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

DataStore::DataStore(List s, NumericVector st, double ppm, int minScans) {
    scans = s;
    scantime = st;
    ghostScan(ppm, minScans);
}

DataStore:: ~DataStore() {
    
}

int DataStore::getTotalScanNumbers() {
    return scans.size();
}

double DataStore::getInitIS2() {
    return initIS2;
}

double DataStore::getInitMZS2() {
    return initMZS2;
}

double DataStore::getInitIS() {
    return initIS;
}

void DataStore::getScan(int scan, std::vector<double> & mzScan, std::vector<double> & intenScan) {
    //pass in as reference and changes scan to scan
    mzScan.clear();
    intenScan.clear();
    
    NumericMatrix ions = getScan(scan);
    
    mzScan = as< vector<double> >(NumericVector(ions(_, 0)));
    intenScan = as< vector<double> >(NumericVector(ions(_, 1)));
}

NumericMatrix DataStore::getScan(int s) {
    return scans[s-1];
}

double DataStore::getScanTime(int s) {
    return scantime[s-1];
}

void DataStore::ghostScan(double ppm, int minScans) {
    //Find most intense ion
    int nScans = scans.size();
    NumericVector bpc(nScans);
    IntegerVector bpcInd(nScans);
    NumericVector bpcMZ(nScans);
    for(int i = 1; i <= nScans; i++) {
        NumericMatrix scan = getScan(i);
        bpcInd[i] = which_max(scan(_, 1));
        bpc[i] = scan(bpcInd[i], 1);
        bpcMZ[i] = scan(bpcInd[i], 0);
    }
    int apexScan = which_max(bpc);
    double apexVal = max(bpc);
    double mzApex = bpcMZ[apexScan];
    
    //Rprintf("apexVal is %f\n", apexVal);
    initIS = sqrt(apexVal);
    
    //Extract feature containing most intense ion
    double mzDev = (mzApex * ppm) / 1000000.0;
    int startScan = apexScan - minScans + 1;
    startScan = startScan < 1 ? 1 : startScan;
    int endScan = apexScan + minScans + 1;
    endScan = endScan > nScans ? nScans : endScan;
    
    list<double> apexFeatInt;
    list<double> apexFeatMZ;

    for(int i = startScan; i <= endScan; i++) {
        NumericMatrix scan = getScan(i);
        NumericVector scanMz = scan(_, 0);
        NumericVector scanInt = scan(_, 1);
        scanInt = scanInt[scanMz >= (mzApex-mzDev) & scanMz <= (mzApex+mzDev)];
        scanMz = scanMz[scanMz >= (mzApex-mzDev) & scanMz <= (mzApex+mzDev)];
        if(scanMz.size() > 1) {
            scanMz = scanMz[which_max(scanInt)];
            scanInt = scanInt[which_max(scanInt)];
        }
        if(scanMz.size() == 1) {
            apexFeatMZ.push_back(scanMz[0]);
            apexFeatInt.push_back(scanInt[0]);
        }
    }
    
    initMZS2 = computeAnySampVar(apexFeatMZ);
    initIS2  = computeAnySampVar(apexFeatInt);
}