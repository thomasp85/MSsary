// Massifquant

//MASSIFQUANT
#include "massifquant/OpOverload.h"
#include "massifquant/Tracker.h"
#include "massifquant/TrMgr.h"
#include "massifquant/DataStore.h"
#include "massifquant/SegProc.h"

// R
#include <Rcpp.h>

const int N_NAMES = 7;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List mqCpp(List scans, NumericVector scantime, double minIntensity, 
           double minCentroids, double consecMissedLim, double ppm, 
           double criticalVal, bool segs, bool scanBack) {

    //store data
    DataStore dkeep(scans, scantime, ppm, minCentroids);
    vector<double> mzScan;
    vector<double> intenScan;
    
    // Get model parameters
    int totalScanNums = dkeep.getTotalScanNumbers();
    double iq =  dkeep.getInitIS2();        // Intensity variance
    double mzq = dkeep.getInitMZS2();       // MZ variance
    double mzr =  sqrt(mzq);                // MZ std. dev
    double ir = dkeep.getInitIS();          // Intensity std. dev
    double * pscantime = scantime.begin();

    // Model breaks down otherwise
    if (mzq == 0) {
        mzq = 1e-6;
        mzr = sqrt(mzq);
    }

    // Show the progress please
    Rprintf("\n Detecting Kalman ROI's ... \n percent finished: ");
    R_FlushConsole();
    R_ProcessEvents();
    
    // Initialize tracker manager
    TrMgr busybody(totalScanNums, minIntensity,
            minCentroids, consecMissedLim,
            ppm, criticalVal, scanBack);
    
    // Start with last scan
    dkeep.getScan(totalScanNums, mzScan, intenScan);
    busybody.setDataScan(mzScan, intenScan);
    busybody.initTrackers(iq, mzq, ir, mzr, totalScanNums);
    
    
    // Begin feature finding - itereate backwards
    double progCount = 0;
    double maxScanNums = double(totalScanNums);
    double progThresh = 10;
    for (int k = totalScanNums - 1; k >= 1; k--) {
        //progress
        double perc  = (progCount/totalScanNums) * 100;
        if (perc > progThresh) {
            Rprintf(" %d  ", int(perc));
            R_FlushConsole();
            R_ProcessEvents();
            progThresh += 10;
        }
        busybody.setCurrScanIdx(k);
        dkeep.getScan(k, mzScan, intenScan);
        busybody.predictScan(mzScan, intenScan);
        busybody.competeAct();
        busybody.manageMissed();
        busybody.manageTracked();
        busybody.initTrackers(iq, mzq, ir, mzr, k);
        progCount++;
    }
    busybody.removeOvertimers();

    Rprintf(" %d\n", 100);
    R_FlushConsole();
    R_ProcessEvents();
    
    // Include segmentation correction if specified
    if (segs) {
        SegProc sproc(busybody.getPicCounts());
        sproc.groupSegments(busybody);
        sproc.collapseSubsets();
        sproc.solderSegs(busybody);
    }
    
    // Create return value stores
    int nRes = busybody.getPicCounts();
    NumericVector mzRes(nRes);
    NumericVector mzminRes(nRes);
    NumericVector mzmaxRes(nRes);
    IntegerVector scminRes(nRes);
    IntegerVector scmaxRes(nRes);
    IntegerVector lengthRes(nRes);
    NumericVector intensityRes(nRes);
    NumericVector maxintRes(nRes);
    List peak(nRes);
    
    // Iterate over trackers and extract features
    for (int i=0;i<busybody.getPicCounts();i++) {
        feature featInfo = busybody.iterOverFeatures(i, pscantime);

        mzRes[i]  = featInfo.mz;
        mzminRes[i] = featInfo.mzmin;
        mzmaxRes[i] = featInfo.mzmax;
        lengthRes[i] = featInfo.length;
        scminRes[i] = featInfo.scmin;
        scmaxRes[i] = featInfo.scmax;
        intensityRes[i] = featInfo.intensity;
        maxintRes[i] = featInfo.maxint;
        peak[i] = NumericVector(featInfo.ions.begin(), featInfo.ions.end());
    }
    return List::create(Named("mzMean")=mzRes, Named("mzMin")=mzminRes, Named("mzMax")=mzmaxRes, Named("length")=lengthRes, Named("scanStart")=scminRes, Named("scanEnd")=scmaxRes, Named("area")=intensityRes, Named("maxHeight")=maxintRes, Named("peak")=peak);
}