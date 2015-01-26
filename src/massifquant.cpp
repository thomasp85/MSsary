// Massifquant

#include <string.h>

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
    const vector<double> retentionTime = as< vector<double> >(scantime);

    // Model breaks down otherwise
    if (mzq == 0) {
        mzq = 1e-6;
        mzr = sqrt(mzq);
    }
    
    // Initialize tracker manager
    TrMgr busybody(totalScanNums, minIntensity,
            minCentroids, consecMissedLim,
            ppm, criticalVal, scanBack);
    
    // Start with last scan
    dkeep.getScan(totalScanNums, mzScan, intenScan);
    busybody.setDataScan(mzScan, intenScan);
    busybody.initTrackers(iq, mzq, ir, mzr, totalScanNums);
    
    // Set up progress bar
    string progress(50, ' ');
    string header = "Finding peaks    ";
    // Begin feature finding - itereate backwards
    double progCount = 0;
    double maxScanNums = double(totalScanNums);
    double progThresh = 1;
    for (int k = totalScanNums - 1; k >= 1; k--) {
        //progress
        double perc  = (progCount/totalScanNums) * 100;
        char buffer[10];
        if (perc > progThresh) {
            progress[int(perc)/2] = '=';
            sprintf(buffer, "%3d", int(perc));
            Rcout << "\r" + header + " |" + progress + "| " + string(buffer) + "%  ";
            R_FlushConsole();
            progThresh += 1;
        }
        R_ProcessEvents();
        
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
    
    // Include segmentation correction if specified
    string footer = "";
    if (segs) {
        SegProc sproc(busybody.getPicCounts());
        sproc.groupSegments(busybody);
        sproc.splitToGroups();
        vector<int> stats = sproc.collapseGroups(busybody);
        char buffer[100];
        sprintf(buffer, "\n%d features collapsed to %d", stats[0], stats[1]);
        footer = footer + string(buffer);
    }
    
    // Create return value stores
    int nRes = busybody.getPicCounts();
    NumericVector mzRes(nRes);
    NumericVector mzminRes(nRes);
    NumericVector mzmaxRes(nRes);
    IntegerVector scstartRes(nRes);
    IntegerVector scendRes(nRes);
    IntegerVector lengthRes(nRes);
    NumericVector intensityRes(nRes);
    NumericVector maxintRes(nRes);
    NumericVector scmaxRes(nRes);
    List peak(nRes);
    NumericVector fwhmRes(nRes);
    
    // Iterate over trackers and extract features
    for (int i=0;i<busybody.getPicCounts();i++) {
        feature featInfo = busybody.iterOverFeatures(i, retentionTime);

        mzRes[i]  = featInfo.mz;
        mzminRes[i] = featInfo.mzmin;
        mzmaxRes[i] = featInfo.mzmax;
        lengthRes[i] = featInfo.length;
        scstartRes[i] = featInfo.scstart;
        scendRes[i] = featInfo.scend;
        intensityRes[i] = featInfo.intensity;
        maxintRes[i] = featInfo.maxint;
        scmaxRes[i] = featInfo.scmax;
        peak[i] = NumericVector(featInfo.ions.begin(), featInfo.ions.end());
        fwhmRes[i] = featInfo.fwhm;
    }
    Rcout << "\rDone              |" + string(50, '=') + "| 100%  " + footer << endl;
    return List::create(Named("mzMean")=mzRes, Named("mzMin")=mzminRes, Named("mzMax")=mzmaxRes, Named("length")=lengthRes, Named("FWHM")=fwhmRes, Named("scanStart")=scstartRes, Named("scanEnd")=scendRes, Named("scanMax")=scmaxRes, Named("area")=intensityRes, Named("maxHeight")=maxintRes, Named("peak")=peak);
}