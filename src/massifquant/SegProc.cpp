//SegProc.cpp

#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>

#include "OpOverload.h"
#include "TrMgr.h"
#include "SegProc.h"

#include <R.h>
#include <Rdefines.h>
#include <Rcpp.h>


using namespace Rcpp;

using namespace std;
SegProc::SegProc(int otn) : origTrNum(otn) {

    edges = vector<int> (origTrNum, -1);
}

SegProc::~SegProc() {
}

void SegProc::groupSegments(TrMgr & busybody) {
    vector<int> picIdx = busybody.getPicIdx();
    list<int> candIdx;
    int ppm = busybody.getPpm();

    vector<int>::iterator it_i;
    
    // Set up progress
    string progress(50, ' ');
    string header = "Grouping peaks   ";
    double progThresh = 1;
    double nPic = double(picIdx.size());
    int i  = -1; //use for segCluster indexing
    
    for (it_i = picIdx.begin(); it_i != picIdx.end(); ++it_i) {
        ++i;
        
        double perc  = (i/nPic) * 100;
        char buffer[10];
        if (perc > progThresh) {
            progress[int(perc)/2] = '=';
            sprintf(buffer, "%3d", int(perc));
            Rcout << "\r" + header + " |" + progress + "| " + string(buffer) + "%  ";
            R_FlushConsole();
            R_ProcessEvents();
            progThresh += 1;
        }
        
        candIdx.clear(); //candidates are different for each tracker
        
        //mz mean and position check
        double mzTol =  busybody.getTracker(*it_i)->getXbar() *  ppm / 1e6;
        std::vector<int>::iterator it_j;
        int j = -1;
        for (it_j = picIdx.begin(); it_j != picIdx.end(); ++it_j) {
            ++j;
            if (*it_i == *it_j) {
                continue;
            }
            double jmeandiff = fabs(busybody.getTracker(*it_i)->getXbar() - busybody.getTracker(*it_j)->getXbar());
            if (jmeandiff > mzTol) {
                continue;
            }
            bool before = busybody.getTracker(*it_j)->getStopScanIdx() < busybody.getTracker(*it_i)->getStartScanIdx();
            int gap = busybody.getTracker(*it_i)->getStartScanIdx() - busybody.getTracker(*it_j)->getStopScanIdx();
            if (before && gap <= MAXGAP) {
                candIdx.push_back(j);
            }
            
        }
        if (candIdx.size() == 0) {
            continue;
        }
        
        // Find best candidate
        int bestMatch = compareCandidates(busybody, i, candIdx);
        edges[i] = bestMatch;
    }
}

void SegProc::splitToGroups() {
    
    vector<int> edgesInverse(edges.size(), -1);
    
    vector<int>::iterator it;
    int i = 0;
    for(it = edges.begin(); it != edges.end(); ++it) {
        if(*it != -1) {
            edgesInverse[*it] = i;
        }
        ++i;
    }
    
    i = 0;
    for(it = edges.begin(); it != edges.end(); ++it) {
        if(*it == -1) {
            int nextPic = edgesInverse[i];
            if(nextPic != -1) {
                vector<int> group;
                group.push_back(i);
                while(nextPic != -1) {
                    group.push_back(nextPic);
                    nextPic = edgesInverse[nextPic];
                }
                edgeGroups.push_back(group);
            }
        }
        ++i;
    }
}

int SegProc::compareCandidates(TrMgr & busybody, const int seed, const std::list<int> candidates) {
    vector<int> picIdx = busybody.getPicIdx();
    int seedInd = picIdx[seed];
    
    int res = -1;
    double bestP = 1;
    int bestInd = -1;
    
    std::list<int>::const_iterator it;
    for (it = candidates.begin(); it != candidates.end(); ++it) {
        int candInd = picIdx[*it];
        
        double varRatio = busybody.getTracker(seedInd)->getS2() / busybody.getTracker(candInd)->getS2();
        if (varRatio < TROBUST1 ||
            varRatio > TROBUST2) {

            ttestEq(busybody.getTracker(seedInd)->getXbar(),
                    busybody.getTracker(candInd)->getXbar(),
                    busybody.getTracker(seedInd)->getTrLen(),
                    busybody.getTracker(candInd)->getTrLen(),
                    busybody.getTracker(seedInd)->getS2(),
                    busybody.getTracker(candInd)->getS2());

        }
        else {
            ttestWelch(busybody.getTracker(seedInd)->getXbar(),
                    busybody.getTracker(candInd)->getXbar(),
                    busybody.getTracker(seedInd)->getTrLen(),
                    busybody.getTracker(candInd)->getTrLen(),
                    busybody.getTracker(seedInd)->getS2(),
                    busybody.getTracker(candInd)->getS2());
        }

        p  = 2*ptest(fabs(t), v, 0, 0);
        
        if(p < bestP) {
            bestP = p;
            bestInd = *it;
        }
    }
    if(bestP < ALPHA) {
        res = bestInd;
    }
    return res;
}

vector<int> SegProc::collapseGroups(TrMgr & busybody) {
    vector<int> picIdx = busybody.getPicIdx();
    vector<int> erasedPicIdx;
    
    int nCollapse = 0;
    
    // Set up progress
    string progress(50, ' ');
    string header = "Collapsing groups";
    double progThresh = 1;
    double nGroups = double(edgeGroups.size());
    int i = 0;
    
    vector< vector<int> >::iterator it_i;
    for(it_i = edgeGroups.begin(); it_i != edgeGroups.end(); ++it_i) {
        
        double perc  = (i/nGroups) * 100;
        char buffer[10];
        if (perc > progThresh) {
            progress[int(perc)/2] = '=';
            sprintf(buffer, "%3d", int(perc));
            Rcout << "\r" + header + " |" + progress + "| " + string(buffer) + "%  ";
            R_FlushConsole();
            R_ProcessEvents();
            progThresh += 1;
        }
        
        nCollapse += it_i->size();
        int master = picIdx[it_i->back()];
        vector<int>::reverse_iterator it_j;
        for(it_j = ++(it_i->rbegin()); it_j != it_i->rend(); ++it_j) {
            int slave = picIdx[*it_j];
            
            list<int> sl = busybody.getTracker(slave)->getScanList();
            list<int> cl = busybody.getTracker(slave)->getCentroidList();
            list<double> ml = busybody.getTracker(slave)->getMzList();
            list<double> il = busybody.getTracker(slave)->getIntensityList();
            busybody.getTracker(master)->appendToTracker(sl, cl, ml, il);
            
            erasedPicIdx.push_back(slave);
        }
        ++i;
    }
    busybody.erasePicElements(erasedPicIdx);
    
    vector<int> stats;
    stats.push_back(nCollapse);
    stats.push_back(edgeGroups.size());
    
    return stats;
}

void SegProc::ttestEq(double xbar1, double xbar2, double n1, double n2, double s12, double s22) {
    v = n1 + n2 - 2;
    double sp2 = ((n1 - 1)*s12 + (n2 - 1)*s22)/v;
    t = (xbar1 - xbar2)/sqrt(sp2*(1/n1 + 1/n2));
}

void SegProc::ttestWelch(double xbar1, double xbar2, double n1, double n2, double s12, double s22) {

    t = (xbar1 - xbar2) / sqrt(s12/n1 + s22/n2);

    double vnumer = (s12 / n1 + s22 / n2 ) * (s12 / n1 + s22 / n2 ); //square the numerator
    double vdenom = (s12 * s12) / (n1*n1 * (n1 - 1)) +
        (s22 * s22) / (n2*n2 * (n2 - 1));

    v = vnumer/vdenom;
}

double SegProc::ptest(double x, double n, int lower_tail, int log_p) {
    /* return  P[ T <= x ]  where
     *  * T ~ t_{n}  (t distrib. with n degrees of freedom).
     *
     *   *  --> ./pnt.c for NON-central
     *    */
    double val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
#endif
    if (n <= 0.0) ML_ERR_return_NAN;

    if(!R_FINITE(x))
        return (x < 0) ? R_DT_0 : R_DT_1;
    if(!R_FINITE(n))
        return R::pnorm(x, 0.0, 1.0, lower_tail, log_p);

#ifdef R_version_le_260
    if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
        /* Approx. from  Abramowitz & Stegun
         * 26.7.8 (p.949) */
        val = 1./(4.*n);
        return pnorm(x*(1. - val)/sqrt(1. + x*x*2.*val), 0.0, 1.0,
                lower_tail, log_p);
    }
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing
     * rather than gaining precision,
     *      * now that pbeta(*, log_p = TRUE)
     *      is much better.
     *           * Note however that a version
     *           of this test *is* needed for
     *           x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
        /* Danger of underflow. So use
         * Abramowitz & Stegun 26.5.4
         *     pbeta(z, a, b) ~ z^a(1-z)^b
         *     / aB(a,b) ~ z^a / aB(a,b),
         *         with z = 1/nx,  a =
         *         n/2,  b= 1/2 :
         *          */
        double lval;
        lval = -0.5*n*(2*log(fabs(x)) - log(n))
            - R::lbeta(0.5*n, 0.5) - log(0.5*n);
        val = log_p ? lval : exp(lval);
    } else {
        val = (n > x * x)
            ? R::pbeta (x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
            : R::pbeta (1. / nx,             n / 2., 0.5, /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if  lower_tail  and  x
     * > 0 (but not both):*/
    if(x <= 0.)
    lower_tail = !lower_tail;

    if(log_p) {
        if(lower_tail) return log1p(-0.5*exp(val));
        else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
        val /= 2.;
        return R_D_Cval(val);
    }
}
