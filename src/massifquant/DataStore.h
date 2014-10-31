#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

class DataStore {
    
    private:
        List scans;
        NumericVector scantime;
        
        double initIS2;
        double initMZS2;
        double initIS;
        
        void ghostScan(double ppm, int minScans);;
    
    public:
        DataStore(List s, NumericVector st, double ppm, int minScans);
        ~DataStore();
        
        int getTotalScanNumbers();
        double getInitIS2();
        double getInitMZS2();
        double getInitIS();
        void getScan(int scan, std::vector<double> & mzScan, std::vector<double> & intenScan);
        NumericMatrix getScan(int s);
        double getScanTime(int s);
};