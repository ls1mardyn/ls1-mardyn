#include "ProfilerPMF.h"

InternalCellProcessor::InternalCellProcessor(const double cr, InternalProfiler* r):CellProcessor{cr,cr},my_profiler{r}{

}

double InternalCellProcessor::DistanceBetweenCOMs(std::array<double,3>& c1, std::array<double,3>& c2){
    double r =0.0;
    std::array<double,3> diff={0.0,0.0,0.0};
    
    for(int i=0;i<diff.size();i++){
        diff[i]=c1[i]-c2[i];
    }

    r = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
    return r;

}

double InternalCellProcessor::PotentialCallBack(double e, double s, double r2){
    double pot=0;

    double inv_r2 = 1.0/r2;
    double lj6 = s*s*inv_r2;
    lj6 = lj6*lj6*lj6;
    double lj12 = lj6*lj6;
    double lj12m6 = lj12 - lj6;
    pot = lj12m6*e*4.0;

    return pot;
}