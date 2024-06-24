#include "WeightFunction.h"

double WeightFunction::WeightValue(const std::array<double,3>& pos, FPRegion& region){
    //Is in FP region
    if(region.isInnerPoint(pos, region._low, region._high)){
        return 1.0;
    }
    //Is in hybrid and need to compute stuff
    if(region.isInnerPoint(pos, region._lowHybrid, region._highHybrid)){
        double weight=0.0;
        double at_width = region._high[0]- region._low[0];
        double hy_width = region._highHybrid[0];
        double x = std::abs(pos[0]-region._center[0]);
        x = M_PI/(2.0*hy_width) * (x - at_width);
        weight = std::cos(x) * std::cos(x);

        return weight;
    }
    //It is in cg and w=0
    return 0.0;
}

