#include "WeightFunction.h"

double WeightFunction::WeightValue(std::array<double,3>& pos, FPRegion& region){
    //Is in FP region
    if(region.isInnerPoint(pos, region._low, region._high)){
        return 1.0;
    }
    //Is in hybrid and need to compute stuff
    if(region.isInnerPoint(pos, region._low, region._high)){

    }
    //It is in cg and w=0
    return 0.0;
}