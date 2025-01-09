#include "Statistics.h"

RegionRDFProfiler::RegionRDFProfiler(ResRegion& r):region{r}{

}


void RegionRDFProfiler::init(int bins){
    total_bins = bins;
    values.resize(bins);
    centers.resize(bins);

    bin_width = region.high[0] -  region.low[0];
    bin_width /= ((double)total_bins);

    InitCenters();

    Log::global_log->info()<<"[AdResSStatistics] Measuring local RDFs with"<<total_bins<<" bins"<<std::endl;


}

void RegionRDFProfiler::MeasureRDF(ParticleContainer* particleContainer){

    RegionParticleIterator it = particleContainer->regionIterator(region.low.data(),region.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);

    for(it;it.isValid();++it){



    }

}


/***************
 * 
 * 
 * CELL PROCESSOR FUNCTIONS
 * bunch of code repetition
 * 
 */