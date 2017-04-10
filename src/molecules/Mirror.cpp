#include "Mirror.h"

using namespace std;
using Log::global_log;

Mirror::Mirror(){}

Mirror::~Mirror(){}

void Mirror::initialize(const std::vector<Component>* /*components*/, double in_yMirr, double in_forceConstant) {
	global_log->info() << "Initializing the mirror function.\n";
	this->_yMirr = in_yMirr;
	this->_forceConstant = in_forceConstant;
}

void Mirror::VelocityChange( ParticleContainer* partContainer, Domain* /*domain*/) {
	double regionLowCorner[3], regionHighCorner[3];
  
	/*!*** Mirror boundary in y-direction on top of the simulation box ****/
	if(partContainer->getBoundingBoxMax(1) > _yMirr){ // if linked cell in the region of the mirror boundary
		for(unsigned d = 0; d < 3; d++){
			regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
			regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
		}

		//perform a check if the region is contained by the particleContainer???
		if(partContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner)){
			#if defined (_OPENMP)
			#pragma omp parallel shared(regionLowCorner, regionHighCorner)
			#endif
			{
				RegionParticleIterator begin = partContainer->iterateRegionBegin(regionLowCorner, regionHighCorner);
				RegionParticleIterator end = partContainer->iterateRegionEnd();

				for(RegionParticleIterator i = begin; i != end; ++i){
					double additionalForce[3];
					additionalForce[0] = 0;
					additionalForce[2] = 0;
					if ((*i).r(1) > _yMirr){
						double distance = (*i).r(1) - _yMirr;
						additionalForce[1] = -_forceConstant * distance;
						(*i).Fljcenteradd(0, additionalForce);
					}
				}
			}
		}
	}
}
