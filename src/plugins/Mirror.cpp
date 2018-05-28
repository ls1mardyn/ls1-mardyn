#include "Mirror.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

Mirror::Mirror()
{
}

Mirror::~Mirror()
{
}

void Mirror::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	global_log->debug() << "Mirror enabled at position: " << _yPos << std::endl;
}

void Mirror::readXML(XMLfileUnits& xmlconfig) {
	_yPos = 0.;
	xmlconfig.getNodeValue("yPos", _yPos);
	global_log->info() << "Mirror: y position = " << _yPos << endl;

	_forceConstant = 100.;
	xmlconfig.getNodeValue("forceConstant", _forceConstant);
	global_log->info() << "Mirror: force constant = " << _forceConstant << endl;

	_direction = 0.;
	xmlconfig.getNodeValue("direction", _direction);
	if(0 == _direction)
		global_log->info() << "Mirror: direction -->|" << endl;
	else if(1 == _direction)
		global_log->info() << "Mirror: direction |<--" << endl;
}

void Mirror::afterForces(
        ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
        unsigned long simstep
)
{
	this->VelocityChange(particleContainer);
}

void Mirror::VelocityChange( ParticleContainer* particleContainer) {
	double regionLowCorner[3], regionHighCorner[3];
  
	/*!*** Mirror boundary in y-direction on top of the simulation box ****/
	if(particleContainer->getBoundingBoxMax(1) > _yPos){ // if linked cell in the region of the mirror boundary
		for(unsigned d = 0; d < 3; d++){
			regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
			regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
		}

		//perform a check if the region is contained by the particleContainer???
		if(particleContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner)){
			#if defined (_OPENMP)
			#pragma omp parallel shared(regionLowCorner, regionHighCorner)
			#endif
			{
				RegionParticleIterator begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner);  //, ParticleIterator::ALL_CELLS);

				for(RegionParticleIterator it = begin; it.hasNext(); it.next()){
					double additionalForce[3];
					additionalForce[0] = 0;
					additionalForce[2] = 0;
					if( (0 == _direction && it->r(1) > _yPos) || (1 == _direction && it->r(1) < _yPos) ) {
						double distance = _yPos - it->r(1);
						additionalForce[1] = _forceConstant * distance;
//						cout << "additionalForce[1]=" << additionalForce[1] << endl;
//						cout << "before: " << (*it) << endl;
//						it->Fljcenteradd(0, additionalForce);  TODO: Can additional force be added on the LJ sites instead of adding to COM of molecule?
						it->Fadd(additionalForce);
//						cout << "after: " << (*it) << endl;
					}
				}
			}
		}
	}
}
