#include "Mirror.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

#include <fstream>
#include <cmath>

using namespace std;
using Log::global_log;

Mirror::Mirror()
{
	_bReflect = true;
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
	if(MD_LEFT_MIRROR == _direction)
		global_log->info() << "Mirror: direction |<--" << endl;
	else if(MD_RIGHT_MIRROR == _direction)
		global_log->info() << "Mirror: direction -->|" << endl;

	// normal distributions
	bool bRet1 = xmlconfig.getNodeValue("norm/vxz", _norm.fname.vxz);
	bool bRet2 = xmlconfig.getNodeValue("norm/vy",  _norm.fname.vy);
	_norm.enabled = bRet1 && bRet2;
	if(true == _norm.enabled) {  // TODO: move this to method: init()? Has to be called before method: afterForces(), within method Simulation::prepare_start()
		global_log->info() << "Mirror uses MB from files." << std::endl;
		this->readNormDistr();
	}
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
					double ry = it->r(1);
					double vy = it->v(1);
					if( (MD_RIGHT_MIRROR == _direction && ry > _yPos) || (MD_LEFT_MIRROR == _direction && ry < _yPos) ) {

						if(true == _bReflect) {

							if( (MD_RIGHT_MIRROR == _direction && vy > 0.) || (MD_LEFT_MIRROR == _direction && vy < 0.) ) {

								if(false == _norm.enabled)
									it->setv(1, -vy);
								else {
									double vx = _norm.vxz.front();
									_norm.vxz.pop_front();
									_norm.vxz.push_back(vx);

									/*double*/ vy = _norm.vy.front();  // double vy already declared
									_norm.vy.pop_front();
									_norm.vy.push_back(vy);

									double vz = _norm.vxz.front();
									_norm.vxz.pop_front();
									_norm.vxz.push_back(vz);

									it->setv(0, vx);
									it->setv(1, vy);
									it->setv(2, vz);
								}
							}
						}
						else {
							double distance = _yPos - ry;
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
}

void Mirror::readNormDistr()
{
	struct {
		std::ifstream vxz;
		std::ifstream vy;
	} ifs;
	ifs.vxz.open(_norm.fname.vxz, std::ios::in);
	ifs.vy.open(_norm.fname.vy, std::ios::in);

	//check to see that the file was opened correctly:
	if (!ifs.vxz.is_open() || !ifs.vy.is_open() ) {
		std::cerr << "There was a problem opening the input file!\n";
		Simulation::exit(-1);//exit or do additional error checking
	}

	double dVal = 0.0;
	//keep storing values from the text file so long as data exists:
	while (ifs.vxz >> dVal) {
		_norm.vxz.push_back(dVal);
	}
	while (ifs.vy >> dVal) {
		if(MD_LEFT_MIRROR == _direction)
			_norm.vy.push_back( abs(dVal) );
		else if (MD_RIGHT_MIRROR == _direction)
			_norm.vy.push_back( abs(dVal) * (-1.) );
		else
			Simulation::exit(-1);
	}
	// close files
	ifs.vxz.close();
	ifs.vy.close();
}







