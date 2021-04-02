#include "MirrorSystem.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;
using Log::global_log;

MirrorSystem::MirrorSystem()
{
	_bDone = false;
}

MirrorSystem::~MirrorSystem()
{
}

void MirrorSystem::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {

}

void MirrorSystem::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
) {

	if(_bDone)
		return;

	global_log->info() << "HELLO beforeEventNewTimestep() ..." << endl;

	Domain* domain = global_simulation->getDomain();

	if(_type == MST_SHIFT) {
		const auto begin = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);

		std::array<double,3> width;
		width.at(0) = _box_new.at(0) - _box_old.at(0);
		width.at(1) = _box_new.at(1) - _box_old.at(1);
		width.at(2) = _box_new.at(2) - _box_old.at(2);

		// shift system
		for(auto it = begin; it.isValid(); ++it) {

			std::array<double,3> oldPos;
			std::array<double,3> newPos;
			oldPos.at(0) = it->r(0);
			oldPos.at(1) = it->r(1);
			oldPos.at(2) = it->r(2);

			newPos.at(0) = oldPos.at(0) + width.at(0)*0.5;
			newPos.at(1) = oldPos.at(1) + width.at(1)*0.5;
			newPos.at(2) = oldPos.at(2) + width.at(2)*0.5;

			cout << domainDecomp->getRank() << ": oldPos = " << oldPos.at(0) << "," << oldPos.at(1) << "," << oldPos.at(2);
			cout << domainDecomp->getRank() << ": newPos = " << newPos.at(0) << "," << newPos.at(1) << "," << newPos.at(2) << endl;

			it->setr(0, newPos.at(0) );
			it->setr(1, newPos.at(1) );
			it->setr(2, newPos.at(2) );
		}
//		particleContainer->update();
		global_log->info() << "System shifted." << endl;
	}
	else if(_type == MST_ENLARGE) {
		// add vapor
		std::vector<std::array<int32_t,3> > mask;
		std::array<int32_t,3> arr;
		// z=0
		// y=0
//		arr.at(0) =  0;
//		arr.at(1) =  0;
//		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  0;
		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  0;
		arr.at(2) =  0;
		mask.push_back(arr);
		// y=1
		arr.at(0) =  0;
		arr.at(1) =  1;
		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  1;
		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  1;
		arr.at(2) =  0;
		mask.push_back(arr);
		// y=-1
		arr.at(0) =  0;
		arr.at(1) = -1;
		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) = -1;
		arr.at(2) =  0;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) = -1;
		arr.at(2) =  0;

		// z=1
		// y=0
		arr.at(0) =  0;
		arr.at(1) =  0;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  0;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  0;
		arr.at(2) =  1;
		mask.push_back(arr);
		// y=1
		arr.at(0) =  0;
		arr.at(1) =  1;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  1;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  1;
		arr.at(2) =  1;
		mask.push_back(arr);
		// y=-1
		arr.at(0) =  0;
		arr.at(1) = -1;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) = -1;
		arr.at(2) =  1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) = -1;
		arr.at(2) =  1;

		// z=-1
		// y=0
		arr.at(0) =  0;
		arr.at(1) =  0;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  0;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  0;
		arr.at(2) = -1;
		mask.push_back(arr);
		// y=1
		arr.at(0) =  0;
		arr.at(1) =  1;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) =  1;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) =  1;
		arr.at(2) = -1;
		mask.push_back(arr);
		// y=-1
		arr.at(0) =  0;
		arr.at(1) = -1;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) =  1;
		arr.at(1) = -1;
		arr.at(2) = -1;
		mask.push_back(arr);
		arr.at(0) = -1;
		arr.at(1) = -1;
		arr.at(2) = -1;

		global_log->info() << "Adding new particles ..." << endl;
		uint64_t numAdded = 0;
		double bbMin[3];
		double bbMax[3];
		domainDecomp->getBoundingBoxMinMax(domain, bbMin, bbMax);
		cout << domainDecomp->getRank() << ": bbMin = " << bbMin[0] << "," << bbMin[1] << "," << bbMin[2] << endl;
		cout << domainDecomp->getRank() << ": bbMax = " << bbMax[0] << "," << bbMax[1] << "," << bbMax[2] << endl;

		for(auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			std::array<double,3> oldPos;
			oldPos.at(0) = it->r(0);
			oldPos.at(1) = it->r(1);
			oldPos.at(2) = it->r(2);

			for(auto m:mask) {
				std::array<double,3> newPos;

				newPos.at(0) = oldPos.at(0) + m.at(0)*_box_old.at(0);
				newPos.at(1) = oldPos.at(1) + m.at(1)*_box_old.at(1);
				newPos.at(2) = oldPos.at(2) + m.at(2)*_box_old.at(2);

//				if( newPos.at(0) > 0.0 && newPos.at(0) < _box_new.at(0) &&
//					newPos.at(1) > 0.0 && newPos.at(1) < _box_new.at(1) &&
//					newPos.at(2) > 0.0 && newPos.at(2) < _box_new.at(2) ) {

				if( newPos.at(0) > bbMin[0] && newPos.at(0) < bbMax[0] &&
					newPos.at(1) > bbMin[1] && newPos.at(1) < bbMax[1] &&
					newPos.at(2) > bbMin[2] && newPos.at(2) < bbMax[2] ) {

					Molecule mol;
					mol.setr(0, newPos.at(0) );
					mol.setr(1, newPos.at(1) );
					mol.setr(2, newPos.at(2) );

					cout << domainDecomp->getRank() << ": newPos = " << newPos.at(0) << "," << newPos.at(1) << "," << newPos.at(2);
					cout << domainDecomp->getRank() << ": oldPos = " << oldPos.at(0) << "," << oldPos.at(1) << "," << oldPos.at(2) << endl;

					particleContainer->addParticle(mol, true, false);
					numAdded++;
				}
			}
		}
		cout << domainDecomp->getRank() << ": Added " << numAdded << " new particles." << endl;
	}
	else if(_type == MST_MIRROR) {

	}
//	particleContainer->update();
	global_simulation->updateParticleContainerAndDecomposition(1.0, false);

	// done
	_bDone = true;
}

void MirrorSystem::readXML(XMLfileUnits& xmlconfig) {

	// type
	uint32_t type = MST_UNKNOWN;
	xmlconfig.getNodeValue("@type", type);
	_type = static_cast<MirrorSystemType>(type);

	// mirror position
	_yPos = 0.;
	xmlconfig.getNodeValue("yPos", _yPos);
	global_log->info() << "MirrorSystem: y position = " << _yPos << endl;

	// old box size
	xmlconfig.getNodeValue("box/old/x", _box_old.at(0));
	xmlconfig.getNodeValue("box/old/y", _box_old.at(1));
	xmlconfig.getNodeValue("box/old/z", _box_old.at(2));
	// old box size
	xmlconfig.getNodeValue("box/new/x", _box_new.at(0));
	xmlconfig.getNodeValue("box/new/y", _box_new.at(1));
	xmlconfig.getNodeValue("box/new/z", _box_new.at(2));
}

void MirrorSystem::beforeForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
) {
}

void MirrorSystem::afterForces(
        ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
        unsigned long simstep
)
{
}







