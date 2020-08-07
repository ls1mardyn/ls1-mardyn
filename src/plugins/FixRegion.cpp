/*
 * FixRegion.cpp
 *
 *  Created on: 9 July 2020
 *      Author: Marx
 */

#include "FixRegion.h"

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void FixRegion::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("xmin", _xMin);
	xmlconfig.getNodeValue("ymin", _yMin);
	xmlconfig.getNodeValue("zmin", _zMin);
	xmlconfig.getNodeValue("xmax", _xMax);
	xmlconfig.getNodeValue("ymax", _yMax);
	xmlconfig.getNodeValue("zmax", _zMax);
}

void FixRegion::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	for (unsigned d = 0; d < 3; d++) {
		_boxLength[d] = domain->getGlobalLength(d);
	}

	// SANITY CHECK
	if (_xMin < 0. || _yMin < 0. || _zMin < 0. || _xMax > _boxLength[0] || _yMax > _boxLength[1] ||
		_zMax > _boxLength[2]) {
		global_log->error() << "[FixRegion] INVALID INPUT!!! DISABLED!" << std::endl;
		global_log->error() << "[FixRegion] HALTING SIMULATION" << std::endl;
		// HALT SIM
		Simulation::exit(1);
		return;
	}

	global_log->info() << "[FixRegion] settings:" << std::endl;
	global_log->info() << "                  xmin: " << _xMin << std::endl;
	global_log->info() << "                  ymin: " << _yMin << std::endl;
	global_log->info() << "                  zmin: " << _zMin << std::endl;
	global_log->info() << "                  xmax: " << _xMax << std::endl;
	global_log->info() << "                  ymax: " << _yMax << std::endl;
	global_log->info() << "                  zmax: " << _zMax << std::endl;

	_molCount = 0;

	std::array<double, 3> min{_xMin, _yMin, _zMin};
	std::array<double, 3> max{_xMax, _yMax, _zMax};
	// ITERATE OVER PARTICLES
	for (auto temporaryMolecule = particleContainer->regionIterator(min.data(), max.data(), ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		 temporaryMolecule.isValid(); ++temporaryMolecule) {
		for (unsigned i = 0; i < 3; i++) {
			temporaryMolecule->setv(i, 0.0);
			temporaryMolecule->setF(i, 0.0);
		}
		_molCount++;
	}
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(_molCount);
	domainDecomp->collCommAllreduceSum();
	_molCount = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
	global_log->info() << _molCount << " molecules are inside a fixed region" << std::endl;
}

void FixRegion::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
							 unsigned long simstep) {}

void FixRegion::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
							unsigned long simstep) {
	std::array<double, 3> min{_xMin, _yMin, _zMin};
	std::array<double, 3> max{_xMax, _yMax, _zMax};
	// ITERATE OVER PARTICLES
	for (auto temporaryMolecule = particleContainer->regionIterator(min.data(), max.data(), ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		 temporaryMolecule.isValid(); ++temporaryMolecule) {
		for (unsigned i = 0; i < 3; i++) {
			temporaryMolecule->setv(i, 0.0);
		}
		temporaryMolecule->clearFM();
	}
}

void FixRegion::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
						unsigned long simstep) {}
