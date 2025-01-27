/*
 * COMaligner.h
 *
 *  Created on: 22 June 2020
 *      Author: Marx
 */

#include "Dropaligner.h"

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void Dropaligner::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("xpos", _xPos);
	xmlconfig.getNodeValue("ypos", _yPos);
	xmlconfig.getNodeValue("zpos", _zPos);
	xmlconfig.getNodeValue("radius", _radius);
	xmlconfig.getNodeValue("interval", _interval);
	xmlconfig.getNodeValue("correctionFactor", _alignmentCorrection);

	// SANITY CHECK
	if (_interval < 1 || _alignmentCorrection < 0 || _alignmentCorrection > 1 || _xPos <= 0. || _yPos <= 0. ||
		_zPos <= 0. || _radius <= 0) {
		Log::global_log->error() << "[Dropaligner] INVALID CONFIGURATION!!! DISABLED!" << std::endl;
		Log::global_log->error() << "[Dropaligner] HALTING SIMULATION" << std::endl;
		_enabled = false;
		// HALT SIM
		Simulation::exit(1);
		return;
	}

	Log::global_log->info() << "[Dropaligner] settings:" << std::endl;
	Log::global_log->info() << "                  xpos: " << _xPos << std::endl;
	Log::global_log->info() << "                  ypos: " << _yPos << std::endl;
	Log::global_log->info() << "                  zpos: " << _zPos << std::endl;
	Log::global_log->info() << "                  radius: " << _radius << std::endl;
	Log::global_log->info() << "                  interval: " << _interval << std::endl;
	Log::global_log->info() << "                  correctionFactor: " << _alignmentCorrection << std::endl;
}

void Dropaligner::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
							   unsigned long simstep) {
	if (_enabled) {
		Log::global_log->debug() << "[Dropaligner] before forces called" << std::endl;

		if ((simstep - 1) % _interval != 0) {
			return;
		}

		// RESET
		for (unsigned d = 0; d < 3; d++) {
			_balance[d] = 0.0;
			_motion[d] = 0.0;
		}
		_mass = 0.;

		// ITERATE OVER PARTICLES
		for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 temporaryMolecule.isValid(); ++temporaryMolecule) {
			double partMass = temporaryMolecule->mass();
			double distanceSquared[3];

			distanceSquared[0] = (_xPos - temporaryMolecule->r(0)) * (_xPos - temporaryMolecule->r(0));
			distanceSquared[1] = (_yPos - temporaryMolecule->r(1)) * (_yPos - temporaryMolecule->r(1));
			distanceSquared[2] = (_zPos - temporaryMolecule->r(2)) * (_zPos - temporaryMolecule->r(2));

			if (distanceSquared[0] + distanceSquared[1] + distanceSquared[2] < _radius * _radius) {
				_mass += partMass;
				for (int d = 0; d < 3; d++) {
					_balance[d] += temporaryMolecule->r(d) * partMass;
				}
			}
		}

		// COMMUNICATION
		auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), _balance[0], _balance[1], _balance[2], _mass);
		collComm.communicate();
		collComm.get(_balance[0], _balance[1], _balance[2], _mass);

		// CALCULATE MOTION

		_motion[0] = -_alignmentCorrection * ((_balance[0] / _mass) - _xPos);
		_motion[1] = -_alignmentCorrection * ((_balance[1] / _mass) - _yPos);
		_motion[2] = -_alignmentCorrection * ((_balance[2] / _mass) - _zPos);

		// MOVE
		for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 temporaryMolecule.isValid(); ++temporaryMolecule) {
			double distanceSquared[3];

			distanceSquared[0] = (_xPos - temporaryMolecule->r(0)) * (_xPos - temporaryMolecule->r(0));
			distanceSquared[1] = (_yPos - temporaryMolecule->r(1)) * (_yPos - temporaryMolecule->r(1));
			distanceSquared[2] = (_zPos - temporaryMolecule->r(2)) * (_zPos - temporaryMolecule->r(2));

			if (distanceSquared[0] + distanceSquared[1] + distanceSquared[2] < _radius * _radius) {
				for (unsigned d = 0; d < 3; d++) {
					temporaryMolecule->move(d, _motion[d]);
				}
			}
		}
	}
}

void Dropaligner::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
						  unsigned long simstep) {}
