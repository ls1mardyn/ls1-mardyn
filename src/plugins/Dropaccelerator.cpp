/*
 * Dropaccelerator.h
 *
 *  Created on: 26 June 2020
 *      Author: Koch
 */

#include "Dropaccelerator.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void Dropaccelerator::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("xposition", _xPosition);
	xmlconfig.getNodeValue("yposition", _yPosition);
	xmlconfig.getNodeValue("zposition", _zPosition);
	xmlconfig.getNodeValue("dropradius", _dropRadius);
	xmlconfig.getNodeValue("velocity", _veloc);
	xmlconfig.getNodeValue("starttime", _startSimStep);
	xmlconfig.getNodeValue("steps", _steps);

	// SANITY CHECK
	if (_interval < 1 || _steps <= 0 || _startSimStep < 0 || _xPosition <= 0. || _yPosition <= 0. || _zPosition <= 0. ||
		_dropRadius <= 0) {
		Log::global_log->error() << "[Dropaccelerator] INVALID CONFIGURATION!!! DISABLED!" << std::endl;
		Log::global_log->error() << "[Dropaccelerator] HALTING SIMULATION" << std::endl;
		_enabled = false;
		// HALT SIM
		Simulation::exit(1);
		return;
	}

	Log::global_log->info() << "[Dropaccelerator] settings:" << std::endl;
	Log::global_log->info() << "                  xposition: " << _xPosition << std::endl;
	Log::global_log->info() << "                  yposition: " << _yPosition << std::endl;
	Log::global_log->info() << "                  zposition: " << _zPosition << std::endl;
	Log::global_log->info() << "                  dropradius: " << _dropRadius << std::endl;
	Log::global_log->info() << "                  velocity: " << _veloc << std::endl;
	Log::global_log->info() << "                  starttime: " << _startSimStep << std::endl;
	Log::global_log->info() << "                  steps: " << _steps << std::endl;
}

void Dropaccelerator::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
								  unsigned long simstep) {
	double corrVeloc = _veloc / _steps;

	if (_enabled) {
		Log::global_log->debug() << "[Dropaccelerator] after forces called" << std::endl;

		if ((simstep - 1) % _interval != 0) {
			return;
		}

		int particlesInDrop = 0;

		// ITERATE OVER PARTICLES AND CHOOSE AND MARK IF MOLECULE IN DROPLET OR NOT
		if (simstep == _startSimStep) {
			// resize first
			_particleIsInDroplet.clear();
			_particleIsInDroplet.resize(global_simulation->getDomain()->getglobalNumMolecules(true, particleContainer, domainDecomp), false);

			for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				 temporaryMolecule.isValid(); ++temporaryMolecule) {
				double distanceSquared[3];

				_particleIsInDroplet[temporaryMolecule->getID()] = false;
				distanceSquared[0] = (_xPosition - temporaryMolecule->r(0)) * (_xPosition - temporaryMolecule->r(0));
				distanceSquared[1] = (_yPosition - temporaryMolecule->r(1)) * (_yPosition - temporaryMolecule->r(1));
				distanceSquared[2] = (_zPosition - temporaryMolecule->r(2)) * (_zPosition - temporaryMolecule->r(2));

				if (distanceSquared[0] + distanceSquared[1] + distanceSquared[2] < _dropRadius * _dropRadius) {
					_particleIsInDroplet[temporaryMolecule->getID()] = true;
					temporaryMolecule->vadd(0, corrVeloc, 0);
					particlesInDrop++;
				}
			}

#ifdef ENABLE_MPI
			MPI_Allreduce(MPI_IN_PLACE, _particleIsInDroplet.data(), _particleIsInDroplet.size(), MPI_CHAR, MPI_LOR,
						  domainDecomp->getCommunicator());
#endif

#ifdef ENABLE_PERSISTENT
			auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), particlesInDrop);
			collComm.persistent();
			collComm.get(particlesInDrop);
#else
			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(particlesInDrop);
			domainDecomp->collCommAllreduceSum();
			particlesInDrop = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();
#endif
		}

		// ITERATE OVER PARTICLES AND ACCELERATE ONLY DROPMOLECULES

		if (simstep > _startSimStep && simstep < (_startSimStep + _steps)) {
			for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				 temporaryMolecule.isValid(); ++temporaryMolecule) {
				if (_particleIsInDroplet[temporaryMolecule->getID()]) {
					temporaryMolecule->vadd(0, corrVeloc, 0);
					particlesInDrop++;
				}
			}

#ifdef ENABLE_PERSISTENT
			auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), particlesInDrop);
			collComm.persistent();
			collComm.get(particlesInDrop);
#else
			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(particlesInDrop);
			domainDecomp->collCommAllreduceSum();
			particlesInDrop = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();
#endif
		}

		// CHECK IF VELOCITY HAS REACHED AND IF NOT ACCELERATE AGAIN

		if (simstep >= (_startSimStep + _steps) && simstep < (_startSimStep + 2 * _steps)) {
			_velocNow = 0;
			particlesInDrop = 0;
			double deltavelocity = 0;

			// GET VELOCITY IN Y-DIRECTION FOR EVERY PARTICLE IN DROPLET

			for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				 temporaryMolecule.isValid(); ++temporaryMolecule) {
				double partVelocity = temporaryMolecule->v(1);

				if (_particleIsInDroplet[temporaryMolecule->getID()]) {
					particlesInDrop += 1;
					_velocNow += partVelocity;
				}
			}

			// COMMUNICATION
#ifdef ENABLE_PERSISTENT
			auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), _velocNow, particlesInDrop);
			collComm.persistent();
			collComm.get(_velocNow, particlesInDrop);
#else
			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendDouble(_velocNow);
			domainDecomp->collCommAllreduceSum();
			_velocNow = domainDecomp->collCommGetDouble();
			domainDecomp->collCommFinalize();

			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(particlesInDrop);
			domainDecomp->collCommAllreduceSum();
			particlesInDrop = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();
#endif
			// CALCULATE AVERAGE SPEED
			_velocNow = _velocNow / particlesInDrop;
			deltavelocity = _veloc - _velocNow;

			// ACCELERATE IF TOO SLOW

			if ((abs(_velocNow) < abs(_veloc)) && (abs(deltavelocity) > abs(corrVeloc))) {
				for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
					 temporaryMolecule.isValid(); ++temporaryMolecule) {
					if (_particleIsInDroplet[temporaryMolecule->getID()]) {
						temporaryMolecule->vadd(0, corrVeloc, 0);
					}
				}
			} else if ((abs(_velocNow) < abs(_veloc)) && (abs(deltavelocity) < abs(corrVeloc))) {
				for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
					 temporaryMolecule.isValid(); ++temporaryMolecule) {
					if (_particleIsInDroplet[temporaryMolecule->getID()]) {
						temporaryMolecule->vadd(0, deltavelocity, 0);
					}
				}
			}
		}
	}
}

void Dropaccelerator::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
							  unsigned long simstep) {}
