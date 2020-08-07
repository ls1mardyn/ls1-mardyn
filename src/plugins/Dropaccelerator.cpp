/*
 * Dropaccelerator.h
 *
 *  Created on: 26 June 2020
 *      Author: Koch
 */

#include "Dropaccelerator.h"

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void Dropaccelerator::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("xposition", _xPosition);
	xmlconfig.getNodeValue("yposition", _yPosition);
	xmlconfig.getNodeValue("zposition", _zPosition);
	xmlconfig.getNodeValue("dropradius", _Dropradius);
	xmlconfig.getNodeValue("velocity", _Veloc);
	xmlconfig.getNodeValue("starttime", _StartSimstep);
	xmlconfig.getNodeValue("steps", _Steps);

	// SANITY CHECK
	if (_INTERVAL < 1 || _Steps <= 0 || _Veloc < 0 || _StartSimstep < 0 || _xPosition <= 0. || _yPosition <= 0. ||
		_zPosition <= 0. || _Dropradius <= 0) {
		global_log->error() << "[Dropaccelerator] INVALID CONFIGURATION!!! DISABLED!" << std::endl;
		global_log->error() << "[Dropaccelerator] HALTING SIMULATION" << std::endl;
		_enabled = false;
		// HALT SIM
		Simulation::exit(1);
		return;
	}

	global_log->info() << "[Dropaccelerator] settings:" << std::endl;
	global_log->info() << "                  xposition: " << _xPosition << std::endl;
	global_log->info() << "                  yposition: " << _yPosition << std::endl;
	global_log->info() << "                  zposition: " << _zPosition << std::endl;
	global_log->info() << "                  dropradius: " << _Dropradius << std::endl;
	global_log->info() << "                  velocity: " << _Veloc << std::endl;
	global_log->info() << "                  starttime: " << _StartSimstep << std::endl;
	global_log->info() << "                  steps: " << _Steps << std::endl;
}

void Dropaccelerator::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
								   unsigned long simstep) {
	double corrVeloc = -(_Veloc / _Steps);

	if (_enabled) {
		global_log->debug() << "[Dropaccelerator] before forces called" << std::endl;

		if ((simstep - 1) % _INTERVAL != 0) {
			return;
		}

		int numberparticles = 0;

		// ITERATE OVER PARTICLES AND CHOOSE AND MARK IF MOLECULE IN DROPLET OR NOT
		if (simstep == _StartSimstep) {
			global_log->info() << "Startsimstep =" << simstep << endl;
			global_log->info() << "corrVeloc = " << corrVeloc << endl;

			// resize first
			_particleIsInDroplet.resize(global_simulation->getDomain()->getglobalNumMolecules(), false);

			for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm) {
				double distanceSquared[3];

				_particleIsInDroplet[tm->getID()] = false;
				distanceSquared[0] = (_xPosition - tm->r(0)) * (_xPosition - tm->r(0));
				distanceSquared[1] = (_yPosition - tm->r(1)) * (_yPosition - tm->r(1));
				distanceSquared[2] = (_zPosition - tm->r(2)) * (_zPosition - tm->r(2));

				if (distanceSquared[0] + distanceSquared[1] + distanceSquared[2] < _Dropradius * _Dropradius) {
					_particleIsInDroplet[tm->getID()] = true;
					tm->vadd(0, corrVeloc, 0);
					numberparticles++;
				}
			}

			// TODO: _particleIsInDroplet has to be communicated!


			// END TODO

			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(numberparticles);
			domainDecomp->collCommAllreduceSum();
			numberparticles = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();

			global_log->info() << "Marked Particles (Start)=" << numberparticles << endl;
		}

		// ITERATE OVER PARTICLES AND ACCELERATE ONLY DROPMOLECULES

		if (simstep > _StartSimstep && simstep < (_StartSimstep + _Steps)) {
			global_log->info() << "Accelerationsimstep = " << simstep << endl;

			for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm) {
				if (_particleIsInDroplet[tm->getID()]) {
					tm->vadd(0, corrVeloc, 0);
					numberparticles++;
				}
			}

			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(numberparticles);
			domainDecomp->collCommAllreduceSum();
			numberparticles = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();

			global_log->info() << "Marked Particles (Acceleration)=" << numberparticles << endl;
		}

		// CHECK IF VELOCITY HAS REACHED AND IF NOT ACCELERATE AGAIN

		if (simstep >= (_StartSimstep + _Steps) && simstep < (_StartSimstep + 2 * _Steps)) {
			_VelocNOW = 0;
			int Particles = 0;
			double deltavelocity;

			global_log->info() << "Postaccelerationsimstep = " << simstep << endl;

			// GET VELOCITY IN Y-DIRECTION FOR EVERY PARTICLE IN DROPLET

			for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm) {
				double partVelocity = tm->v(1);

				if (_particleIsInDroplet[tm->getID()]) {
					Particles += 1;
					_VelocNOW += partVelocity;
				}
			}

			// COMMUNICATION

			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendDouble(_VelocNOW);
			domainDecomp->collCommAllreduceSum();
			_VelocNOW = domainDecomp->collCommGetDouble();
			domainDecomp->collCommFinalize();

			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendInt(Particles);
			domainDecomp->collCommAllreduceSum();
			Particles = domainDecomp->collCommGetInt();
			domainDecomp->collCommFinalize();

			// CALCULATE AVERAGE SPEED
			_VelocNOW = -(_VelocNOW / Particles);

			global_log->info() << "_VelocNOW = " << _VelocNOW << endl;
			global_log->info() << "_VelocNOW-_Veloc = " << _VelocNOW - _Veloc << endl;
			global_log->info() << "corrVeloc = " << corrVeloc << endl;
			global_log->info() << "Marked Particles (Postacceleration)=" << Particles << endl;

			// ACCELERATE IF TOO SLOW

			if (_VelocNOW < _Veloc) {
				global_log->info() << "ACCELERATION with DELTAvelocity" << endl;
				deltavelocity = -(_Veloc - _VelocNOW);

				for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid();
					 ++tm) {
					if (_particleIsInDroplet[tm->getID()]) {
						tm->vadd(0, deltavelocity, 0);
					}
				}
			}
		}
	}
}

void Dropaccelerator::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
							  unsigned long simstep) {}
