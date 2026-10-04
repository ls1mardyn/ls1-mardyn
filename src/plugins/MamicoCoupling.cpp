#ifdef MAMICO_COUPLING

#include "MamicoCoupling.h"

#include "Domain.h"

void MamicoCoupling::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	Log::global_log->info() << "MaMiCo coupling plugin initialized" << std::endl;
}

void MamicoCoupling::beforeForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
								  unsigned long simstep) {
	if (_couplingEnabled) {
		// This object should be set by MaMiCo after the plugins are created in the simulation readxml file.
		// Even though this method is called before the object is set, at this point the coupling switch is always off
		_couplingCellService->processInnerCouplingCellAfterMDTimestep();
		// Before forces are calculated, mamico would like to insert and delete particles using usher.
		// At this point in time, ls1 linked cells should have no halos, and leaving particles should not have been
		// communicated. Particle container is hence updated to make sure halos exist.
		if (_couplingCellService->checkIfDistributeMass(simstep)) {
			global_simulation->updateParticleContainerAndDecomposition(1.0, false);
			// Particle insertion and deletion
			_couplingCellService->distributeMass(simstep);
		}
		// Mamico thermostat
		_couplingCellService->applyTemperatureToMolecules(simstep);
		// Remove halos, as at this point in the code there should be no halos
#ifndef MARDYN_AUTOPAS
		if (_couplingCellService->checkIfDistributeMass(simstep)) {
			particleContainer->deleteOuterParticles();
		}
#endif
		// Unfortunately, at this point, leaving particles should exist on source ranks and not destination
		// Thus this plugin does not entirely preserve the simulation state, since leaving particles are communicated
		// Any plugins after this that work on leaving particles will find zero leaving particles
	}
}

void MamicoCoupling::afterForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
								 unsigned long simstep) {
	if (_couplingEnabled) {
		_couplingCellService->distributeMomentum(simstep);
		_couplingCellService->applyBoundaryForce(simstep);
	}
}

#endif	// MAMICO_COUPLING
