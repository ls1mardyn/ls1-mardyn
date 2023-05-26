/*
 * InMemoryCheckpointing.cpp
 *
 *  Created on: 3 Jul 2018
 *      Author: tchipevn
 */

#include "InMemoryCheckpointing.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

using Log::global_log;

void InMemoryCheckpointing::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 5;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_restartAtIteration = 10;
	xmlconfig.getNodeValue("restartAtIteration", _restartAtIteration);
	global_log->info() << "Restart at iteration (for development purposes): " << _restartAtIteration << std::endl;
}

void InMemoryCheckpointing::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep) {
	if (simstep != _restartAtIteration) {
		return;
	}
	global_log->info() << "InMemoryCheckpointWriter: resetting time to: " << _snapshot.getCurrentTime() << std::endl;
	Domain * domain = global_simulation->getDomain();

	// erase all current molecules
	particleContainer->clear();

	// fill new molecules
	particleContainer->addParticles(const_cast<std::vector<Molecule>& >(_snapshot.getMolecules()));

	// there should be no need to compute the forces again!
	// They are saved in the Molecule objects this time, due to the copy constructor.

	// Note that the forces, rotational moments are usually not saved in checkpoints and have to be recomputed in prepare_start()
	// so eventually, the following calls may be necessary:
	// * Simulation::updateParticleContainerAndDecomposition,
	// * ParticleContainer::traverseCells(cellProcessor)
	// * Simulation::updateForces()


	// set globals
	global_simulation->setSimulationTime(_snapshot.getCurrentTime());
	domain->setGlobalTemperature(_snapshot.getTemperature());

	mardyn_assert(_snapshot.getGlobalNumberOfMolecules() == domain->getglobalNumMolecules(true, &_simulation.getMoleculeContainers(), domainDecomp));
	mardyn_assert(domainDecomp->getRank() == _snapshot.getRank());

}

void InMemoryCheckpointing::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (simstep % _writeFrequency != 0) {
		return;
	}

	// else, write snapshot
	global_log->info() << "InMemoryCheckpointWriter: writing snapshot: " << std::endl;

	// put the molecules in the buffer
	_snapshot.clearMolecules();
	for (auto m = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); m.isValid(); ++m) {
		_snapshot.addMolecule(*m);
	}

	//set time, global number of molecules and temperature
	_snapshot.setCurrentTime(global_simulation->getSimulationTime());
	_snapshot.setGlobalNumberOfMolecules(domain->getglobalNumMolecules(true, &_simulation.getMoleculeContainers(), domainDecomp));
	_snapshot.setTemperature(domain->getGlobalCurrentTemperature());
	_snapshot.setRank(domainDecomp->getRank());

}
