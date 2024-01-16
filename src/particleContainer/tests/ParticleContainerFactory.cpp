/*
 * ParticleContainerFactory.cpp
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#include "particleContainer/tests/ParticleContainerFactory.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"

#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif
#include "Domain.h"

#include "io/ASCIIReader.h"
#include "utils/Logger.h"

#ifdef MARDYN_AUTOPAS

#include "particleContainer/AutoPasContainer.h"
#endif

#include <io/BinaryReader.h>
#include <list>


ParticleContainer* ParticleContainerFactory::createEmptyParticleContainer(Type type) {
	if (type == LinkedCell) {
		double bBoxMin[] = {0.0, 0.0, 0.0, 0.0};
		double bBoxMax[] = {2.0, 2.0, 2.0, 2.0};
		double cutoffRadius = 1.0;
#ifndef MARDYN_AUTOPAS
		LinkedCells* container = new LinkedCells(bBoxMin, bBoxMax, cutoffRadius);
#else
		AutoPasContainer* container = new AutoPasContainer(cutoffRadius);
		container->rebuild(bBoxMin, bBoxMax);
#endif
		return container;

	} else {
		Log::global_log->error() << "ParticleContainerFactory: Unsupported type requested! " << std::endl;
		return nullptr;
	}
}



ParticleContainer* ParticleContainerFactory::createInitializedParticleContainer(
		Type type, Domain* domain, DomainDecompBase* domainDecomposition, double cutoff, const std::string& fileName, bool binary) {
	global_simulation->setcutoffRadius(cutoff);
	global_simulation->setLJCutoff(cutoff);

	std::unique_ptr<InputBase> inputReader;
	if (binary) {
		inputReader.reset(new BinaryReader());
		auto* binaryReader = dynamic_cast<BinaryReader*>(inputReader.get());
		binaryReader->setPhaseSpaceHeaderFile(fileName + ".header.xml");
		binaryReader->setPhaseSpaceFile(fileName + ".dat");
	} else {
		inputReader.reset(new ASCIIReader());
		auto* asciiReader = dynamic_cast<ASCIIReader*>(inputReader.get());
		asciiReader->setPhaseSpaceHeaderFile(fileName);
		asciiReader->setPhaseSpaceFile(fileName);
	}
	inputReader->readPhaseSpaceHeader(domain, 1.0);
	double bBoxMin[3];
	double bBoxMax[3];
	for (int i = 0; i < 3; i++) {
		bBoxMin[i] = domainDecomposition->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomposition->getBoundingBoxMax(i, domain);
	}

	ParticleContainer* moleculeContainer;
	if (type == Type::LinkedCell) {
#ifndef MARDYN_AUTOPAS
		moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, cutoff);
#else
		moleculeContainer = new AutoPasContainer(cutoff);
		moleculeContainer->rebuild(bBoxMin, bBoxMax);
#endif
		#ifdef ENABLE_MPI
		if (auto temp = dynamic_cast<DomainDecomposition *>(domainDecomposition)) {
			temp->initCommunicationPartners(cutoff, domain, moleculeContainer);
		}
		#endif
	} else {
		Log::global_log->error() << "ParticleContainerFactory: Unsupported type requested! " << std::endl;
		return nullptr;
	}

	inputReader->readPhaseSpace(moleculeContainer, domain, domainDecomposition);
	moleculeContainer->deleteOuterParticles();
	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	domain->initParameterStreams(cutoff, cutoff);
	return moleculeContainer;
}
