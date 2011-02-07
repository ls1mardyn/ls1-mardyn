/*
 * ParticleContainerFactory.cpp
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#include "particleContainer/tests/ParticleContainerFactory.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/LinkedCellsCUDA.h"

#include "utils/Logger.h"

using namespace Log;

ParticleContainer* ParticleContainerFactory::createEmptyParticleContainer(type type) {
	if (type == LinkedCell) {
		double bBoxMin[] = {0.0, 0.0, 0.0, 0.0};
		double bBoxMax[] = {2.0, 2.0, 2.0, 2.0};
		double cutoffRadius = 1.0;
		double LJCutoffRadius = 1.0;
		double tersoffCutoffRadius = 1.0;
		double cellsInCutoffRadius = 1.0;

		LinkedCells* container = new LinkedCells(bBoxMin, bBoxMax, cutoffRadius, LJCutoffRadius,
		                                        tersoffCutoffRadius, cellsInCutoffRadius, NULL);
		return container;

	} else {
		global_log->error() << "ParticleContainerFactory: Unsupported type requested! " << std::endl;
		return NULL;
	}
}
