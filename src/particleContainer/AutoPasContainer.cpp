/**
 * @file AutoPasContainer.cpp
 * @author seckler
 * @date 19.09.18
 */

#include "AutoPasContainer.h"

#include "autopas/AutoPas.h"

AutoPasContainer::AutoPasContainer() {
	autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autopasContainer();
}