//
// Created by seckler on 19.09.18.
//

#include "AutoPasContainer.h"

#include "autopas/AutoPas.h"

AutoPasContainer::AutoPasContainer() {
	autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autopasContainer();
}