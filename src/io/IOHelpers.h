#pragma once

#include "molecules/Component.h"
#include "utils/Random.h"

class ParticleContainer;
class DomainDecompBase;

namespace IOHelpers {

void removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components,
					DomainDecompBase* domainDecomp);

void initializeVelocityAccordingToTemperature(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
											  double temperature);

unsigned long makeParticleIdsUniqueAndGetTotalNumParticles(ParticleContainer* particleContainer,
														   DomainDecompBase* domainDecomp);

std::array<double, 3> getRandomVelocityAccordingToTemperature(double temperature, Random& RNG);
}  // namespace IOHelpers