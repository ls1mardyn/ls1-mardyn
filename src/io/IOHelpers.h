#pragma once

#include "molecules/Component.h"
#include "utils/Random.h"

class ParticleContainer;
class DomainDecompBase;

namespace IOHelpers {

/**
 * Removes the momentum from the simulation.
 * Will first sum up the overall momentum (m*v) and mass of the simulation and then apply a velocity shift
 * `-(momentum_sum / mass_sum)` to all particles which removes the momentum from the simulation.
 * @param particleContainer
 * @param components
 * @param domainDecomp
 */
void removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components,
					DomainDecompBase* domainDecomp);

/**
 * Assigns random velocities based on the given temperature to all particles in the given particleContainer.
 *
 * @param particleContainer
 * @param domainDecomp
 * @param temperature
 */
void initializeVelocityAccordingToTemperature(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
											  double temperature);

/**
 * Makes the particle ids of a container unique (throughout all ranks).
 * First an exclusive scan is used to get the number of particles of all processes with lower rank.
 * Each rank then adds that number to the ids of its particles.
 *
 * Assumptions that are made for the particles in particleContainer (for each rank independently, prior to the function
 * call):
 * - The ids of the particle range from 0 to particleContainer->getNumberOfParticles() - 1.
 * - There do not exist two particle with identical id.
 *
 * @param particleContainer
 * @param domainDecomp
 * @return The total number of particles (sum over all processes).
 */
unsigned long makeParticleIdsUniqueAndGetTotalNumParticles(ParticleContainer* particleContainer,
														   DomainDecompBase* domainDecomp);

/**
 * Generates a velocity vector according to the given temperature.
 * The generated velocity vector lies on the surface of a sphere with radius sqrt(3. * temperature).
 * If this function is called multiple times, the generated velocity vectors will be equally distributed on the surface
 * of the sphere.
 * @param temperature
 * @param RNG
 * @return
 */
std::array<double, 3> getRandomVelocityAccordingToTemperature(double temperature, Random& RNG);
}  // namespace IOHelpers
