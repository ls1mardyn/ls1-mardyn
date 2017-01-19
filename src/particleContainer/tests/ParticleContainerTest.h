/*
 * ParticleContainerTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef PARTICLECONTAINERTEST_H_
#define PARTICLECONTAINERTEST_H_

#include "molecules/Component.h"
#include "molecules/Molecule.h"
#include "utils/Testing.h"
#include "utils/TestWithSimulationSetup.h"
#include <vector>

class ParticleContainer;

/**
 * This class is intended as base class for any tests of particleContainers. However
 * it doesn't define a test suite but just provides actual test methods. For each
 * particleContainer there should be a subclass defining the suites.
 *
 * As this class should provide basic testing, it doesn't use any convenience methods
 * for setup.
 */
class ParticleContainerTest: public utils::TestWithSimulationSetup {

public:
	ParticleContainerTest();

	virtual ~ParticleContainerTest();

	/**
	 * Test insertion into a particle container.
	 *
	 * @param container ParticleContainer with corner points [0;0;0] x [10;10;10] and cutoff=2.5
	 */
	void testInsertion(ParticleContainer* container);

	/**
	 * Test the iterator methods of the particle container.
	 *
	 * Precondition: testInsertion has to be successful.
	 *
	 * @param container ParticleContainer with corner points [0;0;0] x [10;10;10] and cutoff=2.5
	 */
	void testMoleculeIteration(ParticleContainer* container);

	/**
	 * Test the update methods of the particle container:
	 * - iterates over the particles
	 * - changes positions and add new particles
	 * - check if update is done the right way
	 * - call deleteOuterparticles
	 *
	 * Precondition: testMoleculeIteration has to be successful.
	 *
	 * @param container ParticleContainer with corner points [0;0;0] x [10;10;10] and cutoff=2.5
	 */
	void testUpdateAndDeleteOuterParticles(ParticleContainer* container);

private:

	/**
	 * initialize the container with 4 molecules at [1/1/1] [2/2/2] [3/3/3] [5.1/5.1/5.1]
	 * so we have 1 pair in the same cell, 1 pair in the adjacent cell, and a single molecule.
	*/
	void setupMolecules(ParticleContainer* container);

	std::vector<Component> _components;
};

#endif /* PARTICLECONTAINERTEST_H_ */
