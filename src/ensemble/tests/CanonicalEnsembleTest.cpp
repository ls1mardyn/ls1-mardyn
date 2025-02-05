/*
 * CanonicalEnsembleTest.cpp
 *
 * @Date: 18.02.2011
 * @Author: eckhardw
 */

#include "CanonicalEnsembleTest.h"
#include "ensemble/CanonicalEnsemble.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "Simulation.h"

#include <vector>

TEST_SUITE_REGISTRATION(CanonicalEnsembleTest);

CanonicalEnsembleTest::CanonicalEnsembleTest() { }

CanonicalEnsembleTest::~CanonicalEnsembleTest() { }


void CanonicalEnsembleTest::UpdateNumMoleculesSequential() {
// remove the ifndef when canonicalensemble can be tested in parallel
#ifndef ENABLE_MPI

	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	// the halo is cleared for freshly initialized particle containers.
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-12x12x12.inp", 1.0);
	CanonicalEnsemble ensemble;

	ensemble.addComponent((*(global_simulation->getEnsemble()->getComponents()))[0]);
	Component* component = ensemble.getComponent(0);


	ensemble.updateGlobalVariable(container, NUM_PARTICLES);
	// has the ensemble counted the right number of particles?
	ASSERT_EQUAL(1728ul, ensemble.N());
	// has the ensemble updated the count of particles per component right?
	ASSERT_EQUAL(1728ul, component->getNumMolecules());

	Molecule molecule(1729, component, 5.5, 5.5, 5.5, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
	container->addParticle(molecule);

	ensemble.updateGlobalVariable(container, NUM_PARTICLES);
	// has the ensemble counted the right number of particles?
	ASSERT_EQUAL(1729ul, ensemble.N());
	// has the ensemble updated the count of particles per component right?
	ASSERT_EQUAL(1729ul, component->getNumMolecules());

	delete _domainDecomposition;
	delete container;
#endif
}


void CanonicalEnsembleTest::UpdateNumMoleculesParallel() {
	//ASSERT_FAIL("NOT YET IMPLEMENTED!");
}
