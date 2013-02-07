/*
 * CanonicalEnsembleTest.cpp
 *
 * @Date: 18.02.2011
 * @Author: eckhardw
 */

#include "CanonicalEnsembleTest.h"
#include "ensemble/CanonicalEnsemble.h"
#include "parallel/DomainDecompDummy.h"
#include "molecules/Molecule.h"
#include "Domain.h"

#include <vector>

using namespace std;

TEST_SUITE_REGISTRATION(CanonicalEnsembleTest);

CanonicalEnsembleTest::CanonicalEnsembleTest() { }

CanonicalEnsembleTest::~CanonicalEnsembleTest() { }


void CanonicalEnsembleTest::UpdateNumMoleculesSequential() {
// remove the ifndef when canonicalensemble can be tested in parallel
#ifndef ENABLE_MPI

	delete _domainDecomposition;
	// will be deleted by tearDown()
	_domainDecomposition = new DomainDecompDummy();

	// the halo is cleared for freshly initialized particle containers.
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::AdaptiveSubCell, "1clj-regular-12x12x12.inp", 1.0);
	vector<Component>& components(_domain->getComponents());
	CanonicalEnsemble ensemble(container, &components);

	ensemble.updateGlobalVariable(NUM_PARTICLES);
	// has the ensemble counted the right number of particles?
	ASSERT_EQUAL(1728ul, ensemble.N());
	// has the ensemble updated the count of particles per component right?
	ASSERT_EQUAL(1728ul, components[0].getNumMolecules());

	Molecule molecule(1729, 0, 5.5, 5.5, 5.5, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &components);
	container->addParticle(molecule);

	ensemble.updateGlobalVariable(NUM_PARTICLES);
	// has the ensemble counted the right number of particles?
	ASSERT_EQUAL(1729ul, ensemble.N());
	// has the ensemble updated the count of particles per component right?
	ASSERT_EQUAL(1729ul, components[0].getNumMolecules());

#endif
}


void CanonicalEnsembleTest::UpdateNumMoleculesParallel() {
	//ASSERT_FAIL("NOT YET IMPLEMENTED!");
}
