/*
 * RDFTest.cpp
 *
 * @Date: 15.02.2011
 * @Author: eckhardw
 */

#include "RDFTest.h"

#include "RDF.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif

#include <sstream>

using namespace std;

TEST_SUITE_REGISTRATION(RDFTest);

RDFTest::RDFTest() {
}

RDFTest::~RDFTest() {
}


void RDFTest::testRDFCountSequential12_LinkedCell() {
	delete _domainDecomposition;
	// will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompDummy();

	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-2x2x3.inp", 1.8);
	testRDFCountSequential12(moleculeContainer);
	delete moleculeContainer;
}


void RDFTest::testRDFCountSequential12_AdaptiveCell() {
	delete _domainDecomposition;
	// will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompDummy();

	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::AdaptiveSubCell, "1clj-regular-2x2x3.inp", 1.8);
	testRDFCountSequential12(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCountSequential12(ParticleContainer* moleculeContainer) {
	ParticlePairs2PotForceAdapter handler(*_domain);
	moleculeContainer->setPairHandler(&handler);

	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	/* The number of pairs counted by the RDF also depends on the particles in the halo.
	 * So count first with the halo being empty, and then being populated. */
	RDF rdf(0.018, 100, 1);
	handler.setRDF(&rdf);
	rdf.tickRDF();
	moleculeContainer->traversePairs();
	rdf.collectRDF(_domainDecomposition);

	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();
	rdf.tickRDF();

	// now the same with halo particles present.
	_domainDecomposition->exchangeMolecules(moleculeContainer, _domain->getComponents(), _domain);
	moleculeContainer->traversePairs();
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(40ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(44ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(4ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(4ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedDistribution[0][0][i]);
		}
	}
}


void RDFTest::testRDFCountLinkedCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-12x12x12.inp", 1.8);
	testRDFCount(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCountAdaptiveCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::AdaptiveSubCell, "1clj-regular-12x12x12.inp", 1.8);
	testRDFCount(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCount(ParticleContainer* moleculeContainer) {
	ParticlePairs2PotForceAdapter handler(*_domain);
	moleculeContainer->setPairHandler(&handler);

	_domainDecomposition->balanceAndExchange(true, moleculeContainer, _domain->getComponents(), _domain);
	moleculeContainer->updateMoleculeCaches();

	RDF rdf(0.018, 100, 1);
	handler.setRDF(&rdf);
	rdf.tickRDF();
	moleculeContainer->traversePairs();
	rdf.collectRDF(_domainDecomposition);

	// assert number of pairs counted
	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();

	rdf.tickRDF();
	moleculeContainer->traversePairs();
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
		stringstream msg;
		msg << "at index " << i;
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL_MSG(msg.str(), 0ul, rdf._globalDistribution[0][0][i]);
		}

		// the accumulated global distribution must be now twice the global distribution
		ASSERT_EQUAL(rdf._globalAccumulatedDistribution[0][0][i], 2 * rdf._globalDistribution[0][0][i]);
	}
}
