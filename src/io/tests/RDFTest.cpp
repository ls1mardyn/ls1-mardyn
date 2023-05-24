/*
 * RDFTest.cpp
 *
 * @Date: 15.02.2011
 * @Author: eckhardw
 */

#include "RDFTest.h"

#include "io/RDF.h"
#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/RDFCellProcessor.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif

#include <sstream>


#if !defined(ENABLE_REDUCED_MEMORY_MODE) && !defined(MARDYN_AUTOPAS)
TEST_SUITE_REGISTRATION(RDFTest);
#else
#pragma message "Compilation info: RDFTest disabled in reduced memory mode and autopas mode"
#endif


RDFTest::RDFTest() {
}

RDFTest::~RDFTest() {
}


void RDFTest::initRDF(RDF &rdf, double intervalLength, unsigned int bins, std::vector<Component>* components) {
	rdf._intervalLength = intervalLength;
	rdf._bins = bins;
	rdf._components = components;
	rdf._writeFrequency = 25000;
	rdf._outputPrefix = "out";
	rdf._readConfig=true;
	rdf.init();
}

void RDFTest::testRDFCountSequential12_LinkedCell() {
	// original pointer will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompBase();

	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-2x2x3.inp", 1.8);
	testRDFCountSequential12(moleculeContainer);
	delete moleculeContainer;
	delete _domainDecomposition;
}

void RDFTest::testRDFCountSequential12(ParticleContainer* moleculeContainer) {

	ParticlePairs2PotForceAdapter handler(*_domain);
	double cutoff = moleculeContainer->getCutoff();
	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	ASSERT_EQUAL((size_t) 1, components->size());

	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	/* The number of pairs counted by the RDF also depends on the particles in the halo.
	 * So count first with the halo being empty, and then being populated. */
	RDF rdf;
	initRDF(rdf, 0.018, 100, components);
	rdf.tickRDF();
	RDFCellProcessor cellProcessor(cutoff, &rdf);
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);

	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._distribution.global[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._distribution.global[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._distribution.global[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._distribution.global[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();
	rdf.tickRDF();

	// now the same with halo particles present.
	_domainDecomposition->exchangeMolecules(moleculeContainer, _domain);
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
//		std::cout << " i=" << i << " global = " <<rdf._distribution.global[0][0][i] << " acc="
//				<< rdf._globalAccumulatedDistribution[0][0][i] << std::endl;
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._distribution.global[0][0][i]);
			ASSERT_EQUAL(40ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._distribution.global[0][0][i]);
			ASSERT_EQUAL(44ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 88) {
			ASSERT_EQUAL(4ul, rdf._distribution.global[0][0][i]);
			ASSERT_EQUAL(4ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._distribution.global[0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._distribution.global[0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedDistribution[0][0][i]);
		}
	}
}


void RDFTest::testRDFCountLinkedCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-12x12x12.inp", 1.8);
	testRDFCount(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCount(ParticleContainer* moleculeContainer) {
	ParticlePairs2PotForceAdapter handler(*_domain);
	double cutoff = moleculeContainer->getCutoff();
	
	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	ASSERT_EQUAL((size_t) 1, components->size());

	moleculeContainer->deleteOuterParticles();
	_domainDecomposition->balanceAndExchange(1.0, false, moleculeContainer, _domain);
	moleculeContainer->updateMoleculeCaches();

	RDF rdf;
	RDFCellProcessor cellProcessor(cutoff, &rdf);
	initRDF(rdf, 0.018, 100, components);

	rdf.tickRDF();
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);

	// assert number of pairs counted
	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._distribution.global[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._distribution.global[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._distribution.global[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._distribution.global[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._distribution.global[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();

	moleculeContainer->deleteOuterParticles();
	_domainDecomposition->balanceAndExchange(1.0, false, moleculeContainer, _domain);
	moleculeContainer->updateMoleculeCaches();

	rdf.tickRDF();
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
		std::stringstream msg;
		msg << "at index " << i;
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._distribution.global[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._distribution.global[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._distribution.global[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._distribution.global[0][0][i]);
		} else {
			ASSERT_EQUAL_MSG(msg.str(), 0ul, rdf._distribution.global[0][0][i]);
		}

		// the accumulated global distribution must be now twice the global distribution
		ASSERT_EQUAL(rdf._globalAccumulatedDistribution[0][0][i], 2 * rdf._distribution.global[0][0][i]);
	}
}



void RDFTest::testSiteSiteRDFLinkedCell() {
	if (_domainDecomposition->getNumProcs() > 8) {
		ASSERT_FAIL("RUN THIS TEST WITH <= 8 PROCESSORS!");
	}
	std::unique_ptr<ParticleContainer> moleculeContainer{initializeFromFile(ParticleContainerFactory::LinkedCell, "2clj-regular.inp", 3.5)};
	testSiteSiteRDF(moleculeContainer.get());
}

void RDFTest::testSiteSiteRDF(ParticleContainer* moleculeContainer) {

	std::cout << "RDFTest::testSiteSiteRDF with " <<_domainDecomposition->getNumProcs()<<" procs" << std::endl;

	ParticlePairs2PotForceAdapter handler(*_domain);
	double cutoff = moleculeContainer->getCutoff();

	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	ASSERT_EQUAL((size_t) 1, components->size());

	_domainDecomposition->balanceAndExchange(1.0, true, moleculeContainer, _domain);
	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	RDF rdf;
	initRDF(rdf, 0.05, 101, components);
	rdf.tickRDF();
	RDFCellProcessor cellProcessor(cutoff, &rdf);
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);

	for (int i = 0; i < 101; i++) {
//		std::cout << "Bin " << i << ": " << rdf._siteDistribution.global[0][0][0][0][i] <<
//					", " << rdf._siteDistribution.global[0][0][0][1][i] <<
//					", " << rdf._siteDistribution.global[0][0][1][0][i] <<
//					", " << rdf._siteDistribution.global[0][0][1][1][i] << std::endl;
		if (i == 20) {
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][1][i]);
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
		} else if (i == 60) {
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][0][i]);
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][1][i]);
		} else if (i == 100) {
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][1][i]);
			ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();
	rdf.tickRDF();

	// test the accumulation of counts...
	moleculeContainer->traverseCells(cellProcessor);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 101; i++) {
			if (i == 20) {
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][1][i]);
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
				// accumulated counts
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
			} else if (i == 60) {
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][0][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][1][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][0][i]);
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][1][i]);
				// accumulated counts
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
			} else if (i == 100) {
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][0][1][i]);
				ASSERT_EQUAL(16ul, rdf._siteDistribution.global[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
				// accumulated counts
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
				ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
			} else {
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][0][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][0][1][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._siteDistribution.global[0][0][1][1][i]);
				// accumulated counts
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
				ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
			}
		}

}

