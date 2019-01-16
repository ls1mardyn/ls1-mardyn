/*
 * DomainDecompBaseTest.cpp
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#include "parallel/tests/DomainDecompBaseTest.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include <set>

TEST_SUITE_REGISTRATION(DomainDecompBaseTest);

DomainDecompBaseTest::DomainDecompBaseTest() = default;

DomainDecompBaseTest::~DomainDecompBaseTest() = default;

void DomainDecompBaseTest::testNoDuplicatedParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
	unsigned long numMols = container->getNumberOfParticles();

	_domainDecomposition->exchangeMolecules(container, _domain);
	container->deleteOuterParticles();

	unsigned long newNumMols = container->getNumberOfParticles();
//	_domain->writeCheckpoint("dump.txt", container, _domainDecomposition);
	ASSERT_EQUAL(numMols, newNumMols);

	delete container;
	delete _domainDecomposition;
}

void DomainDecompBaseTest::testNoDuplicatedParticles() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15.inp", 27.0);
}


void DomainDecompBaseTest::testNoLostParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
	unsigned long numMols = container->getNumberOfParticles();

	double bBoxMin[3];
	double bBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		bBoxMin[dim] = container->getBoundingBoxMin(dim);
		bBoxMax[dim] = container->getBoundingBoxMax(dim);
	}
	std::set<unsigned long> lower[3];  // the id of particles that were close to the lower boundary in the specific dimension are stored here
	std::set<unsigned long> upper[3];  // the id of particles that were close to the upper boundary in the specific dimension are stored here

	for (auto m = container->iterator(); m.isValid(); ++m) {
		for (unsigned short dim = 0; dim < 3; dim++) {
			if (m->r(dim) < bBoxMin[dim] + cutoff / 2.) {
				// we shift particles close to the lower boundary to outside of the lower boundary.
				// in this case they are put to the smallest (in abs values) negative representable number
				// i.e. 2^(-149) = -1.4013e-45 for float resp. 4.94066e-324 for double
				m->setr(dim, std::nexttoward((vcp_real_calc) bBoxMin[dim], bBoxMin[dim] - 1.f));
				lower[dim].insert(m->getID());
			}
			if (m->r(dim) > bBoxMax[dim] - cutoff / 2.) {
				// We shift particles close to the upper boundary to outside of the upper boundary.
				// In this case they are put at minimum to boundingBoxMax, as this is no longer inside of the domain.
				// If the float representation of the maximum is less than the double representation, the next bigger floating point representation is used.
				// Otherwise the maximum is used.
				vcp_real_calc r = (float)bBoxMax[dim] >= bBoxMax[dim] ? bBoxMax[dim] : std::nexttoward((vcp_real_calc) bBoxMax[dim], bBoxMax[dim] + 1.f);
				m->setr(dim, r);
				upper[dim].insert(m->getID());
			}
		}
	}

	container->update();

	_domainDecomposition->exchangeMolecules(container, _domain);
	container->deleteOuterParticles();

	unsigned long newNumMols = container->getNumberOfParticles();
	//_domain->writeCheckpoint("dump.txt", container, _domainDecomposition, false);
	ASSERT_EQUAL(numMols, newNumMols);

	for (auto m = container->iterator(); m.isValid(); ++m) {
		for (int dim = 0; dim < 3; dim++) {
			if (lower[dim].count(m->getID())) {
				// We make sure, that these particles are now at the top part of the domain.
				ASSERT_TRUE(m->r(dim) >= bBoxMax[dim] - cutoff / 2.);
			} else if (upper[dim].count(m->getID())) {
				// We make sure, that these particles are now at the lower part of the domain.
				ASSERT_TRUE(m->r(dim) <= bBoxMin[dim] + cutoff / 2.);
			}
		}
	}

	delete _domainDecomposition;
	delete container;
}

void DomainDecompBaseTest::testNoLostParticles() {
	testNoLostParticlesFilename("H20_NaBr_0.01_T_293.15_DD_2.inp", 3.0);
}

void DomainDecompBaseTest::testExchangeMoleculesSimple() {

	// make sure we have a DomainDecompBase
	_domainDecomposition = new DomainDecompBase();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "LinkedCells.inp", 1.0);

	unsigned int count[8] = {0};
	ASSERT_EQUAL(8ul, container->getNumberOfParticles());

	auto m = container->iterator();
	while ( m.isValid()) {
		count[m->getID()]++;
		++m;
	}

	for (int i = 0; i < 8; i++) {
		ASSERT_EQUAL(1u, count[i]);
		count[i] = 0;
	}

	// after the exchange, every molecule should be replicated 7 times, i.e. there have to be 64 molecules in total
	_domainDecomposition->exchangeMolecules(container, _domain);
	ASSERT_EQUAL(64ul, container->getNumberOfParticles());

	m = container->iterator();
	while(m.isValid()) {
		count[m->getID()]++;
		++m;
	}

	for (int i = 0; i < 8; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
	delete _domainDecomposition;
}

void DomainDecompBaseTest::testExchangeMolecules() {

	// make sure we have a DomainDecompBase
	_domainDecomposition = new DomainDecompBase();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "DomainDecompBase.inp", 21.0);

	unsigned int count[3] = {0};
	ASSERT_EQUAL(3ul, container->getNumberOfParticles());

	auto m = container->iterator();
	while ( m.isValid()) {
		count[m->getID()]++;
		++m;
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(1u, count[i]);
		count[i] = 0;
	}

	// after the exchange, there have to be 21 copies in the halos, i.e. 24 molecules in total
	_domainDecomposition->exchangeMolecules(container, _domain);

	ASSERT_EQUAL(24ul, container->getNumberOfParticles());

	m = container->iterator();
	while(m.isValid()) {
		count[m->getID()]++;
		++m;
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
	delete _domainDecomposition;
}
