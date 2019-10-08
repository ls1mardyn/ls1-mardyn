/*
 * DomainDecompositionTest.cpp
 *
 * @Date: 10.05.2012
 * @Author: eckhardw
 */

#include "DomainDecompositionTest.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecomposition.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"
#include "Domain.h"

TEST_SUITE_REGISTRATION(DomainDecompositionTest);

DomainDecompositionTest::DomainDecompositionTest() = default;

DomainDecompositionTest::~DomainDecompositionTest() = default;

void DomainDecompositionTest::testNoDuplicatedParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown() (delete global_simulation)
	std::unique_ptr<DomainDecomposition> _domainDecomposition {new DomainDecomposition()};

	std::unique_ptr<ParticleContainer> container{
		initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff)};
	auto numMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(numMols);
	_domainDecomposition->collCommAllreduceSum();
	numMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();

	_domainDecomposition->balanceAndExchange(0., true, container.get(), _domain);
	container->deleteOuterParticles();

	auto newNumMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(newNumMols);
	_domainDecomposition->collCommAllreduceSum();
	newNumMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();

	ASSERT_EQUAL(numMols, newNumMols);
}

void DomainDecompositionTest::testNoDuplicatedParticles() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15_DD.inp", 5.0);
}

void DomainDecompositionTest::testNoLostParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown() (delete global_simulation)
	std::unique_ptr<DomainDecomposition> _domainDecomposition{new DomainDecomposition()};

	std::unique_ptr<ParticleContainer> container{
		initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff)};
	auto numMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(numMols);
	_domainDecomposition->collCommAllreduceSum();
	numMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();


	double bBoxMin[3];
	double bBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		bBoxMin[dim] = 0.;
		bBoxMax[dim] = _domain->getGlobalLength(dim);
	}
	std::set<unsigned long> lower[3];  // the id of particles that were close to the lower boundary in the specific dimension are stored here
	std::set<unsigned long> upper[3];  // the id of particles that were close to the upper boundary in the specific dimension are stored here

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
		for (int dim = 0; dim < 3; dim++) {
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

	_domainDecomposition->balanceAndExchange(0., true, container.get(), _domain);
	container->deleteOuterParticles();

	auto newNumMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(newNumMols);
	_domainDecomposition->collCommAllreduceSum();
	newNumMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();

	//_domain->writeCheckpoint("dump.txt", container, _domainDecomposition, false);
	ASSERT_EQUAL(numMols, newNumMols);

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
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
}

void DomainDecompositionTest::testNoLostParticles() {
	testNoLostParticlesFilename("H20_NaBr_0.01_T_293.15_DD_2.inp", 3.0);
}

void DomainDecompositionTest::testExchangeMolecules1Proc() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "DomainDecompBase.inp", 21.0);

	unsigned int count[3] = {0};
	ASSERT_EQUAL(3ul, container->getNumberOfParticles());

	auto m = container->iterator(ParticleIterator::ALL_CELLS);
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

	//ASSERT_EQUAL(24ul, container->getNumberOfParticles());
	size_t molCount = 0;
	for(auto iter = container->iterator(ParticleIterator::ALL_CELLS); iter.isValid(); ++iter) ++molCount;
	ASSERT_EQUAL(24ul, molCount);

	m = container->iterator(ParticleIterator::ALL_CELLS);
	while(m.isValid()) {
		count[m->getID()]++;
		++m;
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
}
