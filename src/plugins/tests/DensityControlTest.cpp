#include "DensityControlTest.h"

#include "utils/Testing.h"
#include "utils/CommVar.h"

TEST_SUITE_REGISTRATION(DensityControlTest);

DensityControlTest::DensityControlTest() {}

DensityControlTest::~DensityControlTest() {}

void DensityControlTest::testDensityControl() {

// So far, Autopas can only handle one component
// Therefore, the following tests are only conducted without Autopas
#ifndef MARDYN_AUTOPAS

	double delta = 1e-6;  // Tolerate deviation between expected and actual value

	// Number of particles per component
	// According to checkpoint file
	std::array<unsigned long, 5> numPrtl_before = {77, 45, 63, 46, 19};   // Before plugin action; total: 250
	// Expected
	std::array<unsigned long, 5> numPrtl_expect = {61, 45, 134, 46, 19};   // After plugin action;  total: 305
	// Counted
	CommVar<std::array<unsigned long, 5>> numPrtl_counted;

	for (unsigned short i = 0; i < 5; i++) {
		numPrtl_counted.global[i] = numPrtl_counted.local[i] = 0ul;
	}

	// Cutoff might be too small, but has no effect anyway
	std::unique_ptr<ParticleContainer> container{initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationMultiComponentMultiPotentials.inp", 6.0)};

	XMLfileUnits inp(getTestDataDirectory()+"/DensityControl.xml");

	std::unique_ptr<DensityControl> plugin {new DensityControl()};

	plugin->init(container.get(), _domainDecomposition, _domain);
	plugin->readXML(inp);
	plugin->beforeForces(container.get(), _domainDecomposition, 0);

	// Get number of particles per component
	for (auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		numPrtl_counted.local[it->componentid()]++;
	}

	// Get global values
#ifdef ENABLE_MPI
	MPI_Allreduce(numPrtl_counted.local.data(), numPrtl_counted.global.data(), numPrtl_counted.local.size(), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
	for (unsigned short i = 0; i < 5; i++) {
		numPrtl_counted.global[i] = numPrtl_counted.local[i];
	}
#endif

	// Compare to expected values
	for (uint32_t cid = 0; cid < 5; cid++) {
		ASSERT_EQUAL_MSG("Total number of particles of component "+std::to_string(cid+1)+" not as expected.", numPrtl_expect[cid], numPrtl_counted.global[cid]);
	}

#endif  //MARDYN_AUTOPAS
}
