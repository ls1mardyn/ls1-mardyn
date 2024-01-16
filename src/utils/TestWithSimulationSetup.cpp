/*
 * TestWithSimulationSetup.cpp
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#include "TestWithSimulationSetup.h"
#include "utils/Logger.h"

#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif



#ifdef UNIT_TESTS

utils::TestWithSimulationSetup::TestWithSimulationSetup()
	: _rank(0), _domain(nullptr), _domainDecomposition(nullptr) { }



utils::TestWithSimulationSetup::~TestWithSimulationSetup() {
}


void utils::TestWithSimulationSetup::setUp() {
	_rank = 0;
	Simulation* sim = new Simulation();
	//deleted the next line, as it is already done in the constructor.
	//sim->initialize(); // this assigns global_simulation.
	#ifdef ENABLE_MPI
		MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_rank) );
	#endif
	_domain = global_simulation->getDomain();
	_domainDecomposition = & global_simulation->domainDecomposition();
}


void utils::TestWithSimulationSetup::tearDown() {
	delete global_simulation;
	global_simulation = nullptr;
	_domain = nullptr;
	_domainDecomposition = nullptr;
}


ParticleContainer* utils::TestWithSimulationSetup::initializeFromFile(
		ParticleContainerFactory::Type type, const std::string& fileName, double cutoff, bool binary) {

	bool checkExist = not binary;  // for binary the prefix is given, so we cannot check the file for existence.
	return ParticleContainerFactory::createInitializedParticleContainer(
			type, _domain, _domainDecomposition, cutoff, getTestDataFilename(fileName, checkExist), binary);
}

#endif /* UNIT_TESTS */
