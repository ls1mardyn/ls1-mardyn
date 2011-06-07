/*
 * TestWithSimulationSetup.cpp
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#include "TestWithSimulationSetup.h"
#include "utils/Logger.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif


using namespace Log;

#ifdef UNIT_TESTS

utils::TestWithSimulationSetup::TestWithSimulationSetup()
	: _rank(0), _domain(NULL), _domainDecomposition(NULL) { }



utils::TestWithSimulationSetup::~TestWithSimulationSetup() {
	if (_domain != NULL) {
		Log::global_log->warning() << "TestCase did not free it' ressources!" << std::endl;
		delete _domain;
	}

	if (_domainDecomposition != NULL) {
		Log::global_log->warning() << "TestCase did not free it' ressources!" << std::endl;
		delete _domainDecomposition;
	}
}


void utils::TestWithSimulationSetup::setUp() {
	_rank = 0;
	#ifdef ENABLE_MPI
		MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_rank) );
	#endif
	_domain = new Domain(_rank, NULL);

	#ifdef ENABLE_MPI
	_domainDecomposition = new DomainDecomposition();
	#else
	_domainDecomposition = new DomainDecompDummy();
	#endif
}


void utils::TestWithSimulationSetup::tearDown() {
	delete _domain;
	_domain = NULL;
	delete _domainDecomposition;
	_domainDecomposition = NULL;
}


ParticleContainer* utils::TestWithSimulationSetup::initializeFromFile(
		ParticleContainerFactory::type type, const char* fileName, double cutoff) {

	return ParticleContainerFactory::createInitializedParticleContainer(
			type, _domain, _domainDecomposition, cutoff, getTestDataFilename(fileName));
}

#endif
