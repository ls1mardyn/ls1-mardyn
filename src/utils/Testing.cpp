/*
 * Testing.cpp
 *
 * @Date: 11.08.2010
 * @Author: eckhardw
 */

#include "utils/Testing.h"
#include "utils/Logger.h"
#include "utils/FileUtils.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"

#ifdef PARALLEL
#include "parallel/DomainDecomposition.h"
#endif

#ifdef UNIT_TESTS
  #ifdef USE_CPPUNIT
    #include<cppunit/ui/text/TestRunner.h>
    #include<cppunit/TestResultCollector.h>
    #include <cppunit/extensions/TestFactoryRegistry.h>
  #else
    #include <tarch/tests/TestCaseRegistry.h>
    #include <tarch/tests/TestCaseCollection.h>
  #endif
#endif

bool runTests() {

#ifdef UNIT_TESTS
	Log::global_log->info() << "Running unit tests!" << std::endl;
#ifdef USE_CPPUNIT
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TextUi::TestRunner runner;
	runner.addTest( registry.makeTest() );
	runner.run();

	const CppUnit::TestResultCollector& collector = runner.result();
	bool testresult = collector.testFailuresTotal() != 0;
	return testresult;
#else
	tarch::tests::TestCaseRegistry& registry = tarch::tests::TestCaseRegistry::getInstance();
	tarch::tests::TestCase& testCases = registry.getTestCaseCollection();
	testCases.run();
	bool testresult = testCases.getNumberOfErrors() != 0;
	return testresult;
#endif
#else
	Log::global_log->error() << std::endl << "Running unit tests demanded, but programme compiled without -DCPPUNIT_TESTS!" << std::endl << std::endl;
	return true;
#endif

}

void setTestDataDirectory(std::string& testDataDirectory) {
#ifdef UNIT_TESTS
	utils::Test::setTestDataDirectory(testDataDirectory);
#endif
}


#ifdef UNIT_TESTS

std::string utils::Test::testDataDirectory("");

utils::Test::Test() : _rank(0), _domain(NULL), _domainDecomposition(NULL) { }

utils::Test::~Test() {
	if (_domain != NULL) {
		Log::global_log->warning() << "TestCase did not free it' ressources!" << std::endl;
		delete _domain;
	}

	if (_domainDecomposition != NULL) {
		Log::global_log->warning() << "TestCase did not free it' ressources!" << std::endl;
		delete _domainDecomposition;
	}
}

void utils::Test::setUp() {
	_rank = 0;
	#ifdef PARALLEL
		MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_rank) );
	#endif
	_domain = new Domain(_rank, NULL);

	#ifdef PARALLEL
	_domainDecomposition = new DomainDecomposition();
	#else
	_domainDecomposition = new DomainDecompDummy();
	#endif
}


void utils::Test::tearDown() {
	delete _domain;
	_domain = NULL;
	delete _domainDecomposition;
	_domainDecomposition = NULL;
}


void utils::Test::setTestDataDirectory(std::string& testDataDir) {
	if (!fileExists(testDataDir.c_str())) {
		Log::global_log->error() << "Directory " << testDataDirectory << " for test input data does not exits!" << std::endl;
		exit(-1);
	}
	testDataDirectory = testDataDir;
}


std::string utils::Test::getTestDataFilename(const std::string& file) {
	std::string fullPath = testDataDirectory + file;

	if (!fileExists(fullPath.c_str())) {
		Log::global_log->error() << "File " << testDataDirectory << " for test input data does not exits!" << std::endl;
		exit(-1);
	}
	return fullPath;
}


ParticleContainer* utils::Test::initializeFromFile(
		ParticleContainerFactory::type type, const char* fileName, double cutoff) {

	return ParticleContainerFactory::createInitializedParticleContainer(
			type, _domain, _domainDecomposition, cutoff, getTestDataFilename(fileName));
}

#endif
