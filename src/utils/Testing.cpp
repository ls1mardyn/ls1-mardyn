/*
 * Testing.cpp
 *
 * @Date: 11.08.2010
 * @Author: eckhardw
 */

#include "utils/Testing.h"
#include "utils/Logger.h"
#include "utils/FileUtils.h"

Log::Logger* test_log;

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


bool runTests(Log::logLevel testLogLevel, std::string& testDataDirectory, const std::string& testcases) {
	Log::logLevel globalLogLevel = Log::global_log->get_log_level();

	test_log = new Log::Logger(testLogLevel);
	if (testLogLevel > Log::Info) {
		Log::global_log->set_log_level(Log::Debug);
	} else {
		Log::global_log->set_log_level(Log::Warning);
	}

	setTestDataDirectory(testDataDirectory);

	bool testresult;

#ifdef UNIT_TESTS
	test_log->info() << "Running unit tests!" << std::endl;
#ifdef USE_CPPUNIT
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TextUi::TestRunner runner;
	runner.addTest( registry.makeTest() );
	runner.run(testcases);

	const CppUnit::TestResultCollector& collector = runner.result();
	testresult = collector.testFailuresTotal() != 0;
#else
	tarch::tests::TestCaseRegistry& registry = tarch::tests::TestCaseRegistry::getInstance();
	tarch::tests::TestCase& testCases = registry.getTestCaseCollection();
	testCases.run();
	testresult = testCases.getNumberOfErrors() != 0;
#endif
#else
	test_log->error() << std::endl << "Running unit tests demanded, but programme compiled without -DCPPUNIT_TESTS!" << std::endl << std::endl;
	testresult = true;
#endif

	Log::global_log->set_log_level(globalLogLevel);
	return testresult;
}

void setTestDataDirectory(std::string& testDataDirectory) {
#ifdef UNIT_TESTS
	utils::Test::setTestDataDirectory(testDataDirectory);
#endif
}


#ifdef UNIT_TESTS

std::string utils::Test::testDataDirectory("");

utils::Test::Test() { }

utils::Test::~Test() { }


void utils::Test::setTestDataDirectory(std::string& testDataDir) {
	if (!fileExists(testDataDir.c_str())) {
		test_log->error() << "Directory " << testDataDirectory << " for test input data does not exits!" << std::endl;
		exit(-1);
	}
	testDataDirectory = testDataDir;
}


std::string utils::Test::getTestDataFilename(const std::string& file) {
	std::string fullPath = testDataDirectory + file;

	if (!fileExists(fullPath.c_str())) {
		test_log->error() << "File " << testDataDirectory << " for test input data does not exits!" << std::endl;
		exit(-1);
	}
	return fullPath;
}

#endif
