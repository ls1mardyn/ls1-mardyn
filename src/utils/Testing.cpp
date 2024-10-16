/*
 * Testing.cpp
 *
 * @Date: 11.08.2010
 * @Author: eckhardw
 */

#include "utils/Testing.h"

#include <memory>

#include "utils/Logger.h"
#include "utils/FileUtils.h"
#include "utils/mardyn_assert.h"

std::shared_ptr<Log::Logger> test_log;

#ifdef UNIT_TESTS
  #ifdef USE_CPPUNIT
    #include<cppunit/ui/text/TestRunner.h>
    #include<cppunit/TestResultCollector.h>
    #include <cppunit/extensions/TestFactoryRegistry.h>
    #include <cppunit/XmlOutputter.h>
	#include <cppunit/TestResult.h>
	#include <cppunit/BriefTestProgressListener.h>
  #else
    #include <tarch/tests/TestCaseRegistry.h>
    #include <tarch/tests/TestCaseCollection.h>
  #endif
#endif /* UNIT_TESTS */


int runTests(Log::logLevel testLogLevel, std::string& testDataDirectory, const std::string& testcases) {
	Log::logLevel globalLogLevel = Log::global_log->get_log_level();

	test_log = std::make_shared<Log::Logger>(testLogLevel);
	test_log->set_do_output(Log::global_log->get_do_output());
	if (testLogLevel > Log::Info) {
		Log::global_log->set_log_level(Log::Debug);
	} else {
		Log::global_log->set_log_level(Log::Warning);
	}

	setTestDataDirectory(testDataDirectory);

	int testresult;

#ifndef UNIT_TESTS
	test_log->error() << "Running unit tests demanded, but program was build without unit test support!" << std::endl;
	testresult = true;

#else /* UNIT_TESTS */
	test_log->info() << "Running unit tests" << std::endl;
#ifdef USE_CPPUNIT
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TextUi::TestRunner runner;
	runner.addTest( registry.makeTest() );
	CppUnit::BriefTestProgressListener progressListener;
	runner.eventManager().addListener( &progressListener );
	runner.run(testcases);

	CppUnit::TestResultCollector& collector = runner.result();
	testresult = collector.testFailuresTotal();

	std::ofstream stream("results.xml");
	CppUnit::XmlOutputter outputter( &collector, stream );
	outputter.write();

#else
	tarch::tests::TestCaseRegistry& registry = tarch::tests::TestCaseRegistry::getInstance();
	tarch::tests::TestCase& testCases = registry.getTestCaseCollection();
	testCases.run();
	testresult = testCases.getNumberOfErrors() != 0;
#endif
#endif /* UNIT_TESTS */

	Log::global_log->set_log_level(globalLogLevel);
	return testresult;
}

void setTestDataDirectory(std::string& testDataDirectory) {
#ifdef UNIT_TESTS
	utils::Test::setTestDataDirectory(testDataDirectory);
#endif /* UNIT_TESTS */
}


#ifdef UNIT_TESTS

std::string utils::Test::testDataDirectory("");

utils::Test::Test() { }

utils::Test::~Test() { }


void utils::Test::setTestDataDirectory(std::string& testDataDir) {
	if (!fileExists(testDataDir.c_str())) {
		std::ostringstream error_message;
		error_message << "Directory '" << testDataDir.c_str() << "' for test input data does not exist!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	testDataDirectory = testDataDir;
}


std::string utils::Test::getTestDataFilename(const std::string& file, bool checkExistence) {
	std::string fullPath = testDataDirectory +"/"+ file;

	if (!fileExists(fullPath.c_str()) and checkExistence) {
		std::ostringstream error_message;
		error_message << "File " << fullPath << " for test input data does not exist!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	return fullPath;
}

#endif /* UNIT_TESTS */
