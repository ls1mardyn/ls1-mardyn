/*
 * Testing.cpp
 *
 * @Date: 11.08.2010
 * @Author: eckhardw
 */

#include "utils/Testing.h"
#include "utils/Logger.h"

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

