/*
 * MarDynTest.cpp
 *
 * @Date: 12.07.2010
 * @Author: eckhardw
 */

#include "tests/MarDynInitOptionsTest.h"
#include "utils/OptionParser.h"

#include <string>

using namespace std;

// declaration of the function to be tested in Mardyn.cpp
optparse::Values& initOptions(int argc, char *argv[],
    optparse::OptionParser& op);

CPPUNIT_TEST_SUITE_REGISTRATION(MarDynInitOptionsTest);

MarDynInitOptionsTest::MarDynInitOptionsTest() {
}

MarDynInitOptionsTest::~MarDynInitOptionsTest() {
}

void MarDynInitOptionsTest::testAllOptions() {
	char* argv[] = { "DummyProgrammeName", "-n", "4", "-t", "-v", "-p",
	    "test_prefix", "--phasespace-file", "phasespace.inp",
	    "--particle-container", "AdaptiveSubCells", "--cutoff-radius", "2.5",
	    "--cells-in-cutoff", "2", "--domain-decomposition", "KDDecomposition",
	    "--timestep-length", "0.00123" };

	int argc = sizeof(argv) / sizeof(char*);
	optparse::OptionParser op;
	optparse::Values options = initOptions(argc, argv, op);

	int n = options.get("timesteps");
	CPPUNIT_ASSERT_ASSERTION_PASS_MESSAGE("timesteps isSetByUser must be true!", options.is_set_by_user("timesteps"));
	CPPUNIT_ASSERT_EQUAL(4, n);

	bool test = options.get("tests");
	CPPUNIT_ASSERT_ASSERTION_PASS_MESSAGE("test option must be set", test);

	string prefix(options.get("outputprefix"));
	CPPUNIT_ASSERT_ASSERTION_PASS_MESSAGE("outputprefix isSetByUser must be true!", options.is_set_by_user("outputprefix"));
	CPPUNIT_ASSERT_EQUAL(string("test_prefix"), prefix);

	bool verbose = options.get("verbose");
	CPPUNIT_ASSERT_ASSERTION_PASS_MESSAGE("verbose option must be set", verbose);

	// TODO: how to query other commandline arguments!?

}

