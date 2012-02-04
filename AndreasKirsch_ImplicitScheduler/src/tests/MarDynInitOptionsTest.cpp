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

TEST_SUITE_REGISTRATION(MarDynInitOptionsTest);

// declaration of the function to be tested in Mardyn.cpp
optparse::Values& initOptions(int argc, const char* const argv[],
    optparse::OptionParser& op);


MarDynInitOptionsTest::MarDynInitOptionsTest() {
}

MarDynInitOptionsTest::~MarDynInitOptionsTest() {
}

void MarDynInitOptionsTest::testAllOptions() {
	const char* const argv[] = { "DummyProgrammeName", "-n", "4", "-t", "-v", "-p",
	    "test_prefix", "--phasespace-file", "phasespace.inp",
	    "--particle-container", "AdaptiveSubCells", "--cutoff-radius", "2.5",
	    "--cells-in-cutoff", "2", "--domain-decomposition", "KDDecomposition",
	    "--timestep-length", "0.00123" };

	int argc = sizeof(argv) / sizeof(char*);
	optparse::OptionParser op;
	optparse::Values options = initOptions(argc, argv, op);

	int n = options.get("timesteps");
	ASSERT_TRUE_MSG("timesteps isSetByUser must be true!", options.is_set_by_user("timesteps"));
	ASSERT_EQUAL(4, n);

	bool test = options.get("tests");
	ASSERT_TRUE_MSG("test option must be set", test);

	string prefix(options.get("outputprefix"));
	ASSERT_TRUE_MSG("outputprefix isSetByUser must be true!", options.is_set_by_user("outputprefix"));
	ASSERT_EQUAL(string("test_prefix"), prefix);

	bool verbose = options.get("verbose");
	ASSERT_TRUE_MSG("verbose option must be set", verbose);

	// TODO: how to query other commandline arguments!?

}

