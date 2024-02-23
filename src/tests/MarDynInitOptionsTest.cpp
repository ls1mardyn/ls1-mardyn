/*
 * MarDynTest.cpp
 *
 * @Date: 12.07.2010
 * @Author: eckhardw
 */

#include "tests/MarDynInitOptionsTest.h"
#include "utils/OptionParser.h"

#include <string>


TEST_SUITE_REGISTRATION(MarDynInitOptionsTest);

// declaration of the function to be tested in Mardyn.cpp
void initOptions(optparse::OptionParser *op);


MarDynInitOptionsTest::MarDynInitOptionsTest() {
}

MarDynInitOptionsTest::~MarDynInitOptionsTest() {
}

void MarDynInitOptionsTest::testAllOptions() {

	optparse::OptionParser op;
	initOptions(&op);

	const char* const argv[] = { "DummyProgrammeName", "-n", "4", "-t", "-v" };
	int argc = sizeof(argv) / sizeof(argv[0]);
	optparse::Values options = op.parse_args(argc, argv);

	ASSERT_TRUE_MSG("timesteps isSetByUser must be true!", options.is_set_by_user("timesteps"));
	ASSERT_EQUAL(4, (int) options.get("timesteps"));

	ASSERT_TRUE_MSG("test option must be set", options.is_set_by_user("tests"));

	ASSERT_TRUE_MSG("verbose option must be set", options.is_set_by_user("verbose"));

	// TODO: how to query other commandline arguments!?

}

