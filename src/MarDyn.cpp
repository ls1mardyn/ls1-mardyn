#include "Simulation.h"

#include "utils/OptionParser.h"
#include <iostream>
#include <iomanip>
#include <ctime>

#ifdef CPPUNIT_TESTS
#include<cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#endif

using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

void runTests();


optparse::Values& initOptions(int argc, char *argv[], optparse::OptionParser& op);


//! @page main
//! In this project, software for molecular dynamics simulation
//! with short-range forces is developed. The aim is to have a parallel code (MPI) 
//! for multi-centered molecules.
//!
//! The role of the main function is to run tests for all classes
//! and to instantiate an object of the Simulation class which
//! is actually responsible for the simulation
//!
int main(int argc, char** argv) {

  cout.precision(6);

  OptionParser op;
  Values options = initOptions(argc, argv, op);
  vector<string> args = op.args();
  unsigned int numargs = args.size();

  if (numargs < 1) {
    op.print_usage();
    exit(1);
  }

  bool tests = options.get("tests");
  if (tests) {
    runTests();
  }

  Simulation simulation(options, args);
  simulation.initialize();

  double runtime = double(clock()) / CLOCKS_PER_SEC;

  simulation.simulate();

  runtime = double(clock()) / CLOCKS_PER_SEC - runtime;

  cout << "main: used " << fixed << setprecision(2) << runtime << " s" << endl;
}


void runTests() {
#ifdef CPPUNIT_TESTS
	cout << "Running unit tests!" << endl;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( registry.makeTest() );
  runner.run();
#else
  cout << endl << "Running unit tests demanded, but programme compiled without -DCPPUNIT_TESTS!" << endl << endl;
#endif
}


Values& initOptions(int argc, char *argv[], OptionParser& op) {

	op = OptionParser()
        						// The last two optional positional arguments are only here for backwards-compatibility
        						.usage("%prog [-n steps] [-p prefix] <configfilename> [<number of timesteps>] [<outputprefix>]\n\nUse option --help to display all available options.")
        						.version("%prog 1.0")
        						.description("MarDyn is a MD simulator. All behavior is controlled via the config file.")
        						// .epilog("background info?")
        						;

	op.add_option("-n", "--steps") .dest("timesteps") .metavar("NUM") .type("int") .set_default(1) .help("number of timesteps to simulate (default: %default)");
	op.add_option("-p", "--outprefix") .dest("outputprefix") .metavar("STR") .help("prefix for output files");
	op.add_option("-v", "--verbose") .action("store_true") .dest("verbose") .metavar("V") .type("bool") .set_default(false) .help("verbose mode: print debugging information (default: %default)");
	op.add_option("-t", "--tests") .action("store_true") .dest("tests") .metavar("T") .type("bool") .set_default(false) .help("unit tests: run built-in unit tests (default: %default)");

	OptionGroup dgroup = OptionGroup(op, "Developer options", "Advanced options for developers and experienced users.");
	dgroup.add_option("--phasespace-file") .metavar("FILE") .help("path to file containing phase space data");
	char const* const pc_choices[] = { "LinkedCells", "AdaptiveSubCells" };
	dgroup.add_option("--particle-container") .choices(&pc_choices[0], &pc_choices[2]) .set_default(pc_choices[0]) .help("container used for locating nearby particles (default: %default)");
	dgroup.add_option("--cutoff-radius") .type("float") .set_default(5.0) .help("radius of sphere around a particle in which forces are considered (default: %default)");
	dgroup.add_option("--cells-in-cutoff") .type("int") .set_default(2) .help("number of cells in cutoff-radius cube (default: %default); only used by LinkedCells particle container");
	char const* const dd_choices[] = { "DomainDecomposition", "KDDecomposition" };
	dgroup.add_option("--domain-decomposition") .choices(&dd_choices[0], &dd_choices[2]) .set_default(dd_choices[0]) .help("domain decomposition strategy for MPI (default: %default)");
	dgroup.add_option("--timestep-length") .type("float") .set_default(0.004) .help("length of one timestep in TODO (default: %default)");
	op.add_option_group(dgroup);

	return op.parse_args(argc, argv);
}


