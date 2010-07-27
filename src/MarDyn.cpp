#include <iostream>
#include <iomanip>
#include <ctime>

#ifdef CPPUNIT_TESTS
#include<cppunit/ui/text/TestRunner.h>
#include<cppunit/TestResultCollector.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#endif

#include "utils/OptionParser.h"
#include "utils/Logger.h"
#include "utils/compile_info.h"
#include "Simulation.h"


using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

//! execute unit tests
//! @return false if no errors occured, true otherwise
bool runTests();


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

#ifdef PARALLEL
	MPI_Init(&argc, &argv);
#endif
	/* Initialize the global log file */
	//string logfileName("MarDyn");
	//global_log = new Log::Logger(Log::All, logfileName);
	global_log = new Log::Logger(Log::Info);
#ifdef PARALLEL
	global_log->set_mpi_output_root(0);
#endif

	char *info_str = new char[MAX_INFO_STRING_LENGTH];
	get_compiler_info(&info_str);
	global_log->info() << "Compiler: " << info_str << endl;
	get_compile_time(&info_str);
	global_log->info() << "Compiled: " << info_str << endl;
#ifdef PARALLEL
	get_mpi_info(&info_str);
	global_log->info() << "MPI library: " << info_str << endl;
#endif
	get_timestamp(&info_str);
	global_log->info() << "Started: " << info_str << endl;
	get_host(&info_str);
	global_log->info() << "Execution host: " << info_str << endl;
#ifdef PARALLEL
	int world_size = 1;
	MPI_Comm_size( MPI_COMM_WORLD, &world_size );
	global_log->info() << "Running with " << world_size << " processes." << endl;
#endif



	cout.precision(6);

	OptionParser op;
	Values options = initOptions(argc, argv, op);
	vector<string> args = op.args();
	unsigned int numargs = args.size();

	bool tests = options.get("tests");
	if (tests) {
		bool testresult = runTests();
		exit(testresult);
	}

	if (numargs < 1) {
		op.print_usage();
		exit(1);
	}

	Simulation simulation(options, args);
	simulation.initialize();

	double runtime = double(clock()) / CLOCKS_PER_SEC;

	simulation.simulate();

	runtime = double(clock()) / CLOCKS_PER_SEC - runtime;

	cout << "main: used " << fixed << setprecision(2) << runtime << " s" << endl;
#ifdef PARALLEL
	MPI_Finalize();
#endif
}

bool runTests() {
#ifdef CPPUNIT_TESTS
	global_log->info() << "Running unit tests!" << endl;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TextUi::TestRunner runner;
	runner.addTest( registry.makeTest() );
	runner.run();

	const CppUnit::TestResultCollector& collector = runner.result();
	bool testresult = collector.testFailuresTotal() != 0;
	return testresult;
#else
	global_log->error() << endl << "Running unit tests demanded, but programme compiled without -DCPPUNIT_TESTS!" << endl << endl;
	return false;
#endif
}


Values& initOptions(int argc, char *argv[], OptionParser& op) {

	cout << "Size of args: " << argc << endl;

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


