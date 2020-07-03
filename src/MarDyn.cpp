#include "MarDyn_version.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "WrapOpenMP.h"

#include "Simulation.h"
#include "utils/compile_info.h"
#include "utils/PrintThreadPinningToCPU.h"
#include "utils/FileUtils.h"
#include "utils/Logger.h"
#include "utils/OptionParser.h"
#include "utils/Testing.h"
#include "utils/Timer.h"
#include "utils/SigsegvHandler.h"

using Log::global_log;
using std::endl;
using optparse::Values;

/**
 * @brief Initialize command line options.
 */
void initOptions(optparse::OptionParser *op) {
	op->usage("%prog [OPTIONS] <configfilename>\n\n"
		"Use option --help to display all available options.\n\n"
		"To execute the built in unit tests run with\n"
		"%prog --tests --test-dir <test input data directory> [<name of testcase>]");
	op->version("%prog " + MARDYN_VERSION);
	op->description("ls1-MarDyn (Large Scale SImulation MoleculAR DYNamics)");
	op->add_option("-a", "--loop-abort-time").dest("loop-abort-time").type("float") .metavar("TIME") .set_default(-1) .help("(optional) max walltime allowed in (s) before stop of main loop (default: %default)");
	op->add_option("-n", "--steps").dest("timesteps").type("int") .metavar("NUM") .set_default(1) .help("number of timesteps to simulate (default: %default)");
	op->add_option("-p", "--outprefix").dest("outputprefix").type("string") .metavar("STR") .set_default("MarDyn") .help("default prefix for output files (default: %default)");
	op->add_option("-v", "--verbose").dest("verbose").type("bool") .action("store_true") .set_default(false) .help("verbose mode: print debugging information (default: %default)");
	op->add_option("--print-meminfo").dest("print-meminfo").type("bool").action("store_true").set_default(false).help("Print memory consumtion info (default: %default)");
	op->add_option("--logfile").dest("logfile").type("string").metavar("PREFIX").set_default("MarDyn").help("enable output to logfile using given prefix for the filename (default: %default)");
	op->add_option("--legacy-cell-processor").dest("legacy-cell-processor").type("bool").action("store_true").set_default(false).help("use legacyCellProcessor (AoS) (default: %default)");
	op->add_option("--final-checkpoint").dest("final-checkpoint").type("int").metavar("(1|0)").set_default(1).help("enable/disable final checkopint (default: %default)");
	op->add_option("--timed-checkpoint").dest("timed-checkpoint").type("float").metavar("TIME").set_default(-1).help("Execution time of the simulation in seconds after which a checkpoint is forced, disable: -1. (default: %default)");
#ifdef ENABLE_SIGHANDLER
	op->add_option("-S", "--sigsegvhandler").dest("sigsegvhandler").type("bool") .action("store_true") .set_default(false) .help("sigsegvhandler: prints stacktrace on sigsegv (default: %default)");
#endif

	op->add_option("-t", "--tests").action("store_true").dest("tests").type("bool").set_default(false).help("unit tests: run built-in unit tests instead of regular simulation");
	op->add_option("-d", "--test-dir").dest("testDataDirectory").type("string").metavar("STR").set_default("").help("unit tests: specify the directory where the in input data required by the tests resides");
}

/**
 * @brief Helper function outputting program build information to given logger
 */
void program_build_info(Log::Logger *log) {
	log->info() << "Compilation info:" << endl;

	char info_str[MAX_INFO_STRING_LENGTH];
	get_compiler_info(info_str);
	log->info() << "	Compiler:	" << info_str << endl;
	get_compile_time(info_str);
	log->info() << "	Compiled on:	" << info_str << endl;
	get_precision_info(info_str);
	log->info() << "	Precision:	" << info_str << endl;
	get_intrinsics_info(info_str);
	log->info() << "	Intrinsics:	" << info_str << endl;
	get_rmm_normal_info(info_str);
	log->info() << "	RMM/normal:	" << info_str << endl;
	get_openmp_info(info_str);
	log->info() << "	OpenMP:		" << info_str << endl;
	get_mpi_info(info_str);
	log->info() << "	MPI:		" << info_str << endl;
}

/**
 * @brief Helper function outputting program invocation information to given logger
 */
void program_execution_info(int argc, char **argv, Log::Logger *log) {
	log->info() << "Execution info:" << endl;

	char info_str[MAX_INFO_STRING_LENGTH];
	get_timestamp(info_str);
	log->info() << "	Started: " << info_str << endl;
	get_host(info_str);
	log->info() << "	Execution host: " << info_str << endl;
	std::stringstream arguments;
	for (int i = 0; i < argc; i++) {
		arguments << " " << argv[i];
	}
	log->info() << "	Started with arguments: " << arguments.str() << endl;

#if defined(_OPENMP)
	int num_threads = mardyn_get_max_threads();
	global_log->info() << "	Running with " << num_threads << " OpenMP threads." << endl;
	// print thread pinning info
	PrintThreadPinningToCPU();
#endif

#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
	global_log->info() << "	Running with " << world_size << " MPI processes." << endl;
#endif
}

/** Run the internal unit tests */
int run_unit_tests(const Values &options, const vector<string> &args) {
	string testcases("");
	if(args.size() == 1) {
		testcases = args[0];
		global_log->info() << "Running unit tests: " << testcases << endl;
	} else {
		global_log->info() << "Running all unit tests!" << endl;
	}
	std::string testDataDirectory(options.get("testDataDirectory"));
	global_log->info() << "Test data directory: " << testDataDirectory << endl;
	Log::logLevel testLogLevel = (options.is_set("verbose") && options.get("verbose")) ? Log::All : Log::Info;
	int testresult = runTests(testLogLevel, testDataDirectory, testcases);
	return testresult;
}

/**
 * The role of the main function is to instantiate an object of the Simulation
 * class which is actually responsible for the simulation and to run tests for
 * all classes.
 */
int main(int argc, char** argv) {
#ifdef ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif

	/* Initialize the global log file */
	global_log = new Log::Logger(Log::Info);
#ifdef ENABLE_MPI
	global_log->set_mpi_output_root(0);
	//global_log->set_mpi_output_all();
#endif

	optparse::OptionParser op;
	initOptions(&op);
	optparse::Values options = op.parse_args(argc, argv);
	vector<string> args = op.args();

	global_log->info() << "Running ls1-MarDyn version " << MARDYN_VERSION << endl;
#ifndef NDEBUG
	global_log->warning() << "This ls1-MarDyn binary is a DEBUG build!" << endl;
#endif

	if( options.is_set_by_user("logfile") ) {
		string logfileNamePrefix(options.get("logfile"));
		global_log->info() << "Using logfile with prefix " << logfileNamePrefix << endl;
		delete global_log;
		global_log = new Log::Logger(Log::Info, logfileNamePrefix);
	}
	if( options.is_set_by_user("verbose") ) {
		global_log->info() << "Enabling verbose log output." << endl;
		global_log->set_log_level(Log::All);
	}
#ifdef ENABLE_SIGHANDLER
	if (options.is_set_by_user("sigsegvhandler")) {
		global_log->info() << "Enabling sigsegvhandler." << endl;
		registerSigsegvHandler();  // from SigsegvHandler.h
	}
#endif
	program_build_info(global_log);
	program_execution_info(argc, argv, global_log);


	/* Run built in tests and exit */
	if (options.is_set_by_user("tests")) {
		int testresult = run_unit_tests(options, args);
		#ifdef ENABLE_MPI
		MPI_Finalize();
		#endif
		exit(testresult); // using exit here should be OK
	}


	/* Set up and run regular Simulation */
	Simulation simulation;

	auto numArgs = args.size();
	if(numArgs != 1) {
		global_log->error() << "Incorrect number of arguments provided." << std::endl;
		op.print_usage();
		Simulation::exit(-1);
	}
	/* First read the given config file if it exists, then overwrite parameters with command line arguments. */
	std::string configFileName(args[0]);
	if( fileExists(configFileName.c_str()) ) {
		global_log->info() << "Config file: " << configFileName << endl;
		simulation.readConfigFile(configFileName);
	} else {
		global_log->error() << "Cannot open config file '" << configFileName << "'" << endl;
		Simulation::exit(-2);
	}

	/* processing command line arguments */
	if ( (int) options.get("legacy-cell-processor") > 0 ) {
		simulation.useLegacyCellProcessor();
		global_log->info() << "--legacy-cell-processor specified, using legacyCellProcessor" << endl;
	}

	if ( (int) options.get("final-checkpoint") > 0 ) {
		simulation.enableFinalCheckpoint();
		global_log->info() << "Final checkpoint enabled" << endl;
	} else {
		simulation.disableFinalCheckpoint();
		global_log->info() << "Final checkpoint disabled." << endl;
	}

	if( options.is_set_by_user("timed-checkpoint") ) {
		double checkpointtime = options.get("timed-checkpoint");
		simulation.setForcedCheckpointTime(checkpointtime);
		global_log->info() << "Enabling checkpoint after execution time: " << checkpointtime << " sec" << endl;
	}

	if (options.is_set_by_user("timesteps")) {
		simulation.setNumTimesteps(options.get("timesteps").operator unsigned long int());
	}
	if (options.is_set_by_user("loop-abort-time")) {
		simulation.setLoopAbortTime(options.get("loop-abort-time").operator double());
	}
	global_log->info() << "Simulating " << simulation.getNumTimesteps() << " steps." << endl;

	if(options.is_set_by_user("print-meminfo")) {
		global_log->info() << "Enabling memory info output" << endl;
		simulation.enableMemoryProfiler();
	}
	size_t lastIndex = configFileName.rfind(".");
	std::string outPrefix = configFileName.substr(0, lastIndex);
	if( options.is_set_by_user("outputprefix") ) {
		outPrefix = options["outputprefix"];
	}
	simulation.setOutputPrefix(outPrefix.c_str());
	global_log->info() << "Default output prefix: " << simulation.getOutputPrefix() << endl;


	simulation.prepare_start();

	Timer sim_timer;
	sim_timer.start();
	simulation.simulate();
	sim_timer.stop();
	double runtime = sim_timer.get_etime();
	//!@todo time only for simulation.simulate not "main"!
	global_log->info() << "main: used " << fixed << setprecision(2) << runtime << " seconds" << endl << fixed << setprecision(5);
	//  FIXME: The statements "<< fixed << setprecision(5)" after endl are so that the next logger timestamp appears as expected. A better solution would be nice, of course.

	// print out total simulation speed
	//!@todo is this correct w.r.t. starting from time > 0 ? We keep changing this...
	const unsigned long numForceCalculations = simulation.getNumTimesteps();
	const double speed = simulation.getTotalNumberOfMolecules() * numForceCalculations / runtime;
	global_log->info() << "Simulation speed: " << scientific << setprecision(6) << speed << " Molecule-updates per second." << endl << fixed << setprecision(5);

	const double iterationsPerSecond = simulation.getNumTimesteps() / runtime;
	global_log->info() << "Iterations per second: " << fixed << setprecision(3) << iterationsPerSecond << endl << fixed << setprecision(5);
	global_log->info() << "Time per iteration: " << fixed << setprecision(3) << 1.0 / iterationsPerSecond << " seconds." << endl << fixed << setprecision(5);

	simulation.finalize();

	delete global_log;

#ifdef ENABLE_MPI
	MPI_Finalize();
#endif
}
