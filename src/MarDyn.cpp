#include "MarDyn_version.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

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
#ifdef MARDYN_AUTOPAS
#include "autopas/Version.h"
#endif

#include "parallel/CollectiveCommunicationPersistent.h"


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
void log_program_build_info() {
	Log::global_log->info() << "Compilation info:" << std::endl;

	char info_str[MAX_INFO_STRING_LENGTH];
	get_compiler_info(info_str);
	Log::global_log->info() << "	Compiler:	" << info_str << std::endl;
	get_compile_time(info_str);
	Log::global_log->info() << "	Compiled on:	" << info_str << std::endl;
	get_precision_info(info_str);
	Log::global_log->info() << "	Precision:	" << info_str << std::endl;
	get_intrinsics_info(info_str);
	Log::global_log->info() << "	Intrinsics:	" << info_str << std::endl;
	get_rmm_normal_info(info_str);
	Log::global_log->info() << "	RMM/normal:	" << info_str << std::endl;
	get_openmp_info(info_str);
	Log::global_log->info() << "	OpenMP:		" << info_str << std::endl;
	get_mpi_info(info_str);
	Log::global_log->info() << "	MPI:		" << info_str << std::endl;
}

/**
 * @brief Helper function outputting program invocation information to given logger
 */
void log_program_execution_info(int argc, char **argv) {
	Log::global_log->info() << "Execution info:" << std::endl;

	char info_str[MAX_INFO_STRING_LENGTH];
	get_timestamp(info_str);
	Log::global_log->info() << "	Started: " << info_str << std::endl;
	get_host(info_str);
	Log::global_log->info() << "	Execution host: " << info_str << std::endl;
	std::stringstream arguments;
	for (int i = 0; i < argc; i++) {
		arguments << " " << argv[i];
	}
	Log::global_log->info() << "	Started with arguments: " << arguments.str() << std::endl;

#if defined(_OPENMP)
	int num_threads = mardyn_get_max_threads();
	Log::global_log->info() << "	Running with " << num_threads << " OpenMP threads." << std::endl;
	// print thread pinning info
	PrintThreadPinningToCPU();
#endif

#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
	Log::global_log->info() << "	Running with " << world_size << " MPI processes." << std::endl;
#endif
}

/** Run the internal unit tests */
int run_unit_tests(const Values &options, const std::vector<std::string> &args) {
	std::string testcases("");
	if(args.size() == 1) {
		testcases = args[0];
		Log::global_log->info() << "Running unit tests: " << testcases << std::endl;
	} else {
		Log::global_log->info() << "Running all unit tests!" << std::endl;
	}
	std::string testDataDirectory(options.get("testDataDirectory"));
	Log::global_log->info() << "Test data directory: " << testDataDirectory << std::endl;
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
	MPI_Env_Wrapper::init_environment(&argc, &argv);
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
	Log::global_log = std::make_unique<Log::Logger>(Log::Info);

	// Open scope to exclude MPI_Init() and MPI_Finalize().
	// This way, all simulation objects are cleaned up before MPI finalizes.
	{

	optparse::OptionParser op;
	initOptions(&op);
	optparse::Values options = op.parse_args(argc, argv);
	std::vector<std::string> args = op.args();

	/* Initialize the global log file */
	if( options.is_set_by_user("logfile") ) {
		// Print to file
		std::string logfileNamePrefix(options.get("logfile"));
		std::string logfile_name = logfileNamePrefix;
#ifdef ENABLE_MPI
		logfile_name += "_R" + std::to_string(world_rank);
#endif
		logfile_name += ".log";
		Log::global_log->info() << "Using logfile " << logfile_name << " for further log output." << std::endl;
		std::shared_ptr<std::ostream> logfile_ptr = std::make_shared<std::ofstream>(logfile_name);
		Log::global_log->set_log_stream(logfile_ptr);
	}

#ifdef ENABLE_MPI
	Log::global_log->set_mpi_output_root(0);
	//global_log->set_mpi_output_all();
#endif

	Log::global_log->info() << "Running ls1-MarDyn version " << MARDYN_VERSION << std::endl;

#ifdef MARDYN_AUTOPAS
	Log::global_log->info() << "Built with AutoPas version " << AutoPas_VERSION << std::endl;
#endif

#ifndef NDEBUG
	Log::global_log->warning() << "This ls1-MarDyn binary is a DEBUG build!" << std::endl;
#endif

	if( options.is_set_by_user("verbose") ) {
		Log::global_log->info() << "Enabling verbose log output." << std::endl;
		Log::global_log->set_log_level(Log::All);
	}
#ifdef ENABLE_SIGHANDLER
	if (options.is_set_by_user("sigsegvhandler")) {
		Log::global_log->info() << "Enabling sigsegvhandler." << std::endl;
		registerSigsegvHandler();  // from SigsegvHandler.h
	}
#endif
	log_program_build_info();
	log_program_execution_info(argc, argv);


	/* Run built in tests and exit */
	if (options.is_set_by_user("tests")) {
		int testresult = run_unit_tests(options, args);
		std::exit(testresult); // using exit here should be OK
	}


	/* Set up and run regular Simulation */
	Simulation simulation;

	auto numArgs = args.size();
	if(numArgs != 1) {
		op.print_usage();
		std::ostringstream error_message;
		error_message << "Incorrect number of arguments provided." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	/* First read the given config file if it exists, then overwrite parameters with command line arguments. */
	std::string configFileName(args[0]);
	if( fileExists(configFileName.c_str()) ) {
		Log::global_log->info() << "Config file: " << configFileName << std::endl;
		simulation.readConfigFile(configFileName);
	} else {
		std::ostringstream error_message;
		error_message << "Cannot open config file '" << configFileName << "'" << std::endl;
		MARDYN_EXIT(error_message.str());
	}

	/* processing command line arguments */
	if ( (int) options.get("legacy-cell-processor") > 0 ) {
		simulation.useLegacyCellProcessor();
		Log::global_log->info() << "--legacy-cell-processor specified, using legacyCellProcessor" << std::endl;
	}

	if ( (int) options.get("final-checkpoint") > 0 ) {
		simulation.enableFinalCheckpoint();
		Log::global_log->info() << "Final checkpoint enabled" << std::endl;
	} else {
		simulation.disableFinalCheckpoint();
		Log::global_log->info() << "Final checkpoint disabled." << std::endl;
	}

	if( options.is_set_by_user("timed-checkpoint") ) {
		double checkpointtime = options.get("timed-checkpoint");
		simulation.setForcedCheckpointTime(checkpointtime);
		Log::global_log->info() << "Enabling checkpoint after execution time: " << checkpointtime << " sec" << std::endl;
	}

	if (options.is_set_by_user("timesteps")) {
		simulation.setNumTimesteps(options.get("timesteps").operator unsigned long int());
	}
	if (options.is_set_by_user("loop-abort-time")) {
		simulation.setLoopAbortTime(options.get("loop-abort-time").operator double());
	}
	Log::global_log->info() << "Simulating " << simulation.getNumTimesteps() << " steps." << std::endl;

	if(options.is_set_by_user("print-meminfo")) {
		Log::global_log->info() << "Enabling memory info output" << std::endl;
		simulation.enableMemoryProfiler();
	}
	size_t lastIndex = configFileName.rfind(".");
	std::string outPrefix = configFileName.substr(0, lastIndex);
	if( options.is_set_by_user("outputprefix") ) {
		outPrefix = options["outputprefix"];
	}
	simulation.setOutputPrefix(outPrefix.c_str());
	Log::global_log->info() << "Default output prefix: " << simulation.getOutputPrefix() << std::endl;


	simulation.prepare_start();

	Timer sim_timer;
	sim_timer.start();
	simulation.simulate();
	sim_timer.stop();
	double runtime = sim_timer.get_etime();
	//!@todo time only for simulation.simulate not "main"!
	Log::global_log->info() << "main: used " << std::fixed << std::setprecision(2) << runtime << " seconds" << std::endl << std::fixed << std::setprecision(5);
	//  FIXME: The statements "<< std::fixed << std::setprecision(5)" after endl are so that the next logger timestamp appears as expected. A better solution would be nice, of course.

	// print out total simulation speed
	const unsigned long numTimesteps = simulation.getNumTimesteps() - simulation.getNumInitTimesteps();
	const double speed = simulation.getTotalNumberOfMolecules() * numTimesteps / runtime;
	Log::global_log->info() << "Simulation speed: " << std::scientific << std::setprecision(6) << speed << " Molecule-updates per second." << std::endl << std::fixed << std::setprecision(5);

	const double iterationsPerSecond = numTimesteps / runtime;
	Log::global_log->info() << "Iterations per second: " << std::fixed << std::setprecision(3) << iterationsPerSecond << std::endl << std::fixed << std::setprecision(5);
	Log::global_log->info() << "Time per iteration: " << std::fixed << std::setprecision(3) << 1.0 / iterationsPerSecond << " seconds." << std::endl << std::fixed << std::setprecision(5);

	double resources = runtime / 3600.0;
#if defined(_OPENMP)
	resources *= mardyn_get_max_threads();
#endif

#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
	resources *= world_size;
#endif
	Log::global_log->info() << "Used resources: " << std::fixed << std::setprecision(3) << resources << " core-hours" << std::endl << std::fixed << std::setprecision(5);

	simulation.finalize();

	} // End of scope to exclude MPI_Init() and MPI_Finalize()
}
