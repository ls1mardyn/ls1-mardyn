#if ENABLE_MPI
#include <mpi.h>
#endif

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "WrapOpenMP.h"

#if ENABLE_MPI
#include "parallel/KDDecomposition.h"
#endif
#include "Simulation.h"
#include "utils/compile_info.h"
#include "utils/FileUtils.h"
#include "utils/Logger.h"
#include "utils/OptionParser.h"
#include "utils/Testing.h"
#include "utils/Timer.h"
#include "utils/SigsegvHandler.h"

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;


optparse::Values& initOptions(int argc, const char* const argv[], optparse::OptionParser& op);

/** Helper function outputting program build information to given logger */
void program_build_info(Log::Logger &log) {
	char info_str[MAX_INFO_STRING_LENGTH];
	get_compiler_info(info_str);
	log << "Compiler: " << info_str << endl;
	get_compile_time(info_str);
	log << "Compiled: " << info_str << endl;
#ifdef ENABLE_MPI
	get_mpi_info(info_str);
	log << "MPI library: " << info_str << endl;
#endif
}

/** Helper function outputting program invocation information to given logger */
void program_execution_info(int argc, char **argv, Log::Logger &log) {
	char info_str[MAX_INFO_STRING_LENGTH];
	get_timestamp(info_str);
	log << "Started: " << info_str << endl;
	get_host(info_str);
	log << "Execution host: " << info_str << endl;
	std::stringstream arguments;
	for (int i = 0; i < argc; i++) {
		arguments << " " << argv[i];
	}
	log << "Started with arguments: " << arguments.str() << endl;
#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
	global_log->info() << "Running with " << world_size << " MPI processes." << endl;
#endif
#if defined(_OPENMP)
	int num_threads = mardyn_get_max_threads();
	global_log->info() << "Running with " << num_threads << " OpenMP threads." << endl;
#endif
}

/** Run the internal unit tests */
int run_unit_tests(Values &options, vector<string> &args) {
	string testcases("");
	if(args.size() == 1) {
		testcases = args[0];
		global_log->info() << "Running unit tests: " << testcases << endl;
	} else {
		global_log->info() << "Running all unit tests!" << endl;
	}
	std::string testDataDirectory(options.get("testDataDirectory"));
	global_log->info() << "Test data directory: " << testDataDirectory << endl;
	Log::logLevel testLogLevel = options.is_set("verbose") && options.get("verbose") ? Log::All : Log::Info;
	int testresult = runTests(testLogLevel, testDataDirectory, testcases);
	return testresult;
}

/** @page main
 * In this project, software for molecular dynamics simulation with short-range
 * forces is developed. The aim is to have a parallel code (MPI) for
 * multi-centered molecules.
 *
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

	OptionParser op;
	Values options = initOptions(argc, argv, op);
	vector<string> args = op.args();

	if( options.is_set_by_user("logfile") ) {
		string logfileName(options.get("logfile"));
		global_log->info() << "Using logfile " << logfileName << endl;
		delete global_log;
		global_log = new Log::Logger(Log::Info, logfileName);
	}
	if( options.is_set_by_user("verbose") ) {
		global_log->info() << "Enabling verbose log output." << endl;
		global_log->set_log_level(Log::All);
	}
	if (options.is_set_by_user("sigsegvhandler")) {
		global_log->info() << "Enabling sigsegvhandler." << endl;
		registerSigsegvHandler();  // from SigsegvHandler.h
	}

	program_build_info(global_log->info());
	program_execution_info(argc, argv, global_log->info());

	if (options.is_set_by_user("tests")) {
		int testresult = run_unit_tests(options, args);
		#ifdef ENABLE_MPI
		MPI_Finalize();
		#endif
		exit(testresult); // using exit here should be OK
	}

	unsigned int numargs = args.size();
	if (numargs < 1) {
		op.print_usage();
		Simulation::exit(-13);
	}

	Simulation simulation;
	simulation.setName(op.prog());

	/** @todo remove unnamed options, present as --steps, --output-prefix below **/
	if( numargs > 2 ) {
		simulation.setOutputPrefix(args[2]);
	}

	/* First read the given config file if it exists, then overwrite parameters with command line arguments. */
	if( fileExists(args[0].c_str()) ) {
		global_log->info() << "Config file: " << args[0] << endl;
		simulation.readConfigFile(args[0]);
	} else {
		global_log->error() << "Cannot open input file '" << args[0] << "'" << endl;
		exit(-54); // Simulation::exit(-54);
	}

	/** @todo remove unnamed options, present as --steps, --output-prefix below **/
	if (numargs > 1) {
		unsigned long steps = 0;
		istringstream(args[1]) >> steps;
		simulation.setNumTimesteps(steps);
	}


	if ( (int) options.get("final-checkpoint") > 0 ) {
		simulation.enableFinalCheckpoint();
		global_log->info() << "Final checkpoint enabled" << endl;
	} else {
		simulation.disableFinalCheckpoint();
		global_log->info() << "Final checkpoint disabled." << endl;
	}

	double time = options.get("timed-checkpoint");
	simulation.setForcedCheckpointTime(time);
	if( options.is_set_by_user("timed-checkpoint") ) {
		global_log->info() << "Enabling checkpoint after execution time: " << time << " sec" << endl;
	}

	if (options.is_set_by_user("timesteps")) {
		simulation.setNumTimesteps(options.get("timesteps").operator unsigned long int());
	}
	global_log->info() << "Simulating " << simulation.getNumTimesteps() << " steps." << endl;
    
	string outprefix(args[numargs-1]);
	size_t foundstr = outprefix.rfind(".");
	if (foundstr!=string::npos && foundstr==outprefix.size()-4) outprefix.erase(foundstr);	// remove .??? suffix
	simulation.setOutputPrefix(outprefix.c_str());
	if( options.is_set_by_user("outputprefix") ) {
		simulation.setOutputPrefix( options["outputprefix"] );
	}
	global_log->info() << "Default output prefix: " << simulation.getOutputPrefix() << endl;


	simulation.prepare_start();

	Timer sim_timer;
	sim_timer.start();
	simulation.simulate();
	sim_timer.stop();
	double runtime = sim_timer.get_etime();
	global_log->info() << "main: used " << fixed << setprecision(2) << runtime << " seconds" << endl;

	// print out total simulation speed
	const unsigned long numForceCalculations = simulation.getNumTimesteps() + 1ul;
	const double speed = simulation.getTotalNumberOfMolecules() * numForceCalculations / runtime;
	global_log->info() << "Simulation speed: " << scientific << speed << " Molecule-updates per second." << endl;

	simulation.finalize();

	delete global_log;

#ifdef ENABLE_MPI
	MPI_Finalize();
#endif
}


Values& initOptions(int argc, const char* const argv[], OptionParser& op) {
	op = OptionParser()
		.usage("%prog [<scenario generator with options> | <configfilename>] <number of timesteps> <outputprefix>\n "
		  "      %prog --tests --test-dir <test input data directory> [<name of testcase>]\n\n"
				"Use option --help to display all available options.")
		.version("%prog 1.1")
		.description("ls1 mardyn (M-olecul-AR DYN-amics)");

	op.add_option("-n", "--steps") .dest("timesteps") .metavar("NUM") .type("int") .set_default(1) .help("number of timesteps to simulate (default: %default)");
	op.add_option("-p", "--outprefix") .dest("outputprefix") .metavar("STR") .type("string") .set_default("MarDyn") .help("default prefix for output files (default: %default)");
	op.add_option("-v", "--verbose") .action("store_true") .dest("verbose") .metavar("V") .type("bool") .set_default(false) .help("verbose mode: print debugging information (default: %default)");
	op.add_option("-S", "--sigsegvhandler") .action("store_true") .dest("sigsegvhandler") .metavar("S") .type("bool") .set_default(false) .help("sigsegvhandler: prints stacktrace on sigsegv(default: %default)");
	op.add_option("--logfile").dest("logfile").type("string").set_default("MarDyn.log").metavar("STRING").help("enable/disable final checkopint (default: %default)");
	op.add_option("--final-checkpoint").dest("final-checkpoint").type("int").set_default(1).metavar("(1|0)").help("enable/disable final checkopint (default: %default)");
	op.add_option("--timed-checkpoint").dest("timed-checkpoint").type("float").set_default(-1).help("Execution time of the simulation in seconds after which a checkpoint is forced.");

    op.add_option("-t", "--tests").action("store_true").dest("tests").metavar("T").type("bool").set_default(false).help("unit tests: run built-in unit tests (default: %default)");
    op.add_option("-d", "--test-dir").dest("testDataDirectory") .metavar("STR") .set_default("") .help("unit tests: specify the directory where the in input data required by the tests resides");

	//OptionGroup dgroup = OptionGroup(op, "Developer options", "Advanced options for developers and experienced users.");
	//dgroup.add_option("--phasespace-file") .metavar("FILE") .help("path to file containing phase space data");
	//char const* const pc_choices[] = { "LinkedCells" };
	//dgroup.add_option("--particle-container") .choices(&pc_choices[0], &pc_choices[1]) .set_default(pc_choices[0]) .help("container used for locating nearby particles (default: %default)");
	//op.add_option_group(dgroup);

	return op.parse_args(argc, argv);
}
