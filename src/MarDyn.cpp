#include <iostream>
#include <iomanip>
#include <ctime>

#include "utils/OptionParser.h"
#include "utils/Logger.h"
#include "utils/compile_info.h"
#include "utils/Testing.h"
#include "utils/FileUtils.h"
#include "Simulation.h"


using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;


optparse::Values& initOptions(int argc, const char* const argv[], optparse::OptionParser& op);


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

#ifdef ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
	/* Initialize the global log file */
	//string logfileName("MarDyn");
	//global_log = new Log::Logger(Log::All, logfileName);
	global_log = new Log::Logger(Log::Info);
	cout.precision(6);
#ifdef ENABLE_MPI
	global_log->set_mpi_output_root(0);
#endif

	/* Print some info about the program itself */
	char info_str[MAX_INFO_STRING_LENGTH];
	get_compiler_info(info_str);
	global_log->info() << "Compiler: " << info_str << endl;
	get_compile_time(info_str);
	global_log->info() << "Compiled: " << info_str << endl;
#ifdef ENABLE_MPI
	get_mpi_info(info_str);
	global_log->info() << "MPI library: " << info_str << endl;
#endif
	get_timestamp(info_str);
	global_log->info() << "Started: " << info_str << endl;
	get_host(info_str);
	global_log->info() << "Execution host: " << info_str << endl;

#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK( MPI_Comm_size( MPI_COMM_WORLD, &world_size ) );
	global_log->info() << "Running with " << world_size << " processes." << endl;
#endif

	OptionParser op;
	Values options = initOptions(argc, argv, op);
	vector<string> args = op.args();
	unsigned int numargs = args.size();
	if (numargs < 1) {
		op.print_usage();
		exit(1);
	}

	if (options.is_set("verbose") && options.get("verbose"))
		global_log->set_log_level(Log::All);

	bool tests(options.get("tests"));
	if (tests) {
		string testcases = "";
		if (numargs == 1 ) {
			testcases = args[0];
			global_log->info() << "Running unit tests: " << testcases << endl;
		} else {
			global_log->info() << "Running all unit tests! (Too many arguments were specified)" << endl;
		}

		std::string testDataDirectory(options.get("testDataDirectory"));
		Log::logLevel testLogLevel = options.is_set("verbose") && options.get("verbose") ? Log::All : Log::Info;
		bool testresult = runTests(testLogLevel, testDataDirectory, testcases);
		if (testresult) {
			#ifdef ENABLE_MPI
			MPI_Finalize();
			#endif
			exit(1);
		} else {
			#ifdef ENABLE_MPI
			MPI_Finalize();
			#endif
			exit(0);
		}
	}

    Simulation simulation;

    /* read the given config file if it exists */
    if( fileExists(args[0].c_str()) ) {
        string inputfilename(args[0]);
        simulation.readConfigFile(inputfilename);
    } else {
        exit(1);
    }

    // set the number of timesteps to be simulated
    if( numargs > 1 ) {
        unsigned long steps = 0;
        istringstream(args[1]) >> steps;
        simulation.setNumTimesteps( steps );
    }
    if (options.is_set_by_user("timesteps")) {
        simulation.setNumTimesteps(options.get("timesteps"));
    }
    global_log->info() << "Simulating " << simulation.getNumTimesteps() << " steps." << endl;

    // set the prefix for output files
    if( numargs > 2 ) {
        simulation.setOutputPrefix( args[2] );
    }
    if( options.is_set_by_user("outputprefix") ) {
        simulation.setOutputPrefix( options["outputprefix"] );
    }
    global_log->info() << "Using output prefix '" << simulation.getOutputPrefix() << "'" << endl;

    // TODO
    //simulation.setCutoffRadius(options.get("cutoff_radius"));


    simulation.prepare_start();

    double runtime = double(clock()) / CLOCKS_PER_SEC;
    simulation.simulate();
    runtime = double(clock()) / CLOCKS_PER_SEC - runtime;
    cout << "main: used " << fixed << setprecision(2) << runtime << " s" << endl;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
}


Values& initOptions(int argc, const char* const argv[], OptionParser& op) {

	op = OptionParser()
		// The last two optional positional arguments are only here for backwards-compatibility
		.usage("%prog [-n steps] [-p prefix] <configfilename> [<number of timesteps>] [<outputprefix>]\n "
		  "      %prog -t -d <test input data directory> [<name of testcase>]\n\n"
				"Use option --help to display all available options.")
		.version("%prog 1.0")
		.description("MarDyn is a MD simulator. All behavior is controlled via the config file.")
		// .epilog("background info?")
		;

	op.add_option("-n", "--steps") .dest("timesteps") .metavar("NUM") .type("int") .set_default(1) .help("number of timesteps to simulate (default: %default)");
	op.add_option("-p", "--outprefix") .dest("outputprefix") .metavar("STR") .help("prefix for output files");
	op.add_option("-v", "--verbose") .action("store_true") .dest("verbose") .metavar("V") .type("bool") .set_default(false) .help("verbose mode: print debugging information (default: %default)");
	op.add_option("-t", "--tests") .action("store_true") .dest("tests") .metavar("T") .type("bool") .set_default(false) .help("unit tests: run built-in unit tests (default: %default)");
	op.add_option("-d", "--test-dir").dest("testDataDirectory") .metavar("STR") .set_default("") .help("unit tests: specify the directory where the in input data required by the tests resides");

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


