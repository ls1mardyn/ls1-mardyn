#include <iostream>
#include <iomanip>
#include <ctime>

#include "utils/OptionParser.h"
#include "utils/Logger.h"
#include "utils/compile_info.h"
#include "utils/Testing.h"
#include "utils/FileUtils.h"
#include "utils/Timer.h"
#include "Simulation.h"

#if ENABLE_MPI
#include "parallel/KDDecomposition2.h"
#include <mpi.h>
#endif

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

    OptionParser op;
    Values options = initOptions(argc, argv, op);
    vector<string> args = op.args();
    unsigned int numargs = args.size();

    if( options.is_set_by_user("verbose") ) {
        global_log->info() << "Enabling verbose log output." << endl;
        global_log->set_log_level(Log::All);
    }

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

	std::stringstream arguments;
	for (int i = 0; i < argc; i++) {
		arguments << " " << argv[i];
	}
	global_log->info() << "Started with arguments: " << arguments.str() << endl;

#ifdef ENABLE_MPI
	int world_size = 1;
	MPI_CHECK( MPI_Comm_size( MPI_COMM_WORLD, &world_size ) );
	global_log->info() << "Running with " << world_size << " processes." << endl;
#endif

	bool tests(options.is_set_by_user("tests"));
	if (tests) {
		string testcases;
		if (numargs == 1 ) {
			testcases = args[0];
			global_log->info() << "Running unit tests: " << testcases << endl;
		} else {
			global_log->info() << "Running all unit tests!" << endl;
		}

		std::string testDataDirectory(options.get("testDataDirectory"));
		Log::logLevel testLogLevel = options.is_set("verbose") && options.get("verbose") ? Log::All : Log::Info;
		bool testresult = runTests(testLogLevel, testDataDirectory, testcases);

		#ifdef ENABLE_MPI
		MPI_Finalize();
		#endif
		if (testresult) {
			exit(1);
		} else {
			exit(0);
		}
	}

	if (numargs < 1) {
		op.print_usage();
		exit(1);
	}

	Simulation simulation;

	/* First read the given config file if it exists, then overwrite parameters with command line arguments. */
	if(args[0] == "mkTcTS") {
		simulation.mkTcTS(options);
		// set the number of timesteps to be simulated
		if (numargs > 2) {
			unsigned long steps = 0;
			istringstream(args[numargs - 2]) >> steps;
			simulation.setNumTimesteps(steps);
		}
		if( numargs > 3 ) {
			simulation.setOutputPrefix( args[numargs - 1] );
		}
	}
	else if(args[0] == "mkesfera") {
		simulation.mkesfera(options);
		// set the number of timesteps to be simulated
		if (numargs > 2) {
			unsigned long steps = 0;
			istringstream(args[numargs - 2]) >> steps;
			simulation.setNumTimesteps(steps);
		}
		if( numargs > 3 ) {
			simulation.setOutputPrefix( args[numargs - 1] );
		}
	}
	else if( fileExists( args[0].c_str()) ) {
		if (numargs > 1) {
			unsigned long steps = 0;
			istringstream(args[1]) >> steps;
			simulation.setNumTimesteps(steps);
		}
		if( numargs > 2 ) {
			simulation.setOutputPrefix( args[2] );
		}
			simulation.readConfigFile( args[0] );
	} else {
		global_log->error() << "Cannot open input file '" << args[0] << "'" << std::endl;
		exit(1);
	}

	if ( (int) options.get("final-checkpoint") > 0 ) {
		simulation.enableFinalCheckpoint();
		global_log->info() << "Final checkpoint enabled" << endl;
	}
	else {
		simulation.disableFinalCheckpoint();
		global_log->info() << "Final checkpoint disbaled." << endl;
	}

	if( options.is_set_by_user("timed-checkpoint") ) {
		double time = options.get("timed-checkpoint");
		global_log->info() << "Enabling checkpoint after execution time: " << time << " sec" << endl;
		simulation.setForcedCheckpointTime(time);
	}

    if (options.is_set_by_user("timesteps")) {
        simulation.setNumTimesteps(options.get("timesteps"));
    }
#if ENABLE_MPI
    if( options.is_set("domain-decomposition") ) {
        string domainDecomp((const char *) options.get("domain-decomposition"));
        global_log->info() << "Setting domain decomposition via command line to " << domainDecomp << endl;
        if ( domainDecomp == "KDDecomposition2" ) {
            DomainDecompBase *domainDecomp;
            int updateFrequency = 100;
            int fullSearchThreshold = 3;
            domainDecomp = (DomainDecompBase*) new KDDecomposition2(simulation.getcutoffRadius(), simulation.getDomain(), updateFrequency, fullSearchThreshold);
            simulation.setDomainDecomposition(domainDecomp);
        }
    }
#endif
    
    global_log->info() << "Simulating " << simulation.getNumTimesteps() << " steps." << endl;

    // set the prefix for output files
    simulation.setOutputPrefix(args[numargs - 1]);
    /*
    if( options.is_set_by_user("outputprefix") ) {
        simulation.setOutputPrefix( options["outputprefix"] );
    }
    */
    global_log->info() << "Using output prefix '" << simulation.getOutputPrefix() << "'" << endl;

    // TODO
    //simulation.setCutoffRadius(options.get("cutoff_radius"));

    simulation.prepare_start();

    Timer sim_timer;
    sim_timer.start();
    simulation.simulate();
    sim_timer.stop();
    double runtime = sim_timer.get_etime();
    global_log->info() << "main: used " << fixed << setprecision(2) << runtime << " s" << endl;

	simulation.finalize();

    delete global_log;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
}


Values& initOptions(int argc, const char* const argv[], OptionParser& op) {

        /*
	op = OptionParser()
		.usage("%prog [-n steps] [-p prefix] <configfilename> [<number of timesteps>] [<outputprefix>]\n "
		  "      %prog -t -d <test input data directory> [<name of testcase>]\n\n"
				"Use option --help to display all available options.")
		.version("%prog 1.0")
		.description("MarDyn is a MD simulator. All behavior is controlled via the config file.")
		// .epilog("background info?")
		;
        */
	op = OptionParser()
		.usage("%prog [<scenario generator with options> | <configfilename>] <number of timesteps> <outputprefix>\n "
		  "      %prog -t -d <test input data directory> [<name of testcase>]\n\n"
				"Use option --help to display all available options.")
		.version("%prog 0.1")
		.description("ls1 mardyn (M-olecul-AR DYN-amics)")
		// .epilog("background info?")
		;

	op.add_option("-n", "--steps") .dest("timesteps") .metavar("NUM") .type("int") .set_default(1) .help("number of timesteps to simulate (default: %default)");
	// op.add_option("-p", "--outprefix") .dest("outputprefix") .metavar("STR") .help("prefix for output files");
	op.add_option("-v", "--verbose") .action("store_true") .dest("verbose") .metavar("V") .type("bool") .set_default(false) .help("verbose mode: print debugging information (default: %default)");
    op.add_option("--final-checkpoint").dest("final-checkpoint").type("int").set_default(1).metavar("(1|0)").help("enable/disable final checkopint (default: %default)");
	op.add_option("--timed-checkpoint").dest("timed-checkpoint").type("float").help("Execution time of the simulation in seconds after which a checkpoint is forced.");
    
    OptionGroup mkTcTS_options = OptionGroup(op, "mkTcTS options", "Options for the mkTcTS scenario generator.");
    mkTcTS_options.add_option("-c").dest("density-1").type("float").help("density 1 (mkTcTS, mkesfera)");
    mkTcTS_options.add_option("-C").dest("density-2").type("float").help("density 2 (mkTcTS, mkesfera)");
    mkTcTS_options.add_option("-N").dest("num-particles").type("int").help("approximate number of fluid molecules (mkTcTS)");
    mkTcTS_options.add_option("-T").dest("temperature").type("float").help("temperature (mkTcTS, mkesfera)");
    
    mkTcTS_options.add_option("-R").dest("cutoff-LJ").type("float").help("Lennard-Jones cutoff (mkTcTS, mkesfera)");
    mkTcTS_options.add_option("-S").dest("shift_LJ").type("bool").help("shift the LJ potential (mkTcTS, mkesfera)").set_default(false);
    mkTcTS_options.add_option("-h").dest("height").type("float").help("height (mkTcTS)");
    mkTcTS_options.add_option("-m").dest("chemical-potential").type("float").help("chemical potential (mkTcTS)");
    mkTcTS_options.add_option("-p").dest("pair-correlation-cutoff").type("float").help("pair correlation cutoff (mkTcTS)");
    op.add_option_group(mkTcTS_options);

	OptionGroup mkesfera_options = OptionGroup(op, "mkesfera options", "Options for the mkesfera scenario generator.");
	mkesfera_options.add_option("-I").dest("droplet-radius-inner").type("float").help("inner droplet radius (mkesfera)");
	mkesfera_options.add_option("-O").dest("droplet-radius-outer").type("float").help("outer droplet radius (mkesfera)");
	op.add_option_group(mkesfera_options);
    
    op.add_option("-t", "--tests").action("store_true").dest("tests").metavar("T").type("bool").set_default(false).help("unit tests: run built-in unit tests (default: %default)");
    op.add_option("-d", "--test-dir").dest("testDataDirectory") .metavar("STR") .set_default("") .help("unit tests: specify the directory where the in input data required by the tests resides");
    
	OptionGroup dgroup = OptionGroup(op, "Developer options", "Advanced options for developers and experienced users.");
	dgroup.add_option("--phasespace-file") .metavar("FILE") .help("path to file containing phase space data");
	char const* const pc_choices[] = { "LinkedCells", "AdaptiveSubCells" };
	dgroup.add_option("--particle-container") .choices(&pc_choices[0], &pc_choices[2]) .set_default(pc_choices[0]) .help("container used for locating nearby particles (default: %default)");
	dgroup.add_option("--cutoff-radius") .type("float") .set_default(5.0) .help("radius of sphere around a particle in which forces are considered (default: %default)");
	dgroup.add_option("--cells-in-cutoff") .type("int") .set_default(2) .help("number of cells in cutoff-radius cube (default: %default); only used by LinkedCells particle container");
	char const* const dd_choices[] = { "DomainDecomposition", "KDDecomposition", "KDDecomposition2" };
	dgroup.add_option("--domain-decomposition").dest("domain-decomposition") .choices(&dd_choices[0], &dd_choices[3]) .set_default(dd_choices[0]) .help("domain decomposition strategy for MPI (default: %default)");
	dgroup.add_option("--timestep-length") .type("float") .set_default(0.004) .help("length of one timestep in TODO (default: %default)");
	op.add_option_group(dgroup);

// 	OptionGroup animake_options = OptionGroup(op, "animake options", "Options for the animake scenario generator.");
// 	animake_options.add_option("-c").dest("fluid-density").help("density (animake)");
// 	animake_options.add_option("-f").dest("fluid").help("fluid type (animake) possible values: CH4 (default), Ar, C2H6, N2, CO2, EOX, JES or VEG");
// 	animake_options.add_option("-g").dest("fluid2").help("second fluid type (animake) possible values: CH4 (default), Ar, C2H6, N2, CO2, EOX, JES or VEG");
// 	animake_options.add_option("-N").dest("num-particles").help("number of molecules (animake)");
// 	animake_options.add_option("-T").dest("temperature").help("temperature (animake)");
// 	animake_options.add_option("-y").dest("length_y").help("box length_y (animake)");
// 	animake_options.add_option("-z").dest("length_z").help("box length_z (animake)");
// 	animake_options.add_option("-Y").dest("reference-mass").help("reference mass [u] (animake)").set_default(1.0);
// 	animake_options.add_option("-s").dest("reference-length").help("reference length [A] (animake)").set_default(1.0);
// 	animake_options.add_option("-W").dest("reference-energy").help("reference energy and temperature [K] (animake)").set_default(1.0);
// 	animake_options.add_option("-x").dest("mole-fraction").help("mole fraction fluid1:fluid2 (animake)");
// 	animake_options.add_option("-m").dest("mu").help("chemical potential (animake)").set_default(0.0);
// 	animake_options.add_option("-J").dest("mixing-eta").help("Lorenz Berthelot mixing parameter 'eta' (animake)").set_default(1.0);
// 	animake_options.add_option("-5").dest("mixing-xi").help("Lorenz Berthelot mixing parameter 'xi' (animake)").set_default(1.0);
// 	op.add_option_group(animake_options);
	
	return op.parse_args(argc, argv);
}


