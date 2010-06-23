/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// Simulation.cpp
#define SIMULATION_SRC
#include "Simulation.h"

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "parallel/DomainDecompBase.h"

#ifdef PARALLEL
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"

#include "utils/xmlfile.h"
#include "io/PartGen.h"
#include "io/ResultWriter.h"
#include "io/XyzWriter.h"
#include "io/PovWriter.h"
#include "io/DecompWriter.h"
#include "io/CheckpointWriter.h"
#include "io/VISWriter.h"
#include "io/InputOldstyle.h"

#include "ensemble/GrandCanonical.h"

#ifdef STEEREO
#include "commands/snapshotCommand.h"
#include "commands/megaMolSnapshotCommand.h"
#include "commands/sendCouplingMDCommand.h"
#include "commands/receiveCouplingMDCommand.h"
#include "commands/getVisDataCommand.h"
#ifdef PARALLEL
#include <steereoMPIIntraCommunicator.h>
#endif //PARALLEL
#include <steerParameterCommand.h>
#include <steereoSocketCommunicator.h>
#endif

#include "utils/OptionParser.h"
#include "utils/Timer.h"
#include "utils/Logger.h"
#include "utils/compile_info.h"

#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

Simulation* global_simulation;

Simulation::Simulation(int *argc, char ***argv) {
#ifdef PARALLEL
  MPI_Init(argc, argv);
#endif

  global_simulation = this;

  /* Initialize the global log file */
  //string logfileName("MarDyn");
  //global_log = new Log::Logger(Log::All, logfileName);
  global_log = new Log::Logger(Log::Info);
#ifdef PARALLEL
  global_log->set_mpi_output_root(0);
#endif

  OptionParser op;
  Values options = initOptions(*argc, *argv, op);
  vector<string> args = op.args();
  unsigned int numargs = args.size();

  if (numargs < 1) {
    op.print_usage();
    exit(1);
  }
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
  int world_size = 1;
#ifdef PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
#endif
  global_log->info() << "Running with " << world_size << " processes." << endl;

#ifndef PARALLEL
  global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
  _domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
  global_log->info() << "Initialization done" << endl;
#endif

  /*
   * default parameters
   */
  _cutoffRadius = 0.0;
  _LJCutoffRadius = 0.0;
  _numberOfTimesteps = 1;
  _outputPrefix = string("mardyn");
  _outputPrefix.append(gettimestring());
  _resultOutputTimesteps = 25;
  _doRecordProfile = false;
  _profileRecordingTimesteps = 7;
  _profileOutputTimesteps = 12500;
  _profileOutputPrefix = "out";
  _collectThermostatDirectedVelocity = 100;
  _zoscillation = false;
  _zoscillator = 512;
  _doRecordRDF = false;
  _RDFOutputTimesteps = 25000;
  _RDFOutputPrefix = "out";
  _initCanonical = 5000;
  _initGrandCanonical = 10000000;
  _initStatistics = 20000;
  h = 0.0;

  _cutoffRadius = options.get("cutoff_radius");

  // store number of timesteps to be simulated
  if (numargs >= 2)
    istringstream(args[1]) >> _numberOfTimesteps;
  if (options.is_set_by_user("timesteps"))
    _numberOfTimesteps = options.get("timesteps");
  global_log->info() << "Simulating " << _numberOfTimesteps << " steps." << endl;

  // store prefix for output files
  if (numargs >= 3)
    _outputPrefix = args[2];
  if (options.is_set_by_user("outputprefix"))
    _outputPrefix = options["outputprefix"];

  if (options.is_set("verbose") && options.get("verbose"))
    global_log->set_log_level(Log::All);

  string inputfilename(args[0]);
  initConfigFile(inputfilename);
}

Simulation::~Simulation() {
#ifdef PARALLEL
  MPI_Finalize();
#endif
}

void Simulation::exit(int exitcode) {
#ifdef PARALLEL
  // terminate all mpi processes and return exitcode
  MPI_Abort( MPI_COMM_WORLD, exitcode );
#else
  // call global exit
  ::exit(exitcode);
#endif
}

const Values& Simulation::initOptions(int argc, char *argv[], OptionParser& op) {

  op = OptionParser()
      // The last two optional positional arguments are only here for backwards-compatibility
      .usage("%prog [-n steps] [-p prefix] <configfilename> [<number of timesteps>] [<outputprefix>]")
      .version("%prog 1.0")
      .description("MarDyn is a MD simulator. All behavior is controlled via the config file.")
      // .epilog("background info?")
  ;

  op.add_option("-n", "--steps") .dest("timesteps") .metavar("NUM") .type("int") .set_default(1) .help("number of timesteps to simulate (default: %default)");
  op.add_option("-p", "--outprefix") .dest("outputprefix") .metavar("STR") .help("prefix for output files");
  op.add_option("-v", "--verbose") .dest("verbose") .metavar("V") .type("bool") .set_default(false) .help("verbose mode: print debugging information (default: %default)");

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

void Simulation::initConfigFile(const string& inputfilename) {
  if (inputfilename.rfind(".xml") == inputfilename.size() - 4) {
    global_log->info() << "command line config file type is XML (*.xml)" << endl;
    initConfigXML(inputfilename);
  }
  else if (inputfilename.rfind(".cfg") == inputfilename.size() - 4) {
    global_log->info() << "command line config file type is oldstyle (*.cfg)" << endl;
    initConfigOldstyle(inputfilename);
  }
  else {
    global_log->info() << "command line config file type is unknown: trying oldstyle" << endl;
    initConfigOldstyle(inputfilename);
  }
}

void Simulation::initConfigXML(const string& inputfilename) {
  global_log->info() << "init XML config file: " << inputfilename << endl;
  XMLfile inp(inputfilename);
  //inp.printXML();
  global_log->debug() << string(inp) << endl;
  if (!inp.changecurrentnode("/mardyn")) {
    global_log->error() << inputfilename << " is not a MarDyn XML input file!" << endl;
    return;
  }
  string version;
  //inp.getNodeValue("/mardyn@version",version);
  inp.getNodeValue("@version", version);
  global_log->debug() << "MarDyn XML config file version " << version << endl;
  if (inp.changecurrentnode("simulation")) {
    string siminpfile, siminptype;
    if (inp.getNodeValue("input", siminpfile)) {
      // input type="oldstyle" to include oldstyle input files for backward compatibility - only temporary!!!
      if (inp.getNodeValue("input@type", siminptype) && siminptype == "oldstyle")
        initConfigOldstyle(siminpfile);
    }
  }
  else {
    global_log->error() << "XML input " << inputfilename << ": no simulation section" << endl;
  }
}

void Simulation::initConfigOldstyle(const string& inputfilename) {
  int ownrank = 0;
#ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &ownrank);
#endif

  global_log->info() << "init oldstyle config file: " << inputfilename << endl;

  // open filestream to the input file
  ifstream inputfilestream(inputfilename.c_str());
  if (!inputfilestream.is_open()) {
    global_log->error() << "Could not open file " << inputfilename << endl;
    exit(1);
  }

  //  std::string inputPath;
  //  unsigned int lastIndex = inputfilename.find_last_of('/',inputfilename.size()-1);
  //  if (lastIndex == string::npos)
  //    inputPath="";
  //  else
  //    inputPath = inputfilename.substr(0, lastIndex+1);


  // used to store one token of the inputfilestream
  string token;

  double timestepLength;
  unsigned cosetid = 0;

  global_log->info() << "Constructing domain ..." << endl;
  _domain = new Domain(ownrank);
  global_log->info() << "Domain construction done." << endl;
  _particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

  // The first line of the config file has to contain the token "MDProjectConfig"
  inputfilestream >> token;
  if ((token != "mardynconfig") && (token != "MDProjectConfig")) {
    global_log->error() << "Not a mardynconfig file! First token: " << token << endl;
    exit(1);
  }

  while (inputfilestream) {
    token.clear();
    inputfilestream >> token;
    global_log->debug() << " [[" << token << "]]" << endl;

    if (token.substr(0, 1) == "#") {
      inputfilestream.ignore(1024, '\n');
      continue;
    }
    if (token == "phaseSpaceFile") {
      string phaseSpaceFileFormat;
      inputfilestream >> phaseSpaceFileFormat;

      if (timestepLength == 0.0) {
        global_log->error() << "timestep missing." << endl;
        exit(1);
      }
      if (phaseSpaceFileFormat == "OldStyle") {
        string phaseSpaceFileName;
        inputfilestream >> phaseSpaceFileName;
        _inputReader = (InputBase*) new InputOldstyle();
        _inputReader->setPhaseSpaceFile(phaseSpaceFileName);
        _inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
        _inputReader->readPhaseSpaceHeader(_domain, timestepLength);
      }
      else if (phaseSpaceFileFormat == "PartGen") {
        string mode;
        string phaseSpaceFileName;
        inputfilestream >> mode >> phaseSpaceFileName;
        _inputReader = (InputBase*) new PartGen();
        _inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
        // PartGen has to modes, a "Homogeneous" mode, where particles
        // with a homogeneous distribution are created, and a "Cluster" mode,
        // where droplets are created. Currently, only the cluster mode is supported,
        // which needs another config line starting with "clusterFile ..." directly
        // after the config line starting with "phaseSpaceFile"
        double gasDensity;
        double fluidDensity;
        double volPercOfFluid;
        string clusterFileName;
        inputfilestream >> token >> gasDensity >> fluidDensity >> volPercOfFluid >> clusterFileName;
        ((PartGen*) _inputReader)->setClusterFile(gasDensity, fluidDensity, volPercOfFluid, clusterFileName);
        ((PartGen*) _inputReader)->readPhaseSpaceHeader(_domain, timestepLength);
      }
      if (this->_LJCutoffRadius == 0.0)
        _LJCutoffRadius = this->_cutoffRadius;
      _domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
    }
    else if (token == "timestepLength") {
      inputfilestream >> timestepLength;
    }
    else if (token == "cutoffRadius") {
      inputfilestream >> _cutoffRadius;
    }
    else if (token == "LJCutoffRadius") {
      inputfilestream >> _LJCutoffRadius;
    }
    else if ((token == "parallelization") || (token == "parallelisation")) {
#ifndef PARALLEL
      global_log->warning() << "Input file demands parallelization, but the current compilation doesn't\n\tsupport parallel execution.\n" << endl;
      inputfilestream >> token;
#else
      inputfilestream >> token;
      if (token=="DomainDecomposition") {
        _domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
      }
      else if(token=="KDDecomposition") {
        _domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, 1.0, 0.0);
      }
#endif
    }
    else if (token == "datastructure") {
      inputfilestream >> token;
      if (token == "LinkedCells") {
        int cellsInCutoffRadius;
        inputfilestream >> cellsInCutoffRadius;
        double bBoxMin[3];
        double bBoxMax[3];
        for (int i = 0; i < 3; i++) {
          bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
          bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
        }
        if (this->_LJCutoffRadius == 0.0)
          _LJCutoffRadius = this->_cutoffRadius;
        _moleculeContainer = new LinkedCells(bBoxMin, bBoxMax,
                                             _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius,
                                             cellsInCutoffRadius, _particlePairsHandler);
      }
      else if (token == "AdaptiveSubCells") {
        double bBoxMin[3];
        double bBoxMax[3];
        for (int i = 0; i < 3; i++) {
          bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
          bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
        }
        //creates a new Adaptive SubCells datastructure
        if (_LJCutoffRadius == 0.0)
          _LJCutoffRadius = _cutoffRadius;
        _moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax,
                                                  _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius,
                                                  _particlePairsHandler);
      }
    }
    else if (token == "output") {
      inputfilestream >> token;
      if (token == "ResultWriter") {
        string outputPathAndPrefix;
        inputfilestream >> outputPathAndPrefix;
        _outputPlugins.push_back(new ResultWriter(outputPathAndPrefix));
        global_log->debug() << "ResultWriter '" << outputPathAndPrefix << "'.\n";
      }
      else if (token == "XyzWriter") {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new XyzWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
        global_log->debug() << "XyzWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
        _outputFrequency = writeFrequency;
      }
      else if (token == "PovWriter") {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new PovWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
        global_log->debug() << "POVWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
      }
      else if (token == "DecompWriter") {
        unsigned long writeFrequency;
        string mode;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> mode >> outputPathAndPrefix;
        _outputPlugins.push_back(new DecompWriter(writeFrequency, mode, outputPathAndPrefix, _numberOfTimesteps, true));
        global_log->debug() << "DecompWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
      }
      else if ((token == "VisittWriter") || (token == "VISWriter")) {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new VISWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
        global_log->debug() << "VISWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
      }
      /*
       else if(token == "VimWriter")
       {
       unsigned long writeFrequency;

       string outputPathAndPrefix;
       inputfilestream >> writeFrequency >> outputPathAndPrefix;
       _outputPlugins.push_back(new VimWriter(_numberOfTimesteps, writeFrequency, outputPathAndPrefix, true));
       if(!ownrank) cout << "Vim " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
       }
       */
    }
    else if (token == "accelerate") {
      cosetid++;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if (token != "comp") {
        global_log->error() << "Expected 'comp' instead of '" << token << "'.\n";
        exit(1);
      }
      int cid = 0;
      while (cid >= 0) {
        inputfilestream >> cid;
        if (cid > 0)
          global_log->info() << "acc. for component " << cid << endl;
        cid--;
        _domain->assignCoset((unsigned) cid, cosetid);
      }
      double v;
      inputfilestream >> v;
      global_log->debug() << "velocity " << v << endl;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if (token != "towards") {
        global_log->error() << "Expected 'towards' instead of '" << token << "'.\n";
        exit(1);
      }
      double dir[3];
      double dirnorm = 0;
      for (unsigned d = 0; d < 3; d++) {
        inputfilestream >> dir[d];
        dirnorm += dir[d] * dir[d];
      }
      dirnorm = 1.0 / sqrt(dirnorm);
      for (unsigned d = 0; d < 3; d++)
        dir[d] *= dirnorm;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if (token != "within") {
        global_log->error() << "Expected 'within' instead of '" << token << "'.\n";
        exit(1);
      }
      double tau;
      inputfilestream >> tau;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if (token != "from") {
        global_log->error() << "Expected 'from' instead of '" << token << "'.\n";
        exit(1);
      }
      double ainit[3];
      for (unsigned d = 0; d < 3; d++)
        inputfilestream >> ainit[3];
      if (timestepLength == 0.0) {
        global_log->error() << "timestep missing." << endl;
        exit(1);
      }
      _domain->specifyComponentSet(cosetid, dir, tau, ainit, timestepLength);
    }
    else if (token == "constantAccelerationTimesteps") {
      unsigned uCAT;
      inputfilestream >> uCAT;
      _domain->setUCAT(uCAT);
    }
    else if (token == "profile") {
      unsigned xun, yun, zun;
      inputfilestream >> xun >> yun >> zun;
      _domain->setupProfile(xun, yun, zun);
      _doRecordProfile = true;
    }
    else if (token == "profileRecordingTimesteps") {
      inputfilestream >> _profileRecordingTimesteps;
    }
    else if (token == "profileOutputTimesteps") {
      inputfilestream >> _profileOutputTimesteps;
    }
    else if (token == "RDF") {
      double interval;
      unsigned bins;
      inputfilestream >> interval >> bins;
      _domain->setupRDF(interval, bins);
      _doRecordRDF = true;
    }
    else if (token == "RDFOutputTimesteps") {
      inputfilestream >> _RDFOutputTimesteps;
    }
    else if (token == "RDFOutputPrefix") {
      inputfilestream >> _RDFOutputPrefix;
    }
    else if (token == "resultOutputTimesteps") {
      inputfilestream >> _resultOutputTimesteps;
    }
    else if (token == "profiledComponent") {
      unsigned cid;
      inputfilestream >> cid;
      cid--;
      _domain->considerComponentInProfile(cid);
    }
    else if (token == "profileOutputPrefix") {
      inputfilestream >> _profileOutputPrefix;
    }
    else if (token == "collectThermostatDirectedVelocity") {
      inputfilestream >> _collectThermostatDirectedVelocity;
    }
    else if (token == "zOscillator") {
      _zoscillation = true;
      inputfilestream >> _zoscillator;
    }
    // chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
    // to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
    else if (token == "chemicalPotential") {
      double imu;
      inputfilestream >> imu;
      inputfilestream >> token;
      if (token != "component") {
        global_log->error() << "Expected 'component' instead of '" << token << "'.\n";
        global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps" << endl;
        exit(1);
      }
      unsigned icid;
      inputfilestream >> icid;
      icid--;
      inputfilestream >> token;
      double x0, y0, z0, x1, y1, z1;
      bool controlVolume = false;
      if (token == "control") {
        controlVolume = true;
        inputfilestream >> x0 >> y0 >> z0;
        inputfilestream >> token;
        if (token != "to") {
          global_log->error() << "Expected 'to' instead of '" << token << "'.\n";
          global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps\n";
          exit(1);
        }
        inputfilestream >> x1 >> y1 >> z1;
        inputfilestream >> token;
      }
      if (token != "conduct") {
        global_log->error() << "Expected 'conduct' instead of '" << token << "'.\n";
        global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps\n";
        exit(1);
      }
      unsigned intest;
      inputfilestream >> intest;
      inputfilestream >> token;
      if (token != "tests") {
        global_log->error() << "Expected 'tests' instead of '" << token << "'.\n";
        global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps" << endl;
        exit(1);
      }
      inputfilestream >> token;
      if (token != "every") {
        global_log->error() << "Expected 'every' instead of '" << token << "'.\n";
        global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps" << endl;
        exit(1);
      }
      unsigned instep;
      inputfilestream >> instep;
      inputfilestream >> token;
      if (token != "steps") {
        global_log->error() << "Expected 'steps' instead of '" << token << "'.\n";
        global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> " << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps" << endl;
        exit(1);
      }
      ChemicalPotential tmu = ChemicalPotential();
      tmu.setMu(icid, imu);
      tmu.setInterval(instep);
      tmu.setInstances(intest);
      if (controlVolume)
        tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
      global_log->info() << setprecision(6)
                         << "chemical Potential " << imu
                         << " component " << icid + 1 << " (internally " << icid << ") conduct "
                         << intest << " tests every " << instep << " steps: ";
      _lmu.push_back(tmu);
      global_log->info() << " pushed back." << endl;
    }
    else if (token == "planckConstant") {
      inputfilestream >> h;
    }
    else if (token == "NVE") {
      _domain->thermostatOff();
    }
    else if (token == "initCanonical") {
      inputfilestream >> _initCanonical;
    }
    else if (token == "initGrandCanonical") {
      inputfilestream >> _initGrandCanonical;
    }
    else if (token == "initStatistics") {
      inputfilestream >> _initStatistics;
    }
  }

  // read particle data
  unsigned long maxid = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);

  if (this->_LJCutoffRadius == 0.0)
    _LJCutoffRadius = this->_cutoffRadius;
  _domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

  // @todo comment
  _integrator = new Leapfrog(timestepLength);

  // test new Decomposition
  _moleculeContainer->update();
  _moleculeContainer->deleteOuterParticles();

  unsigned idi = _lmu.size();
  unsigned j = 0;
  std::list<ChemicalPotential>::iterator cpit;
  for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
    cpit->setIncrement(idi);
    double tmp_molecularMass = _domain->getComponents()[cpit->getComponentID()].m();
    cpit->setSystem(
        _domain->getGlobalLength(0), _domain->getGlobalLength(1),
        _domain->getGlobalLength(2), tmp_molecularMass
    );
    cpit->setGlobalN(
        _domain->N(cpit->getComponentID())
    );
    cpit->setNextID(j + (int) (1.001 * (256 + maxid)));

    cpit->setSubdomain(
        ownrank,
        _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
        _moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
        _moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2)
    );
    /* TODO: thermostat */
    double Tcur = _domain->getCurrentTemperature(0);
    /* FIXME: target temperature from thermostat ID 0 or 1?  */
    double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1)
                                                : _domain->getTargetTemperature(0);
    if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
      Tcur = Ttar;
    cpit->submitTemperature(Tcur);
    if (h != 0.0)
      cpit->setPlanckConstant(h);

    j++;
  }
#ifdef STEEREO
  SteereoLogger::setOutputLevel(2);
  _steer = new SimSteering();
  _steer->setNumberOfQueues(1);
  int portNumber = 44445;
#ifdef PARALLEL
  int ownSize;
  MPI_Comm_size(MPI_COMM_WORLD, &ownSize);
  SteereoMPIIntraCommunicator* mpiIntraComm = new SteereoMPIIntraCommunicator();
  int partNum = 1;
#ifdef STEEREO_PARTITIONS
  partNum = STEEREO_PARTITIONS;
#endif /* STEEREO_PARTITIONS */
  char* partVal = getenv("STEEREO_PARTITIONS");
  if (partVal != NULL) {
    partNum = atoi(partVal);
  }
  global_log->debug() << "going to divide the " << ownSize << " processes into " << partNum << " partitions" << std::endl;
  mpiIntraComm->generateEqually(ownrank, partNum, ownSize);
  global_log->debug() << "equally generated" << std::endl;
  _steer->setIntraCommunicator(mpiIntraComm);
  //if (ownrank == mpiIntraComm->getRoot())
  if (mpiIntraComm->amIRoot()) {
    portNumber += (ownrank * partNum) / ownSize;
#endif /* PARALLEL */
    std::stringstream strstr;
    strstr << portNumber;
    _steer->setCommunicator(new SteereoSocketCommunicator(strstr.str()));
#ifdef PARALLEL
  }
#endif /* PARALLEL */
#endif /* STEEREO */

}

void Simulation::initialize() {
  global_log->info() << "Initializing simulation" << endl;
  // clear halo
  global_log->info() << "Clearing halos" << endl;
  _moleculeContainer->deleteOuterParticles();

  global_log->info() << "Updating domain decomposition" << endl;
  updateParticleContainerAndDecomposition();

  // Force calculation
  global_log->info() << "Performing force calculation" << endl;
  _moleculeContainer->traversePairs();

  // clear halo
  global_log->info() << "Clearing halos" << endl;
  _moleculeContainer->deleteOuterParticles();

  // initialize the radial distribution function
  if (_doRecordRDF)
    _domain->resetRDF();

  if (_domain->isAcceleratingUniformly()) {
    global_log->info() << "Initialising uniform acceleration." << endl;
    unsigned long uCAT = _domain->getUCAT();
    global_log->info() << "uCAT: " << uCAT << " steps." << endl;
    _domain->determineAdditionalAcceleration(
        _domainDecomposition,
        _moleculeContainer,
        uCAT * _integrator->getTimestepLength()
    );
    global_log->info() << "Uniform acceleration initialised." << endl;
  }

  global_log->info() << "Calculating global values" << endl;
  _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
  _domain->calculateVelocitySums(_moleculeContainer);
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, true, 1.0);

  if (_lmu.size() > 0) {
    /* TODO: thermostat */
    double Tcur = _domain->getGlobalCurrentTemperature();
    /* FIXME: target temperature from thermostat ID 0 or 1? */
    double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1)
                                                : _domain->getTargetTemperature(0);
    if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
      Tcur = Ttar;

    list<ChemicalPotential>::iterator cpit;
    if (h == 0.0)
      h = sqrt(6.2831853 * Ttar);
    for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
      cpit->submitTemperature(Tcur);
      cpit->setPlanckConstant(h);
    }
  }

  if (_zoscillation) {
    global_log->debug() << "Initializing z-oscillators" << endl;
    _integrator->init1D(_zoscillator, _moleculeContainer);
  }

  // initialize output
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain);
  }

#ifdef STEEREO
  MegaMolSnapshotCommand::setSimData(this);
  _steer->registerCommand(MegaMolSnapshotCommand::generateNewInstance, "getMegaMolSnapshot");
  _steer->registerCommand(SendCouplingMDCommand::generateNewInstance, "sendCouplingMD");
  _steer->registerCommand(ReceiveCouplingMDCommand::generateNewInstance, "receiveCouplingMD");
  _steer->registerCommand(SnapshotCommand::generateNewInstance, "getSnapshot");
  SnapshotCommand::setSimData(this);
  SteerParameterCommand::registerScalarParameter("temp", _domain, &Domain::getGlobalCurrentTemperature, &Domain::setGlobalTemperature);
  SendCouplingMDCommand::addData(this);
  ReceiveCouplingMDCommand::addData(this);
#ifdef PARALLEL
  int ownrank = _domainDecomposition->getRank();
  global_log->debug() << "_ownrank is " << ownrank << std::endl;
  global_log->debug() << "my local rank is " << _steer->getIntraCommunicator()->getRank() << std::endl;
  global_log->debug() << "my local root rank is " << _steer->getIntraCommunicator()->getRoot() << std::endl;
  //if (ownrank == _steer->getIntraCommunicator()->getRoot())
  if (_steer->getIntraCommunicator()->amIRoot()) {
    global_log->debug() << "going to start listening on rank " << ownrank << std::endl;
#endif
    _steer->startListening();
#ifdef PARALLEL
  }
#endif
#endif

  // activate RDF sampling
  if ((this->_initSimulation > this->_initStatistics) && this->_doRecordRDF)
    this->_domain->tickRDF();

  global_log->info() << "System initialised\n" << endl;
}

void Simulation::simulate() {

  Molecule* tM;

  global_log->info() << "Started simulation" << endl;

  // (universal) constant acceleration (number of) timesteps
  unsigned uCAT = _domain->getUCAT();

  /***************************************************************************/
  /* BEGIN MAIN LOOP                                                         */
  /***************************************************************************/
  // all timers except the ioTimer messure inside the main loop
  Timer loopTimer; // Timer for computation
  Timer perStepIoTimer; // Timer for IO during simulation steps
  Timer ioTimer; // Timer for final IO

  _initSimulation = (unsigned long) (_domain->getCurrentTime() / _integrator->getTimestepLength());

  loopTimer.start();
  for (_simstep = _initSimulation; _simstep <= _numberOfTimesteps; _simstep++) {
    if (_simstep >= _initGrandCanonical) {
      unsigned j = 0;
      list<ChemicalPotential>::iterator cpit;
      for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
        if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
          cpit->prepareTimestep(_moleculeContainer, _domainDecomposition);
        }
        j++;
      }
    }
    global_log->debug() << "timestep " << _simstep << endl;

    _integrator->eventNewTimestep(_moleculeContainer, _domain);

    // activate RDF sampling
    if ((_simstep == this->_initStatistics) && this->_doRecordRDF)
      this->_domain->tickRDF();

    // ensure that all Particles are in the right cells and exchange Particles
    global_log->debug() << "Updating container and decomposition" << endl;
    updateParticleContainerAndDecomposition();

    // Force calculation
    global_log->debug() << "Traversing pairs" << endl;
    _moleculeContainer->traversePairs();

    // test deletions and insertions
    if (_simstep >= _initGrandCanonical) {
      unsigned j = 0;
      list<ChemicalPotential>::iterator cpit;
      for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
        if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
          global_log->debug() << "Grand canonical ensemble(" << j
                              << "): test deletions and insertions" << endl;
          /* TODO: thermostat */
          _moleculeContainer->grandcanonicalStep(&(*cpit), _domain->getGlobalCurrentTemperature());
#ifndef NDEBUG
          /* silly check if random numbers inserted are the same for all processes... */
          cpit->assertSynchronization(_domainDecomposition);
#endif

          int localBalance = _moleculeContainer->localGrandcanonicalBalance();
          int balance = _moleculeContainer->grandcanonicalBalance(_domainDecomposition);
          global_log->debug() << "   b[" << ((balance > 0) ? "+" : "") << balance
                              << "(" << ((localBalance > 0) ? "+" : "")
                              << localBalance << ")" << " / c = "
                              << cpit->getComponentID() << "]   " << endl;
          _domain->Nadd(cpit->getComponentID(), balance, localBalance);
        }

        j++;
      }
    }

    // clear halo
    global_log->debug() << "Delete outer particles" << endl;
    _moleculeContainer->deleteOuterParticles();

    if (_simstep >= _initGrandCanonical) {
      _domain->evaluateRho(_moleculeContainer->getNumberOfParticles(), _domainDecomposition);
    }

    if (!(_simstep % _collectThermostatDirectedVelocity))
      _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
    if (_domain->isAcceleratingUniformly()) {
      if (!(_simstep % uCAT)) {
        global_log->debug() << "Determine the additional acceleration" << endl;
        _domain->determineAdditionalAcceleration(_domainDecomposition,
                                                 _moleculeContainer,
                                                 uCAT * _integrator->getTimestepLength());
      }
      global_log->debug() << "Process the uniform acceleration" << endl;
      _integrator->accelerateUniformly(_moleculeContainer, _domain);
    }

    /*
     * radial distribution function
     */
    if (_simstep >= _initStatistics) {
      if (this->_lmu.size() == 0) {
        this->_domain->record_cv();
      }
      if (this->_doRecordRDF) {
        this->_domain->tickRDF();
        this->_particlePairsHandler->recordRDF();
        this->_moleculeContainer->countParticles(_domain);
      }
    }

    if (_zoscillation) {
      global_log->debug() << "alert z-oscillators" << endl;
      _integrator->zOscillation(_zoscillator, _moleculeContainer);
    }

    // Inform the integrator about the calculated forces
    global_log->debug() << "Inform the integrator" << endl;
    _integrator->eventForcesCalculated(_moleculeContainer, _domain);

    // calculate the global macroscopic values from the local values
    global_log->debug() << "Calculate macroscopic values" << endl;
    _domain->calculateGlobalValues(_domainDecomposition,
                                   _moleculeContainer,
                                   (!(_simstep % _collectThermostatDirectedVelocity)),
                                   Tfactor(_simstep));

    // scale velocity and angular momentum
    if (!_domain->NVE()) {
      global_log->debug() << "Velocity scaling" << endl;
      if (_domain->severalThermostats()) {
        for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
          int thermostat = _domain->getThermostat(tM->componentid());
          if (0 >= thermostat)
            continue;
          if (_domain->thermostatIsUndirected(thermostat)) {
            /* TODO: thermostat */
            tM->scale_v(_domain->getGlobalBetaTrans(thermostat),
                        _domain->getThermostatDirectedVelocity(thermostat, 0),
                        _domain->getThermostatDirectedVelocity(thermostat, 1),
                        _domain->getThermostatDirectedVelocity(thermostat, 2));
          }
          else {
            tM->scale_v(_domain->getGlobalBetaTrans(thermostat));
          }
          tM->scale_D(_domain->getGlobalBetaRot(thermostat));
        }
      }
      else {
        for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
          tM->scale_v(_domain->getGlobalBetaTrans());
          tM->scale_D(_domain->getGlobalBetaRot());
        }
      }
    }

    _domain->advanceTime(_integrator->getTimestepLength());
#ifdef STEEREO
    _steer -> processQueue (0);
#endif
    // measure per timestep IO
    loopTimer.stop();
    perStepIoTimer.start();
    output(_simstep);
    perStepIoTimer.stop();
    loopTimer.start();
  }
  loopTimer.stop();
  /***************************************************************************/
  /* END MAIN LOOP                                                           */
  /***************************************************************************/

  ioTimer.start();
  string cpfile(_outputPrefix + ".restart.xdr");
  _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
  // finish output
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
    (*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
    delete (*outputIter);
  }
  ioTimer.stop();

  global_log->info() << "Computation in main loop took: " << loopTimer.get_etime() << " sec" << endl;
  global_log->info() << "IO in main loop  took:         " << perStepIoTimer.get_etime() << " sec" << endl;
  global_log->info() << "Final IO took:                 " << ioTimer.get_etime() << " sec" << endl;

  delete _domainDecomposition;
  delete _domain;
  delete _particlePairsHandler;
  delete _moleculeContainer;
  delete _integrator;
  delete _inputReader;

}

void Simulation::output(unsigned long simstep) {
  int ownrank = _domainDecomposition->getRank();

  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
    (*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &(_lmu));
  }

  if (_doRecordRDF && !(simstep % _RDFOutputTimesteps)) {
    _domain->collectRDF(_domainDecomposition);
    if (!ownrank) {
      _domain->accumulateRDF();
      for (unsigned i = 0; i < _numberOfComponents; i++) {
        for (unsigned j = i; j < _numberOfComponents; j++) {
          ostringstream osstrm;
          osstrm << _RDFOutputPrefix << "_" << i << "-" << j << ".";
          osstrm.fill('0');
          osstrm.width(9);
          osstrm << right << simstep;
          _domain->outputRDF(osstrm.str().c_str(), i, j);
          osstrm.str("");
          osstrm.clear();
        }
      }
    }
    _domain->resetRDF();
  }

  if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileRecordingTimesteps)) {
    _domain->recordProfile(_moleculeContainer);
  }
  if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileOutputTimesteps)) {
    _domain->collectProfile(_domainDecomposition);
    if (!ownrank) {
      ostringstream osstrm;
      osstrm << _profileOutputPrefix << ".";
      osstrm.fill('0');
      osstrm.width(9);
      osstrm << right << simstep;
      _domain->outputProfile(osstrm.str().c_str());
      osstrm.str("");
      osstrm.clear();
    }
    _domain->resetProfile();
  }

  if ((_outputFrequency > 0) && (0 == simstep % _outputFrequency)) {
    _moleculeContainer->deleteOuterParticles();
    ostringstream osstrm;
    osstrm.str("");
    osstrm << _outputPrefix;
    osstrm << ".restart.xdr";
    _domain->writeCheckpoint(osstrm.str(), _moleculeContainer, _domainDecomposition);
  }

  if (_domain->thermostatWarning())
    global_log->warning() << "Thermostat!" << endl;
  /* TODO: thermostat */
  global_log->info()
      << "Simstep = " << simstep
      << "\tT = " << _domain->getGlobalCurrentTemperature()
      << "\tU_pot = " << _domain->getAverageGlobalUpot()
      << "\tp = " << _domain->getGlobalPressure()
      << endl;
}

void Simulation::updateParticleContainerAndDecomposition() {

  // The particles have moved, so the neighbourhood relations have
  // changed and have to be adjusted
  _moleculeContainer->update();
  //_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
  _domainDecomposition->balanceAndExchange(true, _moleculeContainer, _domain->getComponents(), _domain);
  // The cache of the molecules must be updated/build after the exchange process,
  // as the cache itself isn't transferred
  Molecule* tM;
  for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
    tM->upd_cache();
  }
}

/* FIXME: we shoud provide a more general way of doing this */
double Simulation::Tfactor(unsigned long simstep) {
  double xi = (double) (simstep - _initSimulation) / (double) (_initCanonical - _initSimulation);
  if ((xi < 0.1) || (xi > 0.9)) return 1.0;
  else if (xi < 0.3) return 15.0 * xi - 0.5;
  else if (xi < 0.4) return 10.0 - 20.0 * xi;
  else if (xi < 0.6) return 2.0;
  else return 4 - 10.0 * xi / 3.0;
}

