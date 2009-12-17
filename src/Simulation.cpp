#include "Simulation.h"
#include "Domain.h"
#include "md_io/PartGen.h"
#include "molecules/Molecule.h"
#include "datastructures/LinkedCells.h"
#include "datastructures/AdaptiveSubCells.h"
#include "parallel/DomainDecompBase.h"
#ifdef PARALLEL
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif
#include "datastructures/adapter/ParticlePairs2PotForceAdapter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "md_io/ResultWriter.h"
#include "md_io/XyzWriter.h"
#include "md_io/PovWriter.h"
#include "md_io/DecompWriter.h"
#include "md_io/CheckpointWriter.h"
#include "md_io/VISWriter.h"
#include "md_io/InputOldstyle.h"
#include "ensemble/GrandCanonical.h"

#ifdef STEEREO
#include "commands/snapshotCommand.h"
#include "commands/megaMolSnapshotCommand.h"
#include "commands/sendCouplingMDCommand.h"
#include "commands/receiveCouplingMDCommand.h"
#include <steerParameterCommand.h>
#include <steereoSocketCommunicator.h>
#endif

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include "Timer.h"
#include "Logger.h"

using Log::global_log;

Simulation::Simulation(int *argc, char ***argv)
{
  int ownrank = 0;
#ifdef PARALLEL
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownrank);
#endif

  /* Initialize the global log file */
  string logfileName("MarDyn");
  global_log = new Log::Logger(Log::ALL);
  global_log->set_mpi_output_root(0);

  if (*argc != 4) 
  {
     global_log->error() << "Usage: " << (*argv)[0] << " <configfilename> <number of timesteps> <outputprefix>\n";
     global_log->error() << "Detected argc: " << *argc << endl;
     exit(1);
  }
  
#ifndef PARALLEL
  global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
  _domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
  global_log->info() << "Initialization done" << endl;
#endif

  // open filestream to the input file
  string inputfilename((*argv)[1]);
  ifstream inputfilestream(inputfilename.c_str());
  if ( !inputfilestream.is_open() ) {
    global_log->error() << "Could not open file " << inputfilename << endl;
    exit (1);
  }

  std::string inputPath;
  unsigned int lastIndex = inputfilename.find_last_of('/',inputfilename.size()-1);
  if (lastIndex == string::npos)
    inputPath="";
  else
    inputPath = inputfilename.substr(0, lastIndex+1);
  
  // store number of timesteps to be simulated
  this->_numberOfTimesteps = atol((*argv)[2]);
  global_log->info() << "Simulating " << _numberOfTimesteps << " steps." << endl;

  // store prefix for output files
  _outputPrefix = string((*argv)[3]);

  // used to store one token of the inputfilestream
  string token;

  double timestepLength;
  unsigned cosetid = 0;
  
  /*
   * default parameters
   */
  this->_resultOutputTimesteps = 25;
  this->_doRecordProfile = false;
  this->_profileRecordingTimesteps = 7;
  this->_profileOutputTimesteps = 12500;
  this->_profileOutputPrefix = "out";
  this->_collectThermostatDirectedVelocity = 100;
  this->_zoscillation = false;
  this->_zoscillator = 512;
  this->_doRecordRDF = false;
  this->_RDFOutputTimesteps = 25000;
  this->_RDFOutputPrefix = "out";
  this->_initCanonical = 5000;
  this->_initGrandCanonical = 10000000;
  this->h = 0.0;
  this->_initStatistics = 20000;
  
  global_log->info() << "Constructing domain ..." << endl;
  this->_domain = new Domain(ownrank);
  global_log->info() << "Domain construction done." << endl;
  _particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

  // The first line of the config file has to contain the token "MDProjectConfig"
  inputfilestream >> token;
  if((token != "mardynconfig") && (token != "MDProjectConfig"))
  {
    global_log->error() << "Not a mardynconfig file! " << token << endl;
    exit(1);
  } 
  
  while (inputfilestream) {
    token.clear();
    inputfilestream >> token;
    global_log->debug() << " [[" << token << "]]" << endl;

    if (token.substr(0, 1)=="#"){
      inputfilestream.ignore(1024,'\n');
      continue;
    } 
    if (token=="phaseSpaceFile") {
      string phaseSpaceFileFormat; 
      inputfilestream >> phaseSpaceFileFormat;

      if(timestepLength == 0.0) {
         global_log->error() << "timestep missing." << endl;
         exit(1);
      }
      if(phaseSpaceFileFormat=="OldStyle"){
        string phaseSpaceFileName;
        inputfilestream >> phaseSpaceFileName;
        _inputReader = (InputBase*) new InputOldstyle();
        _inputReader->setPhaseSpaceFile(phaseSpaceFileName);
        _inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
        _inputReader->readPhaseSpaceHeader(_domain, timestepLength, _cutoffRadius);
      }
      else if(phaseSpaceFileFormat=="PartGen"){
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
        ((PartGen*)_inputReader)->setClusterFile(gasDensity, fluidDensity, volPercOfFluid, clusterFileName);
        ((PartGen*)_inputReader)->readPhaseSpaceHeader(_domain, timestepLength, _cutoffRadius);
      }
      _domain->initParameterStreams(_cutoffRadius);
    } 
    else if (token=="timestepLength") {
      inputfilestream >> timestepLength;
    } 
    else if (token=="cutoffRadius") {
      inputfilestream >> _cutoffRadius;
    } 
    else if ((token=="parallelization") || (token == "parallelisation"))
    {
#ifndef PARALLEL
      global_log->warning() << "WARNING: Input file demands parallelization, but the current compilation doesn't\n\tsupport parallel execution.\n" << endl;
      inputfilestream >> token;
#else
      inputfilestream >> token;
      if (token=="DomainDecomposition") {
        _domainDecomposition = (DomainDecompBase*) new DomainDecomposition(argc, argv);
      } 
      else if(token=="KDDecomposition"){
        _domainDecomposition = (DomainDecompBase*) new KDDecomposition(argc, argv, _cutoffRadius, _domain);
      }
#endif
    }
    else if (token=="datastructure") {
      inputfilestream >> token;
      if (token=="LinkedCells") {
        int cellsInCutoffRadius;
        inputfilestream >> cellsInCutoffRadius;
        double bBoxMin[3];
        double bBoxMax[3];
        for (int i=0; i<3; i++) {
          bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
          bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
        }
        _moleculeContainer = new LinkedCells(bBoxMin, bBoxMax,
            _cutoffRadius, _tersoffCutoffRadius, cellsInCutoffRadius, *_particlePairsHandler);
      } 
      else if (token=="AdaptiveSubCells") {
        int cellsInCutoffRadius;
        inputfilestream >> cellsInCutoffRadius;
        double bBoxMin[3];
        double bBoxMax[3];
        for (int i=0; i<3; i++) {
          bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
          bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
        }
        //creates a new Adaptive SubCells datastructure
        _moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax,
              _cutoffRadius, cellsInCutoffRadius, *_particlePairsHandler);
      } 
    } 
    else if (token=="output") {
      inputfilestream >> token;
      if (token=="ResultWriter") {
        string outputPathAndPrefix;
        inputfilestream >> outputPathAndPrefix;
        _outputPlugins.push_back(new ResultWriter(outputPathAndPrefix));
      } 
      else if (token=="XyzWriter") {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new XyzWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
        this->_outputFrequency = writeFrequency;
      } 
      else if (token=="PovWriter") {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new PovWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
      }
      else if (token=="DecompWriter") {
        unsigned long writeFrequency;
        string mode;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> mode >> outputPathAndPrefix;
        _outputPlugins.push_back(new DecompWriter(writeFrequency, mode, outputPathAndPrefix, _numberOfTimesteps, true));
      }
      else if((token == "VisittWriter") || (token == "VISWriter"))
      {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new VISWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
        global_log->debug() << "VisItt " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
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
    else if(token == "accelerate")
    {
      cosetid++;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if(token != "comp")
      {
         global_log->error() << "Expected 'comp' instead of '" << token << "'.\n";
         exit(1);
      }
      int cid = 0;
      while(cid >= 0)
      {
         inputfilestream >> cid;
         if(cid > 0) global_log->info() << "acc. for component " << cid << endl;
         cid--;
         _domain->assignCoset((unsigned)cid, cosetid);
      }
      double v;
      inputfilestream >> v;
      global_log->debug() << "velocity " << v << endl;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if(token != "towards") {
         global_log->error() << "Expected 'towards' instead of '" << token << "'.\n";
         exit(1);
      }
      double dir[3];
      double dirnorm = 0;
      for(unsigned d = 0; d < 3; d++) {
         inputfilestream >> dir[d];
         dirnorm += dir[d]*dir[d];
      }
      dirnorm = 1.0 / sqrt(dirnorm);
      for(unsigned d = 0; d < 3; d++) dir[d] *= dirnorm;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if(token != "within") {
         global_log->error() << "Expected 'within' instead of '" << token << "'.\n";
         exit(1);
      }
      double tau;
      inputfilestream >> tau;
      inputfilestream >> token;
      global_log->debug() << "Found specifier '" << token << "'\n";

      if(token != "from") {
         global_log->error() << "Expected 'from' instead of '" << token << "'.\n";
         exit(1);
      }
      double ainit[3];
      for(unsigned d = 0; d < 3; d++) inputfilestream >> ainit[3];
      if(timestepLength == 0.0) {
         global_log->error() << "timestep missing." << endl;
         exit(1);
      }
      this->_domain->specifyComponentSet(cosetid, dir, tau, ainit, timestepLength);
    }
    else if(token == "constantAccelerationTimesteps")
    {
       unsigned uCAT;
       inputfilestream >> uCAT;
       this->_domain->setUCAT(uCAT);
    }
    else if(token == "profile")
    {
       unsigned xun, yun, zun;
       inputfilestream >> xun >> yun >> zun;
       this->_domain->setupProfile(xun, yun, zun);
       this->_doRecordProfile = true;
    }
    else if(token == "profileRecordingTimesteps")
    {
       inputfilestream >> this->_profileRecordingTimesteps;
    }
    else if(token == "profileOutputTimesteps")
    {
       inputfilestream >> this->_profileOutputTimesteps;
    }
    else if(token == "RDF")
    {
       double interval;
       unsigned bins;
       inputfilestream >> interval >> bins;
       this->_domain->setupRDF(interval, bins);
       this->_doRecordRDF = true;
    }
    else if(token == "RDFOutputTimesteps")
    {
       inputfilestream >> this->_RDFOutputTimesteps;
    }
    else if(token == "RDFOutputPrefix")
    {
       inputfilestream >> this->_RDFOutputPrefix;
    }
    else if(token == "resultOutputTimesteps")
    {
       inputfilestream >> this->_resultOutputTimesteps;
    }
    else if(token == "profiledComponent")
    {
       unsigned cid;
       inputfilestream >> cid;
       cid--;
       this->_domain->considerComponentInProfile(cid);
    }
    else if(token == "profileOutputPrefix")
    {
       inputfilestream >> this->_profileOutputPrefix;
    }
    else if(token == "collectThermostatDirectedVelocity")
    {
       inputfilestream >> this->_collectThermostatDirectedVelocity;
    }
    else if(token == "zOscillator")
    {
       this->_zoscillation = true;
       inputfilestream >> this->_zoscillator;
    }
    // chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
    // to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
    else if(token == "chemicalPotential")
    {
       double imu;
       inputfilestream >> imu;
       inputfilestream >> token;
       if(token != "component")
       {
	  global_log->error() << "Expected 'component' instead of '" << token << "'.\n";
	  global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
	              << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
	              << "conduct <ntest> tests every <nstep> steps" << endl;
          exit(1);
       }
       unsigned icid;
       inputfilestream >> icid;
       icid--;
       inputfilestream >> token;
       double x0, y0, z0, x1, y1, z1;
       bool controlVolume = false;
       if(token == "control")
       {
	  controlVolume = true;
          inputfilestream >> x0 >> y0 >> z0;
          inputfilestream >> token;
          if(token != "to")
          {
             global_log->error() << "Expected 'to' instead of '" << token << "'.\n";
             global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
                         << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
	                 << "conduct <ntest> tests every <nstep> steps\n";
             exit(1);
          }
          inputfilestream >> x1 >> y1 >> z1;
          inputfilestream >> token;
       }
       if(token != "conduct")
       {
          global_log->error() << "Expected 'conduct' instead of '" << token << "'.\n";
          global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
	              << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		      << "conduct <ntest> tests every <nstep> steps\n";
          exit(1);
       }
       unsigned intest;
       inputfilestream >> intest;
       inputfilestream >> token;
       if(token != "tests")
       {
          global_log->error() << "Expected 'tests' instead of '" << token << "'.\n";
          global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
	              << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
	              << "conduct <ntest> tests every <nstep> steps" << endl;
          exit(1);
       }
       inputfilestream >> token;
       if(token != "every")
       {
          global_log->error() << "Expected 'every' instead of '" << token << "'.\n";
          global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
	              << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		      << "conduct <ntest> tests every <nstep> steps" << endl;
          exit(1);
       }
       unsigned instep;
       inputfilestream >> instep;
       inputfilestream >> token;
       if(token != "steps")
       {
          global_log->error() << "Expected 'steps' instead of '" << token << "'.\n";
          global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
	              << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		      << "conduct <ntest> tests every <nstep> steps" << endl;
          exit(1);
       }
       ChemicalPotential tmu = ChemicalPotential();
       tmu.setMu(icid, imu);
       tmu.setInterval(instep); 
       tmu.setInstances(intest);
       if(controlVolume) tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
       global_log->info().precision(6);
       global_log->info() << "chemical Potential " << imu << " component "
	          << icid+1 << " (internally " << icid << ") conduct "
	          << intest << " tests every " << instep << " steps: ";
       this->_lmu.push_back(tmu);
       global_log->info() << " pushed back." << endl;
    }
    else if(token == "planckConstant")
    {
       inputfilestream >> this->h;
    }
    else if(token == "NVE")
    {
       this->_domain->thermostatOff();
    }
    else if(token == "initCanonical")
    {
       inputfilestream >> this->_initCanonical;
    }
    else if(token == "initGrandCanonical")
    {
       inputfilestream >> this->_initGrandCanonical;
    }
    else if(token == "initStatistics")
    {
       inputfilestream >> this->_initStatistics;
    }
  }

  // read particle data
  unsigned long maxid = _inputReader->readPhaseSpace(
     this->_moleculeContainer, this->_cutoffRadius, &this->_lmu,
     this->_domain, this->_domainDecomposition
  );

  _domain->initFarFieldCorr(_cutoffRadius);

  // @todo comment
  _integrator = new Leapfrog(timestepLength);

  // test new Decomposition
  _moleculeContainer->update();
  _moleculeContainer->deleteOuterParticles();

  unsigned idi = this->_lmu.size();
  unsigned j = 0;
  std::list<ChemicalPotential>::iterator cpit;
  for(cpit = this->_lmu.begin(); cpit != this->_lmu.end(); cpit++)
  {
     cpit->setIncrement(idi);
     double tmp_molecularMass = this->_domain->getComponents()[
        cpit->getComponentID()
     ].m();
     cpit->setSystem(
        _domain->getGlobalLength(0), _domain->getGlobalLength(1),
        _domain->getGlobalLength(2), tmp_molecularMass
     );
     cpit->setGlobalN(
        this->_domain->N(cpit->getComponentID())
     );
     cpit->setNextID(j + (int)(1.001 * (256 + maxid)));

     cpit->setSubdomain(
        ownrank, _moleculeContainer->getBoundingBoxMin(0),
        _moleculeContainer->getBoundingBoxMax(0),
        _moleculeContainer->getBoundingBoxMin(1),
        _moleculeContainer->getBoundingBoxMax(1),
        _moleculeContainer->getBoundingBoxMin(2),
        _moleculeContainer->getBoundingBoxMax(2)
     );
     double Tcur = this->_domain->T(0);
     double Ttar = _domain->severalThermostats()? _domain->targetT(1)
                                                : _domain->targetT(0);
     if((Tcur < 0.85*Ttar) || (Tcur > 1.15*Ttar)) Tcur = Ttar;
     cpit->submitTemperature(Tcur);
     if(this->h != 0.0) cpit->setPlanckConstant(this->h);
     
     j++;
  }
#ifdef STEEREO
	_steer = new SimSteering ();
  _steer->setNumberOfQueues (1);

#ifdef PARALLEL
  if (ownrank == 0)
  {
#endif
    _steer->setCommunicator (new SteereoSocketCommunicator ("44445"));
#ifdef PARALLEL
  }
#endif
#endif
}

void Simulation::initialize()
{
  global_log->set_mpi_output_root(0);
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
  if(this->_doRecordRDF) this->_domain->resetRDF();

  if(this->_domain->isAcceleratingUniformly())
  {
     global_log->info() << "Initialising uniform acceleration." << endl;
     unsigned long uCAT = this->_domain->getUCAT();
     global_log->info() << "uCAT: " << uCAT << " steps." << endl;
     this->_domain->determineAdditionalAcceleration( this->_domainDecomposition,
                                                     this->_moleculeContainer,
                                                     uCAT * _integrator->getTimestepLength() );
     global_log->info() << "Uniform acceleration initialised." << endl;
  }
  
  global_log->info() << "Calculating global values" << endl;
  _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
  _domain->calculateVelocitySums(_moleculeContainer);
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, true, 1.0);
  
  if(this->_lmu.size() > 0)
  {
     double Tcur = this->_domain->T(0);
     double Ttar = this->_domain->severalThermostats()? this->_domain->targetT(1)
                                                      : this->_domain->targetT(0);
     if((Tcur < 0.85*Ttar) || (Tcur > 1.15*Ttar)) Tcur = Ttar;
  
     list<ChemicalPotential>::iterator cpit;
     if(this->h == 0.0) this->h = sqrt(6.2831853 * Ttar);
     for(cpit = this->_lmu.begin(); cpit != this->_lmu.end(); cpit++)
     {
        cpit->submitTemperature(Tcur);
        cpit->setPlanckConstant(this->h);
     }
  }

  if(this->_zoscillation)
  {
    global_log->debug() << "Initializing z-oscillators" << endl;
    this->_integrator->init1D(this->_zoscillator, this->_moleculeContainer);
  }
  
  // initialize output 
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain);
  }

#ifdef STEEREO
  MegaMolSnapshotCommand::setSimData (this);
  _steer->registerCommand (MegaMolSnapshotCommand::generateNewInstance, "getMegaMolSnapshot");
  _steer->registerCommand (SendCouplingMDCommand::generateNewInstance, "sendCouplingMD");
  _steer->registerCommand (ReceiveCouplingMDCommand::generateNewInstance, "receiveCouplingMD");
  _steer->registerCommand (SnapshotCommand::generateNewInstance, "getSnapshot");
  SnapshotCommand::setSimData (this);
  SteerParameterCommand::registerScalarParameter ("temp", _domain, &Domain::getGlobalTemperature,
																									&Domain::setGlobalTemperature);
  SendCouplingMDCommand::setData(this);
  ReceiveCouplingMDCommand::setData(this);
#ifdef PARALLEL
  std::cout << "_ownrank is " << ownrank << std::endl;
  if (this->_domain->_ownrank == 0)
  {
#endif
    _steer->startListening ();
#ifdef PARALLEL
  }
#endif
#endif
  global_log->info() << "System initialised\n" << endl;
}

void Simulation::simulate()
{

  Molecule* tM;

  global_log->set_mpi_output_root(0);
  global_log->info() << "Started simulation" << endl;
  

  // (universal) constant acceleration (number of) timesteps
  unsigned uCAT = this->_domain->getUCAT();

  /***************************************************************************/
  /* BEGIN MAIN LOOP                                                         */
  /***************************************************************************/
  Timer loopTimer;
  loopTimer.start();
  
  this->_initSimulation = (unsigned long)(this->_domain->getCurrentTime()
                             / _integrator->getTimestepLength());

  for (unsigned long simstep=this->_initSimulation; simstep<=_numberOfTimesteps; simstep++)
  {
    if(simstep >= this->_initGrandCanonical)
    {
       unsigned j = 0;
       list<ChemicalPotential>::iterator cpit;
       for( cpit = this->_lmu.begin();
	    cpit != this->_lmu.end();
	    cpit++ )
       {
	  if(!((simstep + 2*j + 3) % cpit->getInterval()))
	  {
	     cpit->prepareTimestep(
	        this->_moleculeContainer, this->_domainDecomposition
             );
	  }
	  j++;
       }
    }
    global_log->debug() << "timestep " << simstep << endl;
     
    _integrator->eventNewTimestep(_moleculeContainer, _domain);

    // ensure that all Particles are in the right cells and exchange Particles
    global_log->debug() << "Updating container and decomposition" << endl;
    updateParticleContainerAndDecomposition();

    // Force calculation
    global_log->debug() << "Traversing pairs" << endl;
    _moleculeContainer->traversePairs();
    
    /*
     * test deletions and insertions
     */
    if(simstep >= _initGrandCanonical)
    {
       unsigned j = 0;
       list<ChemicalPotential>::iterator cpit;
       for( cpit = this->_lmu.begin();
	    cpit != this->_lmu.end();
	    cpit++ )
       {
	  if(!((simstep + 2*j + 3) % cpit->getInterval()))
	  {     
             global_log->debug() << "Grand canonical ensemble(" << j
		         << "): test deletions and insertions" << endl;
             this->_moleculeContainer->grandcanonicalStep(
	        &(*cpit), this->_domain->T(0)
             );
	     cpit->assertSynchronization(this->_domainDecomposition);
	     
	     int localBalance = _moleculeContainer->localGrandcanonicalBalance();
	     int balance = _moleculeContainer->grandcanonicalBalance(
	        this->_domainDecomposition
	     );
             global_log->debug() << "   b[" << ((balance > 0)? "+": "") << balance
                         << "(" << ((localBalance > 0)? "+": "")
		         << localBalance << ")" << " / c = "
		         << cpit->getComponentID() << "]   " << endl;
             this->_domain->Nadd(
	        cpit->getComponentID(), balance, localBalance
	     );
	  }
	  
	  j++;
       }
    }

    // clear halo
    global_log->debug() << "Delete outer particles" << endl;
    _moleculeContainer->deleteOuterParticles();
    
    /*
     * radial distribution function
     */
    if((simstep >= this->_initStatistics) && (this->_doRecordRDF))
    {
       this->_particlePairsHandler->recordRDF();
       this->_moleculeContainer->countParticles(this->_domain);
       this->_domain->tickRDF();
    }

    if(simstep >= _initGrandCanonical)
    {
       this->_domain->evaluateRho
       (
          this->_moleculeContainer->getNumberOfParticles(),
          this->_domainDecomposition
       );
    }
    
    if(!(simstep % _collectThermostatDirectedVelocity))
       _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
    if(this->_domain->isAcceleratingUniformly())
    {
      if(!(simstep % uCAT))
      {
        global_log->debug() << "Determine the additional acceleration" << endl;
        this->_domain->determineAdditionalAcceleration( this->_domainDecomposition,
                                                        this->_moleculeContainer,
                                                        uCAT * _integrator->getTimestepLength() );
      }
      global_log->debug() << "Process the uniform acceleration" << endl;
      this->_integrator->accelerateUniformly(this->_moleculeContainer, this->_domain);
    }

    if(this->_zoscillation)
    {
      global_log->debug() << "alert z-oscillators" << endl;
      this->_integrator->zOscillation(this->_zoscillator, this->_moleculeContainer);
    }

    // Inform the integrator about the calculated forces
    global_log->debug() << "Inform the integrator" << endl;
    _integrator->eventForcesCalculated(_moleculeContainer, _domain);

    // calculate the global macroscopic values from the local values
    global_log->debug() << "Calculate macroscopic values" << endl;
    _domain->calculateGlobalValues( _domainDecomposition, _moleculeContainer,
                                    (!(simstep % _collectThermostatDirectedVelocity)),
                                    Tfactor(simstep) );

    // scale velocity and angular momentum
    if(!this->_domain->NVE())
    {
       global_log->debug() << "Velocity scaling" << endl;
       if(this->_domain->severalThermostats())
       {
          for( tM = _moleculeContainer->begin();
               tM != _moleculeContainer->end();
               tM = _moleculeContainer->next() )
          {
             int thermostat = this->_domain->getThermostat(tM->componentid());
             if(0 >= thermostat) continue;
             if(this->_domain->thermostatIsUndirected(thermostat))
             {
                tM->scale_v( _domain->getGlobalBetaTrans(thermostat),
                             _domain->thermostatv(thermostat, 0),
                             _domain->thermostatv(thermostat, 1),
                             _domain->thermostatv(thermostat, 2)  );
             }
             else
             {
                tM->scale_v(_domain->getGlobalBetaTrans(thermostat));
             }
             tM->scale_D(_domain->getGlobalBetaRot(thermostat));  
          }
       }
       else
       {
          for( tM = _moleculeContainer->begin();
               tM != _moleculeContainer->end();
               tM = _moleculeContainer->next() )
          {
             tM->scale_v(_domain->getGlobalBetaTrans());
             tM->scale_D(_domain->getGlobalBetaRot());
	  }
       }
    }

    _domain->advanceTime(_integrator->getTimestepLength());
#ifdef STEEREO
		_steer -> processQueue (0);
#endif
    output(simstep);
  }
  loopTimer.stop();
  /***************************************************************************/
  /* END MAIN LOOP                                                           */
  /***************************************************************************/

  Timer ioTimer;
  ioTimer.start();
  string cpfile(_outputPrefix+".restart.xdr");
  _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
  // finish output
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
      (*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
    delete (*outputIter);
  }
  ioTimer.stop();

  global_log->info() << "Simulation loop took: " << loopTimer.get_etime() << " sec" << endl;
  global_log->info() << "Simulation io took:   " << ioTimer.get_etime() << " sec" << endl;

  delete _domainDecomposition;
  delete _domain;
  delete _particlePairsHandler;
  delete _moleculeContainer;
  delete _integrator;
  delete _inputReader;

}

void Simulation::output(unsigned long simstep)
{
  int ownrank = _domainDecomposition->getRank();
  global_log->set_mpi_output_root(0);

  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++)
  {
    (*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &(this->_lmu));
  }
  
  if(this->_doRecordRDF && !(simstep % this->_RDFOutputTimesteps))
  {
    this->_domain->collectRDF(this->_domainDecomposition);
    if(!ownrank)
    {
      this->_domain->accumulateRDF();
      for(unsigned i=0; i < this->_numberOfComponents; i++)
      {
        for(unsigned j=i; j < this->_numberOfComponents; j++)
        {
          ostringstream osstrm;
          osstrm << this->_RDFOutputPrefix << "_" << i << "-" << j << ".";
          osstrm.fill('0');
          osstrm.width(9);
          osstrm << right << simstep;
          this->_domain->outputRDF(osstrm.str().c_str(), i, j);
          osstrm.str(""); osstrm.clear();
        }
      }
    }
    this->_domain->resetRDF();
  }
  
  if((simstep >= this->_initStatistics) && this->_doRecordProfile
                                        && !(simstep % this->_profileRecordingTimesteps))
  {
    this->_domain->recordProfile(this->_moleculeContainer);
  }
  if((simstep >= this->_initStatistics) && this->_doRecordProfile 
                                        && !(simstep % this->_profileOutputTimesteps))
  {
    this->_domain->collectProfile(this->_domainDecomposition);
    if(!ownrank)
    {
      ostringstream osstrm;
      osstrm << this->_profileOutputPrefix << ".";
      osstrm.fill('0');
      osstrm.width(9);
      osstrm << right << simstep;
      this->_domain->outputProfile(osstrm.str().c_str());
      osstrm.str(""); osstrm.clear();
    }
    this->_domain->resetProfile();
  }
  
  if(!(simstep % this->_outputFrequency))
  {
     _moleculeContainer->deleteOuterParticles();
     ostringstream osstrm;
     osstrm.str("");
     osstrm << this->_outputPrefix;
     // osstrm.width(3);
     // osstrm.fill('0');
     // osstrm << this->_domain->ownrank();
     osstrm << ".restart.xdr";
     this->_domain->writeCheckpoint(
        osstrm.str(), _moleculeContainer, _domainDecomposition
     );
  }

  if(_domain->thermostatWarning())
      global_log->warning() << "Thermostat!" << endl;
  global_log->info() << "Simstep = "<< simstep
      << "\tT = " << _domain->T(0)
      << "\tU_pot = " << _domain->getAverageGlobalUpot() 
      << "\tp = " << _domain->getGlobalPressure() << endl;

#if 0 
  // CHECKFORREMOVAL
  if(ownrank==0)
  {
     double a = _domain->getDirectedVelocity(1);
     if(a > 0)
     {
        cout << "\t\t" << _domain->getDirectedVelocity(1, 2)
             << "\t" << _domain->getUniformAcceleration(1, 2)
             << "\t" << _domain->getUniformAcceleration(1)
             << "\t" << _domain->getCosetN(1);
     }
     cout << "\t";
     list<ChemicalPotential>::iterator cpit;
     for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++)
     {
        cout << (this->_domain->thermostatWarning()? "*": "")
             << cpit->getGlobalN() << "  ";
     }
     cout << (this->_domain->thermostatWarning()? "  *": "");
  }
#endif
}

void Simulation::updateParticleContainerAndDecomposition()
{

  // The particles have moved, so the neighbourhood relations have
  // changed and have to be adjusted
  _moleculeContainer->update();
  //_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
  _domainDecomposition->balanceAndExchange(true, _moleculeContainer, _domain->getComponents(), _domain, this->_cutoffRadius);
  // The cache of the molecules must be updated/build after the exchange process,
  // as the cache itself isn't transferred 
  Molecule* tM;
  for (tM = _moleculeContainer->begin(); tM
       != _moleculeContainer->end(); tM = _moleculeContainer->next()){
    tM->upd_cache();
  }
}

double Simulation::Tfactor(unsigned long simstep)
{
   double xi = (double)(simstep - this->_initSimulation) / (double)(this->_initCanonical - this->_initSimulation);
   if((xi < 0.1) || (xi > 0.9)) return 1.0;
   else if(xi < 0.3) return 15.0*xi - 0.5;
   else if(xi < 0.4) return 10.0 - 20.0*xi;
   else if(xi < 0.6) return 2.0;
   else return 4 - 10.0*xi/3.0;
}

