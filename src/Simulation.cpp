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

#include <fstream>
#include <iostream>
#include <iterator>

Simulation::Simulation(int *argc, char ***argv)
{
  int ownrank = 0;
#ifdef PARALLEL
  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownrank);
#endif
  if (*argc != 4) {
    if (ownrank == 0) {
      cout << "Usage: " << *argv[0] << " <configfilename> <number of timesteps> <outputprefix>" << endl;
    }
    exit(1);
  }

#ifndef PARALLEL
  _domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
#endif

  // open filestream to the input file
  string inputfilename((*argv)[1]);
  ifstream inputfilestream(inputfilename.c_str());
  if ( !inputfilestream.is_open() ) {
    cerr << "ERROR: Could not open file " << inputfilename << endl;
    exit (1);
  }

  std::string inputPath;
  unsigned int lastIndex = inputfilename.find_last_of('/',inputfilename.size()-1);
  if (lastIndex == string::npos) {
      inputPath="";
  }    
  else {
    inputPath = inputfilename.substr(0, lastIndex+1);
  }    
  

  // store number of timesteps to be simulated
  _numberOfTimesteps = atol((*argv)[2]);

  // store prefix for output files
  _outputPrefix = string((*argv)[3]);

  // used to store one token of the inputfilestream
  string token;

  double timestepLength;

  _domain = new Domain(ownrank);
  _particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

  // The first line of the config file has to contain the token "MDProjectConfig"
  inputfilestream >> token;
  if (token != "MDProjectConfig") {
    cerr << "Not a MDProject config file! " << token << endl;
    exit(1);
  }
  
  while (inputfilestream) {
    token.clear();
    inputfilestream >> token;

    if (token.substr(0, 1)=="#"){
      inputfilestream.ignore(1024,'\n');
    } 
    else if (token=="phaseSpaceFile") {
      string phaseSpaceFileFormat; 
      inputfilestream >> phaseSpaceFileFormat;

      if(phaseSpaceFileFormat=="OldStyle"){
        string phaseSpaceFileName;
        inputfilestream >> phaseSpaceFileName;
        _inputReader = (InputBase*) new InputOldstyle();
        _inputReader->setPhaseSpaceFile(phaseSpaceFileName);
        _inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
        _inputReader->readPhaseSpaceHeader(_domain);
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
        ((PartGen*)_inputReader)->readPhaseSpaceHeader(_domain);
      }
      _domain->initParameterStreams(_cutoffRadius);
    } 
    else if (token=="timestepLength") {
      inputfilestream >> timestepLength;
    } 
    else if (token=="cutoffRadius") {
      inputfilestream >> _cutoffRadius;
    } 
    else if (token=="parallelisation") {
#ifndef PARALLEL
      cerr << "\nWARING: Input file demands parallelisation, but the current compilation doesn't\n\tsupport parallel execution.\n" << endl;
      inputfilestream >> token;
      //exit(1);
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
    }
  }
  
  // read particle data
  _inputReader->readPhaseSpace(_moleculeContainer, _domain, _domainDecomposition);

  _domain->initFarFieldCorr(_cutoffRadius);
  

  // @todo comment
  _integrator = new Leapfrog(timestepLength);

  // test new Decomposition
  _moleculeContainer->update();
  _moleculeContainer->deleteOuterParticles();

}

void Simulation::initialize()
{
  // clear halo
  _moleculeContainer->deleteOuterParticles();

  updateParticleContainerAndDecomposition();
  // Force calculation

  _moleculeContainer->traversePairs();

  // clear halo
  _moleculeContainer->deleteOuterParticles();

  _domain->calculateVelocitySums(_moleculeContainer);

  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer);

  // initialize output 
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain);
  }
}

void Simulation::simulate()
{

  Molecule* tempMolecule;

  initialize();

  // MAIN LOOP
  for (unsigned long simstep=1; simstep<=_numberOfTimesteps; simstep++) {
    _integrator->eventNewTimestep(_moleculeContainer, _domain);

    // ensure that all Particles are in the right cells and exchange Particles
    updateParticleContainerAndDecomposition();

    // Force calculation
    _moleculeContainer->traversePairs();

    // clear halo
    _moleculeContainer->deleteOuterParticles();

    // Inform the integrator about the calculated forces
    _integrator->eventForcesCalculated(_moleculeContainer, _domain);

    // calculate the global macroscopic values from the local values
    _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer);

    // scale velocity and angular momentum
    // @todo why here? what about preF
    for (tempMolecule = _moleculeContainer->begin(); tempMolecule
         != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
      tempMolecule->scale_v(_domain->getGlobalBetaTrans());
      tempMolecule->scale_D(_domain->getGlobalBetaRot());
    }

    _domain->advanceTime(_integrator->getTimestepLength());
    output(simstep);
  }
  string cpfile(_outputPrefix+".restart.inp");
  _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);

  // finish output
  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
      (*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
    delete (*outputIter);
  }

  delete _domainDecomposition;
  delete _domain;
  delete _particlePairsHandler;
  delete _moleculeContainer;
  delete _integrator;
  delete _inputReader;

}

void Simulation::output(int simstep)
{
  int ownrank = _domainDecomposition->getRank();

  std::list<OutputBase*>::iterator outputIter;
  for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
    (*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep);
  }

  if (ownrank==0)
    cout << simstep << "\t" << _domain->getAverageGlobalUpot() << "\t"
         << _domain->getGlobalPressure() << "\t" << endl; //double(runclock)/CLOCKS_PER_SEC << " s" << endl;

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
  Molecule* tempMolecule;
  for (tempMolecule = _moleculeContainer->begin(); tempMolecule
       != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
    tempMolecule->upd_cache();
  }
}
