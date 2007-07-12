#include "Simulation.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "datastructures/LinkedCells.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"
#ifdef PARALLEL
#include "parallel/DomainDecomposition.h"
#endif
#include "datastructures/adapter/ParticlePairs2PotForceAdapter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "md_io/OutputBase.h"
#include "md_io/ResultWriter.h"
#include "md_io/XyzWriter.h"

#include <fstream>
#include <iostream>

utils::Log Simulation::_log("Simulation");



Simulation::Simulation(int *argc, char ***argv){
  #ifdef PARALLEL
    _domainDecomposition = (parallel::DomainDecompBase*) new parallel::DomainDecomposition(argc, argv);
  #else
    _domainDecomposition = (parallel::DomainDecompBase*) new parallel::DomainDecompDummy();
  #endif
  int ownrank = _domainDecomposition->getRank();

  // two parameters are necessary: config-file and number of timesteps
  // if any other number of parameters is provied, print usage message and exit
  if(*argc != 4) {
    if(ownrank == 0) {
      cout << "Usage: " << *argv[0] << " <configfilename> <number of timesteps> <outputprefix>" << endl;
    }
    delete _domainDecomposition;
    exit(1);
  }
  
  // open filestream to the input file
  string inputfilename((*argv)[1]);  
  ifstream inputfilestream(inputfilename.c_str());
  
  // store number of timesteps to be simulated
  _numberOfTimesteps = atol((*argv)[2]);

  // store prefix for output files
  _outputPrefix = string((*argv)[3]);

  // used to store one token of the inputfilestream
  string token;
  string phaseSpaceFileName;
  double timestepLength;
  
  _domain = new Domain(ownrank);
  _particlePairsHandler = new datastructures::ParticlePairs2PotForceAdapter(*_domain);
  
  // The first line of the config file has to contain the token "MDProjectConfig"
  inputfilestream >> token;
  if(token != "MDProjectConfig") {
    cerr << "Not a MDProject config file! " << token << endl;
    exit(1);
  }   
  while(inputfilestream){
    token.clear();
    inputfilestream >> token;

    if(token.substr(0,1)=="#"){
      inputfilestream.ignore(INT_MAX,'\n');
    }
    else if(token=="phaseSpaceFile"){      
      inputfilestream >> phaseSpaceFileName;
      _domain->setPhaseSpaceFile(phaseSpaceFileName);
      _domain->readPhaseSpaceHeader();
      _domain->initParameterStreams(_cutoffRadius);
    }
    else if(token=="timestepLength") {
      inputfilestream >> timestepLength;
    }
    else if(token=="cutoffRadius")   {
      inputfilestream >> _cutoffRadius;
    } 
    else if(token=="datastructure")  {
      inputfilestream >> token;
      if(token=="LinkedCells"){

        int cellsInCutoffRadius;
        inputfilestream >> cellsInCutoffRadius;
        double bBoxMin[3];
        double bBoxMax[3];
        for(int i=0;i<3;i++) {
          bBoxMin[i] = _domainDecomposition->getCoords(i)*_domain->getGlobalLength(i)/_domainDecomposition->getGridSize(i);
          bBoxMax[i] = (_domainDecomposition->getCoords(i)+1)*_domain->getGlobalLength(i)/_domainDecomposition->getGridSize(i);
        }
        
        _moleculeContainer = new datastructures::LinkedCells<Molecule>(bBoxMin, bBoxMax,
                                                             _cutoffRadius, cellsInCutoffRadius, *_particlePairsHandler);
        
      }
    }
    else if(token=="output")  {
      inputfilestream >> token;
      if(token=="ResultWriter"){
        string outputPathAndPrefix;
        inputfilestream >> outputPathAndPrefix;
        _outputPlugins.push_back(new md_io::ResultWriter(outputPathAndPrefix));
      }
      else if(token=="XyzWriter"){
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new md_io::XyzWriter(_numberOfTimesteps, writeFrequency,outputPathAndPrefix));
      }
    }
  }
  
  
  _domain->readPhaseSpaceData(_moleculeContainer);
  _domain->initFarFieldCorr(_cutoffRadius);
  
  // @todo comment
  _integrator = new integrators::Leapfrog(timestepLength);
  

}


void Simulation::initialize(){

  // clear halo
  _moleculeContainer->deleteOuterParticles();
  
  updateParticleContainerAndDecomposition();

  // Force calculation
  _moleculeContainer->traversePairs();
 
  // clear halo
  _moleculeContainer->deleteOuterParticles();

  //! @todo calculation of macroscopic values, so that the output in
  //!       step 1 is correct. This doesn't work yet as for the methode
  //!       _domain->calculateGlobalValues(...), the iterator has
  //!       to be executed before (sets summv2 and sumIw2)
  _domain->calculateVelocitySums(_moleculeContainer);
  
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer);

  // initialize output
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain); 
  }

}


void Simulation::simulate(){
  
  Molecule* tempMolecule;
  int ownrank = _domainDecomposition->getRank();
  if(ownrank==0) _log.info("simulation(...)", "Starting Simulation: ");
  //double runclock = clock();

  initialize();
  // MAIN LOOP
  for(unsigned long simstep=1; simstep<=_numberOfTimesteps; simstep++){
    
    _integrator->eventNewTimestep(_moleculeContainer, _domain);
    
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
    for(tempMolecule = _moleculeContainer->begin(); tempMolecule != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
      tempMolecule->scale_v(_domain->getGlobalBetaTrans());
      tempMolecule->scale_D(_domain->getGlobalBetaRot());      
    }
    
    _domain->advanceTime(_integrator->getTimestepLength());
    
    output(simstep);
  }
  
  string cpfile(_outputPrefix+".restart.inp");
  _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
  
  
  // finish output
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain); 
  }
  
  
  delete _domainDecomposition;
  delete _domain;
  // wait for all processes to reach the end of the program
  //_DomainDecomposition->barrier();
}

void Simulation::output(int simstep){
  int ownrank = _domainDecomposition->getRank();
  
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep); 
  }

  if(ownrank==0) cout << simstep << "\t" << _domain->getAverageGlobalUpot() << "\t" << _domain->getGlobalPressure()  << "\t" << endl; //double(runclock)/CLOCKS_PER_SEC << " s" << endl;

}

void Simulation::updateParticleContainerAndDecomposition(){
  
  _domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);

  // The cache of the molecules must be updated/build after the exchange process,
  // as the cache itself isn't transferred
  Molecule* tempMolecule;
  for(tempMolecule = _moleculeContainer->begin(); tempMolecule != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
    tempMolecule->upd_cache();
  }

  // The particles have moved, so the neighbourhood relations have
  // changed and have to be adjusted
  _moleculeContainer->update();
  
}
