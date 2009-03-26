#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <string>
#include <list>

using namespace std;

class Domain;
class ParticleContainer;
class ParticlePairsHandler;
class Integrator;
class OutputBase;
class DomainDecompBase;
class InputBase;


using namespace std;

//! @brief controls the whole simulation process
//! @author Martin Buchholz
//!
//! Some of the simulation parameters are provided in a config file.
//! The order of the parameters in the config file is important.
//! Thats's because e.g. the datastructure can only be built after the 
//! phasespace has been read. \n
//! The config file starts with a line containing the token "MDProjectConfig"
//! followed by the following parameters (possibly mixed with comment lines starting with "#"):
//! - timestepLength: Uses by the Iterator to calculate new velocities and positions
//! - cutoffRadius: Determines the maximum distance for which the force between two 
//!                   molecules still has to be calculated
//! - phaseSpaceFile: Full path to the file containing the phase space
//! - parallelisation: Parallelisation scheme to be used
//!                    - DomainDecomposition: standard spacial domain decomposition into
//!                                           cuboid regions of equal size
//! - datastructure: Datastructure to be used (e.g. Linked Cells) followed by
//!                    the parameters for the datastructures
//!                  The datastructure LinkedCells needs one additional parameter, 
//!                  which is the number of cells in the cutoff radius (equals to the 
//!                  cutoff radius divided by the cell length). 
//!  
//! Example for the config file: \n
//! 
//! \pre
//! MDProjectConfig 
//! timestepLength 0.00005 
//! cutoffRadius 3.0 
//! phaseSpaceFile /home/buchholm/MDProject/phasespace/phasespace.dat 
//! # datastructure followed by the parameters for the datastructure 
//! # for LinkedCells, the cellsInCutoffRadius has to be provided 
//! datastructure LinkedCells 2 
//! parallelisation DomainDecomposition
//! \endpre
class Simulation{
 public:
  //! @brief instantiate simulation object
  //!
  //! The Constructor opens the file with the given filename and reads in
  //! all parameters for the simulaion and initializes the following member variables:
  //! - timestepLength: 
  //! - cutoffRadius
  //! - phaseSpace
  //! - moleculeContainer
  //! @param argc Pointer to the number of arguments passed to the programm. Needed for MPI
  //! @param argv the arguments, also needed for MPI
  Simulation(int *argc, char ***argv);

  //! @brief calculate all values for the starting timepoint
  //! 
  //! After the input file has been read in, only the information which is
  //! directly in the inputfile is known at time step zero.
  //! This includes e.g. Molecule position and velocities, but not forces
  //! which act on the molecules or macroscopic values at the potential.
  //! This method has to ensure, that all values are available at
  //! time step zero
  void initialize();

  //! @brief Controls the main loop of the simulation.
  //!
  //! precondition for this method is that initialize() has been called.
  //! The main loop calls all methods which have to be called in each
  //! iteration step, e.g. initializes the molecule Container, calculates
  //! the forces, ... \n
  //! For the integration, a seperate integration object is used.
  //! This method is written in a way that should make it possible to use
  //! different integrators without having to change anything.
  //! Whenever something is done which might make it necessary for an integrator
  //! to do something, the integrator is informed about that using
  //! the integrators corresponding method (see integrator documentation)
  //! It follows a corse outline of what has to be done:
  //! - Do all output that has to be done in each time step
  //! - Inform the iterator that all information is available (new timestep)
  //! - As the position of the molecules changed, the domain decomposition and
  //!     the datastructure have to be updated
  //! - calculate forces (and afterwards delete halo)
  //! - Inform the iterator that forces for the new positions have been calculated 
  //! - The iterator should now have finished everything to be done in this time step,
  //!     so the macroscopic values can be calculated
  //! - velocity and angular momentum have to be scaled
  void simulate();
 
  //! @brief output results 
  //! @param simstep timestep of the output
  //! @todo comment
  void output(int simstep);
    
  //! The following things have to be done here:
  //! - bring all molecules to the corresponding processes (including copies for halo)
  //! - update the caches of the molecules
  //! - update the ParticleContainer
  void updateParticleContainerAndDecomposition();
      
    
 private:
    
  //! maximum distance at which the forces between two molecules still have to be calculated.
  double _cutoffRadius;

  //! Number of discrete time steps for the simulation        
  unsigned long _numberOfTimesteps;

  //! Incremental output flag NEW
  bool _increment;

  //! Datastructure for finding neighbours efficiently
  ParticleContainer* _moleculeContainer;
    
  //! Handler describing what action is to be done for each particle pair
  ParticlePairsHandler* _particlePairsHandler;
    
  //! module which handles the domain decomposition
  DomainDecompBase* _domainDecomposition;

  //! numerical solver for the particles equations of motion
  Integrator* _integrator;
    
  //! all macroscopic (local and global) information
  Domain* _domain;

  //! responsible for reading in the phasespace (header+data)
  InputBase* _inputReader;

  //! prefix for the names of all output files
  string _outputPrefix;
    
  //! frequency of the checkpoint writer
  long _outputFrequency;
    
  //! list of output plugins to use
  std::list<OutputBase*> _outputPlugins;
    
    
};
#endif /*SIMULATION_H_*/
