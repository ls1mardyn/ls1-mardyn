/*************************************************************************
 * Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <string>
#include <list>
//#include <vector>
//#include <sstream>

#include "ensemble/GrandCanonical.h"
#ifdef STEEREO
  #include <simSteering.h>
#endif

class Domain;
class ParticleContainer;
class ParticlePairsHandler;
class Integrator;
class OutputBase;
class DomainDecompBase;
class InputBase;
class Molecule;
namespace optparse {
  class OptionParser;
  class Values;
}

//! @brief controls the whole simulation process
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
//!
//! Some of the simulation parameters are provided in a config file.
//! The order of the parameters in the config file is important.
//! Thats's because e.g. the datastructure can only be built after the 
//! phasespace has been read. \n
//!
//! The config file usually has the file ending ".txt" and
//! starts with a line containing the token "mardynconfig"
//! followed by the following parameters, among others
//! (possibly mixed with comment lines starting with "#"):
//! - timestepLength: Uses by the Iterator to calculate new velocities and positions
//! - cutoffRadius: Determines the maximum distance for which the force between two 
//!                   molecules still has to be calculated
//! - phaseSpaceFile: Full path to the XDR file containing the phase space
//! - parallelization: Parallelisation scheme to be used
//!                    - DomainDecomposition: standard spacial domain decomposition into
//!                                           cuboid regions of equal size
//! - datastructure: Datastructure to be used (e.g. Linked Cells) followed by
//!                    the parameters for the datastructures
//!                  The datastructure LinkedCells needs one additional parameter, 
//!                  which is the number of cells in the cutoff radius (equals to the 
//!                  cutoff radius divided by the cell length). 
//!  
//! Example for a config file: \n
//! 
//! \pre
//! mardynconfig 
//! timestepLength 0.00005 
//! cutoffRadius 3.0 
//! phaseSpaceFile OldStype phasespace.xdr
//! # datastructure followed by the parameters for the datastructure 
//! # for LinkedCells, the cellsInCutoffRadius has to be provided 
//! datastructure LinkedCells 1 
//! parallelization DomainDecomposition
//! \endpre
class Simulation{
 public:
  //! @brief instantiate simulation object
  //!
  //! The Constructor processes the command line arguments
  //! @param argc Pointer to the number of arguments passed to the programm. Needed for MPI
  //! @param argv Pointer to the list of arguments, also needed for MPI
  Simulation(int *argc, char ***argv);

  //! @brief destruct simulation object
  //!
  ~Simulation();

  //! @brief Terminate simulation with a given exit code.
  //!
  //! The exit method takes care over the right way to terminate the application in a correct way
  //! for the different parallelization schemes. e.g. terminating other processes in MPI parallel
  //! execution mode.
  int exit( int exitcode );

  const optparse::Values& initOptions(int argc, char *argv[], optparse::OptionParser& op);

  //! @brief process configuration file
  //! 
  //! calls initConfigXML or initConfigOldStyle
  //! @param inputfilename filename of the input file
  void initConfigFile(const std::string& inputfilename);
  void initConfigFile(const char* inputfilename) { initConfigFile(std::string(inputfilename)); }

  //! @brief process XML configuration file (*.xml)
  //! 
  //! Opens the XML file with the given filename and reads in all parameters
  //! for the simulaion and initializes the following member variables:
  //! - timestepLength: 
  //! - cutoffRadius
  //! - phaseSpace
  //! - moleculeContainer
  //! @param inputfilename filename of the XML input file
  void initConfigXML(const std::string& inputfilename);
  void initConfigXML(const char* inputfilename) { initConfigXML(std::string(inputfilename)); }

  //! @brief process oldstyle configuration file (*.cfg)
  //! 
  //! Opens the file with the given filename and reads in all parameters
  //! for the simulaion and initializes the following member variables:
  //! - timestepLength: 
  //! - cutoffRadius
  //! - phaseSpace
  //! - moleculeContainer
  void initConfigOldstyle(const std::string& inputfilename);
  void initConfigOldstyle(const char* inputfilename) { initConfigOldstyle(std::string(inputfilename)); }

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
  //! It follows a coarse outline of what has to be done:
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
  //! 
  //! The present method serves as a redirection to the actual output.
  //! That includes a) particular output objects included in _outputPlugins,
  //! b) conventional output methods, i.e. outputRDF for the radial distribution
  //! function and writeCheckpoint for XDR (configuration) files.
  void output(unsigned long simstep);

  //! The following things have to be done here:
  //! - bring all molecules to the corresponding processes (including copies for halo)
  //! - update the caches of the molecules
  //! - update the ParticleContainer
  void updateParticleContainerAndDecomposition();

  void setDomainDecomposition (DomainDecompBase* ddb) {_domainDecomposition = ddb;};
  Domain* getDomain () {return _domain;};
  ParticleContainer* getMolecules () {return _moleculeContainer;};
  unsigned long getSimStep () {return _simstep;};
  double getcutoffRadius() const {return _cutoffRadius;}

  //! @brief temperature increase factor during automatic equilibration
  //! @param current simulation time step
  //!
  //! The present version of Mardyn provides the option of
  //! equilibrating the system at an increased temperature. This has,
  //! of course, the advantage that the equilibration is accelerated.
  //! In a special case that is of particular interest to a relevant
  //! subset of the developers, namely MD simulation of nucleation,
  //! it also reduces the amount of clusters formed during
  //! equilibration. Then, the actual nucleation process can be
  //! investigated for the equilibrated system.
  //!
  //! Equilibration at an increased temperature is nothing new, in
  //! principle it could be considered a variant of simulated
  //! annealing. For the purpose of accelerating the relaxation,
  //! starting from an unrealistic initial state, it was e.g.
  //! proposed in the thesis of M. Kreitmeir.
  //!
  //! The key parameter to regulating the equilibration is
  //! this->_initCanonical. The system is gradually heated and cooled
  //! again, and from this->_initCanonical on the temperature is
  //! constant (i.e. this method returns 1.0).
  //!
  //! Thus, if temperature increase is to be avoided, the user should
  //! set the "initCanonical" parameter to "0".
  double Tfactor(unsigned long simstep);
    
 private:
    
  //! maximum distance at which the forces between two molecules still have to be calculated.
  double _cutoffRadius;

  //! LJ cutoff (may be smaller than the RDF/electrostatics cutoff)
  double _LJCutoffRadius;

  //! external cutoff radius for the Tersoff potential
  double _tersoffCutoffRadius;
  
  //! flag specifying whether planar interface profiles are recorded
  bool _doRecordProfile;
  //! Interval between two evaluations of the profile.
  //! This means that only 1 / _profileRecordingTimesteps of the
  //! internally available data are actually used, so if precision is
  //! a concern, set the value to 1. On the other hand, the program
  //! may be accelerated somewhat by increasing the interval.
  unsigned _profileRecordingTimesteps;
  //! Aggregation interval for the profile data, i.e. if _profileRecordingTimesteps
  //! is 100 and _profileOutputTimesteps is 20 000, this means that
  //! the profiles found in the output are averages over 200 configurations.
  unsigned _profileOutputTimesteps;
  //! Although the meaning of this should be obvious, it may be noted
  //! that the time step and "rhpry" (density), "vzpry" (z-velocity),
  //! and Tpry (kinetic energy) will be attached to the prefix for
  //! the different profiles.
  std::string _profileOutputPrefix;
  //! flag specifying whether the radial distribution function is recorded
  bool _doRecordRDF;
  //! aggregation interval for the RDF data
  unsigned _RDFOutputTimesteps;
  //! the timestep, the respective component IDs, and "rdf" are
  //! appended to this prefix
  std::string _RDFOutputPrefix;
  //! controls the interval between two lines in the most important
  //! output file generated by Mardyn, i.e. the res file.
  unsigned _resultOutputTimesteps;
  
  //! A thermostat can be specified to account for the directed
  //! motion, which means that only the undirected kinetic energy is
  //! maintained constant. In many cases, this is an essential to
  //! simulating the system in the desired way, most notably for
  //! applying MD to nanoscopic fluid dynamics.
  //!
  //! Obviously, that requires calculating the directed velocity.
  //! Since the velocity is often fairly constant over time, it is
  //! not necessary to redetermine it every time the thermostat is
  //! applied. The property this->_collectThermostatDirectedVelocity
  //! regulates how often (in number of time steps) the directed
  //! velocity is evaluated.
  unsigned _collectThermostatDirectedVelocity;

  //! Sometimes during equilibration, a solid wall surrounded by
  //! liquid may experience a stress or an excessive pressure, which
  //! could damage its structure. With the flag this->_zoscillation,
  //! the z coordinate for some of the solid atoms (i.e. those which
  //! include Tersoff sites) is fixed, so that no motion in z
  //! direction can occur for the wall structure.
  bool _zoscillation;
  //! The fixed z coordinate applies to 1 out of this->_zoscillator
  //! atoms for all solid components.
  unsigned _zoscillator;
  
  //! Number of discrete time steps for the simulation        
  unsigned long _numberOfTimesteps;
  
  //! Incremental output flag NEW
  bool _increment;

  //! initial number of steps
  unsigned long _initSimulation;
  //! step number for the end of the configurational equilibration
  unsigned long _initCanonical;
  //! step number for activation of the muVT ensemble
  unsigned long _initGrandCanonical;
  //! step number for activation of all sorts of statistics
  unsigned long _initStatistics;

  // this should be obvious, right?
  unsigned _numberOfComponents;
    
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
  std::string _outputPrefix;
    
  //! frequency of the checkpoint writer
  long _outputFrequency;
  // unsigned _restartOutputInterval;
    
  //! list of output plugins to use
  std::list<OutputBase*> _outputPlugins;
  
  //! List of ChemicalPotential objects, each of which describes a
  //! particular control volume for the grand canonical ensemble with
  //! respect to one of the simulated components.
  //!
  //! It may at first be unclear why one could want to specify
  //! several grand canonical ensembles, which are then stored in a
  //! list. However, note that for every component a distinct
  //! chemical potential can be specified, and this is of course
  //! essential in certain cases. Also, different chemical potentials
  //! can be specified for different control volumes to induce a
  //! gradient of the chemical potential.
  std::list<ChemicalPotential> _lmu;
    
  //! This is Planck's constant. (Required for the Metropolis
  //! criterion which is used for the grand canonical ensemble).
  //! Of course, what is actually specified here is not the value
  //! of h in and of itself, since that is a universal constant, but
  //! the value of h AS EXPRESSED IN REDUCED UNITS, i.e. for the
  //! internal use of the program.
  double h;

#ifdef STEEREO
  SimSteering* _steer;
#endif

  //! simulation time step
  unsigned long _simstep;
};
#endif /*SIMULATION_H_*/

