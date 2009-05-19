#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <fstream>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"

class ParticleContainer;
class DomainDecompBase; 

using namespace std;

//! @brief This class is used to read in the phasespace and to handle macroscopic values
//! @author Martin Bernreuther, Martin Buchholz
//!
//! The main purpose of this class is responsible for all macroscopic value (Potential,...).
//! It is important to differentiate between local and global version of those values
//! Macroscopic values are values that aggregate some property of a set of molecules.
//! As this program is designed to run on a parallel computer, there are typically 
//! several processes. Each process has an instance of this class, but with a different 
//! subset of molecules. Therefore, also the macroscopic values are only representative 
//! for the "local" domain of that process. Thus, local member variables that represent 
//! "local" macroscopic values, begin with _local...
//! 
//! At some points of the simulation, macroscopic values for the whole set of molecules 
//! have to be calculated. Those values are stored in member variables beginning with _global. 
class Domain{
 public:
  //! The constructor sets _localRank to rank and initializes all member variables 
  Domain(int rank);

  //! @brief reads in the data of all molecules
  //! 
  //! The Molecule Data starts in a new line with the string "MoleculeFormat"
  //! followed by whitespace and a string representing the format.
  //! In the standard case (format ICRVQD), the following values are provided for each molecule:
  //! \li id of the molecule (int)
  //! \li id of the component of the molecule (int)
  //! \li Coordinates: x, y, z (all double)
  //! \li velocities: vx, vy, vz (all double)
  //! \li Orientation (quaternion): q0, q1, q2, q3 (all double)
  //! \li Angular Momentum: Dx, Dy, Dz (all double)
  //! 
  //! An example can be seen in the documentation of this class
  //! @param particleContainer Here the Molecules from the input file are stored 
  //void readPhaseSpaceData(ParticleContainer* particleContainer);
  //void readPhaseSpaceData(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, ParticleGenerator *partGen);

  //! @brief writes a checkpoint file that can be used to continue the simulation
  //!
  //! The format of the checkpointfile written by this method is the same as the format
  //! of the input file. 
  //! @param filename Name of the checkpointfile (including path)
  //! @param particleContainer The molecules that have to be written to the file are stored here
  //! @param domainDecomp In the parallel version, the file has to be written by more than one process.
  //!                     Methods to achieve this are available in domainDecomp 
  void writeCheckpoint(string filename, ParticleContainer* particleContainer,
		       DomainDecompBase* domainDecomp) const;

  //! @brief initialize far field correction parameters
  //! 
  //! By limiting the calculation to pairs of particles which have
  //! less distance than the given cutoff radius, an error is made
  //! By calculating approximations for the neglected pairs, the
  //! error can be reduced
  //! @param cutoffRadius cutoff radius
  //! @todo How does the Correction work? Give reference to some paper,
  //!       documentation in the implementation
  void initFarFieldCorr(double cutoffRadius); // CHECKED

  //! @brief initialize parameter streams
  //!
  //! This method should only be called, after the component information
  //! and all molecule data have been read in
  //! @param cutoffRadius cutoff radius
  void initParameterStreams(double cutoffRadius);

  //! @brief set the potential of the local process
  void setLocalUpot(double Upot); 

  //! @brief get the potential of the local process
  double getLocalUpot() const;

  //! @brief set the virial of the local process
  void setLocalVirial(double Virial); 

  //! @brief get the virial of the local process
  double getLocalVirial() const;   
    
  //! @brief get thermostat scaling for translations 
  double getGlobalBetaTrans() const;

  //! @brief get thermostat scaling for rotations 
  double getGlobalBetaRot() const;

  //! @brief return the length of the domain
  //!
  //! @param index dimension for which the length should be returned
  double getGlobalLength(int index) const;
    
  //! @brief set the length of the domain
  //!
  //! @param index dimension for which the length should be set
  //! @param index value which should be set
  void setGlobalLength(int index, double length);

  //! @brief get the global temperature
  double getGlobalTemperature() const;

  //! @brief set the global temperature
  void setGlobalTemperature(double temp);

  //! @brief get the mixcoeff
  vector<double> & getmixcoeff();

  //! @brief get the epsilonRF
  double getepsilonRF() const;

  //! @brief set the epsilonRF
  void setepsilonRF(double erf);

  //! @brief get globalNumMolecules
  unsigned long getglobalNumMolecules() const;

  //! @brief set globalNumMolecules
  void setglobalNumMolecules(unsigned long glnummol);

  //! @brief get the global pressure
  //!
  //! @todo provide justification for the formula
  double getGlobalPressure() const;

  //! @brief get the global average potential per particle
  //!
  //! Before this method is called, it has to be sure that the
  //! global potential has been calculated (method calculateGlobalValues)
  double getAverageGlobalUpot() const;

  //! @brief get the global average virial per particle
  //!
  //! Before this method is called, it has to be sure that the
  //! global virial has been calculated (method calculateGlobalValues)
  double getAverageGlobalVirial() const;

  //! @brief sets _localSummv2 to the given value
  void setLocalSummv2(double summv2);
    
  //! @brief sets _localSumIw2 to the given value
  void setLocalSumIw2(double sumIw2);
	
  //! @brief get local rank 
  int getlocalRank();
    
  //! @brief get inpversion
  unsigned long getinpversion();

  //! @brief set inpversion
  void setinpversion(unsigned long inpv);

  //! @brief get globalRho
  double getglobalRho();

  //! @brief set globalRho
  void setglobalRho(double grho);

  //! @brief get globalRotDOF
  unsigned long getglobalRotDOF();

  //! @brief set globalRotDOF
  void setglobalRotDOF(unsigned long grotdof);

  //! @brief get the current time
  double getCurrentTime();
    
  //! @brief get the current time
  void setCurrentTime(double curtime);
    
  //! @brief advance the current time by timestep
  void advanceTime(double timestep);

  //! @brief get a reference to the vector of Components
  vector<Component>& getComponents();
    
  //! @brief get the parameter streams
  Comp2Param& getComp2Params();
       
  //! @brief calculate the global macroscopic values
  //!
  //! @param domainDecomp domain decomposition
  //! @param particleContainer particle Container
  //! @todo more detailled description
  void calculateGlobalValues(DomainDecompBase* domainDecomp,
			     ParticleContainer* particleContainer);
    
  //! @brief calculate _localSummv2 and _localSumIw2
  //! @todo more detailled description
  void calculateVelocitySums(ParticleContainer* partCont);
    
 private:
 
  //! rank of the local process
  int _localRank;

  //! @brief Version of the input file
  //!
  //! even though it is desirable, that the format of the input file
  //! doesn't change, is sometimes does change. When that happens,
  //! the code which reads in the input file (parser) has to be changed as well.
  //! old versions of the input file then can't be read any more.
  //! So whenever the parser is changed, _inpversion is set to the
  //! date of the change (YYYYMMDD) (hard-coded). Only input files
  //! with the same version are sure to be processed correctly
  unsigned long _inpversion;
    
  //! Potential of the local process
  double _localUpot;
  //! Virial of the local process
  double _localVirial;
  //! global Potential
  double _globalUpot;
  //! global virial
  double _globalVirial;
  //! @todo WHAT IS THIS? 
  double _globalBetaTrans;
  //! @todo WHAT IS THIS?
  double _globalBetaRot;
  //! global density
  double _globalRho;
  //! global Temperature
  double _globalTemperature;
  //! global Number of rotational degrees of freedom
  unsigned long _globalRotDOF;
  //! global Number of Molecules
  //! @todo redundancy?
  unsigned long _globalNumMolecules;
  //! side length of the cubic simulation box
  double _globalLength[3];

  //! local sum (over all molecules) of the mass multiplied with the squared velocity
  //! @todo rename this 
  double _localSummv2;
  //! local sum (over all molecules) of the moment of inertia 
  //! multiplied with the squared  rotational velocity
  //! @todo rename this
  double _localSumIw2; 
    
  //! reaction field
  //! @todo WHAT IS THIS?
  //! @todo local or global?
  double _epsilonRF;

  //! Potential correction for the error made by the cutoff
  //! @todo local or global?
  double _UpotCorr;
  //! Virial correction for the error made by the cutoff
  //! @todo local or global?
  double _VirialCorr;     

  //! @todo comment
  double _currentTime;

  //! Components resp. molecule types
  vector<Component> _components;
  //! parameter streams for each possible pair of molecule-types
  Comp2Param _comp2params;
  //! modified Lorentz-Berthelot mixing rule parameters
  //! @todo more explanation
  vector<double> _mixcoeff;

};


#endif /*DOMAIN_H_*/
