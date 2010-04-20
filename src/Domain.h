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

#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <fstream>
#include <map>
#include <queue>
#include <list>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"

class Molecule;
class ParticleContainer;
class DomainDecompBase; 

using namespace std;

//! @brief This class is used to read in the phasespace and to handle macroscopic values
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
//!
//! This class is responsible for all macroscopic values.
//! It is important to differentiate between local and global version of those values
//! Macroscopic values are values that aggregate some property of a set of molecules.
//! As this program is designed to run on a parallel computer, there are typically 
//! several processes. Each process has an instance of this class, but with a different 
//! subset of molecules. Therefore, also the macroscopic values are only representative 
//! for the "local" domain of that process. 
//!
//! member variables that represent "local" macroscopic values begin with _local
//!
//! properties of the global system begin with _global if only the rank 0 process
//! obtains the correct value (e.g., as a sum over corresponding local properties)
//! and with _universal if the global value must be communicated to all processes.
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
		       DomainDecompBase* domainDecomp);

  //! @brief initialize far field correction parameters
  //! 
  //! By limiting the calculation to pairs of particles which have
  //! less distance than the given cutoff radius, an error is made
  //! By calculating approximations for the neglected pairs, the
  //! error can be reduced.
  //! The way the cutoff correction works is essentially explained
  //! in the lecture notes "Molekulare Thermodynamik" by Vrabec for
  //! the LJ potential, and it is analogous for electrostatics.
  //! Probably all books on molecular simulation explain that method,
  //! e.g. have a look at Allen/Tildesley.
  //!
  //! The idea is that for radii above the cutoff radius, the
  //! radial distribution function is assumed to be exactly one, so
  //! that the correction for the internal energy as well as the
  //! virial can be obtained immediately from an integral over the
  //! pair potential and its derivative, respectively.
  //! For quadrupoles, the long-range correction turns out to be zero
  //! so that it is not calculated here. Molecules with point charges
  //! and a zero net charge interact over long distances according to
  //! their dipole moments, which is why effective dipoles are
  //! calculated for the far field correction.
  //!
  //! If ions appear, so that the net charge of a molecule is distinct
  //! from zero, the far field terms become much more complicated, and
  //! THAT IS NOT IMPLEMENTED AT PRESENT.
  //!
  //! @param cutoffRadius cutoff radius for electrostatics
  //! @param cutoffRadiusLJ cutoff radius for the LJ potential
  void initFarFieldCorr(double cutoffRadius, double cutoffRadiusLJ); 

  //! @brief initialize parameter streams
  //!
  //! This method should only be called, after the component information
  //! and all molecule data have been read in
  //! @param cutoffRadius cutoff radius
  void initParameterStreams(double cutoffRadius, double cutoffRadiusLJ);

  //! @brief set the potential of the local process
  void setLocalUpot(double Upot); 

  //! @brief get the potential of the local process
  double getLocalUpot() const;

  //! @brief set the virial of the local process
  void setLocalVirial(double Virial); 

  //! @brief get the virial of the local process
  double getLocalVirial() const;   
    
  //! @brief get thermostat scaling for translations 
  double getGlobalBetaTrans();
  double getGlobalBetaTrans(int thermostat);

  //! @brief get thermostat scaling for rotations 
  double getGlobalBetaRot();
  double getGlobalBetaRot(int thermostat);

  //! @brief return the length of the domain
  //!
  //! @param index dimension for which the length should be returned
  double getGlobalLength(int index) const;
    
  //! @brief set the length of the domain
  //!
  //! @param index dimension for which the length should be set
  //! @param index value which should be set
  void setGlobalLength(int index, double length);

  //! @brief get the global temperature for the whole system (i.e. thermostat ID 0)
  double getGlobalCurrentTemperature() { return this->getCurrentTemperature(0); }
  double getCurrentTemperature(int thermostat) { return this->_globalTemperatureMap[thermostat]; }
  double getTargetTemperature(int thermostat) { return this->_universalTargetTemperature[thermostat]; }

  //! @brief set the global temperature
  void setGlobalTemperature(double T);
  void setTargetTemperature(int thermostat, double T);

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
  double getGlobalPressure();

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

  //! @brief sets _localSummv2 to the given value (for thermostat ID 0)
  void setLocalSummv2(double summv2);
  //! @brief sets _localSummv2 to the given value
  void setLocalSummv2(double summv2, int thermostat);
    
  //! @brief sets _localSumIw2 to the given value (for thermostat ID 0)
  void setLocalSumIw2(double sumIw2);
  //! @brief sets _localSumIw2 to the given value
  void setLocalSumIw2(double sumIw2, int thermostat);
	
    //! @brief sets _localThermostatN and _localRotationalDOF for thermostat
    void setLocalNrotDOF(int thermostat, unsigned long N, unsigned long rotDOF)
    {
       this->_localThermostatN[thermostat] = N;
       this->_localRotationalDOF[thermostat] = rotDOF;
    }
    unsigned getComponentRotDOF(int cid) { return this->_components[cid].rot_dof(); }
  
  //! @brief get local rank 
  int getlocalRank();
    
  //! @brief get input version
  unsigned long getinpversion();

  //! @brief set input version
  void setinpversion(unsigned long inputVersion);

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

  //! @brief get a reference to the vector of components
  vector<Component>& getComponents();

  //! @brief add a component to the vector of components
  void addComponent(Component component);

  //! @brief get the parameter streams
  Comp2Param& getComp2Params();
       
  //! @brief calculate the global macroscopic values
  //!
  //! @param domainDecomp domain decomposition
  //! @param particleContainer particle Container
  //! @param collectThermostatVelocities flag stating whether the directed velocity should be collected for the corresponding thermostats
  //! @param Tfactor temporary factor applied to the temperature during equilibration
  //!
  //! Essentially, this method calculates all thermophysical values
  //! that require communication between the subdomains of the system.
  //!
  //! In particular, it determines:
  //!   - the potential energy and the virial
  //!   - if (collectThermostatVelocities == true), the directed velocity associated with the appropriately marked thermostats
  //!   - the kinetic energy and number of DOF associated with each thermostat
  //!   - on that basis, the correction factors for translation and rotation, also applying some heuristics in case of an explosion
  //!
  //! The heuristics applied in case of an explosion were already commented
  //! on, in the appropriate place (Domain.cpp), but apparently it should be
  //! repeated here (according to the requirements communicated by TUM):
  //! Sometimes an "explosion" occurs, for instance, by inserting a particle
  //! at an unfavorable position, or due to imprecise integration of the
  //! equations of motion. The thermostat would then remove the kinetic energy from
  //! the rest of the system, effectively destroying a simulation run that could
  //! be saved by a more intelligent version, which is provided by the present
  //! heuristics. It is essential that such an intervention does not occur
  //! regularly, therefore it is limited to three times every 3000 time steps.
  //! That limit is implemented using the properties _universalSelectiveThermostatCounter,
  //! _universalSelectiveThermostatWarning, and _universalSelectiveThermostatError.
  void calculateGlobalValues(
     DomainDecompBase* domainDecomp, ParticleContainer* particleContainer,
     bool collectThermostatVelocities, double Tfactor
  );
  //! @brief calls this->calculateGlobalValues with Tfactor = 1 and without velocity collection
  void calculateGlobalValues(DomainDecompBase* domainDecomp,
			     ParticleContainer* particleContainer)
  {
     this->calculateGlobalValues(domainDecomp, particleContainer, false, 1.0);
  }
  
  //! @brief calculate _localSummv2 and _localSumIw2
  //!
  //! The present method calculates the translational and rotational
  //! kinetic energies associated with the thermostats of the system,
  //! together with the corresponding number of DOF. If a thermostat
  //! is marked as applying to the undirected motion only, this is
  //! explicitly considered by subtracting the average directed
  //! motion corresponding to the thermostat from the motion of each
  //! individual particle, for the purpose of aggregating the kinetic
  //! energy of the undirected motion.
  void calculateVelocitySums(ParticleContainer* partCont);
  
    //! Calculates the LOCAL directed kinetic energy associated
    //! with the appropriately marked thermostats
    void calculateThermostatDirectedVelocity(ParticleContainer* partCont);
    //! A thermostat is referred to as undirected here if it
    //! explicitly excludes the kinetic energy associated with
    //! the directed motion of the respective components.
    bool thermostatIsUndirected(int thermostat) { return this->_universalUndirectedThermostat[thermostat]; }
    //! @brief returns the directed velocity associated with a thermostat
    //! @param th ID of the thermostat
    //! @param d coordinate 0 (v_x), 1 (v_y), or 2 (v_z).
    //!
    //! It should be obvious that this method only returns sensible values
    //! for thermostats marked as "undirected", because otherwise the
    //! directed velocity is not explicitly computed.
    double getThermostatDirectedVelocity(int thermostat, int d) { return this->_universalThermostatDirectedVelocity[d][thermostat]; }

    int ownrank() { return this->_localRank; }

    //! @brief returns whether there are several distinct thermostats in the system
    bool severalThermostats() { return this->_universalComponentwiseThermostat; }
    //! @brief thermostat to be applied to component cid
    //! @param cid ID of the respective component
    int getThermostat(int cid) { return this->_universalThermostatID[cid]; }
    //! @brief disables the componentwise thermostat (so that a single thermostat is applied to all DOF)
    void disableComponentwiseThermostat() { this->_universalComponentwiseThermostat = false; }
    //! @brief enables the componentwise thermostat
    void enableComponentwiseThermostat();
    //! @brief returns the ID of the "last" thermostat in the system
    //!
    //! The idea is that ID -1 refers to DOF without thermostats,
    //! while ID 0 refers to the entire system (or to the unique thermostat if the componentwise
    //! thermostat is disabled). In case of componentwise thermostats being applied,
    //! these have the IDs from 1 to this->maxThermostat().
    unsigned maxThermostat()
    {
       return (_universalComponentwiseThermostat)? (_universalThermostatN.size() - 2): 0;
    }
    //! @brief associates a component with a thermostat
    //! @param cid internal ID of the component
    //! @param th internal ID of the thermostat
    void setComponentThermostat(int cid, int thermostat)
    {
       if((0 > cid) || (0 >= thermostat)) exit(787);
       this->_universalThermostatID[cid] = thermostat;
       this->_universalThermostatN[thermostat] = 0;
    }
    //! @brief enables the "undirected" flag for the specified thermostat
    //! @param tst internal ID of the respective thermostat
    void enableUndirectedThermostat(int thermostat);

    //! @brief assigns a coset ID to a component (ID)
    void assignCoset(unsigned cid, unsigned cosetid) { _universalComponentSetID[cid] = cosetid; }
    //! @brief sets the information on the acceleration model for one coset
    void specifyComponentSet(unsigned cosetid, double v[3], double tau, double ainit[3], double timestep);
    //! @brief sets the number of timesteps between two updates of the uniform acceleration
    void setUCAT(unsigned uCAT) { this->_universalConstantAccelerationTimesteps = uCAT; }
    //! @brief returns the number of timesteps between two updates of the uniform acceleration
    unsigned getUCAT() { return this->_universalConstantAccelerationTimesteps; }
    //! @brief are there any cosets?
    bool isAcceleratingUniformly() { return ( (this->_universalTau.size() > 0)
                                              && (this->_universalConstantAccelerationTimesteps > 0) ); }
    //! @brief updates the intensity and direction of the uniform acceleration
    void determineAdditionalAcceleration
    (
        DomainDecompBase* domainDecomp,
        ParticleContainer* molCont, double dtConstantAcc
    );
    //! @brief returns the acceleration map (necessary for passing data to the integrator)
    map<unsigned, double>* getUAA() { return this->_universalAdditionalAcceleration; }
    //! @brief returns the cosetid of a component (0 for unaccelerated components)
    unsigned getComponentSet(unsigned cid)
    {
       if(_universalComponentSetID.find(cid) == _universalComponentSetID.end()) return 0;
       else return this->_universalComponentSetID[cid];
    }
    //! @brief returns the absolute directed velocity for a component set
    //! @param cosetid ID of the component set
    double getDirectedVelocity(unsigned cosetid);
    //! @brief returns the directed velocity for a component set
    //! @param cosetid ID of the component set
    //! @param d x direction (0), y direction (1), or z direction (2)
    double getDirectedVelocity(unsigned cosetid, unsigned d);
    //! @brief returns the absolute external acceleration for a component set
    //! @param cosetid ID of the component set
    double getUniformAcceleration(unsigned cosetid);
    //! @brief returns the external acceleration for a component set
    //! @param cosetid ID of the component set
    //! @param d x direction (0), y direction (1), or z direction (2)
    double getUniformAcceleration(unsigned cosetid, unsigned d);
    //! @brief returns the difference between the desired velocity and the global average velocity
    double getMissingVelocity(unsigned cosetid, unsigned d);
    //! @brief total number of particles that belong to the specified component set
    double getCosetN(unsigned cosetid) { return this->_globalN[cosetid]; }

    void setupProfile(unsigned xun, unsigned yun, unsigned zun);
    void considerComponentInProfile(int cid); 
    void recordProfile(ParticleContainer* molCont);
    void collectProfile(DomainDecompBase* domainDecomp);
    void outputProfile(const char* prefix);
    void resetProfile();

    unsigned maxCoset() { return this->_universalTau.size(); }

    double N() {return this->_globalNumMolecules;}
    double N(unsigned cid) { return this->_components[cid].numMolecules(); }
    void Nadd(unsigned cid, int N, int localN);

   double getGlobalLength(int d) { return _globalLength[d]; }
   double getGlobalVolume() { return (_globalLength[0] *  _globalLength[1] *  _globalLength[2]); }

   void setupRDF(double interval, unsigned bins);
   void resetRDF();
   void collectRDF(DomainDecompBase* domainDecomp);
   void outputRDF(const char* prefix, unsigned i, unsigned j);
   void accumulateRDF();
   void tickRDF() { this->_universalRDFTimesteps++; }
   inline void observeRDF(unsigned i) { this->_localCtr[i] ++; }
   inline void observeRDF(double dd, unsigned i, unsigned j)
   {
      if(this->_universalRDFTimesteps < 0) return;
      if(dd > this->ddmax) return;
      if(i > j) { this->observeRDF(dd, j, i); return; }
      unsigned l = (unsigned)floor(sqrt(dd)/this->_universalInterval);
      this->_localDistribution[i][j-i][l] ++;
   }
   void thermostatOff() { this->_universalNVE = true; }
   void thermostatOn() { this->_universalNVE = false; }
   bool NVE() { return this->_universalNVE; }
   bool thermostatWarning() { return (this->_universalSelectiveThermostatWarning > 0); }

   void evaluateRho(unsigned long localN, DomainDecompBase* comm);

   void init_cv(unsigned N, double U, double UU)
   {
      this->_globalUSteps = N;
      this->_globalSigmaU = U;
      this->_globalSigmaUU = UU;
   }
   void record_cv();
   double cv();

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
  //! global density
  double _globalRho;
  //! global Number of Molecules
  //! @todo redundancy?
  unsigned long _globalNumMolecules;
  //! side length of the cubic simulation box
  double _globalLength[3];

    //! does a componentwise thermostat apply?
    bool _universalComponentwiseThermostat;
    //! thermostat IDs. negative: no thermostat, 0: global, positive: componentwise
    //! in the case of a componentwise thermostat, all components are assigned
    //! a thermostat ID different from zero.
    map<int, int> _universalThermostatID;
    //! _localThermostatN[0] and _universalThermostatN[0] are always the total number
    //! of particles in the subdomain and, respectively, the entire domain
    map<int, unsigned long> _localThermostatN;
    map<int, unsigned long> _universalThermostatN;
    map<int, unsigned long> _localRotationalDOF;
    map<int, unsigned long> _universalRotationalDOF;
    //! _globalTemperatureMap[0] is always the temperature of the whole system,
    //! including components to which no thermostat is applied.
    //! The _globalTemperatureMap stores actual CURRENT temperatures, whereas
    //! the temperature objective of the thermostat is stored in _universalTargetTemperature
    map<int, double> _globalTemperatureMap;
    map<int, double> _universalTargetTemperature;
    map<int, double> _universalBTrans;
    map<int, double> _universalBRot;
    //! should the directed movement be subtracted when calculating the temperature?
    map<int, bool> _universalUndirectedThermostat;
    //! stores the velocity of the directed movement
    map<int, double> _universalThermostatDirectedVelocity[3];
    map<int, double> _localThermostatDirectedVelocity[3];
    bool _universalNVE;

    /// calculate new value of the uniform acceleration each # timesteps
    unsigned _universalConstantAccelerationTimesteps;
    /// assigns a component set ID to some of the components
    map<unsigned, unsigned> _universalComponentSetID;
    /// local number of molecules that belong to a given component set ID
    map<unsigned, unsigned long> _localN;
    /// global number of molecules that belong to a given component set ID
    map<unsigned, unsigned long> _globalN;
    /// local sum of the velocity vectors corresponding to a given component set ID
    map<unsigned, long double> _localVelocitySum[3];
    /// global sum of the velocity vectors corresponding to a given component set ID
    map<unsigned, long double> _globalVelocitySum[3];
    /// uniform acceleration
    map<unsigned, double> _universalAdditionalAcceleration[3];
    /// target average velocity for the molecules of a coset
    map<unsigned, double> _globalTargetVelocity[3];
    /// delay variable tau of a coset
    map<unsigned, double> _universalTau;
    /// queue of previously recorded velocity sums
    map<unsigned, deque<long double> > _globalPriorVelocitySums[3];
    /// number of items in the velocity queue
    map<unsigned, unsigned> _globalVelocityQueuelength;

    //! computation of the isochoric heat capacity
    unsigned _globalUSteps;
    double _globalSigmaU;
    double _globalSigmaUU;    

    //! 1 / dimension of a profile cuboid
    double _universalInvProfileUnit[3];
    //! number of successive profile cuboids in x/y/z direction
    unsigned _universalNProfileUnits[3];
    //! local N profile map
    map<unsigned, long double> _localNProfile;
    //! global N profile map
    map<unsigned, double> _globalNProfile;
    //! local directed velocity profile map
    map<unsigned, long double> _localvProfile[3];
    //! global directed velocity  profile map
    map<unsigned, double> _globalvProfile[3];
    //! local kinetic energy profile map
    map<unsigned, long double> _localKineticProfile;
    //! global kinetic energy profile map
    map<unsigned, double> _globalKineticProfile;
    //! local counter w. r. t. degrees of freedom
    map<unsigned, long double> _localDOFProfile;
    //! global counter w. r. t. degrees of freedom
    map<unsigned, double> _globalDOFProfile;
    //! how many _evaluated_ timesteps are currently accumulated in the profile?
    unsigned _globalAccumulatedDatasets;
    //! which components should be considered?
    map<unsigned, bool> _universalProfiledComponents;

    bool _doCollectRDF;
    double _universalInterval;
    unsigned _universalBins;
    int _universalRDFTimesteps, _universalAccumulatedTimesteps;
    double ddmax;
    unsigned long *_localCtr, *_globalCtr, *_globalAccumulatedCtr;
    unsigned long ***_localDistribution, ***_globalDistribution, ***_globalAccumulatedDistribution;

    int _universalSelectiveThermostatCounter;
    int _universalSelectiveThermostatWarning;
    int _universalSelectiveThermostatError;

    //! local sum (over all molecules) of the mass multiplied with the squared velocity
    map<int, double> _local2KETrans;
    //! local sum (over all molecules) of the moment of inertia 
    //! multiplied with the squared  rotational velocity
    map<int, double> _local2KERot; 
    
  //! reaction field
  //!
  //! This is neither "local" nor "global" but a parameter of the reaction field method.
  //! (So one might regard it as "global" formally.)
  //! For a description of the reaction field method cf. the dissertation of J. Stoll.
  //! It was introduced by J. A. Barker and R. O. Watts, Mol. Phys. 26: 789-792 (1973).
  double _epsilonRF;

  //! Global potential correction for the error made by the cutoff
  double _UpotCorr;
  //! Global virial correction for the error made by the cutoff
  double _VirialCorr;     

  //! Contains the time t in reduced units
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
