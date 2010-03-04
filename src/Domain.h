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
//! @author Martin Bernreuther, Martin Buchholz, et al. (2009)
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

  //! @brief get the global temperature
  double getGlobalTemperature() { return this->T(0); }
  double T(int thermostat) { return this->_globalTemperatureMap[thermostat]; }
  double targetT(int thermostat) { return this->_universalTargetTemperature[thermostat]; }

  //! @brief set the global temperature
  void setGlobalTemperature(double temp);
  void setTargetT(int thermostat, double T);

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

  //! @brief sets _localSummv2 to the given value
  void setLocalSummv2(double summv2);
  void setLocalSummv2(double summv2, int thermostat);
    
  //! @brief sets _localSumIw2 to the given value
  void setLocalSumIw2(double sumIw2);
  void setLocalSumIw2(double sumIw2, int thermostat);
	
    //! @brief sets _localThermostatN(i) and _localRotationalDOF(i)
    void setLocalNrotDOF(int i, unsigned long N, unsigned long rotDOF)
    {
       this->_localThermostatN[i] = N;
       this->_localRotationalDOF[i] = rotDOF;
    }
    unsigned getComponentRotDOF(int cid) { return this->_components[cid].rot_dof(); }
  
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
  //! @todo more detailled description
  void calculateGlobalValues(DomainDecompBase* domainDecomp,
			     ParticleContainer* particleContainer)
  {
     this->calculateGlobalValues(domainDecomp, particleContainer, false, 1.0);
  }
  void calculateGlobalValues(
     DomainDecompBase* domainDecomp, ParticleContainer* particleContainer,
     bool collectThermostatVelocities, double Tfactor
  );
  
  //! @brief calculate _localSummv2 and _localSumIw2
  //! @todo more detailled description
  void calculateVelocitySums(ParticleContainer* partCont);
  
    void calculateThermostatDirectedVelocity(ParticleContainer* partCont);
    bool thermostatIsUndirected(int th) { return this->_universalUndirectedThermostat[th]; }
    double thermostatv(int th, int d) { return this->_universalThermostatDirectedVelocity[d][th]; }

    int ownrank() { return this->_localRank; }

    bool severalThermostats() { return this->_universalComponentwiseThermostat; }
    int getThermostat(int cid) { return this->_universalThermostatID[cid]; }
    void disableCT() { this->_universalComponentwiseThermostat = false; }
    void enableCT();
    unsigned maxThermostat()
    {
       return (_universalComponentwiseThermostat)? (_universalThermostatN.size() - 2): 0;
    }
    void setComponentThermostat(int cid, int th)
    {
       if((0 > cid) || (0 >= th)) exit(787);
       this->_universalThermostatID[cid] = th;
       this->_universalThermostatN[th] = 0;
    }
    void enableUndirectedThermostat(int tst);

    /// assigns a coset ID to a component (ID)
    void assignCoset(unsigned cid, unsigned cosetid) { _universalComponentSetID[cid] = cosetid; }
    /// sets the information on the acceleration model for one coset
    void specifyComponentSet(unsigned cosetid, double v[3], double tau, double ainit[3], double timestep);
    /// sets the number of timesteps between two updates of the uniform acceleration
    void setUCAT(unsigned uCAT) { this->_universalConstantAccelerationTimesteps = uCAT; }
    /// returns the number of timesteps between two updates of the uniform acceleration
    unsigned getUCAT() { return this->_universalConstantAccelerationTimesteps; }
    /// are there any cosets?
    bool isAcceleratingUniformly() { return ( (this->_universalTau.size() > 0)
                                              && (this->_universalConstantAccelerationTimesteps > 0) ); }
    /// updates the intensity and direction of the uniform acceleration
    void determineAdditionalAcceleration
    (
        DomainDecompBase* domainDecomp,
        ParticleContainer* molCont, double dtConstantAcc
    );
    /// returns the acceleration map (necessary for passing data to the integrator)
    map<unsigned, double>* getUAA() { return this->_universalAdditionalAcceleration; }
    /// returns the cosetid of a component (0 for unaccelerated components)
    unsigned getComponentSet(unsigned cid)
    {
       if(_universalComponentSetID.find(cid) == _universalComponentSetID.end()) return 0;
       else return this->_universalComponentSetID[cid];
    }
    double getDirectedVelocity(unsigned cosetid);
    double getDirectedVelocity(unsigned cosetid, unsigned d);
    double getUniformAcceleration(unsigned cosetid);
    double getUniformAcceleration(unsigned cosetid, unsigned d);
    /// returns the difference between the desired velocity and the global average velocity
    double getMissingVelocity(unsigned cosetid, unsigned d);
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
      if(dd > this->ddmax) return;
      if(i > j) { this->observeRDF(dd, j, i); return; }
      unsigned l = (unsigned)(sqrt(dd)/this->_universalInterval);
      this->_localDistribution[i][j-i][l] ++;
   }
   void thermostatOff() { this->_universalNVE = true; }
   void thermostatOn() { this->_universalNVE = false; }
   bool NVE() { return this->_universalNVE; }
   bool thermostatWarning() { return (this->_universalSelectiveThermostatWarning > 0); }

   void evaluateRho(unsigned long localN, DomainDecompBase* comm);

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
    unsigned _universalRDFTimesteps, _universalAccumulatedTimesteps;
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
