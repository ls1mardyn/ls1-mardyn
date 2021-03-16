
#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <map>
#include <array>
#include <cstdint>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include "utils/CommVar.h"
/* 
 * TODO add comments for variables 
 */
#define CHECKPOINT_FILE_VERSION 20160512  /**< checkpoint file version */

#define MIN_BETA 0.9  /**< minimal scaling factor before an explosion is detected */
#define KINLIMIT_PER_T 10.0

#include "molecules/MoleculeForwardDeclaration.h"
class ParticleContainer;
class DomainDecompBase;
class XMLfileUnits;

//! @brief This class is used to read in the phasespace and to handle macroscopic values
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2011)
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
class Domain {

private:
	Domain();
	Domain(Domain &domain);

	Domain& operator=(Domain &domain);
	
public:
	//! The constructor sets _localRank to rank and initializes all member variables
	Domain(int rank);

	void readXML(XMLfileUnits& xmlconfig);
	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param currentTime The current time to be printed.
	//! @param useBinaryFormat indicates wheter binary I/O is used or not
	void writeCheckpoint( std::string filename, ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, double currentTime, bool useBinaryFormat = false);

	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the header file (including path). Dependent on the I/O format, this could
	//! 				be the same file as the data file
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param currentTime The current time to be printed.
	void writeCheckpointHeader(std::string filename,
			ParticleContainer* particleContainer,
			const DomainDecompBase* domainDecomp, double currentTime);

	void writeCheckpointHeaderXML(std::string filename,
			ParticleContainer* particleContainer,
			const DomainDecompBase* domainDecomp, double currentTime);

	//! @brief initialize far field correction parameters
	//!
	//! By limiting the calculation to pairs of particles which have
	//! less distance than the given cutoff radius, an error is made
	//! By calculating approximations for the neglected pairs, the
	//! error can be reduced.
	//! The way the cutoff correction works is probably explained
	//! by all books on molecular simulation, e.g. Allen/Tildesley [1].
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
	//! [1] Computer Simulation of Liquids (1989), Clarendon, Oxford.
	//!
	//! @param cutoffRadius cutoff radius for electrostatics
	//! @param cutoffRadiusLJ cutoff radius for the LJ potential
//	void initFarFieldCorr(double cutoffRadius, double cutoffRadiusLJ);

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
	
	//! @brief set the fluid and fluid-solid potential of the local process
	void setLocalUpotCompSpecific(double UpotCspec);

	//! @brief set the number of fluid phase components (specified in the config-file)
	void setNumFluidComponents(unsigned nc);
	
	//! @brief get the numbr of fluid molecules as specified in the config file (*_1R.cfg)
	unsigned getNumFluidComponents();
	
	//! @brief get the fluid and fluid-solid potential of the local process
	double getLocalUpotCompSpecific();

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
	//! @param d  dimension for which the length should be returned
	double getGlobalLength(int d) const { return _globalLength[d]; }

	//! @brief set the length of the domain
	//!
	//! @param index dimension for which the length should be set
	//! @param length value which should be set
	void setGlobalLength(int index, double length);

	//! @brief get the global temperature for the whole system (i.e. thermostat ID 0)
	double getGlobalCurrentTemperature() { return getCurrentTemperature(0); }
	double getCurrentTemperature(int thermostatID) { return _globalTemperatureMap[thermostatID]; }
	double getTargetTemperature(int thermostatID) { return _universalTargetTemperature[thermostatID]; }

	//! @brief set the global temperature
	void setGlobalTemperature(double T) { setTargetTemperature(0, T); }
	void setTargetTemperature(int thermostatID, double T);

	//! @brief get the mixcoeff
	std::vector<double> & getmixcoeff();

	//! @brief get the epsilonRF
	double getepsilonRF() const;

	//! @brief set the epsilonRF
	void setepsilonRF(double erf);

	//! @brief get globalNumMolecules
	unsigned long getglobalNumMolecules() const;

	//! @brief set globalNumMolecules
	void setglobalNumMolecules(unsigned long glnummol);
	
	//! @brief set numMolecules in FixRegion 

    void setNumFixRegion(unsigned long nf);
        
    //! @brief set numMolecules in FixRegion 
    unsigned long getNumFixRegion();

	//! @brief update globalNumMolecules
	void updateglobalNumMolecules(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

	//! @brief get local/global max. moleculeID
	CommVar<uint64_t> getMaxMoleculeID() const;

	//! @brief update max. moleculeID
	void updateMaxMoleculeID(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

	//! @brief get the global pressure
	double getGlobalPressure();

	//! @brief get the global average potential per particle
	//!
	//! Before this method is called, it has to be sure that the
	//! global potential has been calculated (method calculateGlobalValues)
	double getAverageGlobalUpot() const;
        double getGlobalUpot() const;

	//! by Stefan Becker: return the average global potential of the fluid-fluid and fluid-solid interaction (but NOT solid-solid interaction)
	double getAverageGlobalUpotCSpec();

	//! by Stefan Becker: determine and return the totel number of fluid molecules
	//! this method assumes all molecules with a component-ID less than _numFluidComponent to be fluid molecules 
	unsigned long getNumFluidMolecules();

	//! @brief get the global average virial per particle
	//!
	//! Before this method is called, it has to be sure that the
	//! global virial has been calculated (method calculateGlobalValues)
	double getAverageGlobalVirial() const;

	//! @brief sets _localSummv2 to the given value
	void setLocalSummv2(double summv2, int thermostat);

	//! @brief sets _localSumIw2 to the given value
	void setLocalSumIw2(double sumIw2, int thermostat);

	//! @brief sets _localThermostatN and _localRotationalDOF for thermostat
	void setLocalNrotDOF(int thermostat, unsigned long N, unsigned long rotDOF) {
		this->_localThermostatN[thermostat] = N;
		this->_localRotationalDOF[thermostat] = rotDOF;
	}
	
	//! @brief get globalRho
	double getglobalRho();

	//! @brief set globalRho
	void setglobalRho(double grho);

	//! @brief get globalRotDOF
	unsigned long getglobalRotDOF();

	//! @brief set globalRotDOF
	void setglobalRotDOF(unsigned long grotdof);

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

	/* FIXME: alternatively: default values for function parameters */
	//! @brief calls this->calculateGlobalValues with Tfactor = 1 and without velocity collection
	void calculateGlobalValues(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer) {
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
	double getThermostatDirectedVelocity(int thermostat, int d) { return this->_universalThermostatDirectedVelocity[thermostat][d]; }

	//! @brief returns whether there are several distinct thermostats in the system
	bool severalThermostats() { return this->_componentwiseThermostat; }
	//! @brief thermostat to be applied to component cid
	//! @param cid ID of the respective component
	int getThermostat(int cid) { return this->_componentToThermostatIdMap[cid]; }
	//! @brief disables the componentwise thermostat (so that a single thermostat is applied to all DOF)
	void disableComponentwiseThermostat() { this->_componentwiseThermostat = false; }
	//! @brief enables the componentwise thermostat
	void enableComponentwiseThermostat();
	//! @brief returns the ID of the "last" thermostat in the system
	//!
	//! The idea is that ID -1 refers to DOF without thermostats,
	//! while ID 0 refers to the entire system (or to the unique thermostat if the componentwise
	//! thermostat is disabled). In case of componentwise thermostats being applied,
	//! these have the IDs from 1 to this->maxThermostat().
	unsigned maxThermostat() {
		return (_componentwiseThermostat)? (_universalThermostatN.size() - 2): 0;
	}
	//! @brief associates a component with a thermostat
	//! @param cid internal ID of the component
	//! @param th internal ID of the thermostat
	void setComponentThermostat(int cid, int thermostat) {
		if((0 > cid) || (0 >= thermostat)) Simulation::exit(787);
		this->_componentToThermostatIdMap[cid] = thermostat;
		this->_universalThermostatN[thermostat] = 0;
	}
	//! @brief enables the "undirected" flag for the specified thermostat
	//! @param tst internal ID of the respective thermostat
	void enableUndirectedThermostat(int thermostat);

	unsigned long N() {return _globalNumMolecules;}

	void Nadd(unsigned cid, int N, int localN);

	double getGlobalVolume() const { return (_globalLength[0] *  _globalLength[1] *  _globalLength[2]); }

	void thermostatOff() { this->_universalNVE = true; }
	void thermostatOn() { this->_universalNVE = false; }
	bool NVE() { return this->_universalNVE; }
	bool thermostatWarning() { return (this->_universalSelectiveThermostatWarning > 0); }

	void evaluateRho(unsigned long localN, DomainDecompBase* comm);
	void submitDU(unsigned cid, double DU, double* r);
	void setLambda(double lambda) { this->_universalLambda = lambda; }
        void setDensityCoefficient(float coeff) { _globalDecisiveDensity = coeff; }
        void setProfiledComponentMass(double m) { _universalProfiledComponentMass = m; }

	void init_cv(unsigned N, double U, double UU) {
		this->_globalUSteps = N;
		this->_globalSigmaU = U;
		this->_globalSigmaUU = UU;
	}
	void record_cv();
	double cv();

    // by Stefan Becker <stefan.becker@mv.uni-kl.de>
	/* method returning the sigma parameter of a component 
	=> needed in the output of the MmspdWriter (specifying the particles' radii in a movie) */
	double getSigma(unsigned cid, unsigned nthSigma);
	// needed for the MmspdWriter (MegaMol)
	unsigned getNumberOfComponents();
	
	void setUpotCorr(double upotcorr){ _UpotCorr = upotcorr; }
	void setVirialCorr(double virialcorr){ _VirialCorr = virialcorr; }

    // explosion heuristics, NOTE: turn off when using slab thermostat
    void SetExplosionHeuristics(bool bVal) { _bDoExplosionHeuristics = bVal; }

private:

	//! rank of the local process
	int _localRank;

	//! Potential of the local process
	double _localUpot;
	//! by Stefan Becker: component specific potential of the local process (fluid-fluid and fluid-solid, but not solid-solid)
	double _localUpotCspecif;
	//! by Stefan Becker: number of fluid components specified in the config file 
	//! --> used to determine the average component specific potential energy
	unsigned _numFluidComponent;

	//! Virial of the local process
	double _localVirial;
	//! global Potential
	double _globalUpot;
	//! global component specific potential (fluid-fluid and fluid-solid but NOT solid-solid)
	double _globalUpotCspecif;
	//! global virial
	double _globalVirial;
	//! global density
	double _globalRho;
	//! global Number of Molecules
	//! @todo redundancy?
	unsigned long _globalNumMolecules;
	CommVar<uint64_t> _maxMoleculeID;
	//! side length of the cubic simulation box
	double _globalLength[3];

	//! does a componentwise thermostat apply?
	bool _componentwiseThermostat;
	//! thermostat IDs. negative: no thermostat, 0: global, positive: componentwise
	//! in the case of a componentwise thermostat, all components are assigned
	//! a thermostat ID different from zero.
	std::map<int, int> _componentToThermostatIdMap;

	//! _localThermostatN[0] and _universalThermostatN[0] are always the total number
	//! of particles in the subdomain and, respectively, the entire domain
	std::map<int, unsigned long> _localThermostatN;
	std::map<int, unsigned long> _universalThermostatN;
	std::map<int, unsigned long> _localRotationalDOF;
	std::map<int, unsigned long> _universalRotationalDOF;
	//! _globalTemperatureMap[0] is always the temperature of the whole system,
	//! including components to which no thermostat is applied.
	//! The _globalTemperatureMap stores actual CURRENT temperatures, whereas
	//! the temperature objective of the thermostat is stored in _universalTargetTemperature
	std::map<int, double> _globalTemperatureMap;
	std::map<int, double> _universalTargetTemperature;
	std::map<int, double> _universalBTrans;
	std::map<int, double> _universalBRot;
	//! should the directed movement be subtracted when calculating the temperature?
	std::map<int, bool> _universalUndirectedThermostat;
	//! stores the velocity of the directed movement
	std::map<int, std::array<double, 3> > _universalThermostatDirectedVelocity;
	std::map<int, std::array<double, 3> > _localThermostatDirectedVelocity;

	/* FIXME: This info should go into an ensemble class */
	bool _universalNVE;

	//! computation of the isochoric heat capacity
	unsigned _globalUSteps;
	double _globalSigmaU;
	double _globalSigmaUU;
	//! which components should be considered?
	std::map<unsigned, bool> _universalProfiledComponents;
        double _universalProfiledComponentMass;  // set from outside
        double _universalLambda;  // set from outside
        float _globalDecisiveDensity;  // set from outside

	int _universalSelectiveThermostatCounter;
	int _universalSelectiveThermostatWarning;
	int _universalSelectiveThermostatError;

	//! local sum (over all molecules) of the mass multiplied with the squared velocity
	std::map<int, double> _local2KETrans;
	//! local sum (over all molecules) of the moment of inertia
	//! multiplied with the squared  rotational velocity
	std::map<int, double> _local2KERot; 
	
	// number of molecules in FixRegion
	unsigned long _numMoleculesFix;

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


	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;
	//! modified Lorentz-Berthelot mixing rule parameters
	//! @todo more explanation
	std::vector<double> _mixcoeff;

    // explosion heuristics, NOTE: turn off when using slab thermostat
    bool _bDoExplosionHeuristics;
};


#endif /*DOMAIN_H_*/
