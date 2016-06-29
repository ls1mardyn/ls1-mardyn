
#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <map>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"

#include "ensemble/GrandCanonical.h"
/* 
 * TODO add comments for variables 
 */
#define CHECKPOINT_FILE_VERSION 20140131  /**< checkpoint file version */

#define MIN_BETA 0.9  /**< minimal scaling factor before an explosion is detected */
#define KINLIMIT_PER_T 10.0

class Molecule;
class ParticleContainer;
class DomainDecompBase; 
class PressureGradient;
class XMLfileUnits;

using namespace std;

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
	Domain(int rank, PressureGradient* pg);

	void readXML(XMLfileUnits& xmlconfig);
	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	void writeCheckpoint( std::string filename, ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp );

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
        void setLocalVirial(double VIx, double VIy, double VIz, double VIILL, double VIILM);

	//! @brief get the virial of the local process
	double getLocalVirial() const;

	//! @brief get thermostat scaling for translations
	double getGlobalBetaTrans();
	double getGlobalBetaTrans(int thermostat);
	
	//! @brief get thermostat scaling for translations for one-dimesional thermostat
	double getGlobalAlphaTrans();
	double getGlobalAlphaTrans(int thermostat);


	//! @brief get thermostat scaling for rotations
	double getGlobalBetaRot();
	double getGlobalBetaRot(int thermostat);
	
	//! @brief get thermostat scaling for rotations
	long double getThT_heatFlux(int thermostat) { return _universalThT_heatFlux[thermostat]; }
	void setThT_heatFlux(int thermostat, double value) { _universalThT_heatFlux[thermostat] = value; }
	
	double getGlobalVirial(int d) { return this->_globalVirialI[d]; }
        double getGlobalVirialIILL() { return this->_globalVirialIILL; }
        double getGlobalVirialIILM() { return this->_globalVirialIILM; }
	
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
	// index = 0: startTime; index = 1: endTime
	unsigned long get_thermostatTimeSlot(int thermostat, int index) {return this->_thermostatTimeSlot[index][thermostat]; }

	//! @brief set the global temperature
	void setGlobalTemperature(double T);
	void setTargetTemperature(int thermostat, double T);
	void set1DimThermostat(int thermostat, int dimension);
	void setThermostatTimeSlot(int thermostat_id, unsigned long startTime, unsigned long endTime);

	//! @brief get the mixcoeff
	std::vector<double> & getmixcoeff();

	//! @brief get the epsilonRF
	double getepsilonRF() const;

	//! @brief set the epsilonRF
	void setepsilonRF(double erf);

	//! @brief get globalNumMolecules
	unsigned long getglobalNumMolecules() const;
	unsigned long getglobalOrigNumMolecules() const { return _globalOrigNumMolecules; }

	//! @brief set globalNumMolecules
	void setglobalNumMolecules(unsigned long glnummol);
	void setglobalOrigNumMolecules(unsigned long glnummol) { _globalOrigNumMolecules = glnummol; }

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
	void setLocalSummv2(double summv2, int thermostat);
	
	//! @brief sets localSummv2_1Dim to the given value
	void setLocalSummv2_1Dim(double summv2_1Dim, int thermostat) { this->_local2KETrans_1Dim[thermostat] = summv2_1Dim; }

	//! @brief sets _localSumIw2 to the given value
	void setLocalSumIw2(double sumIw2, int thermostat);

	//! @brief sets _localThermostatN and _localRotationalDOF for thermostat
	void setLocalNrotDOF(int thermostat, unsigned long N, unsigned long rotDOF) {
		this->_localThermostatN[thermostat] = N;
		this->_localRotationalDOF[thermostat] = rotDOF;
	}
	
	int getUniversalNProfileUnits_Stress(int d) {return this->_universalNProfileUnits_Stress[d]; }
	double getUniversalInvProfileUnit_Stress(int d) {return this->_universalInvProfileUnit_Stress[d]; }
	int getUniversalNProfileUnitsConfinement(int d) {return this->_universalNProfileUnitsConfinement[d]; }
	double getUniversalInvProfileUnitConfinement(int d) {return this->_universalInvProfileUnitConfinement[d]; }
	
	int getUniversalNProfileUnitsStressConfinement(int d) {return this->_universalNProfileUnitsStressConfinement[d]; }
	double getUniversalInvProfileUnitStressConfinement(int d) {return this->_universalInvProfileUnitStressConfinement[d]; }

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
	//! @brief calls this->calculate_maxSlabDist2GlobalValues with Tfactor = 1 and without velocity collection
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
	double getThermostatDirectedVelocity(int thermostat, int d) { return this->_universalThermostatDirectedVelocity[d][thermostat]; }
	
	void enable1DimThermostat(int thermostat) { this->_scale_v_1Dim[thermostat] = true; }
	std::map<int, int> getDim() { return this->_dimToThermostat; }
	void setDim(int thermostat, int dim) { this->_dimToThermostat[thermostat] = dim; }
	//! one-dimensional thermostat applied?
	bool isScaling1Dim(int thermostat) { return this->_scale_v_1Dim[thermostat]; }
	bool getAlphaTransCorrection(int thermostat) { return this->_alphaTransCorrection[thermostat]; }
	
	void enableStressCalculation(int component) { this->_stressCalc[component] = true; }
	bool isStressCalculating(int component) {return this->_stressCalc[component]; }
	
	bool isBulkPressure(int component) {return this->_bulkComponent[component]; }
	double getBulkBoundary(int dim) {return this->_bulkCorner[dim]; }
	
	// barostat
	bool isBarostat(int component) {return this->_barostatComponent[component]; }
	double getControl_top(int d) {return this->_control_top[d]; }
	double getControl_bottom(int d) {return this->_control_bottom[d]; }
	
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
		if((0 > cid) || (0 >= thermostat)) exit(787);
		this->_componentToThermostatIdMap[cid] = thermostat;
		this->_universalThermostatN[thermostat] = 0;
	}
	//! @brief enables the "undirected" flag for the specified thermostat
	//! @param tst internal ID of the respective thermostat
	void enableUndirectedThermostat(int thermostat);

	PressureGradient* getPG() { return this->_universalPG; }

	void setupProfile(unsigned xun, unsigned yun, unsigned zun);
	void setupSlabProfile(unsigned xun, unsigned yun, unsigned zun);
	void setupStressProfile(unsigned xun, unsigned yun, unsigned zun);
	void setupBulkPressure(double xmin, double xmax, double ymin, double ymax, unsigned cid);
	void setupConfinementProperties(double wallThickness, double horDist, double vertDist, double radius2, int cid, double xmax, double ymax, double zmax, unsigned long upperID, unsigned long lowerID);
	void setupConfinementProfile(unsigned xun,unsigned yun, double correlationLength);
	void considerComponentInProfile(int cid) { this->_universalProfiledComponents[cid] = true; }
	void considerComponentInProfileSlab(int cid) { this->_universalProfiledComponentsSlab[cid] = true; }
	void recordProfile(ParticleContainer* molCont);
	void recordSlabProfile(ParticleContainer* molCont);
	void recordStressProfile(ParticleContainer* molCont);
	void recordBulkPressure(ParticleContainer* molCont);
	void recordConfinementProperties(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned long simstep, unsigned long initStatistics);
	void collectProfile(DomainDecompBase* domainDecomp);
	void collectSlabProfile(DomainDecompBase* domainDecomp);
	void collectStressProfile(DomainDecompBase* domainDecomp);
	void collectBulkPressure(DomainDecompBase* dode);
	void collectConfinementProperties(DomainDecompBase* dode);
	void outputProfile(const char* prefix);
	void outputSlabProfile(const char* prefix);
	void outputStressProfile(const char* prefix);
	void outputConfinementProperties(const char* prefix, PressureGradient* pg);
	void resetProfile();
	void resetSlabProfile();
        void resetStressProfile();
	void resetBulkPressure();
	void resetConfinementProperties();
	long double **allocStressMatrix(size_t rows, size_t cols);
	void dellocStressMatrix(long double **matrix, size_t rows, size_t cols);
	long double **_universalStress;
	long double **_localStress;
	long double **_localStressConfinement;
	long double **_globalStressConfinement;
	
	double getGlobalBulkPressure() {return this->_universalBulkPressure; }
	void setGlobalBulkPressure(double bulkPressure) { this->_universalBulkPressure = bulkPressure; }
	double getGlobalBulkDensity() {return this->_universalBulkDensity;}
	void setGlobalBulkDensity(double bulkDensity) { this->_universalBulkDensity = bulkDensity; }
	double getGlobalBulkTemperature() {return this->_universalBulkTemp; }
	void setGlobalBulkTemperature(double bulkTemp) { this->_universalBulkTemp = bulkTemp; }
	
	long double getPressureVirial() {return this->_universalPressureVirial; }
	long double getPressureKin() {return this->_universalPressureKin; }
	double getBulkVolume() {return this->_bulkVolume; }
	long double getBulkN() {return this->_universalNBulk; }
	unsigned getAccumulatedDatasetsBulkPressure() {return this->_globalAccumulatedDatasets_BulkPressure; }
	void setStressXX(long double stress) {this->_stressXX = stress; }
	void setStressYY(long double stress) {this->_stressYY = stress; }
	long double getStressXX() {return this->_stressXX; }
	long double getStressYY() {return this->_stressYY; }
	
	// barostat
	long double getPressureVirial_barostat() {return this->_universalPressureVirial_barostat; }
	long double getPressureKin_barostat() {return this->_universalPressureKin_barostat; }
	long double getN_barostat() {return this->_universalN_barostat; }
	unsigned getAccumulatedDatasets_barostat() {return this->_globalAccumulatedDatasets_barostat; }
	void setupBarostat(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, unsigned cid);
	void recordBarostat(ParticleContainer* molCont);
	void collectBarostat(DomainDecompBase* dode);
	void resetBarostat();
	void setDifferentBarostatInterval(bool boolean) {this->_differentBarostatInterval = boolean; }
	bool getDifferentBarostatInterval() {return this->_differentBarostatInterval; }
	
	// Cylindrical Density Profile
	void confinementDensity(double radius1, double radius2, double xCentre, double yCentre);
	bool isCylindrical();
	long int unID(double qx, double qy, double qz);
	void outputCylindricalProfile(const char* prefix);	
	
	unsigned long N() {return _globalNumMolecules;}

	void Nadd(unsigned cid, int N, int localN);

	double getGlobalLength(int d) { return _globalLength[d]; }
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
	
	
	// Measurements in the confinement
	//
	unsigned long get_universalNConfinement() { return this->_universalNConfinement; }
	long double get_universalPressureVirial_Confinement() { return this->_universalPressureVirial_Confinement; }
	long double get_universalPressureKin_Confinement() { return this->_universalPressureKin_Confinement; }
	long double get_universalDOFProfile_Confinement() { return this->_universalDOFProfile_Confinement; }
	long double get_universalKineticProfile_Confinement() { return this->_universalKineticProfile_Confinement; }
	long double get_globalForce_Confinement(int d) { return this->_globalForce_Confinement[d]; }
	long double get_dMax() { return this->_dMax; }
	double get_dAveraged() { return this->_dAveraged; }
	double get_confinedVolume() { return this->_confinedVolume; }
	double get_confinedSolidSurface() { return this->_confinedSolidSurface; }
	bool isConfinement(int component) { return this->_confinementComponent[component]; }
	double get_confinementMidPoint(int d) { return this->_confinementMidPoint[d]; }
	void collectForcesOnComponentConfinement(ParticleContainer* molCont);
	unsigned getRank() { return this->_localRank; }
	double getConfinementEdge(int d) { return this->_confinementEdge[d]; }
	
	void specifyOutputFormat(std::string outputFormat) { this->_outputFormat = outputFormat; }
	std::string getOutputFormat() { return this->_outputFormat; }
	
	void specifyStressCalcMethodConfinement(std::string stressCalc) { this->_stressCalcMethodConfinement = stressCalc; }
	std::string getStressCalcMethodConfinement() { return this->_stressCalcMethodConfinement; }
	
	unsigned long getSimstep();
	unsigned long getInitStatistics();
	double getTimestepLength();
	unsigned getStressRecordTimeStep();
	unsigned getConfinementRecordTimeStep();
	double getCutoffRadius();
	unsigned getBarostatTimeInit();
	unsigned getBarostatTimeEnd();
	bool isShearRate();
private:

	//! rank of the local process
	int _localRank;

	//! Potential of the local process
	double _localUpot;
	//! Virial of the local process
	double _localVirial;
	double _localVirialI[3];
        double _localVirialIILL;
        double _localVirialIILM;
	//! global Potential
	double _globalUpot;
	//! global virial
	double _globalVirial;
	double _globalVirialI[3];
        double _globalVirialIILL;
        double _globalVirialIILM;
	//! global density
	double _globalRho;
	//! global Number of Molecules
	//! @todo redundancy?
	unsigned long _globalNumMolecules;
	unsigned long _globalOrigNumMolecules;
	//! side length of the cubic simulation box
	double _globalLength[3];

	//! does a componentwise thermostat apply?
	bool _componentwiseThermostat;
	//! thermostat IDs. negative: no thermostat, 0: global, positive: componentwise
	//! in the case of a componentwise thermostat, all components are assigned
	//! a thermostat ID different from zero.
	std::map<int, int> _componentToThermostatIdMap;
	
	//! one-dimensional thermostat (velocity scaling)
	std::map < int, int > _dimToThermostat;

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
	std::map<int, double> _universalATrans;
	std::map<int, double> _universalBRot;
	std::map<int, bool> _alphaTransCorrection;
	std::map<int, long double> _universalThT_heatFlux;
	//! should the directed movement be subtracted when calculating the temperature?
	std::map<int, bool> _universalUndirectedThermostat;
	//! stores the velocity of the directed movement
	std::map<int, double> _universalThermostatDirectedVelocity[3];
	std::map<int, double> _localThermostatDirectedVelocity[3];
	std::map<int, unsigned long>_thermostatTimeSlot[2];

	/* FIXME: This info should go into an ensemble class */
	bool _universalNVE;

	PressureGradient* _universalPG;

	//! computation of the isochoric heat capacity
	unsigned _globalUSteps;
	double _globalSigmaU;
	double _globalSigmaUU;

	//! 1 / dimension of a profile bin
	double _universalInvProfileUnit[3];
	//! number of successive profile bins in x/y/z direction
	unsigned _universalNProfileUnits[3];
	//! local N profile map
	std::map<unsigned, long double> _localNProfile;
	//! global N profile map
	std::map<unsigned, double> _universalNProfile;
	//! local directed velocity profile map
	std::map<unsigned, long double> _localvProfile[3];
	//! global directed velocity  profile map
	std::map<unsigned, double> _universalvProfile[3];
	//! local kinetic energy profile map
	std::map<unsigned, long double> _localKineticProfile;
	//! global kinetic energy profile map
	std::map<unsigned, double> _universalKineticProfile;
	//! local counter w. r. t. degrees of freedom
	std::map<unsigned, long double> _localDOFProfile;
	//! global counter w. r. t. degrees of freedom
	std::map<unsigned, double> _universalDOFProfile;
	//! how many _evaluated_ timesteps are currently accumulated in the profile?
	unsigned _globalAccumulatedDatasets;
	//! which components should be considered?
	std::map<unsigned, bool> _universalProfiledComponents;
        double _universalProfiledComponentMass;  // set from outside
        
        /*
	 * Slab Profile
        */
	//! 1 / dimension of a profile bin
	double _universalInvProfileUnitSlab[3];
	//! number of successive profile bins in x/y/z direction
	unsigned _universalNProfileUnitsSlab[3];
	//! local N profile map
	std::map<unsigned, long double> _localNProfileSlab;
	//! global N profile map
	std::map<unsigned, double> _universalNProfileSlab;
	//! local directed velocity profile map
	std::map<unsigned, long double> _localvProfileSlab[3];
	//! global directed velocity  profile map
	std::map<unsigned, double> _universalvProfileSlab[3];
	//! local kinetic energy profile map
	std::map<unsigned, long double> _localKineticProfileSlab;
	//! global kinetic energy profile map
	std::map<unsigned, double> _universalKineticProfileSlab;
	//! local counter w. r. t. degrees of freedom
	std::map<unsigned, long double> _localDOFProfileSlab;
	//! global counter w. r. t. degrees of freedom
	std::map<unsigned, double> _universalDOFProfileSlab;
	//! how many _evaluated_ timesteps are currently accumulated in the profile?
	unsigned _globalAccumulatedDatasetsSlab;
	//! which components should be considered?
	std::map<unsigned, bool> _universalProfiledComponentsSlab;
	double _maxSlabDist2;
	double _universalCenterZ;
	
	//! 1 / dimension of a profile cuboid for stresses in solids
	double _universalInvProfileUnit_Stress[3];
	//! number of successive profile cuboids in x/y/z direction for stress in solids
	unsigned _universalNProfileUnits_Stress[3];
	//! local N profile map
	std::map<unsigned, long double> _localNProfile_Stress;
	//! local N profile map of molecules that are nor part of the profiled component but part of the profile cuboid
	std::map<unsigned, long double> _localNProfileResidual_Stress;
	//! global N profile map
	std::map<unsigned, long double> _universalNProfile_Stress;
	//! global N profile map
	std::map<unsigned, long double> _universalNProfileResidual_Stress;
	//! how many _evaluated_ timesteps are currently accumulated in the profile?
	unsigned _globalAccumulatedDatasets_Stress;
	double _maxSlabDist2_Stress;
	double _universalCenterZ_Stress;
	std::map<unsigned, bool> _stressProfiledComponents;
	
	std::map<unsigned, double> _bulkCorner;
	std::map<unsigned, bool> _bulkComponent;
	long double _localPressureKin, _localPressureVirial, _universalPressureKin, _universalPressureVirial;
	long double _localNBulk, _universalNBulk;
	double _bulkVolume;
	double _universalBulkPressure, _universalBulkDensity, _universalBulkTemp, _stressXX, _stressYY;
	unsigned _globalAccumulatedDatasets_BulkPressure;
	
	//barostat
	std::map<unsigned, bool> _barostatComponent;
	double _control_top[3], _control_bottom[3];
	long double _localPressureKin_barostat, _localPressureVirial_barostat, _universalPressureKin_barostat, _universalPressureVirial_barostat;
	long double _localN_barostat, _universalN_barostat;
	unsigned _globalAccumulatedDatasets_barostat;
	bool _differentBarostatInterval;
	
	// Cylindrical Density Profile
	bool _universalCylindricalGeometry;
	double _lowerAsperityRadius, _upperAsperityRadius, _universalR2max, _universalR2min, _universalCentre[3];	
	
	// Measurements in the confinement
	std::map<int, unsigned long> _localNConfinement;
	std::map<int, unsigned long> _localWallNConfinement;
	std::map<int, long double> _localPressureVirial_Confinement;
	std::map<int, long double> _localPressureKin_Confinement;
	std::map<int, long double> _localVirialForce_Confinement;
	std::map<int, long double> _localVirialKin_Confinement;
	std::map<int, long double> _localDOFProfile_Confinement;
	std::map<int, long double> _localKineticProfile_Confinement;
	std::map<int, long double> _localvProfile_Confinement[3];
	std::map<int, unsigned long> _globalNConfinement;
	std::map<int, unsigned long> _globalWallNConfinement;
	std::map<int, long double> _globalVirialForce_Confinement;
	std::map<int, long double> _globalVirialKin_Confinement;
	std::map<int, long double> _globalPressureVirial_Confinement;
	std::map<int, long double> _globalPressureKin_Confinement;	
	std::map<int, long double> _globalDOFProfile_Confinement;
	std::map<int, long double> _globalKineticProfile_Confinement;
	std::map<int, long double> _globalvProfile_Confinement[3];
	std::map<int, long double> _localForceConfinement[3];
	std::map<int, long double> _globalForceConfinement[3];
	std::map<int, long double> _localFluidForce_Confinement[3];
	std::map<int, long double> _globalFluidForce_Confinement[3];
	std::map<int, long double> _confinementEdge;
	std::map<int, long double> _confinementMidPoint;
	std::map<int, bool> _confinementComponent;
	std::map<int, unsigned long> _confinementMidPointID;
	long double _dMax;
	double _dAveraged;
	std::map<int, unsigned> _dBinFailCount;
	unsigned _dBinFailureCount;
	double _confinedVolume;
	double _confinedSolidSurface;
	std::map<int, unsigned> _universalNProfileUnitsConfinement;
	std::map<int, double> _universalInvProfileUnitConfinement;
	std::map<int, unsigned> _universalNProfileUnitsStressConfinement;
	std::map<int, double> _universalInvProfileUnitStressConfinement;
	std::map<unsigned, double> _lowUpper;
	std::map<unsigned, double> _highLower;
	std::map<int, long double> _dBin;
	unsigned long _globalAccumulatedDatasets_ConfinementProperties;
	unsigned long _universalNConfinement;
	long double _universalPressureVirial_Confinement;
	long double _universalPressureKin_Confinement;
	long double _universalDOFProfile_Confinement;
	long double _universalKineticProfile_Confinement;
	long double _globalForce_Confinement[3];
	std::map<int, unsigned long>_localFluidicArea;
	std::map<int, unsigned long>_globalFluidicArea;
	std::map<int, std::map<unsigned long, long double> > _localDiffusiveHeatflux;
	std::map<int, std::map<unsigned long, long double> > _globalDiffusiveHeatflux;
	std::map<int, std::map<unsigned long, long double> > _localConvectiveKinHeatflux;
	std::map<int, std::map<unsigned long, long double> > _globalConvectiveKinHeatflux;
	std::map<int, std::map<unsigned long, long double> > _localConvectivePotHeatflux;
	std::map<int, std::map<unsigned long, long double> > _globalConvectivePotHeatflux;
	std::map<int, std::map<unsigned long, long double> > _localTotalHeatflux;
	std::map<int, std::map<unsigned long, long double> > _globalTotalHeatflux;
	std::map<int, std::map<unsigned long, long double> > _localDiffusiveMovement;
	std::map<int, std::map<unsigned long, long double> > _globalDiffusiveMovement;
	bool _isConfinement;
	
	
        std::map<unsigned, double> _universalTProfile; 
        std::map<unsigned, long double> _localWidomProfile;  // submit individually
        std::map<unsigned, double> _globalWidomProfile;
        std::map<unsigned, long double> _localWidomProfileTloc;  // submit individually
        std::map<unsigned, double> _globalWidomProfileTloc;
        std::map<unsigned, long double> _localWidomInstances;  // submit individually
        std::map<unsigned, double> _globalWidomInstances;
        std::map<unsigned, long double> _localWidomInstancesTloc;  // submit individually
        std::map<unsigned, double> _globalWidomInstancesTloc;
        double _universalLambda;  // set from outside
        float _globalDecisiveDensity;  // set from outside

	int _universalSelectiveThermostatCounter;
	int _universalSelectiveThermostatWarning;
	int _universalSelectiveThermostatError;

	//! local sum (over all molecules) of the mass multiplied with the squared velocity
	std::map<int, double> _local2KETrans;
	
	//! local sum (over all molecules for one direction) of the mass multiplied with the squared velocity
	std::map<int, double> _local2KETrans_1Dim;
	
	std::map<int, bool> _scale_v_1Dim;
	
	std::map<int, bool> _stressCalc;
	
	//! local sum (over all molecules) of the moment of inertia
	//! multiplied with the squared  rotational velocity
	std::map<int, double> _local2KERot; 

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

	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;
	//! modified Lorentz-Berthelot mixing rule parameters
	//! @todo more explanation
	std::vector<double> _mixcoeff;

	std::string _outputFormat;
	string _matlab;
	string _vtk;
	string _all;
	std::string _stressCalcMethodConfinement;
};


#endif /*DOMAIN_H_*/
