/*
 * DistControl.h
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#ifndef DISTCONTROL_H_
#define DISTCONTROL_H_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdint>
#include "plugins/PluginBase.h"
#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "utils/CommVar.h"
#include "molecules/MoleculeForwardDeclaration.h"


class XMLfileUnits;
class Domain;
class ParticleContainer;
class DomainDecompBase;

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

enum DistControlUpdateMethods
{
	DCUM_UNKNOWN = 0,
	DCUM_DENSITY_PROFILE = 1,
	DCUM_DENSITY_PROFILE_DERIVATION = 2,
//	DCUM_FORCE_PROFILE = 3,
};

enum DistControlInitMethods
{
	DCIM_UNKNOWN = 0,
	DCIM_START_CONFIGURATION = 1,
	DCIM_MIDPOINT_VALUES = 2,
	DCIM_READ_FROM_FILE = 3,
};

class DistControl : public PluginBase, public ControlInstance, public SubjectBase
{
public:
	DistControl();
	virtual ~DistControl();

	std::string getShortName() override {return "DiC";}

	/** @brief Read in XML configuration for DistControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="DistControl">
			<control>
				<update>5000</update>   <!-- update frequency -->
			</control>
			<filenames>
				<control>DistControl.dat</control>         <!-- log file of updated positions -->
				<profiles>DistControlProfiles</profiles>   <!-- file prefix for density profiles to determine interface position(s) -->
			</filenames>
			<subdivision type="width">   <!-- type="number|width" for subdivision of domain into bins for density profile sampling -->
				<width>FLOAT</width>     <!-- bin width -->
				<number>1</number>       <!-- number of bins -->
			</subdivision>
			<init type="startconfig">    <!-- type="startconfig|values|file" for init positions-->
				<values> <left>FLOAT</left> <right>FLOAT</right> </values>   <!-- specify init values for left and right interface -->
				<file>../path/to/file/DistControl.dat</file>                 <!-- read values from specified file -->
				<simstep>INT</simstep>                                       <!-- read from file in line simstep=INT -->
			</init>
			<method type="denderiv">   <!-- type="density"|denderiv" method to determine interface positions-->
				<componentID>INT</componentID>                            <!-- target component, 0:all components -->
				<neighbourvals algorithm="smooth">INT</neighbourvals>     <!-- neighbour values used to smooth the profile -->
				<neighbourvals algorithm="derivate">INT</neighbourvals>   <!-- neighbour values used to calculate derivation profile by linear regression -->
				<density>FLOAT</density>                                  <!-- vapor density to identify vapor phase -->
			</method>
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain) override;

	void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep, bool signalled) override {};

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep) override;

	void siteWiseForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep) override {};
	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep) override {}

	void endStep(ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep) override {};

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("DistControl");}
	static PluginBase* createInstance() {return new DistControl();}

	// set subdivision
	void SetSubdivision(const uint32_t& numBins) {_binParams.count = numBins; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivision(const double& dBinWidth) {_binParams.width = dBinWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void PrepareSubdivision();  // need to be called before PrepareDataStructures()

	// data structures
	void PrepareDataStructures();

	// init
	void InitPositions(const double& dInterfaceMidLeft, const double& dInterfaceMidRight);

	double GetInterfaceMidLeft() {return _dInterfaceMidLeft;}
	double GetInterfaceMidRight() {return _dInterfaceMidRight;}
	unsigned int GetUpdateFreq() {return _controlFreqs.update;}
	unsigned int GetWriteFreqProfiles() {return _controlFreqs.write.profiles;}

	void Init(ParticleContainer* particleContainer);
	void WriteHeader();
	void WriteData(const uint64_t& simstep);
	void WriteDataProfiles(const uint64_t& simstep);

	// place method inside loop over molecule container
	void SampleProfiles(Molecule* mol);

	void UpdatePositionsInit(ParticleContainer* particleContainer);  // call in Simulation::prepare_start()
	void UpdatePositions(const uint64_t& simstep);

	// SubjectBase methods
	void registerObserver(ObserverBase* observer) override;
	void unregisterObserver(ObserverBase* observer) override;
	void informObserver() override;

private:
	// place methods after the loop
	void CalcProfiles();
	void EstimateInterfaceMidpoint();  // called by UpdatePositions
	void EstimateInterfaceMidpointsByForce();
	void ResetLocalValues();

	// data structures
	void InitDataStructures();

	// processing profiles
	void SmoothProfile(double* dData, double* dSmoothData, const uint64_t& nNumVals, const uint32_t& nNeighbourVals);
	void SmoothProfiles(const uint32_t& nNeighbourVals);
	void DerivateProfile(double* dDataX, double* dDataY, double* dDerivDataY, const uint64_t& nNumVals, const uint32_t& nNeighbourVals);
	void DerivateProfiles(const uint32_t& nNeighbourVals);


private:
	double _dInterfaceMidLeft;
	double _dInterfaceMidRight;

	uint16_t _nNumComponents;
	uint16_t _nTargetCompID;
	uint64_t _nNumValuesScalar;
	std::vector<uint64_t> _nOffsets;

	CommVar<std::vector<uint64_t> > _nNumMolecules;
	CommVar<std::vector<uint64_t> > _dForceSum;
	std::vector<double> _dMidpointPositions;
	std::vector<double> _dDensityProfile;
	std::vector<double> _dDensityProfileSmoothed;
	std::vector<double> _dDensityProfileSmoothedDerivation;
	std::vector<double> _dForceProfile;
	std::vector<double> _dForceProfileSmoothed;

	// update method
	int _nMethod;
	double _dVaporDensity;
	uint16_t _nNeighbourValsSmooth;
	uint16_t _nNeighbourValsDerivate;

	int _nMethodInit;
	std::string _strFilenameInit;
	uint64_t _nRestartTimestep;

	// write data
	std::string _strFilename;
	std::string _strFilenameProfilesPrefix;

	// observer
	std::vector<ObserverBase*> _observer;

	int _nSubdivisionOpt;

	struct BinParamsType
	{
		uint32_t count;
		double width;
		double invWidth;
		double volume;
	} _binParams;

	struct ControlFreqType
	{
		uint32_t update;
		struct WriteFreqType
		{
			uint32_t data;
			uint32_t profiles;
		} write;
	} _controlFreqs;

};  // class DistControl


#endif /* DISTCONTROL_H_ */
