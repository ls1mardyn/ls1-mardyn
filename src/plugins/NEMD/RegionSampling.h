/*
 * RegionSampling.h
 *
 *  Created on: 18.03.2015
 *      Author: mheinen
 */

#ifndef REGIONSAMPLING_H_
#define REGIONSAMPLING_H_

#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "plugins/PluginBase.h"

#include <vector>
#include <array>
#include <cstdint>

enum RegionSamplingDimensions
{
	RS_DIMENSION_X = 0,
	RS_DIMENSION_Y = 1,
	RS_DIMENSION_Z = 2,
};

enum RegionSamplingFileTypes : uint8_t
{
	RSFT_ASCII = 1,
	RSFT_BINARY = 2
};

struct ComponentSpecificParamsVDF
{
	bool bSamplingEnabled;
	double dVeloMax;
	double dVelocityClassWidth;
	double dInvVelocityClassWidth;
	uint32_t numVelocityClasses;
	uint32_t numVals;
	uint32_t nOffsetDataStructure;
	std::vector<double> dDiscreteVelocityValues;
};

class DistControl;
class XMLfileUnits;
class Domain;
class DomainDecompBase;
class RegionSampling;

class SampleRegion : public CuboidRegionObs
{
public:
	SampleRegion(RegionSampling* parent, double dLowerCorner[3], double dUpperCorner[3] );
	virtual ~SampleRegion();

	void readXML(XMLfileUnits& xmlconfig);

	// set parameters
	void setParamProfiles( unsigned long initSamplingProfiles, unsigned long writeFrequencyProfiles, unsigned long stopSamplingProfiles)
	{
		_initSamplingProfiles   = initSamplingProfiles;
		_writeFrequencyProfiles = writeFrequencyProfiles;
		_stopSamplingProfiles   = stopSamplingProfiles;
		_SamplingEnabledProfiles = true;
		_dInvertNumSamplesProfiles = 1. / ( (double)(_writeFrequencyProfiles) );  // needed for e.g. density profiles
	}
	void setParamProfiles()
	{
		_SamplingEnabledProfiles = true;
		_dInvertNumSamplesProfiles = 1. / ( (double)(_writeFrequencyProfiles) );  // needed for e.g. density profiles
	}

	// subdivision
	void setSubdivisionProfiles(const unsigned int& nNumSlabs) {_nNumBinsProfiles = nNumSlabs; _nSubdivisionOptProfiles = SDOPT_BY_NUM_SLABS;}
	void setSubdivisionProfiles(const double& dSlabWidth) {_dBinWidthProfilesInit = dSlabWidth; _nSubdivisionOptProfiles = SDOPT_BY_SLAB_WIDTH;}
	void setSubdivisionVDF(const unsigned int& nNumSlabs) {_numBinsVDF = nNumSlabs; _nSubdivisionOptVDF = SDOPT_BY_NUM_SLABS;}
	void setSubdivisionVDF(const double& dSlabWidth) {_dBinWidthVDFInit = dSlabWidth; _nSubdivisionOptVDF = SDOPT_BY_SLAB_WIDTH;}
	void prepareSubdivisionProfiles();  // need to be called before data structure allocation
	void prepareSubdivisionVDF();  		// need to be called before data structure allocation
	void prepareSubdivisionFieldYR();  	// need to be called before data structure allocation


	// init sampling
	void initSamplingProfiles(int nDimension);
	void initSamplingVDF(int nDimension);
	void initSamplingFieldYR(int nDimension);
	void doDiscretisationProfiles(int nDimension);
	void doDiscretisationVDF(int nDimension);
	void doDiscretisationFieldYR(int nDimension);

	// molecule container loop methods
	void sampleProfiles(Molecule* molecule, int nDimension, unsigned long simstep);
	void sampleVDF(Molecule* molecule, int nDimension, unsigned long simstep);
	void sampleFieldYR(Molecule* molecule, unsigned long simstep);

	// calc global values
	void calcGlobalValuesProfiles(DomainDecompBase* domainDecomp, Domain* domain);
	void calcGlobalValuesVDF();
	void calcGlobalValuesFieldYR(DomainDecompBase* domainDecomp, Domain* domain);

	// output
	void writeDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);
	void writeDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep);
	void writeDataFieldYR(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);

	void updateSlabParameters();

private:
	// reset local values
	void resetLocalValuesProfiles();
	void resetOutputDataProfiles();
	void resetLocalValuesVDF();
	void resetLocalValuesFieldYR();

	void initComponentSpecificParamsVDF();
	void showComponentSpecificParamsVDF();

	static void get_v(double* q, Molecule* mol);
	static void get_F(double* q, Molecule* mol);
	static void get_v2(double& q, Molecule* mol);
	static void get_F2(double& q, Molecule* mol);
	void(*_fptr)(double*, Molecule*);
	void(*_f2ptr)(double&, Molecule*);

	// observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	DistControl* getDistControl();

private:
	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

	// instances / ID
	static unsigned short _nStaticID;

	uint32_t _numComponents;

	// ******************
	// sampling variables
	// ******************

	// --- temperature / density profiles, SamplingZone ---

	// control
	bool _bDiscretisationDoneProfiles;
	bool _SamplingEnabledProfiles;
	int _nSubdivisionOptProfiles;

	// parameters
	unsigned long _initSamplingProfiles;
	unsigned long _writeFrequencyProfiles;
	unsigned long _stopSamplingProfiles;
	unsigned int _nNumBinsProfiles;

	double  _dBinWidthProfiles;
	double  _dBinWidthProfilesInit;
	double  _dBinVolumeProfiles;
	std::vector<double> _dBinMidpointsProfiles;

	// offsets
	// TODO: Use only 1 offset array: for scalar one can take the offsets of dimension x, right???
	std::array<std::vector<unsigned long>, 3> _nOffsetScalar; //                                  [direction all|+|-][component]
	std::array<std::array<std::vector<unsigned long>, 3>, 3> _nOffsetVector;  // [dimension x|y|z][direction all|+|-][component]

	unsigned long _nNumValsScalar;
	unsigned long _nNumValsVector;

	std::vector<double> _vecMass;
	double _dInvertNumSamplesProfiles;
	double _dInvertBinVolumeProfiles;
	double _dInvertBinVolSamplesProfiles;

	// Scalar quantities
	// [direction all|+|-][component][position]
	std::vector<unsigned long> _nNumMoleculesLocal;
	std::vector<unsigned long> _nNumMoleculesGlobal;
	std::vector<unsigned long> _nRotDOFLocal;
	std::vector<unsigned long> _nRotDOFGlobal;
	std::vector<double> _d2EkinRotLocal;
	std::vector<double> _d2EkinRotGlobal;
	std::vector<double> _dEpotLocal;
	std::vector<double> _dEpotGlobal;

	// output profiles
	std::vector<double> _dDensity;
	std::vector<double> _d2EkinTotal;
	std::vector<double> _d2EkinTrans;
	std::vector<double> _d2EkinDrift;
	std::vector<double> _d2EkinRot;
	std::vector<double> _d2EkinT;
	std::vector<double> _dEpot;
	std::vector<double> _dTemperature;
	std::vector<double> _dTemperatureTrans;
	std::vector<double> _dTemperatureRot;

	// Vector quantities
	// [dimension x|y|z][direction all|+|-][component][position]
	std::vector<double> _dVelocityLocal;
	std::vector<double> _dVelocityGlobal;
	std::vector<double> _dSquaredVelocityLocal;
	std::vector<double> _dSquaredVelocityGlobal;
	std::vector<double> _dForceLocal;
	std::vector<double> _dForceGlobal;
	std::vector<double> _dHeatfluxLocal1;
	std::vector<double> _dHeatfluxLocal2;
	std::vector<double> _dHeatfluxLocal3;
	std::vector<double> _dHeatfluxGlobal1;
	std::vector<double> _dHeatfluxGlobal2;
	std::vector<double> _dHeatfluxGlobal3;
	std::vector<double> _dVirialLocal;
	std::vector<double> _dVirialGlobal;

	// output profiles
	std::vector<double> _dForce;
	std::vector<double> _dHeatflux1;
	std::vector<double> _dHeatflux2;
	std::vector<double> _dHeatflux3;
	std::vector<double> _dDriftVelocity;
	std::vector<double> _d2EkinTransComp;
	std::vector<double> _d2EkinDriftComp;
	std::vector<double> _dTemperatureComp;
	std::vector<double> _dVirial;

	// --- VDF ---

	// parameters
	unsigned long _initSamplingVDF;
	unsigned long _writeFrequencyVDF;
	unsigned long _stopSamplingVDF;

	// control
	bool _bDiscretisationDoneVDF;
	bool _SamplingEnabledVDF;
	int _nSubdivisionOptVDF;

	double _dBinWidthVDF;
	double _dInvBinWidthVDF;
	double _dBinWidthVDFInit;
	std::vector<double> _dBinMidpointsVDF;

	// component specific parameters
	std::vector<ComponentSpecificParamsVDF> _vecComponentSpecificParamsVDF;

	// data structures
	uint32_t _numBinsVDF;
	uint32_t _numValsVDF;

	// local
	std::vector<uint64_t> _VDF_pjy_abs_local;
	std::vector<uint64_t> _VDF_pjy_pvx_local;
	std::vector<uint64_t> _VDF_pjy_pvy_local;
	std::vector<uint64_t> _VDF_pjy_pvz_local;
	std::vector<uint64_t> _VDF_pjy_nvx_local;
	std::vector<uint64_t> _VDF_pjy_nvz_local;

	std::vector<uint64_t> _VDF_njy_abs_local;
	std::vector<uint64_t> _VDF_njy_pvx_local;
	std::vector<uint64_t> _VDF_njy_pvz_local;
	std::vector<uint64_t> _VDF_njy_nvx_local;
	std::vector<uint64_t> _VDF_njy_nvy_local;
	std::vector<uint64_t> _VDF_njy_nvz_local;

	// global
	std::vector<uint64_t> _VDF_pjy_abs_global;
	std::vector<uint64_t> _VDF_pjy_pvx_global;
	std::vector<uint64_t> _VDF_pjy_pvy_global;
	std::vector<uint64_t> _VDF_pjy_pvz_global;
	std::vector<uint64_t> _VDF_pjy_nvx_global;
	std::vector<uint64_t> _VDF_pjy_nvz_global;

	std::vector<uint64_t> _VDF_njy_abs_global;
	std::vector<uint64_t> _VDF_njy_pvx_global;
	std::vector<uint64_t> _VDF_njy_pvz_global;
	std::vector<uint64_t> _VDF_njy_nvx_global;
	std::vector<uint64_t> _VDF_njy_nvy_global;
	std::vector<uint64_t> _VDF_njy_nvz_global;

	std::array<std::array<uint64_t*,4>,3> _dataPtrs;
	std::string _fnamePrefixVDF;

	// --- fieldYR ---

	// control
	bool _bDiscretisationDoneFieldYR;
	bool _SamplingEnabledFieldYR;
	int _nSubdivisionOptFieldYR_Y;
	int _nSubdivisionOptFieldYR_R;

	// parameters
	uint64_t _initSamplingFieldYR;
	uint64_t _writeFrequencyFieldYR;
	uint64_t _stopSamplingFieldYR;
	uint32_t _nNumBinsFieldYR;
	uint32_t _nNumShellsFieldYR;
	double  _dBinWidthInitFieldYR;
	double  _dShellWidthInitFieldYR;

	double  _dBinWidthFieldYR;
	double  _dShellWidthFieldYR;
	double  _dShellWidthSquaredFieldYR;
	std::vector<double> _dBinMidpointsFieldYR;
	std::vector<double> _dShellMidpointsFieldYR;
	std::vector<double> _dShellVolumesFieldYR;
	std::vector<double> _dInvShellVolumesFieldYR;
	// equal shell volumes
	double _dShellVolumeFieldYR;
	double _dInvShellVolumeFieldYR;

	// offsets
	std::array<std::array<std::vector<std::vector<uint64_t>>, 3>, 3> _nOffsetFieldYR; // [dimension x|y|z][component][section][positionR], dimension = 0 for scalar values, section: 0: all, 1: upper, 2: lower section
	uint64_t _nNumValsFieldYR;

	// Scalar quantities
	// [component][section][positionR][positionY]
	std::vector<uint64_t> _nNumMoleculesFieldYRLocal;
	std::vector<uint64_t> _nNumMoleculesFieldYRGlobal;

	// output profiles
	std::vector<double> _dDensityFieldYR;

	// output file
	std::string _strFilePrefixFieldYR;
	uint8_t _nFileTypeFieldYR;
	
	bool _boolSingleComp;
};

class XMLfileUnits;
class RegionSampling : public ControlInstance, public PluginBase
{
public:
	RegionSampling();
	virtual ~RegionSampling();

	std::string getShortName() override {return "ReS";}

	/** @brief Read in XML configuration for DistControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
	<plugin name="RegionSampling">
		<region>
			<coords>   <!-- lc and uc: lower and upper corner of cuboid sampling region -->
				<lcx>FLOAT</lcx> <lcy refcoordsID="0">FLOAT</lcy> <lcz>FLOAT</lcz>
				<ucx>FLOAT</ucx> <ucy refcoordsID="0">FLOAT</ucy> <ucz>FLOAT</ucz>
			</coords>

			<sampling type="profiles">   <!-- Sampling profiles of various scalar and vector quantities, e.g. temperature, density, force, hydrodynamic velocity -->
				<control>
					<start>INT</start>           <!-- start time step -->
					<frequency>INT</frequency>   <!-- frequency of writing profiles -->
					<stop>INT</stop>             <!-- stop time step -->
				</control>
				<subdivision type="width">       <!-- type="number | width" => subdivision of region into bins -->
					<width>FLOAT</width>         <!-- bin width -->
					<number>INT</number>         <!-- number of bins -->
				</subdivision>
			</sampling>

			<sampling type="VDF" single_component="1">                <!-- Sampling of velocity distribution functions (VDF); Do not differ between components, if single_component is 1-->
				<control>
					<start>INT</start>           <!-- start time step -->
					<frequency>INT</frequency>   <!-- frequency of writing profiles -->
					<stop>INT</stop>             <!-- stop time step -->
				</control>
				<subdivision type="width">       <!-- type="number | width" => subdivision of region into bins -->
					<width>FLOAT</width>         <!-- bin width -->
					<number>INT</number>         <!-- number of bins -->
				</subdivision>
				<discretizations>
					<discretization cid="INT">          <!-- discretization of the velocity into discrete classes, for component cid="INT" -->
						<numclasses>INT</numclasses>    <!-- number of velocity classes -->
						<maxvalue>FLOAT</maxvalue>      <!-- maximum velocity that will be sampled -->
					</discretization>
				</discretizations>
			</sampling>

			<sampling type="fieldYR">             <!-- Sampling of a density field \rho(y,r) in cylinder shells, depending on a Cartesian coordinate y and radius r -->
				<outputfile type="binary">        <!-- type="ASCII | binary" of output files -->
					<prefix>STRING</prefix>       <!-- file prefix of output files -->
				</outputfile>
				<control>
					<start>INT</start>           <!-- start time step -->
					<frequency>INT</frequency>   <!-- frequency of writing profiles -->
					<stop>INT</stop>             <!-- stop time step -->
				</control>
				<subdivision dim="y" type="width">   <!-- dim="y": subdivision of region into bins (y axis), type="number | width" -->
					<width>FLOAT</width>             <!-- bin width -->
					<number>INT</number>             <!-- number of bins -->
				</subdivision>
				<subdivision dim="r" type="width">   <!-- dim="r": subdivision of region into cylinder shells around y axis, type="number | width" -->
					<width>FLOAT</width>             <!-- width of the inner most shell, width of outer shells will be calculated such that the shell volume is const. -->
					<number>INT</number>             <!-- number of shells -->
				</subdivision>
			</sampling>
		</region>

		<region>
			...
		</region>
	</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override {}

	void endStep(
			ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep) override;

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("RegionSampling");}
	static PluginBase* createInstance() {return new RegionSampling();}

protected:
	void addRegion(SampleRegion* region);
	int getNumRegions() {return _vecSampleRegions.size();}
	SampleRegion* getSampleRegion(unsigned short nRegionID) {return _vecSampleRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

	// sample profiles and vdf
	void doSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
	// write out profiles and vdf
	void writeData(DomainDecompBase* domainDecomp, unsigned long simstep);
	void prepareRegionSubdivisions();  // need to be called before allocating the data structures

private:
	std::vector<SampleRegion*> _vecSampleRegions;
	
	unsigned long _initSamplingProfiles;
	unsigned long _writeFrequencyProfiles;
	unsigned long _initSamplingVDF;
	unsigned long _writeFrequencyVDF;
};

#endif /* REGIONSAMPLING_H_ */

