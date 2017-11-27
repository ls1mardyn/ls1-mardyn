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
	void SetParamProfiles( unsigned long initSamplingProfiles, unsigned long writeFrequencyProfiles, unsigned long stopSamplingProfiles)
	{
		_initSamplingProfiles   = initSamplingProfiles;
		_writeFrequencyProfiles = writeFrequencyProfiles;
		_stopSamplingProfiles   = stopSamplingProfiles;
		_SamplingEnabledProfiles = true;
		_dInvertNumSamplesProfiles = 1. / ( (double)(_writeFrequencyProfiles) );  // needed for e.g. density profiles
	}
	void SetParamProfiles()
	{
		_SamplingEnabledProfiles = true;
		_dInvertNumSamplesProfiles = 1. / ( (double)(_writeFrequencyProfiles) );  // needed for e.g. density profiles
	}

	// subdivision
	void SetSubdivisionProfiles(const unsigned int& nNumSlabs) {_nNumBinsProfiles = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionProfiles(const double& dSlabWidth) {_dBinWidthProfilesInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void SetSubdivisionVDF(const unsigned int& nNumSlabs) {_numBinsVDF = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionVDF(const double& dSlabWidth) {_dBinWidthVDFInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void PrepareSubdivisionProfiles();  // need to be called before data structure allocation
	void PrepareSubdivisionVDF();  		// need to be called before data structure allocation
	void PrepareSubdivisionFieldYR();  	// need to be called before data structure allocation


	// init sampling
	void InitSamplingProfiles(int nDimension);
	void InitSamplingVDF(int nDimension);
	void InitSamplingFieldYR(int nDimension);
	void DoDiscretisationProfiles(int nDimension);
	void DoDiscretisationVDF(int nDimension);
	void DoDiscretisationFieldYR(int nDimension);

	// molecule container loop methods
	void SampleProfiles(Molecule* molecule, int nDimension);
	void SampleVDF(Molecule* molecule, int nDimension);
	void SampleFieldYR(Molecule* molecule);

	// calc global values
	void CalcGlobalValuesProfiles(DomainDecompBase* domainDecomp, Domain* domain);
	void CalcGlobalValuesVDF();
	void CalcGlobalValuesFieldYR(DomainDecompBase* domainDecomp, Domain* domain);

	// output
	void WriteDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);
	void WriteDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep);
	void WriteDataFieldYR(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);

	void UpdateSlabParameters();

private:
	// reset local values
	void ResetLocalValuesProfiles();
	void ResetOutputDataProfiles();
	void ResetLocalValuesVDF();
	void ResetLocalValuesFieldYR();

	void InitComponentSpecificParamsVDF();
	void ShowComponentSpecificParamsVDF();
private:
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

	// parameters
	unsigned long _initSamplingProfiles;
	unsigned long _writeFrequencyProfiles;
	unsigned long _stopSamplingProfiles;
	unsigned int _nNumBinsProfiles;

	double  _dBinWidthProfiles;
	double  _dBinWidthProfilesInit;
	double  _dBinVolumeProfiles;
	double* _dBinMidpointsProfiles;

	// offsets
	// TODO: Use only 1 offset array: for scalar one can take the offsets of dimension x, right???
	unsigned long**  _nOffsetScalar;  //                  [direction all|+|-][component]
	unsigned long*** _nOffsetVector;  // [dimension x|y|z][direction all|+|-][component]

	unsigned long _nNumValsScalar;
	unsigned long _nNumValsVector;

	std::vector<double> _vecMass;
	double _dInvertNumSamplesProfiles;
	double _dInvertBinVolumeProfiles;
	double _dInvertBinVolSamplesProfiles;

	// Scalar quantities
	// [direction all|+|-][component][position]
	unsigned long* _nNumMoleculesLocal;
	unsigned long* _nNumMoleculesGlobal;
	unsigned long* _nRotDOFLocal;
	unsigned long* _nRotDOFGlobal;
	double* _d2EkinRotLocal;
	double* _d2EkinRotGlobal;

	// output profiles
	double* _dDensity;
	double* _d2EkinTotal;
	double* _d2EkinTrans;
	double* _d2EkinDrift;
	double* _d2EkinRot;
	double* _d2EkinT;
	double* _dTemperature;
	double* _dTemperatureTrans;
	double* _dTemperatureRot;

	// Vector quantities
	// [dimension x|y|z][direction all|+|-][component][position]
	double* _dVelocityLocal;
	double* _dVelocityGlobal;
	double* _dSquaredVelocityLocal;
	double* _dSquaredVelocityGlobal;
	double* _dForceLocal;
	double* _dForceGlobal;

	// output profiles
	double* _dForce;
	double* _dDriftVelocity;
	double* _d2EkinTransComp;
	double* _d2EkinDriftComp;
	double* _dTemperatureComp;

	// --- VDF ---

	// parameters
	unsigned long _initSamplingVDF;
	unsigned long _writeFrequencyVDF;
	unsigned long _stopSamplingVDF;

	// control
	bool _bDiscretisationDoneVDF;
	bool _SamplingEnabledVDF;

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
	uint64_t* _VDF_pjy_abs_local;
	uint64_t* _VDF_pjy_pvx_local;
	uint64_t* _VDF_pjy_pvy_local;
	uint64_t* _VDF_pjy_pvz_local;
	uint64_t* _VDF_pjy_nvx_local;
	uint64_t* _VDF_pjy_nvz_local;

	uint64_t* _VDF_njy_abs_local;
	uint64_t* _VDF_njy_pvx_local;
	uint64_t* _VDF_njy_pvz_local;
	uint64_t* _VDF_njy_nvx_local;
	uint64_t* _VDF_njy_nvy_local;
	uint64_t* _VDF_njy_nvz_local;

	// global
	uint64_t* _VDF_pjy_abs_global;
	uint64_t* _VDF_pjy_pvx_global;
	uint64_t* _VDF_pjy_pvy_global;
	uint64_t* _VDF_pjy_pvz_global;
	uint64_t* _VDF_pjy_nvx_global;
	uint64_t* _VDF_pjy_nvz_global;

	uint64_t* _VDF_njy_abs_global;
	uint64_t* _VDF_njy_pvx_global;
	uint64_t* _VDF_njy_pvz_global;
	uint64_t* _VDF_njy_nvx_global;
	uint64_t* _VDF_njy_nvy_global;
	uint64_t* _VDF_njy_nvz_global;

	std::array<std::array<uint64_t*,4>,3> _dataPtrs;

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
	double* _dBinMidpointsFieldYR;
	double* _dShellMidpointsFieldYR;
	double* _dShellVolumesFieldYR;
	double* _dInvShellVolumesFieldYR;
	// equal shell volumes
	double _dShellVolumeFieldYR;
	double _dInvShellVolumeFieldYR;

	// offsets
	uint64_t**** _nOffsetFieldYR;  // [dimension x|y|z][component][section][positionR], dimension = 0 for scalar values, section: 0: all, 1: upper, 2: lower section
	uint64_t _nNumValsFieldYR;

	// Scalar quantities
	// [component][section][positionR][positionY]
	uint64_t* _nNumMoleculesFieldYRLocal;
	uint64_t* _nNumMoleculesFieldYRGlobal;

	// output profiles
	double* _dDensityFieldYR;

	// output file
	std::string _strFilePrefixFieldYR;
	uint8_t _nFileTypeFieldYR;
};

class XMLfileUnits;
class RegionSampling : public ControlInstance
{
public:
	RegionSampling(Domain* domain, DomainDecompBase* domainDecomp);
	virtual ~RegionSampling();

	std::string GetShortName() {return "ReS";}
	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(SampleRegion* region);
	int GetNumRegions() {return _vecSampleRegions.size();}
	SampleRegion* GetSampleRegion(unsigned short nRegionID) {return _vecSampleRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

	void PrepareRegionSubdivisions();  // need to be called before method Init() that allocates the data structures
	void Init();
	void DoSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);

	void WriteData(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);

private:
	std::vector<SampleRegion*> _vecSampleRegions;

	unsigned long _initSamplingProfiles;
	unsigned long _writeFrequencyProfiles;
	unsigned long _initSamplingVDF;
	unsigned long _writeFrequencyVDF;
};

#endif /* REGIONSAMPLING_H_ */

