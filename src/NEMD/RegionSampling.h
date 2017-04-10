/*
 * RegionSampling.h
 *
 *  Created on: 18.03.2015
 *      Author: mheinen
 */

#ifndef REGIONSAMPLING_H_
#define REGIONSAMPLING_H_

enum RegionSamplingDimensions
{
	RS_DIMENSION_X = 0,
	RS_DIMENSION_Y = 1,
	RS_DIMENSION_Z = 2,
};

#include <vector>
#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "molecules/MoleculeForwardDeclaration.h"

class Domain;
class DomainDecompBase;
class RegionSampling;

class SampleRegion : public CuboidRegionObs
{
public:
	SampleRegion(ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3] );
	~SampleRegion();

	// set parameters
	void SetParamProfiles( unsigned long initSamplingProfiles, unsigned long writeFrequencyProfiles)
	{
		_initSamplingProfiles   = initSamplingProfiles;
		_writeFrequencyProfiles = writeFrequencyProfiles;
		_SamplingEnabledProfiles = true;
		_dInvertNumSamplesProfiles = 1. / ( (double)(_writeFrequencyProfiles) );  // needed for e.g. density profiles
	}
	void SetParamVDF( unsigned long initSamplingVDF, unsigned long writeFrequencyVDF,
					  unsigned int nNumDiscreteStepsVDF, double dVeloMax)
	{
		_initSamplingVDF       = initSamplingVDF;
		_writeFrequencyVDF     = writeFrequencyVDF;
		_nNumDiscreteStepsVDF  = nNumDiscreteStepsVDF;
		 _dVeloMax             = dVeloMax;
		_SamplingEnabledVDF = true;
	}

	// subdivision
	void SetSubdivisionProfiles(const unsigned int& nNumSlabs) {_nNumBinsProfiles = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionProfiles(const double& dSlabWidth) {_dBinWidthProfilesInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void SetSubdivisionVDF(const unsigned int& nNumSlabs) {_nNumBinsVDF = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionVDF(const double& dSlabWidth) {_dBinWidthVDFInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void PrepareSubdivisionProfiles();  // need to be called before data structure allocation
	void PrepareSubdivisionVDF();  		// need to be called before data structure allocation


	// init sampling
	void InitSamplingProfiles(int nDimension);
	void InitSamplingVDF(int nDimension);
	void DoDiscretisationProfiles(int nDimension);
	void DoDiscretisationVDF(int nDimension);

	// molecule container loop methods
	void SampleProfiles(Molecule* molecule, int nDimension);
	void SampleVDF(Molecule* molecule, int nDimension);

	// calc global values
	void CalcGlobalValuesProfiles(DomainDecompBase* domainDecomp, Domain* domain);
	void CalcGlobalValuesVDF();

	// output
	void WriteDataProfilesOld(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);
	void WriteDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);
	void WriteDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep);

	void UpdateSlabParameters();

private:
	// reset local values
	void ResetLocalValuesProfiles();
	void ResetOutputDataProfiles();
	void ResetLocalValuesVDF();

private:

	// instances / ID
	static unsigned short _nStaticID;

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

	unsigned int _nNumComponents;
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
	unsigned int _nNumBinsVDF;
	unsigned int _nNumDiscreteStepsVDF;
	double _dVeloMax;

	// control
	bool _bDiscretisationDoneVDF;
	bool _SamplingEnabledVDF;

	double  _dBinWidthVDF;
	double  _dBinWidthVDFInit;
	double* _dBinMidpointsVDF;
	double* _dDiscreteVelocityValues;

	// local
	unsigned long** _veloDistrMatrixLocal_py_abs;
	unsigned long** _veloDistrMatrixLocal_py_pvx;
	unsigned long** _veloDistrMatrixLocal_py_pvy;
	unsigned long** _veloDistrMatrixLocal_py_pvz;
	unsigned long** _veloDistrMatrixLocal_py_nvx;
	unsigned long** _veloDistrMatrixLocal_py_nvy;
	unsigned long** _veloDistrMatrixLocal_py_nvz;

	unsigned long** _veloDistrMatrixLocal_ny_abs;
	unsigned long** _veloDistrMatrixLocal_ny_pvx;
	unsigned long** _veloDistrMatrixLocal_ny_pvy;
	unsigned long** _veloDistrMatrixLocal_ny_pvz;
	unsigned long** _veloDistrMatrixLocal_ny_nvx;
	unsigned long** _veloDistrMatrixLocal_ny_nvy;
	unsigned long** _veloDistrMatrixLocal_ny_nvz;

	// global
	unsigned long** _veloDistrMatrixGlobal_py_abs;
	unsigned long** _veloDistrMatrixGlobal_py_pvx;
	unsigned long** _veloDistrMatrixGlobal_py_pvy;
	unsigned long** _veloDistrMatrixGlobal_py_pvz;
	unsigned long** _veloDistrMatrixGlobal_py_nvx;
	unsigned long** _veloDistrMatrixGlobal_py_nvy;
	unsigned long** _veloDistrMatrixGlobal_py_nvz;

	unsigned long** _veloDistrMatrixGlobal_ny_abs;
	unsigned long** _veloDistrMatrixGlobal_ny_pvx;
	unsigned long** _veloDistrMatrixGlobal_ny_pvy;
	unsigned long** _veloDistrMatrixGlobal_ny_pvz;
	unsigned long** _veloDistrMatrixGlobal_ny_nvx;
	unsigned long** _veloDistrMatrixGlobal_ny_nvy;
	unsigned long** _veloDistrMatrixGlobal_ny_nvz;
};


class RegionSampling : public ControlInstance
{
public:
	RegionSampling(Domain* domain, DomainDecompBase* domainDecomp);
	~RegionSampling();

	std::string GetShortName() {return "ReS";}
	void AddRegion(SampleRegion* region);
	int GetNumRegions() {return _vecSampleRegions.size();}
	SampleRegion* GetSampleRegion(unsigned short nRegionID) {return _vecSampleRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

	void PrepareRegionSubdivisions();  // need to be called before method Init() that allocates the data structures
	void Init();
	void DoSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);

	void WriteData(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);

	unsigned short GetNumComponents() {return _nNumComponents;}

private:
	std::vector<SampleRegion*> _vecSampleRegions;

	unsigned long _initSamplingProfiles;
	unsigned long _writeFrequencyProfiles;
	unsigned long _initSamplingVDF;
	unsigned long _writeFrequencyVDF;

	unsigned short _nNumComponents;
};

#endif /* REGIONSAMPLING_H_ */

