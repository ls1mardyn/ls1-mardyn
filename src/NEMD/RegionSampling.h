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

class Domain;
class DomainDecompBase;
class Molecule;
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
    }
    void SetParamVDF( unsigned long initSamplingVDF, unsigned long writeFrequencyVDF,
                      unsigned int nNumDiscreteStepsVDF, double dVeloMax)
    {
        _initSamplingVDF       = initSamplingVDF;
        _writeFrequencyVDF     = writeFrequencyVDF;
        _nNumDiscreteStepsVDF  = nNumDiscreteStepsVDF;
         _dVeloMax             = dVeloMax;
    }

    // subdivision
	void SetSubdivisionProfiles(const unsigned int& nNumSlabs) {_nNumShellsProfiles = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionProfiles(const double& dSlabWidth) {_dShellWidthProfilesInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
	void SetSubdivisionVDF(const unsigned int& nNumSlabs) {_nNumShellsVDF = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
	void SetSubdivisionVDF(const double& dSlabWidth) {_dShellWidthVDFInit = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
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
    void WriteDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);
    void WriteDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep);

    void UpdateSlabParameters();

private:
    // reset local values
    void ResetLocalValuesProfiles();
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

    // parameters
    unsigned long _initSamplingProfiles;
    unsigned long _writeFrequencyProfiles;
    unsigned int _nNumShellsProfiles;

    double _dShellWidthProfiles;
    double _dShellWidthProfilesInit;
    double _dShellVolumeProfiles;
    double* _dShellMidpointsProfiles;

    unsigned long* _nNumMoleculesSumLocal;
    unsigned long* _nNumMoleculesSumCumulativeLocal;
    unsigned long* _nNumMoleculesSumGlobal;
    unsigned long* _nNumMoleculesSumCumulativeGlobal;

    // counting j+, j-
    unsigned long* _nNumMoleculesPlusSumCumulativeLocal;
    unsigned long* _nNumMoleculesMinusSumCumulativeLocal;
    unsigned long* _nNumMoleculesPlusSumCumulativeGlobal;
    unsigned long* _nNumMoleculesMinusSumCumulativeGlobal;

    // [position][x,y,z-component]
    double** _dVelocityComponentSumsLocal;
    double** _dVelocityComponentSumsCumulativeLocal;
    double** _dSquaredVelocityComponentSumsLocal;
    double** _dSquaredVelocityComponentSumsCumulativeLocal;
    double** _dVelocityComponentPlusSumsCumulativeLocal;
    double** _dVelocityComponentMinusSumsCumulativeLocal;

    double** _dVelocityComponentSumsGlobal;
    double** _dVelocityComponentSumsCumulativeGlobal;
    double** _dSquaredVelocityComponentSumsGlobal;
    double** _dSquaredVelocityComponentSumsCumulativeGlobal;
    double** _dVelocityComponentPlusSumsCumulativeGlobal;
    double** _dVelocityComponentMinusSumsCumulativeGlobal;

    // output values
    double** _dDriftVelocityGlobal;
    double** _dDriftVelocityAverageGlobal;
    double** _dTemperatureComponentGlobal;
    double** _dTemperatureComponentAverageGlobal;
    double*  _dDensityGlobal;
    double*  _dDensityAverageGlobal;
    double** _dDriftVelocityPlusAverageGlobal;
    double** _dDriftVelocityMinusAverageGlobal;

    // --- temperature profiles ---


    // --- VDF ---

    // parameters
    unsigned long _initSamplingVDF;
    unsigned long _writeFrequencyVDF;
    unsigned int _nNumShellsVDF;
    unsigned int _nNumDiscreteStepsVDF;
    double _dVeloMax;

    // control
    bool _bDiscretisationDoneVDF;
    bool _SamplingEnabledVDF;

    double  _dShellWidthVDF;
    double  _dShellWidthVDFInit;
    double* _dShellMidpointsVDF;
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

    // --- VDF ---


    // --- componentwise temperature ---

    unsigned long** _nNumMoleculesCompLocal;
    unsigned long** _nNumMoleculesCompGlobal;

    unsigned long** _nRotDOFCompLocal;
    unsigned long** _nRotDOFCompGlobal;

    double** _d2EkinTransCompLocal;
    double** _d2EkinTransCompGlobal;

    double** _d2EkinRotCompLocal;
    double** _d2EkinRotCompGlobal;

    double** _dTemperatureCompGlobal;
    double** _dDensityCompGlobal;


    // --- componentwise; x,y,z ; j+/j-; slabwise; rho, vx,vy,vz; Fx,Fy,Fz ---

    // [component][position]
    unsigned long** _nNumMoleculesCompLocal_py;
    unsigned long** _nNumMoleculesCompLocal_ny;
    unsigned long** _nNumMoleculesCompGlobal_py;
    unsigned long** _nNumMoleculesCompGlobal_ny;

    // [component][position]
    double** _dDensityCompGlobal_py;
    double** _dDensityCompGlobal_ny;

    // [component][vx,vy,vz][position]
    double*** _dVelocityCompLocal_py;
    double*** _dVelocityCompLocal_ny;
    double*** _dVelocityCompGlobal_py;
    double*** _dVelocityCompGlobal_ny;

    // [component][fx,fy,fz][position]
    double*** _dForceCompLocal_py;
    double*** _dForceCompLocal_ny;
    double*** _dForceCompGlobal_py;
    double*** _dForceCompGlobal_ny;
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
