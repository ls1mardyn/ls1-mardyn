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

class Domain;
class DomainDecompBase;
class Molecule;
class RegionSampling;
class SampleRegion
{
public:
    SampleRegion(RegionSampling* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned short nID);
    ~SampleRegion();

    // standard region methods
    double* GetLowerCorner() {return _dLowerCorner;}
    double* GetUpperCorner() {return _dUpperCorner;}
    void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal;}
    void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal;}
    double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
    void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
    unsigned short GetID() {return _nID;}

    // set parameters
    void SetParamProfiles( unsigned long initSamplingProfiles, unsigned long writeFrequencyProfiles,
                           unsigned int nNumShellsProfiles)
    {
        _initSamplingProfiles   = initSamplingProfiles;
        _writeFrequencyProfiles = writeFrequencyProfiles;
        _nNumShellsProfiles     = nNumShellsProfiles;
    }
    void SetParamVDF( unsigned long initSamplingVDF, unsigned long writeFrequencyVDF,
                      unsigned int nNumShellsVDF, unsigned int nNumDiscreteStepsVDF, double dVeloMax)
    {
        _initSamplingVDF       = initSamplingVDF;
        _writeFrequencyVDF     = writeFrequencyVDF;
        _nNumShellsVDF         = nNumShellsVDF;
        _nNumDiscreteStepsVDF  = nNumDiscreteStepsVDF;
         _dVeloMax             = dVeloMax;
    }

    // init sampling
    void InitSamplingProfiles(int nDimension, Domain* domain);
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

private:
    // reset local values
    void ResetLocalValuesProfiles();
    void ResetLocalValuesVDF();

private:

    double _dLowerCorner[3];
    double _dUpperCorner[3];
    unsigned short _nID;

    RegionSampling* _parent;


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


class RegionSampling
{
public:
    RegionSampling(Domain* domain);
    ~RegionSampling();

    void AddRegion(double dLowerCorner[3], double dUpperCorner[3] );
    int GetNumRegions() {return _vecSampleRegions.size();}
    SampleRegion* GetSampleRegion(unsigned short nRegionID) {return &(_vecSampleRegions.at(nRegionID-1) ); }  // vector index starts with 0, region index with 1

    void Init(Domain* domain);
    void DoSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);

    void WriteData(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain);

    unsigned short GetNumComponents() {return _nNumComponents;}

private:
    std::vector<SampleRegion> _vecSampleRegions;

    unsigned long _initSamplingProfiles;
    unsigned long _writeFrequencyProfiles;
    unsigned long _initSamplingVDF;
    unsigned long _writeFrequencyVDF;

    Domain* _domain;
    unsigned short _nNumComponents;
};


#endif /* REGIONSAMPLING_H_ */
