/*
 * RegionSampling.cpp
 *
 *  Created on: 18.03.2015
 *      Author: mheinen
 */

#include "RegionSampling.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;


SampleRegion::SampleRegion(RegionSampling* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned short nID)
{
    // set lower and upper corner
    for(unsigned short d=0; d<3; ++d)
    {
        _dLowerCorner[d] = dLowerCorner[d];
        _dUpperCorner[d] = dUpperCorner[d];
    }

    // set ID
    _nID = nID;

    // init control
    _bDiscretisationDoneProfiles = false;
    _bDiscretisationDoneVDF = false;

    // store parent pointer
    _parent = parent;
}


SampleRegion::~SampleRegion()
{

}

void SampleRegion::InitSamplingProfiles(int nDimension, Domain* domain)
{
    // shell width
    double dNumShellsTemperature = (double) _nNumShellsProfiles;
    _dShellWidthProfiles = this->GetWidth(nDimension) / dNumShellsTemperature;

    // shell volume
    double dArea;

    switch(nDimension)
    {
    case RS_DIMENSION_X:
        dArea = domain->getGlobalLength(1) * domain->getGlobalLength(2);
        break;

    case RS_DIMENSION_Y:
        dArea = domain->getGlobalLength(0) * domain->getGlobalLength(2);
        break;

    case RS_DIMENSION_Z:
        dArea = domain->getGlobalLength(0) * domain->getGlobalLength(1);
        break;
    default:
        dArea = domain->getGlobalLength(0) * domain->getGlobalLength(2);
    }
    _dShellVolumeProfiles = _dShellWidthProfiles * dArea;


    // discrete values: shell midpoints, velocity values
    _dShellMidpointsProfiles = new double[_nNumShellsProfiles];

    // number of molecules
    _nNumMoleculesSumLocal            = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesSumCumulativeLocal  = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesSumGlobal           = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesSumCumulativeGlobal = new unsigned long[_nNumShellsProfiles];

    // j+, j-
    _nNumMoleculesPlusSumCumulativeLocal = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesPlusSumCumulativeGlobal = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesMinusSumCumulativeLocal = new unsigned long[_nNumShellsProfiles];
    _nNumMoleculesMinusSumCumulativeGlobal = new unsigned long[_nNumShellsProfiles];

    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
        _dShellMidpointsProfiles[s] = 0.;

        // number of molecules
        _nNumMoleculesSumLocal[s] = 0;
        _nNumMoleculesSumCumulativeLocal[s] = 0;
        _nNumMoleculesSumGlobal[s] = 0;
        _nNumMoleculesSumCumulativeGlobal[s] = 0;

        // j+, j-
        _nNumMoleculesPlusSumCumulativeLocal[s] = 0;
        _nNumMoleculesPlusSumCumulativeGlobal[s] = 0;
        _nNumMoleculesMinusSumCumulativeLocal[s] = 0;
        _nNumMoleculesMinusSumCumulativeGlobal[s] = 0;
    }

    // output values

    // local
    _dVelocityComponentSumsLocal                  = new double*[_nNumShellsProfiles];
    _dVelocityComponentSumsCumulativeLocal        = new double*[_nNumShellsProfiles];
    _dSquaredVelocityComponentSumsLocal           = new double*[_nNumShellsProfiles];
    _dSquaredVelocityComponentSumsCumulativeLocal = new double*[_nNumShellsProfiles];

    // global
    _dVelocityComponentSumsGlobal                  = new double*[_nNumShellsProfiles];
    _dVelocityComponentSumsCumulativeGlobal        = new double*[_nNumShellsProfiles];
    _dSquaredVelocityComponentSumsGlobal           = new double*[_nNumShellsProfiles];
    _dSquaredVelocityComponentSumsCumulativeGlobal = new double*[_nNumShellsProfiles];

    _dDriftVelocityGlobal               = new double*[_nNumShellsProfiles];
    _dDriftVelocityAverageGlobal        = new double*[_nNumShellsProfiles];
    _dTemperatureComponentGlobal        = new double*[_nNumShellsProfiles];
    _dTemperatureComponentAverageGlobal = new double*[_nNumShellsProfiles];

    _dDensityGlobal        = new double[_nNumShellsProfiles];
    _dDensityAverageGlobal = new double[_nNumShellsProfiles];


    // j+, j-
    _dVelocityComponentPlusSumsCumulativeLocal = new double*[_nNumShellsProfiles];
    _dVelocityComponentMinusSumsCumulativeLocal = new double*[_nNumShellsProfiles];
    _dVelocityComponentPlusSumsCumulativeGlobal = new double*[_nNumShellsProfiles];
    _dVelocityComponentMinusSumsCumulativeGlobal = new double*[_nNumShellsProfiles];

    _dDriftVelocityPlusAverageGlobal = new double*[_nNumShellsProfiles];
    _dDriftVelocityMinusAverageGlobal = new double*[_nNumShellsProfiles];


    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
        // local
        _dVelocityComponentSumsLocal[s]                  = new double[3];
        _dVelocityComponentSumsCumulativeLocal[s]        = new double[3];
        _dSquaredVelocityComponentSumsLocal[s]           = new double[3];
        _dSquaredVelocityComponentSumsCumulativeLocal[s] = new double[3];

        // global
        _dVelocityComponentSumsGlobal[s]                  = new double[3];
        _dVelocityComponentSumsCumulativeGlobal[s]        = new double[3];
        _dSquaredVelocityComponentSumsGlobal[s]           = new double[3];
        _dSquaredVelocityComponentSumsCumulativeGlobal[s] = new double[3];

        _dDriftVelocityGlobal[s]               = new double[3];
        _dDriftVelocityAverageGlobal[s]        = new double[3];
        _dTemperatureComponentGlobal[s]        = new double[3];
        _dTemperatureComponentAverageGlobal[s] = new double[3];

        // j+, j-
        _dVelocityComponentPlusSumsCumulativeLocal[s] = new double[3];
        _dVelocityComponentMinusSumsCumulativeLocal[s] = new double[3];
        _dVelocityComponentPlusSumsCumulativeGlobal[s] = new double[3];
        _dVelocityComponentMinusSumsCumulativeGlobal[s] = new double[3];

        _dDriftVelocityPlusAverageGlobal[s] = new double[3];
        _dDriftVelocityMinusAverageGlobal[s] = new double[3];
    }


    // init values
    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
        // init density
        _dDensityGlobal[s] = 0.;
        _dDensityAverageGlobal[s] = 0.;

        for(unsigned int d = 0; d < 3; ++d)  // 3 diminsions: x, y, z
        {
            // local
            _dVelocityComponentSumsLocal[s][d]                  = 0.;
            _dVelocityComponentSumsCumulativeLocal[s][d]        = 0.;
            _dSquaredVelocityComponentSumsLocal[s][d]           = 0.;
            _dSquaredVelocityComponentSumsCumulativeLocal[s][d] = 0.;

            // global
            _dVelocityComponentSumsGlobal[s][d]                  = 0.;
            _dVelocityComponentSumsCumulativeGlobal[s][d]        = 0.;
            _dSquaredVelocityComponentSumsGlobal[s][d]           = 0.;
            _dSquaredVelocityComponentSumsCumulativeGlobal[s][d] = 0.;

            _dDriftVelocityGlobal[s][d]               = 0.;
            _dDriftVelocityAverageGlobal[s][d]        = 0.;
            _dTemperatureComponentGlobal[s][d]        = 0.;
            _dTemperatureComponentAverageGlobal[s][d] = 0.;

            // j+, j-
            _dVelocityComponentPlusSumsCumulativeLocal[s][d] = 0.;
            _dVelocityComponentMinusSumsCumulativeLocal[s][d] = 0.;
            _dVelocityComponentPlusSumsCumulativeGlobal[s][d] = 0.;
            _dVelocityComponentMinusSumsCumulativeGlobal[s][d] = 0.;

            _dDriftVelocityPlusAverageGlobal[s][d]  = 0.;
            _dDriftVelocityMinusAverageGlobal[s][d] = 0.;
        }
    }

//    cout << "_initSamplingProfiles = " << _initSamplingProfiles << endl;
//    cout << "_writeFrequencyProfiles = " << _writeFrequencyProfiles << endl;
//    cout << "_nNumShellsProfiles = " << _nNumShellsProfiles << endl;


    // componentwise temperature
    unsigned int nNumComponents;
    nNumComponents = domain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

    _nNumMoleculesCompLocal  = new unsigned long* [nNumComponents];
    _nNumMoleculesCompGlobal = new unsigned long* [nNumComponents];
    _nRotDOFCompLocal  = new unsigned long* [nNumComponents];
    _nRotDOFCompGlobal = new unsigned long* [nNumComponents];

    _d2EkinTransCompLocal = new double* [nNumComponents];
    _d2EkinTransCompGlobal = new double* [nNumComponents];
    _d2EkinRotCompLocal = new double* [nNumComponents];
    _d2EkinRotCompGlobal = new double* [nNumComponents];

    _dTemperatureCompGlobal = new double* [nNumComponents];
    _dDensityCompGlobal = new double* [nNumComponents];

    for(unsigned short c=0; c < nNumComponents; ++c)
    {
        _nNumMoleculesCompLocal[c]  = new unsigned long[_nNumShellsProfiles];
        _nNumMoleculesCompGlobal[c] = new unsigned long[_nNumShellsProfiles];
        _nRotDOFCompLocal[c]  = new unsigned long[_nNumShellsProfiles];
        _nRotDOFCompGlobal[c] = new unsigned long[_nNumShellsProfiles];

        _d2EkinTransCompLocal[c]  = new double[_nNumShellsProfiles];
        _d2EkinTransCompGlobal[c] = new double[_nNumShellsProfiles];
        _d2EkinRotCompLocal[c]  = new double[_nNumShellsProfiles];
        _d2EkinRotCompGlobal[c] = new double[_nNumShellsProfiles];

        _dTemperatureCompGlobal[c] = new double[_nNumShellsProfiles];
        _dDensityCompGlobal[c] = new double[_nNumShellsProfiles];
    }

    // init
    for(unsigned short c=0; c < nNumComponents; ++c)
    {
        for(unsigned short s=0; s < _nNumShellsProfiles; ++s)
        {
            _nNumMoleculesCompLocal[c][s]  = 0;
            _nNumMoleculesCompGlobal[c][s]  = 0;
            _nRotDOFCompLocal[c][s]  = 0;
            _nRotDOFCompGlobal[c][s]  = 0;

            _d2EkinTransCompLocal[c][s]  = 0.;
            _d2EkinTransCompGlobal[c][s]  = 0.;
            _d2EkinRotCompLocal[c][s]  = 0.;
            _d2EkinRotCompGlobal[c][s]  = 0.;

            _dTemperatureCompGlobal[c][s]  = 0.;
            _dDensityCompGlobal[c][s]  = 0.;
        }
    }

    this->DoDiscretisationProfiles(RS_DIMENSION_Y);
}


void SampleRegion::InitSamplingVDF(int nDimension)
{
    // shell width
    double dNumShellsVDF = (double) _nNumShellsVDF;
    _dShellWidthVDF = this->GetWidth(nDimension) / dNumShellsVDF;

    // discrete values: shell midpoints, velocity values
    _dShellMidpointsVDF = new double[_nNumShellsVDF];

    for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
    {
        _dShellMidpointsVDF[s] = 0.;
    }

    _dDiscreteVelocityValues = new double[_nNumDiscreteStepsVDF];

    for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
    {
        _dDiscreteVelocityValues[v] = 0.;
    }

    // local
    _veloDistrMatrixLocal_py_abs = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_pvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_pvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_pvz = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_nvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_nvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_py_nvz = new unsigned long*[_nNumShellsVDF];

    _veloDistrMatrixLocal_ny_abs = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_pvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_pvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_pvz = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_nvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_nvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixLocal_ny_nvz = new unsigned long*[_nNumShellsVDF];

    // global
    _veloDistrMatrixGlobal_py_abs = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_pvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_pvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_pvz = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_nvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_nvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_py_nvz = new unsigned long*[_nNumShellsVDF];

    _veloDistrMatrixGlobal_ny_abs = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_pvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_pvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_pvz = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_nvx = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_nvy = new unsigned long*[_nNumShellsVDF];
    _veloDistrMatrixGlobal_ny_nvz = new unsigned long*[_nNumShellsVDF];


    for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
    {
        // local
        _veloDistrMatrixLocal_py_abs[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_pvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_pvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_pvz[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_nvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_nvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_py_nvz[s] = new unsigned long[_nNumDiscreteStepsVDF];

        _veloDistrMatrixLocal_ny_abs[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_pvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_pvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_pvz[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_nvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_nvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixLocal_ny_nvz[s] = new unsigned long[_nNumDiscreteStepsVDF];

        // global
        _veloDistrMatrixGlobal_py_abs[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_pvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_pvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_pvz[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_nvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_nvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_py_nvz[s] = new unsigned long[_nNumDiscreteStepsVDF];

        _veloDistrMatrixGlobal_ny_abs[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_pvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_pvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_pvz[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_nvx[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_nvy[s] = new unsigned long[_nNumDiscreteStepsVDF];
        _veloDistrMatrixGlobal_ny_nvz[s] = new unsigned long[_nNumDiscreteStepsVDF];

    }

    // init values
    for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
    {
        for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
        {
            // local
            _veloDistrMatrixLocal_py_abs[s][v] = 0;
            _veloDistrMatrixLocal_py_pvx[s][v] = 0;
            _veloDistrMatrixLocal_py_pvy[s][v] = 0;
            _veloDistrMatrixLocal_py_pvz[s][v] = 0;
            _veloDistrMatrixLocal_py_nvx[s][v] = 0;
            _veloDistrMatrixLocal_py_nvy[s][v] = 0;
            _veloDistrMatrixLocal_py_nvz[s][v] = 0;

            _veloDistrMatrixLocal_ny_abs[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvx[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvy[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvz[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvx[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvy[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvz[s][v] = 0;

            // global
            _veloDistrMatrixGlobal_py_abs[s][v] = 0;
            _veloDistrMatrixGlobal_py_pvx[s][v] = 0;
            _veloDistrMatrixGlobal_py_pvy[s][v] = 0;
            _veloDistrMatrixGlobal_py_pvz[s][v] = 0;
            _veloDistrMatrixGlobal_py_nvx[s][v] = 0;
            _veloDistrMatrixGlobal_py_nvy[s][v] = 0;
            _veloDistrMatrixGlobal_py_nvz[s][v] = 0;

            _veloDistrMatrixGlobal_ny_abs[s][v] = 0;
            _veloDistrMatrixGlobal_ny_pvx[s][v] = 0;
            _veloDistrMatrixGlobal_ny_pvy[s][v] = 0;
            _veloDistrMatrixGlobal_ny_pvz[s][v] = 0;
            _veloDistrMatrixGlobal_ny_nvx[s][v] = 0;
            _veloDistrMatrixGlobal_ny_nvy[s][v] = 0;
            _veloDistrMatrixGlobal_ny_nvz[s][v] = 0;

        }
    }
/*
    cout << "_initSamplingVDF = " << _initSamplingVDF << endl;
    cout << "_writeFrequencyVDF = " << _writeFrequencyVDF << endl;
    cout << "_nNumShellsVDF = " << _nNumShellsVDF << endl;
    cout << "_nNumDiscreteStepsVDF = " << _nNumDiscreteStepsVDF << endl;
*/

    // discrete velocity values
    this->DoDiscretisationVDF(RS_DIMENSION_Y);
}


void SampleRegion::DoDiscretisationProfiles(int nDimension)
{
    if(_bDiscretisationDoneProfiles == true)  // if allready done -> return
        return;

    double* dLowerCorner = this->GetLowerCorner();

    // calc shell midpoints
    for(unsigned int s = 0; s < _nNumShellsProfiles; s++)
    {
        _dShellMidpointsProfiles[s] = (s + 0.5) * _dShellWidthProfiles + dLowerCorner[nDimension];
    }

    _bDiscretisationDoneProfiles = true;
}


void SampleRegion::DoDiscretisationVDF(int nDimension)
{
    if(_bDiscretisationDoneVDF == true)  // if allready done -> return
        return;

//    double dVeloMax = _dVeloMax * 1.1;  // velocity discretisation, highest value with safety factor
    double dNumDiscreteStepsVDF = (double) _nNumDiscreteStepsVDF;
    double dDeltaVelo = _dVeloMax / dNumDiscreteStepsVDF;

    double* dLowerCorner = this->GetLowerCorner();

    // calc shell midpoints
    for(unsigned int s = 0; s < _nNumShellsVDF; s++)
    {
        _dShellMidpointsVDF[s] = (s + 0.5) * _dShellWidthVDF + dLowerCorner[nDimension];
    }


    // calc discrete velocity values
    for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
    {
        _dDiscreteVelocityValues[v] = dDeltaVelo * (v + 0.5);
    }

    _bDiscretisationDoneVDF = true;
}



void SampleRegion::SampleProfiles(Molecule* molecule, int nDimension)
{
    unsigned int nPosIndex;
    unsigned int nIndexMax = _nNumShellsProfiles - 1;

    // do not reset profile matrices here!!!
    // BUT: reset profile before calling this function!!!

    double v[3];

    // calc position index
    double* dLowerCorner = this->GetLowerCorner();
    double dPosRelative = molecule->r(nDimension) - dLowerCorner[nDimension];

    nPosIndex = (unsigned int) floor(dPosRelative / _dShellWidthProfiles);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
        return;

    v[0] = molecule->v(0);
    v[1] = molecule->v(1);
    v[2] = molecule->v(2);

    // component velocity sums (drift velocities)
    // [position][x,y,z-component]
    _dVelocityComponentSumsLocal[nPosIndex][0] += v[0];
    _dVelocityComponentSumsLocal[nPosIndex][1] += v[1];
    _dVelocityComponentSumsLocal[nPosIndex][2] += v[2];

    _dVelocityComponentSumsCumulativeLocal[nPosIndex][0] += v[0];
    _dVelocityComponentSumsCumulativeLocal[nPosIndex][1] += v[1];
    _dVelocityComponentSumsCumulativeLocal[nPosIndex][2] += v[2];


    // TODO: perhabs calc m * v_drift instead of v_drift
    // important when we have mixture with different molecule masses

    // squared component velocity sums (need to calculate kinetic energies)
    // [position][x,y,z-component]
    _dSquaredVelocityComponentSumsLocal[nPosIndex][0] += v[0]*v[0];
    _dSquaredVelocityComponentSumsLocal[nPosIndex][1] += v[1]*v[1];
    _dSquaredVelocityComponentSumsLocal[nPosIndex][2] += v[2]*v[2];

    _dSquaredVelocityComponentSumsCumulativeLocal[nPosIndex][0] += v[0]*v[0];
    _dSquaredVelocityComponentSumsCumulativeLocal[nPosIndex][1] += v[1]*v[1];
    _dSquaredVelocityComponentSumsCumulativeLocal[nPosIndex][2] += v[2]*v[2];

    // sum up molecules taken into account
    _nNumMoleculesSumLocal[nPosIndex]++;
    _nNumMoleculesSumCumulativeLocal[nPosIndex]++;

    // j+, j-
    bool bVelocityIsPlus = (v[1] > 0);

/*
    cout << "v[1] = " << v[1] << endl;
    cout << "bVelocityIsPlus = " << bVelocityIsPlus << endl;
    cout << "!bVelocityIsPlus = " << !bVelocityIsPlus << endl;
    cout << "bVelocityIsPlus * v[1] = " << bVelocityIsPlus * v[1] << endl;
    cout << "!bVelocityIsPlus * v[1] = " << !bVelocityIsPlus * v[1] << endl;
*/

    _nNumMoleculesPlusSumCumulativeLocal[nPosIndex]  +=  bVelocityIsPlus;
    _nNumMoleculesMinusSumCumulativeLocal[nPosIndex] += !bVelocityIsPlus;

    // y-direction
    _dVelocityComponentPlusSumsCumulativeLocal[nPosIndex][1]  +=  bVelocityIsPlus * v[1];
    _dVelocityComponentMinusSumsCumulativeLocal[nPosIndex][1] += !bVelocityIsPlus * v[1];


    // componentwise temperature
    unsigned int cid = molecule->componentid()+1;  // starts with 0

    _nNumMoleculesCompLocal[cid][nPosIndex]++;
    _nRotDOFCompLocal[cid][nPosIndex] += molecule->component()->getRotationalDegreesOfFreedom();

    molecule->calculate_mv2_Iw2(_d2EkinTransCompLocal[cid][nPosIndex], _d2EkinRotCompLocal[cid][nPosIndex]);

    // total
    _nNumMoleculesCompLocal[0][nPosIndex]++;
    _nRotDOFCompLocal[0][nPosIndex] += molecule->component()->getRotationalDegreesOfFreedom();

    molecule->calculate_mv2_Iw2(_d2EkinTransCompLocal[0][nPosIndex], _d2EkinRotCompLocal[0][nPosIndex]);

}


void SampleRegion::SampleVDF(Molecule* molecule, int nDimension)
{
    // return if discretisation is not done yet
    // reason: v_max has to be determined first (when not set manually)
    if(false == _bDiscretisationDoneVDF)
        return;

    unsigned int nPosIndex;
    unsigned int nVelocityIndex;
    unsigned int nIndexMax     = _nNumShellsVDF - 1;
    unsigned int nIndexMaxVelo = _nNumDiscreteStepsVDF - 1;

    double dVelocity;
//    double dMaxVelo = _dVeloMax * 1.1;  // velocity discretisation, highest value with safety factor

    // do not reset profile matrices here!!

    // calc position index
    double* dLowerCorner = this->GetLowerCorner();
    double dPosRelative = molecule->r(nDimension) - dLowerCorner[nDimension];

    nPosIndex = (unsigned int) floor(dPosRelative / _dShellWidthVDF);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)
        return;

    dVelocity = sqrt( molecule->v2() );
    nVelocityIndex = (unsigned int) floor(dVelocity / _dVeloMax * _nNumDiscreteStepsVDF);


    // respect finite resolution of velocity
    if(nVelocityIndex > nIndexMaxVelo)
    {
//        global_log->info() << "Molecule with: v > v_max detected: v = " << dVelocity << ", v_max = " << dMaxVelo << endl;
        return;
    }

    // calculate velocity vector indices for velocity components
    double v_dim[3];
    unsigned int naVelocityIndex[3];

    for(unsigned int d=0; d<3; ++d)
    {
        v_dim[d] = molecule->v(d);
        naVelocityIndex[d] = (unsigned int) floor( fabs( v_dim[d] ) / _dVeloMax * _nNumDiscreteStepsVDF);
    }

    int pv[3];
    int nv[3];

    for(unsigned int d=0; d<3; ++d)
    {
        if(v_dim[d] > 0.)
        {
            pv[d] = 1;
            nv[d] = 0;
        }
        else
        {
            pv[d] = 0;
            nv[d] = 1;
        }
    }

    // positive y-direction
    if(v_dim[1] > 0.)
    {
        // absolute
        _veloDistrMatrixLocal_py_abs[nPosIndex][nVelocityIndex]++;

        _veloDistrMatrixLocal_py_pvx[nPosIndex][naVelocityIndex[0] ] += pv[0];
        _veloDistrMatrixLocal_py_pvy[nPosIndex][naVelocityIndex[1] ] += pv[1];
        _veloDistrMatrixLocal_py_pvz[nPosIndex][naVelocityIndex[2] ] += pv[2];

        _veloDistrMatrixLocal_py_nvx[nPosIndex][naVelocityIndex[0] ] += nv[0];
        _veloDistrMatrixLocal_py_nvy[nPosIndex][naVelocityIndex[1] ] += nv[1];
        _veloDistrMatrixLocal_py_nvz[nPosIndex][naVelocityIndex[2] ] += nv[2];
    }
    // negative y-direction
    else
    {
        // absolute
        _veloDistrMatrixLocal_ny_abs[nPosIndex][nVelocityIndex]++;

        _veloDistrMatrixLocal_ny_pvx[nPosIndex][naVelocityIndex[0] ] += pv[0];
        _veloDistrMatrixLocal_ny_pvy[nPosIndex][naVelocityIndex[1] ] += pv[1];
        _veloDistrMatrixLocal_ny_pvz[nPosIndex][naVelocityIndex[2] ] += pv[2];

        _veloDistrMatrixLocal_ny_nvx[nPosIndex][naVelocityIndex[0] ] += nv[0];
        _veloDistrMatrixLocal_ny_nvy[nPosIndex][naVelocityIndex[1] ] += nv[1];
        _veloDistrMatrixLocal_ny_nvz[nPosIndex][naVelocityIndex[2] ] += nv[2];
    }
}


void SampleRegion::CalcGlobalValuesProfiles(DomainDecompBase* domainDecomp, Domain* domain)
{
    int ownRank = domainDecomp->getRank();
//  int numprocs = domainDecomp->getNumProcs();

    double dInvertDOF;
    double dInvertDOF_plus, dInvertDOF_minus;
    double dInvertShellVolume = 1. / _dShellVolumeProfiles;

    double dInvertNumSamples = (double) (_writeFrequencyProfiles);  // TODO: perhabs in future different from writeFrequency
    dInvertNumSamples = 1. / dInvertNumSamples;


    // <<< important !! >>>
    // reset global values before reduce operation, cause this function can be called more than once
    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
        _nNumMoleculesSumGlobal[s] = 0;
        _nNumMoleculesSumCumulativeGlobal[s] = 0;

        for(unsigned short d = 0; d<3; ++d)
        {
            _dVelocityComponentSumsGlobal[s][d] = 0.;
            _dVelocityComponentSumsCumulativeGlobal[s][d] = 0.;
            _dSquaredVelocityComponentSumsGlobal[s][d] = 0.;
            _dSquaredVelocityComponentSumsCumulativeGlobal[s][d] = 0.;
        }
    }


    // perform reduce operation, process further calculations
#ifdef ENABLE_MPI

    MPI_Reduce( _nNumMoleculesSumLocal,           _nNumMoleculesSumGlobal,           _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( _nNumMoleculesSumCumulativeLocal, _nNumMoleculesSumCumulativeGlobal, _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( _nNumMoleculesPlusSumCumulativeLocal, _nNumMoleculesPlusSumCumulativeGlobal, _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( _nNumMoleculesMinusSumCumulativeLocal, _nNumMoleculesMinusSumCumulativeGlobal, _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

#endif


    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
    #ifdef ENABLE_MPI

        MPI_Reduce( _dVelocityComponentSumsLocal[s],                  _dVelocityComponentSumsGlobal[s],                  3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _dVelocityComponentSumsCumulativeLocal[s],        _dVelocityComponentSumsCumulativeGlobal[s],        3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _dSquaredVelocityComponentSumsLocal[s],           _dSquaredVelocityComponentSumsGlobal[s],           3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _dSquaredVelocityComponentSumsCumulativeLocal[s], _dSquaredVelocityComponentSumsCumulativeGlobal[s], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _dVelocityComponentPlusSumsCumulativeLocal[s],        _dVelocityComponentPlusSumsCumulativeGlobal[s],        3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _dVelocityComponentMinusSumsCumulativeLocal[s],        _dVelocityComponentMinusSumsCumulativeGlobal[s],        3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    #endif

        if(ownRank != 0)
            continue;

        // actual values
        dInvertDOF = (double) (_nNumMoleculesSumGlobal[s]);
        dInvertDOF = 1. / dInvertDOF;

        _dDriftVelocityGlobal[s][0] = _dVelocityComponentSumsGlobal[s][0] * dInvertDOF;
        _dDriftVelocityGlobal[s][1] = _dVelocityComponentSumsGlobal[s][1] * dInvertDOF;
        _dDriftVelocityGlobal[s][2] = _dVelocityComponentSumsGlobal[s][2] * dInvertDOF;

        // TODO: has to be multiplied by molecule mass
        _dTemperatureComponentGlobal[s][0] = _dSquaredVelocityComponentSumsGlobal[s][0] * dInvertDOF;
        _dTemperatureComponentGlobal[s][1] = _dSquaredVelocityComponentSumsGlobal[s][1] * dInvertDOF;
        _dTemperatureComponentGlobal[s][2] = _dSquaredVelocityComponentSumsGlobal[s][2] * dInvertDOF;

        _dDensityGlobal[s] = _nNumMoleculesSumGlobal[s] * dInvertShellVolume;

        // average values
        dInvertDOF = (double) (_nNumMoleculesSumCumulativeGlobal[s]);
        dInvertDOF = 1. / dInvertDOF;

        dInvertDOF_plus = (double) (_nNumMoleculesPlusSumCumulativeGlobal[s]);
        dInvertDOF_plus = 1. / dInvertDOF_plus;
        dInvertDOF_minus = (double) (_nNumMoleculesMinusSumCumulativeGlobal[s]);
        dInvertDOF_minus = 1. / dInvertDOF_minus;

//        if(ownRank == 0)
//        {
//            cout << "[" << ownRank << "]: Region: " << this->GetID() << endl;
//            cout << "[" << ownRank << "]: _nNumMoleculesSumCumulativeLocal[" << s << "] = " << _nNumMoleculesSumCumulativeLocal[s] << endl;
//            cout << "[" << ownRank << "]: _nNumMoleculesSumCumulativeGlobal[" << s << "] = " << _nNumMoleculesSumCumulativeGlobal[s] << endl;
//            cout << "[" << ownRank << "]: dInvertDOF = " << dInvertDOF << endl;
//            cout << "[" << ownRank << "]: _dVelocityComponentSumsCumulativeGlobal[" << s << "][0] = " << _dVelocityComponentSumsCumulativeGlobal[s][0] << endl;
//        }

        _dDriftVelocityAverageGlobal[s][0] = _dVelocityComponentSumsCumulativeGlobal[s][0] * dInvertDOF;
        _dDriftVelocityAverageGlobal[s][1] = _dVelocityComponentSumsCumulativeGlobal[s][1] * dInvertDOF;
        _dDriftVelocityAverageGlobal[s][2] = _dVelocityComponentSumsCumulativeGlobal[s][2] * dInvertDOF;

        // j+, j-
        _dDriftVelocityPlusAverageGlobal[s][0] = _dVelocityComponentPlusSumsCumulativeGlobal[s][0] * dInvertDOF_plus;
        _dDriftVelocityPlusAverageGlobal[s][1] = _dVelocityComponentPlusSumsCumulativeGlobal[s][1] * dInvertDOF_plus;
        _dDriftVelocityPlusAverageGlobal[s][2] = _dVelocityComponentPlusSumsCumulativeGlobal[s][2] * dInvertDOF_plus;

        _dDriftVelocityMinusAverageGlobal[s][0] = _dVelocityComponentMinusSumsCumulativeGlobal[s][0] * dInvertDOF_minus;
        _dDriftVelocityMinusAverageGlobal[s][1] = _dVelocityComponentMinusSumsCumulativeGlobal[s][1] * dInvertDOF_minus;
        _dDriftVelocityMinusAverageGlobal[s][2] = _dVelocityComponentMinusSumsCumulativeGlobal[s][2] * dInvertDOF_minus;

//        if(ownRank == 0)
//        {
//            cout << "[" << ownRank << "]: _dDriftVelocityAverageGlobal[" << s << "][0] = " << _dDriftVelocityAverageGlobal[s][0] << endl;
//        }

        // TODO: has to be multiplied by molecule mass
        _dTemperatureComponentAverageGlobal[s][0] = _dSquaredVelocityComponentSumsCumulativeGlobal[s][0] * dInvertDOF;
        _dTemperatureComponentAverageGlobal[s][1] = _dSquaredVelocityComponentSumsCumulativeGlobal[s][1] * dInvertDOF;
        _dTemperatureComponentAverageGlobal[s][2] = _dSquaredVelocityComponentSumsCumulativeGlobal[s][2] * dInvertDOF;

        _dDensityAverageGlobal[s] = _nNumMoleculesSumCumulativeGlobal[s] * dInvertShellVolume * dInvertNumSamples;
    }


    // componentwise temperature
    unsigned short nNumComponents = domain->getNumberOfComponents()+1;

#ifdef ENABLE_MPI
    for(unsigned int c = 0; c < nNumComponents; ++c)
    {
        MPI_Reduce( _nNumMoleculesCompLocal[c], _nNumMoleculesCompGlobal[c], _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _nRotDOFCompLocal[c],       _nRotDOFCompGlobal[c],       _nNumShellsProfiles, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _d2EkinTransCompLocal[c],   _d2EkinTransCompGlobal[c],   _nNumShellsProfiles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( _d2EkinRotCompLocal[c],     _d2EkinRotCompGlobal[c],     _nNumShellsProfiles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif



    for(unsigned int c = 0; c < nNumComponents; ++c)
    {
        for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
        {
            double dInvertDOF = (double) (1. / (3. * _nNumMoleculesCompGlobal[c][s] + _nRotDOFCompGlobal[c][s] ) );
            _dTemperatureCompGlobal[c][s] = (_d2EkinTransCompGlobal[c][s] + _d2EkinRotCompGlobal[c][s] ) * dInvertDOF;

            // density componentwise
            _dDensityCompGlobal[c][s] = _nNumMoleculesCompGlobal[c][s] * dInvertShellVolume * dInvertNumSamples;
        }
    }

}



void SampleRegion::CalcGlobalValuesVDF()
{
    #ifdef ENABLE_MPI

        for(unsigned int s = 0; s < _nNumShellsVDF; s++)
        {
            // positive y-direction

            // absolute
            MPI_Reduce( _veloDistrMatrixLocal_py_abs[s], _veloDistrMatrixGlobal_py_abs[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

            MPI_Reduce( _veloDistrMatrixLocal_py_pvx[s], _veloDistrMatrixGlobal_py_pvx[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_py_pvy[s], _veloDistrMatrixGlobal_py_pvy[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_py_pvz[s], _veloDistrMatrixGlobal_py_pvz[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

            MPI_Reduce( _veloDistrMatrixLocal_py_nvx[s], _veloDistrMatrixGlobal_py_nvx[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_py_nvy[s], _veloDistrMatrixGlobal_py_nvy[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_py_nvz[s], _veloDistrMatrixGlobal_py_nvz[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


            // negative y-direction

            // absolute
            MPI_Reduce( _veloDistrMatrixLocal_ny_abs[s], _veloDistrMatrixGlobal_ny_abs[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

            MPI_Reduce( _veloDistrMatrixLocal_ny_pvx[s], _veloDistrMatrixGlobal_ny_pvx[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_ny_pvy[s], _veloDistrMatrixGlobal_ny_pvy[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_ny_pvz[s], _veloDistrMatrixGlobal_ny_pvz[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

            MPI_Reduce( _veloDistrMatrixLocal_ny_nvx[s], _veloDistrMatrixGlobal_ny_nvx[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_ny_nvy[s], _veloDistrMatrixGlobal_ny_nvy[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce( _veloDistrMatrixLocal_ny_nvz[s], _veloDistrMatrixGlobal_ny_nvz[s], _nNumDiscreteStepsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        }
    #else
        for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
        {
            for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
            {
                // positive y-direction
                _veloDistrMatrixGlobal_py_abs[s][v] = _veloDistrMatrixLocal_py_abs[s][v];
                _veloDistrMatrixGlobal_py_pvx[s][v] = _veloDistrMatrixLocal_py_pvx[s][v];
                _veloDistrMatrixGlobal_py_pvy[s][v] = _veloDistrMatrixLocal_py_pvy[s][v];
                _veloDistrMatrixGlobal_py_pvz[s][v] = _veloDistrMatrixLocal_py_pvz[s][v];
                _veloDistrMatrixGlobal_py_nvx[s][v] = _veloDistrMatrixLocal_py_nvx[s][v];
                _veloDistrMatrixGlobal_py_nvy[s][v] = _veloDistrMatrixLocal_py_nvy[s][v];
                _veloDistrMatrixGlobal_py_nvz[s][v] = _veloDistrMatrixLocal_py_nvz[s][v];

                // negative y-direction
                _veloDistrMatrixGlobal_ny_abs[s][v] = _veloDistrMatrixLocal_ny_abs[s][v];
                _veloDistrMatrixGlobal_ny_pvx[s][v] = _veloDistrMatrixLocal_ny_pvx[s][v];
                _veloDistrMatrixGlobal_ny_pvy[s][v] = _veloDistrMatrixLocal_ny_pvy[s][v];
                _veloDistrMatrixGlobal_ny_pvz[s][v] = _veloDistrMatrixLocal_ny_pvz[s][v];
                _veloDistrMatrixGlobal_ny_nvx[s][v] = _veloDistrMatrixLocal_ny_nvx[s][v];
                _veloDistrMatrixGlobal_ny_nvy[s][v] = _veloDistrMatrixLocal_ny_nvy[s][v];
                _veloDistrMatrixGlobal_ny_nvz[s][v] = _veloDistrMatrixLocal_ny_nvz[s][v];
            }
        }
    #endif

}




void SampleRegion::WriteDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain)
{
    // sampling starts after initial timestep (_initSamplingVDF) and with respect to write frequency (_writeFrequencyVDF)
    if( simstep <= _initSamplingProfiles )
        return;

    if ( (simstep - _initSamplingProfiles) % _writeFrequencyProfiles != 0 )
        return;


    // calc global values
    this->CalcGlobalValuesProfiles(domainDecomp, domain);

    // reset local values
    this->ResetLocalValuesProfiles();


    // writing .dat-files
    std::stringstream outputstream;
    std::stringstream filenamestream;
    filenamestream << "T-rho_region" << this->GetID() << "_TS" << simstep << ".dat";

    #ifdef ENABLE_MPI
        int rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank== 0)
        {
    #endif



            // header
            //outputstream << "           pos         v_d,x         v_d,y         v_d,z            Tx            Ty            Tz           rho        v_d,y+        v_d,y-";
            outputstream << "           pos         v_d,x         v_d,y         v_d,z            Tx            Ty            Tz";
            outputstream << "           rho        v_d,y+        v_d,y-          DOF+          DOF-       DOF_ges";
            outputstream << endl;

            // data
            double dVal;

            for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
            {
                outputstream << std::setw(14) << std::setprecision(6) << _dShellMidpointsProfiles[s];

                // drift x
                dVal = _dDriftVelocityAverageGlobal[s][0];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // drift y
                dVal = _dDriftVelocityAverageGlobal[s][1];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // drift z
                dVal = _dDriftVelocityAverageGlobal[s][2];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;


                // temperature Tx
                dVal = _dTemperatureComponentAverageGlobal[s][0];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // temperature Ty
                dVal = _dTemperatureComponentAverageGlobal[s][1];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // temperature Tz
                dVal = _dTemperatureComponentAverageGlobal[s][2];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;


                // density rho
                dVal = _dDensityAverageGlobal[s];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;


                // drift y+
                dVal = _dDriftVelocityPlusAverageGlobal[s][1];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // drift y-
                dVal = _dDriftVelocityMinusAverageGlobal[s][1];

                if( abs(dVal) < 0.00001 )
                    outputstream << std::setw(14) << fixed << std::setprecision(10) << dVal;
                else
                    outputstream << std::setw(14) << std::setprecision(6) << dVal;

                // DOF+, DOF-, DOF_ges
                outputstream << std::setw(14) << _nNumMoleculesPlusSumCumulativeGlobal[s];
                outputstream << std::setw(14) << _nNumMoleculesMinusSumCumulativeGlobal[s];
                outputstream << std::setw(14) << _nNumMoleculesSumCumulativeGlobal[s];

                outputstream << endl;
            }

            // Datei zum schreiben öffnen, daten schreiben
            ofstream fileout(filenamestream.str().c_str(), ios::out);
            fileout << outputstream.str();
            fileout.close();

            // global_log->info() << "files closed." << endl;

    #ifdef ENABLE_MPI
        }
    #endif


    // writing .dat-files
    std::stringstream outputstream_comp;
    std::stringstream filenamestream_comp;
    filenamestream_comp << "T-rho_comp_region" << this->GetID() << "_TS" << simstep << ".dat";

    #ifdef ENABLE_MPI
        rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank== 0)
        {
    #endif

            unsigned int nNumComponents;
            nNumComponents = domain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

            // header
            outputstream_comp << "           pos";

            // temperature
            for(unsigned short c = 0; c < nNumComponents; ++c)
            {
                outputstream_comp << "          T[" << c << "]";
            }

            // density
            for(unsigned short c = 0; c < nNumComponents; ++c)
            {
                outputstream_comp << "        rho[" << c << "]";
            }

            outputstream_comp << endl;

            // data
            double dVal;

            for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
            {

                outputstream_comp << std::setw(14) << std::setprecision(6) << _dShellMidpointsProfiles[s];

                for(unsigned short c = 0; c < nNumComponents; ++c)
                {
                    // temperature componentwise
                    dVal = _dTemperatureCompGlobal[c][s];

                    if( abs(dVal) < 0.00001 )
                        outputstream_comp << std::setw(14) << fixed << std::setprecision(10) << dVal;
                    else
                        outputstream_comp << std::setw(14) << std::setprecision(6) << dVal;
                }

                for(unsigned short c = 0; c < nNumComponents; ++c)
                {
                    // temperature componentwise
                    dVal = _dDensityCompGlobal[c][s];

                    if( abs(dVal) < 0.00001 )
                        outputstream_comp << std::setw(14) << fixed << std::setprecision(10) << dVal;
                    else
                        outputstream_comp << std::setw(14) << std::setprecision(6) << dVal;
                }

                outputstream_comp << endl;
            }

            // Datei zum schreiben öffnen, daten schreiben
            ofstream fileout(filenamestream_comp.str().c_str(), ios::out);
            fileout << outputstream_comp.str();
            fileout.close();

            // global_log->info() << "files closed." << endl;

    #ifdef ENABLE_MPI
        }
    #endif
}






void SampleRegion::WriteDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    // sampling starts after initial timestep (_initSamplingVDF) and with respect to write frequency (_writeFrequencyVDF)
    if( simstep <= _initSamplingVDF )
        return;

    if ( (simstep - _initSamplingVDF) % _writeFrequencyVDF != 0 )
        return;


    // calc global values
    this->CalcGlobalValuesVDF();  // calculate global velocity distribution sums


    // reset local values
    this->ResetLocalValuesVDF();

    // writing .rpf-files
    std::stringstream outputstreamVelo_py_abs;
    std::stringstream outputstreamVelo_py_pvx;
    std::stringstream outputstreamVelo_py_pvy;
    std::stringstream outputstreamVelo_py_pvz;
    std::stringstream outputstreamVelo_py_nvx;
    std::stringstream outputstreamVelo_py_nvy;
    std::stringstream outputstreamVelo_py_nvz;

    std::stringstream outputstreamVelo_ny_abs;
    std::stringstream outputstreamVelo_ny_pvx;
    std::stringstream outputstreamVelo_ny_pvy;
    std::stringstream outputstreamVelo_ny_pvz;
    std::stringstream outputstreamVelo_ny_nvx;
    std::stringstream outputstreamVelo_ny_nvy;
    std::stringstream outputstreamVelo_ny_nvz;

    std::stringstream filenamestreamVelo_py_abs;
    std::stringstream filenamestreamVelo_py_pvx;
    std::stringstream filenamestreamVelo_py_pvy;
    std::stringstream filenamestreamVelo_py_pvz;
    std::stringstream filenamestreamVelo_py_nvx;
    std::stringstream filenamestreamVelo_py_nvy;
    std::stringstream filenamestreamVelo_py_nvz;

    std::stringstream filenamestreamVelo_ny_abs;
    std::stringstream filenamestreamVelo_ny_pvx;
    std::stringstream filenamestreamVelo_ny_pvy;
    std::stringstream filenamestreamVelo_ny_pvz;
    std::stringstream filenamestreamVelo_ny_nvx;
    std::stringstream filenamestreamVelo_ny_nvy;
    std::stringstream filenamestreamVelo_ny_nvz;

    filenamestreamVelo_py_abs << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfabs";
    filenamestreamVelo_py_pvx << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfpvx";
    filenamestreamVelo_py_pvy << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfpvy";
    filenamestreamVelo_py_pvz << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfpvz";
    filenamestreamVelo_py_nvx << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfnvx";
    filenamestreamVelo_py_nvy << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfnvy";
    filenamestreamVelo_py_nvz << "VDF_region" << this->GetID() << "_TS" << simstep << "_py.vdfnvz";

    filenamestreamVelo_ny_abs << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfabs";
    filenamestreamVelo_ny_pvx << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfpvx";
    filenamestreamVelo_ny_pvy << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfpvy";
    filenamestreamVelo_ny_pvz << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfpvz";
    filenamestreamVelo_ny_nvx << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfnvx";
    filenamestreamVelo_ny_nvy << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfnvy";
    filenamestreamVelo_ny_nvz << "VDF_region" << this->GetID() << "_TS" << simstep << "_ny.vdfnvz";


    #ifdef ENABLE_MPI
        int rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank== 0)
        {
    #endif

            // header
            outputstreamVelo_py_abs << "v/y         ";
            outputstreamVelo_py_pvx << "v/y         ";
            outputstreamVelo_py_pvy << "v/y         ";
            outputstreamVelo_py_pvz << "v/y         ";
            outputstreamVelo_py_nvx << "v/y         ";
            outputstreamVelo_py_nvy << "v/y         ";
            outputstreamVelo_py_nvz << "v/y         ";

            outputstreamVelo_ny_abs << "v/y         ";
            outputstreamVelo_ny_pvx << "v/y         ";
            outputstreamVelo_ny_pvy << "v/y         ";
            outputstreamVelo_ny_pvz << "v/y         ";
            outputstreamVelo_ny_nvx << "v/y         ";
            outputstreamVelo_ny_nvy << "v/y         ";
            outputstreamVelo_ny_nvz << "v/y         ";

            // first line - discrete radius values
            for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
            {
                outputstreamVelo_py_abs << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvx << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvy << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvz << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvx << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvy << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvz << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);

                outputstreamVelo_ny_abs << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvx << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvy << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvz << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvx << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvy << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvz << _dShellMidpointsVDF[s] << std::setw(12) << std::setprecision(6);
            }
            outputstreamVelo_py_abs << endl;
            outputstreamVelo_py_pvx << endl;
            outputstreamVelo_py_pvy << endl;
            outputstreamVelo_py_pvz << endl;
            outputstreamVelo_py_nvx << endl;
            outputstreamVelo_py_nvy << endl;
            outputstreamVelo_py_nvz << endl;

            outputstreamVelo_ny_abs << endl;
            outputstreamVelo_ny_pvx << endl;
            outputstreamVelo_ny_pvy << endl;
            outputstreamVelo_ny_pvz << endl;
            outputstreamVelo_ny_nvx << endl;
            outputstreamVelo_ny_nvy << endl;
            outputstreamVelo_ny_nvz << endl;

            // global_log->info() << "header done." << endl;

            // velocity distribution matrix
            for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
            {
                outputstreamVelo_py_abs << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvx << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvy << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_pvz << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvx << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvy << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_py_nvz << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);

                outputstreamVelo_ny_abs << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvx << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvy << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_pvz << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvx << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvy << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);
                outputstreamVelo_ny_nvz << _dDiscreteVelocityValues[v] << std::setw(12) << std::setprecision(6);

                for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
                {
                    outputstreamVelo_py_abs << _veloDistrMatrixGlobal_py_abs[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_pvx << _veloDistrMatrixGlobal_py_pvx[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_pvy << _veloDistrMatrixGlobal_py_pvy[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_pvz << _veloDistrMatrixGlobal_py_pvz[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_nvx << _veloDistrMatrixGlobal_py_nvx[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_nvy << _veloDistrMatrixGlobal_py_nvy[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_py_nvz << _veloDistrMatrixGlobal_py_nvz[s][v] << std::setw(12) << std::setprecision(6);

                    outputstreamVelo_ny_abs << _veloDistrMatrixGlobal_ny_abs[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_pvx << _veloDistrMatrixGlobal_ny_pvx[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_pvy << _veloDistrMatrixGlobal_ny_pvy[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_pvz << _veloDistrMatrixGlobal_ny_pvz[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_nvx << _veloDistrMatrixGlobal_ny_nvx[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_nvy << _veloDistrMatrixGlobal_ny_nvy[s][v] << std::setw(12) << std::setprecision(6);
                    outputstreamVelo_ny_nvz << _veloDistrMatrixGlobal_ny_nvz[s][v] << std::setw(12) << std::setprecision(6);
                }

                outputstreamVelo_py_abs << endl;
                outputstreamVelo_py_pvx << endl;
                outputstreamVelo_py_pvy << endl;
                outputstreamVelo_py_pvz << endl;
                outputstreamVelo_py_nvx << endl;
                outputstreamVelo_py_nvy << endl;
                outputstreamVelo_py_nvz << endl;

                outputstreamVelo_ny_abs << endl;
                outputstreamVelo_ny_pvx << endl;
                outputstreamVelo_ny_pvy << endl;
                outputstreamVelo_ny_pvz << endl;
                outputstreamVelo_ny_nvx << endl;
                outputstreamVelo_ny_nvy << endl;
                outputstreamVelo_ny_nvz << endl;
            }

            // Datei zum schreiben öffnen, daten schreiben
            ofstream fileoutVelo_py_abs(filenamestreamVelo_py_abs.str().c_str(), ios::out);
            fileoutVelo_py_abs << outputstreamVelo_py_abs.str();
            fileoutVelo_py_abs.close();

            ofstream fileoutVelo_py_pvx(filenamestreamVelo_py_pvx.str().c_str(), ios::out);
            fileoutVelo_py_pvx << outputstreamVelo_py_pvx.str();
            fileoutVelo_py_pvx.close();

            ofstream fileoutVelo_py_pvy(filenamestreamVelo_py_pvy.str().c_str(), ios::out);
            fileoutVelo_py_pvy << outputstreamVelo_py_pvy.str();
            fileoutVelo_py_pvy.close();

            ofstream fileoutVelo_py_pvz(filenamestreamVelo_py_pvz.str().c_str(), ios::out);
            fileoutVelo_py_pvz << outputstreamVelo_py_pvz.str();
            fileoutVelo_py_pvz.close();

            ofstream fileoutVelo_py_nvx(filenamestreamVelo_py_nvx.str().c_str(), ios::out);
            fileoutVelo_py_nvx << outputstreamVelo_py_nvx.str();
            fileoutVelo_py_nvx.close();

            ofstream fileoutVelo_py_nvy(filenamestreamVelo_py_nvy.str().c_str(), ios::out);
            fileoutVelo_py_nvy << outputstreamVelo_py_nvy.str();
            fileoutVelo_py_nvy.close();

            ofstream fileoutVelo_py_nvz(filenamestreamVelo_py_nvz.str().c_str(), ios::out);
            fileoutVelo_py_nvz << outputstreamVelo_py_nvz.str();
            fileoutVelo_py_nvz.close();


            ofstream fileoutVelo_ny_abs(filenamestreamVelo_ny_abs.str().c_str(), ios::out);
            fileoutVelo_ny_abs << outputstreamVelo_ny_abs.str();
            fileoutVelo_ny_abs.close();

            ofstream fileoutVelo_ny_pvx(filenamestreamVelo_ny_pvx.str().c_str(), ios::out);
            fileoutVelo_ny_pvx << outputstreamVelo_ny_pvx.str();
            fileoutVelo_ny_pvx.close();

            ofstream fileoutVelo_ny_pvy(filenamestreamVelo_ny_pvy.str().c_str(), ios::out);
            fileoutVelo_ny_pvy << outputstreamVelo_ny_pvy.str();
            fileoutVelo_ny_pvy.close();

            ofstream fileoutVelo_ny_pvz(filenamestreamVelo_ny_pvz.str().c_str(), ios::out);
            fileoutVelo_ny_pvz << outputstreamVelo_ny_pvz.str();
            fileoutVelo_ny_pvz.close();

            ofstream fileoutVelo_ny_nvx(filenamestreamVelo_ny_nvx.str().c_str(), ios::out);
            fileoutVelo_ny_nvx << outputstreamVelo_ny_nvx.str();
            fileoutVelo_ny_nvx.close();

            ofstream fileoutVelo_ny_nvy(filenamestreamVelo_ny_nvy.str().c_str(), ios::out);
            fileoutVelo_ny_nvy << outputstreamVelo_ny_nvy.str();
            fileoutVelo_ny_nvy.close();

            ofstream fileoutVelo_ny_nvz(filenamestreamVelo_ny_nvz.str().c_str(), ios::out);
            fileoutVelo_ny_nvz << outputstreamVelo_ny_nvz.str();
            fileoutVelo_ny_nvz.close();

            // global_log->info() << "files closed." << endl;

    #ifdef ENABLE_MPI
        }
    #endif
}



// private methods

void SampleRegion::ResetLocalValuesVDF()
{
    // reset values
    for(unsigned int s = 0; s < _nNumShellsVDF; ++s)
    {
        // reset local velocity profile arrays
        for(unsigned int v = 0; v < _nNumDiscreteStepsVDF; ++v)
        {
            _veloDistrMatrixLocal_py_abs[s][v] = 0;
            _veloDistrMatrixLocal_py_pvx[s][v] = 0;
            _veloDistrMatrixLocal_py_pvy[s][v] = 0;
            _veloDistrMatrixLocal_py_pvz[s][v] = 0;
            _veloDistrMatrixLocal_py_nvx[s][v] = 0;
            _veloDistrMatrixLocal_py_nvy[s][v] = 0;
            _veloDistrMatrixLocal_py_nvz[s][v] = 0;

            _veloDistrMatrixLocal_ny_abs[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvx[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvy[s][v] = 0;
            _veloDistrMatrixLocal_ny_pvz[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvx[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvy[s][v] = 0;
            _veloDistrMatrixLocal_ny_nvz[s][v] = 0;
        }
    }
}


void SampleRegion::ResetLocalValuesProfiles()
{

    // reset cumulative data structures TODO: <-- should be placed elsewhere
    for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
    {
        // reset local num molecules
        _nNumMoleculesSumCumulativeLocal[s] = 0;

        // j+, j-
        _nNumMoleculesPlusSumCumulativeLocal[s] = 0;
        _nNumMoleculesMinusSumCumulativeLocal[s] = 0;

        // reset velocity sums
        for(unsigned int d = 0; d < 3; ++d)
        {
            // [position][x,y,z-component]
            _dVelocityComponentSumsCumulativeLocal[s][d] = 0.;
            _dSquaredVelocityComponentSumsCumulativeLocal[s][d] = 0.;

            // j+, j-
            _dVelocityComponentPlusSumsCumulativeLocal[s][d] = 0.;
            _dVelocityComponentMinusSumsCumulativeLocal[s][d] = 0.;
        }
    }


    // componentwise temperature / density
    for(unsigned int c = 0; c < _parent->GetNumComponents(); ++c)
    {
        for(unsigned int s = 0; s < _nNumShellsProfiles; ++s)
        {
            _nNumMoleculesCompLocal[c][s] = 0;
            _nRotDOFCompLocal[c][s] = 0;
            _d2EkinTransCompLocal[c][s] = 0.;
            _d2EkinRotCompLocal[c][s] = 0.;
        }
    }

}




// class RegionSampling

RegionSampling::RegionSampling(Domain* domain)
{
    //store domain pointer
    _domain = domain;

    // number of components
    _nNumComponents = domain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

}


RegionSampling::~RegionSampling()
{
}

void RegionSampling::AddRegion(double dLowerCorner[3], double dUpperCorner[3] )
{
    unsigned short nID = this->GetNumRegions() + 1;

    _vecSampleRegions.push_back(SampleRegion(this, dLowerCorner, dUpperCorner, nID) );
}

void RegionSampling::Init(Domain* domain)
{
    // init data structures
    std::vector<SampleRegion>::iterator it;

    for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
    {
        (*it).InitSamplingProfiles(RS_DIMENSION_Y, domain);
        (*it).InitSamplingVDF(RS_DIMENSION_Y);
    }
}

void RegionSampling::DoSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
    // sample profiles and vdf
    std::vector<SampleRegion>::iterator it;

    for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
    {
        (*it).SampleProfiles(mol, RS_DIMENSION_Y);
        (*it).SampleVDF(mol, RS_DIMENSION_Y);
    }
}

void RegionSampling::WriteData(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain)
{
    // write out profiles and vdf
    std::vector<SampleRegion>::iterator it;

    for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
    {
        (*it).WriteDataProfiles(domainDecomp, simstep, domain);
        (*it).WriteDataVDF(domainDecomp, simstep);
    }
}














