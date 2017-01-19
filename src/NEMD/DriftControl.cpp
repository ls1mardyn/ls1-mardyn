/*
 * MEX.cpp
 *
 *  Created on: 16.02.2015
 *      Author: mheinen
 */

#include "NEMD/DriftControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
//#include <cmath>
//#include <list>
//#include <map>
#include "Domain.h"
//#include "NEMD/ParticleInsertion.h"
//#include "utils/Random.h"
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
//#include <iterator>  // std::advance

using namespace std;


// class ControlRegion

ControlRegion::ControlRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, double dDirection[3], const double& dTargetVal)
{
    // region span
    for(unsigned short d=0; d<3; ++d)
    {
        _dLowerCorner[d] = dLowerCorner[d];
        _dUpperCorner[d] = dUpperCorner[d];
    }

    // normalize direction vector
    double dLength = 0.;

    for(unsigned short d=0; d<3; ++d)
    {
       dLength += dDirection[d] * dDirection[d];
    }
    dLength = sqrt(dLength);

    // drift velocity target direction, target vector
    for(unsigned short d=0; d<3; ++d)
    {
        _dDriftDirection[d] = dDirection[d] / dLength;
        _dDriftVelocityTargetVector[d] = _dDriftDirection[d] * dTargetVal;
    }

//    cout << "_dDriftDirection = " << _dDriftDirection[0] << ", " << _dDriftDirection[1] << ", " << _dDriftDirection[2] << endl;
//    cout << "_dDriftVelocityTargetVector = " << _dDriftVelocityTargetVector[0] << ", " << _dDriftVelocityTargetVector[1] << ", " << _dDriftVelocityTargetVector[2] << endl;

    // drift velocity target value
    _dDriftVelocityTargetVal = dTargetVal;

    // vector to add to each molecule to obtain target drift velocity of system inside control region
    for(unsigned short d=0; d<3; ++d)
    {
        _dAddVector[d] = 0.;
    }


    // target component ID
    _nTargetComponentID = nTargetComponentID;
}


ControlRegion::~ControlRegion()
{

}


void ControlRegion::CalcGlobalValues(DomainDecompBase* domainDecomp)
{
#ifdef ENABLE_MPI

    MPI_Allreduce( _dDriftVelocityLocal, _dDriftVelocityGlobal, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _dEkinLocal, _dEkinGlobal, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

    for(unsigned short d=0; d<3; ++d)
    {
        _dDriftVelocityGlobal[d] = _dDriftVelocityLocal[d];
    }
#endif

//    cout << "_dDriftVelocityLocal[0] = " << _dDriftVelocityLocal[0] << endl;
//    cout << "_dDriftVelocityLocal[1] = " << _dDriftVelocityLocal[1] << endl;
//    cout << "_dDriftVelocityLocal[2] = " << _dDriftVelocityLocal[2] << endl;
//
//    cout << "_dDriftVelocityGlobal[0] = " << _dDriftVelocityGlobal[0] << endl;
//    cout << "_dDriftVelocityGlobal[1] = " << _dDriftVelocityGlobal[1] << endl;
//    cout << "_dDriftVelocityGlobal[2] = " << _dDriftVelocityGlobal[2] << endl;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendUnsLong(_nNumMoleculesLocal);
    domainDecomp->collCommAllreduceSum();
    _nNumMoleculesGlobal = domainDecomp->collCommGetUnsLong();
    domainDecomp->collCommFinalize();

//    cout << "_nNumMoleculesLocal = " << _nNumMoleculesLocal << endl;
//    cout << "_nNumMoleculesGlobal = " << _nNumMoleculesGlobal << endl;

    // velocity sum has to be divided by number of molecules taken into account
    if(_nNumMoleculesGlobal > 0)
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dDriftVelocityGlobal[d] = _dDriftVelocityGlobal[d] / (double) (_nNumMoleculesGlobal);
        }
    }

    // calculate vector to add to each molecule to obtain target drift velocity of system inside control region
    if(0. == _dDriftVelocityTargetVal)
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dAddVector[d] = _dDriftVelocityGlobal[d] * -1.;
        }
    }
    else
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dAddVector[d] = _dDriftVelocityTargetVector[d] - _dDriftVelocityGlobal[d];
        }
    }

//    cout << "_dDriftVelocityGlobal = " << _dDriftVelocityGlobal[0] << ", " << _dDriftVelocityGlobal[1] << ", " << _dDriftVelocityGlobal[2] << endl;
//    cout << "_dDriftVelocityTargetVector = " << _dDriftVelocityTargetVector[0] << ", " << _dDriftVelocityTargetVector[1] << ", " << _dDriftVelocityTargetVector[2] << endl;
//    cout << "_dAddVector = " << _dAddVector[0] << ", " << _dAddVector[1] << ", " << _dAddVector[2] << endl;

}

void ControlRegion::MeasureDrift(Molecule* mol, DomainDecompBase* domainDecomp)
{
    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // sum up drift and kin. energy
    for(unsigned short d = 0; d<3; ++d)
    {
        double v = mol->v(d);

        _dDriftVelocityLocal[d] += v;
        _dEkinLocal[d] += v*v;  // * mol->mass();
    }

    // count num molecules
    _nNumMoleculesLocal++;
}

void ControlRegion::CalcScaleFactor()
{
    double dEkin_new = 0;
    double dEkin_old = 0;
    double dDeltaEkin = 0;

    // sum up x, y, z dimension of kin. energy values
    for(unsigned short d = 0; d<3; ++d)
    {
        dEkin_old += _dEkinGlobal[d];

        // (v^2_target - v^2_drift) * N * m yields the kinetic energy added or removed from system
        // dDeltaEkin is positiv, when energy added
        // dDeltaEkin is negativ, when energy removed
        dDeltaEkin += _dDriftVelocityTargetVector[d] * _dDriftVelocityTargetVector[d] - _dDriftVelocityGlobal[d] * _dDriftVelocityGlobal[d];
    }

    dDeltaEkin = dDeltaEkin * _nNumMoleculesGlobal;
    dEkin_new = dEkin_old +dDeltaEkin;

//    cout << "dDeltaEkin = " << dDeltaEkin << endl;
//    cout << "dEkin_new = " << dEkin_new << endl;
//    cout << "dEkin_old = " << dEkin_old << endl;

    double f = 1.;

    if(dEkin_new > 0.)
        f = dEkin_old / dEkin_new;

    _dScaleFactor = sqrt(f);

    // no scaling
    _dScaleFactor = 1.;
}

void ControlRegion::ControlDrift(Molecule* mol)
{
    if(_nNumMoleculesGlobal < 1)
        return;

    // check componentID
        if(mol->componentid()+1 != _nTargetComponentID)  // program intern componentID starts with 0
            return;

    // check if molecule is inside
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // change velocity vector
    for(unsigned short d = 0; d<3; ++d)
    {
        double v = mol->v(d);

//        cout << "v[" << d << "] = " << v << endl;
//
//        cout << "_dAddVector[" << d << "] = " << _dAddVector[d] << endl;

        // control drift
        mol->setv( d, (v + _dAddVector[d]) * _dScaleFactor);

//        cout << "v[" << d << "] = " << mol->v(d) << endl;
    }
}

void ControlRegion::ResetLocalValues()
{
    // reset local values
    for(unsigned short d = 0; d<3; ++d)
    {
        _dDriftVelocityLocal[d] = 0.;
        _dEkinLocal[d] = 0.;
    }

    _nNumMoleculesLocal = 0;

    // scale factor
    _dScaleFactor = 1.;
}



// class DriftControl

DriftControl::DriftControl(unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;
}

DriftControl::~DriftControl()
{

}

void DriftControl::AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, double dDirection[3], const double& dTargetVal)
{
    _vecControlRegions.push_back(ControlRegion(dLowerCorner, dUpperCorner, nTargetComponentID, dDirection, dTargetVal) );
}

void DriftControl::MeasureDrift(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<ControlRegion>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).MeasureDrift(mol, domainDecomp);
    }
}

void DriftControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<ControlRegion>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).CalcGlobalValues(domainDecomp);
    }
}

void DriftControl::CalcScaleFactors(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc scale factors for control regions
    std::vector<ControlRegion>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).CalcScaleFactor();
    }
}

void DriftControl::ControlDrift(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<ControlRegion>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ControlDrift(mol);
    }
}

void DriftControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<ControlRegion>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ResetLocalValues();
    }
}




