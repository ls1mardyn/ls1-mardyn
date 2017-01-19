/*
 * DriftControl.h
 *
 *  Created on: 22.02.2015
 *      Author: mheinen
 */

#ifndef DRIFTCONTROL_H_
#define DRIFTCONTROL_H_

#include <vector>
#include <string>
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"


class ControlRegion
{
public:
    ControlRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, double dDirection[3], const double& dTargetVal);
    ~ControlRegion();

    double* GetLowerCorner() {return _dLowerCorner;}
    double* GetUpperCorner() {return _dUpperCorner;}
    void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal;}
    void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal;}
    double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
    void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
    unsigned long GetDriftVelocityGlobal(unsigned short nDim) {return _dDriftVelocityGlobal[nDim];}
    void CalcGlobalValues(DomainDecompBase* domainDecomp);
    void MeasureDrift(Molecule* mol, DomainDecompBase* domainDecomp);

    // After drift is measured by calling MeasureDrift() the scale factor has to be calculated
    // by calling CalcScaleFactor() to be able to rescale kin. energy of control region to its origin value
    void CalcScaleFactor();  // needed to keep kin. energy in control region

    void ControlDrift(Molecule* mol);

    void ResetLocalValues();

private:
    double _dLowerCorner[3];
    double _dUpperCorner[3];
    double _dDriftVelocityLocal[3];
    double _dDriftVelocityGlobal[3];
    double _dEkinLocal[3];
    double _dEkinGlobal[3];
    double _dScaleFactor;
    unsigned long _nNumMoleculesLocal;
    unsigned long _nNumMoleculesGlobal;
    unsigned short _nRegionID;

    double _dDriftDirection[3];
    double _dDriftVelocityTargetVal;
    double _dDriftVelocityTargetVector[3];
    double _dAddVector[3];

    unsigned int _nTargetComponentID;
};


class DriftControl
{
public:
    DriftControl(unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~DriftControl();

    void AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, double dDirection[3], const double& dTargetVal);
    int GetNumRegions() {return _vecControlRegions.size();}
    ControlRegion* GetControlRegion(unsigned short nRegionID) {return &(_vecControlRegions.at(nRegionID-1) ); }  // vector index starts with 0, region index with 1

    void Init(unsigned long simstep);
    void MeasureDrift(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
    void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
    void CalcScaleFactors(unsigned long simstep);
    void ControlDrift(Molecule* mol, unsigned long simstep);

    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

private:
    std::vector<ControlRegion> _vecControlRegions;
    unsigned long _nControlFreq;
    unsigned long _nStart;
    unsigned long _nStop;
};


#endif /* DRIFTCONTROL_H_ */
