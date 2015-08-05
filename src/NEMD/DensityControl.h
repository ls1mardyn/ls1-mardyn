/*
 * DensityControl.h
 *
 *  Created on: 29.05.2015
 *      Author: mheinen
 */

#ifndef DENSITYCONTROL_H_
#define DENSITYCONTROL_H_

#include <vector>
#include <string>
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"

class Simulation;
class Domain;
class DomainDecompBase;
class DensityControl;

class ControlRegionD
{
public:
    ControlRegionD(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity);
    ~ControlRegionD();

    void InitMPI();
    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    double* GetLowerCorner() {return _dLowerCorner;}
    double* GetUpperCorner() {return _dUpperCorner;}
    void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal;}
    void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal;}
    double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
    void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
    double GetDensityGlobal() {return _dDensityGlobal;}

    void CalcGlobalValues(DomainDecompBase* domainDecomp);
    void UpdateGlobalDensity(DomainDecompBase* domainDecomp, bool bDeleteMolecule);
    void MeasureDensity(Molecule* mol, DomainDecompBase* domainDecomp);

    void ControlDensity(DomainDecompBase* domainDecomp, Molecule* mol, Simulation* simulation, bool& bDeleteMolecule);

    void ResetLocalValues();
    void UpdateVolume() {_dInvertVolume = 1. / (this->GetWidth(0) * this->GetWidth(1) * this->GetWidth(2) );}

private:
    double _dLowerCorner[3];
    double _dUpperCorner[3];
    double _dInvertVolume;
    unsigned long _nNumMoleculesLocal;
    unsigned long _nNumMoleculesGlobal;
    unsigned short _nRegionID;

    double _dTargetDensity;
    double _dDensityGlobal;
    unsigned int _nTargetComponentID;  // inert gas

    DensityControl* _parent;

    int* _ranks;
    bool _bProcessIsRelevant;
    MPI_Comm _newcomm;
};


class DensityControl
{
public:
    DensityControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~DensityControl();

    void AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity);
    int GetNumRegions() {return _vecControlRegions.size();}
    ControlRegionD* GetControlRegion(unsigned short nRegionID) {return &(_vecControlRegions.at(nRegionID-1) ); }  // vector index starts with 0, region index with 1

    void Init(unsigned long simstep);
    void MeasureDensity(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
    void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
    void UpdateGlobalDensities(DomainDecompBase* domainDecomp, unsigned long simstep, bool bDeleteMolecule );
    void ControlDensity(DomainDecompBase* domainDecomp, Molecule* mol, Simulation* simulation, unsigned long simstep, bool& bDeleteMolecule);

    unsigned long GetControlFreq() {return _nControlFreq;}
    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    Domain* GetDomain() {return _domain;}
    DomainDecompBase* GetDomainDecomposition() {return _domainDecomp;}

    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

private:
    std::vector<ControlRegionD> _vecControlRegions;
    unsigned long _nControlFreq;
    unsigned long _nStart;
    unsigned long _nStop;

    Domain* _domain;
    DomainDecompBase* _domainDecomp;

    bool _bProcessIsRelevant;
};


#endif /* DENSITYCONTROL_H_ */
