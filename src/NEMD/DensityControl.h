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
#include "utils/ObserverBase.h"
#include "utils/Region.h"

class Simulation;
class Domain;
class DomainDecompBase;
class DensityControl;

namespace dec
{

class ControlRegion : public CuboidRegionObs
{
public:
    ControlRegion(ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity);
    ~ControlRegion();

    void CheckBounds();
    void Init();
    void InitMPI();
    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    double GetDensityGlobal() {return _dDensityGlobal;}

    void CalcGlobalValues();
    void UpdateGlobalDensity(bool bDeleteMolecule);
    void MeasureDensity(Molecule* mol);

    void ControlDensity(Molecule* mol, Simulation* simulation, bool& bDeleteMolecule);

    void ResetLocalValues();
    void UpdateVolume() {_dVolume = this->GetWidth(0) * this->GetWidth(1) * this->GetWidth(2); _dInvertVolume = 1./_dVolume;}

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(unsigned long simstep);

private:
    double _dVolume;
    double _dInvertVolume;
    unsigned long _nNumMoleculesLocal;
    unsigned long _nNumMoleculesGlobal;

    double _dTargetDensity;
    double _dDensityGlobal;
    unsigned int _nTargetComponentID;  // inert gas

    int* _ranks;
    bool _bProcessIsRelevant;
    MPI_Comm _newcomm;

    // deleted molecules data
    unsigned long _nDeletedNumMoleculesLocal;
    unsigned long _nDeletedNumMoleculesGlobal;
    double _dDeletedEkinLocal[3];
    double _dDeletedEkinGlobal[3];
    double _dDeletedVelocityLocal[3];
    double _dDeletedVelocityGlobal[3];
    double _dDeletedVelocitySquaredLocal[3];
    double _dDeletedVelocitySquaredGlobal[3];

    // instances / ID
    static unsigned short _nStaticID;
};

}


class DensityControl : public ControlInstance
{
public:
    DensityControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~DensityControl();

    std::string GetShortName() {return "DeC";}
    void AddRegion(dec::ControlRegion* region);
    int GetNumRegions() {return _vecControlRegions.size();}
    dec::ControlRegion* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

    void CheckRegionBounds();
    void Init(unsigned long simstep);
    void InitMPI();
    void MeasureDensity(Molecule* mol, unsigned long simstep);
    void CalcGlobalValues(unsigned long simstep);
    void UpdateGlobalDensities(unsigned long simstep, bool bDeleteMolecule );
    void ControlDensity(Molecule* mol, Simulation* simulation, unsigned long simstep, bool& bDeleteMolecule);

    unsigned long GetControlFreq() {return _nControlFreq;}
    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(unsigned long simstep);

private:
    std::vector<dec::ControlRegion*> _vecControlRegions;
    unsigned long _nControlFreq;
    unsigned long _nStart;
    unsigned long _nStop;

    bool _bProcessIsRelevant;

    unsigned int _nWriteFreqDeleted;
};


#endif /* DENSITYCONTROL_H_ */
