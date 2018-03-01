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
#include "utils/ObserverBase.h"
#include "utils/Region.h"

namespace drc
{

class ControlRegion : public CuboidRegionObs
{
public:
	ControlRegion(ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3] );
	ControlRegion(ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3],
			unsigned int nTargetComponentID, const double dDirection[3], const double& dTargetVal);
	virtual ~ControlRegion();

	void readXML(XMLfileUnits& xmlconfig);
    unsigned long GetDriftVelocityGlobal(unsigned short nDim) {return _dDriftVelocityGlobal[nDim];}
    void CalcGlobalValues(DomainDecompBase* domainDecomp);
    void MeasureDrift(Molecule* mol);

    // After drift is measured by calling MeasureDrift() the scale factor has to be calculated
    // by calling CalcScaleFactor() to be able to rescale kin. energy of control region to its origin value
    void CalcScaleFactor();  // needed to keep kin. energy in control region

    void ControlDrift(Molecule* mol);
    void ResetLocalValues();

private:
	void PrepareDriftVector();
    void PrepareDriftVector(const double dDirection[3], const double& dTargetVal);

private:
    double _dDriftVelocityLocal[3];
    double _dDriftVelocityGlobal[3];
    double _dEkinLocal[3];
    double _dEkinGlobal[3];
    double _dScaleFactor;
    unsigned long _nNumMoleculesLocal;
    unsigned long _nNumMoleculesGlobal;

    double _dDriftDirection[3];
    double _dDriftVelocityTargetVal;
    double _dDriftVelocityTargetVector[3];
    double _dAddVector[3];

    unsigned int _nTargetComponentID;

    // instances / ID
    static unsigned short _nStaticID;
};

}

class XMLfileUnits;
class DriftControl : public ControlInstance
{
public:
	DriftControl(Domain* domain, DomainDecompBase* domainDecomp);
    DriftControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
	virtual ~DriftControl();

	void readXML(XMLfileUnits& xmlconfig);
    std::string GetShortName() {return "DrC";}
    void AddRegion(drc::ControlRegion* region);
    int GetNumRegions() {return _vecControlRegions.size();}
    drc::ControlRegion* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

    void Init(unsigned long simstep);
    void MeasureDrift(Molecule* mol, unsigned long simstep);
    void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
    void CalcScaleFactors(unsigned long simstep);
    void ControlDrift(Molecule* mol, unsigned long simstep);

    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

private:
    std::vector<drc::ControlRegion*> _vecControlRegions;
	unsigned long _nStart;
	unsigned long _nControlFreq;
	unsigned long _nStop;
};


#endif /* DRIFTCONTROL_H_ */
