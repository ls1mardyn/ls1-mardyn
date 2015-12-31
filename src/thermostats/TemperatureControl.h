/*
 * TemperatureControl.h
 *
 *  Created on: 27.05.2015
 *      Author: mheinen
 */

#ifndef TEMPERATURECONTROL_H_
#define TEMPERATURECONTROL_H_

#include <vector>
#include <string>

#include "molecules/Molecule.h"

class DomainDecompBase;
class ParticleContainer;
class AccumulatorBase;
class TemperatureControl;
class ControlRegionT
{
public:
    ControlRegionT( TemperatureControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp,
                    double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections, unsigned short nRegionID, unsigned int nNumSlabsDeltaEkin );
    ~ControlRegionT();

    void Init();
    unsigned short GetID() {return _nRegionID;}
    double* GetLowerCorner() {return _dLowerCorner;}
    double* GetUpperCorner() {return _dUpperCorner;}
    void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal; this->UpdateSlabParameters();}
    void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal; this->UpdateSlabParameters();}
    double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
    void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
    void CalcGlobalValues(DomainDecompBase* domainDecomp);
    void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp);

    void ControlTemperature(Molecule* mol);

    void ResetLocalValues();

    void UpdateSlabParameters();

    // write out data --> heat supply
    void CalcGlobalValuesDeltaEkin();
    void ResetValuesDeltaEkin();
    void WriteHeaderDeltaEkin(DomainDecompBase* domainDecomp, Domain* domain);
    void WriteDataDeltaEkin(DomainDecompBase* domainDecomp, unsigned long simstep);

private:
    TemperatureControl* _parent;

    double _dLowerCorner[3];
    double _dUpperCorner[3];

    unsigned int _nNumSlabs;
    unsigned int _nNumSlabsReserve;
    double _dSlabWidth;
    double _dSlabWidthInit;

    unsigned long* _nNumMoleculesLocal;
    unsigned long* _nNumMoleculesGlobal;
    unsigned long* _nRotDOFLocal;
    unsigned long* _nRotDOFGlobal;

    double* _d2EkinTransLocal;
    double* _d2EkinTransGlobal;
    double* _d2EkinRotLocal;
    double* _d2EkinRotGlobal;

    double* _dBetaTransGlobal;
    double* _dBetaRotGlobal;

    double _dTargetTemperature;
    double _dTemperatureExponent;
    unsigned int _nTargetComponentID;
    unsigned short _nNumThermostatedTransDirections;

    unsigned short _nRegionID;

    AccumulatorBase* _accumulator;

    // heat supply
    unsigned int _nNumSlabsDeltaEkin;
    double _dSlabWidthDeltaEkin;

    unsigned long* _nNumMoleculesSumLocal;
    unsigned long* _nNumMoleculesSumGlobal;

    double* _dDelta2EkinTransSumLocal;
    double* _dDelta2EkinTransSumGlobal;

};


class Domain;
class DomainDecompBase;
class TemperatureControl
{
public:
    TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~TemperatureControl();

    void AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp, double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections);
    int GetNumRegions() {return _vecControlRegions.size();}
    ControlRegionT* GetControlRegion(unsigned short nRegionID) {return &(_vecControlRegions.at(nRegionID-1) ); }  // vector index starts with 0, region index with 1

    void Init(unsigned long simstep);
    void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
    void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
    void ControlTemperature(Molecule* mol, unsigned long simstep);

    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    // loops over molecule container
    void DoLoopsOverMolecules(DomainDecompBase*, ParticleContainer* particleContainer, unsigned long simstep);

    Domain* GetDomain() {return _domain;}
    DomainDecompBase* GetDomainDecomposition() {return _domainDecomp;}

    // heat supply
    void SetDeltaEkinParameters( unsigned int nWriteFreqDeltaEkin, unsigned int nNumSlabsDeltaEkin)
    {
        _nWriteFreqDeltaEkin = nWriteFreqDeltaEkin; _nNumSlabsDeltaEkin = nNumSlabsDeltaEkin;
    }
    void WriteDataDeltaEkin(DomainDecompBase* domainDecomp, unsigned long simstep);

private:
    std::vector<ControlRegionT> _vecControlRegions;
    unsigned long _nControlFreq;
    unsigned long _nStart;
    unsigned long _nStop;

    // heat supply
    unsigned int _nWriteFreqDeltaEkin;
    unsigned int _nNumSlabsDeltaEkin;

    Domain* _domain;
    DomainDecompBase* _domainDecomp;
};


// Accumulate kinetic energy dependent on which translatoric directions should be thermostated

class AccumulatorBase
{
protected:
    AccumulatorBase() {}
    virtual ~AccumulatorBase() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol) = 0;
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum) = 0;
};

class AccumulatorX : public AccumulatorBase
{
public:
    AccumulatorX() {}
    virtual ~AccumulatorX() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double m  = mol->mass();

        return m * vx*vx;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(0, mol->v(0) * vcorr);
    }
};

class AccumulatorY : public AccumulatorBase
{
public:
    AccumulatorY() {}
    virtual ~AccumulatorY() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vy = mol->v(1);
        double m  = mol->mass();

        return m * vy*vy;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(1, mol->v(1) * vcorr);
    }
};

class AccumulatorZ : public AccumulatorBase
{
public:
    AccumulatorZ() {}
    virtual ~AccumulatorZ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * vz*vz;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(2, mol->v(2) * vcorr);
    }
};

class AccumulatorXY : public AccumulatorBase
{
public:
    AccumulatorXY() {}
    virtual ~AccumulatorXY() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vy = mol->v(1);
        double m  = mol->mass();

        return m * (vx*vx + vy*vy);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(0, mol->v(0) * vcorr);
        mol->setv(1, mol->v(1) * vcorr);
    }
};

class AccumulatorXZ : public AccumulatorBase
{
public:
    AccumulatorXZ() {}
    virtual ~AccumulatorXZ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vx*vx + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(0, mol->v(0) * vcorr);
        mol->setv(2, mol->v(2) * vcorr);
    }
};

class AccumulatorYZ : public AccumulatorBase
{
public:
    AccumulatorYZ() {}
    virtual ~AccumulatorYZ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vy = mol->v(1);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vy*vy + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(1, mol->v(1) * vcorr);
        mol->setv(2, mol->v(2) * vcorr);
    }
};

class AccumulatorXYZ : public AccumulatorBase
{
public:
    AccumulatorXYZ() {}
    virtual ~AccumulatorXYZ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vy = mol->v(1);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vx*vx + vy*vy + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr)
    {
        mol->setv(0, mol->v(0) * vcorr);
        mol->setv(1, mol->v(1) * vcorr);
        mol->setv(2, mol->v(2) * vcorr);
    }
};


// Accumulators with calculation of heat supply

class AccumulatorXQ : public AccumulatorBase
{
public:
    AccumulatorXQ() {}
    virtual ~AccumulatorXQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double m  = mol->mass();

        return m * vx*vx;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[0] = mol->v(0);
        v_new[0] = v_old[0] * vcorr;

        // set new velocities
        mol->setv(0, v_new[0]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[0] * v_new[0]) - (v_old[0] * v_old[0]) );
    }
};

class AccumulatorYQ : public AccumulatorBase
{
public:
    AccumulatorYQ() {}
    virtual ~AccumulatorYQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vy = mol->v(1);
        double m  = mol->mass();

        return m * vy*vy;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[1] = mol->v(1);
        v_new[1] = v_old[1] * vcorr;

        // set new velocities
        mol->setv(1, v_new[1]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[1] * v_new[1]) - (v_old[1] * v_old[1]) );
    }
};

class AccumulatorZQ : public AccumulatorBase
{
public:
    AccumulatorZQ() {}
    virtual ~AccumulatorZQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * vz*vz;
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[2] = mol->v(2);
        v_new[2] = v_old[2] * vcorr;

        // set new velocities
        mol->setv(2, v_new[2]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[2] * v_new[2]) - (v_old[2] * v_old[2]) );
    }
};

class AccumulatorXYQ : public AccumulatorBase
{
public:
    AccumulatorXYQ() {}
    virtual ~AccumulatorXYQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vy = mol->v(1);
        double m  = mol->mass();

        return m * (vx*vx + vy*vy);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[1] = mol->v(1);
        v_new[0] = v_old[0] * vcorr;
        v_new[1] = v_old[1] * vcorr;

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(1, v_new[1]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[0] * v_new[0]) - (v_old[0] * v_old[0]) );
        dDelta2EkinTransSum += m * ( (v_new[1] * v_new[1]) - (v_old[1] * v_old[1]) );
    }
};

class AccumulatorXZQ : public AccumulatorBase
{
public:
    AccumulatorXZQ() {}
    virtual ~AccumulatorXZQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vx*vx + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[2] = mol->v(2);
        v_new[0] = v_old[0] * vcorr;
        v_new[2] = v_old[2] * vcorr;

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(2, v_new[2]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[0] * v_new[0]) - (v_old[0] * v_old[0]) );
        dDelta2EkinTransSum += m * ( (v_new[2] * v_new[2]) - (v_old[2] * v_old[2]) );
    }
};

class AccumulatorYZQ : public AccumulatorBase
{
public:
    AccumulatorYZQ() {}
    virtual ~AccumulatorYZQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vy = mol->v(1);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vy*vy + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[1] = mol->v(1);
        v_old[2] = mol->v(2);
        v_new[1] = v_old[1] * vcorr;
        v_new[2] = v_old[2] * vcorr;

        // set new velocities
        mol->setv(1, v_new[1]);
        mol->setv(2, v_new[2]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[1] * v_new[1]) - (v_old[1] * v_old[1]) );
        dDelta2EkinTransSum += m * ( (v_new[2] * v_new[2]) - (v_old[2] * v_old[2]) );
    }
};

class AccumulatorXYZQ : public AccumulatorBase
{
public:
    AccumulatorXYZQ() {}
    virtual ~AccumulatorXYZQ() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol)
    {
        double vx = mol->v(0);
        double vy = mol->v(1);
        double vz = mol->v(2);
        double m  = mol->mass();

        return m * (vx*vx + vy*vy + vz*vz);
    }
    virtual void ScaleVelocityComponents(Molecule* mol, double vcorr, double& dDelta2EkinTransSum)
    {
        double m  = mol->mass();
        double v_old[3];
        double v_new[3];

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[1] = mol->v(1);
        v_old[2] = mol->v(2);
        v_new[0] = v_old[0] * vcorr;
        v_new[1] = v_old[1] * vcorr;
        v_new[2] = v_old[2] * vcorr;

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(1, v_new[1]);
        mol->setv(2, v_new[2]);

        // accumulate added kinetic energy (heat)
        dDelta2EkinTransSum += m * ( (v_new[0] * v_new[0]) - (v_old[0] * v_old[0]) );
        dDelta2EkinTransSum += m * ( (v_new[1] * v_new[1]) - (v_old[1] * v_old[1]) );
        dDelta2EkinTransSum += m * ( (v_new[2] * v_new[2]) - (v_old[2] * v_old[2]) );
    }
};


#endif /* TEMPERATURECONTROL_H_ */
