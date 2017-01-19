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
#include "utils/ObserverBase.h"
#include "utils/Region.h"

enum TemperatureControlTypes
{
	TCT_UNKNOWN = 0,
	TCT_CONSTANT_TEMPERATURE = 1,
	TCT_TEMPERATURE_ADJUST = 2,
	TCT_TEMPERATURE_GRADIENT = 3,
	TCT_TEMPERATURE_GRADIENT_LOWER = 4,
	TCT_TEMPERATURE_GRADIENT_RAISE = 5
};

class DomainDecompBase;
class ParticleContainer;
class AccumulatorBase;
class TemperatureControl;

namespace tec
{

class ControlRegion : public CuboidRegionObs
{
public:
    ControlRegion( ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nComp,
                   double* dTargetTemperature, double dTemperatureExponent, std::string strTransDirections,
                   int nTemperatureControlType, unsigned long nStartAdjust, unsigned long nStopAdjust, unsigned long nAdjustFreq );
    virtual ~ControlRegion();

    void SetSubdivision(unsigned int nNumSlabs) {_nNumSlabs = nNumSlabs;}
    void SetSubdivision(double dSlabWidth) {_dSlabWidthInit = dSlabWidth;}
    void PrepareSubdivision();  // need to be called before InitDataStructures()
    void PrepareDataStructures();
    void CalcGlobalValues(unsigned long simstep);
    void MeasureKineticEnergy(Molecule* mol);
    void ControlTemperature(Molecule* mol);
    void ResetLocalValues();
    void UpdateSlabParameters();

    // write out data --> heat supply
    void CalcGlobalValuesDeltaEkin();
    void ResetValuesDeltaEkin();
    void WriteHeaderDeltaEkin();
    void WriteDataDeltaEkin(unsigned long simstep);

    // set adjust parameters
    void SetTemperatureControlAdjustParameters(unsigned long nStartAdjust, unsigned long nStopAdjust, unsigned long nAdjustFreq)
    {
    	_nStartAdjust = nStartAdjust;
    	_nStopAdjust  = nStopAdjust;
    	_nAdjustFreq  = nAdjustFreq;
    }

	// private methods
private:
    void InitDataStructurePointers();
    void AdjustTemperatureGradient();
    void AllocateDataStructuresT();
	void InitDataStructuresT();
	void AllocateDataStructuresDEkin();
	void InitDataStructuresDEkin();

private:
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

    double _dTargetTemperature[2];
    double* _dTargetTemperatureVec;
    double _dTemperatureExponent;
    unsigned int _nTargetComponentID;
    unsigned short _nNumThermostatedTransDirections;

    double _dTargetTemperatureActual;
    double _dDeltaTemperatureAdjust;
    unsigned long _nStartAdjust;
    unsigned long _nStopAdjust;
    unsigned long _nAdjustFreq;
    unsigned int _nTemperatureControlType;

    AccumulatorBase* _accumulator;

    // heat supply
    unsigned int _nNumSlabsDeltaEkin;
    unsigned int _nNumSlabsDEkinReserve;
    double _dSlabWidthDeltaEkin;

    unsigned long* _nNumMoleculesSumLocal;
    unsigned long* _nNumMoleculesSumGlobal;

    double* _dDelta2EkinTransSumLocal;
    double* _dDelta2EkinTransSumGlobal;

    // instances / ID
    static unsigned short _nStaticID;
};

}

class Domain;
class DomainDecompBase;
class TemperatureControl : public ControlInstance
{
public:
    TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~TemperatureControl();

    std::string GetShortName() {return "TeC";}
    void AddRegion(tec::ControlRegion* region);
    int GetNumRegions() {return _vecControlRegions.size();}
    tec::ControlRegion* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

    void PrepareRegionSubdivisions();
    void PrepareRegionDataStructures();
    void InitControl(unsigned long simstep);
    void MeasureKineticEnergy(Molecule* mol, unsigned long simstep);
    void CalcGlobalValues(unsigned long simstep);
    void ControlTemperature(Molecule* mol, unsigned long simstep);

    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    // loops over molecule container
    void DoLoopsOverMolecules(ParticleContainer* particleContainer, unsigned long simstep);

    // heat supply
    void SetDeltaEkinParameters( unsigned int nWriteFreqDeltaEkin, unsigned int nNumSlabsDeltaEkin)
    {
        _nWriteFreqDeltaEkin = nWriteFreqDeltaEkin; _nNumSlabsDeltaEkin = nNumSlabsDeltaEkin; _bWriteDataDeltaEkin = true;
    }
    void WriteDataDeltaEkin(unsigned long simstep);

private:
    std::vector<tec::ControlRegion*> _vecControlRegions;
    unsigned long _nControlFreq;
    unsigned long _nStart;
    unsigned long _nStop;

    // heat supply
    bool _bWriteDataDeltaEkin;
    unsigned int _nWriteFreqDeltaEkin;
    unsigned int _nNumSlabsDeltaEkin;
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
