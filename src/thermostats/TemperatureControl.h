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
#include <map>
#include <cstdint>

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

struct TimestepControl {
	uint64_t start, stop, freq;
};

template <typename T>
class TargetCtrl {
	public:
	T start, actual, stop;
};

struct AdjustCtrl {
	double delta;
	TimestepControl control;
};

struct TemperatureCtrl {
	TargetCtrl<double> temperature;
	double exponent;
	uint32_t cID;
	uint16_t numDirections, controlType;
	AdjustCtrl adjust;
	double* vec;
};

class XMLfileUnits;
class DomainDecompBase;
class ParticleContainer;
class AccumulatorBase;
class TemperatureControl;

namespace tec
{

class ControlRegion : public CuboidRegionObs
{
public:
    ControlRegion( TemperatureControl* const parent, double dLowerCorner[3], double dUpperCorner[3] );
    virtual ~ControlRegion();

	void readXML(XMLfileUnits& xmlconfig);
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

	virtual void Print(std::ostream& os)
	{
		os << "----------------------------------------------------------------" << std::endl;
		os << "ID: " << _nID << std::endl;
		os << "width: " << this->GetWidth(0) << " " << this->GetWidth(1) << " " << this->GetWidth(2) << std::endl;
		double lc[3];
		double uc[3];
		this->GetLowerCorner(lc);
		this->GetUpperCorner(uc);
		os << "lowerCorner: " << lc[0] << " " << lc[1] << " " << lc[2] << std::endl;
		os << "upperCorner: " << uc[0] << " " << uc[1] << " " << uc[2] << std::endl;
		os << "target T: " << _target.temperature.start << ", " << _target.temperature.stop << std::endl;
		os << "target cid: " << _target.cID << std::endl;
		os << "----------------------------------------------------------------" << std::endl;
	}

	// private methods
private:
    void InitDataStructurePointers();
    void AdjustTemperatureGradient();
    void AllocateDataStructuresT();
	void InitDataStructuresT();
	void AllocateDataStructuresDEkin();
	void InitDataStructuresDEkin();
	void PrepareControlType();
	void PrepareAccumulator(const std::string& strTransDirections);

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

    // heat supply
    unsigned int _nNumSlabsDeltaEkin;
    unsigned int _nNumSlabsDEkinReserve;
    double _dSlabWidthDeltaEkin;

    unsigned long* _nNumMoleculesSumLocal;
    unsigned long* _nNumMoleculesSumGlobal;

    double* _dDelta2EkinTransSumLocal;
    double* _dDelta2EkinTransSumGlobal;
	double* _dDelta2EkinRotSumLocal;
	double* _dDelta2EkinRotSumGlobal;

    // instances / ID
    static unsigned short _nStaticID;

	TemperatureCtrl _target;
	AccumulatorBase* _accumulator;
};

} // namespace tec

struct paramLineHeat
{
	unsigned int nWriteFreq;
	std::string slabsKey;
	std::string slabsVal;
	unsigned int nWriteFreqRegions;
};

class Domain;
class DomainDecompBase;
class TemperatureControl : public ControlInstance
{
public:
	TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp);
    TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
    ~TemperatureControl();

    std::string GetShortName() {return "TeC";}
    void AddRegion(tec::ControlRegion* region);
	void readXML(XMLfileUnits& xmlconfig);
    int GetNumRegions() {return _vecControlRegions.size();}
    tec::ControlRegion* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

    void PrepareRegionSubdivisions();
    void PrepareRegionDataStructures();
    void InitControl(unsigned long simstep);
    void MeasureKineticEnergy(Molecule* mol, unsigned long simstep);
    void CalcGlobalValues(unsigned long simstep);
    void ControlTemperature(Molecule* mol, unsigned long simstep);

    unsigned long GetStart() {return _control.start;}
    unsigned long GetStop()  {return _control.stop;}

    // loops over molecule container
    void DoLoopsOverMolecules(ParticleContainer* particleContainer, unsigned long simstep);

    // heat supply
	void SetSubdivisionHeat(const uint32_t& nSubdivisionHeatType);
	void SetDeltaEkinParameters(const paramLineHeat &paramLine);
	void PrepareDatastructuresHeat();
	void SampleDeltaEkin(Molecule* mol, const double &dAdded2EkinTrans, const double &dAdded2EkinRot);
	void CalcGlobalValuesDeltaEkin();
	void ResetLocalValuesDeltaEkin();
    void WriteDataDeltaEkin(unsigned long simstep);

private:
    std::vector<tec::ControlRegion*> _vecControlRegions;
	TimestepControl _control;

    // heat supply
    bool _bWriteDataDeltaEkin;
    unsigned int _nWriteFreqDeltaEkin;
	unsigned int _nWriteFreqDeltaEkinRegions;
    unsigned int _nNumSlabsDeltaEkin;
	double _dSlabWidthDeltaEkin;

	std::vector<double> _dSlabCoords;
	std::vector<std::vector<unsigned long> > _nNumMoleculesSumLocal;
	std::vector<std::vector<unsigned long> > _nNumMoleculesSumGlobal;
	std::vector<std::vector<double> > _dDelta2EkinTransSumLocal;
	std::vector<std::vector<double> > _dDelta2EkinTransSumGlobal;
	std::vector<std::vector<double> > _dDelta2EkinRotSumLocal;
	std::vector<std::vector<double> > _dDelta2EkinRotSumGlobal;
};

// Accumulate kinetic energy dependent on which translatoric directions should be thermostated

class AccumulatorBase
{
protected:
    AccumulatorBase() {}
public:
    virtual ~AccumulatorBase() {}

public:
    virtual double CalcKineticEnergyContribution(Molecule* mol) = 0;
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot) = 0;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[0] = mol->v(0);
        v_new[0] = v_old[0] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(0, v_new[0]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[1] = mol->v(1);
        v_new[1] = v_old[1] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(1, v_new[1]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[2] = mol->v(2);
        v_new[2] = v_old[2] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(2, v_new[2]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[1] = mol->v(1);
        v_new[0] = v_old[0] * vcorr;
        v_new[1] = v_old[1] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(1, v_new[1]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[2] = mol->v(2);
        v_new[0] = v_old[0] * vcorr;
        v_new[2] = v_old[2] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(2, v_new[2]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[1] = mol->v(1);
        v_old[2] = mol->v(2);
        v_new[1] = v_old[1] * vcorr;
        v_new[2] = v_old[2] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(1, v_new[1]);
        mol->setv(2, v_new[2]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
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
	virtual void ScaleVelocityComponents(Molecule* mol, const double &vcorr, const double &Dcorr,
		double& dAdded2EkinTrans, double& dAdded2EkinRot)
    {
        double v_old[3];
        double v_new[3];
		double d2EkinTrans_old, d2EkinRot_old;
		d2EkinTrans_old = d2EkinRot_old = 0.;

        // calc new velocities
        v_old[0] = mol->v(0);
        v_old[1] = mol->v(1);
        v_old[2] = mol->v(2);
        v_new[0] = v_old[0] * vcorr;
        v_new[1] = v_old[1] * vcorr;
        v_new[2] = v_old[2] * vcorr;

		// store old kinetic energy
		d2EkinTrans_old = mol->U2_trans();
		d2EkinRot_old   = mol->U2_rot();

        // set new velocities
        mol->setv(0, v_new[0]);
        mol->setv(1, v_new[1]);
        mol->setv(2, v_new[2]);

		// set new angular momentum
		mol->scale_D(Dcorr);

		// calc added kinetic energy
		dAdded2EkinTrans = mol->U2_trans() - d2EkinTrans_old;
		dAdded2EkinRot   = mol->U2_rot()   - d2EkinRot_old;
    }
};


#endif /* TEMPERATURECONTROL_H_ */
