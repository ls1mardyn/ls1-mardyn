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
class ControlRegionT
{
public:
	ControlRegionT(double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp,
			double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections,
			unsigned long nWriteFreqBeta, std::string strFilenamePrefix);
	~ControlRegionT();

	unsigned int GetID(){return _nID;}
	void Init();

	double* GetLowerCorner() {return _dLowerCorner;}
	double* GetUpperCorner() {return _dUpperCorner;}
	void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal;}
	void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal;}
	double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
	void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
	void CalcGlobalValues(DomainDecompBase* domainDecomp);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp);

	void ControlTemperature(Molecule* mol);

	void ResetLocalValues();

	// beta log file
	void InitBetaLogfile();
	void WriteBetaLogfile(unsigned long simstep);

private:
	unsigned short _nID;

	// instances / ID
	static unsigned short _nStaticID;

	double _dLowerCorner[3];
	double _dUpperCorner[3];

	unsigned int _nNumSlabs;
	double _dSlabWidth;

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

	std::string _strFilenamePrefixBetaLog;
	unsigned long _nWriteFreqBeta;
	unsigned long _numSampledConfigs;
	double _dBetaTransSumGlobal;
	double _dBetaRotSumGlobal;
};


class Domain;
class XMLfileUnits;
class TemperatureControl
{
public:
	TemperatureControl();
	TemperatureControl(unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
	~TemperatureControl();

	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp,
			double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections,
			unsigned long nWriteFreqBeta, std::string strFilenamePrefix);
	int GetNumRegions() {return _vecControlRegions.size();}
	ControlRegionT* GetControlRegion(unsigned short nRegionID) {return &(_vecControlRegions.at(nRegionID-1) ); }  // vector index starts with 0, region index with 1

	void Init(unsigned long simstep);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
	void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
	void ControlTemperature(Molecule* mol, unsigned long simstep);

	unsigned long GetStart() {return _nStart;}
	unsigned long GetStop()  {return _nStop;}

	// beta log file
	void InitBetaLogfiles();
	void WriteBetaLogfiles(unsigned long simstep);

	// loops over molecule container
	void DoLoopsOverMolecules(DomainDecompBase*, ParticleContainer* particleContainer, unsigned long simstep);

private:
	std::vector<ControlRegionT> _vecControlRegions;
	unsigned long _nControlFreq;
	unsigned long _nStart;
	unsigned long _nStop;
};


// Accumulate kinetic energy dependent on which translatoric directions should be thermostated

class AccumulatorBase
{
protected:
	AccumulatorBase() {}
	virtual ~AccumulatorBase() {}

public:
	virtual double CalcKineticEnergyContribution(Molecule* mol) = 0;
	virtual void ScaleVelocityComponents(Molecule* mol, double vcorr) = 0;
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

#endif /* TEMPERATURECONTROL_H_ */
