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
#include <cstdint>

#include "molecules/Molecule.h"
#include "ThermostatVariables.h"
#include "utils/Random.h"
#include "integrators/Integrator.h"
#include "plugins/NEMD/DistControl.h"
#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "utils/CommVar.h"

class DistControl;
class XMLfileUnits;
class DomainDecompBase;
class ParticleContainer;
class Accumulator;
class TemperatureControl;
class ControlRegionT : public CuboidRegionObs
{
public:
	ControlRegionT(TemperatureControl* const parent);
	~ControlRegionT();

	void readXML(XMLfileUnits& xmlconfig);
	unsigned int GetID(){return _nID;}
	void VelocityScalingInit(XMLfileUnits &xmlconfig, std::string strDirections);
	void CalcGlobalValues(DomainDecompBase* domainDecomp);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp);
	void ControlTemperature(Molecule* mol);
	void ResetLocalValues();

	// beta log file
	void InitBetaLogfile();
	void WriteBetaLogfile(unsigned long simstep);

	enum LocalControlMethod {
		VelocityScaling,
		Andersen,
	};
	LocalControlMethod _localMethod;

	void registerAsObserver();

	// measure added kin. energy
	void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	// create accumulator object dependent on which translatoric directions should be thermostated (xyz)
	Accumulator* CreateAccumulatorInstance(std::string strTransDirections);

	// observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	DistControl* getDistControl();

	// instances / ID
	static unsigned short _nStaticID;

	unsigned int _nNumSlabs;
	double _dSlabWidth;

	std::vector<LocalAndGlobalThermostatVariables> _thermVars;

	double _dTargetTemperature;
	double _dTemperatureExponent;
	unsigned int _nTargetComponentID;
	unsigned short _nNumThermostatedTransDirections;

	Accumulator* _accumulator;

	std::string _strFilenamePrefixBetaLog;
	unsigned long _nWriteFreqBeta;
	unsigned long _numSampledConfigs;
	double _dBetaTransSumGlobal;
	double _dBetaRotSumGlobal;
	double _nuAndersen;
	double _timestep;
	double _nuDt;
	Random _rand;

	bool _bIsObserver;

	struct AddedEkin {
		uint32_t writeFreq;
		CommVar<std::vector<double> > data;  // \Delta E_kin * 2/m = v^2_2 - v^2_1
	} _addedEkin;
};


class Domain;
class TemperatureControl : public ControlInstance
{
public:
	TemperatureControl();
	~TemperatureControl();

	std::string getShortName() override {return "TeC";}
	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(ControlRegionT* region);
	int GetNumRegions() {return _vecControlRegions.size();}
	ControlRegionT* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1
	void prepare_start();

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
	void VelocityScalingPreparation(DomainDecompBase *, ParticleContainer *, unsigned long simstep);

	// measure added kin. energy
	void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	std::vector<ControlRegionT*> _vecControlRegions;
	unsigned long _nControlFreq;
	unsigned long _nStart;
	unsigned long _nStop;

	enum ControlMethod {
		VelocityScaling,
		Andersen,
		Mixed
	};
	ControlMethod _method = VelocityScaling;
};


// Accumulate kinetic energy dependent on which translatoric directions should be thermostated

class Accumulator
{
private:
	bool _accumulateX, _accumulateY, _accumulateZ;

public:
	Accumulator(bool accX, bool accY, bool accZ) :
			_accumulateX(accX), _accumulateY(accY), _accumulateZ(accZ) {
	}

	double CalcKineticEnergyContribution(Molecule* mol) {
		double vx = _accumulateX ? mol->v(0) : 0.0;
		double vy = _accumulateY ? mol->v(1) : 0.0;
		double vz = _accumulateZ ? mol->v(2) : 0.0;
		double m  = mol->mass();

		return m * (vx*vx + vy*vy + vz*vz);
	}
	void ScaleVelocityComponents(Molecule* mol, double vcorr) {
		if (_accumulateX) mol->setv(0, mol->v(0) * vcorr);
		if (_accumulateY) mol->setv(1, mol->v(1) * vcorr);
		if (_accumulateZ) mol->setv(2, mol->v(2) * vcorr);
	}
};

#endif /* TEMPERATURECONTROL_H_ */
