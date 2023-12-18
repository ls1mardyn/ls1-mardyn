/*
 * SphericalTemperatureControl.h
 *
 *  Created on: 13.11.2023
 *      Author: jniemann, based on TemperatureControl.cpp by mheinen
 */

#ifndef SPHERICALTEMPERATURECONTROL_H_
#define SPHERICALTEMPERATURECONTROL_H_

#include <cstdint>
#include <string>
#include <vector>

#include "ThermostatVariables.h"
#include "Accumulator.h"
#include "integrators/Integrator.h"
#include "molecules/Molecule.h"
#include "plugins/NEMD/DistControl.h"
#include "utils/CommVar.h"
#include "utils/ObserverBase.h"
#include "utils/Random.h"
#include "utils/Region.h"

class DistControl;
class XMLfileUnits;
class DomainDecompBase;
class ParticleContainer;
// class Accumulator;
class SphericalTemperatureControl;

class AbstrControlRegionT {
public:
	AbstrControlRegionT();
	~AbstrControlRegionT();

	/** @brief Read in XML configuration for TemperatureControl and all its included objects.
	*TODO: write something here.
	*/


	virtual void readXML(XMLfileUnits& xmlconfig) = 0;
	
	virtual unsigned int GetID() = 0;
	void VelocityScalingInit(XMLfileUnits& xmlconfig);
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

	virtual void registerAsObserver() = 0;
	virtual void update(SubjectBase* subject) = 0;

	// measure added kin. energy
	void InitAddedEkin();
	void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

protected:
	// observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	DistControl* getDistControl();

	// instances / ID
	static unsigned short _nStaticID;

private:
	// create accumulator object dependent on which translatoric directions should be thermostated (xyz)
	Accumulator* CreateAccumulatorInstance();



	// double radiusBondary;

	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<std::vector<LocalThermostatVariables>> _localThermVarsThreadBuffer;
	std::vector<GlobalThermostatVariables> _globalThermVars;

	double _dTargetTemperature;
	double _dTemperatureExponent;
	unsigned int _nTargetComponentID;

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
		CommVar<std::vector<double>> data;  // \Delta E_kin * 2/m = v^2_2 - v^2_1
	} _addedEkin;

	/**
	 * Thread buffer for the local kinetic energy.
	 */
	std::vector<std::vector<double>> _addedEkinLocalThreadBuffer;

	struct Ramp {
		bool enabled;
		float start, end, delta, slope;
		struct Update {
			uint64_t start, stop, delta, elapsed;
			uint32_t freq;
		} update;
	} _ramp;
};



class SphericalControlRegionT : public AbstrControlRegionT, public SphericalRegionObs {
public:
	SphericalControlRegionT(SphericalTemperatureControl* const parent);
	~SphericalControlRegionT();

	/** @brief Read in XML configuration for TemperatureControl and all its included objects.
	*TODO: write something here.
	*/


	void readXML(XMLfileUnits& xmlconfig);
	unsigned int GetID()  { return _nID; };
	// void VelocityScalingInit(XMLfileUnits& xmlconfig);
	// void CalcGlobalValues(DomainDecompBase* domainDecomp);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp);
	void ControlTemperature(Molecule* mol);
	// void ResetLocalValues();

	// // beta log file
	// void InitBetaLogfile();
	// void WriteBetaLogfile(unsigned long simstep);

	enum LocalControlMethod {
		VelocityScaling,
		Andersen,
	};
	LocalControlMethod _localMethod;

	void registerAsObserver();
	void update(SubjectBase* subject) override;

	// // measure added kin. energy
	// void InitAddedEkin();
	// void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	// // create accumulator object dependent on which translatoric directions should be thermostated (xyz)
	// Accumulator* CreateAccumulatorInstance();

	// // observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	// DistControl* getDistControl();

	// // instances / ID
	// static unsigned short _nStaticID;

	// double radiusBondary;

	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<std::vector<LocalThermostatVariables>> _localThermVarsThreadBuffer;
	std::vector<GlobalThermostatVariables> _globalThermVars;

	double _dTargetTemperature;
	double _dTemperatureExponent;
	unsigned int _nTargetComponentID;

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
		CommVar<std::vector<double>> data;  // \Delta E_kin * 2/m = v^2_2 - v^2_1
	} _addedEkin;

	/**
	 * Thread buffer for the local kinetic energy.
	 */
	std::vector<std::vector<double>> _addedEkinLocalThreadBuffer;
};


class SphereComplementControlRegionT : public AbstrControlRegionT, public SphereComplementRegionObs {
public: 
	SphereComplementControlRegionT(SphericalTemperatureControl* const parent);
	~SphereComplementControlRegionT();

	void readXML(XMLfileUnits& xmlconfig);
	unsigned int GetID()  { return SphericalRegion::_nID; };
	// void VelocityScalingInit(XMLfileUnits& xmlconfig);
	// void CalcGlobalValues(DomainDecompBase* domainDecomp);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp);
	void ControlTemperature(Molecule* mol);
	// void ResetLocalValues();

	// // beta log file
	// void InitBetaLogfile();
	// void WriteBetaLogfile(unsigned long simstep);

	enum LocalControlMethod {
		VelocityScaling,
		Andersen,
	};
	LocalControlMethod _localMethod;

	void registerAsObserver();
	void update(SubjectBase* subject) override;

	// // measure added kin. energy
	// void InitAddedEkin();
	// void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	// // create accumulator object dependent on which translatoric directions should be thermostated (xyz)
	// Accumulator* CreateAccumulatorInstance();

	// // observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	// DistControl* getDistControl();

	// // instances / ID
	// static unsigned short _nStaticID;

	// double radiusBondary;

	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<std::vector<LocalThermostatVariables>> _localThermVarsThreadBuffer;
	std::vector<GlobalThermostatVariables> _globalThermVars;

	double _dTargetTemperature;
	double _dTemperatureExponent;
	unsigned int _nTargetComponentID;

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
		CommVar<std::vector<double>> data;  // \Delta E_kin * 2/m = v^2_2 - v^2_1
	} _addedEkin;

	/**
	 * Thread buffer for the local kinetic energy.
	 */
	std::vector<std::vector<double>> _addedEkinLocalThreadBuffer;
};









class Domain;
class SphericalTemperatureControl : public ControlInstance {
public:
	SphericalTemperatureControl();
	~SphericalTemperatureControl();

	std::string getShortName() override { return "TeC"; }
	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(AbstrControlRegionT* region);
	int GetNumRegions() { return _vecControlRegions.size(); }
	AbstrControlRegionT* GetControlRegion(unsigned short nRegionID) {
		return _vecControlRegions.at(nRegionID - 1);
	}  // vector index starts with 0, region index with 1
	void prepare_start();

	void Init(unsigned long simstep);
	void MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep);
	void CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep);
	void ControlTemperature(Molecule* mol, unsigned long simstep);

	unsigned long GetStart() { return _nStart; }
	unsigned long GetStop() { return _nStop; }

	// beta log file
	void InitBetaLogfiles();
	void WriteBetaLogfiles(unsigned long simstep);

	// loops over molecule container
	void DoLoopsOverMolecules(DomainDecompBase*, ParticleContainer* particleContainer, unsigned long simstep);
	void VelocityScalingPreparation(DomainDecompBase*, ParticleContainer*, unsigned long simstep);

	// measure added kin. energy
	void InitAddedEkin();
	void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	std::vector<AbstrControlRegionT*> _vecControlRegions;
	unsigned long _nControlFreq;
	unsigned long _nStart;
	unsigned long _nStop;

	enum ControlMethod { VelocityScaling, Andersen, Mixed };
	ControlMethod _method = VelocityScaling;
};




// // // Accumulate kinetic energy dependent on which translatoric directions should be thermostated
// class Accumulator {
// private:
// 	bool _accumulateX, _accumulateY, _accumulateZ;

// public:
// 	Accumulator(bool accX, bool accY, bool accZ) : _accumulateX(accX), _accumulateY(accY), _accumulateZ(accZ) {}

// 	double CalcKineticEnergyContribution(Molecule* mol) {
// 		double vx = _accumulateX ? mol->v(0) : 0.0;
// 		double vy = _accumulateY ? mol->v(1) : 0.0;
// 		double vz = _accumulateZ ? mol->v(2) : 0.0;
// 		double m = mol->mass();

// 		return m * (vx * vx + vy * vy + vz * vz);
// 	}
// 	void ScaleVelocityComponents(Molecule* mol, double vcorr) {
// 		// global_log->info() << "[SphericalTemperatureControl]: Accumulator::ScaleVelocityComponents(Molecule " << mol->getID()<< ", vcorr "<<vcorr<<"); "<< endl;
// 		for(int d=0; d<3; d++){
// 			// global_log->info() << "[SphericalTemperatureControl]: before correction: " << mol->getID()<< "->v("<<d<<") == "<<mol->v(d)<<" "<< endl;
// 		}
// 		if (_accumulateX) mol->setv(0, mol->v(0) * vcorr);
// 		if (_accumulateY) mol->setv(1, mol->v(1) * vcorr);
// 		if (_accumulateZ) mol->setv(2, mol->v(2) * vcorr);
// 		for(int d=0; d<3; d++){
// 			// global_log->info() << "[SphericalTemperatureControl]: after correction: " << mol->getID()<< "->v("<<d<<") == "<<mol->v(d)<<" "<< endl;
// 		}
// 	}
// };

#endif /* SPHERICALTEMPERATURECONTROL_H_ */
