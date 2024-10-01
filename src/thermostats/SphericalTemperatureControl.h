/*
 * SphericalTemperatureControl.h
 *
 *  Created on: 13.11.2023
 *      Author: jniemann, based on TemperatureControl plugin by mheinen
 * 
 * Description: In vle of droplets, the Temperature may develop differently in both phases, if only a global thermostat is used, 
 * 				due to systematic differences in how well the integrators function in both phases.
 * 				To correct this, a thermostat with the ability to adjust the Temperature inside and outside a spherical Region serperatly from one another is needed.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "ThermostatVariables.h"
#include "thermostats/Accumulator.h"
#include "integrators/Integrator.h"
#include "molecules/Molecule.h"
#include "plugins/NEMD/DistControl.h"
#include "utils/CommVar.h"
#include "utils/ObserverBase.h"
#include "utils/Random.h"
#include "utils/Region.h"


/*  XML setup:
<thermostat type="SphericalTemperatureControl">
	<control>
		<start>INT</start>   						<!-- Simstep to start sampling; default 0 -->
		<frequency>INT</frequency>					<!-- Sampling Frequency; default 1 -->
		<stop>INT</stop>							<!-- Simstep to stop sampling; default MAX -->
	</control>
	<regions>
		<region>
			<coords>
				<lcx>DOUBLE</lcx> <lcy>DOUBLE</lcy> <lcz>DOUBLE</lcz>			<!-- lower Corner of Box (outer region) -->
				<ucx>box</ucx> <ucy>box</ucy> <ucz>box</ucz>					<!-- upper Corner of Box (outer region) -->
				<ctrx>DOUBLE</ctrx> <ctry>DOUBLE</ctry> <ctrz>DOUBLE</ctrz>		<!-- Center of Spherical Region -->
				<radius>DOUBLE</radius>											<!-- Radius that seperates the two regions -->
			</coords>
			<target>
				<temperature>DOUBLE</temperature>
				<component>0</component>
			</target>
			<settings>
				<numslabs>1</numslabs>							<!-- Legacy; please ignore / don't change. values different from 1 will not have an effect. -->
				<exponent>DOUBLE</exponent>
				<directions>xyz</directions>
			</settings>
			<writefreq>INT</writefreq>							<!-- Write Frequency; default 1000 -->
			<fileprefix>STRING</fileprefix>
		</region>
	</regions>
</thermostat>
 */

class DistControl;
class XMLfileUnits;
class DomainDecompBase;
class ParticleContainer;
class SphericalTemperatureControl;

class AbstrControlRegionT {
public:
	AbstrControlRegionT();
	virtual ~AbstrControlRegionT();

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

	virtual bool ContainsMolecule(Molecule* mol) = 0;

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

	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<LocalThermostatVariables> _localThermVarsThreadBuffer;
	GlobalThermostatVariables _globalThermVars;

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
		CommVar<double> data;  // \Delta E_kin * 2/m = v^2_2 - v^2_1
	} _addedEkin;

	/**
	 * Thread buffer for the local kinetic energy.
	 */
	std::vector<double> _addedEkinLocalThreadBuffer;

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
	~SphericalControlRegionT() override;

	bool ContainsMolecule(Molecule* mol) override;

	void readXML(XMLfileUnits& xmlconfig) override;
	unsigned int GetID() override { return _nID; };

	enum LocalControlMethod {
		VelocityScaling,
		Andersen,
	};
	LocalControlMethod _localMethod;

	void registerAsObserver() override;
	void update(SubjectBase* subject) override;

private:
	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<LocalThermostatVariables> _localThermVarsThreadBuffer;
	GlobalThermostatVariables _globalThermVars;

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
	~SphereComplementControlRegionT() override;

	bool ContainsMolecule(Molecule* mol) override;
	void readXML(XMLfileUnits& xmlconfig) override;
	unsigned int GetID() override { return SphericalRegion::_nID; };

	enum LocalControlMethod {
		VelocityScaling,
		Andersen,
	};
	LocalControlMethod _localMethod;

	void registerAsObserver() override;
	void update(SubjectBase* subject) override;

private:
	/**
	 * Thread buffer for the local thermostat variables.
	 */
	std::vector<LocalThermostatVariables> _localThermVarsThreadBuffer;
	GlobalThermostatVariables _globalThermVars;

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
	unsigned long _nControlFreq {1};
	unsigned long _nStart {0};
	unsigned long _nStop {std::numeric_limits<unsigned long>::max()};;

	enum ControlMethod { VelocityScaling, Andersen, Mixed };
	ControlMethod _method = VelocityScaling;
};
