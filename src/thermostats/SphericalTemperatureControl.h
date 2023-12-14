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

#include "thermostats/TemperatureControl.h" //defines the accumulator class
#include "ThermostatVariables.h"
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
class SphericalControlRegionT : public SphericalRegionObs {
public:
	SphericalControlRegionT(SphericalTemperatureControl* const parent);
	~SphericalControlRegionT();

	/** @brief Read in XML configuration for TemperatureControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<thermostats>
			<thermostat type="TemperatureControl">
				<control>
					<start>UNSIGNED_LONG</start>           <!-- Timestep turning thermostat ON -->
					<frequency>UNSIGNED_LONG</frequency>   <!-- Thermosatt is active every <frequency>-th time step -->
					<stop>UNSIGNED_LONG</stop>             <!-- Timestep turning thermostat OFF -->
				</control>
				<regions>
					<region>
						<coords>
							<lcx>DOUBLE</lcx> <lcy>DOUBLE</lcy> <lcz>DOUBLE</lcz>
							<ucx>DOUBLE</ucx> <ucy>DOUBLE</ucy> <ucz>DOUBLE</ucz>
						</coords>
						<target>
							<temperature>DOUBLE</temperature>   <!-- target temperature -->
							<ramp>                    <!-- ramp temperature from <start> to <end> value -->
								<start>0.70</start>   <!-- start temperature -->
								<end>0.80</end>       <!-- end temperature -->
								<update>
									<start>UNSIGNED_LONG</start>   <!-- Timestep of ramping start -->
									<stop>UNSIGNED_LONG</stop>     <!-- Timestep of ramping stop -->
									<freq>UNSIGNED_LONG</freq>     <!-- adjust target temperature every <freq>-th time
	 step -->
								</update>
							</ramp>
							<component>1</component>   <!-- target component -->
						</target>
						<settings>
							<numslabs>UNSIGNED_LONG</numslabs>   <!-- Divide region into <numslabs> slabs -->
							<exponent>DOUBLE</exponent>          <!-- Damping exponent of thermostat -->
							<directions>xyz</directions>         <!-- Translational degrees of freedom to be considered
	 for thermostating: x|y|z|xy|xz|yz|xyz -->
						</settings>
						<writefreq>5000</writefreq>         <!-- Log thermostat scaling factors betaTrans and betaRot
	 --> <fileprefix>betalog</fileprefix>    <!-- Prefix of log file -->
					</region>
					</region>   <!-- Add as much regions as you want! -->
						<!-- ...  -->
					</region>
				</regions>
			</thermostat>
		</thermostats>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	unsigned int GetID() { return _nID; }
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

	void registerAsObserver();
	void update(SubjectBase* subject) override;

	// measure added kin. energy
	void InitAddedEkin();
	void writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep);

private:
	// create accumulator object dependent on which translatoric directions should be thermostated (xyz)
	Accumulator* CreateAccumulatorInstance();

	// observer mechanism: update region coords dependent on the interface position, determined by plugin DistControl
	DistControl* getDistControl();

	// instances / ID
	static unsigned short _nStaticID;

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

class Domain;
class SphericalTemperatureControl : public ControlInstance {
public:
	SphericalTemperatureControl();
	~SphericalTemperatureControl();

	std::string getShortName() override { return "TeC"; }
	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(SphericalControlRegionT* region);
	int GetNumRegions() { return _vecControlRegions.size(); }
	SphericalControlRegionT* GetControlRegion(unsigned short nRegionID) {
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
	std::vector<SphericalControlRegionT*> _vecControlRegions;
	unsigned long _nControlFreq;
	unsigned long _nStart;
	unsigned long _nStop;

	enum ControlMethod { VelocityScaling, Andersen, Mixed };
	ControlMethod _method = VelocityScaling;
};

// // Accumulate kinetic energy dependent on which translatoric directions should be thermostated

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
// 		if (_accumulateX) mol->setv(0, mol->v(0) * vcorr);
// 		if (_accumulateY) mol->setv(1, mol->v(1) * vcorr);
// 		if (_accumulateZ) mol->setv(2, mol->v(2) * vcorr);
// 	}
// };

#endif /* SPHERICALTEMPERATURECONTROL_H_ */
