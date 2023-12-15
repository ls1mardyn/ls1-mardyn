/*
 * SphericalTemperatureControl.cpp
 *
 *  Created on: 13.11.2023
 *      Author: jniemann, based on TemperatureControl.cpp by mheinen
 */

#include "thermostats/SphericalTemperatureControl.h"
#include "Domain.h"
#include "WrapOpenMP.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/FileUtils.h"
#include "utils/xmlfileUnits.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// init static ID --> instance counting
unsigned short SphericalControlRegionT::_nStaticID = 0;

bool _bDebugOutput = false;

// class SphericalControlRegionT
SphericalControlRegionT::SphericalControlRegionT(SphericalTemperatureControl* const parent)
	: SphericalRegionObs(parent),
	  _localMethod(VelocityScaling),
	  _dTargetTemperature(0.0),
	  _dTemperatureExponent(0.0),
	  _nTargetComponentID(0),
	  _accumulator(nullptr),
	  _strFilenamePrefixBetaLog("beta_log"),
	  _nWriteFreqBeta(1000),
	  _numSampledConfigs(0),
	  _dBetaTransSumGlobal(0.0),
	  _dBetaRotSumGlobal(0.0),
	  _nuAndersen(0.0),
	  _timestep(0.0),
	  _nuDt(0.0),
	  _rand(Random()),
	  _bIsObserver(false) {
	// ID
	_nID = ++_nStaticID;
}

SphericalControlRegionT::~SphericalControlRegionT() { delete _accumulator; }

Accumulator* SphericalControlRegionT::CreateAccumulatorInstance() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: CreateAccumulatorInstance has been called "<< endl;
	Accumulator* accumulator = new Accumulator(true, true, true);
	return accumulator;
}

void SphericalControlRegionT::readXML(XMLfileUnits& xmlconfig) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: readXML has been called "<< endl;

	Domain* domain = global_simulation->getDomain();
	double ctr[3];
	double rad;
	std::string strVal[3];

	// coordinates
	xmlconfig.getNodeValue("coords/ctrx", ctr[0]);
	xmlconfig.getNodeValue("coords/ctry", ctr[1]);
	xmlconfig.getNodeValue("coords/ctrz", ctr[2]);
	xmlconfig.getNodeValue("coords/radius", rad);
	// read upper corner

	if(_bDebugOutput) global_log->info() << "SphericalTemperatureControl: Center:" << ctr[0] << ", " << ctr[1] << ", " << ctr[2] << ", radius: "<< rad<< endl;

	for (uint8_t d = 0; d < 3; ++d) {
		_dCenter[d] = ctr[d];
	}
	_dRadius = rad;
	_dRadiusSquared = pow(rad, 2);

	// observer mechanism
				std::vector<uint32_t> refCoordsID(6, 0);
				xmlconfig.getNodeValue("coords/lcx@refcoordsID", refCoordsID.at(0));
				xmlconfig.getNodeValue("coords/lcy@refcoordsID", refCoordsID.at(1));
				xmlconfig.getNodeValue("coords/lcz@refcoordsID", refCoordsID.at(2));
				xmlconfig.getNodeValue("coords/ucx@refcoordsID", refCoordsID.at(3));
				xmlconfig.getNodeValue("coords/ucy@refcoordsID", refCoordsID.at(4));
				xmlconfig.getNodeValue("coords/ucz@refcoordsID", refCoordsID.at(5));

				_bIsObserver = (std::accumulate(refCoordsID.begin(), refCoordsID.end(), 0u)) > 0;
				if (true == _bIsObserver) this->PrepareAsObserver(refCoordsID);
				// Registration as observer has to be done later by method prepare_start() when DistControl plugin is present.

	// target values
	xmlconfig.getNodeValue("target/temperature", _dTargetTemperature);
	xmlconfig.getNodeValue("target/component", _nTargetComponentID);


	// ControlMethod "VelocityScaling/Andersen/Mixed"
	std::string methods = "";
	xmlconfig.getNodeValue("method", methods);
	if (methods != "") {
		if (methods == "VelocityScaling") {
			_localMethod = VelocityScaling;

			// init data structures
			this->VelocityScalingInit(xmlconfig);
		} else if (methods == "Andersen") {
			_localMethod = Andersen;
			xmlconfig.getNodeValue("settings/nu", _nuAndersen);
			_timestep = global_simulation->getIntegrator()->getTimestepLength();
			_nuDt = _nuAndersen * _timestep;
		} else {
			global_log->error() << "[TemperatureControl] REGION: Invalid 'method' param: " << methods << std::endl;
			Simulation::exit(-1);
		}
		global_log->info() << "[TemperatureControl] REGION 'method' param: " << methods << std::endl;
	}
	//
	else {
		_localMethod = VelocityScaling;
		global_log->info() << "[TemperatureControl] REGION: no method specified, selecting VelocityScaling"
						   << std::endl;

		// init data structures
		this->VelocityScalingInit(xmlconfig);
	}

	// measure added kin. energy
	_addedEkin.writeFreq = 1000;
	xmlconfig.getNodeValue("added_ekin/writefreq", _addedEkin.writeFreq);
}

void SphericalControlRegionT::VelocityScalingInit(XMLfileUnits& xmlconfig) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: VelocityScalingInit has been called "<< endl;
	// settings
	xmlconfig.getNodeValue("settings/exponent", _dTemperatureExponent);

	// create accumulator instance
	_accumulator = this->CreateAccumulatorInstance();

	// write control for beta_trans and beta_rot log file
	_nWriteFreqBeta = 1000;
	_strFilenamePrefixBetaLog = "beta_log";
	xmlconfig.getNodeValue("writefreq", _nWriteFreqBeta);
	xmlconfig.getNodeValue("fileprefix", _strFilenamePrefixBetaLog);
	if (_nWriteFreqBeta == 0) {
		global_log->warning()
			<< "Temperature Control: write Frequency was specified to be zero. This is NOT allowed. Reset it to 1000."
			<< std::endl;
		_nWriteFreqBeta = 1000;
	}
	this->InitBetaLogfile();
	_localThermVarsThreadBuffer.resize(mardyn_get_max_threads());
	for (auto& localThermVars : _localThermVarsThreadBuffer) {
		localThermVars.resize(1); //_nNumSlabs==1
	}
	_globalThermVars.resize(1); //_nNumSlabs==1

	// init data structure for measure of added kin. energy
	_addedEkin.data.local.resize(1); //_nNumSlabs==1
	_addedEkin.data.global.resize(1); //_nNumSlabs==1

	_addedEkinLocalThreadBuffer.resize(mardyn_get_max_threads());
	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		addedEkinLocal.resize(1); //_nNumSlabs==1
	}

	std::vector<double>& v = _addedEkin.data.local;
	std::fill(v.begin(), v.end(), 0.);

	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		std::fill(addedEkinLocal.begin(), addedEkinLocal.end(), 0.);
	}
}

void SphericalControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: CalcGlobalValues has been called "<< endl;
	if (_localMethod != VelocityScaling) return;
	domainDecomp->collCommInit(4);
	for (unsigned s = 0; s < 1; ++s) { //_nNumSlabs==1
		unsigned long numMolecules{}, numRotationalDOF{};
		double ekinTrans{}, ekinRot{};
		for (const auto& threadBuffers : _localThermVarsThreadBuffer) {
			const auto& localVar = threadBuffers[s];
			numMolecules += localVar._numMolecules;
			numRotationalDOF += localVar._numRotationalDOF;
			ekinTrans += localVar._ekinTrans;
			ekinRot += localVar._ekinRot;
		}
		domainDecomp->collCommAppendUnsLong(numMolecules);
		domainDecomp->collCommAppendUnsLong(numRotationalDOF);
		domainDecomp->collCommAppendDouble(ekinRot);
		domainDecomp->collCommAppendDouble(ekinTrans);
	}
	domainDecomp->collCommAllreduceSum();
	for (unsigned s = 0; s < 1; ++s) { //_nNumSlabs==1
		GlobalThermostatVariables& globalTV = _globalThermVars[s];  // do not forget &
		globalTV._numMolecules = domainDecomp->collCommGetUnsLong();
		globalTV._numRotationalDOF = domainDecomp->collCommGetUnsLong();
		globalTV._ekinRot = domainDecomp->collCommGetDouble();
		globalTV._ekinTrans = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	// Adjust target temperature
	uint64_t simstep = _simulation.getSimulationStep();

	// calc betaTrans, betaRot, and their sum
	double dBetaTransSumSlabs = 0.;
	double dBetaRotSumSlabs = 0.;

	for (unsigned int s = 0; s < 1; ++s) { //_nNumSlabs==1
		GlobalThermostatVariables& globalTV = _globalThermVars[s];  // do not forget &
		if (globalTV._numMolecules < 1)
			globalTV._betaTrans = 1.;
		else
			globalTV._betaTrans = pow(
				3 * globalTV._numMolecules * _dTargetTemperature / globalTV._ekinTrans,
				_dTemperatureExponent);

		if (globalTV._numRotationalDOF < 1)
			globalTV._betaRot = 1.;
		else
			globalTV._betaRot =
				pow(globalTV._numRotationalDOF * _dTargetTemperature / globalTV._ekinRot, _dTemperatureExponent);

		// calc sums over all slabs
		dBetaTransSumSlabs += globalTV._betaTrans;
		dBetaRotSumSlabs += globalTV._betaRot;
	}
	// calc ensemble average of beta_trans, beta_rot
	_dBetaTransSumGlobal += dBetaTransSumSlabs;
	_dBetaRotSumGlobal += dBetaRotSumSlabs;
	_numSampledConfigs++;
}

void SphericalControlRegionT::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* /*domainDecomp*/) {
	if (_localMethod != VelocityScaling) return;
	// check componentID
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule inside control region
	double distanceFromCenterSquared = 0;
	for(int d = 0; d<3; d++){
		distanceFromCenterSquared += mol->r(d)-_dCenter[d];
	}
	if(distanceFromCenterSquared > _dRadiusSquared) return;

	// sum up transl. kinetic energy (2x)
	/*
		double vx = mol->v(0);
	//    double vy = mol->v(1);
		double vz = mol->v(2);
		double m  = mol->mass();

	//    _d2EkinTransLocal += m*(vx*vx + vy*vy);
		_d2EkinTransLocal[nPosIndex] += m*(vx*vx + vz*vz);
	*/

	auto myThreadNum = mardyn_get_thread_num();
	double nPosIndex = 0;
	LocalThermostatVariables& localTV = _localThermVarsThreadBuffer[myThreadNum].at(nPosIndex);  // do not forget &  
	localTV._ekinTrans += _accumulator->CalcKineticEnergyContribution(mol);

	// sum up rot. kinetic energy (2x)
	double dDummy = 0.;
	double ekinRot = 0.;

	mol->calculate_mv2_Iw2(dDummy, ekinRot);
	localTV._ekinRot += ekinRot;

	// count num molecules
	localTV._numMolecules++;

	// count rotational DOF
	localTV._numRotationalDOF += mol->component()->getRotationalDegreesOfFreedom();
}

void SphericalControlRegionT::ControlTemperature(Molecule* mol) {
	// check componentID
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule is inside
	double distanceFromCenterSquared = 0;
	for(int d = 0; d<3; d++){
		distanceFromCenterSquared += mol->r(d)-_dCenter[d];
	}
	if(distanceFromCenterSquared > _dRadiusSquared) return;

	// check for method
	if (_localMethod == VelocityScaling) {
		
		double nPosIndex = 0;
		GlobalThermostatVariables& globalTV = _globalThermVars[nPosIndex];  // do not forget &
		if (globalTV._numMolecules < 1) return;
		// scale velocity
		double vcorr = 2. - 1. / globalTV._betaTrans;
		double Dcorr = 2. - 1. / globalTV._betaRot;

		// measure added kin. energy
		double v2_old = mol->v2();

		_accumulator->ScaleVelocityComponents(mol, vcorr);

		// measure added kin. energy
		double v2_new = mol->v2();
		int mythread = mardyn_get_thread_num();
		_addedEkinLocalThreadBuffer[mythread].at(nPosIndex) += (v2_new - v2_old);

		mol->scale_D(Dcorr);
	} else if (_localMethod == Andersen) {
		double stdDevTrans, stdDevRot;
		if (_rand.rnd() < _nuDt) {
			stdDevTrans = sqrt(_dTargetTemperature / mol->mass());
			for (unsigned short d = 0; d < 3; d++) {
				stdDevRot = sqrt(_dTargetTemperature * mol->getI(d));
				mol->setv(d, _rand.gaussDeviate(stdDevTrans));
				mol->setD(d, _rand.gaussDeviate(stdDevRot));
			}
		}
	} else {
		global_log->error() << "[TemperatureControl] Invalid localMethod param: " << _localMethod << std::endl;
		Simulation::exit(-1);
	}
}

void SphericalControlRegionT::ResetLocalValues() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::ResetLocalValues has been called "<< endl;
	// reset local values
	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		for (unsigned int s = 0; s < 1; ++s) {//_nNumSlabs==1
			_localThermVarsThreadBuffer[thread][s].clear();
		}
	}
}

void SphericalControlRegionT::InitBetaLogfile() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::InitBetaLogfile has been called "<< endl;
	if (_localMethod == VelocityScaling) {
#ifdef ENABLE_MPI
		DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition());
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank != 0) return;
#endif
		// touch file
		const std::string fname = _strFilenamePrefixBetaLog + "_reg" + std::to_string(this->GetID()) + ".dat";
		std::ofstream ofs;
		ofs.open(fname, std::ios::out);
		ofs << setw(12) << "simstep"
			<< setw(24) << "dBetaTrans"
			<< setw(24) << "dBetaRot"
			<< std::endl;
		ofs.close();
	}
}

void SphericalControlRegionT::WriteBetaLogfile(unsigned long simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::WriteBetaLogfile has been called "<< endl;
	if (_localMethod != VelocityScaling) {
		return;
	}
	if (0 != (simstep % _nWriteFreqBeta)) {
		return;
	}

	DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0) return;
#endif

	// double dInvNumConfigsSlabs = 1. / (double)(_numSampledConfigs*_nNumSlabs);
	double dInvNumConfigsSlabs = 1. / (double)(_numSampledConfigs); //_nNumSlabs==1
	double dBetaTrans = _dBetaTransSumGlobal * dInvNumConfigsSlabs;
	double dBetaRot = _dBetaRotSumGlobal * dInvNumConfigsSlabs;

	// writing to file
	const std::string fname = _strFilenamePrefixBetaLog + "_reg" + std::to_string(this->GetID()) + ".dat";
	std::ofstream ofs;
	ofs.open(fname, std::ios::app);
	ofs << setw(12) << simstep
		<< FORMAT_SCI_MAX_DIGITS << dBetaTrans
		<< FORMAT_SCI_MAX_DIGITS << dBetaRot
		<< std::endl;
	ofs.close();

	// reset averaged values
	_numSampledConfigs = 0;
	_dBetaTransSumGlobal = 0.;
	_dBetaRotSumGlobal = 0.;
}

DistControl* SphericalControlRegionT::getDistControl() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::getDistControl has been called "<< endl;
	DistControl* distControl = nullptr;
	std::list<PluginBase*>& plugins = *(global_simulation->getPluginList());
	for (auto&& pit : plugins) {
		std::string name = pit->getPluginName();
		if (name == "DistControl") {
			distControl = dynamic_cast<DistControl*>(pit);
		}
	}
	return distControl;
}
void SphericalControlRegionT::registerAsObserver() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::registerAsObserver has been called "<< endl;
	if (true == _bIsObserver) {
		DistControl* distControl = this->getDistControl();
		if (distControl != nullptr)
			distControl->registerObserver(this);
		else {
			global_log->error() << "TemperatureControl->region[" << this->GetID()
								<< "]: Initialization of plugin DistControl is needed before! Program exit..." << endl;
			Simulation::exit(-1);
		}
	}
}

void SphericalControlRegionT::update(SubjectBase* subject) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::update  has been called "<< endl;
	SphericalRegionObs::update(subject);
}

void SphericalControlRegionT::InitAddedEkin() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::InitAddedEkin has been called "<< endl;
	if (_localMethod == VelocityScaling) {
#ifdef ENABLE_MPI
		DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition());
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank != 0) return;
#endif
		// touch file
		const std::string fname = "addedEkin_reg" + std::to_string(this->GetID()) + "_cid" + std::to_string(_nTargetComponentID) + ".dat";
		std::ofstream ofs;
		ofs.open(fname, std::ios::out);
		ofs << setw(12) << "simstep";
		for (int i = 0; i < 1; ++i) { //_nNumSlabs==1
			std::string s = "bin" + std::to_string(i+1);
			ofs << setw(24) << s;
		}
		ofs << std::endl;
		ofs.close();
	}
}

void SphericalControlRegionT::writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::writeAddedEkin has been called "<< endl;
	if (_localMethod != VelocityScaling) return;

	if (simstep % _addedEkin.writeFreq != 0) return;

	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		mardyn_assert(_addedEkin.data.local.size() == 1);//_nNumSlabs==1
		for (size_t slabID = 0; slabID < 1; ++slabID) { //_nNumSlabs==1
			_addedEkin.data.local[slabID] += _addedEkinLocalThreadBuffer[thread][slabID];
		}
	}
	// calc global values
#ifdef ENABLE_MPI
	MPI_Reduce(_addedEkin.data.local.data(), _addedEkin.data.global.data(), _addedEkin.data.local.size(), MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);
#else
	std::memcpy(_addedEkin.data.global.data(), _addedEkin.data.local.data(), _addedEkin.data.local.size());
#endif

	// reset local values
	std::vector<double>& vl = _addedEkin.data.local;
	std::fill(vl.begin(), vl.end(), 0.);
	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		std::fill(addedEkinLocal.begin(), addedEkinLocal.end(), 0.);
	}

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0) return;
#endif

	// dekin = d2ekin * 0.5
	std::vector<double>& vg = _addedEkin.data.global;
	for (double& it : vg) {
		it *= 0.5;
	}

	// writing .dat-files
	const std::string fname = "addedEkin_reg" + std::to_string(this->GetID()) + "_cid" + std::to_string(_nTargetComponentID) + ".dat";
	std::ofstream ofs;
	ofs.open(fname, std::ios::app);
	
	ofs << setw(12) << simstep;
	for (double& it : vg) {
		ofs << FORMAT_SCI_MAX_DIGITS << it;
	}
	ofs << std::endl;
	ofs.close();
}

// class SphericalTemperatureControl
SphericalTemperatureControl::SphericalTemperatureControl() {}

SphericalTemperatureControl::~SphericalTemperatureControl() {
	for (auto region : _vecControlRegions) {
		delete region;
	}
}

void SphericalTemperatureControl::readXML(XMLfileUnits& xmlconfig) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::readXML has been called "<< endl;
	// control
	xmlconfig.getNodeValue("control/start", _nStart);
	xmlconfig.getNodeValue("control/frequency", _nControlFreq);
	xmlconfig.getNodeValue("control/stop", _nStop);
	global_log->info() << "Start control from simstep: " << _nStart << endl;
	global_log->info() << "Control with frequency: " << _nControlFreq << endl;
	global_log->info() << "Stop control at simstep: " << _nStop << endl;

	// turn on/off explosion heuristics
	// domain->setExplosionHeuristics(bUseExplosionHeuristics);

	// add regions
	uint32_t numRegions = 0;
	XMLfile::Query query = xmlconfig.query("regions/region");
	numRegions = query.card();
	global_log->info() << "Number of control regions: " << numRegions << endl;
	if (numRegions < 1) {
		global_log->warning() << "No region parameters specified." << endl;
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputRegionIter;
	for (outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++) {
		xmlconfig.changecurrentnode(outputRegionIter);
		SphericalControlRegionT* region = new SphericalControlRegionT(this);
		region->readXML(xmlconfig);
		this->AddRegion(region);
	}

	bool Vel = false;
	bool And = false;
	// check for mixed mode
	for (auto&& reg : _vecControlRegions) {
		if (reg->_localMethod == SphericalControlRegionT::LocalControlMethod::VelocityScaling)
			Vel = true;
		else if (reg->_localMethod == SphericalControlRegionT::LocalControlMethod::Andersen)
			And = true;
	}
	if (Vel && And) {
		_method = Mixed;
		global_log->info() << "[TemperatureControl] Mixed methods across regions\n";
	} else if (!Vel && And) {
		_method = Andersen;
		global_log->info() << "[TemperatureControl] Andersen in all regions\n";
	} else {
		_method = VelocityScaling;
		global_log->info() << "[TemperatureControl] VelocityControl in all regions\n";
	}
}

void SphericalTemperatureControl::prepare_start() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::prepare_start has been called "<< endl;
	for (auto&& reg : _vecControlRegions) {
		reg->registerAsObserver();
	}
	this->InitBetaLogfiles();
	this->InitAddedEkin();
}

void SphericalTemperatureControl::AddRegion(SphericalControlRegionT* region) { _vecControlRegions.push_back(region); }

void SphericalTemperatureControl::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->MeasureKineticEnergy(mol, domainDecomp);
}

void SphericalTemperatureControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::CalcGlobalValues has been called "<< endl;
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->CalcGlobalValues(domainDecomp);
}

void SphericalTemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) {
		reg->ControlTemperature(mol);
	}
}

void SphericalTemperatureControl::Init(unsigned long simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::Init has been called "<< endl;
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->ResetLocalValues();
}

void SphericalTemperatureControl::InitBetaLogfiles() {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::InitBetaLogfiles has been called "<< endl;
	for (auto&& reg : _vecControlRegions) reg->InitBetaLogfile();
}

void SphericalTemperatureControl::WriteBetaLogfiles(unsigned long simstep) {
	for (auto&& reg : _vecControlRegions) reg->WriteBetaLogfile(simstep);
}

void SphericalTemperatureControl::InitAddedEkin() {
	for (auto&& reg : _vecControlRegions) reg->InitAddedEkin();
}

void SphericalTemperatureControl::writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep) {
	for (auto&& reg : _vecControlRegions) reg->writeAddedEkin(domainDecomp, simstep);
}

/**
 * @brief Decide which ControlMethod to use
 *
 * @param domainDecomposition
 * @param particleContainer
 * @param simstep
 */
void SphericalTemperatureControl::DoLoopsOverMolecules(DomainDecompBase* domainDecomposition,
											  ParticleContainer* particleContainer, const unsigned long simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::DoLoopsOverMolecules has been called "<< endl;
	if (_method == VelocityScaling || _method == Mixed) {
		this->VelocityScalingPreparation(domainDecomposition, particleContainer, simstep);
		global_log->debug() << "[TemperatureControl] VelocityScalingPreparation\n";
	}

	// iterate over all molecules. ControlTemperature depends on _localMethod for Region molecule is in
#if defined(_OPENMP)
// gcc 7 and 8 implicitly declare const variables shared regardless of the default.
// Additionally, explicitly marking them as shared produces an error.
// All other compilers need the explicit declaration when default is none.
// Therefore set default to shared here...
#pragma omp parallel default(shared) shared(particleContainer)
#endif
	for (auto tM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
		// control temperature
		this->ControlTemperature(&(*tM), simstep);
	}

	// measure added kin. energy
	this->writeAddedEkin(domainDecomposition, simstep);
}

/**
 * @brief Prepare for VelocityScaling control method
 *
 * @param domainDecomposition
 * @param particleContainer
 * @param simstep
 */
void SphericalTemperatureControl::VelocityScalingPreparation(DomainDecompBase* domainDecomposition,
													ParticleContainer* particleContainer, const unsigned long simstep) {
	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::VelocityScalingPreparation has been called "<< endl;
	// respect start/stop
	if (this->GetStart() <= simstep && this->GetStop() > simstep) {
		// init temperature control
		this->Init(simstep);
#if defined(_OPENMP)
#pragma omp parallel
#endif
		for (auto tM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
			// measure kinetic energy
			this->MeasureKineticEnergy(&(*tM), domainDecomposition, simstep);
		}

		// calc global values
		this->CalcGlobalValues(domainDecomposition, simstep);

		// write beta_trans, beta_rot log-files
		this->WriteBetaLogfiles(simstep);
	}
}
