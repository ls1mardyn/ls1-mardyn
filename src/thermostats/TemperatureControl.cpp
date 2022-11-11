/*
 * TemperatureControl.cpp
 *
 *  Created on: 27.05.2015
 *      Author: mheinen
 */

#include "thermostats/TemperatureControl.h"
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
unsigned short ControlRegionT::_nStaticID = 0;

// class ControlRegionT

ControlRegionT::ControlRegionT(TemperatureControl* const parent)
	: CuboidRegionObs(parent),
	  _localMethod(VelocityScaling),
	  _nNumSlabs(1),
	  _dSlabWidth(0.0),
	  _dTargetTemperature(0.0),
	  _dTemperatureExponent(0.0),
	  _nTargetComponentID(0),
	  _nNumThermostatedTransDirections(0),
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

ControlRegionT::~ControlRegionT() { delete _accumulator; }

Accumulator* ControlRegionT::CreateAccumulatorInstance(const std::string& strTransDirections) {
	Accumulator* accumulator;

	if (strTransDirections == "x") {
		accumulator = new Accumulator(true, false, false);
		_nNumThermostatedTransDirections = 1;
	} else if (strTransDirections == "y") {
		accumulator = new Accumulator(false, true, false);
		_nNumThermostatedTransDirections = 1;
	} else if (strTransDirections == "z") {
		accumulator = new Accumulator(false, false, true);
		_nNumThermostatedTransDirections = 1;
	} else if (strTransDirections == "xy") {
		accumulator = new Accumulator(true, true, false);
		_nNumThermostatedTransDirections = 2;
	} else if (strTransDirections == "xz") {
		accumulator = new Accumulator(true, false, true);
		_nNumThermostatedTransDirections = 2;
	} else if (strTransDirections == "yz") {
		accumulator = new Accumulator(false, true, true);
		_nNumThermostatedTransDirections = 2;
	} else if (strTransDirections == "xyz") {
		accumulator = new Accumulator(true, true, true);
		_nNumThermostatedTransDirections = 3;
	} else
		accumulator = nullptr;

	return accumulator;
}

void ControlRegionT::readXML(XMLfileUnits& xmlconfig) {
	Domain* domain = global_simulation->getDomain();
	double lc[3];
	double uc[3];
	std::string strVal[3];
	std::string strDirections;

	// coordinates
	xmlconfig.getNodeValue("coords/lcx", lc[0]);
	xmlconfig.getNodeValue("coords/lcy", lc[1]);
	xmlconfig.getNodeValue("coords/lcz", lc[2]);
	xmlconfig.getNodeValue("coords/ucx", strVal[0]);
	xmlconfig.getNodeValue("coords/ucy", strVal[1]);
	xmlconfig.getNodeValue("coords/ucz", strVal[2]);
	// read upper corner
	for (uint8_t d = 0; d < 3; ++d) uc[d] = (strVal[d] == "box") ? domain->getGlobalLength(d) : atof(strVal[d].c_str());

#ifndef NDEBUG
	global_log->info() << "TemperatureControl: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;
#endif

	for (uint8_t d = 0; d < 3; ++d) {
		_dLowerCorner[d] = lc[d];
		_dUpperCorner[d] = uc[d];
	}

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

	// temperature ramp
	_ramp.enabled = false;
	_ramp.update.elapsed = 0;
	bool bRet = true;
	bRet = bRet && xmlconfig.getNodeValue("target/ramp/start", _ramp.start);
	bRet = bRet && xmlconfig.getNodeValue("target/ramp/end", _ramp.end);
	bRet = bRet && xmlconfig.getNodeValue("target/ramp/update/start", _ramp.update.start);
	bRet = bRet && xmlconfig.getNodeValue("target/ramp/update/stop", _ramp.update.stop);
	bRet = bRet && xmlconfig.getNodeValue("target/ramp/update/freq", _ramp.update.freq);
	if (bRet) {
		_ramp.enabled = true;
		global_log->info() << "[TemperatureControl] REGION " << _nStaticID
						   << ": Temperature ramp enabled with start=" << _ramp.start << ", end=" << _ramp.end
						   << ", update.start=" << _ramp.update.start << ", update.stop=" << _ramp.update.stop
						   << ", update.freq=" << _ramp.update.freq << std::endl;
		// calc remaining values
		_ramp.delta = _ramp.end - _ramp.start;
		_ramp.update.delta = _ramp.update.stop - _ramp.update.start;
		_ramp.slope = _ramp.delta / _ramp.update.delta;
	}

	// ControlMethod "VelocityScaling/Andersen/Mixed"
	std::string methods = "";
	xmlconfig.getNodeValue("method", methods);
	if (methods != "") {
		if (methods == "VelocityScaling") {
			_localMethod = VelocityScaling;

			// init data structures
			this->VelocityScalingInit(xmlconfig, strDirections);
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
		this->VelocityScalingInit(xmlconfig, strDirections);
	}

	// measure added kin. energy
	_addedEkin.writeFreq = 1000;
	xmlconfig.getNodeValue("added_ekin/writefreq", _addedEkin.writeFreq);
}

void ControlRegionT::VelocityScalingInit(XMLfileUnits& xmlconfig, std::string strDirections) {
	// settings
	xmlconfig.getNodeValue("settings/numslabs", _nNumSlabs);
	if (_nNumSlabs < 1) {
		global_log->fatal() << "TemperatureControl: need at least one slab! (settings/numslabs)";
		Simulation::exit(932);
	}
	xmlconfig.getNodeValue("settings/exponent", _dTemperatureExponent);
	xmlconfig.getNodeValue("settings/directions", strDirections);
	// calc slab width
	_dSlabWidth = this->GetWidth(1) / ((double)(_nNumSlabs));
	// create accumulator instance
	_accumulator = this->CreateAccumulatorInstance(strDirections);

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
		localThermVars.resize(_nNumSlabs);
	}
	_globalThermVars.resize(_nNumSlabs);

	// init data structure for measure of added kin. energy
	_addedEkin.data.local.resize(_nNumSlabs);
	_addedEkin.data.global.resize(_nNumSlabs);

	_addedEkinLocalThreadBuffer.resize(mardyn_get_max_threads());
	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		addedEkinLocal.resize(_nNumSlabs);
	}

	std::vector<double>& v = _addedEkin.data.local;
	std::fill(v.begin(), v.end(), 0.);

	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		std::fill(addedEkinLocal.begin(), addedEkinLocal.end(), 0.);
	}
}

void ControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp) {
	if (_localMethod != VelocityScaling) return;
	domainDecomp->collCommInit(_nNumSlabs * 4);
	for (unsigned s = 0; s < _nNumSlabs; ++s) {
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
	for (unsigned s = 0; s < _nNumSlabs; ++s) {
		GlobalThermostatVariables& globalTV = _globalThermVars[s];  // do not forget &
		globalTV._numMolecules = domainDecomp->collCommGetUnsLong();
		globalTV._numRotationalDOF = domainDecomp->collCommGetUnsLong();
		globalTV._ekinRot = domainDecomp->collCommGetDouble();
		globalTV._ekinTrans = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	// Adjust target temperature
	uint64_t simstep = _simulation.getSimulationStep();

	if (_ramp.enabled && simstep > _ramp.update.start && simstep <= _ramp.update.stop) {
		if (_ramp.update.elapsed == 0)
			_dTargetTemperature = _ramp.start;
		else if (_ramp.update.stop == simstep)
			_dTargetTemperature = _ramp.end;
		else if (_ramp.update.elapsed % _ramp.update.freq == 0)
			_dTargetTemperature = _ramp.start + _ramp.update.elapsed * _ramp.slope;
		_ramp.update.elapsed++;
	}

	// calc betaTrans, betaRot, and their sum
	double dBetaTransSumSlabs = 0.;
	double dBetaRotSumSlabs = 0.;

	for (unsigned int s = 0; s < _nNumSlabs; ++s) {
		GlobalThermostatVariables& globalTV = _globalThermVars[s];  // do not forget &
		if (globalTV._numMolecules < 1)
			globalTV._betaTrans = 1.;
		else
			globalTV._betaTrans = pow(
				_nNumThermostatedTransDirections * globalTV._numMolecules * _dTargetTemperature / globalTV._ekinTrans,
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

void ControlRegionT::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* /*domainDecomp*/) {
	if (_localMethod != VelocityScaling) return;
	// check componentID
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule inside control region
	for (unsigned short d = 0; d < 3; ++d) {
		double dPos = mol->r(d);

		if (dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d]) return;
	}

	unsigned int nPosIndex;
	unsigned int nIndexMax = _nNumSlabs - 1;

	// calc position index
	std::array<double,3> dLowerCorner = this->GetLowerCorner();
	double dPosRelative = mol->r(1) - dLowerCorner[1];

	nPosIndex = (unsigned int)floor(dPosRelative / _dSlabWidth);

	// ignore outer (halo) molecules
	if (nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
		return;

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

void ControlRegionT::ControlTemperature(Molecule* mol) {
	// check componentID
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule is inside
	for (unsigned short d = 0; d < 3; ++d) {
		double dPos = mol->r(d);

		if (dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d]) return;
	}

	// check for method
	if (_localMethod == VelocityScaling) {
		unsigned int nPosIndex;
		unsigned int nIndexMax = _nNumSlabs - 1;

		// calc position index
		std::array<double,3> dLowerCorner = this->GetLowerCorner();
		double dPosRelative = mol->r(1) - dLowerCorner[1];

		nPosIndex = (unsigned int)floor(dPosRelative / _dSlabWidth);

		// ignore outer (halo) molecules
		if (nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
			return;
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

void ControlRegionT::ResetLocalValues() {
	// reset local values
	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		for (unsigned int s = 0; s < _nNumSlabs; ++s) {
			_localThermVarsThreadBuffer[thread][s].clear();
		}
	}
}

void ControlRegionT::InitBetaLogfile() {
	if (_localMethod == VelocityScaling) {
#ifdef ENABLE_MPI
		DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition());
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank != 0) return;
#endif

		std::stringstream filenamestream;
		filenamestream << _strFilenamePrefixBetaLog << "_reg" << this->GetID() << ".dat";
		std::stringstream outputstream;
		outputstream.write(reinterpret_cast<const char*>(&_nWriteFreqBeta), 8);

		std::ofstream fileout(filenamestream.str().c_str(), std::ios::out | std::ios::binary);
		fileout << outputstream.str();
		fileout.close();
	}
}

void ControlRegionT::WriteBetaLogfile(unsigned long simstep) {
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

	std::stringstream filenamestream;
	filenamestream << _strFilenamePrefixBetaLog << "_reg" << this->GetID() << ".dat";
	std::stringstream outputstream;
	double dInvNumConfigs = 1. / (double)(_numSampledConfigs);
	double dBetaTrans = _dBetaTransSumGlobal * dInvNumConfigs;
	double dBetaRot = _dBetaRotSumGlobal * dInvNumConfigs;
	outputstream.write(reinterpret_cast<const char*>(&dBetaTrans), 8);
	outputstream.write(reinterpret_cast<const char*>(&dBetaRot), 8);

	ofstream fileout(filenamestream.str().c_str(), std::ios::app | std::ios::binary);
	fileout << outputstream.str();
	fileout.close();

	// reset averaged values
	_numSampledConfigs = 0;
	_dBetaTransSumGlobal = 0.;
	_dBetaRotSumGlobal = 0.;
}

DistControl* ControlRegionT::getDistControl() {
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
void ControlRegionT::registerAsObserver() {
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

void ControlRegionT::update(SubjectBase* subject) {
	CuboidRegionObs::update(subject);
	// update slab width
	_dSlabWidth = this->GetWidth(1) / ((double)(_nNumSlabs));
}

void ControlRegionT::writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep) {
	if (_localMethod != VelocityScaling) return;

	if (simstep % _addedEkin.writeFreq != 0) return;

	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		mardyn_assert(_addedEkin.data.local.size() == _nNumSlabs);
		for (size_t slabID = 0; slabID < _nNumSlabs; ++slabID) {
			_addedEkin.data.local[0] += _addedEkinLocalThreadBuffer[thread][slabID];
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
	std::stringstream filenamestream;
	filenamestream << "addedEkin_reg" << this->GetID() << "_cid" << _nTargetComponentID << ".dat";

	std::stringstream outputstream;
	outputstream.write(reinterpret_cast<const char*>(_addedEkin.data.global.data()), 8 * _addedEkin.data.global.size());

	ofstream fileout(filenamestream.str().c_str(), std::ios::app | std::ios::binary);
	fileout << outputstream.str();
	fileout.close();
}

// class TemperatureControl
TemperatureControl::TemperatureControl() {}

TemperatureControl::~TemperatureControl() {
	for (auto region : _vecControlRegions) {
		delete region;
	}
}

void TemperatureControl::readXML(XMLfileUnits& xmlconfig) {
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
		ControlRegionT* region = new ControlRegionT(this);
		region->readXML(xmlconfig);
		this->AddRegion(region);
	}

	bool Vel = false;
	bool And = false;
	// check for mixed mode
	for (auto&& reg : _vecControlRegions) {
		if (reg->_localMethod == ControlRegionT::LocalControlMethod::VelocityScaling)
			Vel = true;
		else if (reg->_localMethod == ControlRegionT::LocalControlMethod::Andersen)
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

void TemperatureControl::prepare_start() {
	for (auto&& reg : _vecControlRegions) {
		reg->registerAsObserver();
	}
	this->InitBetaLogfiles();
}

void TemperatureControl::AddRegion(ControlRegionT* region) { _vecControlRegions.push_back(region); }

void TemperatureControl::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->MeasureKineticEnergy(mol, domainDecomp);
}

void TemperatureControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->CalcGlobalValues(domainDecomp);
}

void TemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) {
		reg->ControlTemperature(mol);
	}
}

void TemperatureControl::Init(unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->ResetLocalValues();
}

void TemperatureControl::InitBetaLogfiles() {
	for (auto&& reg : _vecControlRegions) reg->InitBetaLogfile();
}

void TemperatureControl::WriteBetaLogfiles(unsigned long simstep) {
	for (auto&& reg : _vecControlRegions) reg->WriteBetaLogfile(simstep);
}

void TemperatureControl::writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep) {
	for (auto&& reg : _vecControlRegions) reg->writeAddedEkin(domainDecomp, simstep);
}

/**
 * @brief Decide which ControlMethod to use
 *
 * @param domainDecomposition
 * @param particleContainer
 * @param simstep
 */
void TemperatureControl::DoLoopsOverMolecules(DomainDecompBase* domainDecomposition,
											  ParticleContainer* particleContainer, const unsigned long simstep) {
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
void TemperatureControl::VelocityScalingPreparation(DomainDecompBase* domainDecomposition,
													ParticleContainer* particleContainer, const unsigned long simstep) {
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
