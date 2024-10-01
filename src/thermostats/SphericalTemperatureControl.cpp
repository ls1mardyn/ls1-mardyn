/*
 * SphericalTemperatureControl.cpp
 *
 *  Created on: 13.11.2023
 *      Author: jniemann, based on TemperatureControl plugin by mheinen
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

using Log::global_log;

// init static ID --> instance counting
unsigned short AbstrControlRegionT::_nStaticID = 0;

// class AbstrControlRegionT
AbstrControlRegionT::AbstrControlRegionT(){}
AbstrControlRegionT::~AbstrControlRegionT() { delete _accumulator; }

Accumulator* AbstrControlRegionT::CreateAccumulatorInstance() {	
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::CreateAccumulatorInstance has been called "<< std::endl;
#endif
	Accumulator* accumulator = new Accumulator(true, true, true);
	return accumulator;
}

void AbstrControlRegionT::readXML(XMLfileUnits& xmlconfig) {
    #ifdef DEBUG
	global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::readXML has been called "<< std::endl;
    #endif
	//  individual data (radius, center, corners, ... to be read in readXML-functions of specific ControlRegionT-classes)!

	// target values
	xmlconfig.getNodeValue("target/temperature", _dTargetTemperature);
	xmlconfig.getNodeValue("target/component", _nTargetComponentID);

	// ControlMethod "VelocityScaling/Andersen"
	std::string methods;
	xmlconfig.getNodeValue("method", methods);
	if (methods == "") {
		_localMethod = VelocityScaling;
		global_log->info() << "[TemperatureControl] REGION: no method specified, selecting VelocityScaling"
						<< std::endl;
		// init data structures
		this->VelocityScalingInit(xmlconfig);
	} else if (methods == "VelocityScaling") {
		_localMethod = VelocityScaling;
		// init data structures
		this->VelocityScalingInit(xmlconfig);
		global_log->info() << "[TemperatureControl] REGION 'method' param: " << methods << std::endl;
	} else if (methods == "Andersen") {
		_localMethod = Andersen;
		xmlconfig.getNodeValue("settings/nu", _nuAndersen);
		_timestep = global_simulation->getIntegrator()->getTimestepLength();
		_nuDt = _nuAndersen * _timestep;
		global_log->info() << "[TemperatureControl] REGION 'method' param: " << methods << std::endl;
	} else {
		global_log->error() << "[TemperatureControl] REGION: Invalid 'method' param: " << methods << std::endl;
		Simulation::exit(-1);
	}

	// measure added kin. energy
	_addedEkin.writeFreq = 1000;
	xmlconfig.getNodeValue("added_ekin/writefreq", _addedEkin.writeFreq);
}

void AbstrControlRegionT::VelocityScalingInit(XMLfileUnits& xmlconfig) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::VelocityScalingInit has been called "<< std::endl;
#endif
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

	_addedEkinLocalThreadBuffer.resize(mardyn_get_max_threads());

	_addedEkin.data.local = 0.;
	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		addedEkinLocal = 0.;
	}
}

void AbstrControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::CalcGlobalValues has been called "<< std::endl;
#endif
	if (_localMethod != VelocityScaling) return;
	domainDecomp->collCommInit(4);
	
	unsigned long numMolecules{}, numRotationalDOF{};
	double ekinTrans{}, ekinRot{};
	for (const auto& localVar : _localThermVarsThreadBuffer) {
		numMolecules += localVar._numMolecules;
		numRotationalDOF += localVar._numRotationalDOF;
		ekinTrans += localVar._ekinTrans;
		ekinRot += localVar._ekinRot;
	}

	domainDecomp->collCommAppendUnsLong(numMolecules);
	domainDecomp->collCommAppendUnsLong(numRotationalDOF);
	domainDecomp->collCommAppendDouble(ekinRot);
	domainDecomp->collCommAppendDouble(ekinTrans);
	domainDecomp->collCommAllreduceSum();
	_globalThermVars._numMolecules = domainDecomp->collCommGetUnsLong();
	_globalThermVars._numRotationalDOF = domainDecomp->collCommGetUnsLong();
	_globalThermVars._ekinRot = domainDecomp->collCommGetDouble();
	_globalThermVars._ekinTrans = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

	// Adjust target temperature
	uint64_t simstep = _simulation.getSimulationStep();


	if (_globalThermVars._numMolecules < 1)
		_globalThermVars._betaTrans = 1.;
	else
		_globalThermVars._betaTrans = pow(
			3 * _globalThermVars._numMolecules * _dTargetTemperature / _globalThermVars._ekinTrans,
			_dTemperatureExponent);

	if (_globalThermVars._numRotationalDOF < 1)
		_globalThermVars._betaRot = 1.;
	else
		_globalThermVars._betaRot = pow(_globalThermVars._numRotationalDOF * _dTargetTemperature / _globalThermVars._ekinRot, _dTemperatureExponent);

	// calc ensemble average of beta_trans, beta_rot
	_dBetaTransSumGlobal += _globalThermVars._betaTrans;
	_dBetaRotSumGlobal += _globalThermVars._betaRot;
	_numSampledConfigs++;
}

void AbstrControlRegionT::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* /*domainDecomp*/) {
	if (_localMethod != VelocityScaling) return;
	// check componentID
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
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
	if(!this->ContainsMolecule(mol)) return;

	auto myThreadNum = mardyn_get_thread_num();
	LocalThermostatVariables& localTV = _localThermVarsThreadBuffer[myThreadNum];  // do not forget &  
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

void AbstrControlRegionT::ControlTemperature(Molecule* mol) {
	if (mol->componentid() + 1 != _nTargetComponentID &&
		0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	if(!this->ContainsMolecule(mol)) return;

	// check for method
	if (_localMethod == VelocityScaling) {
		
		GlobalThermostatVariables& globalTV = _globalThermVars;  // do not forget &
		if (_globalThermVars._numMolecules < 1) return;
		// scale velocity
		double vcorr = 2. - 1. / _globalThermVars._betaTrans;
		double Dcorr = 2. - 1. / _globalThermVars._betaRot;

		// measure added kin. energy
		double v2_old = mol->v2();

		// if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: calling  _accumulator->ScaleVelocityComponents(mol, vcorr); "<< std::endl;
		// if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: mol =  "<< mol->getID()<<", vcorr = "<< vcorr << " " << std::endl;
		_accumulator->ScaleVelocityComponents(mol, vcorr);

		// measure added kin. energy
		double v2_new = mol->v2();
		int mythread = mardyn_get_thread_num();
		_addedEkinLocalThreadBuffer[mythread] += (v2_new - v2_old);

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

void AbstrControlRegionT::ResetLocalValues() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::ResetLocalValues has been called "<< std::endl;
#endif
	// reset local values
	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		_localThermVarsThreadBuffer[thread].clear();
	}
}

void AbstrControlRegionT::InitBetaLogfile() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::InitBetaLogfile has been called "<< std::endl;
#endif
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
		ofs << std::setw(12) << "simstep"
			<< std::setw(24) << "dBetaTrans"
			<< std::setw(24) << "dBetaRot"
			<< std::endl;
		ofs.close();
	}
}

void AbstrControlRegionT::WriteBetaLogfile(unsigned long simstep) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::WriteBetaLogfile has been called "<< std::endl;
#endif
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

	double dBetaTrans = _dBetaTransSumGlobal  / (double)(_numSampledConfigs);
	double dBetaRot = _dBetaRotSumGlobal  / (double)(_numSampledConfigs);

	// writing to file
	const std::string fname = _strFilenamePrefixBetaLog + "_reg" + std::to_string(this->GetID()) + ".dat";
	std::ofstream ofs;
	ofs.open(fname, std::ios::app);
	ofs << std::setw(12) << simstep
		<< FORMAT_SCI_MAX_DIGITS << dBetaTrans
		<< FORMAT_SCI_MAX_DIGITS << dBetaRot
		<< std::endl;
	ofs.close();

	// reset averaged values
	_numSampledConfigs = 0;
	_dBetaTransSumGlobal = 0.;
	_dBetaRotSumGlobal = 0.;
}

DistControl* AbstrControlRegionT::getDistControl() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::getDistControl has been called "<< std::endl;
#endif
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

// void AbstrControlRegionT::update(SubjectBase* subject) {
// 	if(_bDebugOutput) global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::update  has been called "<< std::endl;
// 	SphericalRegionObs::update(subject);
// }

void AbstrControlRegionT::InitAddedEkin() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::InitAddedEkin has been called "<< std::endl;
#endif
	if (_localMethod == VelocityScaling) {
#ifdef ENABLE_MPI
		DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition());
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank != 0) return;
#endif
		// touch file
		const std::string fname = "addedEkin_reg" + std::to_string(this->GetID()) + "_cid" + std::to_string(_nTargetComponentID) + ".dat";
		std::ofstream ofs(fname, std::ios::out);
		ofs << std::setw(12) << "simstep";
		ofs << std::setw(24) << "regionData";
		ofs << std::endl;
		ofs.close();
	}
}

void AbstrControlRegionT::writeAddedEkin(DomainDecompBase* domainDecomp, const uint64_t& simstep) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: AbstrControlRegionT::writeAddedEkin has been called "<< std::endl;
#endif
	if (_localMethod != VelocityScaling) return;

	if (simstep % _addedEkin.writeFreq != 0) return;

	for (int thread = 0; thread < mardyn_get_max_threads(); ++thread) {
		_addedEkin.data.local += _addedEkinLocalThreadBuffer[thread] ;
	}
	// calc global values
#ifdef ENABLE_MPI
	MPI_Reduce(_addedEkin.data.local.data(), _addedEkin.data.global.data(), _addedEkin.data.local.size(), MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);
#else
	_addedEkin.data.global = _addedEkin.data.local;
#endif

	// reset local values
	_addedEkin.data.local = 0.;
	for (auto& addedEkinLocal : _addedEkinLocalThreadBuffer) {
		addedEkinLocal = 0.;
	}

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0) return;
#endif

	// dekin = d2ekin * 0.5
	_addedEkin.data.global *= 0.5;

	// writing .dat-files
	const std::string fname = "addedEkin_reg" + std::to_string(this->GetID()) + "_cid" + std::to_string(_nTargetComponentID) + ".dat";
	std::ofstream ofs;
	ofs.open(fname, std::ios::app);
	
	ofs << std::setw(12) << simstep;
	ofs << FORMAT_SCI_MAX_DIGITS << _addedEkin.data.global;
	ofs << std::endl;
	ofs.close();
}






// class SphericalControlRegionT
SphericalControlRegionT::SphericalControlRegionT(SphericalTemperatureControl* const parent)
	: AbstrControlRegionT(), SphericalRegionObs(parent),
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

void SphericalControlRegionT::readXML(XMLfileUnits& xmlconfig) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: readXML has been called "<< std::endl;
#endif

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

#ifdef DEBUG
    global_log->info() << "SphericalTemperatureControl: Center:" << ctr[0] << ", " << ctr[1] << ", " << ctr[2] << ", radius: "<< rad<< std::endl;
#endif

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


	AbstrControlRegionT::readXML(xmlconfig);
}


bool SphericalControlRegionT::ContainsMolecule(Molecule* mol) {
	// check if molecule inside control region
	double distanceFromCenterSquared = 0;
	for(int d = 0; d<3; d++){
		distanceFromCenterSquared += std::pow(mol->r(d)-_dCenter[d],2);
	}
	if(distanceFromCenterSquared > _dRadiusSquared) return false;
	else return true;
}


void SphericalControlRegionT::registerAsObserver() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::registerAsObserver has been called "<< std::endl;
#endif
	
	if (true == _bIsObserver) {
		DistControl* distControl = this->getDistControl();
		if (distControl != nullptr)
			distControl->registerObserver(this);
		else {
			global_log->error() << "TemperatureControl->region[" << this->GetID()
								<< "]: Initialization of plugin DistControl is needed before! Program exit..." << std::endl;
			Simulation::exit(-1);
		}
	}
}

void SphericalControlRegionT::update(SubjectBase* subject) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::update  has been called "<< std::endl;
#endif
	SphericalRegionObs::update(subject);
}
// class SphericalControlRegionT





// class SphereComplementControlRegionT
SphereComplementControlRegionT::SphereComplementControlRegionT(SphericalTemperatureControl* const parent)
	: AbstrControlRegionT(), SphereComplementRegionObs(parent), 
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
	SphericalRegion::_nID = ++_nStaticID;
}

SphereComplementControlRegionT::~SphereComplementControlRegionT() { delete _accumulator; }


void SphereComplementControlRegionT::readXML(XMLfileUnits& xmlconfig) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphereComplementControlRegionT::readXML has been called "<< std::endl;
#endif

	Domain* domain = global_simulation->getDomain();
	double lc[3] = {0.,0.,0.}; //default
	double uc[3];
	double ctr[3];
	double rad;
	std::string strVal[3] = {"box", "box", "box"}; //default

	// coordinates
	xmlconfig.getNodeValue("coords/lcx", lc[0]);
	xmlconfig.getNodeValue("coords/lcy", lc[1]);
	xmlconfig.getNodeValue("coords/lcz", lc[2]);
	xmlconfig.getNodeValue("coords/ucx", strVal[0]);
	xmlconfig.getNodeValue("coords/ucy", strVal[1]);
	xmlconfig.getNodeValue("coords/ucz", strVal[2]);
	xmlconfig.getNodeValue("coords/ctrx", ctr[0]);
	xmlconfig.getNodeValue("coords/ctry", ctr[1]);
	xmlconfig.getNodeValue("coords/ctrz", ctr[2]);
	xmlconfig.getNodeValue("coords/radius", rad);
	// read upper corner
	for (uint8_t d = 0; d < 3; ++d) uc[d] = (strVal[d] == "box") ? domain->getGlobalLength(d) : atof(strVal[d].c_str());
#ifndef NDEBUG
	global_log->info() << "TemperatureControl: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << std::endl;
#endif
#ifdef DEBUG
    global_log->info() << "SphericalTemperatureControl: Center:" << ctr[0] << ", " << ctr[1] << ", " << ctr[2] << ", radius: "<< rad<< std::endl;
#endif

	for (uint8_t d = 0; d < 3; ++d) {
		_dLowerCorner[d] = lc[d];
		_dUpperCorner[d] = uc[d];
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


	AbstrControlRegionT::readXML(xmlconfig);
}


bool SphereComplementControlRegionT::ContainsMolecule(Molecule* mol){
	// check if molecule inside control region
	double distanceFromCenterSquared = 0;
	for (unsigned short d = 0; d < 3; ++d) {
		double dPos = mol->r(d);
		if (dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d]) return false;
		distanceFromCenterSquared += std::pow(dPos-_dCenter[d],2);
	}
	if(distanceFromCenterSquared < _dRadiusSquared) return false;
	else return true;
}



void SphereComplementControlRegionT::registerAsObserver() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::registerAsObserver has been called "<< std::endl;
#endif
	
	if (true == _bIsObserver) {
		DistControl* distControl = this->getDistControl();
		if (distControl != nullptr)
			distControl->registerObserver(this);
		else {
			global_log->error() << "TemperatureControl->region[" << this->GetID()
								<< "]: Initialization of plugin DistControl is needed before! Program exit..." << std::endl;
			Simulation::exit(-1);
		}
	}
}

void SphereComplementControlRegionT::update(SubjectBase* subject) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalControlRegionT::update  has been called "<< std::endl;
#endif
	SphereComplementRegionObs::update(subject);
}

// class SphereComplementControlRegionT









// class SphericalTemperatureControl
SphericalTemperatureControl::SphericalTemperatureControl() {}

SphericalTemperatureControl::~SphericalTemperatureControl() {
	for (auto region : _vecControlRegions) {
		delete region;
	}
}

void SphericalTemperatureControl::readXML(XMLfileUnits& xmlconfig) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::readXML has been called "<< std::endl;
#endif
	// control
	xmlconfig.getNodeValue("control/start", _nStart);
	xmlconfig.getNodeValue("control/frequency", _nControlFreq);
	xmlconfig.getNodeValue("control/stop", _nStop);
	global_log->info() << "Start control from simstep: " << _nStart << std::endl;
	global_log->info() << "Control with frequency: " << _nControlFreq << std::endl;
	global_log->info() << "Stop control at simstep: " << _nStop << std::endl;

	// turn on/off explosion heuristics
	// domain->setExplosionHeuristics(bUseExplosionHeuristics);

	// add regions
	uint32_t numRegions = 0;
	XMLfile::Query query = xmlconfig.query("regions/region");
	numRegions = query.card();
	global_log->info() << "Number of control regions: " << numRegions << std::endl;
	if (numRegions < 1) {
		global_log->warning() << "No region parameters specified." << std::endl;
	}
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputRegionIter;
	for (outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++) {
		xmlconfig.changecurrentnode(outputRegionIter);
		SphericalControlRegionT* sphericalRegion = new SphericalControlRegionT(this);
		SphereComplementControlRegionT* complementRegion = new SphereComplementControlRegionT(this);
		sphericalRegion->readXML(xmlconfig);
		complementRegion->readXML(xmlconfig);
		this->AddRegion(sphericalRegion);
		this->AddRegion(complementRegion);
	}

	bool Vel = false;
	bool And = false;
	// check for mixed mode
	for (auto&& reg : _vecControlRegions) {
		if (reg->_localMethod == AbstrControlRegionT::LocalControlMethod::VelocityScaling)
			Vel = true;
		else if (reg->_localMethod == AbstrControlRegionT::LocalControlMethod::Andersen)
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
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::prepare_start has been called "<< std::endl;
#endif
	for (auto&& reg : _vecControlRegions) {
		reg->registerAsObserver();
	}
	this->InitBetaLogfiles();
	this->InitAddedEkin();
}

void SphericalTemperatureControl::AddRegion(AbstrControlRegionT* region) { _vecControlRegions.push_back(region); }

void SphericalTemperatureControl::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->MeasureKineticEnergy(mol, domainDecomp);
}

void SphericalTemperatureControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::CalcGlobalValues has been called "<< std::endl;
#endif
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->CalcGlobalValues(domainDecomp);
}

void SphericalTemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::ControlTemperature has been called "<< std::endl;
#endif
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) {
		reg->ControlTemperature(mol);
	}
	
}

void SphericalTemperatureControl::Init(unsigned long simstep) {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::Init has been called "<< std::endl;
#endif
	if (simstep % _nControlFreq != 0) return;
	if (simstep <= this->GetStart() || simstep > this->GetStop()) return;

	for (auto&& reg : _vecControlRegions) reg->ResetLocalValues();
}

void SphericalTemperatureControl::InitBetaLogfiles() {
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::InitBetaLogfiles has been called "<< std::endl;
#endif
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
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::DoLoopsOverMolecules has been called "<< std::endl;
#endif
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
#ifdef DEBUG
    global_log->info() << "[SphericalTemperatureControl]: SphericalTemperatureControl::VelocityScalingPreparation has been called "<< std::endl;
#endif
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