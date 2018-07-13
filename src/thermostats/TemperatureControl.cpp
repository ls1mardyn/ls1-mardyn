/*
 * TemperatureControl.cpp
 *
 *  Created on: 27.05.2015
 *      Author: mheinen
 */

#include "thermostats/TemperatureControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "utils/xmlfileUnits.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>

using namespace std;

// init static ID --> instance counting
unsigned short ControlRegionT::_nStaticID = 0;

// class ControlRegionT

ControlRegionT::ControlRegionT()
		:
		_nID(0),
		_dLowerCorner{0.,0.,0.},
		_dUpperCorner{0.,0.,0.},
		_nNumSlabs(0),
		_dSlabWidth(0.0),
		_thermVars(),
		_dTargetTemperature(0.0),
		_dTemperatureExponent(0.0),
		_nTargetComponentID(0),
		_nNumThermostatedTransDirections(0),
		_nRegionID(0),
		_accumulator(nullptr),
		_strFilenamePrefixBetaLog("beta_log"),
		_nWriteFreqBeta(1000),
		_numSampledConfigs(0),
		_dBetaTransSumGlobal(0.0),
		_dBetaRotSumGlobal(0.0)
{
	// ID
	_nID = ++_nStaticID;
}


ControlRegionT::~ControlRegionT()
{
}

Accumulator* ControlRegionT::CreateAccumulatorInstance(std::string strTransDirections)
{
	Accumulator* accumulator;

	if(strTransDirections == "x")
	{
		accumulator = new Accumulator(true, false, false);
		_nNumThermostatedTransDirections = 1;
	}
	else if(strTransDirections == "y")
	{
		accumulator = new Accumulator(false, true, false);
		_nNumThermostatedTransDirections = 1;
	}
	else if(strTransDirections == "z")
	{
		accumulator = new Accumulator(false, false, true);
		_nNumThermostatedTransDirections = 1;
	}
	else if(strTransDirections == "xy")
	{
		accumulator = new Accumulator(true, true, false);
		_nNumThermostatedTransDirections = 2;
	}
	else if(strTransDirections == "xz")
	{
		accumulator = new Accumulator(true, false, true);
		_nNumThermostatedTransDirections = 2;
	}
	else if(strTransDirections == "yz")
	{
		accumulator = new Accumulator(false, true, true);
		_nNumThermostatedTransDirections = 2;
	}
	else if(strTransDirections == "xyz")
	{
		accumulator = new Accumulator(true, true, true);
		_nNumThermostatedTransDirections = 3;
	}
	else
		accumulator = NULL;

	return accumulator;
}

void ControlRegionT::readXML(XMLfileUnits& xmlconfig)
{
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
	for(uint8_t d=0; d<3; ++d)
		uc[d] = (strVal[d] == "box") ? domain->getGlobalLength(d) : atof(strVal[d].c_str() );

#ifndef NDEBUG
	global_log->info() << "TemperatureControl: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;
#endif

	for(uint8_t d=0; d<3; ++d)
	{
		_dLowerCorner[d] = lc[d];
		_dUpperCorner[d] = uc[d];
	}

	// target values
	xmlconfig.getNodeValue("target/temperature", _dTargetTemperature);
	xmlconfig.getNodeValue("target/component", _nTargetComponentID);

	// ControlMethod "VelocityScaling/Andersen/Mixed"
	std::string methods = "";
	xmlconfig.getNodeValue("method", methods);
	if(methods != "") {
		if(methods == "VelocityScaling") {
			_localMethod = VelocityScaling;

			// init data structures
			this->VelocityScalingInit(xmlconfig, strDirections);
		}
		else if(methods == "Andersen") {
			_localMethod = Andersen;
			xmlconfig.getNodeValue("settings/nu", _nuAndersen);
			_timestep = global_simulation->getIntegrator()->getTimestepLength();
			_nuDt = _nuAndersen*_timestep;
		}
		else {
			global_log -> error() << "[TemperatureControl] REGION: Invalid 'method' param: " << methods << std::endl;
			Simulation::exit(-1);
		}
		global_log -> info() << "[TemperatureControl] REGION 'method' param: " << methods << std::endl;
	}
	//
	else {
		_localMethod = VelocityScaling;
		global_log -> info() << "[TemperatureControl] REGION: no method specified, selecting VelocityScaling" << std::endl;

		// init data structures
		this->VelocityScalingInit(xmlconfig, strDirections);
	}
}

void ControlRegionT::VelocityScalingInit(XMLfileUnits &xmlconfig, std::string strDirections)
{
	// settings
	xmlconfig.getNodeValue("settings/numslabs", _nNumSlabs);
	xmlconfig.getNodeValue("settings/exponent", _dTemperatureExponent);
	xmlconfig.getNodeValue("settings/directions", strDirections);
	// calc slab width
	_dSlabWidth = this->GetWidth(1) / ( (double)(_nNumSlabs) );
	// create accumulator instance
	_accumulator = this->CreateAccumulatorInstance(strDirections);

	// write control for beta_trans and beta_rot log file
	_nWriteFreqBeta = 1000;
	_strFilenamePrefixBetaLog = "beta_log";
	xmlconfig.getNodeValue("writefreq",  _nWriteFreqBeta);
	xmlconfig.getNodeValue("fileprefix", _strFilenamePrefixBetaLog);
	if(_nWriteFreqBeta==0) {
		global_log->warning() << "Temperature Control: write Frequency was specified to be zero. This is NOT allowed. Reset it to 1000." << std::endl;
		_nWriteFreqBeta = 1000;
	}
	this->InitBetaLogfile();
	_thermVars.resize(_nNumSlabs);
}

void ControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp )
{
	if(_localMethod != VelocityScaling)
		return;
	domainDecomp->collCommInit(_nNumSlabs * 4);
	for (unsigned s = 0; s < _nNumSlabs; ++s) {
		LocalThermostatVariables & localTV = _thermVars[s]._local; // do not forget &
		domainDecomp->collCommAppendUnsLong(localTV._numMolecules);
		domainDecomp->collCommAppendUnsLong(localTV._numRotationalDOF);
		domainDecomp->collCommAppendDouble(localTV._ekinRot);
		domainDecomp->collCommAppendDouble(localTV._ekinTrans);
	}
	domainDecomp->collCommAllreduceSum();
	for (unsigned s = 0; s < _nNumSlabs; ++s) {
		GlobalThermostatVariables & globalTV = _thermVars[s]._global;  // do not forget &
		globalTV._numMolecules = domainDecomp->collCommGetUnsLong();
		globalTV._numRotationalDOF = domainDecomp->collCommGetUnsLong();
		globalTV._ekinRot = domainDecomp->collCommGetDouble();
		globalTV._ekinTrans = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	// calc betaTrans, betaRot, and their sum
	double dBetaTransSumSlabs = 0.;
	double dBetaRotSumSlabs = 0.;

	for(unsigned int s = 0; s<_nNumSlabs; ++s)
	{
		GlobalThermostatVariables & globalTV = _thermVars[s]._global;  // do not forget &
		if( globalTV._numMolecules < 1 )
			globalTV._betaTrans = 1.;
		else
			globalTV._betaTrans = pow(_nNumThermostatedTransDirections * globalTV._numMolecules * _dTargetTemperature / globalTV._ekinTrans, _dTemperatureExponent);

		if( globalTV._numRotationalDOF < 1 )
			globalTV._betaRot = 1.;
		else
			globalTV._betaRot = pow( globalTV._numRotationalDOF * _dTargetTemperature / globalTV._ekinRot, _dTemperatureExponent);

		// calc sums over all slabs
		dBetaTransSumSlabs += globalTV._betaTrans;
		dBetaRotSumSlabs   += globalTV._betaRot;
	}
	// calc ensemble average of beta_trans, beta_rot
	_dBetaTransSumGlobal += dBetaTransSumSlabs;
	_dBetaRotSumGlobal   += dBetaRotSumSlabs;
	_numSampledConfigs++;

//    cout << "_nNumMoleculesGlobal = " << _nNumMoleculesGlobal << endl;
//    cout << "_dBetaTransGlobal = " << _dBetaTransGlobal << endl;
//    cout << "_dTargetTemperature = " << _dTargetTemperature << endl;
//    cout << "_d2EkinRotGlobal = " << _d2EkinRotGlobal << endl;
//
//    cout << "_nRotDOFGlobal = " << _nRotDOFGlobal << endl;
//    cout << "_dBetaRotGlobal = " << _dBetaRotGlobal << endl;
//    cout << "_d2EkinRotGlobal = " << _d2EkinRotGlobal << endl;

}

void ControlRegionT::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* /*domainDecomp*/)
{
	if(_localMethod != VelocityScaling)
		return;
	// check componentID
	if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule inside control region
	for(unsigned short d = 0; d<3; ++d)
	{
		double dPos = mol->r(d);

		if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
			return;
	}

	unsigned int nPosIndex;
	unsigned int nIndexMax = _nNumSlabs - 1;

	// calc position index
	double* dLowerCorner = this->GetLowerCorner();
	double dPosRelative = mol->r(1) - dLowerCorner[1];

	nPosIndex = (unsigned int) floor(dPosRelative / _dSlabWidth);

	// ignore outer (halo) molecules
	if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
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

	LocalThermostatVariables & localTV = _thermVars[nPosIndex]._local;  // do not forget &
	localTV._ekinTrans += _accumulator->CalcKineticEnergyContribution(mol);

	// sum up rot. kinetic energy (2x)
	double dDummy = 0.;

	mol->calculate_mv2_Iw2(dDummy, localTV._ekinRot );

	// count num molecules
	localTV._numMolecules++;

	// count rotational DOF
	localTV._numRotationalDOF += mol->component()->getRotationalDegreesOfFreedom();
}


void ControlRegionT::ControlTemperature(Molecule* mol)
{
	// check componentID
	if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

	// check if molecule is inside
	for(unsigned short d = 0; d<3; ++d)
	{
		double dPos = mol->r(d);

		if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
			return;
	}

	// check for method
	if(_localMethod == VelocityScaling) {
		unsigned int nPosIndex;
		unsigned int nIndexMax = _nNumSlabs - 1;

		// calc position index
		double *dLowerCorner = this->GetLowerCorner();
		double dPosRelative = mol->r(1) - dLowerCorner[1];

		nPosIndex = (unsigned int) floor(dPosRelative / _dSlabWidth);

		// ignore outer (halo) molecules
		if (nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
			return;

		GlobalThermostatVariables &globalTV = _thermVars[nPosIndex]._global;  // do not forget &
		if (globalTV._numMolecules < 1)
			return;

		// scale velocity
		double vcorr = 2. - 1. / globalTV._betaTrans;
		double Dcorr = 2. - 1. / globalTV._betaRot;

/*
	mol->setv(0, mol->v(0) * vcorr);
//    mol->setv(1, mol->v(1) * vcorr);
	mol->setv(2, mol->v(2) * vcorr);
*/

		_accumulator->ScaleVelocityComponents(mol, vcorr);

		mol->scale_D(Dcorr);
	}
	else if(_localMethod == Andersen) {
		double stdDevTrans, stdDevRot;
		if(_rand.rnd() < _nuDt) {
			stdDevTrans = sqrt(_dTargetTemperature/mol->mass());
			for(unsigned short d = 0; d < 3; d++) {
				stdDevRot = sqrt(_dTargetTemperature*mol->getI(d));
				mol->setv(d, _rand.gaussDeviate(stdDevTrans));
				mol->setD(d, _rand.gaussDeviate(stdDevRot));
			}
		}
	}
	else {
		global_log -> error() << "[TemperatureControl] Invalid localMethod param: " << _localMethod << std::endl;
		Simulation::exit(-1);
	}
}

void ControlRegionT::ResetLocalValues()
{
	// reset local values
	for(unsigned int s = 0; s<_nNumSlabs; ++s)
	{
		_thermVars[s]._local.clear();
	}
}

void ControlRegionT::InitBetaLogfile()
{
    if(_localMethod == VelocityScaling) {
        DomainDecompBase *domainDecomp = &(global_simulation->domainDecomposition());

#ifdef ENABLE_MPI
        int rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank!= 0)
            return;
#endif

        std::stringstream filenamestream;
        filenamestream << _strFilenamePrefixBetaLog << "_reg" << this->GetID() << ".dat";
        std::stringstream outputstream;
        outputstream.write(reinterpret_cast<const char *>(&_nWriteFreqBeta), 8);

        ofstream fileout(filenamestream.str().c_str(), std::ios::out | std::ios::binary);
        fileout << outputstream.str();
        fileout.close();
    }
}

void ControlRegionT::WriteBetaLogfile(unsigned long simstep)
{
    if(_localMethod != VelocityScaling)
        return;
	if(0 != (simstep % _nWriteFreqBeta) )
		return;

	DomainDecompBase* domainDecomp = &(global_simulation->domainDecomposition() );

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank!= 0)
		return;
#endif

	std::stringstream filenamestream;
	filenamestream << _strFilenamePrefixBetaLog << "_reg" << this->GetID() << ".dat";
	std::stringstream outputstream;
	double dInvNumConfigs = 1. / (double)(_numSampledConfigs);
	double dBetaTrans = _dBetaTransSumGlobal * dInvNumConfigs;
	double dBetaRot   = _dBetaRotSumGlobal   * dInvNumConfigs;
	outputstream.write(reinterpret_cast<const char*>(&dBetaTrans), 8);
	outputstream.write(reinterpret_cast<const char*>(&dBetaRot),   8);

	ofstream fileout(filenamestream.str().c_str(), std::ios::app | std::ios::binary);
	fileout << outputstream.str();
	fileout.close();

	// reset averaged values
	_numSampledConfigs = 0;
	_dBetaTransSumGlobal = 0.;
	_dBetaRotSumGlobal = 0.;
}

// class TemperatureControl
TemperatureControl::TemperatureControl()
{
}

TemperatureControl::~TemperatureControl()
{

}

void TemperatureControl::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/start", _nStart);
	xmlconfig.getNodeValue("control/frequency", _nControlFreq);
	xmlconfig.getNodeValue("control/stop", _nStop);
	global_log->info() << "Start control from simstep: " << _nStart << endl;
	global_log->info() << "Control with frequency: " << _nControlFreq << endl;
	global_log->info() << "Stop control at simstep: " << _nStop << endl;

	// turn on/off explosion heuristics
	// domain->SetExplosionHeuristics(bUseExplosionHeuristics);

	// add regions
	uint32_t numRegions = 0;
	XMLfile::Query query = xmlconfig.query("regions/region");
	numRegions = query.card();
	global_log->info() << "Number of control regions: " << numRegions << endl;
	if(numRegions < 1) {
		global_log->warning() << "No region parameters specified." << endl;
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputRegionIter;
	for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
	{
		xmlconfig.changecurrentnode(outputRegionIter);
		ControlRegionT* region = new ControlRegionT();
		region->readXML(xmlconfig);
		this->AddRegion(region);
	}

	bool Vel = false;
	bool And = false;
	// check for mixed mode
	for(auto&& reg : _vecControlRegions){
		if(reg->_localMethod == ControlRegionT::LocalControlMethod::VelocityScaling)
			Vel = true;
		else if(reg->_localMethod == ControlRegionT::LocalControlMethod::Andersen)
			And = true;
	}
	if(Vel && And){
		_method = Mixed;
		global_log -> info() << "[TemperatureControl] Mixed methods across regions\n";
	}
	else if(!Vel && And){
		_method = Andersen;
		global_log -> info() << "[TemperatureControl] Andersen in all regions\n";
	}
	else{
		_method = VelocityScaling;
		global_log -> info() << "[TemperatureControl] VelocityControl in all regions\n";
	}
}

void TemperatureControl::AddRegion(ControlRegionT* region)
{
	_vecControlRegions.push_back(region);
}

void TemperatureControl::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if(simstep % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->MeasureKineticEnergy(mol, domainDecomp);
}

void TemperatureControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if(simstep % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->CalcGlobalValues(domainDecomp);
}


void TemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep)
{
	if(simstep % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->ControlTemperature(mol);
}

void TemperatureControl::Init(unsigned long simstep)
{
	if(simstep % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->ResetLocalValues();
}

void TemperatureControl::InitBetaLogfiles()
{
	for(auto&& reg : _vecControlRegions)
		reg->InitBetaLogfile();
}

void TemperatureControl::WriteBetaLogfiles(unsigned long simstep)
{
	for(auto&& reg : _vecControlRegions)
		reg->WriteBetaLogfile(simstep);
}

/**
 * @brief Decide which ControlMethod to use
 *
 * @param domainDecomposition
 * @param particleContainer
 * @param simstep
 */
void TemperatureControl::DoLoopsOverMolecules(DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned long simstep)
{
	if(_method == VelocityScaling || _method == Mixed){
		this->VelocityScalingPreparation(domainDecomposition, particleContainer, simstep);
		global_log->debug() << "[TemperatureControl] VelocityScalingPreparation\n";
	}

	// iterate over all molecules. ControlTemperature depends on _localMethod for Region molecule is in
	ParticleIterator tM;
	for( tM  = particleContainer->iterator();
		 tM.hasNext();
		 tM.next())
	{
		// control temperature
		this->ControlTemperature(&(*tM), simstep);
	}
}

/**
 * @brief Prepare for VelocityScaling control method
 *
 * @param domainDecomposition
 * @param particleContainer
 * @param simstep
 */
void TemperatureControl::VelocityScalingPreparation(DomainDecompBase *domainDecomposition,
													ParticleContainer *particleContainer, unsigned long simstep) {
	// respect start/stop
	if(this->GetStart() <= simstep && this->GetStop() > simstep)
	{
//		global_log->info() << "Thermostat ON!" << endl;

		ParticleIterator tM;

		// init temperature control
		this->Init(simstep);

		for( tM  = particleContainer->iterator();
			 tM.hasNext();
			 tM.next())
		{
			// measure kinetic energy
			this->MeasureKineticEnergy(&(*tM), domainDecomposition, simstep);

//          cout << "id = " << tM->id() << ", (vx,vy,vz) = " << tM->v(0) << ", " << tM->v(1) << ", " << tM->v(2) << endl;
		}

		// calc global values
		this->CalcGlobalValues(domainDecomposition, simstep);

		// write beta_trans, beta_rot log-files
		this->WriteBetaLogfiles(simstep);
	}
}



