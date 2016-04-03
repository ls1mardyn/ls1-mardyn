#include "Leapfrog.h"

#include <map>

#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/PressureGradient.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"


using namespace std;
using Log::global_log;

Leapfrog::Leapfrog(double timestepLength) :	Integrator(timestepLength) {
	init();
}

void Leapfrog::init() {
	// set starting state
	_state = STATE_POST_FORCE_CALCULATION;
}

Leapfrog::~Leapfrog() {}

void Leapfrog::readXML(XMLfileUnits& xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << endl;
	assert(_timestepLength > 0);
}

void Leapfrog::eventForcesCalculated(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		transition2to3(molCont, domain);
	}
}

void Leapfrog::eventNewTimestep(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		transition3to1(molCont, domain);
		transition1to2(molCont, domain);
	}
}

void Leapfrog::transition1to2(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_NEW_TIMESTEP) {
		Molecule* tempMolecule;
		double vcorr = 2. - 1. / domain->getGlobalBetaTrans();
		double Dcorr = 2. - 1. / domain->getGlobalBetaRot();
		unsigned molID;
		unsigned int cid;
		int thermostat;
		string moved ("moved");
		string fixed ("fixed");
		unsigned cid_moved =  domain->getPG()->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
		unsigned cid_fixed =  domain->getPG()->getCidMovement(fixed, domain->getNumberOfComponents()) - 1;

		
		for (tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()) {
		  molID = tempMolecule->id();
		  cid = tempMolecule->componentid();
		  thermostat = domain->getThermostat(cid);
		  
		  if((molID%250 == 0 && cid == cid_fixed) || (_simulation.getSimulationStep() <= _simulation.getInitStatistics() && (molID%500 == 0 && cid == cid_moved))){}	// each 250th Molecules in Component CID=3 and each 500th in CID=2 (for t<=_initStatistics) are fixed to allow Temperatures T > 0K!
		  else {					
		    if(domain->isScaling1Dim(thermostat) && domain->getAlphaTransCorrection(thermostat) == false){
			vcorr = 2. - 1. / domain->getGlobalAlphaTrans();
		    }
		    tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr, domain);
		  }
		}

		this->_state = STATE_PRE_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition1to2(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition2to3(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		Molecule* tM;
		map<int, unsigned long> N;
		map<int, unsigned long> rotDOF;
		map<int, double> summv2;
		map<int, double> sumIw2;
		map<int, double> summv2_1Dim;
		int molID;
		unsigned int cid;
		
		string moved ("moved");
		string fixed ("fixed");
		unsigned cid_moved =  domain->getPG()->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
		unsigned cid_fixed =  domain->getPG()->getCidMovement(fixed, domain->getNumberOfComponents()) - 1;

		double dt_half = 0.5 * this->_timestepLength;
		if (domain->severalThermostats()) {
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				molID = tM->id();
				cid = tM->componentid();
				if((molID%250 == 0 && cid == cid_fixed) || (_simulation.getSimulationStep() <= _simulation.getInitStatistics() && (molID%500 == 0 && cid == cid_moved))){}	// each 250th Molecules in Component CID=3 and each 500th in CID=2 (for t<=_initStatistics) are fixed to allow Temperatures T > 0K!
				else {
				  int thermostat = domain->getThermostat(cid);
				  tM->upd_postF(dt_half, summv2[thermostat], summv2_1Dim[thermostat], sumIw2[thermostat], domain);
				  
				  N[thermostat]++;
				  rotDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
				}
			}
		}
		else {
			unsigned long Ngt = 0;
			unsigned long rotDOFgt = 0;
			double summv2gt = 0.0;
			double sumIw2gt = 0.0;
			double summv2_1Dimgt = 0.0;
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				tM->upd_postF(dt_half, summv2gt, summv2_1Dimgt, sumIw2gt, domain);
				assert(summv2gt >= 0.0);
				Ngt++;
				rotDOFgt += tM->component()->getRotationalDegreesOfFreedom();
			}
			N[0] = Ngt;
			rotDOF[0] = rotDOFgt;
			summv2[0] = summv2gt;
			sumIw2[0] = sumIw2gt;
			summv2_1Dim[0] = summv2_1Dimgt;
		} 
		for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
			domain->setLocalSummv2(thermit->second, thermit->first);
			domain->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
			domain->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
			if(domain->isScaling1Dim(thermit->first)){
			  domain->setLocalSummv2_1Dim(summv2_1Dim[thermit->first], thermit->first);
			}
		}

		this->_state = STATE_POST_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition2to3(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition3to1(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		this->_state = STATE_NEW_TIMESTEP;
	}
	else {
		global_log->error() << "Leapfrog::transition3to1(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::accelerateUniformly(ParticleContainer* molCont, Domain* domain) {
	map<unsigned, double>* additionalAcceleration = domain->getPG()->getUAA();
	vector<Component> comp = *(_simulation.getEnsemble()->components());
	vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	for (compit = comp.begin(); compit != comp.end(); compit++) {
		unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
		if (cosetid != 0)
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = _timestepLength * additionalAcceleration[d][cosetid];
		else
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = 0;
	}

	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::accelerateInstantaneously(DomainDecompBase* domainDecomp, ParticleContainer* molCont, Domain* domain) {
	//vector<Component> comp = *(_simulation.getEnsemble()->components());
	//vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	string moved ("moved");
	unsigned cid_moved =  domain->getPG()->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
	
	//calculates globalN and globalVelocitySum for getMissingVelocity(); necessary!
	domain->getPG()->prepare_getMissingVelocity(domainDecomp, molCont, cid_moved, domain->getNumberOfComponents());
	// Just the velocity in x-direction (horizontally) is regulated to a constant value; in y and z the velocity is unrestricted 			
	for (int d = 0; d < 3; d++)
	  componentwiseVelocityDelta[d][cid_moved] = domain->getPG()->getMissingVelocity(cid_moved, d);
	
	for(int d = 0; d < 3; d++)
	  domain->getPG()->addGlobalVelSumBeforeAcc(d, cid_moved, domain->getPG()->getGlobalVelSum(d, cid_moved) / domain->getPG()->getGlobalN(cid_moved));
	  
	for (Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
	// Control of velocity after artificial acceleration
	domain->getPG()->prepare_getMissingVelocity(domainDecomp, molCont, cid_moved, domain->getNumberOfComponents());

	for(int d = 0; d < 3; d++)
	  domain->getPG()->addGlobalVelSumAfterAcc(d, cid_moved, domain->getPG()->getGlobalVelSum(d,cid_moved) / domain->getPG()->getGlobalN(cid_moved));
}

void Leapfrog::init1D(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->setXY();
}

void Leapfrog::zOscillation(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		/* TODO: use cid instead of complicated id + tersoff */
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->resetXY();
}
