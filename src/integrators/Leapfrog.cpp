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
		unsigned int cid;
		int thermostat;
		unsigned molID;

		string moved ("moved");
		string fixed ("fixed");
		unsigned cid_moved =  domain->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
		unsigned cid_fixed =  domain->getCidMovement(fixed, domain->getNumberOfComponents()) - 1;

		for (tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()) {
		  molID = tempMolecule->id();	
		  cid = tempMolecule->componentid();
		  thermostat = domain->getThermostat(cid);
		  if(domain->isThermostatLayer() == true || domain->isThermostatWallLayer() == true){
			thermostat = domain->moleculeInLayer(tempMolecule->r(0), tempMolecule->r(1), tempMolecule->r(2), tempMolecule->componentid());
		  }
		  if (domain->isFixedEvery() == true && molID%domain->getFixEvery() == 0 && cid == cid_fixed){
				for(unsigned short d=0;d<3;++d) tempMolecule->setv(d, 0.0);	
		  }else if (domain->isFixedRegion() == true && cid == cid_fixed
				&& tempMolecule->r(0) > domain->getFixedRegion(0) && tempMolecule->r(0) < domain->getFixedRegion(1) 
				&& tempMolecule->r(1) > domain->getFixedRegion(2) && tempMolecule->r(1) < domain->getFixedRegion(3)
				&& tempMolecule->r(2) > domain->getFixedRegion(4) && tempMolecule->r(2) < domain->getFixedRegion(5)){ 
				for(unsigned short d=0;d<3;++d) tempMolecule->setv(d, 0.0);	
		  }else if ((domain->isFixedEvery() == false && domain->isFixedRegion() == false && molID%500 == 0 && cid == cid_fixed) 
				|| (_simulation.getSimulationStep() <= _simulation.getInitStatistics() && (molID%500 == 0 && cid == cid_moved))){
				for(unsigned short d=0;d<3;++d) tempMolecule->setv(d, 0.0);		
		  }else {				
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
		unsigned int cid;
		int molID;
		
		string moved ("moved");
		string fixed ("fixed");
		unsigned cid_moved =  domain->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
		unsigned cid_fixed =  domain->getCidMovement(fixed, domain->getNumberOfComponents()) - 1;

		double dt_half = 0.5 * this->_timestepLength;
		if (domain->severalThermostats()) {
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				molID = tM->id();
				cid = tM->componentid();
				int thermostat = domain->getThermostat(cid);
				if (domain->isFixedEvery() == true && molID%domain->getFixEvery() == 0 && cid == cid_fixed){ 
					for(unsigned short d=0;d<3;++d) tM->setv(d, 0.0); 
				}else if (domain->isFixedRegion() == true && cid == cid_fixed
					&& tM->r(0) > domain->getFixedRegion(0) && tM->r(0) < domain->getFixedRegion(1) 
					&& tM->r(1) > domain->getFixedRegion(2) && tM->r(1) < domain->getFixedRegion(3)
					&& tM->r(2) > domain->getFixedRegion(4) && tM->r(2) < domain->getFixedRegion(5)){ 
					for(unsigned short d=0;d<3;++d) tM->setv(d, 0.0); 
				}else if ((domain->isFixedEvery() == false && domain->isFixedRegion() == false && molID%500 == 0 
					&& cid == cid_fixed) || (_simulation.getSimulationStep() <= _simulation.getInitStatistics() 
					&& (molID%500 == 0 && cid == cid_moved))){
					for(unsigned short d=0;d<3;++d) tM->setv(d, 0.0);
				}else {
					if(domain->isThermostatLayer() == true || domain->isThermostatWallLayer() == true){
						thermostat = domain->moleculeInLayer(tM->r(0), tM->r(1), tM->r(2), tM->componentid());
					}
					tM->upd_postF(dt_half, summv2[thermostat], summv2_1Dim[thermostat], sumIw2[thermostat], domain);
				 
					N[thermostat]++;
					rotDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
					if (domain->isSpringDamped() == true && tM->componentid() == cid_moved &&  domain->getSpringConst() == 0 
						&& tM->r(1) > domain->getAverageY()-1.5){
						N[thermostat]--;
						rotDOF[thermostat] -= tM->component()->getRotationalDegreesOfFreedom();
					}
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
				molID = tM->id();
				cid = tM->componentid();
				if (domain->isFixedEvery() == true && molID%domain->getFixEvery() == 0 && cid == cid_fixed){ }
				else if (domain->isFixedRegion() == true && cid == cid_fixed
					&& tM->r(0) > domain->getFixedRegion(0) && tM->r(0) < domain->getFixedRegion(1) 
					&& tM->r(1) > domain->getFixedRegion(2) && tM->r(1) < domain->getFixedRegion(3)
					&& tM->r(2) > domain->getFixedRegion(4) && tM->r(2) < domain->getFixedRegion(5)){ }
				else if ((domain->isFixedEvery() == false && domain->isFixedRegion() == false && molID%500 == 0 
					&& cid == cid_fixed) || (_simulation.getSimulationStep() <= _simulation.getInitStatistics() 
					&& (molID%500 == 0 && cid == cid_moved))){ }
				else {
					tM->upd_postF(dt_half, summv2gt, summv2_1Dimgt, sumIw2gt, domain);
					assert(summv2gt >= 0.0);
					Ngt++;
					rotDOFgt += tM->component()->getRotationalDegreesOfFreedom();
					if (domain->isSpringDamped() == true && tM->componentid() == cid_moved &&  domain->getSpringConst() == 0 
						&& tM->r(1) > domain->getAverageY()-1.5){
						Ngt--;
						rotDOFgt -= tM->component()->getRotationalDegreesOfFreedom();
					}
				}
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

// accelerates the component that is designated to be "moved" instantaneously (without ramp profile)
// and in each time step to a certain velocity 
void Leapfrog::accelerateInstantaneously(DomainDecompBase* domainDecomp, ParticleContainer* molCont, Domain* domain) {
	map<unsigned, double> componentwiseVelocityDelta[3];
	string moved ("moved");
	unsigned cid_moved =  domain->getCidMovement(moved, domain->getNumberOfComponents()) - 1;
	//calculates globalN and globalVelocitySum for getMissingVelocity(); necessary!
	domain->getPG()->prepare_getMissingVelocity(domainDecomp, molCont, cid_moved, domain->getNumberOfComponents(), _simulation.getDirectedVelocityTime());
	
	// time span in which the acceleration or target velocity is gradually increased to its final value
	int influencingTime = _simulation.getSimulationStep() -_simulation.getInitStatistics() - 25000;
	// asymptotic value for the influenceFactor
	double influenceFactor = 0.1265769815;
	
	// during influencingTime span the influence of the acceleration delta_v is slowly decreased from 1 to 0.1265769815
	// ---> at the beginning the influence of the acceleration is a lot more intense than later ---> increasing stability!
	if(influencingTime > 0 && abs(influencingTime) < _simulation.getDirectedVelocityTime())
		influenceFactor = 1.0 - exp((-1)*(double)_simulation.getDirectedVelocityTime()/((double)influencingTime*exp(2)));
	else if(influencingTime <= 0)
		influenceFactor = 1.0;
        
	// returns the value for delta_v = v_target - 2*v_directed,currently + v_directed,average
	for (int d = 0; d < 3; d++)
	  componentwiseVelocityDelta[d][cid_moved] = domain->getPG()->getMissingVelocity(cid_moved, d, _simulation.getSimulationStep(), _simulation.getInitStatistics());
	
	for(int d = 0; d < 3; d++)
	  domain->getPG()->addGlobalVelSumBeforeAcc(d, cid_moved, domain->getPG()->getGlobalVelSum(d, cid_moved) / domain->getPG()->getGlobalN(cid_moved));
	
	// acceleration by adding delta_v to all molecules of the component
	for (Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(influenceFactor*componentwiseVelocityDelta[0][cid],
		              influenceFactor*componentwiseVelocityDelta[1][cid],
		              influenceFactor*componentwiseVelocityDelta[2][cid]);
	}
	// Control of velocity after artificial acceleration
	domain->getPG()->prepare_getMissingVelocity(domainDecomp, molCont, cid_moved, domain->getNumberOfComponents(), _simulation.getDirectedVelocityTime());

	for(int d = 0; d < 3; d++)
	  domain->getPG()->addGlobalVelSumAfterAcc(d, cid_moved, domain->getPG()->getGlobalVelSum(d,cid_moved) / domain->getPG()->getGlobalN(cid_moved));
}

// keeps the velocity of the system at its box margins in y-direction at v=0, while it accelerates the system
// in the middle at y=0.5*y_max to a certain velocity to maintain the target shear rate
void Leapfrog::shearRate(DomainDecompBase* domainDecomp, ParticleContainer* molCont, Domain* domain) {
	double shearVelocityDelta, shearVelocityTarget, directedVel, directedVelAverage;
	double shearRateBox[4];
	// geometrical dimensions of the box where shear is applied
	// x_min
	shearRateBox[0] = domain->getPG()->getShearRateBox(0);
	// x_max
	shearRateBox[1] = domain->getPG()->getShearRateBox(1);
	// y_min
	shearRateBox[2] = domain->getPG()->getShearRateBox(2);
	// y_max
	shearRateBox[3] = domain->getPG()->getShearRateBox(3);
	// half width of the stripe in which the target velocity is applied
	double shearWidth = domain->getPG()->getShearWidth();
	double shearYmax = shearRateBox[3] - shearRateBox[2];
	double shearRate = domain->getPG()->getShearRate();
	double slowAcceleration = 1.0;
	unsigned cid = domain->getPG()->getShearComp();
	unsigned yun;
	// time span in which the acceleration or target velocity is gradually increased to its final value
	int influencingTime = _simulation.getSimulationStep() -_simulation.getInitStatistics() - domain->getPG()->getShearRampTime();
	// asymptotic value for the influenceFactor
	double influenceFactor = 0.1265769815;
	string Confinement ("Confinement");
	string Default ("Default");
	string noDirVel ("noDirVel");
	double mv2 = 0.0;
	double Iw2 = 0.0;
	
	// slowAcceleration = [0;1]
	if((_simulation.getSimulationStep()-_simulation.getInitStatistics()) < domain->getPG()->getShearRampTime())
			slowAcceleration = (double)(_simulation.getSimulationStep()-_simulation.getInitStatistics())/domain->getPG()->getShearRampTime();
	// during influencingTime span the influence of the acceleration delta_v is slowly decreased from 1 to 0.1265769815
	// ---> at the beginning the influence of the acceleration is a lot more intense than later ---> increasing stability!
	if(influencingTime > 0 && abs(influencingTime) < _simulation.getDirectedVelocityTime())
		influenceFactor = 1.0 - exp((-1)*(double)_simulation.getDirectedVelocityTime()/((double)influencingTime*exp(2)));
	else if(influencingTime <= 0)
		influenceFactor = 1.0;
	//calculates _globalShearN and _globalShearVelocitySum for getDirectedShearVel(); necessary!
	domain->getPG()->prepareShearRate(molCont, domainDecomp, _simulation.getDirectedVelocityTime());

	for (Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		if(thismol->componentid() == cid && thismol->r(0) >= shearRateBox[0] && thismol->r(0) <= shearRateBox[1] && thismol->r(1) >= shearRateBox[2] && thismol->r(1) <= shearRateBox[2] + shearWidth){
			yun = 0;
		}else if(thismol->componentid() == cid && thismol->r(0) >= shearRateBox[0] && thismol->r(0) <= shearRateBox[1] && thismol->r(1) >= 0.5 * (shearRateBox[3]-shearRateBox[2]) - shearWidth && thismol->r(1) < 0.5 * (shearRateBox[3]-shearRateBox[2])){
			yun = 1;
		}else if(thismol->componentid() == cid && thismol->r(0) >= shearRateBox[0] && thismol->r(0) <= shearRateBox[1] && thismol->r(1) >= 0.5 * (shearRateBox[3]-shearRateBox[2]) && thismol->r(1) <= 0.5 * (shearRateBox[3]-shearRateBox[2]) + shearWidth){
			yun = 2;	
		}else if(thismol->componentid() == cid && thismol->r(0) >= shearRateBox[0] && thismol->r(0) <= shearRateBox[1] && thismol->r(1) >= shearRateBox[3] - shearWidth && thismol->r(1) <= shearRateBox[3]){
			yun = 3;
		}else
			yun = 4;
		if((_simulation.isShearRate() && yun < 4) || (_simulation.isShearForce() && (yun == 0 || yun == 3))){
		  directedVel = domain->getPG()->getDirectedShearVel(yun);
		  directedVelAverage = domain->getPG()->getDirectedShearVelAverage(yun);
		  if(_simulation.isShearRate()){
			shearVelocityTarget = shearRate * shearYmax/2 - fabs((shearYmax/2 - thismol->r(1)) * shearRate);
			shearVelocityTarget = slowAcceleration*shearVelocityTarget;
		  }else
			shearVelocityTarget = 0.0;  
		  
		  // returns the value for delta_v = v_target - 2*v_directed,currently + v_directed,average
		  shearVelocityDelta = shearVelocityTarget - 2*directedVel + 1*directedVelAverage;
		  
		  // damping of the influence as 1 - exp(A/t)
		  if(influencingTime > 0)
				shearVelocityDelta = influenceFactor*shearVelocityDelta;
		  
		  if(_simulation.isShearRate()){
			// adaption of the directed velocities
			thismol->setDirectedVelocity(0, shearVelocityTarget); 
			thismol->setDirectedVelocitySlab(0, shearVelocityTarget);
			thismol->setDirectedVelocityStress(0, shearVelocityTarget); 
			thismol->setDirectedVelocityConfinement(0, shearVelocityTarget);
			thismol->setDirectedVelocity(1, 0.0); 
			thismol->setDirectedVelocitySlab(1, 0.0);
			thismol->setDirectedVelocityStress(1, 0.0); 
			thismol->setDirectedVelocityConfinement(1, 0.0);
			thismol->setDirectedVelocity(2, 0.0); 
			thismol->setDirectedVelocitySlab(2, 0.0);
			thismol->setDirectedVelocityStress(2, 0.0); 
			thismol->setDirectedVelocityConfinement(2, 0.0);
		  
			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2, noDirVel);
			
			// calculation of the energetical influence of the shear
			domain->addPreShearEnergyConf(mv2);

			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2, Default);
			// calculation of the energetical influence of the shear
			domain->addPreShearEnergyDefault(mv2);
		  }
		  // acceleration!
		  thismol->vadd(shearVelocityDelta,0.0,0.0);
		  domain->addShearVelDelta(shearVelocityDelta);
		  
		  if(_simulation.isShearRate()){
			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2, noDirVel);
			// calculation of the energetical influence of the shear
			domain->addPostShearEnergyConf(mv2);

			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2, Default);
			// calculation of the energetical influence of the shear
			domain->addPostShearEnergyDefault(mv2);
		  }
		}
	}
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
