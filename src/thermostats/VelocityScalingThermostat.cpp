#include "thermostats/VelocityScalingThermostat.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"

#include <iomanip>

using namespace std;
using Log::global_log;

VelocityScalingThermostat::VelocityScalingThermostat() : _globalBetaTrans(1), _globalBetaRot(1), _globalVelocity(NULL), _componentwise(false) {
	_globalVelocity = new double[3];
}

VelocityScalingThermostat::~VelocityScalingThermostat() {
	delete[] _globalVelocity;
}

void VelocityScalingThermostat::setBetaTrans(int componentId, double beta) {
	_componentBetaTrans[componentId] = beta;
}

void VelocityScalingThermostat::setBetaRot(int componentId, double beta) {
	_componentBetaRot[componentId] = beta;
}

void VelocityScalingThermostat::setGlobalVelocity(double v[3]) {
	for(int d = 0; d < 3; d++) {
		_globalVelocity[d] = v[d];
	}
}
void VelocityScalingThermostat::setVelocity(int componentId, double v[3]) {
	if( _componentVelocity.find(componentId) == _componentVelocity.end() ) {
		_componentVelocity[componentId] = new double[3];
	}
	for(int d = 0; d < 3; d++) {
		_componentVelocity[componentId][d] = v[d];
	}
}

void VelocityScalingThermostat::apply(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
	Molecule *molecule;
	
	if(_componentwise ) {
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			int thermostatId;
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			int cid = molecule->componentid();
			double alphaTrans = _globalAlphaTrans;
			thermostatId = _simulation.getDomain()->getThermostat(cid);
			betaTrans = _componentBetaTrans[thermostatId];
			betaRot   = _componentBetaRot[thermostatId];
			alphaTrans = _componentAlphaTrans[thermostatId];
			map<int, int> Dim = _simulation.getDomain()->getDim();
			map<int, double*>::iterator vIter;
			if( (vIter = _componentVelocity.find(thermostatId)) != _componentVelocity.end()) {
			      	if(_simulation.getDomain()->isScaling1Dim(thermostatId) && _simulation.getDomain()->getAlphaTransCorrection(thermostatId) == false){
				      molecule->scale_v_1Dim(alphaTrans, Dim[thermostatId]);
				}
				else{
				      molecule->scale_v(betaTrans);
				}
			}
			else {
				double *v = vIter->second;
				molecule->vsub(v[0], v[1], v[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(v[0], v[1], v[2]);
			}
			molecule->scale_D(betaRot);
		}
		// Calculate directed velocities in control volumes of dimension [1 * 1 * Z_max]
		calculateDirectedVelocities(moleculeContainer, dode);
		
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {	
			// virial calculation (kinetic and force part) for stresses in solids
			if(_simulation.getDomain()->isStressCalculating(molecule->componentid()) && !(_simulation.getDomain()->getSimstep() % _simulation.getDomain()->getStressRecordTimeStep())){
			  for(int d = 0; d < 3; d++){
			    long double v_d = molecule->v(d) - molecule->getDirectedVelocity(d);
			    for(int e = 0; e < 3; e++){
			      long double v_e = molecule->v(e) - molecule->getDirectedVelocity(d);
			      long double kinTrans = v_d * v_e * molecule->mass();
			      molecule->addVirialKin(d, e, kinTrans);
			    }
			  }
			}
			if(_simulation.getDomain()->isBulkPressure(molecule->componentid())){
			  if((molecule->r(0) >= _simulation.getDomain()->getBulkBoundary(0)-_simulation.getLJCutoff() && molecule->r(0) <= _simulation.getDomain()->getBulkBoundary(1)+_simulation.getLJCutoff() && molecule->r(1) >= _simulation.getDomain()->getBulkBoundary(2)-_simulation.getLJCutoff() && molecule->r(1) <= _simulation.getDomain()->getBulkBoundary(3)+_simulation.getLJCutoff())){
			    for (unsigned short d = 0; d < 3; ++d){
				long double v_d = molecule->v(d) - molecule->getDirectedVelocity(d);
				long double kinTrans = v_d * v_d * molecule->mass();
				molecule->addPressureKin(d, kinTrans);
			    }
			  }
			}
			if(_simulation.getDomain()->isBarostat(molecule->componentid()) && ((_simulation.getDomain()->getSimstep() >= _simulation.getBarostatTimeInit() && _simulation.getDomain()->getSimstep() <= _simulation.getBarostatTimeEnd()))){
			  if((molecule->r(0) >= _simulation.getDomain()->getControl_bottom(0) && molecule->r(0) <= _simulation.getDomain()->getControl_top(0) && molecule->r(1) >= _simulation.getDomain()->getControl_bottom(1) && molecule->r(1) <= _simulation.getDomain()->getControl_top(1))){
			    for (unsigned short d = 0; d < 3; ++d){
				long double v_d = molecule->v(d) - molecule->getDirectedVelocity(d);
				long double kinTrans = v_d * v_d * molecule->mass();
				molecule->addPressureKin_barostat(d, kinTrans);
			    }
			  }
			}
			if(_simulation.getDomain()->isConfinement(molecule->componentid())){
			  if(molecule->r(0) >= _simulation.getDomain()->getConfinementEdge(0) && molecule->r(0) <= _simulation.getDomain()->getConfinementEdge(1) && molecule->r(1) >= _simulation.getDomain()->get_confinementMidPoint(3)-_simulation.getDomain()->getConfinementEdge(5)-0.755 && molecule->r(1) <= _simulation.getDomain()->get_confinementMidPoint(1)+_simulation.getDomain()->getConfinementEdge(5)+0.755 && !(_simulation.getDomain()->getSimstep() % _simulation.getDomain()->getConfinementRecordTimeStep())){
			    for (unsigned short d = 0; d < 3; ++d){
				long double v_d = molecule->v(d) - molecule->getDirectedVelocity(d);
				long double kinTrans = v_d * v_d * molecule->mass();
				molecule->addPressureKinConfinement(d, kinTrans);
				for(int e = 0; e < 3; e++){
				  long double v_e = molecule->v(e) - molecule->getDirectedVelocity(d);
				  long double kinTrans = v_d * v_e * molecule->mass();
				  molecule->addVirialKinConfinement(d, e, kinTrans);
				}
			    }
			  }
			}

		}
	}
	else {
		double betaTrans = _globalBetaTrans;
		double betaRot = _globalBetaRot;
		global_log->debug() << "Beta rot: " << betaRot << endl;
		global_log->debug() << "Beta trans: " << betaTrans << endl;
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			if(_globalVelocity == NULL) {
				molecule->scale_v(betaTrans);
			}
			else {
				molecule->vsub(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
			}
			molecule->scale_D(betaRot);
		}
	}
}

// Calculates directed velocities for control volumes of dimension [1 * 1 * Z_max] averaged for each timestep
void VelocityScalingThermostat::calculateDirectedVelocities(ParticleContainer *moleculeContainer, DomainDecompBase *dode) {
   double xMax = _simulation.getDomain()->getGlobalLength(0);
   double yMax = _simulation.getDomain()->getGlobalLength(1);
   double roundX, roundY;
   unsigned CID;

   // initialize molecule counting; number of molecules of a certain component-id in each control volume
   for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
     for (double xCoord = 0.5; xCoord < xMax; xCoord++)
       for (double yCoord = 0.5; yCoord < yMax; yCoord++){
	  _countedMolecules[cid][xCoord][yCoord] = 0;
	  for (unsigned d = 0; d < 3; d++){
	    _directedVelocity[cid][xCoord][yCoord][d] = 0.0;
	    _universalDirectedVelocity[cid][xCoord][yCoord][d] = 0.0;
	  }
       }
   }
   Molecule *molecule;

   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
      roundX = floor(molecule->r(0)) + 0.5;
      roundY = floor(molecule->r(1)) + 0.5;
      CID = molecule->componentid();
      _countedMolecules[CID][roundX][roundY]++;
      
      for (unsigned d = 0; d < 3; d++)
	_directedVelocity[CID][roundX][roundY][d] += molecule->v(d);
   }
     
   //Parallelization
   for (unsigned cid = 0; cid < _simulation.getDomain()->getNumberOfComponents(); cid++){
      dode->collCommInit(4*round(xMax)*round(yMax));
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    dode->collCommAppendUnsLong(_countedMolecules[cid][x][y]);
	    for(int d = 0; d < 3; d++)
		dode->collCommAppendDouble(_directedVelocity[cid][x][y][d]);
	  }
      dode->collCommAllreduceSum();
	for(double x = 0.5; x < xMax; x++)
	  for(double y = 0.5; y < yMax; y++){
	    _universalCountedMolecules[cid][x][y] = dode->collCommGetUnsLong();
	    for(int d = 0; d < 3; d++){
		_universalDirectedVelocity[cid][x][y][d] = dode->collCommGetDouble();
		if (_universalCountedMolecules[cid][x][y] != 0)
		  _universalDirectedVelocity[cid][x][y][d] = _universalDirectedVelocity[cid][x][y][d]/_universalCountedMolecules[cid][x][y];
		else
		  _universalDirectedVelocity[cid][x][y][d] = 0.0;
	    }
	  }
      dode->collCommFinalize();
   }
   
   
   // attachment of the _universalDirectedVelocity to each molecule in dependence on the control volume part of
   for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()){
     double x = floor(molecule->r(0)) + 0.5;
     double y = floor(molecule->r(1)) + 0.5;
     unsigned cid = molecule->componentid();
     for (int d=0; d<3; d++){
	molecule->setDirectedVelocity(d, _universalDirectedVelocity[cid][x][y][d]); 
     }
   }
}
